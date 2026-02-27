#!/usr/bin/env python3
"""
Cell Segmentation using Edge-based Seeded Watershed

This script segments plant cells from CFW (cell wall) channel using edge detection
and seeded watershed algorithm.

Author: smFISH Analysis Pipeline
Version: 1.0.0
"""

import argparse
import json
import numpy as np
from skimage import io, filters, morphology, feature, measure, segmentation
from scipy import ndimage
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================================
# PARAMETER DEFINITIONS WITH DETAILED EXPLANATIONS
# ============================================================================

PARAM_DEFINITIONS = {
    'pixel_size_xy': {
        'type': float,
        'default': 0.44,
        'unit': 'μm/pixel',
        'range': (0.01, 10.0),
        'description': 'Physical size of one pixel in XY plane',
        'how_to_determine': '''
            1. Check microscope metadata (usually in CZI file)
            2. Use objective magnification and camera pixel size
            3. Formula: pixel_size = (camera_pixel_size) / (magnification * tube_lens_factor)
            4. Typical values:
               - 100x objective: 0.06-0.10 μm/pixel
               - 63x objective: 0.10-0.16 μm/pixel
               - 40x objective: 0.16-0.25 μm/pixel
               - 20x objective: 0.30-0.50 μm/pixel
               - 10x objective: 0.60-1.00 μm/pixel
        ''',
        'impact': 'Affects all size-based parameters (min/max cell size, nucleus diameter)'
    },
    
    'pixel_size_z': {
        'type': float,
        'default': 0.50,
        'unit': 'μm/pixel',
        'range': (0.01, 10.0),
        'description': 'Physical size of one Z-step',
        'how_to_determine': '''
            1. Check Z-stack settings in microscope software
            2. Usually 0.3-1.0 μm for confocal microscopy
            3. Should be ≤ 2 * (pixel_size_xy) for isotropic resolution
            4. Typical values:
               - High resolution: 0.2-0.3 μm
               - Standard: 0.4-0.6 μm
               - Fast acquisition: 0.8-1.0 μm
        ''',
        'impact': 'Affects 3D volume calculations and Z-direction smoothing'
    },
    
    'canny_sigma': {
        'type': float,
        'default': 2.0,
        'unit': 'pixels',
        'range': (0.5, 10.0),
        'description': 'Gaussian smoothing sigma before edge detection',
        'how_to_determine': '''
            1. Start with 2.0 for typical images
            2. Increase (3-5) if image is noisy → smoother edges
            3. Decrease (1-1.5) if cell walls are thin → sharper edges
            4. Rule of thumb: sigma ≈ (cell_wall_thickness_in_pixels) / 2
            5. Visualize: check edge detection output, edges should be continuous
        ''',
        'impact': 'Higher values → smoother but less detailed edges; Lower values → more detailed but noisier edges'
    },
    
    'canny_low': {
        'type': float,
        'default': 10,
        'unit': 'intensity',
        'range': (1, 100),
        'description': 'Lower threshold for Canny edge detection',
        'how_to_determine': '''
            1. Should be lower than canny_high (typically 1/3 to 1/2 of canny_high)
            2. Start with 10 for 8-bit images (0-255 range)
            3. Increase if too many false edges (background noise)
            4. Decrease if missing weak cell wall edges
            5. Check histogram: set to ~5-10% of peak intensity
        ''',
        'impact': 'Lower values → more edges detected (including weak ones); Higher values → fewer edges (only strong ones)'
    },
    
    'canny_high': {
        'type': float,
        'default': 30,
        'unit': 'intensity',
        'range': (10, 200),
        'description': 'Upper threshold for Canny edge detection',
        'how_to_determine': '''
            1. Should be 2-3x of canny_low
            2. Start with 30 for 8-bit images
            3. Increase if edges are broken/incomplete
            4. Decrease if too many edges (over-segmentation)
            5. Check histogram: set to ~15-30% of peak intensity
        ''',
        'impact': 'Higher values → only strong edges; Lower values → more edges including weaker ones'
    },
    
    'min_distance': {
        'type': int,
        'default': 5,
        'unit': 'pixels',
        'range': (1, 50),
        'description': 'Minimum distance between cell centers (seed points)',
        'how_to_determine': '''
            1. Should be smaller than typical cell diameter
            2. Formula: min_distance ≈ (min_cell_diameter_μm) / (2 * pixel_size_xy)
            3. For cells of 10-20 μm diameter with 0.44 μm/pixel:
               min_distance ≈ 10 / (2 * 0.44) ≈ 11 pixels
            4. Start with 5-10 and adjust:
               - Too small → over-segmentation (one cell split into multiple)
               - Too large → under-segmentation (multiple cells merged)
        ''',
        'impact': 'Smaller values → more seed points → more cells detected; Larger values → fewer seed points → cells may be merged'
    },
    
    'min_cell_size': {
        'type': int,
        'default': 2000,
        'unit': 'voxels (3D pixels)',
        'range': (100, 100000),
        'description': 'Minimum cell volume to keep',
        'how_to_determine': '''
            1. Calculate expected cell volume:
               volume_voxels = (cell_diameter_μm)³ / (pixel_size_xy² * pixel_size_z)
            2. Example: 10 μm diameter cell
               volume = 10³ / (0.44² * 0.5) ≈ 10,300 voxels
            3. Set min_cell_size to 20-30% of expected volume to filter debris
            4. Check output: if losing real cells, decrease this value
        ''',
        'impact': 'Filters out small objects (debris, noise). Too high → lose small cells; Too low → keep noise'
    },
    
    'max_cell_size': {
        'type': int,
        'default': 50000,
        'unit': 'voxels (3D pixels)',
        'range': (1000, 1000000),
        'description': 'Maximum cell volume to keep',
        'how_to_determine': '''
            1. Calculate expected maximum cell volume (see min_cell_size)
            2. Set to 3-5x of typical cell volume to filter merged cells
            3. Example: if typical cell is 10,000 voxels, set max to 30,000-50,000
            4. Check output: if losing large cells, increase this value
        ''',
        'impact': 'Filters out very large objects (merged cells, artifacts). Too low → lose large cells; Too high → keep merged cells'
    },
    
    'roi_y_min': {
        'type': int,
        'default': 256,
        'unit': 'pixels',
        'range': (0, 4096),
        'description': 'Minimum Y coordinate for region of interest (optional)',
        'how_to_determine': '''
            1. Use 0 to process entire image
            2. Set to non-zero to focus on specific region (e.g., lower half)
            3. Visualize image first to identify cell-containing regions
            4. Example: if cells are in lower half of 512px image, set roi_y_min=256
        ''',
        'impact': 'Restricts analysis to specific region, ignoring areas without cells'
    }
}


def print_parameter_guide():
    """Print detailed parameter guide"""
    print("\n" + "="*80)
    print("PARAMETER GUIDE FOR CELL SEGMENTATION")
    print("="*80)
    
    for param_name, param_info in PARAM_DEFINITIONS.items():
        print(f"\n{param_name.upper()}")
        print("-" * 40)
        print(f"Type: {param_info['type'].__name__}")
        print(f"Default: {param_info['default']} {param_info['unit']}")
        print(f"Range: {param_info['range']}")
        print(f"Description: {param_info['description']}")
        print(f"\nHow to determine:")
        print(param_info['how_to_determine'].strip())
        print(f"\nImpact: {param_info['impact']}")
    
    print("\n" + "="*80)


def segment_cells(
    cfw_image_path,
    output_dir,
    pixel_size_xy=0.44,
    pixel_size_z=0.50,
    canny_sigma=2.0,
    canny_low=10,
    canny_high=30,
    min_distance=5,
    min_cell_size=2000,
    max_cell_size=50000,
    roi_y_min=256
):
    """
    Segment cells from CFW channel using edge-based seeded watershed.
    
    See PARAM_DEFINITIONS for detailed parameter explanations.
    """
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Loading CFW image from {cfw_image_path}...")
    cfw_image = io.imread(cfw_image_path)
    
    # Apply ROI if specified
    if roi_y_min > 0:
        print(f"Applying ROI: Y >= {roi_y_min}")
        cfw_image = cfw_image[:, roi_y_min:, :]
    
    print(f"Image shape: {cfw_image.shape}")
    print(f"Image dtype: {cfw_image.dtype}, range: {cfw_image.min()}-{cfw_image.max()}")
    
    # Edge detection
    print(f"\nPerforming edge detection (sigma={canny_sigma}, low={canny_low}, high={canny_high})...")
    edges_3d = np.zeros_like(cfw_image, dtype=bool)
    for z in range(cfw_image.shape[0]):
        edges_3d[z] = feature.canny(
            cfw_image[z],
            sigma=canny_sigma,
            low_threshold=canny_low,
            high_threshold=canny_high
        )
    
    # Cell interior (inverse of edges)
    cell_interior = ~edges_3d
    
    # Distance transform
    print("Computing distance transform...")
    distance = ndimage.distance_transform_edt(cell_interior)
    
    # Find local maxima as seeds
    print(f"Finding seed points (min_distance={min_distance})...")
    from skimage.feature import peak_local_max
    local_max_coords = peak_local_max(
        distance,
        min_distance=min_distance,
        threshold_abs=1,
        exclude_border=False
    )
    
    print(f"Found {len(local_max_coords)} seed points")
    
    # Create markers
    markers = np.zeros(distance.shape, dtype=int)
    for i, coord in enumerate(local_max_coords):
        markers[tuple(coord)] = i + 1
    
    # Watershed segmentation
    print("Performing watershed segmentation...")
    labels = segmentation.watershed(-distance, markers, mask=cell_interior)
    
    # Filter by size
    print(f"Filtering cells by size ({min_cell_size}-{max_cell_size} voxels)...")
    props = measure.regionprops(labels)
    filtered_labels = np.zeros_like(labels)
    cell_id = 1
    cell_info = []
    
    for prop in props:
        if min_cell_size <= prop.area <= max_cell_size:
            filtered_labels[labels == prop.label] = cell_id
            
            # Calculate physical dimensions
            bbox = prop.bbox
            size_z = (bbox[3] - bbox[0]) * pixel_size_z
            size_y = (bbox[4] - bbox[1]) * pixel_size_xy
            size_x = (bbox[5] - bbox[2]) * pixel_size_xy
            
            cell_info.append({
                'cell_id': cell_id,
                'volume_voxels': int(prop.area),
                'centroid_z': float(prop.centroid[0]),
                'centroid_y': float(prop.centroid[1]),
                'centroid_x': float(prop.centroid[2]),
                'size_z_um': round(size_z, 2),
                'size_y_um': round(size_y, 2),
                'size_x_um': round(size_x, 2)
            })
            cell_id += 1
    
    print(f"Kept {len(cell_info)} cells after filtering")
    
    # Save results
    print("\nSaving results...")
    io.imsave(str(output_dir / 'cell_labels.tif'), filtered_labels.astype(np.uint16), check_contrast=False)
    
    with open(output_dir / 'cell_info.json', 'w') as f:
        json.dump(cell_info, f, indent=2)
    
    # Generate report
    generate_report(cfw_image, filtered_labels, cell_info, output_dir)
    
    print(f"\n✓ Cell segmentation complete!")
    print(f"  Output directory: {output_dir}")
    print(f"  Cells detected: {len(cell_info)}")
    
    return filtered_labels, cell_info


def generate_report(cfw_image, labels, cell_info, output_dir):
    """Generate visualization report"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    
    # MIP of CFW
    axes[0, 0].imshow(cfw_image.max(axis=0), cmap='gray')
    axes[0, 0].set_title(f'CFW Signal (MIP)')
    axes[0, 0].axis('off')
    
    # MIP of labels
    label_mip = labels.max(axis=0)
    axes[0, 1].imshow(label_mip, cmap='nipy_spectral')
    axes[0, 1].set_title(f'Cell Labels (MIP, N={len(cell_info)})')
    axes[0, 1].axis('off')
    
    # Z=20 slice
    z_slice = min(20, cfw_image.shape[0] - 1)
    axes[1, 0].imshow(cfw_image[z_slice], cmap='gray')
    axes[1, 0].set_title(f'CFW Signal (Z={z_slice})')
    axes[1, 0].axis('off')
    
    axes[1, 1].imshow(labels[z_slice], cmap='nipy_spectral')
    axes[1, 1].set_title(f'Cell Labels (Z={z_slice})')
    axes[1, 1].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'segmentation_report.png', dpi=150, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Segment plant cells from CFW channel',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic usage with default parameters
  python 02_segment_cells.py input.tif output_dir/
  
  # Show parameter guide
  python 02_segment_cells.py --help-params
  
  # Custom parameters for noisy images
  python 02_segment_cells.py input.tif output_dir/ \\
    --canny-sigma 3.0 --canny-low 15 --canny-high 40
  
  # Focus on lower half of image
  python 02_segment_cells.py input.tif output_dir/ \\
    --roi-y-min 256
        '''
    )
    
    parser.add_argument('cfw_image', help='Path to CFW channel image (TIFF)')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('--help-params', action='store_true', help='Show detailed parameter guide')
    
    # Add all parameters with detailed help
    for param_name, param_info in PARAM_DEFINITIONS.items():
        arg_name = '--' + param_name.replace('_', '-')
        parser.add_argument(
            arg_name,
            type=param_info['type'],
            default=param_info['default'],
            help=f"{param_info['description']} (default: {param_info['default']} {param_info['unit']})"
        )
    
    args = parser.parse_args()
    
    if args.help_params:
        print_parameter_guide()
        return
    
    segment_cells(
        args.cfw_image,
        args.output_dir,
        pixel_size_xy=args.pixel_size_xy,
        pixel_size_z=args.pixel_size_z,
        canny_sigma=args.canny_sigma,
        canny_low=args.canny_low,
        canny_high=args.canny_high,
        min_distance=args.min_distance,
        min_cell_size=args.min_cell_size,
        max_cell_size=args.max_cell_size,
        roi_y_min=args.roi_y_min
    )


if __name__ == '__main__':
    main()
