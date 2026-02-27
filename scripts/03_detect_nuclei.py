#!/usr/bin/env python3
"""
Nucleus Detection with Morphological Constraints

This script detects spherical nuclei (DAPI signal) within segmented cells,
handling weak signals and filtering by sphericity.

Author: smFISH Analysis Pipeline
Version: 1.0.0
"""

import argparse
import json
import numpy as np
from skimage import io, filters, measure, morphology
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
        'how_to_determine': 'Same as in cell segmentation (see 02_segment_cells.py --help-params)',
        'impact': 'Must match the value used in cell segmentation'
    },
    
    'pixel_size_z': {
        'type': float,
        'default': 0.50,
        'unit': 'μm/pixel',
        'range': (0.01, 10.0),
        'description': 'Physical size of one Z-step',
        'how_to_determine': 'Same as in cell segmentation',
        'impact': 'Must match the value used in cell segmentation'
    },
    
    'nucleus_diameter_min': {
        'type': float,
        'default': 5.0,
        'unit': 'μm',
        'range': (1.0, 50.0),
        'description': 'Minimum nucleus diameter',
        'how_to_determine': '''
            1. Measure nucleus diameter in your images using ImageJ/Fiji
            2. For plant cells: typically 5-10 μm
            3. For animal cells: typically 8-15 μm
            4. Set to ~80% of smallest expected nucleus diameter
            5. Too small → detect noise as nuclei
            6. Too large → miss small nuclei
        ''',
        'impact': 'Filters out small objects. Decrease if missing small nuclei; increase if detecting noise'
    },
    
    'nucleus_diameter_max': {
        'type': float,
        'default': 10.0,
        'unit': 'μm',
        'range': (2.0, 100.0),
        'description': 'Maximum nucleus diameter',
        'how_to_determine': '''
            1. Measure largest nucleus in your images
            2. Set to ~120% of largest expected nucleus diameter
            3. Too small → miss large nuclei
            4. Too large → detect merged nuclei or artifacts
        ''',
        'impact': 'Filters out large objects. Increase if missing large nuclei; decrease if detecting merged nuclei'
    },
    
    'min_intensity': {
        'type': float,
        'default': 10,
        'unit': 'intensity',
        'range': (0, 255),
        'description': 'Minimum intensity threshold for nucleus signal',
        'how_to_determine': '''
            1. Open image in ImageJ/Fiji
            2. Draw ROI around a nucleus
            3. Measure mean intensity
            4. Set min_intensity to 50-70% of mean nucleus intensity
            5. For weak DAPI signals: use 10-20
            6. For strong DAPI signals: use 30-50
            7. Check histogram: set to background + 2*std
        ''',
        'impact': 'Lower values → detect weaker signals (more nuclei, more noise); Higher values → only strong signals'
    },
    
    'max_intensity': {
        'type': float,
        'default': 40,
        'unit': 'intensity',
        'range': (10, 255),
        'description': 'Maximum intensity threshold for nucleus signal',
        'how_to_determine': '''
            1. Measure mean intensity of nuclei (see min_intensity)
            2. Set to 120-150% of mean nucleus intensity
            3. For weak DAPI signals: use 30-50
            4. For strong DAPI signals: use 80-150
            5. Should be > min_intensity (typically 2-4x)
        ''',
        'impact': 'Higher values → include brighter regions; Lower values → only dim signals (useful for weak DAPI)'
    },
    
    'min_sphericity': {
        'type': float,
        'default': 0.5,
        'unit': 'dimensionless',
        'range': (0.0, 1.0),
        'description': 'Minimum sphericity to consider as nucleus',
        'how_to_determine': '''
            1. Sphericity = (surface area of sphere with same volume) / (actual surface area)
            2. Perfect sphere: sphericity = 1.0
            3. Elongated object: sphericity < 0.5
            4. For plant cell nuclei: use 0.4-0.6 (slightly elongated)
            5. For animal cell nuclei: use 0.6-0.8 (more spherical)
            6. Start with 0.5 and adjust:
               - Too high → miss elongated nuclei
               - Too low → detect non-spherical artifacts
        ''',
        'impact': 'Higher values → only spherical objects; Lower values → allow elongated objects'
    },
    
    'gaussian_sigma': {
        'type': float,
        'default': 2.0,
        'unit': 'pixels',
        'range': (0.5, 10.0),
        'description': 'Gaussian smoothing sigma for noise reduction',
        'how_to_determine': '''
            1. Start with 2.0 for typical images
            2. Increase (3-5) if image is very noisy
            3. Decrease (1-1.5) if nuclei are very small
            4. Rule of thumb: sigma ≈ (nucleus_diameter_pixels) / 5
            5. For nucleus diameter 10 μm with 0.44 μm/pixel:
               sigma ≈ (10/0.44) / 5 ≈ 4.5 pixels
        ''',
        'impact': 'Higher values → smoother, less noise, but may blur small nuclei; Lower values → preserve details but more noise'
    },
    
    'process_n_cells': {
        'type': int,
        'default': 0,
        'unit': 'cells',
        'range': (0, 10000),
        'description': 'Number of cells to process (0 = all)',
        'how_to_determine': '''
            1. Use 0 to process all cells (recommended for final analysis)
            2. Use 5-10 for quick testing and parameter tuning
            3. Use 50-100 for validation before full analysis
        ''',
        'impact': 'Only affects processing time, not results quality'
    }
}


def print_parameter_guide():
    """Print detailed parameter guide"""
    print("\n" + "="*80)
    print("PARAMETER GUIDE FOR NUCLEUS DETECTION")
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
    print("\nQUICK START GUIDE:")
    print("="*80)
    print("""
1. First run with default parameters to see results
2. If nuclei are missed:
   - Decrease min_intensity (e.g., from 10 to 5)
   - Decrease min_sphericity (e.g., from 0.5 to 0.4)
   - Adjust nucleus_diameter_min/max to match your data
3. If too many false positives:
   - Increase min_intensity
   - Increase min_sphericity
   - Decrease nucleus_diameter_max
4. If nuclei are detected but shapes are wrong:
   - Adjust gaussian_sigma (increase for smoother shapes)
    """)


def detect_nuclei(
    dapi_image_path,
    cell_labels_path,
    cell_info_path,
    output_dir,
    pixel_size_xy=0.44,
    pixel_size_z=0.50,
    nucleus_diameter_min=5.0,
    nucleus_diameter_max=10.0,
    min_intensity=10,
    max_intensity=40,
    min_sphericity=0.5,
    gaussian_sigma=2.0,
    process_n_cells=0
):
    """
    Detect spherical nuclei within segmented cells.
    
    See PARAM_DEFINITIONS for detailed parameter explanations.
    """
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print(f"Loading DAPI image from {dapi_image_path}...")
    dapi_image = io.imread(dapi_image_path)
    
    print(f"Loading cell labels from {cell_labels_path}...")
    cell_labels = io.imread(cell_labels_path)
    
    print(f"Loading cell info from {cell_info_path}...")
    with open(cell_info_path, 'r') as f:
        cell_info = json.load(f)
    
    print(f"DAPI image shape: {dapi_image.shape}")
    print(f"DAPI dtype: {dapi_image.dtype}, range: {dapi_image.min()}-{dapi_image.max()}")
    print(f"Total cells: {len(cell_info)}")
    
    # Calculate volume thresholds in voxels
    min_volume_voxels = int((4/3) * np.pi * ((nucleus_diameter_min / 2) / pixel_size_xy)**2 * ((nucleus_diameter_min / 2) / pixel_size_z))
    max_volume_voxels = int((4/3) * np.pi * ((nucleus_diameter_max / 2) / pixel_size_xy)**2 * ((nucleus_diameter_max / 2) / pixel_size_z))
    
    print(f"\nNucleus volume range: {min_volume_voxels}-{max_volume_voxels} voxels")
    print(f"Intensity range: {min_intensity}-{max_intensity}")
    print(f"Min sphericity: {min_sphericity}")
    
    # Smooth DAPI image
    print(f"\nSmoothing DAPI image (sigma={gaussian_sigma})...")
    dapi_smooth = filters.gaussian(dapi_image, sigma=(gaussian_sigma, gaussian_sigma, gaussian_sigma))
    
    # Process cells
    cells_to_process = cell_info[:process_n_cells] if process_n_cells > 0 else cell_info
    print(f"\nProcessing {len(cells_to_process)} cells...")
    
    nucleus_mask = np.zeros_like(cell_labels, dtype=bool)
    nucleus_info = []
    
    for i, cell in enumerate(cells_to_process):
        cell_id = cell['cell_id']
        
        if (i + 1) % 10 == 0:
            print(f"  Processing cell {i+1}/{len(cells_to_process)}...")
        
        # Extract cell mask
        cell_mask = (cell_labels == cell_id)
        
        # Extract DAPI signal in cell
        dapi_in_cell = dapi_smooth * cell_mask
        
        # Threshold
        nucleus_candidate = (dapi_in_cell >= min_intensity) & (dapi_in_cell <= max_intensity) & cell_mask
        
        if nucleus_candidate.sum() == 0:
            continue
        
        # Label connected components
        nucleus_labels = measure.label(nucleus_candidate)
        nucleus_props = measure.regionprops(nucleus_labels, intensity_image=dapi_smooth)
        
        # Filter by volume and sphericity
        for prop in nucleus_props:
            volume = prop.area
            
            if volume < min_volume_voxels or volume > max_volume_voxels:
                continue
            
            # Calculate sphericity
            surface_area = prop.area_filled * 6  # Approximation
            sphere_surface_area = np.pi ** (1/3) * (6 * volume) ** (2/3)
            sphericity = sphere_surface_area / surface_area if surface_area > 0 else 0
            
            if sphericity < min_sphericity:
                continue
            
            # This is a valid nucleus
            nucleus_mask[nucleus_labels == prop.label] = True
            
            # Calculate radius
            radius_voxels = (3 * volume / (4 * np.pi)) ** (1/3)
            radius_um = radius_voxels * ((pixel_size_xy**2 * pixel_size_z) ** (1/3))
            
            nucleus_info.append({
                'cell_id': cell_id,
                'nucleus_id': len(nucleus_info) + 1,
                'volume_voxels': int(volume),
                'centroid_z': float(prop.centroid[0]),
                'centroid_y': float(prop.centroid[1]),
                'centroid_x': float(prop.centroid[2]),
                'radius_um': round(radius_um, 2),
                'sphericity': round(sphericity, 2),
                'mean_intensity': round(prop.intensity_mean, 1),
                'max_intensity': int(prop.intensity_max)
            })
    
    print(f"\n✓ Detected {len(nucleus_info)} nuclei in {len(cells_to_process)} cells")
    print(f"  Detection rate: {len(nucleus_info)}/{len(cells_to_process)} = {100*len(nucleus_info)/len(cells_to_process):.1f}%")
    
    # Save results
    print("\nSaving results...")
    io.imsave(str(output_dir / 'nucleus_mask.tif'), nucleus_mask.astype(np.uint8) * 255, check_contrast=False)
    
    with open(output_dir / 'nucleus_info.json', 'w') as f:
        json.dump(nucleus_info, f, indent=2)
    
    # Generate report
    generate_report(dapi_image, cell_labels, nucleus_mask, nucleus_info, output_dir)
    
    print(f"\n✓ Nucleus detection complete!")
    print(f"  Output directory: {output_dir}")
    print(f"  Nuclei detected: {len(nucleus_info)}")
    
    return nucleus_mask, nucleus_info


def generate_report(dapi_image, cell_labels, nucleus_mask, nucleus_info, output_dir):
    """Generate visualization report"""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # DAPI MIP
    axes[0, 0].imshow(dapi_image.max(axis=0), cmap='gray')
    axes[0, 0].set_title('DAPI Signal (MIP)')
    axes[0, 0].axis('off')
    
    # Nucleus mask MIP
    axes[0, 1].imshow(nucleus_mask.max(axis=0), cmap='Blues')
    axes[0, 1].set_title(f'Nucleus Mask (MIP, N={len(nucleus_info)})')
    axes[0, 1].axis('off')
    
    # Overlay MIP
    axes[0, 2].imshow(dapi_image.max(axis=0), cmap='gray')
    axes[0, 2].imshow(nucleus_mask.max(axis=0), cmap='Reds', alpha=0.5)
    axes[0, 2].set_title('Overlay (MIP)')
    axes[0, 2].axis('off')
    
    # Z=20 slice
    z_slice = min(20, dapi_image.shape[0] - 1)
    axes[1, 0].imshow(dapi_image[z_slice], cmap='gray')
    axes[1, 0].set_title(f'DAPI Signal (Z={z_slice})')
    axes[1, 0].axis('off')
    
    axes[1, 1].imshow(nucleus_mask[z_slice], cmap='Blues')
    axes[1, 1].set_title(f'Nucleus Mask (Z={z_slice})')
    axes[1, 1].axis('off')
    
    axes[1, 2].imshow(dapi_image[z_slice], cmap='gray')
    axes[1, 2].imshow(nucleus_mask[z_slice], cmap='Reds', alpha=0.5)
    axes[1, 2].set_title(f'Overlay (Z={z_slice})')
    axes[1, 2].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'nucleus_detection_report.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Statistics
    if nucleus_info:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        
        radii = [n['radius_um'] for n in nucleus_info]
        sphericities = [n['sphericity'] for n in nucleus_info]
        
        axes[0].hist(radii, bins=20, edgecolor='black')
        axes[0].set_xlabel('Nucleus Radius (μm)')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Nucleus Size Distribution')
        axes[0].axvline(np.mean(radii), color='red', linestyle='--', label=f'Mean: {np.mean(radii):.2f} μm')
        axes[0].legend()
        
        axes[1].hist(sphericities, bins=20, edgecolor='black')
        axes[1].set_xlabel('Sphericity')
        axes[1].set_ylabel('Count')
        axes[1].set_title('Nucleus Sphericity Distribution')
        axes[1].axvline(np.mean(sphericities), color='red', linestyle='--', label=f'Mean: {np.mean(sphericities):.2f}')
        axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / 'nucleus_statistics.png', dpi=150, bbox_inches='tight')
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Detect spherical nuclei in segmented cells',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic usage with default parameters
  python 03_detect_nuclei.py dapi.tif cell_labels.tif cell_info.json output_dir/
  
  # Show parameter guide
  python 03_detect_nuclei.py --help-params
  
  # For weak DAPI signals
  python 03_detect_nuclei.py dapi.tif cell_labels.tif cell_info.json output_dir/ \\
    --min-intensity 5 --max-intensity 30 --min-sphericity 0.4
  
  # Quick test on first 10 cells
  python 03_detect_nuclei.py dapi.tif cell_labels.tif cell_info.json output_dir/ \\
    --process-n-cells 10
        '''
    )
    
    parser.add_argument('dapi_image', help='Path to DAPI channel image (TIFF)')
    parser.add_argument('cell_labels', help='Path to cell labels image (TIFF)')
    parser.add_argument('cell_info', help='Path to cell info JSON file')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('--help-params', action='store_true', help='Show detailed parameter guide')
    
    # Add all parameters
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
    
    detect_nuclei(
        args.dapi_image,
        args.cell_labels,
        args.cell_info,
        args.output_dir,
        pixel_size_xy=args.pixel_size_xy,
        pixel_size_z=args.pixel_size_z,
        nucleus_diameter_min=args.nucleus_diameter_min,
        nucleus_diameter_max=args.nucleus_diameter_max,
        min_intensity=args.min_intensity,
        max_intensity=args.max_intensity,
        min_sphericity=args.min_sphericity,
        gaussian_sigma=args.gaussian_sigma,
        process_n_cells=args.process_n_cells
    )


if __name__ == '__main__':
    main()
