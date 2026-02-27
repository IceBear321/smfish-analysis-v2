# Parameter Tuning Guide

This guide helps you optimize parameters for your specific microscopy data.

## Table of Contents

1. [Quick Start Workflow](#quick-start-workflow)
2. [Cell Segmentation Parameters](#cell-segmentation-parameters)
3. [Nucleus Detection Parameters](#nucleus-detection-parameters)
4. [RNA Detection Parameters](#rna-detection-parameters)
5. [Common Issues and Solutions](#common-issues-and-solutions)
6. [Parameter Dependencies](#parameter-dependencies)

---

## Quick Start Workflow

### Step 1: Determine Physical Dimensions

**Before anything else, you MUST know your pixel sizes:**

```bash
# Check CZI metadata
python -c "import czifile; czi = czifile.CziFile('your_file.czi'); print(czi.metadata)"
```

Look for `<Distance Id="X">` and `<Distance Id="Z">` tags.

**Typical values by objective:**
- 100x: 0.06-0.10 μm/pixel
- 63x: 0.10-0.16 μm/pixel
- 40x: 0.16-0.25 μm/pixel
- 20x: 0.30-0.50 μm/pixel (← your case: 0.44 μm/pixel)
- 10x: 0.60-1.00 μm/pixel

### Step 2: Test with Default Parameters

```bash
# Run full pipeline with defaults
python scripts/02_segment_cells.py channel_cfw.tif output/segmentation/
python scripts/03_detect_nuclei.py channel_dapi.tif output/segmentation/cell_labels.tif \\
    output/segmentation/cell_info.json output/nuclei/
```

### Step 3: Evaluate Results

Open the generated reports:
- `segmentation_report.png`: Check if cells are correctly segmented
- `nucleus_detection_report.png`: Check if nuclei are detected

### Step 4: Adjust Parameters

See sections below for specific parameter tuning.

---

## Cell Segmentation Parameters

### Parameter Priority (tune in this order)

1. **`pixel_size_xy` and `pixel_size_z`** (CRITICAL)
   - Must be correct, or all size-based filtering will fail
   - Check microscope metadata

2. **`canny_sigma`** (affects edge quality)
   - Default: 2.0
   - **Problem**: Edges are broken/incomplete
     - **Solution**: Increase to 3.0-5.0 (smoother edges)
   - **Problem**: Edges are too thick, cells merge
     - **Solution**: Decrease to 1.0-1.5 (sharper edges)

3. **`canny_low` and `canny_high`** (affects edge sensitivity)
   - Default: 10, 30
   - **Problem**: Too many edges (background noise)
     - **Solution**: Increase both (e.g., 15, 40)
   - **Problem**: Missing cell walls
     - **Solution**: Decrease both (e.g., 5, 20)
   - **Rule**: `canny_high` should be 2-3x of `canny_low`

4. **`min_distance`** (affects cell separation)
   - Default: 5 pixels
   - **Problem**: One cell split into multiple
     - **Solution**: Increase (e.g., 10-15)
   - **Problem**: Multiple cells merged into one
     - **Solution**: Decrease (e.g., 3-5)
   - **Formula**: `min_distance ≈ min_cell_diameter / (2 * pixel_size_xy)`

5. **`min_cell_size` and `max_cell_size`** (filters debris and merged cells)
   - Calculate expected cell volume:
     ```
     For 10 μm diameter cell with 0.44 μm/pixel XY and 0.5 μm/pixel Z:
     volume = (10³) / (0.44² * 0.5) ≈ 10,300 voxels
     ```
   - Set `min_cell_size` to 20-30% of expected (e.g., 2000-3000)
   - Set `max_cell_size` to 300-500% of expected (e.g., 30,000-50,000)

### Example Configurations

**For noisy images:**
```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-sigma 3.0 \\
    --canny-low 15 \\
    --canny-high 40
```

**For small cells (5-10 μm):**
```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-distance 3 \\
    --min-cell-size 500 \\
    --max-cell-size 10000
```

**For large cells (20-30 μm):**
```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-distance 15 \\
    --min-cell-size 10000 \\
    --max-cell-size 100000
```

---

## Nucleus Detection Parameters

### Parameter Priority (tune in this order)

1. **`nucleus_diameter_min` and `nucleus_diameter_max`** (CRITICAL)
   - Measure nucleus diameter in ImageJ/Fiji first!
   - For plant cells: typically 5-10 μm
   - Set min to 80% of smallest, max to 120% of largest

2. **`min_intensity` and `max_intensity`** (affects detection sensitivity)
   - Default: 10, 40 (for weak DAPI)
   - **How to determine**:
     ```
     1. Open image in ImageJ/Fiji
     2. Draw ROI around a nucleus
     3. Analyze → Measure (Ctrl+M)
     4. Note "Mean" intensity
     5. Set min_intensity = 0.5 * mean
     6. Set max_intensity = 1.5 * mean
     ```
   - **Problem**: No nuclei detected
     - **Solution**: Decrease `min_intensity` (e.g., 5)
   - **Problem**: Too many false positives
     - **Solution**: Increase `min_intensity` (e.g., 20)

3. **`min_sphericity`** (filters non-spherical objects)
   - Default: 0.5
   - **Problem**: Missing elongated nuclei
     - **Solution**: Decrease to 0.4
   - **Problem**: Detecting non-nuclear structures
     - **Solution**: Increase to 0.6-0.7

4. **`gaussian_sigma`** (smoothing)
   - Default: 2.0
   - **Problem**: Nuclei have rough edges
     - **Solution**: Increase to 3.0-4.0
   - **Problem**: Small nuclei disappear
     - **Solution**: Decrease to 1.0-1.5

### Example Configurations

**For weak DAPI signals:**
```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-intensity 5 \\
    --max-intensity 30 \\
    --min-sphericity 0.4 \\
    --gaussian-sigma 2.5
```

**For strong DAPI signals:**
```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-intensity 30 \\
    --max-intensity 100 \\
    --min-sphericity 0.6
```

**For small nuclei (3-5 μm):**
```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --nucleus-diameter-min 3.0 \\
    --nucleus-diameter-max 6.0 \\
    --gaussian-sigma 1.5
```

---

## RNA Detection Parameters

### Key Parameters

1. **`min_sigma` and `max_sigma`** (blob size range)
   - Default: 1.0, 3.0 pixels
   - Determines the size of RNA spots to detect
   - **Formula**: `sigma ≈ RNA_spot_diameter / (2 * pixel_size_xy)`
   - For 0.5-1.5 μm RNA spots with 0.44 μm/pixel:
     - `min_sigma ≈ 0.5 / (2 * 0.44) ≈ 0.6`
     - `max_sigma ≈ 1.5 / (2 * 0.44) ≈ 1.7`

2. **`threshold`** (detection sensitivity)
   - Default: 0.01 (relative to image max)
   - **Problem**: Too few RNA detected
     - **Solution**: Decrease to 0.005
   - **Problem**: Too many false positives
     - **Solution**: Increase to 0.02-0.05

3. **`overlap`** (how much blobs can overlap)
   - Default: 0.5 (50%)
   - Increase to 0.7 if RNA spots are very close
   - Decrease to 0.3 if detecting merged spots

---

## Common Issues and Solutions

### Issue 1: No cells detected

**Possible causes:**
1. Wrong pixel size → all cells filtered out by size
2. `min_distance` too large → no seed points found
3. `canny_low/high` too high → no edges detected

**Solutions:**
```bash
# Check image first
python -c "from skimage import io; img = io.imread('input.tif'); print(f'Shape: {img.shape}, Range: {img.min()}-{img.max()}')"

# Try with relaxed parameters
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-low 5 \\
    --canny-high 20 \\
    --min-distance 3 \\
    --min-cell-size 500
```

### Issue 2: Cells are merged

**Cause**: `min_distance` too large or edges incomplete

**Solution:**
```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-distance 3 \\
    --canny-sigma 2.5  # Smoother edges to connect gaps
```

### Issue 3: Cells are over-segmented

**Cause**: `min_distance` too small

**Solution:**
```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-distance 10
```

### Issue 4: No nuclei detected

**Possible causes:**
1. DAPI signal too weak
2. `min_intensity` too high
3. Nucleus size range wrong

**Solutions:**
```bash
# Check DAPI intensity
python -c "from skimage import io; import numpy as np; img = io.imread('dapi.tif'); print(f'Mean: {img.mean():.1f}, Std: {img.std():.1f}, Max: {img.max()}')"

# Try with lower threshold
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-intensity 5 \\
    --max-intensity 50 \\
    --min-sphericity 0.3
```

### Issue 5: Too many false nucleus detections

**Cause**: `min_intensity` too low or `min_sphericity` too low

**Solution:**
```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-intensity 15 \\
    --min-sphericity 0.6
```

---

## Parameter Dependencies

### Cell Segmentation → Nucleus Detection

The following parameters MUST match:
- `pixel_size_xy`
- `pixel_size_z`

### Nucleus Detection → RNA Detection

The nucleus mask is used to classify RNA as nuclear vs cytoplasmic.

### All Steps → 3D Visualization

Physical dimensions affect the aspect ratio of 3D visualizations.

---

## Iterative Tuning Workflow

1. **Start with defaults**
   ```bash
   ./run_pipeline.sh input.czi output/
   ```

2. **Check segmentation report**
   - If cells look good → proceed to step 4
   - If not → adjust cell segmentation parameters (step 3)

3. **Tune cell segmentation**
   ```bash
   python scripts/02_segment_cells.py --help-params  # Read guide
   # Adjust parameters based on issues observed
   python scripts/02_segment_cells.py input.tif output/ --canny-sigma 3.0 ...
   ```

4. **Check nucleus detection report**
   - If nuclei look good → proceed to step 6
   - If not → adjust nucleus detection parameters (step 5)

5. **Tune nucleus detection**
   ```bash
   python scripts/03_detect_nuclei.py --help-params  # Read guide
   # Adjust parameters
   python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ --min-intensity 5 ...
   ```

6. **Run RNA detection and 3D visualization**
   ```bash
   python scripts/04_detect_rna.py ...
   python scripts/05_visualize_3d.py ...
   ```

7. **Iterate** until results are satisfactory

---

## Tips for Success

1. **Always visualize intermediate results**
   - Don't skip the report images
   - Use ImageJ/Fiji to inspect TIFF files

2. **Start with a small test region**
   - Use `--roi-y-min` to focus on a specific area
   - Use `--process-n-cells 10` for quick testing

3. **Document your parameters**
   - Save successful parameter sets in `configs/` directory
   - Use descriptive names (e.g., `weak_dapi_config.txt`)

4. **Measure ground truth**
   - Manually measure cell/nucleus sizes in ImageJ
   - Use these measurements to set size thresholds

5. **Check image quality first**
   - Noisy images → increase smoothing (`canny_sigma`, `gaussian_sigma`)
   - Blurry images → decrease smoothing

6. **One parameter at a time**
   - Change one parameter, check result, repeat
   - Don't change multiple parameters simultaneously

---

## Getting Help

If you're stuck:

1. Run with `--help-params` to see detailed parameter explanations
2. Check `docs/TROUBLESHOOTING.md` for common issues
3. See `examples/` for configuration examples
4. Open an issue on GitHub with:
   - Your microscope setup (objective, pixel size)
   - Example image (small region)
   - Parameters you tried
   - Error messages or unexpected results
