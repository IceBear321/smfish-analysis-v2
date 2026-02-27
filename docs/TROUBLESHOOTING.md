# Troubleshooting Guide

This guide helps you diagnose and fix common issues in the smFISH analysis pipeline.

## Table of Contents

1. [Installation Issues](#installation-issues)
2. [Cell Segmentation Issues](#cell-segmentation-issues)
3. [Nucleus Detection Issues](#nucleus-detection-issues)
4. [RNA Detection Issues](#rna-detection-issues)
5. [Visualization Issues](#visualization-issues)
6. [Performance Issues](#performance-issues)
7. [File Format Issues](#file-format-issues)

---

## Installation Issues

### Issue: `ModuleNotFoundError: No module named 'aicspylibczi'`

**Cause**: Required packages not installed

**Solution**:
```bash
pip install -r requirements.txt
```

### Issue: `ImportError: libGL.so.1: cannot open shared object file`

**Cause**: Missing system libraries for image processing

**Solution** (Ubuntu/Debian):
```bash
sudo apt-get update
sudo apt-get install -y libgl1-mesa-glx libglib2.0-0
```

**Solution** (CentOS/RHEL):
```bash
sudo yum install -y mesa-libGL glib2
```

### Issue: `pip install` fails with permission error

**Cause**: No write permission to system Python

**Solution**: Use virtual environment
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

## Cell Segmentation Issues

### Issue: No cells detected

**Symptoms**:
- `cell_labels.tif` is all zeros
- `cell_info.json` is empty
- Report shows "0 cells detected"

**Possible Causes & Solutions**:

#### 1. Wrong pixel size

**Check**:
```bash
# Your pixel size in parameters
grep "pixel_size" your_command.sh

# Actual pixel size in CZI metadata
python -c "import czifile; czi = czifile.CziFile('input.czi'); print(czi.metadata)" | grep Distance
```

**Solution**: Correct the `--pixel-size-xy` and `--pixel-size-z` parameters

#### 2. All cells filtered out by size

**Check**: Look at the console output for "Filtered by size" messages

**Solution**: Decrease `--min-cell-size` or increase `--max-cell-size`
```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-cell-size 500 \\
    --max-cell-size 100000
```

#### 3. No edges detected

**Check**: Open `output/segmentation/segmentation_report.png` and look at the "Edges" panel

**Solution**: Decrease Canny thresholds
```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-low 5 \\
    --canny-high 15
```

#### 4. No seed points for watershed

**Check**: Look for "Found X seed points" in console output

**Solution**: Decrease `--min-distance`
```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-distance 3
```

---

### Issue: Cells are merged (under-segmentation)

**Symptoms**:
- Multiple cells labeled as one
- Cell count is too low
- Cell sizes are too large

**Solutions**:

#### 1. Decrease min_distance

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-distance 3  # Default is 5
```

**Rationale**: More seed points → better separation

#### 2. Improve edge detection

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-sigma 2.5 \\  # Smoother edges
    --canny-low 8 \\
    --canny-high 25
```

**Rationale**: Better edges → better watershed boundaries

#### 3. Filter out merged cells

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --max-cell-size 30000  # Remove abnormally large cells
```

---

### Issue: Cells are over-segmented

**Symptoms**:
- One cell split into multiple pieces
- Cell count is too high
- Cell sizes are too small

**Solutions**:

#### 1. Increase min_distance

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-distance 10  # Default is 5
```

**Rationale**: Fewer seed points → less splitting

#### 2. Increase smoothing

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-sigma 3.0  # More smoothing
```

**Rationale**: Smoother edges → less fragmentation

#### 3. Filter out small fragments

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-cell-size 3000  # Remove small pieces
```

---

### Issue: Edges are incomplete

**Symptoms**:
- Cell walls have gaps
- Cells leak into each other

**Solution**: Increase smoothing to connect gaps
```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-sigma 3.0  # Higher smoothing
```

---

### Issue: Too much background detected as cells

**Symptoms**:
- Debris labeled as cells
- Non-cellular regions segmented

**Solutions**:

#### 1. Increase edge thresholds

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-low 15 \\
    --canny-high 40
```

#### 2. Use ROI to focus on cell region

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --roi-y-min 256  # Skip upper half
```

#### 3. Increase minimum cell size

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --min-cell-size 3000
```

---

## Nucleus Detection Issues

### Issue: No nuclei detected

**Symptoms**:
- `nucleus_mask.tif` is all zeros
- `nucleus_info.json` shows no nuclei
- Report shows "0 nuclei detected"

**Possible Causes & Solutions**:

#### 1. DAPI signal too weak

**Check**: Measure DAPI intensity in ImageJ
```
1. Open DAPI channel in ImageJ
2. Draw ROI around a nucleus
3. Analyze → Measure
4. Note "Mean" intensity
```

**Solution**: Lower `--min-intensity`
```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-intensity 5  # Default is 10
```

#### 2. Nucleus size range wrong

**Check**: Measure nucleus diameter in ImageJ

**Solution**: Adjust diameter range
```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --nucleus-diameter-min 3.0 \\  # Smaller
    --nucleus-diameter-max 12.0    # Larger
```

#### 3. Nuclei filtered by sphericity

**Check**: Look for "Filtered by sphericity" in console output

**Solution**: Relax sphericity requirement
```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-sphericity 0.3  # Default is 0.5
```

#### 4. No cells detected in previous step

**Check**: Verify `cell_labels.tif` is not empty

**Solution**: Fix cell segmentation first (see above)

---

### Issue: Too many false positive nuclei

**Symptoms**:
- Non-nuclear structures detected
- More nuclei than cells
- Nuclei in background regions

**Solutions**:

#### 1. Increase intensity threshold

```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-intensity 20  # Higher threshold
```

#### 2. Increase sphericity requirement

```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-sphericity 0.6  # Stricter shape filtering
```

#### 3. Narrow size range

```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --nucleus-diameter-min 6.0 \\
    --nucleus-diameter-max 9.0  # Tighter range
```

---

### Issue: Nuclei are fragmented

**Symptoms**:
- One nucleus split into multiple pieces
- Nucleus count is too high

**Solution**: Increase smoothing
```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --gaussian-sigma 3.0  # More smoothing
```

---

### Issue: Nuclei have irregular shapes

**Symptoms**:
- Nuclei are not spherical
- Nuclei have rough edges

**Solutions**:

#### 1. Increase smoothing

```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --gaussian-sigma 3.5
```

#### 2. Relax sphericity (if nuclei are naturally elongated)

```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-sphericity 0.4
```

---

### Issue: Nucleus extends outside cell boundary

**Symptoms**:
- Nucleus mask overlaps with adjacent cells
- Nucleus is larger than cell

**Cause**: This is expected behavior in the current algorithm (Version 1)

**Note**: The algorithm uses smoothed masks for RNA filtering but does not constrain nucleus within cell boundary. This is by design to avoid losing RNA molecules near boundaries.

**If you need strict containment**: Modify `03_detect_nuclei.py` to add boundary constraint (see code comments)

---

## RNA Detection Issues

### Issue: No RNA detected

**Symptoms**:
- `rna_molecules.json` is empty
- RNA count is 0 for all cells

**Possible Causes & Solutions**:

#### 1. Threshold too high

**Solution**: Decrease threshold
```bash
python scripts/04_detect_rna.py cy3.tif labels.tif nucleus.tif output/ \\
    --threshold 0.005  # Default is 0.01
```

#### 2. Sigma range wrong

**Check**: Measure RNA spot size in ImageJ

**Solution**: Adjust sigma range
```bash
python scripts/04_detect_rna.py cy3.tif labels.tif nucleus.tif output/ \\
    --min-sigma 0.5 \\
    --max-sigma 4.0
```

#### 3. No cells detected

**Check**: Verify previous steps completed successfully

---

### Issue: Too many false positive RNA spots

**Symptoms**:
- RNA count is abnormally high (>1000 per cell)
- Background noise detected as RNA

**Solutions**:

#### 1. Increase threshold

```bash
python scripts/04_detect_rna.py cy3.tif labels.tif nucleus.tif output/ \\
    --threshold 0.02  # Higher threshold
```

#### 2. Narrow sigma range

```bash
python scripts/04_detect_rna.py cy3.tif labels.tif nucleus.tif output/ \\
    --min-sigma 1.5 \\
    --max-sigma 2.5  # Tighter range
```

#### 3. Decrease overlap

```bash
python scripts/04_detect_rna.py cy3.tif labels.tif nucleus.tif output/ \\
    --overlap 0.3  # Less overlap allowed
```

---

### Issue: RNA spots are merged

**Symptoms**:
- Large blobs detected instead of individual spots
- RNA count is too low

**Solution**: Decrease overlap threshold
```bash
python scripts/04_detect_rna.py cy3.tif labels.tif nucleus.tif output/ \\
    --overlap 0.3  # Default is 0.5
```

---

## Visualization Issues

### Issue: `ValueError: No cells to visualize`

**Cause**: Specified cell IDs not found

**Solution**: Check available cell IDs
```bash
# List all detected cells
python -c "
import json
with open('output/segmentation/cell_info.json') as f:
    cells = json.load(f)
print('Available cell IDs:', list(cells.keys()))
"

# Then visualize existing cells
python scripts/05_visualize_3d.py ... --cell-ids 1 2 3
```

---

### Issue: 3D visualization is distorted

**Symptoms**:
- Cells appear stretched or squashed
- Aspect ratio is wrong

**Cause**: Incorrect pixel sizes

**Solution**: Verify and correct pixel sizes
```bash
python scripts/05_visualize_3d.py ... \\
    --pixel-size-xy 0.44 \\
    --pixel-size-z 0.50
```

---

### Issue: Visualization is too slow

**Symptoms**:
- Browser freezes
- HTML file is very large

**Solutions**:

#### 1. Visualize fewer cells

```bash
python scripts/05_visualize_3d.py ... \\
    --cell-ids 1 2  # Only 2 cells instead of all
```

#### 2. Reduce smoothing (fewer vertices)

```bash
python scripts/05_visualize_3d.py ... \\
    --cell-smooth 1.0 \\  # Less smoothing
    --nucleus-smooth 0.5
```

#### 3. Make RNA spots smaller

```bash
python scripts/05_visualize_3d.py ... \\
    --rna-size 1.0  # Smaller spots
```

---

### Issue: RNA spots not visible in 3D

**Cause**: RNA spots too small or wrong color

**Solution**: Increase RNA size
```bash
python scripts/05_visualize_3d.py ... \\
    --rna-size 3.0  # Larger spots
```

---

### Issue: Surfaces are too rough

**Cause**: Not enough smoothing

**Solution**: Increase smoothing
```bash
python scripts/05_visualize_3d.py ... \\
    --cell-smooth 2.0 \\
    --nucleus-smooth 1.5
```

---

### Issue: Surfaces are too smooth (loss of detail)

**Cause**: Too much smoothing

**Solution**: Decrease smoothing
```bash
python scripts/05_visualize_3d.py ... \\
    --cell-smooth 0.5 \\
    --nucleus-smooth 0.5
```

---

## Performance Issues

### Issue: Processing is very slow

**Solutions**:

#### 1. Process subset of data

```bash
# Focus on specific region
python scripts/02_segment_cells.py input.tif output/ \\
    --roi-y-min 256 \\
    --roi-y-max 512

# Process only a few cells
python scripts/04_detect_rna.py ... \\
    --process-n-cells 10
```

#### 2. Reduce image size

```bash
# Downsample image (if resolution allows)
python -c "
from skimage import io, transform
img = io.imread('input.tif')
img_small = transform.downscale_local_mean(img, (1, 2, 2))  # Downsample XY by 2x
io.imsave('input_small.tif', img_small.astype(img.dtype))
"
```

#### 3. Use fewer Z-slices

```bash
# Extract every other slice
python -c "
from skimage import io
img = io.imread('input.tif')
img_subset = img[::2, :, :]  # Every 2nd slice
io.imsave('input_subset.tif', img_subset)
"
```

---

### Issue: Out of memory error

**Solutions**:

#### 1. Process in chunks

```bash
# Process cells in batches
python scripts/04_detect_rna.py ... --process-n-cells 20
# Then process next 20, etc.
```

#### 2. Reduce image bit depth

```bash
# Convert 16-bit to 8-bit (if dynamic range allows)
python -c "
from skimage import io, exposure
img = io.imread('input.tif')
img_8bit = exposure.rescale_intensity(img, out_range='uint8').astype('uint8')
io.imsave('input_8bit.tif', img_8bit)
"
```

#### 3. Close other applications

Free up RAM by closing unnecessary programs

---

## File Format Issues

### Issue: `ValueError: Cannot read CZI file`

**Cause**: Corrupted or unsupported CZI format

**Solutions**:

#### 1. Try with different CZI reader

```bash
# Use bioformats instead
pip install python-bioformats
```

#### 2. Convert to TIFF in ImageJ/Fiji

```
1. Open CZI in Fiji
2. File → Save As → TIFF
3. Use TIFF files as input
```

---

### Issue: `TypeError: Image must be 2D or 3D`

**Cause**: Wrong image dimensions (e.g., 4D or 5D)

**Solution**: Extract specific channel/timepoint
```bash
python -c "
from skimage import io
img = io.imread('input.tif')
print(f'Shape: {img.shape}')
# If 4D (T, Z, Y, X), extract first timepoint
if len(img.shape) == 4:
    img_3d = img[0, :, :, :]
    io.imsave('input_3d.tif', img_3d)
"
```

---

### Issue: `OSError: Unable to open file`

**Cause**: File path is wrong or file doesn't exist

**Solution**: Check file path
```bash
# Verify file exists
ls -lh /path/to/your/file.tif

# Use absolute path
python scripts/02_segment_cells.py /absolute/path/to/input.tif output/
```

---

## Getting More Help

### Enable debug output

```bash
# Add --verbose flag (if available)
python scripts/02_segment_cells.py input.tif output/ --verbose
```

### Check intermediate files

```bash
# Inspect intermediate results
ls -lh output/segmentation/
# Open TIFF files in ImageJ/Fiji
```

### Report an issue

If you're still stuck, open an issue on GitHub with:

1. **System information**:
   ```bash
   python --version
   pip list | grep -E "scikit-image|scipy|numpy"
   ```

2. **Command you ran**:
   ```bash
   python scripts/02_segment_cells.py input.tif output/ --canny-sigma 2.0 ...
   ```

3. **Error message** (full traceback)

4. **Example data** (small region if possible):
   ```bash
   # Extract small region for sharing
   python -c "
   from skimage import io
   img = io.imread('input.tif')
   small = img[:10, 100:200, 100:200]  # Small region
   io.imsave('input_small.tif', small)
   "
   ```

5. **What you expected vs what you got**

6. **Microscope setup** (objective, pixel size, sample type)

---

## Quick Diagnostic Checklist

Before asking for help, check:

- [ ] Pixel sizes are correct
- [ ] Input files exist and are readable
- [ ] Image dimensions are correct (3D: Z, Y, X)
- [ ] Image intensity range is reasonable (not all zeros)
- [ ] Previous pipeline steps completed successfully
- [ ] Output directories have write permission
- [ ] Sufficient disk space available
- [ ] Sufficient RAM available
- [ ] Python packages are up to date

```bash
# Quick diagnostic script
python -c "
import sys
from skimage import io
import numpy as np

print('Python version:', sys.version)
print('NumPy version:', np.__version__)

img = io.imread('input.tif')
print(f'Image shape: {img.shape}')
print(f'Image dtype: {img.dtype}')
print(f'Image range: {img.min()}-{img.max()}')
print(f'Image mean: {img.mean():.1f}')
print(f'Image size: {img.nbytes / 1e6:.1f} MB')
"
```
