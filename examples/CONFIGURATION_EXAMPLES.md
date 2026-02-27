# Configuration Examples for Different Scenarios

This document provides ready-to-use parameter configurations for common microscopy scenarios.

## Table of Contents

1. [Plant Cells with Weak DAPI](#plant-cells-with-weak-dapi) ← **Your case (stack_oil_2.czi)**
2. [Plant Cells with Strong DAPI](#plant-cells-with-strong-dapi)
3. [Animal Cells (Mammalian)](#animal-cells-mammalian)
4. [Yeast Cells](#yeast-cells)
5. [High-Resolution Confocal](#high-resolution-confocal)
6. [Low-Resolution Widefield](#low-resolution-widefield)
7. [Noisy Images](#noisy-images)
8. [Dense Cell Populations](#dense-cell-populations)

---

## Plant Cells with Weak DAPI

**Scenario**: Plant root/leaf cells with CFW cell wall staining and weak DAPI nuclear staining

**Microscope setup**:
- Objective: 20x
- Pixel size: 0.44 μm/pixel (XY), 0.50 μm/pixel (Z)
- Cell size: 10-20 μm
- Nucleus size: 5-10 μm
- DAPI intensity: 10-40 (8-bit)

### Cell Segmentation

```bash
python scripts/02_segment_cells.py channel_cfw.tif output/segmentation/ \\
    --pixel-size-xy 0.44 \\
    --pixel-size-z 0.50 \\
    --canny-sigma 2.0 \\
    --canny-low 10 \\
    --canny-high 30 \\
    --min-distance 5 \\
    --min-cell-size 2000 \\
    --max-cell-size 50000 \\
    --roi-y-min 256  # Focus on lower half where cells are
```

**Rationale**:
- `canny_sigma=2.0`: Moderate smoothing for incomplete cell walls
- `min_distance=5`: Small value for densely packed cells
- `roi_y_min=256`: Skip upper half with no cells/debris

### Nucleus Detection

```bash
python scripts/03_detect_nuclei.py channel_dapi.tif \\
    output/segmentation/cell_labels.tif \\
    output/segmentation/cell_info.json \\
    output/nuclei/ \\
    --pixel-size-xy 0.44 \\
    --pixel-size-z 0.50 \\
    --nucleus-diameter-min 5.0 \\
    --nucleus-diameter-max 10.0 \\
    --min-intensity 10 \\
    --max-intensity 40 \\
    --min-sphericity 0.5 \\
    --gaussian-sigma 2.0
```

**Rationale**:
- `min_intensity=10`: Low threshold for weak DAPI
- `max_intensity=40`: Upper bound to avoid saturated pixels
- `min_sphericity=0.5`: Moderate requirement (plant nuclei can be slightly elongated)

### RNA Detection

```bash
python scripts/04_detect_rna.py channel_cy3.tif \\
    output/segmentation/cell_labels.tif \\
    output/nuclei/nucleus_mask.tif \\
    output/rna/ \\
    --pixel-size-xy 0.44 \\
    --pixel-size-z 0.50 \\
    --min-sigma 1.0 \\
    --max-sigma 3.0 \\
    --threshold 0.01 \\
    --overlap 0.5
```

### 3D Visualization

```bash
python scripts/05_visualize_3d.py \\
    output/segmentation/cell_labels.tif \\
    output/nuclei/nucleus_info.json \\
    output/rna/rna_molecules.json \\
    output/visualization/ \\
    --cell-ids 23 14 \\  # Visualize specific cells
    --cell-smooth 1.5 \\
    --nucleus-smooth 1.0 \\
    --rna-size 2.0
```

---

## Plant Cells with Strong DAPI

**Scenario**: Plant cells with good DAPI staining (e.g., meristematic cells)

**Differences from weak DAPI**:
- Higher DAPI intensity: 30-150
- More spherical nuclei

### Nucleus Detection (adjusted)

```bash
python scripts/03_detect_nuclei.py channel_dapi.tif \\
    output/segmentation/cell_labels.tif \\
    output/segmentation/cell_info.json \\
    output/nuclei/ \\
    --pixel-size-xy 0.44 \\
    --pixel-size-z 0.50 \\
    --nucleus-diameter-min 5.0 \\
    --nucleus-diameter-max 10.0 \\
    --min-intensity 30 \\      # Higher threshold
    --max-intensity 150 \\     # Higher upper bound
    --min-sphericity 0.6 \\    # Stricter sphericity
    --gaussian-sigma 2.0
```

---

## Animal Cells (Mammalian)

**Scenario**: HeLa, HEK293, or other mammalian cell lines

**Microscope setup**:
- Objective: 63x
- Pixel size: 0.10 μm/pixel (XY), 0.30 μm/pixel (Z)
- Cell size: 15-25 μm
- Nucleus size: 8-15 μm
- DAPI intensity: 50-200 (strong)

### Cell Segmentation

```bash
python scripts/02_segment_cells.py channel_membrane.tif output/segmentation/ \\
    --pixel-size-xy 0.10 \\
    --pixel-size-z 0.30 \\
    --canny-sigma 1.5 \\
    --canny-low 15 \\
    --canny-high 40 \\
    --min-distance 80 \\       # Larger in pixels (8 μm / 0.1 μm/pixel)
    --min-cell-size 50000 \\   # (15 μm)³ / (0.1² * 0.3) ≈ 112,500 voxels
    --max-cell-size 500000
```

### Nucleus Detection

```bash
python scripts/03_detect_nuclei.py channel_dapi.tif \\
    output/segmentation/cell_labels.tif \\
    output/segmentation/cell_info.json \\
    output/nuclei/ \\
    --pixel-size-xy 0.10 \\
    --pixel-size-z 0.30 \\
    --nucleus-diameter-min 8.0 \\
    --nucleus-diameter-max 15.0 \\
    --min-intensity 50 \\      # Strong DAPI
    --max-intensity 200 \\
    --min-sphericity 0.7 \\    # Mammalian nuclei are more spherical
    --gaussian-sigma 3.0
```

---

## Yeast Cells

**Scenario**: S. cerevisiae or S. pombe

**Microscope setup**:
- Objective: 100x
- Pixel size: 0.065 μm/pixel (XY), 0.20 μm/pixel (Z)
- Cell size: 3-5 μm
- Nucleus size: 1-2 μm

### Cell Segmentation

```bash
python scripts/02_segment_cells.py channel_brightfield.tif output/segmentation/ \\
    --pixel-size-xy 0.065 \\
    --pixel-size-z 0.20 \\
    --canny-sigma 1.0 \\       # Small sigma for small cells
    --canny-low 10 \\
    --canny-high 25 \\
    --min-distance 30 \\       # (3 μm / 0.065 μm/pixel) ≈ 46 pixels
    --min-cell-size 500 \\     # Small cells
    --max-cell-size 5000
```

### Nucleus Detection

```bash
python scripts/03_detect_nuclei.py channel_dapi.tif \\
    output/segmentation/cell_labels.tif \\
    output/segmentation/cell_info.json \\
    output/nuclei/ \\
    --pixel-size-xy 0.065 \\
    --pixel-size-z 0.20 \\
    --nucleus-diameter-min 1.0 \\
    --nucleus-diameter-max 2.5 \\
    --min-intensity 20 \\
    --max-intensity 100 \\
    --min-sphericity 0.6 \\
    --gaussian-sigma 1.5       # Small sigma for small nuclei
```

---

## High-Resolution Confocal

**Scenario**: High-quality confocal images with good SNR

**Characteristics**:
- Sharp edges
- Low noise
- Good contrast

### Recommended Parameters

```bash
# Cell segmentation
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-sigma 1.0 \\       # Low smoothing (sharp edges)
    --canny-low 15 \\          # Higher threshold (good contrast)
    --canny-high 40 \\
    --min-distance 5

# Nucleus detection
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --gaussian-sigma 1.5 \\    # Low smoothing (low noise)
    --min-sphericity 0.6       # Stricter (good shape preservation)
```

---

## Low-Resolution Widefield

**Scenario**: Widefield microscopy with lower SNR

**Characteristics**:
- Blurry edges
- Higher noise
- Lower contrast

### Recommended Parameters

```bash
# Cell segmentation
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-sigma 3.0 \\       # High smoothing (reduce noise)
    --canny-low 5 \\           # Lower threshold (low contrast)
    --canny-high 15 \\
    --min-distance 10

# Nucleus detection
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --gaussian-sigma 3.0 \\    # High smoothing (reduce noise)
    --min-sphericity 0.4 \\    # Relaxed (blurry shapes)
    --min-intensity 5          # Lower threshold
```

---

## Noisy Images

**Scenario**: Images with high background noise or autofluorescence

### Strategies

1. **Increase smoothing**:
   ```bash
   --canny-sigma 4.0
   --gaussian-sigma 3.5
   ```

2. **Increase thresholds**:
   ```bash
   --canny-low 20
   --canny-high 50
   --min-intensity 20
   ```

3. **Stricter shape filtering**:
   ```bash
   --min-sphericity 0.6
   ```

### Full Example

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-sigma 4.0 \\
    --canny-low 20 \\
    --canny-high 50

python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --gaussian-sigma 3.5 \\
    --min-intensity 20 \\
    --min-sphericity 0.6
```

---

## Dense Cell Populations

**Scenario**: Cells are tightly packed with minimal gaps

### Challenges

- Cells touching each other
- Difficult to separate individual cells
- Watershed may over-segment

### Recommended Parameters

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-sigma 2.5 \\       # Moderate smoothing to connect gaps
    --min-distance 8 \\        # Larger to prevent over-segmentation
    --min-cell-size 3000       # Filter out small fragments
```

**Tips**:
- Increase `min_distance` until cells are no longer over-segmented
- Use `canny_sigma` to smooth edges and connect small gaps
- May need to manually correct some merged cells

---

## Creating Your Own Configuration

### Step 1: Measure Your Data

```bash
# Check image properties
python -c "
from skimage import io
import numpy as np
img = io.imread('your_image.tif')
print(f'Shape: {img.shape}')
print(f'Dtype: {img.dtype}')
print(f'Range: {img.min()}-{img.max()}')
print(f'Mean: {img.mean():.1f}')
print(f'Std: {img.std():.1f}')
"
```

### Step 2: Measure Cell/Nucleus Sizes in ImageJ

1. Open image in ImageJ/Fiji
2. Draw line across cell/nucleus
3. Analyze → Measure (Ctrl+M)
4. Note the diameter in μm

### Step 3: Calculate Parameters

```python
# Cell size in voxels
pixel_size_xy = 0.44  # μm/pixel
pixel_size_z = 0.50   # μm/pixel
cell_diameter_um = 15  # measured

cell_diameter_pixels = cell_diameter_um / pixel_size_xy
cell_volume_voxels = (cell_diameter_um ** 3) / (pixel_size_xy ** 2 * pixel_size_z)

print(f"Cell diameter: {cell_diameter_pixels:.1f} pixels")
print(f"Cell volume: {cell_volume_voxels:.0f} voxels")

# Set parameters
min_distance = int(cell_diameter_pixels / 2)
min_cell_size = int(cell_volume_voxels * 0.2)  # 20% of expected
max_cell_size = int(cell_volume_voxels * 5.0)  # 500% of expected
```

### Step 4: Save Configuration

Create a file `configs/my_config.txt`:

```bash
# My microscopy configuration
# Date: 2024-01-01
# Sample: Plant root cells
# Objective: 20x

# Physical dimensions
--pixel-size-xy 0.44
--pixel-size-z 0.50

# Cell segmentation
--canny-sigma 2.0
--canny-low 10
--canny-high 30
--min-distance 5
--min-cell-size 2000
--max-cell-size 50000

# Nucleus detection
--nucleus-diameter-min 5.0
--nucleus-diameter-max 10.0
--min-intensity 10
--max-intensity 40
--min-sphericity 0.5
--gaussian-sigma 2.0
```

### Step 5: Use Configuration

```bash
# Load parameters from file
python scripts/02_segment_cells.py input.tif output/ @configs/my_config.txt

# Or use xargs
cat configs/my_config.txt | xargs python scripts/02_segment_cells.py input.tif output/
```

---

## Troubleshooting Configurations

### Configuration doesn't work?

1. **Check pixel sizes first**
   - Wrong pixel size → all size-based filtering fails
   - Verify in microscope metadata

2. **Visualize intermediate results**
   - Check `*_report.png` files
   - Identify which step is failing

3. **Adjust one parameter at a time**
   - Start with default config
   - Change one parameter
   - Check result
   - Repeat

4. **Measure ground truth**
   - Manually measure cell/nucleus sizes
   - Compare with detected sizes
   - Adjust size thresholds accordingly

5. **Check image quality**
   - Noisy → increase smoothing
   - Blurry → decrease smoothing
   - Low contrast → decrease thresholds

---

## Need Help?

1. See `docs/PARAMETER_TUNING_GUIDE.md` for detailed parameter explanations
2. See `docs/TROUBLESHOOTING.md` for common issues
3. Run scripts with `--help-params` for built-in guidance
4. Open an issue on GitHub with your configuration and results
