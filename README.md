# smFISH Multi-Channel Analysis Pipeline

A comprehensive Python pipeline for analyzing single-molecule fluorescence in situ hybridization (smFISH) microscopy data with multi-channel 3D imaging.

## Features

- **Multi-channel signal separation**: Separates CFW (cell wall) and DAPI (nucleus) from mixed channels using shape-based features
- **Robust cell segmentation**: Handles incomplete cell walls using edge detection + watershed algorithm
- **Sensitive nucleus detection**: Detects weak DAPI signals using 3D Gaussian filtering + sphericity filtering
- **Accurate RNA quantification**: Detects single RNA molecules using LoG blob detection
- **Beautiful 3D visualization**: Creates interactive 3D visualizations with smooth surfaces using Marching Cubes algorithm
- **Comprehensive documentation**: Detailed parameter guides, configuration examples, and troubleshooting

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Documentation](#documentation)
- [Examples](#examples)
- [Citation](#citation)
- [License](#license)

## Installation

### Requirements

- Python 3.8 or higher
- 4GB+ RAM (8GB+ recommended for large images)
- Modern web browser (for 3D visualization)

### Install from source

```bash
# Clone repository
git clone https://github.com/yourusername/smfish-analysis.git
cd smfish-analysis

# Create virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\\Scripts\\activate

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

Core packages:
- `aicspylibczi`: CZI file reading
- `scikit-image`: Image processing
- `scipy`: Scientific computing
- `numpy`: Numerical computing
- `pandas`: Data manipulation
- `plotly`: Interactive 3D visualization

See `requirements.txt` for complete list.

## Quick Start

### 1. Prepare your data

Organize your CZI file:
```
my_project/
├── input.czi          # Multi-channel 3D microscopy image
└── output/            # Output directory (will be created)
```

### 2. Run the full pipeline

```bash
# Make script executable
chmod +x run_pipeline.sh

# Run with default parameters
./run_pipeline.sh input.czi output/
```

This will:
1. Extract channels from CZI
2. Separate CFW and DAPI signals
3. Segment cells
4. Detect nuclei
5. Detect RNA molecules
6. Generate 3D visualizations

### 3. View results

Open the generated HTML files in your browser:
- `output/visualization/cell_*.html`: Interactive 3D visualizations
- `output/segmentation/segmentation_report.png`: Cell segmentation quality check
- `output/nuclei/nucleus_detection_report.png`: Nucleus detection quality check
- `output/rna/rna_statistics.csv`: RNA counts per cell

## Pipeline Overview

### Step 1: Channel Extraction

Extract individual channels from multi-channel CZI file.

```bash
python scripts/01_extract_channels.py input.czi output/channels/
```

**Output**:
- `channel_cy3.tif`: RNA signal (Cy3)
- `channel_mixed.tif`: CFW + DAPI mixed signal
- `channel_cfw.tif`: CFW (cell wall) signal (separated)
- `channel_dapi.tif`: DAPI (nucleus) signal (separated)

### Step 2: Cell Segmentation

Segment individual cells from CFW channel using edge detection + watershed.

```bash
python scripts/02_segment_cells.py \\
    output/channels/channel_cfw.tif \\
    output/segmentation/ \\
    --pixel-size-xy 0.44 \\
    --pixel-size-z 0.50
```

**Key parameters**:
- `--pixel-size-xy`, `--pixel-size-z`: Physical pixel dimensions (μm/pixel)
- `--canny-sigma`: Edge smoothing (default: 2.0)
- `--min-distance`: Minimum distance between cell centers (default: 5 pixels)
- `--min-cell-size`, `--max-cell-size`: Cell size filters (voxels)

**Output**:
- `cell_labels.tif`: Labeled cell mask
- `cell_info.json`: Cell metadata (positions, volumes)
- `segmentation_report.png`: Quality check visualization

### Step 3: Nucleus Detection

Detect nuclei within each cell from DAPI channel.

```bash
python scripts/03_detect_nuclei.py \\
    output/channels/channel_dapi.tif \\
    output/segmentation/cell_labels.tif \\
    output/segmentation/cell_info.json \\
    output/nuclei/ \\
    --pixel-size-xy 0.44 \\
    --pixel-size-z 0.50 \\
    --nucleus-diameter-min 5.0 \\
    --nucleus-diameter-max 10.0
```

**Key parameters**:
- `--nucleus-diameter-min`, `--nucleus-diameter-max`: Expected nucleus size (μm)
- `--min-intensity`, `--max-intensity`: DAPI intensity range
- `--min-sphericity`: Shape filter (0-1, higher = more spherical)
- `--gaussian-sigma`: Smoothing strength (default: 2.0)

**Output**:
- `nucleus_mask.tif`: Binary nucleus mask
- `nucleus_info.json`: Nucleus metadata (positions, volumes)
- `nucleus_detection_report.png`: Quality check visualization

### Step 4: RNA Detection

Detect individual RNA molecules from Cy3 channel using LoG blob detection.

```bash
python scripts/04_detect_rna.py \\
    output/channels/channel_cy3.tif \\
    output/segmentation/cell_labels.tif \\
    output/nuclei/nucleus_mask.tif \\
    output/rna/ \\
    --pixel-size-xy 0.44 \\
    --pixel-size-z 0.50
```

**Key parameters**:
- `--min-sigma`, `--max-sigma`: RNA spot size range (pixels)
- `--threshold`: Detection sensitivity (0.001-0.1)
- `--overlap`: Maximum blob overlap (0-1)

**Output**:
- `rna_molecules.json`: RNA coordinates and assignments
- `rna_statistics.csv`: RNA counts per cell (nuclear vs cytoplasmic)

### Step 5: 3D Visualization

Generate interactive 3D visualizations of cells, nuclei, and RNA.

```bash
python scripts/05_visualize_3d.py \\
    output/segmentation/cell_labels.tif \\
    output/nuclei/nucleus_info.json \\
    output/rna/rna_molecules.json \\
    output/visualization/ \\
    --cell-ids 1 2 3 \\
    --pixel-size-xy 0.44 \\
    --pixel-size-z 0.50
```

**Key parameters**:
- `--cell-ids`: Which cells to visualize (space-separated)
- `--cell-smooth`, `--nucleus-smooth`: Surface smoothing (default: 1.5, 1.0)
- `--rna-size`: RNA spot size in visualization (default: 2.0)

**Output**:
- `cell_*.html`: Interactive 3D visualizations (one per cell)

## Documentation

### For Users

- **[Parameter Tuning Guide](docs/PARAMETER_TUNING_GUIDE.md)**: How to optimize parameters for your data
- **[Configuration Examples](examples/CONFIGURATION_EXAMPLES.md)**: Ready-to-use configs for different scenarios
- **[Troubleshooting](docs/TROUBLESHOOTING.md)**: Common issues and solutions

### For Developers

- **[Algorithm Details](docs/ALGORITHM.md)**: Technical description of algorithms
- **[API Reference](docs/API.md)**: Function and class documentation

### Quick Help

Every script has built-in parameter help:

```bash
# General help
python scripts/02_segment_cells.py --help

# Detailed parameter explanations
python scripts/02_segment_cells.py --help-params
```

## Examples

### Example 1: Plant cells with weak DAPI

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --pixel-size-xy 0.44 --pixel-size-z 0.50 \\
    --canny-sigma 2.0 --min-distance 5

python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --pixel-size-xy 0.44 --pixel-size-z 0.50 \\
    --nucleus-diameter-min 5.0 --nucleus-diameter-max 10.0 \\
    --min-intensity 10 --max-intensity 40 \\
    --min-sphericity 0.5
```

### Example 2: Mammalian cells with strong DAPI

```bash
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --pixel-size-xy 0.10 --pixel-size-z 0.30 \\
    --nucleus-diameter-min 8.0 --nucleus-diameter-max 15.0 \\
    --min-intensity 50 --max-intensity 200 \\
    --min-sphericity 0.7
```

### Example 3: Noisy images

```bash
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-sigma 4.0 \\     # More smoothing
    --canny-low 20 \\        # Higher threshold
    --canny-high 50

python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --gaussian-sigma 3.5 \\  # More smoothing
    --min-intensity 20 \\    # Higher threshold
    --min-sphericity 0.6     # Stricter shape filter
```

See [Configuration Examples](examples/CONFIGURATION_EXAMPLES.md) for more scenarios.

## Parameter Tuning Workflow

1. **Check your pixel sizes** (CRITICAL!)
   ```bash
   python -c "import czifile; czi = czifile.CziFile('input.czi'); print(czi.metadata)" | grep Distance
   ```

2. **Run with default parameters**
   ```bash
   ./run_pipeline.sh input.czi output/
   ```

3. **Check quality reports**
   - `segmentation_report.png`: Are cells correctly segmented?
   - `nucleus_detection_report.png`: Are nuclei correctly detected?

4. **Adjust parameters if needed**
   - See [Parameter Tuning Guide](docs/PARAMETER_TUNING_GUIDE.md)
   - Use `--help-params` for detailed explanations

5. **Iterate until satisfied**

## Common Issues

### No cells detected?

```bash
# Try with relaxed parameters
python scripts/02_segment_cells.py input.tif output/ \\
    --canny-low 5 --canny-high 20 \\
    --min-distance 3 \\
    --min-cell-size 500
```

### No nuclei detected?

```bash
# Try with lower threshold
python scripts/03_detect_nuclei.py dapi.tif labels.tif info.json output/ \\
    --min-intensity 5 \\
    --min-sphericity 0.3
```

### Cells are merged?

```bash
# Decrease min_distance for better separation
python scripts/02_segment_cells.py input.tif output/ \\
    --min-distance 3
```

See [Troubleshooting Guide](docs/TROUBLESHOOTING.md) for more solutions.

## Performance Tips

### For large images

```bash
# Process only a region of interest
python scripts/02_segment_cells.py input.tif output/ \\
    --roi-y-min 256 --roi-y-max 512

# Process only a subset of cells
python scripts/04_detect_rna.py ... --process-n-cells 10
```

### For faster visualization

```bash
# Visualize fewer cells
python scripts/05_visualize_3d.py ... --cell-ids 1 2

# Use less smoothing (fewer vertices)
python scripts/05_visualize_3d.py ... \\
    --cell-smooth 1.0 --nucleus-smooth 0.5
```

## Project Structure

```
smfish-analysis/
├── scripts/                    # Analysis scripts
│   ├── 01_extract_channels.py
│   ├── 02_segment_cells.py
│   ├── 03_detect_nuclei.py
│   ├── 04_detect_rna.py
│   └── 05_visualize_3d.py
├── docs/                       # Documentation
│   ├── PARAMETER_TUNING_GUIDE.md
│   ├── TROUBLESHOOTING.md
│   ├── ALGORITHM.md
│   └── API.md
├── examples/                   # Configuration examples
│   └── CONFIGURATION_EXAMPLES.md
├── tests/                      # Unit tests
│   └── test_pipeline.py
├── requirements.txt            # Python dependencies
├── run_pipeline.sh            # Full pipeline script
├── LICENSE                     # MIT License
└── README.md                   # This file
```

## Output Structure

After running the pipeline:

```
output/
├── channels/                   # Extracted channels
│   ├── channel_cy3.tif
│   ├── channel_cfw.tif
│   └── channel_dapi.tif
├── segmentation/              # Cell segmentation results
│   ├── cell_labels.tif
│   ├── cell_info.json
│   └── segmentation_report.png
├── nuclei/                    # Nucleus detection results
│   ├── nucleus_mask.tif
│   ├── nucleus_info.json
│   └── nucleus_detection_report.png
├── rna/                       # RNA detection results
│   ├── rna_molecules.json
│   └── rna_statistics.csv
└── visualization/             # 3D visualizations
    ├── cell_1.html
    ├── cell_2.html
    └── ...
```

## Algorithm Overview

### Cell Segmentation
1. **Edge detection**: Sobel + Canny edge detector
2. **Distance transform**: Find cell centers
3. **Watershed**: Separate touching cells
4. **Size filtering**: Remove debris and merged cells

### Nucleus Detection
1. **Gaussian smoothing**: Reduce noise
2. **Thresholding**: Separate nuclei from background
3. **Connected components**: Label individual nuclei
4. **Sphericity filtering**: Remove non-spherical objects
5. **Size filtering**: Keep only nuclei within expected size range

### RNA Detection
1. **LoG blob detection**: Detect bright spots
2. **Cell assignment**: Assign each RNA to a cell
3. **Nuclear/cytoplasmic classification**: Based on nucleus mask

### 3D Visualization
1. **Marching Cubes**: Generate smooth 3D meshes
2. **Gaussian smoothing**: Further smooth surfaces
3. **Plotly rendering**: Interactive 3D visualization

See [Algorithm Details](docs/ALGORITHM.md) for more information.

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{smfish_analysis,
  title = {smFISH Multi-Channel Analysis Pipeline},
  author = {Zhihe Cai},
  year = {2024},
  url = {https://github.com/IceBear321/smfish-analysis-v2}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **Documentation**: See `docs/` directory
- **Issues**: Open an issue on GitHub
- **Discussions**: Use GitHub Discussions for questions

## Acknowledgments

- **aicspylibczi**: For CZI file reading
- **scikit-image**: For image processing algorithms
- **plotly**: For interactive 3D visualization
- **Community**: Thanks to all contributors and users

---

**Questions?** Check the [documentation](docs/) or open an issue!
