#!/bin/bash

# smFISH Multi-Channel Analysis Pipeline
# Full pipeline script

set -e  # Exit on error

# Check arguments
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input.czi> <output_dir> [options]"
    echo ""
    echo "Example:"
    echo "  $0 data/stack_oil_2.czi output/"
    echo ""
    echo "Optional environment variables:"
    echo "  PIXEL_SIZE_XY=0.44    # μm/pixel (default: 0.44)"
    echo "  PIXEL_SIZE_Z=0.50     # μm/pixel (default: 0.50)"
    echo "  ROI_Y_MIN=256         # Focus on lower region (optional)"
    exit 1
fi

INPUT_CZI="$1"
OUTPUT_DIR="$2"
shift 2

# Default parameters (can be overridden by environment variables)
PIXEL_SIZE_XY="${PIXEL_SIZE_XY:-0.44}"
PIXEL_SIZE_Z="${PIXEL_SIZE_Z:-0.50}"

echo "=========================================="
echo "smFISH Multi-Channel Analysis Pipeline"
echo "=========================================="
echo "Input: $INPUT_CZI"
echo "Output: $OUTPUT_DIR"
echo "Pixel size: ${PIXEL_SIZE_XY} x ${PIXEL_SIZE_XY} x ${PIXEL_SIZE_Z} μm"
echo ""

# Create output directories
mkdir -p "$OUTPUT_DIR"/{channels,segmentation,nuclei,rna,visualization}

# Step 1: Extract channels
echo "Step 1/5: Extracting channels..."
python scripts/01_extract_channels.py "$INPUT_CZI" "$OUTPUT_DIR/channels/"
echo "✓ Channels extracted"
echo ""

# Step 2: Segment cells
echo "Step 2/5: Segmenting cells..."
CMD="python scripts/02_segment_cells.py \
    $OUTPUT_DIR/channels/channel_cfw.tif \
    $OUTPUT_DIR/segmentation/ \
    --pixel-size-xy $PIXEL_SIZE_XY \
    --pixel-size-z $PIXEL_SIZE_Z"

# Add ROI if specified
if [ -n "$ROI_Y_MIN" ]; then
    CMD="$CMD --roi-y-min $ROI_Y_MIN"
fi

# Add any additional arguments
CMD="$CMD $@"

eval $CMD
echo "✓ Cells segmented"
echo ""

# Step 3: Detect nuclei
echo "Step 3/5: Detecting nuclei..."
python scripts/03_detect_nuclei.py \
    "$OUTPUT_DIR/channels/channel_dapi.tif" \
    "$OUTPUT_DIR/segmentation/cell_labels.tif" \
    "$OUTPUT_DIR/segmentation/cell_info.json" \
    "$OUTPUT_DIR/nuclei/" \
    --pixel-size-xy $PIXEL_SIZE_XY \
    --pixel-size-z $PIXEL_SIZE_Z \
    --nucleus-diameter-min 5.0 \
    --nucleus-diameter-max 10.0 \
    --min-intensity 10 \
    --max-intensity 40
echo "✓ Nuclei detected"
echo ""

# Step 4: Detect RNA
echo "Step 4/5: Detecting RNA molecules..."
python scripts/04_detect_rna.py \
    "$OUTPUT_DIR/channels/channel_cy3.tif" \
    "$OUTPUT_DIR/segmentation/cell_labels.tif" \
    "$OUTPUT_DIR/nuclei/nucleus_mask.tif" \
    "$OUTPUT_DIR/rna/" \
    --pixel-size-xy $PIXEL_SIZE_XY \
    --pixel-size-z $PIXEL_SIZE_Z
echo "✓ RNA detected"
echo ""

# Step 5: Generate 3D visualizations (first 5 cells)
echo "Step 5/5: Generating 3D visualizations..."
# Get first 5 cell IDs from cell_info.json
CELL_IDS=$(python -c "
import json
with open('$OUTPUT_DIR/segmentation/cell_info.json') as f:
    cells = json.load(f)
    ids = sorted([int(k) for k in cells.keys()])[:5]
    print(' '.join(map(str, ids)))
")

if [ -n "$CELL_IDS" ]; then
    python scripts/05_visualize_3d.py \
        "$OUTPUT_DIR/segmentation/cell_labels.tif" \
        "$OUTPUT_DIR/nuclei/nucleus_info.json" \
        "$OUTPUT_DIR/rna/rna_molecules.json" \
        "$OUTPUT_DIR/visualization/" \
        --cell-ids $CELL_IDS \
        --pixel-size-xy $PIXEL_SIZE_XY \
        --pixel-size-z $PIXEL_SIZE_Z
    echo "✓ 3D visualizations generated"
else
    echo "⚠ No cells to visualize"
fi
echo ""

echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="
echo "Results:"
echo "  - Channels: $OUTPUT_DIR/channels/"
echo "  - Segmentation: $OUTPUT_DIR/segmentation/"
echo "  - Nuclei: $OUTPUT_DIR/nuclei/"
echo "  - RNA: $OUTPUT_DIR/rna/"
echo "  - Visualization: $OUTPUT_DIR/visualization/"
echo ""
echo "Next steps:"
echo "  1. Check quality reports:"
echo "     - $OUTPUT_DIR/segmentation/segmentation_report.png"
echo "     - $OUTPUT_DIR/nuclei/nucleus_detection_report.png"
echo "  2. View 3D visualizations in browser:"
echo "     - $OUTPUT_DIR/visualization/cell_*.html"
echo "  3. Check RNA statistics:"
echo "     - $OUTPUT_DIR/rna/rna_statistics.csv"
echo ""
echo "For parameter tuning, see docs/PARAMETER_TUNING_GUIDE.md"
