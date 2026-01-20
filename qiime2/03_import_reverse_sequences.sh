#!/bin/bash
# =============================================================================
# 03_import_reverse_sequences.sh
# Import Reverse FASTA Sequences into QIIME2 (Preprocessing Required)
# =============================================================================
#
# Description:
#   Import reverse FASTA sequences from eDNA Explorer into QIIME2.
#
#   IMPORTANT: Reverse FASTA headers use "_paired_R_" while BIOM feature IDs
#   use "_paired_F_". This script automatically transforms the headers before
#   import.
#
# Prerequisites:
#   - QIIME2 environment activated
#   - Python 3 available (for preprocessing)
#
# Input:
#   - Reverse FASTA file (e.g., "fasta/16S_Bacteria_paired_R.fasta")
#
# Output:
#   - QIIME2 artifact: rep-seqs-reverse.qza
#   - Visualization: rep-seqs-reverse.qzv
#   - Transformed FASTA: (temporary, deleted after import)
#
# Usage:
#   ./03_import_reverse_sequences.sh
#   ./03_import_reverse_sequences.sh fasta/my-marker_paired_R.fasta
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Input FASTA file (can be overridden via command line argument)
REVERSE_FASTA="${1:-fasta/16S_Bacteria_paired_R.fasta}"

# Output directory
OUTPUT_DIR="qiime2_output"

# Output file names
OUTPUT_QZA="${OUTPUT_DIR}/rep-seqs-reverse.qza"
OUTPUT_QZV="${OUTPUT_DIR}/rep-seqs-reverse.qzv"

# Temporary transformed FASTA
TEMP_FASTA="${OUTPUT_DIR}/.reverse_transformed_temp.fasta"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

echo "=== QIIME2 Reverse Sequence Import ==="
echo ""

# Check QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "Error: QIIME2 not found. Please activate your QIIME2 environment:"
    echo "  conda activate qiime2-2024.10"
    exit 1
fi

# Check input file exists
if [[ ! -f "$REVERSE_FASTA" ]]; then
    echo "Error: FASTA file not found: $REVERSE_FASTA"
    echo ""
    echo "Usage: $0 [fasta_file]"
    echo "Example: $0 fasta/16S_Bacteria_paired_R.fasta"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# Preprocess: Transform Headers
# -----------------------------------------------------------------------------

echo "Input:  $REVERSE_FASTA"
echo "Output: $OUTPUT_QZA"
echo ""

# Count sequences in input
SEQ_COUNT=$(grep -c "^>" "$REVERSE_FASTA" || echo "0")
echo "Sequences in input: $SEQ_COUNT"

echo ""
echo "Step 1: Transforming FASTA headers..."
echo "  _paired_R_ → _paired_F_"

# Show example transformation
EXAMPLE_HEADER=$(head -1 "$REVERSE_FASTA")
EXAMPLE_TRANSFORMED=$(echo "$EXAMPLE_HEADER" | sed 's/_paired_R_/_paired_F_/g')
echo ""
echo "  Example:"
echo "    Before: $EXAMPLE_HEADER"
echo "    After:  $EXAMPLE_TRANSFORMED"
echo ""

# Transform headers using awk (handles large files efficiently)
awk '{
    if (/^>/) {
        gsub(/_paired_R_/, "_paired_F_")
    }
    print
}' "$REVERSE_FASTA" > "$TEMP_FASTA"

echo "✓ Headers transformed"

# -----------------------------------------------------------------------------
# Import Transformed Sequences
# -----------------------------------------------------------------------------

echo ""
echo "Step 2: Importing into QIIME2..."

qiime tools import \
    --input-path "$TEMP_FASTA" \
    --output-path "$OUTPUT_QZA" \
    --type 'FeatureData[Sequence]'

echo "✓ Sequences imported: $OUTPUT_QZA"

# -----------------------------------------------------------------------------
# Cleanup Temporary File
# -----------------------------------------------------------------------------

rm -f "$TEMP_FASTA"
echo "✓ Temporary file cleaned up"

# -----------------------------------------------------------------------------
# Generate Sequence Visualization
# -----------------------------------------------------------------------------

echo ""
echo "Step 3: Generating sequence summary..."

qiime feature-table tabulate-seqs \
    --i-data "$OUTPUT_QZA" \
    --o-visualization "$OUTPUT_QZV"

echo "✓ Sequence visualization: $OUTPUT_QZV"

# -----------------------------------------------------------------------------
# Display Summary Info
# -----------------------------------------------------------------------------

echo ""
echo "=== Import Complete ==="
echo ""
echo "Artifact info:"
qiime tools peek "$OUTPUT_QZA"

echo ""
echo "View sequences with:"
echo "  qiime tools view $OUTPUT_QZV"
echo ""
echo "Note: Reverse sequence headers have been transformed to match"
echo "      BIOM feature IDs (using _paired_F_ format)."
