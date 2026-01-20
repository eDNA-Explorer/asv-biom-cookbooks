#!/bin/bash
# =============================================================================
# 02_import_forward_sequences.sh
# Import Forward FASTA Sequences into QIIME2
# =============================================================================
#
# Description:
#   Import forward FASTA sequences from eDNA Explorer into QIIME2.
#   Forward FASTA headers match BIOM feature IDs directly - no preprocessing
#   required.
#
# Prerequisites:
#   - QIIME2 environment activated
#
# Input:
#   - Forward FASTA file (e.g., "fasta/16S_Bacteria_paired_F.fasta")
#
# Output:
#   - QIIME2 artifact: rep-seqs-forward.qza
#   - Visualization: rep-seqs-forward.qzv
#
# Usage:
#   ./02_import_forward_sequences.sh
#   ./02_import_forward_sequences.sh fasta/my-marker_paired_F.fasta
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Input FASTA file (can be overridden via command line argument)
FORWARD_FASTA="${1:-fasta/16S_Bacteria_paired_F.fasta}"

# Output directory
OUTPUT_DIR="qiime2_output"

# Output file names
OUTPUT_QZA="${OUTPUT_DIR}/rep-seqs-forward.qza"
OUTPUT_QZV="${OUTPUT_DIR}/rep-seqs-forward.qzv"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

echo "=== QIIME2 Forward Sequence Import ==="
echo ""

# Check QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "Error: QIIME2 not found. Please activate your QIIME2 environment:"
    echo "  conda activate qiime2-2024.10"
    exit 1
fi

# Check input file exists
if [[ ! -f "$FORWARD_FASTA" ]]; then
    echo "Error: FASTA file not found: $FORWARD_FASTA"
    echo ""
    echo "Usage: $0 [fasta_file]"
    echo "Example: $0 fasta/16S_Bacteria_paired_F.fasta"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# Import Forward Sequences
# -----------------------------------------------------------------------------

echo "Input:  $FORWARD_FASTA"
echo "Output: $OUTPUT_QZA"
echo ""

# Count sequences in input
SEQ_COUNT=$(grep -c "^>" "$FORWARD_FASTA" || echo "0")
echo "Sequences in input: $SEQ_COUNT"

echo ""
echo "Importing forward sequences..."
qiime tools import \
    --input-path "$FORWARD_FASTA" \
    --output-path "$OUTPUT_QZA" \
    --type 'FeatureData[Sequence]'

echo "✓ Sequences imported: $OUTPUT_QZA"

# -----------------------------------------------------------------------------
# Generate Sequence Visualization
# -----------------------------------------------------------------------------

echo ""
echo "Generating sequence summary..."

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
echo "Feature ID Format:"
echo "  Forward FASTA headers match BIOM feature IDs directly."
echo "  Example: >16S_Bacteria_paired_F_0"
