#!/bin/bash
# =============================================================================
# 01_extract_bundle.sh
# Extract eDNA Explorer BIOM Bundle and Decompress FASTA Files
# =============================================================================
#
# Description:
#   Extracts a tar bundle from eDNA Explorer and decompresses zstd-compressed
#   FASTA files for fast sequence lookups.
#
# Prerequisites:
#   - zstd installed (conda install -c conda-forge zstd)
#
# Input:
#   - Bundle tar file (e.g., "16S_Bacteria-biom.tar")
#
# Output:
#   - Extracted directory with decompressed FASTA files
#
# Usage:
#   ./01_extract_bundle.sh 16S_Bacteria-biom.tar
#   ./01_extract_bundle.sh 16S_Bacteria-biom.tar output_dir
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

BUNDLE_FILE="${1:-}"
OUTPUT_DIR="${2:-}"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

echo "=== Bundle Extraction ==="
echo ""

if [[ -z "$BUNDLE_FILE" ]]; then
    echo "Error: No bundle file specified."
    echo ""
    echo "Usage: $0 <bundle.tar> [output_dir]"
    echo "Example: $0 16S_Bacteria-biom.tar"
    exit 1
fi

if [[ ! -f "$BUNDLE_FILE" ]]; then
    echo "Error: Bundle file not found: $BUNDLE_FILE"
    exit 1
fi

# Check zstd is available
if ! command -v zstd &> /dev/null; then
    echo "Error: zstd not found. Please install it:"
    echo "  conda install -c conda-forge zstd"
    exit 1
fi

# Determine output directory from bundle name if not specified
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="${BUNDLE_FILE%.tar}"
    OUTPUT_DIR="${OUTPUT_DIR%-biom}"
fi

# -----------------------------------------------------------------------------
# Extract Bundle
# -----------------------------------------------------------------------------

echo "Input:  $BUNDLE_FILE"
echo "Output: $OUTPUT_DIR/"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Extracting bundle..."
tar -xf "$BUNDLE_FILE" -C "$OUTPUT_DIR"

echo "Extracted files:"
ls -la "$OUTPUT_DIR/"

# -----------------------------------------------------------------------------
# Decompress FASTA Files
# -----------------------------------------------------------------------------

echo ""
echo "Decompressing FASTA files..."

ZST_COUNT=0
for zst_file in "$OUTPUT_DIR"/*.fasta.zst; do
    if [[ -f "$zst_file" ]]; then
        echo "  Decompressing: $(basename "$zst_file")"
        zstd -d --rm "$zst_file"
        ((ZST_COUNT++))
    fi
done

if [[ $ZST_COUNT -eq 0 ]]; then
    echo "  No .fasta.zst files found (already decompressed?)"
else
    echo "  Decompressed $ZST_COUNT FASTA files"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

echo ""
echo "=== Extraction Complete ==="
echo ""
echo "Directory contents:"
ls -lh "$OUTPUT_DIR/"

echo ""
echo "Next steps:"
echo "  1. Export taxonomy for fast lookup:"
echo "     ./02_export_taxonomy.sh $OUTPUT_DIR/*-taxa.biom"
echo ""
echo "  2. Find sequences:"
echo "     ./03_find_sequence.sh $OUTPUT_DIR/*.fasta <ASV_ID>"
echo ""
echo "  3. Find counts:"
echo "     ./05_find_counts.sh $OUTPUT_DIR/*-asv.biom <ASV_ID>"
