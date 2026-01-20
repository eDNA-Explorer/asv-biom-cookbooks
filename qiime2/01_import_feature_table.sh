#!/bin/bash
# =============================================================================
# 01_import_feature_table.sh
# Import eDNA Explorer BIOM File as QIIME2 FeatureTable
# =============================================================================
#
# Description:
#   Import a BIOM v2.1.0 file from eDNA Explorer into QIIME2 as a
#   FeatureTable[Frequency] artifact.
#
# Prerequisites:
#   - QIIME2 environment activated (conda activate qiime2-2024.10)
#
# Input:
#   - BIOM file (e.g., "biom/16S_Bacteria-asv.biom")
#
# Output:
#   - QIIME2 artifact: feature-table.qza
#   - Visualization: feature-table.qzv (optional)
#
# Usage:
#   ./01_import_feature_table.sh
#   ./01_import_feature_table.sh biom/my-marker-asv.biom
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Input BIOM file (can be overridden via command line argument)
BIOM_FILE="${1:-biom/16S_Bacteria-asv.biom}"

# Output directory
OUTPUT_DIR="qiime2_output"

# Output file names
OUTPUT_QZA="${OUTPUT_DIR}/feature-table.qza"
OUTPUT_QZV="${OUTPUT_DIR}/feature-table.qzv"

# Sample metadata (optional, for visualization)
METADATA_FILE="metadata/sample_metadata.tsv"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

echo "=== QIIME2 Feature Table Import ==="
echo ""

# Check QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "Error: QIIME2 not found. Please activate your QIIME2 environment:"
    echo "  conda activate qiime2-2024.10"
    exit 1
fi

# Check input file exists
if [[ ! -f "$BIOM_FILE" ]]; then
    echo "Error: BIOM file not found: $BIOM_FILE"
    echo ""
    echo "Usage: $0 [biom_file]"
    echo "Example: $0 biom/16S_Bacteria-asv.biom"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# Import Feature Table
# -----------------------------------------------------------------------------

echo "Input:  $BIOM_FILE"
echo "Output: $OUTPUT_QZA"
echo ""

echo "Importing BIOM file..."
qiime tools import \
    --input-path "$BIOM_FILE" \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path "$OUTPUT_QZA"

echo "✓ Feature table imported: $OUTPUT_QZA"

# -----------------------------------------------------------------------------
# Generate Summary Visualization
# -----------------------------------------------------------------------------

echo ""
echo "Generating feature table summary..."

if [[ -f "$METADATA_FILE" ]]; then
    echo "  Including sample metadata from: $METADATA_FILE"
    qiime feature-table summarize \
        --i-table "$OUTPUT_QZA" \
        --m-sample-metadata-file "$METADATA_FILE" \
        --o-visualization "$OUTPUT_QZV"
else
    echo "  No metadata file found, generating basic summary"
    qiime feature-table summarize \
        --i-table "$OUTPUT_QZA" \
        --o-visualization "$OUTPUT_QZV"
fi

echo "✓ Summary visualization: $OUTPUT_QZV"

# -----------------------------------------------------------------------------
# Display Summary Info
# -----------------------------------------------------------------------------

echo ""
echo "=== Import Complete ==="
echo ""
echo "Artifact info:"
qiime tools peek "$OUTPUT_QZA"

echo ""
echo "View the summary visualization with:"
echo "  qiime tools view $OUTPUT_QZV"
echo ""
echo "Next steps:"
echo "  1. Import sequences: ./02_import_forward_sequences.sh"
echo "  2. Import taxonomy: ./04_import_taxonomy.sh"
echo "  3. Run diversity analysis: ./05_diversity_analysis.sh"
