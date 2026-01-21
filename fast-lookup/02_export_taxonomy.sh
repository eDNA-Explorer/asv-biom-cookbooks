#!/bin/bash
# =============================================================================
# 02_export_taxonomy.sh
# Export Taxonomy from BIOM to TSV for Fast Grep Lookups
# =============================================================================
#
# Description:
#   Exports observation metadata (taxonomy) from a BIOM file to a simple TSV
#   format that can be quickly searched with ripgrep.
#
# Prerequisites:
#   - biom-format installed (conda install -c bioconda biom-format)
#
# Input:
#   - Taxa BIOM file (e.g., "16S_Bacteria-taxa.biom")
#
# Output:
#   - taxonomy.tsv with columns: FeatureID, Taxon, Confidence
#
# Usage:
#   ./02_export_taxonomy.sh 16S_Bacteria-taxa.biom
#   ./02_export_taxonomy.sh 16S_Bacteria-taxa.biom custom_output.tsv
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

BIOM_FILE="${1:-}"
OUTPUT_FILE="${2:-taxonomy.tsv}"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

echo "=== Taxonomy Export ==="
echo ""

if [[ -z "$BIOM_FILE" ]]; then
    echo "Error: No BIOM file specified."
    echo ""
    echo "Usage: $0 <taxa.biom> [output.tsv]"
    echo "Example: $0 16S_Bacteria-taxa.biom"
    exit 1
fi

if [[ ! -f "$BIOM_FILE" ]]; then
    echo "Error: BIOM file not found: $BIOM_FILE"
    exit 1
fi

# Check biom is available
if ! command -v biom &> /dev/null; then
    echo "Error: biom-format not found. Please install it:"
    echo "  conda install -c bioconda biom-format"
    exit 1
fi

# -----------------------------------------------------------------------------
# Export Taxonomy
# -----------------------------------------------------------------------------

echo "Input:  $BIOM_FILE"
echo "Output: $OUTPUT_FILE"
echo ""

echo "Exporting observation metadata..."

# Export to TSV with observation metadata
biom convert \
    -i "$BIOM_FILE" \
    -o "$OUTPUT_FILE" \
    --to-tsv \
    --header-key taxonomy

# Clean up the header (remove the "# Constructed from biom file" line)
if [[ "$(head -1 "$OUTPUT_FILE")" == "# Constructed from biom file" ]]; then
    tail -n +2 "$OUTPUT_FILE" > "${OUTPUT_FILE}.tmp"
    mv "${OUTPUT_FILE}.tmp" "$OUTPUT_FILE"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

LINE_COUNT=$(wc -l < "$OUTPUT_FILE" | tr -d ' ')
echo ""
echo "=== Export Complete ==="
echo ""
echo "Exported $LINE_COUNT features to: $OUTPUT_FILE"
echo ""
echo "Preview (first 5 lines):"
head -5 "$OUTPUT_FILE"
echo ""
echo "Search examples:"
echo "  # Find taxonomy for specific ASV"
echo "  rg '^16S_Bacteria_paired_F_123\t' $OUTPUT_FILE"
echo ""
echo "  # Find all Escherichia"
echo "  rg -i 'g__Escherichia' $OUTPUT_FILE"
echo ""
echo "  # Count features per phylum"
echo "  rg -o 'p__[^;]+' $OUTPUT_FILE | sort | uniq -c | sort -rn"
