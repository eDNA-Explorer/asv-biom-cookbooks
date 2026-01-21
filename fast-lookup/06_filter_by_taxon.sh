#!/bin/bash
# =============================================================================
# 06_filter_by_taxon.sh
# Filter BIOM by Taxonomic Group (e.g., All Escherichia)
# =============================================================================
#
# Description:
#   Find all ASVs matching a taxonomic term and extract their counts.
#   Replaces: grep "g__Escherichia" giant.tsv | cut columns...
#
# Workflow:
#   1. Search taxonomy.tsv for matching ASV IDs
#   2. Subset BIOM to those features
#   3. Optionally convert to TSV
#
# Prerequisites:
#   - ripgrep, biom-format installed
#   - taxonomy.tsv exported (run 02_export_taxonomy.sh first)
#
# Input:
#   - Taxonomy TSV file
#   - BIOM file
#   - Taxonomic search term (e.g., "g__Escherichia", "p__Proteobacteria")
#
# Output:
#   - Subset BIOM file (and optional TSV)
#
# Usage:
#   ./06_filter_by_taxon.sh "g__Escherichia" feature-table.biom taxonomy.tsv
#   ./06_filter_by_taxon.sh "p__Proteobacteria" feature-table.biom taxonomy.tsv --to-tsv
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

TAXON="${1:-}"
BIOM_FILE="${2:-}"
TAXONOMY_FILE="${3:-taxonomy.tsv}"
TO_TSV="${4:-}"

# Output file names (derived from taxon)
SAFE_TAXON=$(echo "$TAXON" | sed 's/[^a-zA-Z0-9_]/_/g')
OUTPUT_IDS="${SAFE_TAXON}_ids.txt"
OUTPUT_BIOM="${SAFE_TAXON}.biom"
OUTPUT_TSV="${SAFE_TAXON}.tsv"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ -z "$TAXON" ]] || [[ -z "$BIOM_FILE" ]]; then
    echo "Usage: $0 <taxon> <biom_file> [taxonomy.tsv] [--to-tsv]"
    echo ""
    echo "Examples:"
    echo "  # Get all Escherichia ASVs"
    echo "  $0 'g__Escherichia' feature-table.biom taxonomy.tsv"
    echo ""
    echo "  # Get all Proteobacteria and export to TSV"
    echo "  $0 'p__Proteobacteria' feature-table.biom taxonomy.tsv --to-tsv"
    echo ""
    echo "  # Search by family"
    echo "  $0 'f__Enterobacteriaceae' feature-table.biom taxonomy.tsv"
    echo ""
    echo "Taxonomy prefixes:"
    echo "  k__ = Kingdom    p__ = Phylum     c__ = Class"
    echo "  o__ = Order      f__ = Family     g__ = Genus     s__ = Species"
    exit 1
fi

if [[ ! -f "$TAXONOMY_FILE" ]]; then
    echo "Error: Taxonomy file not found: $TAXONOMY_FILE" >&2
    echo "Run 02_export_taxonomy.sh first to create taxonomy.tsv" >&2
    exit 1
fi

if [[ ! -f "$BIOM_FILE" ]]; then
    echo "Error: BIOM file not found: $BIOM_FILE" >&2
    exit 1
fi

# Check tools
if ! command -v rg &> /dev/null; then
    echo "Error: ripgrep not found. Please install it:" >&2
    echo "  conda install -c conda-forge ripgrep" >&2
    exit 1
fi

if ! command -v biom &> /dev/null; then
    echo "Error: biom-format not found. Please install it:" >&2
    echo "  conda install -c bioconda biom-format" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Step 1: Find Matching ASV IDs
# -----------------------------------------------------------------------------

echo "=== Filter by Taxon: $TAXON ==="
echo ""

echo "Step 1: Finding matching ASV IDs..."
rg "$TAXON" "$TAXONOMY_FILE" | cut -f1 > "$OUTPUT_IDS"

MATCH_COUNT=$(wc -l < "$OUTPUT_IDS" | tr -d ' ')

if [[ "$MATCH_COUNT" -eq 0 ]]; then
    echo "No ASVs found matching: $TAXON" >&2
    rm -f "$OUTPUT_IDS"
    exit 1
fi

echo "  Found $MATCH_COUNT ASVs matching '$TAXON'"
echo "  Saved to: $OUTPUT_IDS"

# -----------------------------------------------------------------------------
# Step 2: Subset BIOM
# -----------------------------------------------------------------------------

echo ""
echo "Step 2: Subsetting BIOM file..."

biom subset-table \
    -i "$BIOM_FILE" \
    -o "$OUTPUT_BIOM" \
    -a observation \
    --ids-fp "$OUTPUT_IDS"

echo "  Created: $OUTPUT_BIOM"

# Show summary
echo ""
echo "Step 3: Summary..."
biom summarize-table -i "$OUTPUT_BIOM"

# -----------------------------------------------------------------------------
# Step 3 (Optional): Convert to TSV
# -----------------------------------------------------------------------------

if [[ "$TO_TSV" == "--to-tsv" ]] || [[ "$4" == "--to-tsv" ]]; then
    echo ""
    echo "Step 4: Converting to TSV..."
    biom convert -i "$OUTPUT_BIOM" -o "$OUTPUT_TSV" --to-tsv
    echo "  Created: $OUTPUT_TSV"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

echo ""
echo "=== Complete ==="
echo ""
echo "Output files:"
echo "  Feature IDs: $OUTPUT_IDS"
echo "  BIOM file:   $OUTPUT_BIOM"
if [[ -f "$OUTPUT_TSV" ]]; then
    echo "  TSV file:    $OUTPUT_TSV"
fi
echo ""
echo "Quick commands:"
echo "  # View BIOM summary"
echo "  biom summarize-table -i $OUTPUT_BIOM"
echo ""
echo "  # Convert to TSV (if not already)"
echo "  biom convert -i $OUTPUT_BIOM -o $OUTPUT_TSV --to-tsv"
echo ""
echo "  # Get sequences for these ASVs"
echo "  seqkit grep -n -f $OUTPUT_IDS rep-seqs.fasta"
