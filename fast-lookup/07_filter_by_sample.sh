#!/bin/bash
# =============================================================================
# 07_filter_by_sample.sh
# Filter BIOM by Sample ID - "What's in this sample?"
# =============================================================================
#
# Description:
#   Extract all features (ASVs) present in a specific sample, optionally
#   sorted by abundance (top hits).
#
# Workflow:
#   1. Subset BIOM to single sample
#   2. Show summary statistics
#   3. Optionally convert to TSV and show top hits
#
# Prerequisites:
#   - biom-format installed (conda install -c bioconda biom-format)
#
# Input:
#   - BIOM file (e.g., "feature-table.biom")
#   - Sample ID (e.g., "SampleA" or "clxyz123abc")
#
# Output:
#   - Subset BIOM file
#   - Optional TSV with top hits
#
# Usage:
#   ./07_filter_by_sample.sh feature-table.biom SampleA
#   ./07_filter_by_sample.sh feature-table.biom SampleA --top 20
#   ./07_filter_by_sample.sh feature-table.biom SampleA --to-tsv
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

BIOM_FILE="${1:-}"
SAMPLE_ID="${2:-}"
shift 2 2>/dev/null || true

# Parse optional arguments
TOP_N=""
TO_TSV=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --top)
            TOP_N="$2"
            shift 2
            ;;
        --to-tsv)
            TO_TSV=true
            shift
            ;;
        *)
            shift
            ;;
    esac
done

# Output file names
SAFE_SAMPLE=$(echo "$SAMPLE_ID" | sed 's/[^a-zA-Z0-9_-]/_/g')
OUTPUT_BIOM="${SAFE_SAMPLE}.biom"
OUTPUT_TSV="${SAFE_SAMPLE}.tsv"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ -z "$BIOM_FILE" ]] || [[ -z "$SAMPLE_ID" ]]; then
    echo "Usage: $0 <biom_file> <sample_id> [--top N] [--to-tsv]"
    echo ""
    echo "Examples:"
    echo "  # Get all features in a sample"
    echo "  $0 feature-table.biom SampleA"
    echo ""
    echo "  # Get top 20 features by abundance"
    echo "  $0 feature-table.biom SampleA --top 20"
    echo ""
    echo "  # Export to TSV"
    echo "  $0 feature-table.biom SampleA --to-tsv"
    echo ""
    echo "  # Combine: top 10 as TSV"
    echo "  $0 feature-table.biom SampleA --top 10 --to-tsv"
    exit 1
fi

if [[ ! -f "$BIOM_FILE" ]]; then
    echo "Error: BIOM file not found: $BIOM_FILE" >&2
    exit 1
fi

# Check biom is available
if ! command -v biom &> /dev/null; then
    echo "Error: biom-format not found. Please install it:" >&2
    echo "  conda install -c bioconda biom-format" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Step 1: Subset BIOM to Sample
# -----------------------------------------------------------------------------

echo "=== Sample Contents: $SAMPLE_ID ==="
echo ""

echo "Step 1: Subsetting to sample..."
biom subset-table \
    -i "$BIOM_FILE" \
    -o "$OUTPUT_BIOM" \
    -a sample \
    --ids "$SAMPLE_ID" 2>/dev/null

if [[ ! -s "$OUTPUT_BIOM" ]]; then
    echo "Error: Sample not found in BIOM file: $SAMPLE_ID" >&2
    exit 1
fi

echo "  Created: $OUTPUT_BIOM"

# -----------------------------------------------------------------------------
# Step 2: Summary
# -----------------------------------------------------------------------------

echo ""
echo "Step 2: Summary..."
biom summarize-table -i "$OUTPUT_BIOM"

# -----------------------------------------------------------------------------
# Step 3: Convert and Show Top Hits
# -----------------------------------------------------------------------------

if [[ "$TO_TSV" == true ]] || [[ -n "$TOP_N" ]]; then
    echo ""
    echo "Step 3: Converting to TSV..."
    biom convert -i "$OUTPUT_BIOM" -o "$OUTPUT_TSV" --to-tsv
    echo "  Created: $OUTPUT_TSV"

    if [[ -n "$TOP_N" ]]; then
        echo ""
        echo "Top $TOP_N features by abundance:"
        echo "─────────────────────────────────────────"
        echo "Feature_ID                                Count"
        echo "─────────────────────────────────────────"
        # Skip header lines, extract feature ID and count, sort by count descending
        tail -n +3 "$OUTPUT_TSV" | \
            awk -F'\t' '{print $1"\t"$2}' | \
            sort -t$'\t' -k2 -rn | \
            head -n "$TOP_N" | \
            awk -F'\t' '{printf "%-40s %s\n", $1, $2}'
    fi
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

echo ""
echo "=== Complete ==="
echo ""
echo "Output files:"
echo "  BIOM file: $OUTPUT_BIOM"
if [[ -f "$OUTPUT_TSV" ]]; then
    echo "  TSV file:  $OUTPUT_TSV"
fi
echo ""
echo "Quick commands:"
echo "  # View full TSV"
echo "  biom convert -i $OUTPUT_BIOM -o - --to-tsv | less"
echo ""
echo "  # Get sequences for features in this sample"
echo "  biom convert -i $OUTPUT_BIOM -o - --to-tsv | tail -n +3 | cut -f1 > ${SAFE_SAMPLE}_features.txt"
echo "  seqkit grep -n -f ${SAFE_SAMPLE}_features.txt rep-seqs.fasta"
