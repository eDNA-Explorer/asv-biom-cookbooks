#!/bin/bash
# =============================================================================
# 05_find_counts.sh
# Find Sample Counts for an ASV Using biom subset-table
# =============================================================================
#
# Description:
#   Quickly retrieve the count data for a specific ASV across all samples.
#   Uses biom subset-table to avoid loading the entire table into memory.
#   Replaces: grep "ASV123" giant.tsv (to get counts per sample)
#
# Prerequisites:
#   - biom-format installed (conda install -c bioconda biom-format)
#
# Input:
#   - BIOM file (e.g., "16S_Bacteria-asv.biom")
#   - ASV ID (e.g., "16S_Bacteria_paired_F_123")
#
# Output:
#   - TSV showing sample counts for the ASV
#
# Usage:
#   ./05_find_counts.sh feature-table.biom 16S_Bacteria_paired_F_123
#   ./05_find_counts.sh feature-table.biom 16S_Bacteria_paired_F_123 --keep-zeros
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

BIOM_FILE="${1:-}"
ASV_ID="${2:-}"
KEEP_ZEROS="${3:-}"

# Temporary file for subset
TMP_BIOM=$(mktemp /tmp/biom_subset.XXXXXX.biom)
trap "rm -f $TMP_BIOM" EXIT

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ -z "$BIOM_FILE" ]] || [[ -z "$ASV_ID" ]]; then
    echo "Usage: $0 <biom_file> <asv_id> [--keep-zeros]"
    echo ""
    echo "Examples:"
    echo "  $0 feature-table.biom 16S_Bacteria_paired_F_123"
    echo "  $0 feature-table.biom 16S_Bacteria_paired_F_123 --keep-zeros"
    echo ""
    echo "Options:"
    echo "  --keep-zeros    Include samples with zero counts (default: non-zero only)"
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
# Subset and Convert
# -----------------------------------------------------------------------------

# Subset BIOM to single ASV
biom subset-table \
    -i "$BIOM_FILE" \
    -o "$TMP_BIOM" \
    -a observation \
    --ids "$ASV_ID" 2>/dev/null

# Check if subset succeeded
if [[ ! -s "$TMP_BIOM" ]]; then
    echo "Error: ASV not found in BIOM file: $ASV_ID" >&2
    exit 1
fi

# Convert to TSV and display
echo "# Counts for: $ASV_ID"
echo "# Sample_ID	Count"

# Convert and format output
biom convert -i "$TMP_BIOM" -o /dev/stdout --to-tsv 2>/dev/null | \
    tail -n +2 | \  # Skip header comments
    if [[ "$KEEP_ZEROS" != "--keep-zeros" ]]; then
        # Filter out zero counts and transpose for readability
        awk -F'\t' 'NR==1 {for(i=2;i<=NF;i++) samples[i]=$i; next} {for(i=2;i<=NF;i++) if($i>0) print samples[i]"\t"$i}'
    else
        # Keep all samples, transpose for readability
        awk -F'\t' 'NR==1 {for(i=2;i<=NF;i++) samples[i]=$i; next} {for(i=2;i<=NF;i++) print samples[i]"\t"$i}'
    fi
