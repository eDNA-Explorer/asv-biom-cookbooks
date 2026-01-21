#!/bin/bash
# =============================================================================
# complete_workflow.sh
# Complete ASV Lookup - Sequence, Taxonomy, and Counts
# =============================================================================
#
# Description:
#   Retrieves all information for a specific ASV in one command:
#   - DNA sequence from FASTA
#   - Taxonomy assignment from taxonomy TSV
#   - Sample counts from BIOM
#
#   Replaces: grep "ASV123" giant.tsv (one-stop lookup)
#
# Prerequisites:
#   - seqkit, ripgrep, biom-format installed
#   - Taxonomy TSV exported (run 02_export_taxonomy.sh first)
#
# Usage:
#   ./complete_workflow.sh <asv_id> [options]
#
# Options:
#   --fasta <file>      FASTA file (default: auto-detect)
#   --taxonomy <file>   Taxonomy TSV (default: taxonomy.tsv)
#   --biom <file>       BIOM file (default: auto-detect *-asv.biom)
#   --dir <directory>   Directory containing all files (default: .)
#
# Examples:
#   ./complete_workflow.sh 16S_Bacteria_paired_F_123
#   ./complete_workflow.sh 16S_Bacteria_paired_F_123 --dir extracted_bundle/
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Parse Arguments
# -----------------------------------------------------------------------------

ASV_ID=""
FASTA_FILE=""
TAXONOMY_FILE="taxonomy.tsv"
BIOM_FILE=""
DATA_DIR="."

while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta)
            FASTA_FILE="$2"
            shift 2
            ;;
        --taxonomy)
            TAXONOMY_FILE="$2"
            shift 2
            ;;
        --biom)
            BIOM_FILE="$2"
            shift 2
            ;;
        --dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 <asv_id> [--fasta <file>] [--taxonomy <file>] [--biom <file>] [--dir <directory>]"
            echo ""
            echo "Retrieves sequence, taxonomy, and counts for an ASV."
            echo ""
            echo "Options:"
            echo "  --fasta <file>      FASTA file (default: auto-detect)"
            echo "  --taxonomy <file>   Taxonomy TSV (default: taxonomy.tsv)"
            echo "  --biom <file>       BIOM file (default: auto-detect *-asv.biom)"
            echo "  --dir <directory>   Directory containing files (default: .)"
            exit 0
            ;;
        -*)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
        *)
            ASV_ID="$1"
            shift
            ;;
    esac
done

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ -z "$ASV_ID" ]]; then
    echo "Error: No ASV ID specified." >&2
    echo "Usage: $0 <asv_id> [options]" >&2
    exit 1
fi

# Check tools
for cmd in seqkit rg biom; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd not found. Run: conda activate fast-lookup" >&2
        exit 1
    fi
done

# Auto-detect files if not specified
if [[ -z "$FASTA_FILE" ]]; then
    FASTA_FILE=$(find "$DATA_DIR" -maxdepth 1 -name "*_F.fasta" -o -name "*_paired_F.fasta" 2>/dev/null | head -1)
fi

if [[ ! -f "$TAXONOMY_FILE" ]] && [[ -f "$DATA_DIR/taxonomy.tsv" ]]; then
    TAXONOMY_FILE="$DATA_DIR/taxonomy.tsv"
fi

if [[ -z "$BIOM_FILE" ]]; then
    BIOM_FILE=$(find "$DATA_DIR" -maxdepth 1 -name "*-asv.biom" 2>/dev/null | head -1)
fi

# -----------------------------------------------------------------------------
# Lookup Functions
# -----------------------------------------------------------------------------

print_separator() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
}

lookup_sequence() {
    echo "SEQUENCE"
    echo "--------"
    if [[ -z "$FASTA_FILE" ]] || [[ ! -f "$FASTA_FILE" ]]; then
        echo "(No FASTA file found)"
        return
    fi
    echo "Source: $FASTA_FILE"
    echo ""

    RESULT=$(seqkit grep -n -r -p "^${ASV_ID}$" "$FASTA_FILE" 2>/dev/null || true)
    if [[ -n "$RESULT" ]]; then
        echo "$RESULT"
    else
        echo "(Not found)"
    fi
}

lookup_taxonomy() {
    echo "TAXONOMY"
    echo "--------"
    if [[ ! -f "$TAXONOMY_FILE" ]]; then
        echo "(No taxonomy.tsv found - run 02_export_taxonomy.sh first)"
        return
    fi
    echo "Source: $TAXONOMY_FILE"
    echo ""

    RESULT=$(rg -m 1 "^${ASV_ID}\t" "$TAXONOMY_FILE" 2>/dev/null || true)
    if [[ -n "$RESULT" ]]; then
        # Parse and format taxonomy
        TAXON=$(echo "$RESULT" | cut -f2)
        echo "Taxonomy: $TAXON"
    else
        echo "(Not found)"
    fi
}

lookup_counts() {
    echo "SAMPLE COUNTS"
    echo "-------------"
    if [[ -z "$BIOM_FILE" ]] || [[ ! -f "$BIOM_FILE" ]]; then
        echo "(No BIOM file found)"
        return
    fi
    echo "Source: $BIOM_FILE"
    echo ""

    # Create temp file for subset
    TMP_BIOM=$(mktemp /tmp/biom_subset.XXXXXX.biom)
    trap "rm -f $TMP_BIOM" RETURN

    # Subset to single ASV
    if ! biom subset-table -i "$BIOM_FILE" -o "$TMP_BIOM" -a observation --ids "$ASV_ID" 2>/dev/null; then
        echo "(ASV not found in BIOM)"
        return
    fi

    # Convert and show non-zero counts
    echo "Sample_ID                         Count"
    echo "--------------------------------------"
    biom convert -i "$TMP_BIOM" -o /dev/stdout --to-tsv 2>/dev/null | \
        tail -n +2 | \
        awk -F'\t' 'NR==1 {for(i=2;i<=NF;i++) samples[i]=$i; next} {for(i=2;i<=NF;i++) if($i>0) printf "%-35s %s\n", samples[i], $i}'

    rm -f "$TMP_BIOM"
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

echo ""
echo "ASV LOOKUP: $ASV_ID"
echo "================================================================================"

print_separator
lookup_sequence

print_separator
lookup_taxonomy

print_separator
lookup_counts

echo ""
echo "================================================================================"
