#!/bin/bash
# =============================================================================
# 03b_find_sequence_indexed.sh
# Find ASV Sequence by ID Using Indexed FASTA (O(1) Lookup)
# =============================================================================
#
# Description:
#   Instantly retrieve DNA sequences using a pre-built FASTA index.
#   This is O(1) constant-time lookup vs O(n) linear scan with seqkit grep.
#
#   Use this for:
#   - Large FASTA files (100k+ sequences)
#   - Repeated lookups
#   - Batch operations
#   - Scripts/pipelines
#
#   Use 03_find_sequence.sh (seqkit grep) for:
#   - One-off searches
#   - Pattern/regex matching
#   - Files without an index
#
# Prerequisites:
#   - samtools installed (conda install -c bioconda samtools)
#   - FASTA index created: ./00_index_fasta.sh rep-seqs.fasta
#
# Input:
#   - Indexed FASTA file (with .fai file present)
#   - ASV ID(s) or file of IDs
#
# Output:
#   - FASTA record(s)
#
# Usage:
#   ./03b_find_sequence_indexed.sh rep-seqs.fasta ASV_001
#   ./03b_find_sequence_indexed.sh rep-seqs.fasta ASV_001 ASV_002 ASV_003
#   ./03b_find_sequence_indexed.sh rep-seqs.fasta -f asv_ids.txt
#   ./03b_find_sequence_indexed.sh rep-seqs.fasta ASV_001:1-100  # subsequence
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

FASTA_FILE="${1:-}"
shift 2>/dev/null || true  # Remove first arg, continue if none

# Remaining args are either IDs, -f flag, or region specs
IDS=()
FROM_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--file)
            FROM_FILE="${2:-}"
            shift 2
            ;;
        --seq-only)
            SEQ_ONLY=true
            shift
            ;;
        *)
            IDS+=("$1")
            shift
            ;;
    esac
done

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ -z "$FASTA_FILE" ]]; then
    echo "Usage: $0 <fasta_file> <asv_id> [asv_id2 ...] [options]"
    echo "       $0 <fasta_file> -f <ids_file>"
    echo ""
    echo "Examples:"
    echo "  $0 rep-seqs.fasta ASV_001                    # Single lookup"
    echo "  $0 rep-seqs.fasta ASV_001 ASV_002 ASV_003    # Multiple IDs"
    echo "  $0 rep-seqs.fasta -f asv_ids.txt             # Batch from file"
    echo "  $0 rep-seqs.fasta ASV_001:1-100              # Subsequence (bases 1-100)"
    echo "  $0 rep-seqs.fasta ASV_001 --seq-only         # Sequence only (no header)"
    echo ""
    echo "Options:"
    echo "  -f, --file <file>   Read IDs from file (one per line)"
    echo "  --seq-only          Output sequence only (no FASTA header)"
    echo ""
    echo "Prerequisites:"
    echo "  Create index first: ./00_index_fasta.sh rep-seqs.fasta"
    exit 1
fi

if [[ ! -f "$FASTA_FILE" ]]; then
    echo "Error: FASTA file not found: $FASTA_FILE" >&2
    exit 1
fi

# Check for index
INDEX_FILE="${FASTA_FILE}.fai"
if [[ ! -f "$INDEX_FILE" ]]; then
    echo "Error: Index file not found: $INDEX_FILE" >&2
    echo ""
    echo "Create the index first:" >&2
    echo "  ./00_index_fasta.sh $FASTA_FILE" >&2
    echo ""
    echo "Or use the non-indexed version (slower for large files):" >&2
    echo "  ./03_find_sequence.sh $FASTA_FILE <asv_id>" >&2
    exit 1
fi

# Check samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found. Please install it:" >&2
    echo "  conda install -c bioconda samtools" >&2
    exit 1
fi

# Need either IDs or file
if [[ ${#IDS[@]} -eq 0 ]] && [[ -z "$FROM_FILE" ]]; then
    echo "Error: No ASV IDs provided" >&2
    echo "Usage: $0 <fasta_file> <asv_id> [asv_id2 ...]" >&2
    echo "       $0 <fasta_file> -f <ids_file>" >&2
    exit 1
fi

# Validate file if provided
if [[ -n "$FROM_FILE" ]] && [[ ! -f "$FROM_FILE" ]]; then
    echo "Error: IDs file not found: $FROM_FILE" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Lookup
# -----------------------------------------------------------------------------

run_faidx() {
    if [[ "$SEQ_ONLY" == true ]]; then
        # Strip headers, output sequence only
        samtools faidx "$@" 2>/dev/null | grep -v "^>" | tr -d '\n'
        echo  # Add final newline
    else
        samtools faidx "$@" 2>/dev/null
    fi
}

if [[ -n "$FROM_FILE" ]]; then
    # Batch lookup from file
    # Remove empty lines from input file
    run_faidx "$FASTA_FILE" -r <(grep -v '^$' "$FROM_FILE")
else
    # Lookup specific IDs
    run_faidx "$FASTA_FILE" "${IDS[@]}"
fi

EXIT_CODE=$?

if [[ $EXIT_CODE -ne 0 ]]; then
    echo "Warning: Some sequences may not have been found" >&2
fi

exit $EXIT_CODE
