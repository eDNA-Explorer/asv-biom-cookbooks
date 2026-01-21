#!/bin/bash
# =============================================================================
# 03_find_sequence.sh
# Find ASV Sequence by ID Using seqkit
# =============================================================================
#
# Description:
#   Quickly retrieve the DNA sequence for a specific ASV from a FASTA file.
#   Replaces: grep "ASV123" giant.tsv (to get sequence)
#
# Prerequisites:
#   - seqkit installed (conda install -c bioconda seqkit)
#
# Input:
#   - FASTA file (e.g., "16S_Bacteria_paired_F.fasta")
#   - ASV ID (e.g., "16S_Bacteria_paired_F_123")
#
# Output:
#   - FASTA record or sequence only
#
# Usage:
#   ./03_find_sequence.sh rep-seqs.fasta 16S_Bacteria_paired_F_123
#   ./03_find_sequence.sh rep-seqs.fasta 16S_Bacteria_paired_F_123 --seq-only
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

FASTA_FILE="${1:-}"
ASV_ID="${2:-}"
SEQ_ONLY="${3:-}"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ -z "$FASTA_FILE" ]] || [[ -z "$ASV_ID" ]]; then
    echo "Usage: $0 <fasta_file> <asv_id> [--seq-only]"
    echo ""
    echo "Examples:"
    echo "  $0 rep-seqs.fasta 16S_Bacteria_paired_F_123"
    echo "  $0 rep-seqs.fasta 16S_Bacteria_paired_F_123 --seq-only"
    echo ""
    echo "Options:"
    echo "  --seq-only    Output sequence only (no header)"
    exit 1
fi

if [[ ! -f "$FASTA_FILE" ]]; then
    # Check for compressed version
    if [[ -f "${FASTA_FILE}.zst" ]]; then
        echo "Error: FASTA file is compressed. Decompress first:" >&2
        echo "  zstd -d ${FASTA_FILE}.zst" >&2
        exit 1
    fi
    echo "Error: FASTA file not found: $FASTA_FILE" >&2
    exit 1
fi

# Check seqkit is available
if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit not found. Please install it:" >&2
    echo "  conda install -c bioconda seqkit" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Search
# -----------------------------------------------------------------------------

# Use regex anchor for exact match
PATTERN="^${ASV_ID}$"

if [[ "$SEQ_ONLY" == "--seq-only" ]]; then
    # Output sequence only (no header, single line)
    seqkit grep -n -r -p "$PATTERN" "$FASTA_FILE" 2>/dev/null | seqkit seq -s -w 0
else
    # Output full FASTA record
    seqkit grep -n -r -p "$PATTERN" "$FASTA_FILE" 2>/dev/null
fi

# Check if anything was found
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo "No sequence found for: $ASV_ID" >&2
    exit 1
fi
