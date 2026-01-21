#!/bin/bash
# =============================================================================
# 00_index_fasta.sh
# Create FASTA Index for O(1) Random Access Lookups
# =============================================================================
#
# Description:
#   Creates an index file (.fai) that enables instant O(1) sequence lookups
#   instead of O(n) linear scanning. Essential for large FASTA files.
#
# Prerequisites:
#   - samtools installed (conda install -c bioconda samtools)
#
# Input:
#   - FASTA file (uncompressed or bgzip-compressed)
#
# Output:
#   - Index file: <fasta>.fai
#   - For bgzip files: also creates <fasta>.gzi
#
# Usage:
#   ./00_index_fasta.sh rep-seqs.fasta
#   ./00_index_fasta.sh rep-seqs.fasta.gz      # bgzip-compressed
#   ./00_index_fasta.sh rep-seqs.fasta --bgzip # compress then index
#
# When to use indexed lookups:
#   - FASTA is large (100k+ sequences, >100 MB)
#   - Repeated lookups on the same file
#   - Batch operations ("give me these 500 ASVs")
#   - Scripts/pipelines with multiple queries
#
# When NOT to use (stick with seqkit grep):
#   - One-off ad-hoc searches
#   - Pattern/regex searches in sequence content
#   - Tiny FASTA files (<10k sequences)
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

FASTA_FILE="${1:-}"
OPTION="${2:-}"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ -z "$FASTA_FILE" ]]; then
    echo "Usage: $0 <fasta_file> [--bgzip]"
    echo ""
    echo "Examples:"
    echo "  $0 rep-seqs.fasta           # Index uncompressed FASTA"
    echo "  $0 rep-seqs.fasta.gz        # Index bgzip-compressed FASTA"
    echo "  $0 rep-seqs.fasta --bgzip   # Compress with bgzip, then index"
    echo ""
    echo "Options:"
    echo "  --bgzip    Compress the FASTA with bgzip before indexing"
    echo "             (saves disk space while keeping random access)"
    echo ""
    echo "Output files:"
    echo "  <fasta>.fai    Index file (enables O(1) lookups)"
    echo "  <fasta>.gzi    Block index (for bgzip files only)"
    exit 1
fi

# Check samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found. Please install it:" >&2
    echo "  conda install -c bioconda samtools" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Handle --bgzip option
# -----------------------------------------------------------------------------

if [[ "$OPTION" == "--bgzip" ]]; then
    if [[ ! -f "$FASTA_FILE" ]]; then
        echo "Error: FASTA file not found: $FASTA_FILE" >&2
        exit 1
    fi

    if [[ "$FASTA_FILE" == *.gz ]]; then
        echo "Error: File is already compressed: $FASTA_FILE" >&2
        echo "If it's gzip (not bgzip), decompress first: gzip -d $FASTA_FILE" >&2
        exit 1
    fi

    # Check bgzip is available
    if ! command -v bgzip &> /dev/null; then
        echo "Error: bgzip not found. Please install samtools/htslib:" >&2
        echo "  conda install -c bioconda samtools" >&2
        exit 1
    fi

    echo "Compressing with bgzip: $FASTA_FILE"
    bgzip "$FASTA_FILE"
    FASTA_FILE="${FASTA_FILE}.gz"
    echo "Created: $FASTA_FILE"
fi

# -----------------------------------------------------------------------------
# Validate FASTA file exists
# -----------------------------------------------------------------------------

if [[ ! -f "$FASTA_FILE" ]]; then
    # Check for compressed version
    if [[ -f "${FASTA_FILE}.zst" ]]; then
        echo "Error: FASTA is zstd-compressed. Decompress first:" >&2
        echo "  zstd -d ${FASTA_FILE}.zst" >&2
        echo ""
        echo "Then optionally recompress with bgzip for indexed access:" >&2
        echo "  $0 ${FASTA_FILE} --bgzip" >&2
        exit 1
    fi
    echo "Error: FASTA file not found: $FASTA_FILE" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Create Index
# -----------------------------------------------------------------------------

echo "Creating index for: $FASTA_FILE"

# Run samtools faidx
if samtools faidx "$FASTA_FILE" 2>&1; then
    echo ""
    echo "Index created successfully!"
    echo ""

    # Show what was created
    if [[ "$FASTA_FILE" == *.gz ]]; then
        echo "Files created:"
        echo "  ${FASTA_FILE}.fai  (sequence index)"
        echo "  ${FASTA_FILE}.gzi  (block index for bgzip)"
    else
        echo "File created:"
        echo "  ${FASTA_FILE}.fai"
    fi

    # Show index stats
    echo ""
    echo "Index stats:"
    NUM_SEQS=$(wc -l < "${FASTA_FILE}.fai")
    INDEX_SIZE=$(ls -lh "${FASTA_FILE}.fai" | awk '{print $5}')
    echo "  Sequences indexed: $NUM_SEQS"
    echo "  Index size: $INDEX_SIZE"

    echo ""
    echo "You can now use fast O(1) lookups:"
    echo "  samtools faidx $FASTA_FILE ASV_ID"
    echo "  ./03b_find_sequence_indexed.sh $FASTA_FILE ASV_ID"
else
    echo ""
    echo "Error: Index creation failed." >&2
    echo ""
    echo "Common issues:" >&2
    echo "  1. Inconsistent line lengths within sequences" >&2
    echo "     Fix: seqkit seq -w 60 input.fasta > reformatted.fasta" >&2
    echo ""
    echo "  2. File is gzip (not bgzip) compressed" >&2
    echo "     Fix: gzip -d file.fasta.gz && bgzip file.fasta" >&2
    exit 1
fi
