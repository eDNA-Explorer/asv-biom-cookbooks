#!/bin/bash
# =============================================================================
# preprocess_reverse_fasta.sh
# Transform reverse FASTA headers to match BIOM feature IDs
# =============================================================================
#
# Description:
#   eDNA Explorer BIOM files use forward read names as feature IDs:
#       16S_Bacteria_paired_F_0
#
#   Forward FASTA headers match directly:
#       >16S_Bacteria_paired_F_0
#
#   Reverse FASTA headers use _paired_R_:
#       >16S_Bacteria_paired_R_0
#
#   This script transforms reverse FASTA headers to match BIOM feature IDs
#   by replacing _paired_R_ with _paired_F_.
#
# Prerequisites:
#   - sed or awk (included in most Unix systems)
#
# Usage:
#   ./preprocess_reverse_fasta.sh <input.fasta> <output.fasta>
#
# Example:
#   ./preprocess_reverse_fasta.sh \
#       fasta/16S_Bacteria_paired_R.fasta \
#       fasta/16S_Bacteria_paired_R_transformed.fasta
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Help
# -----------------------------------------------------------------------------

show_help() {
    cat << EOF
Reverse FASTA Preprocessor

Transforms FASTA headers from _paired_R_ to _paired_F_ format
to match eDNA Explorer BIOM feature IDs.

Usage:
    $0 <input.fasta> <output.fasta>

Arguments:
    input.fasta     Input reverse FASTA file
    output.fasta    Output transformed FASTA file

Example:
    $0 fasta/16S_Bacteria_paired_R.fasta fasta/16S_Bacteria_paired_R_transformed.fasta

Transformation:
    Before: >16S_Bacteria_paired_R_0
    After:  >16S_Bacteria_paired_F_0

EOF
}

# -----------------------------------------------------------------------------
# Parse Arguments
# -----------------------------------------------------------------------------

if [[ $# -eq 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    show_help
    exit 0
fi

if [[ $# -ne 2 ]]; then
    echo "Error: Expected 2 arguments (input and output paths)"
    echo ""
    echo "Usage: $0 <input.fasta> <output.fasta>"
    exit 1
fi

INPUT_FASTA="$1"
OUTPUT_FASTA="$2"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "Error: Input file not found: $INPUT_FASTA"
    exit 1
fi

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_FASTA")
if [[ ! -d "$OUTPUT_DIR" ]] && [[ "$OUTPUT_DIR" != "." ]]; then
    mkdir -p "$OUTPUT_DIR"
fi

# -----------------------------------------------------------------------------
# Process
# -----------------------------------------------------------------------------

echo "============================================================"
echo "  Reverse FASTA Preprocessor (Bash)"
echo "============================================================"
echo ""
echo "Input:  $INPUT_FASTA"
echo "Output: $OUTPUT_FASTA"
echo ""

# Show example transformation
FIRST_HEADER=$(head -1 "$INPUT_FASTA")
if [[ "$FIRST_HEADER" == ">"* ]]; then
    TRANSFORMED=$(echo "$FIRST_HEADER" | sed 's/_paired_R_/_paired_F_/g')
    echo "Transformation:"
    echo "  Before: $FIRST_HEADER"
    echo "  After:  $TRANSFORMED"
    echo ""
fi

# Count sequences before processing
SEQ_COUNT=$(grep -c "^>" "$INPUT_FASTA" || echo "0")
echo "Processing $SEQ_COUNT sequences..."

# Method 1: Using sed (simplest, works on most systems)
# Only transforms lines starting with > (headers)
sed 's/_paired_R_/_paired_F_/g' "$INPUT_FASTA" > "$OUTPUT_FASTA"

# Method 2: Using awk (more control, same result)
# Uncomment to use instead of sed:
# awk '{
#     if (/^>/) {
#         gsub(/_paired_R_/, "_paired_F_")
#     }
#     print
# }' "$INPUT_FASTA" > "$OUTPUT_FASTA"

# Method 3: Using perl (if available)
# Uncomment to use instead of sed:
# perl -pe 's/_paired_R_/_paired_F_/g if /^>/' "$INPUT_FASTA" > "$OUTPUT_FASTA"

# Verify output
OUTPUT_SEQ_COUNT=$(grep -c "^>" "$OUTPUT_FASTA" || echo "0")

if [[ "$SEQ_COUNT" -ne "$OUTPUT_SEQ_COUNT" ]]; then
    echo "Warning: Sequence count mismatch!"
    echo "  Input:  $SEQ_COUNT"
    echo "  Output: $OUTPUT_SEQ_COUNT"
fi

# Count how many were actually transformed
TRANSFORMED_COUNT=$(diff <(grep "^>" "$INPUT_FASTA") <(grep "^>" "$OUTPUT_FASTA") | grep "^<" | wc -l | tr -d ' ')

echo ""
echo "Statistics:"
echo "  Total sequences:  $SEQ_COUNT"
echo "  Transformed:      $TRANSFORMED_COUNT"
echo "  Unchanged:        $((SEQ_COUNT - TRANSFORMED_COUNT))"
echo ""
echo "Output written to: $OUTPUT_FASTA"

if [[ $((SEQ_COUNT - TRANSFORMED_COUNT)) -gt 0 ]]; then
    echo ""
    echo "Note: $((SEQ_COUNT - TRANSFORMED_COUNT)) sequence(s) did not contain '_paired_R_'"
    echo "      These headers were written unchanged."
fi
