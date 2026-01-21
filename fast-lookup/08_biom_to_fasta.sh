#!/bin/bash
# =============================================================================
# 08_biom_to_fasta.sh
# Extract FASTA Sequences from BIOM Filter Results
# =============================================================================
#
# Description:
#   Combines BIOM filtering with indexed FASTA lookup to extract sequences
#   for ASVs matching specific criteria - WITHOUT materializing giant TSV files.
#
#   Workflow: BIOM filter → Extract IDs → Indexed FASTA lookup
#
# Prerequisites:
#   - samtools installed (conda install -c bioconda samtools)
#   - biom-format installed (conda install -c bioconda biom-format)
#   - FASTA index created: ./00_index_fasta.sh rep-seqs.fasta
#
# Input:
#   - BIOM file (feature table)
#   - Indexed FASTA file
#   - Filter criteria (taxonomy term, sample ID, or ID file)
#
# Output:
#   - FASTA file with matching sequences
#
# Usage:
#   # By taxonomy (requires taxonomy.tsv from 02_export_taxonomy.sh)
#   ./08_biom_to_fasta.sh --taxon "g__Bacteroides" table.biom rep-seqs.fasta taxonomy.tsv
#
#   # By sample (ASVs present in a sample)
#   ./08_biom_to_fasta.sh --sample "SampleA" table.biom rep-seqs.fasta
#
#   # By ID file (pre-filtered list)
#   ./08_biom_to_fasta.sh --ids asv_ids.txt rep-seqs.fasta
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Parse Arguments
# -----------------------------------------------------------------------------

MODE=""
FILTER_VALUE=""
BIOM_FILE=""
FASTA_FILE=""
TAXONOMY_FILE=""
OUTPUT_FILE=""

usage() {
    echo "Usage: $0 <mode> <filter> <biom_file> <fasta_file> [taxonomy_file] [-o output.fasta]"
    echo ""
    echo "Modes:"
    echo "  --taxon <term>    Filter by taxonomy (e.g., 'g__Bacteroides')"
    echo "  --sample <id>     Filter by sample ID (ASVs present in sample)"
    echo "  --ids <file>      Use pre-filtered ID file (one ID per line)"
    echo ""
    echo "Examples:"
    echo "  # Get all Bacteroides sequences"
    echo "  $0 --taxon 'g__Bacteroides' table.biom rep-seqs.fasta taxonomy.tsv"
    echo ""
    echo "  # Get sequences for ASVs in SampleA"
    echo "  $0 --sample 'SampleA' table.biom rep-seqs.fasta"
    echo ""
    echo "  # Get sequences from pre-filtered ID list"
    echo "  $0 --ids my_asvs.txt rep-seqs.fasta"
    echo ""
    echo "  # Specify output file"
    echo "  $0 --taxon 'g__Bacteroides' table.biom rep-seqs.fasta taxonomy.tsv -o bacteroides.fasta"
    echo ""
    echo "Prerequisites:"
    echo "  - Create FASTA index: ./00_index_fasta.sh rep-seqs.fasta"
    echo "  - For --taxon mode: ./02_export_taxonomy.sh to create taxonomy.tsv"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --taxon)
            MODE="taxon"
            FILTER_VALUE="$2"
            shift 2
            ;;
        --sample)
            MODE="sample"
            FILTER_VALUE="$2"
            shift 2
            ;;
        --ids)
            MODE="ids"
            FILTER_VALUE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            # Positional arguments
            if [[ -z "$BIOM_FILE" ]] && [[ "$1" == *.biom ]]; then
                BIOM_FILE="$1"
            elif [[ -z "$FASTA_FILE" ]] && [[ "$1" == *.fasta* ]]; then
                FASTA_FILE="$1"
            elif [[ -z "$TAXONOMY_FILE" ]] && [[ "$1" == *.tsv ]]; then
                TAXONOMY_FILE="$1"
            elif [[ -z "$FASTA_FILE" ]] && [[ -f "$1" ]]; then
                # Could be FASTA without .fasta extension or IDs file
                FASTA_FILE="$1"
            else
                echo "Warning: Unknown argument: $1" >&2
            fi
            shift
            ;;
    esac
done

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

if [[ -z "$MODE" ]]; then
    echo "Error: No mode specified (--taxon, --sample, or --ids)" >&2
    usage
fi

if [[ -z "$FILTER_VALUE" ]]; then
    echo "Error: No filter value provided" >&2
    usage
fi

# Mode-specific validation
case "$MODE" in
    taxon)
        if [[ -z "$BIOM_FILE" ]] || [[ -z "$FASTA_FILE" ]] || [[ -z "$TAXONOMY_FILE" ]]; then
            echo "Error: --taxon mode requires: <biom_file> <fasta_file> <taxonomy_file>" >&2
            usage
        fi
        if [[ ! -f "$TAXONOMY_FILE" ]]; then
            echo "Error: Taxonomy file not found: $TAXONOMY_FILE" >&2
            echo "Create it with: ./02_export_taxonomy.sh" >&2
            exit 1
        fi
        ;;
    sample)
        if [[ -z "$BIOM_FILE" ]] || [[ -z "$FASTA_FILE" ]]; then
            echo "Error: --sample mode requires: <biom_file> <fasta_file>" >&2
            usage
        fi
        ;;
    ids)
        if [[ -z "$FASTA_FILE" ]]; then
            # In ids mode, FILTER_VALUE is the ids file, next arg is fasta
            echo "Error: --ids mode requires: <ids_file> <fasta_file>" >&2
            usage
        fi
        if [[ ! -f "$FILTER_VALUE" ]]; then
            echo "Error: IDs file not found: $FILTER_VALUE" >&2
            exit 1
        fi
        ;;
esac

# Check FASTA exists
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "Error: FASTA file not found: $FASTA_FILE" >&2
    exit 1
fi

# Check for FASTA index
if [[ ! -f "${FASTA_FILE}.fai" ]]; then
    echo "Error: FASTA index not found: ${FASTA_FILE}.fai" >&2
    echo "Create it with: ./00_index_fasta.sh $FASTA_FILE" >&2
    exit 1
fi

# Check BIOM exists (if needed)
if [[ -n "$BIOM_FILE" ]] && [[ ! -f "$BIOM_FILE" ]]; then
    echo "Error: BIOM file not found: $BIOM_FILE" >&2
    exit 1
fi

# Check tools
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found" >&2
    exit 1
fi

if [[ "$MODE" != "ids" ]] && ! command -v biom &> /dev/null; then
    echo "Error: biom-format not found" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Create temp directory for intermediate files
# -----------------------------------------------------------------------------

TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

IDS_FILE="$TEMP_DIR/ids.txt"

# -----------------------------------------------------------------------------
# Extract IDs based on mode
# -----------------------------------------------------------------------------

case "$MODE" in
    taxon)
        echo "Filtering by taxonomy: $FILTER_VALUE" >&2

        # Get ASV IDs matching the taxonomy term
        rg "$FILTER_VALUE" "$TAXONOMY_FILE" | cut -f1 > "$IDS_FILE"

        NUM_IDS=$(wc -l < "$IDS_FILE" | tr -d ' ')
        if [[ "$NUM_IDS" -eq 0 ]]; then
            echo "Error: No ASVs found matching '$FILTER_VALUE'" >&2
            exit 1
        fi
        echo "Found $NUM_IDS ASVs matching taxonomy" >&2
        ;;

    sample)
        echo "Filtering by sample: $FILTER_VALUE" >&2

        # Subset BIOM to single sample, then get feature IDs
        TEMP_BIOM="$TEMP_DIR/sample.biom"
        biom subset-table -i "$BIOM_FILE" -o "$TEMP_BIOM" \
            --axis sample --ids "$FILTER_VALUE" 2>/dev/null

        # Extract feature IDs (observations with non-zero counts)
        biom convert -i "$TEMP_BIOM" -o "$TEMP_DIR/sample.tsv" --to-tsv 2>/dev/null
        tail -n +3 "$TEMP_DIR/sample.tsv" | cut -f1 > "$IDS_FILE"

        NUM_IDS=$(wc -l < "$IDS_FILE" | tr -d ' ')
        if [[ "$NUM_IDS" -eq 0 ]]; then
            echo "Error: No ASVs found in sample '$FILTER_VALUE'" >&2
            exit 1
        fi
        echo "Found $NUM_IDS ASVs in sample" >&2
        ;;

    ids)
        echo "Using ID file: $FILTER_VALUE" >&2

        # Copy and clean the IDs file (remove empty lines)
        grep -v '^$' "$FILTER_VALUE" > "$IDS_FILE"

        NUM_IDS=$(wc -l < "$IDS_FILE" | tr -d ' ')
        echo "Found $NUM_IDS ASVs in file" >&2
        ;;
esac

# -----------------------------------------------------------------------------
# Fetch sequences using indexed FASTA
# -----------------------------------------------------------------------------

echo "Fetching sequences from indexed FASTA..." >&2

if [[ -n "$OUTPUT_FILE" ]]; then
    samtools faidx "$FASTA_FILE" -r "$IDS_FILE" > "$OUTPUT_FILE"
    echo "Output written to: $OUTPUT_FILE" >&2
else
    # Output to stdout
    samtools faidx "$FASTA_FILE" -r "$IDS_FILE"
fi

echo "Done! Retrieved sequences for $NUM_IDS ASVs" >&2
