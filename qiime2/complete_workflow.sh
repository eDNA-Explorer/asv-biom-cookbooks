#!/bin/bash
# =============================================================================
# complete_workflow.sh
# Complete QIIME2 Workflow for eDNA Explorer Data
# =============================================================================
#
# Description:
#   End-to-end workflow for importing and analyzing eDNA Explorer BIOM files
#   in QIIME2. This script:
#   1. Imports feature table from BIOM
#   2. Imports forward sequences
#   3. Imports reverse sequences (with preprocessing)
#   4. Imports taxonomy
#   5. Generates taxonomic bar plots
#   6. Runs diversity analysis
#
# Prerequisites:
#   - QIIME2 environment activated
#   - Input files in expected locations (see Configuration)
#
# Input Directory Structure:
#   your-data/
#   ├── biom/
#   │   ├── {MARKER}-asv.biom
#   │   └── {MARKER}-taxa.biom
#   ├── fasta/
#   │   ├── {MARKER}_paired_F.fasta
#   │   └── {MARKER}_paired_R.fasta
#   └── metadata/
#       └── sample_metadata.tsv
#
# Output:
#   qiime2_output/
#   ├── feature-table.qza
#   ├── rep-seqs-forward.qza
#   ├── rep-seqs-reverse.qza
#   ├── taxonomy.qza
#   ├── taxa-barplot.qzv
#   ├── diversity/
#   │   ├── core-metrics/
#   │   └── *.qzv
#   └── *.qzv (various visualizations)
#
# Usage:
#   ./complete_workflow.sh
#   ./complete_workflow.sh 16S_Bacteria  # Specify marker name
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Marker name (can be overridden via command line argument)
MARKER="${1:-16S_Bacteria}"

# Input files (derived from marker name)
ASV_BIOM="biom/${MARKER}-asv.biom"
TAXA_BIOM="biom/${MARKER}-taxa.biom"
FORWARD_FASTA="fasta/${MARKER}_paired_F.fasta"
REVERSE_FASTA="fasta/${MARKER}_paired_R.fasta"
METADATA_FILE="metadata/sample_metadata.tsv"

# Output directory
OUTPUT_DIR="qiime2_output"

# Diversity analysis parameters
SAMPLING_DEPTH="auto"  # Set to a number to override
GROUP_COLUMN="country" # Metadata column for comparisons

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

echo "============================================================"
echo "  QIIME2 Complete Workflow for eDNA Explorer"
echo "============================================================"
echo ""
echo "Marker: $MARKER"
echo ""

# Check QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "Error: QIIME2 not found. Please activate your QIIME2 environment:"
    echo "  conda activate qiime2-2024.10"
    exit 1
fi

echo "QIIME2 version: $(qiime --version | head -1)"
echo ""

# Check input files
echo "Checking input files..."
MISSING_FILES=0

for FILE in "$ASV_BIOM" "$FORWARD_FASTA" "$METADATA_FILE"; do
    if [[ -f "$FILE" ]]; then
        echo "  ✓ $FILE"
    else
        echo "  ✗ $FILE (not found)"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
done

# Optional files
for FILE in "$TAXA_BIOM" "$REVERSE_FASTA"; do
    if [[ -f "$FILE" ]]; then
        echo "  ✓ $FILE"
    else
        echo "  ○ $FILE (optional, not found)"
    fi
done

if [[ $MISSING_FILES -gt 0 ]]; then
    echo ""
    echo "Error: $MISSING_FILES required file(s) not found."
    echo "Please check your file paths and try again."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""

# -----------------------------------------------------------------------------
# Step 1: Import Feature Table
# -----------------------------------------------------------------------------

echo "============================================================"
echo "  Step 1: Import Feature Table"
echo "============================================================"
echo ""

qiime tools import \
    --input-path "$ASV_BIOM" \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path "${OUTPUT_DIR}/feature-table.qza"

echo "✓ Feature table imported"

# Generate summary
qiime feature-table summarize \
    --i-table "${OUTPUT_DIR}/feature-table.qza" \
    --m-sample-metadata-file "$METADATA_FILE" \
    --o-visualization "${OUTPUT_DIR}/feature-table.qzv"

echo "✓ Feature table summary generated"
echo ""

# -----------------------------------------------------------------------------
# Step 2: Import Forward Sequences
# -----------------------------------------------------------------------------

echo "============================================================"
echo "  Step 2: Import Forward Sequences"
echo "============================================================"
echo ""

qiime tools import \
    --input-path "$FORWARD_FASTA" \
    --output-path "${OUTPUT_DIR}/rep-seqs-forward.qza" \
    --type 'FeatureData[Sequence]'

echo "✓ Forward sequences imported"

qiime feature-table tabulate-seqs \
    --i-data "${OUTPUT_DIR}/rep-seqs-forward.qza" \
    --o-visualization "${OUTPUT_DIR}/rep-seqs-forward.qzv"

echo "✓ Forward sequence summary generated"
echo ""

# -----------------------------------------------------------------------------
# Step 3: Import Reverse Sequences (Optional)
# -----------------------------------------------------------------------------

echo "============================================================"
echo "  Step 3: Import Reverse Sequences"
echo "============================================================"
echo ""

if [[ -f "$REVERSE_FASTA" ]]; then
    echo "Preprocessing reverse FASTA (transforming _paired_R_ → _paired_F_)..."

    # Transform headers
    TEMP_FASTA="${OUTPUT_DIR}/.reverse_temp.fasta"
    awk '{
        if (/^>/) {
            gsub(/_paired_R_/, "_paired_F_")
        }
        print
    }' "$REVERSE_FASTA" > "$TEMP_FASTA"

    qiime tools import \
        --input-path "$TEMP_FASTA" \
        --output-path "${OUTPUT_DIR}/rep-seqs-reverse.qza" \
        --type 'FeatureData[Sequence]'

    rm -f "$TEMP_FASTA"

    echo "✓ Reverse sequences imported"

    qiime feature-table tabulate-seqs \
        --i-data "${OUTPUT_DIR}/rep-seqs-reverse.qza" \
        --o-visualization "${OUTPUT_DIR}/rep-seqs-reverse.qzv"

    echo "✓ Reverse sequence summary generated"
else
    echo "○ Skipping (reverse FASTA not found)"
fi
echo ""

# -----------------------------------------------------------------------------
# Step 4: Import Taxonomy
# -----------------------------------------------------------------------------

echo "============================================================"
echo "  Step 4: Import Taxonomy"
echo "============================================================"
echo ""

# Use taxa BIOM if available, otherwise extract from ASV BIOM
TAXONOMY_SOURCE="$TAXA_BIOM"
if [[ ! -f "$TAXONOMY_SOURCE" ]]; then
    TAXONOMY_SOURCE="$ASV_BIOM"
fi

echo "Extracting taxonomy from: $TAXONOMY_SOURCE"

# Extract taxonomy using Python
TAXONOMY_TSV="${OUTPUT_DIR}/taxonomy.tsv"

python3 << EOF
import biom

table = biom.load_table("$TAXONOMY_SOURCE")

with open("$TAXONOMY_TSV", 'w') as f:
    f.write("Feature ID\tTaxon\n")

    for feature_id in table.ids(axis='observation'):
        metadata = table.metadata(feature_id, axis='observation')

        if metadata is None:
            taxonomy = "Unassigned"
        elif 'taxonomy' in metadata:
            tax_list = metadata['taxonomy']
            if isinstance(tax_list, list):
                taxonomy = "; ".join(str(t) for t in tax_list if t)
            else:
                taxonomy = str(tax_list)
        elif 'taxonomic_path' in metadata:
            taxonomy = str(metadata['taxonomic_path'])
        else:
            taxonomy = "Unassigned"

        f.write(f"{feature_id}\t{taxonomy}\n")

print(f"Extracted taxonomy for {len(table.ids(axis='observation'))} features")
EOF

qiime tools import \
    --input-path "$TAXONOMY_TSV" \
    --type 'FeatureData[Taxonomy]' \
    --input-format TSVTaxonomyFormat \
    --output-path "${OUTPUT_DIR}/taxonomy.qza"

echo "✓ Taxonomy imported"
echo ""

# -----------------------------------------------------------------------------
# Step 5: Taxonomic Bar Plot
# -----------------------------------------------------------------------------

echo "============================================================"
echo "  Step 5: Generate Taxonomic Bar Plot"
echo "============================================================"
echo ""

qiime taxa barplot \
    --i-table "${OUTPUT_DIR}/feature-table.qza" \
    --i-taxonomy "${OUTPUT_DIR}/taxonomy.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --o-visualization "${OUTPUT_DIR}/taxa-barplot.qzv"

echo "✓ Taxonomic bar plot generated"
echo ""

# -----------------------------------------------------------------------------
# Step 6: Diversity Analysis
# -----------------------------------------------------------------------------

echo "============================================================"
echo "  Step 6: Diversity Analysis"
echo "============================================================"
echo ""

DIVERSITY_DIR="${OUTPUT_DIR}/diversity"
mkdir -p "$DIVERSITY_DIR"

# Determine sampling depth
if [[ "$SAMPLING_DEPTH" == "auto" ]]; then
    echo "Determining sampling depth..."

    TEMP_DIR=$(mktemp -d)
    qiime tools export \
        --input-path "${OUTPUT_DIR}/feature-table.qza" \
        --output-path "$TEMP_DIR" 2>/dev/null

    SAMPLING_DEPTH=$(python3 << EOF
import biom
import numpy as np
table = biom.load_table("$TEMP_DIR/feature-table.biom")
depths = table.sum(axis='sample')
# Use 10th percentile as sampling depth
print(int(np.percentile(depths, 10)))
EOF
)

    rm -rf "$TEMP_DIR"
    echo "  Using sampling depth: $SAMPLING_DEPTH"
fi

# Alpha rarefaction
echo ""
echo "Running alpha rarefaction..."
qiime diversity alpha-rarefaction \
    --i-table "${OUTPUT_DIR}/feature-table.qza" \
    --p-max-depth 10000 \
    --m-metadata-file "$METADATA_FILE" \
    --o-visualization "${DIVERSITY_DIR}/alpha-rarefaction.qzv"

echo "✓ Alpha rarefaction complete"

# Core diversity metrics
echo ""
echo "Computing core diversity metrics..."
qiime diversity core-metrics \
    --i-table "${OUTPUT_DIR}/feature-table.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --output-dir "${DIVERSITY_DIR}/core-metrics"

echo "✓ Core metrics computed"

# Alpha significance
echo ""
echo "Testing alpha diversity significance..."
qiime diversity alpha-group-significance \
    --i-alpha-diversity "${DIVERSITY_DIR}/core-metrics/shannon_vector.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --o-visualization "${DIVERSITY_DIR}/shannon-significance.qzv" \
    2>/dev/null || echo "  Warning: Could not compute Shannon significance"

# Beta significance (PERMANOVA)
echo ""
echo "Running PERMANOVA on Bray-Curtis distances..."
qiime diversity beta-group-significance \
    --i-distance-matrix "${DIVERSITY_DIR}/core-metrics/bray_curtis_distance_matrix.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --m-metadata-column "$GROUP_COLUMN" \
    --p-method permanova \
    --p-permutations 999 \
    --o-visualization "${DIVERSITY_DIR}/bray-curtis-permanova.qzv" \
    2>/dev/null || echo "  Warning: Could not compute PERMANOVA"

echo "✓ Diversity analysis complete"
echo ""

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

echo "============================================================"
echo "  Workflow Complete!"
echo "============================================================"
echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Generated artifacts:"
echo ""
echo "  Feature Table:"
echo "    ${OUTPUT_DIR}/feature-table.qza"
echo "    ${OUTPUT_DIR}/feature-table.qzv"
echo ""
echo "  Sequences:"
echo "    ${OUTPUT_DIR}/rep-seqs-forward.qza"
if [[ -f "${OUTPUT_DIR}/rep-seqs-reverse.qza" ]]; then
echo "    ${OUTPUT_DIR}/rep-seqs-reverse.qza"
fi
echo ""
echo "  Taxonomy:"
echo "    ${OUTPUT_DIR}/taxonomy.qza"
echo "    ${OUTPUT_DIR}/taxa-barplot.qzv"
echo ""
echo "  Diversity:"
echo "    ${DIVERSITY_DIR}/alpha-rarefaction.qzv"
echo "    ${DIVERSITY_DIR}/core-metrics/"
echo "    ${DIVERSITY_DIR}/shannon-significance.qzv"
echo "    ${DIVERSITY_DIR}/bray-curtis-permanova.qzv"
echo ""
echo "View visualizations:"
echo "  qiime tools view ${OUTPUT_DIR}/*.qzv"
echo "  qiime tools view ${DIVERSITY_DIR}/*.qzv"
echo ""
