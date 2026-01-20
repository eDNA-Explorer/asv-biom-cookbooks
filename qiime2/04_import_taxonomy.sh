#!/bin/bash
# =============================================================================
# 04_import_taxonomy.sh
# Import Taxonomy from eDNA Explorer BIOM into QIIME2
# =============================================================================
#
# Description:
#   Extract and import taxonomy from an eDNA Explorer BIOM file into QIIME2.
#   This creates a FeatureData[Taxonomy] artifact that can be used for
#   taxonomic analysis and visualization.
#
# Prerequisites:
#   - QIIME2 environment activated
#   - Python 3 with biom-format package (usually included with QIIME2)
#
# Input:
#   - BIOM file with embedded taxonomy (e.g., "biom/16S_Bacteria-taxa.biom")
#
# Output:
#   - QIIME2 artifact: taxonomy.qza
#   - Visualization: taxonomy.qzv
#   - TSV file: taxonomy.tsv (intermediate)
#
# Usage:
#   ./04_import_taxonomy.sh
#   ./04_import_taxonomy.sh biom/my-marker-taxa.biom
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Input BIOM file (can be overridden via command line argument)
# The taxa BIOM file has taxonomy in an easier-to-extract format
BIOM_FILE="${1:-biom/16S_Bacteria-taxa.biom}"

# Output directory
OUTPUT_DIR="qiime2_output"

# Output file names
TAXONOMY_TSV="${OUTPUT_DIR}/taxonomy.tsv"
OUTPUT_QZA="${OUTPUT_DIR}/taxonomy.qza"
OUTPUT_QZV="${OUTPUT_DIR}/taxonomy.qzv"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

echo "=== QIIME2 Taxonomy Import ==="
echo ""

# Check QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "Error: QIIME2 not found. Please activate your QIIME2 environment:"
    echo "  conda activate qiime2-2024.10"
    exit 1
fi

# Check input file exists
if [[ ! -f "$BIOM_FILE" ]]; then
    echo "Error: BIOM file not found: $BIOM_FILE"
    echo ""
    echo "Usage: $0 [biom_file]"
    echo "Example: $0 biom/16S_Bacteria-taxa.biom"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# Extract Taxonomy from BIOM
# -----------------------------------------------------------------------------

echo "Input:  $BIOM_FILE"
echo "Output: $OUTPUT_QZA"
echo ""

echo "Step 1: Extracting taxonomy from BIOM file..."

# Use Python to extract taxonomy from BIOM file
python3 << EOF
import sys
try:
    import biom
except ImportError:
    print("Error: biom-format package not found.")
    print("Install with: pip install biom-format")
    sys.exit(1)

# Load BIOM file
print(f"  Loading: $BIOM_FILE")
table = biom.load_table("$BIOM_FILE")

# Check for observation metadata
if table.metadata(axis='observation') is None:
    print("Error: No observation metadata found in BIOM file.")
    print("This file may not contain taxonomy information.")
    sys.exit(1)

# Extract taxonomy
print(f"  Found {len(table.ids(axis='observation'))} features")

# Write TSV file
with open("$TAXONOMY_TSV", 'w') as f:
    f.write("Feature ID\tTaxon\n")

    for feature_id in table.ids(axis='observation'):
        metadata = table.metadata(feature_id, axis='observation')

        if metadata is None:
            taxonomy = "Unassigned"
        elif 'taxonomy' in metadata:
            # Taxonomy as list: ["k__Bacteria", "p__Proteobacteria", ...]
            tax_list = metadata['taxonomy']
            if isinstance(tax_list, list):
                taxonomy = "; ".join(str(t) for t in tax_list if t)
            else:
                taxonomy = str(tax_list)
        elif 'taxonomic_path' in metadata:
            # Taxonomy as string
            taxonomy = str(metadata['taxonomic_path'])
        else:
            taxonomy = "Unassigned"

        f.write(f"{feature_id}\t{taxonomy}\n")

print(f"  Taxonomy extracted to: $TAXONOMY_TSV")
EOF

if [[ ! -f "$TAXONOMY_TSV" ]]; then
    echo "Error: Failed to extract taxonomy"
    exit 1
fi

# Count taxa
TAXA_COUNT=$(tail -n +2 "$TAXONOMY_TSV" | wc -l | tr -d ' ')
echo "✓ Extracted taxonomy for $TAXA_COUNT features"

# Show sample
echo ""
echo "Sample taxonomy entries:"
head -5 "$TAXONOMY_TSV" | column -t -s $'\t'

# -----------------------------------------------------------------------------
# Import Taxonomy into QIIME2
# -----------------------------------------------------------------------------

echo ""
echo "Step 2: Importing taxonomy into QIIME2..."

qiime tools import \
    --input-path "$TAXONOMY_TSV" \
    --type 'FeatureData[Taxonomy]' \
    --input-format TSVTaxonomyFormat \
    --output-path "$OUTPUT_QZA"

echo "✓ Taxonomy imported: $OUTPUT_QZA"

# -----------------------------------------------------------------------------
# Generate Taxonomy Visualization
# -----------------------------------------------------------------------------

echo ""
echo "Step 3: Generating taxonomy visualization..."

qiime metadata tabulate \
    --m-input-file "$OUTPUT_QZA" \
    --o-visualization "$OUTPUT_QZV"

echo "✓ Taxonomy visualization: $OUTPUT_QZV"

# -----------------------------------------------------------------------------
# Display Summary Info
# -----------------------------------------------------------------------------

echo ""
echo "=== Import Complete ==="
echo ""
echo "Artifact info:"
qiime tools peek "$OUTPUT_QZA"

echo ""
echo "View taxonomy with:"
echo "  qiime tools view $OUTPUT_QZV"
echo ""
echo "Use in analysis:"
echo "  qiime taxa barplot \\"
echo "      --i-table feature-table.qza \\"
echo "      --i-taxonomy $OUTPUT_QZA \\"
echo "      --m-metadata-file metadata/sample_metadata.tsv \\"
echo "      --o-visualization taxa-barplot.qzv"
