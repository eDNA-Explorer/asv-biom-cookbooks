#!/bin/bash
# =============================================================================
# 05_diversity_analysis.sh
# Run Diversity Analysis on eDNA Explorer Data in QIIME2
# =============================================================================
#
# Description:
#   Perform alpha and beta diversity analysis on eDNA Explorer data that has
#   been imported into QIIME2. This includes:
#   - Alpha rarefaction curves
#   - Core diversity metrics (Shannon, observed features, Bray-Curtis, Jaccard)
#   - Alpha diversity significance testing
#   - Beta diversity group significance (PERMANOVA)
#
# Prerequisites:
#   - QIIME2 environment activated
#   - Feature table imported (feature-table.qza)
#   - Sample metadata file (TSV format)
#
# Input:
#   - Feature table artifact (qiime2_output/feature-table.qza)
#   - Sample metadata (metadata/sample_metadata.tsv)
#
# Output:
#   - Alpha rarefaction visualization
#   - Core diversity metrics directory
#   - Significance test results
#
# Usage:
#   ./05_diversity_analysis.sh
#   ./05_diversity_analysis.sh 5000  # Specify sampling depth
#
# =============================================================================

set -e  # Exit on error

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Input files
FEATURE_TABLE="qiime2_output/feature-table.qza"
METADATA_FILE="metadata/sample_metadata.tsv"

# Output directory
OUTPUT_DIR="qiime2_output/diversity"

# Sampling depth for rarefaction (can be overridden via command line)
# Set to "auto" to determine from data
SAMPLING_DEPTH="${1:-auto}"

# Metadata column for group comparisons
GROUP_COLUMN="country"  # Change to match your metadata

# Maximum depth for rarefaction curves
MAX_DEPTH="${2:-10000}"

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

echo "=== QIIME2 Diversity Analysis ==="
echo ""

# Check QIIME2 is available
if ! command -v qiime &> /dev/null; then
    echo "Error: QIIME2 not found. Please activate your QIIME2 environment:"
    echo "  conda activate qiime2-2024.10"
    exit 1
fi

# Check feature table exists
if [[ ! -f "$FEATURE_TABLE" ]]; then
    echo "Error: Feature table not found: $FEATURE_TABLE"
    echo "Run 01_import_feature_table.sh first."
    exit 1
fi

# Check metadata exists
if [[ ! -f "$METADATA_FILE" ]]; then
    echo "Error: Metadata file not found: $METADATA_FILE"
    echo "Please provide a sample metadata TSV file."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# Determine Sampling Depth
# -----------------------------------------------------------------------------

echo "Feature table: $FEATURE_TABLE"
echo "Metadata file: $METADATA_FILE"
echo ""

if [[ "$SAMPLING_DEPTH" == "auto" ]]; then
    echo "Determining appropriate sampling depth..."

    # Export feature table to get sample depths
    TEMP_DIR=$(mktemp -d)
    qiime tools export \
        --input-path "$FEATURE_TABLE" \
        --output-path "$TEMP_DIR" 2>/dev/null

    # Get sample depths using biom
    DEPTHS=$(python3 << EOF
import biom
import numpy as np
table = biom.load_table("$TEMP_DIR/feature-table.biom")
depths = table.sum(axis='sample')
print(f"min:{int(min(depths))}")
print(f"median:{int(np.median(depths))}")
print(f"recommended:{int(np.percentile(depths, 10))}")
EOF
)

    MIN_DEPTH=$(echo "$DEPTHS" | grep "min:" | cut -d: -f2)
    MEDIAN_DEPTH=$(echo "$DEPTHS" | grep "median:" | cut -d: -f2)
    RECOMMENDED=$(echo "$DEPTHS" | grep "recommended:" | cut -d: -f2)

    rm -rf "$TEMP_DIR"

    echo "  Minimum sample depth: $MIN_DEPTH"
    echo "  Median sample depth: $MEDIAN_DEPTH"
    echo "  Recommended depth (10th percentile): $RECOMMENDED"

    # Use 10th percentile as sampling depth
    SAMPLING_DEPTH=$RECOMMENDED

    if [[ $SAMPLING_DEPTH -lt 1000 ]]; then
        echo ""
        echo "Warning: Recommended depth ($SAMPLING_DEPTH) is very low."
        echo "Consider filtering low-depth samples first."
    fi
fi

echo ""
echo "Using sampling depth: $SAMPLING_DEPTH"
echo ""

# -----------------------------------------------------------------------------
# Alpha Rarefaction
# -----------------------------------------------------------------------------

echo "=== Step 1: Alpha Rarefaction Curves ==="
echo ""

qiime diversity alpha-rarefaction \
    --i-table "$FEATURE_TABLE" \
    --p-max-depth "$MAX_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --o-visualization "${OUTPUT_DIR}/alpha-rarefaction.qzv"

echo "✓ Alpha rarefaction: ${OUTPUT_DIR}/alpha-rarefaction.qzv"

# -----------------------------------------------------------------------------
# Core Diversity Metrics
# -----------------------------------------------------------------------------

echo ""
echo "=== Step 2: Core Diversity Metrics ==="
echo ""

# Core metrics computes:
# - Alpha: Shannon, observed features, evenness, Faith's PD (if tree provided)
# - Beta: Bray-Curtis, Jaccard, unweighted/weighted UniFrac (if tree provided)

qiime diversity core-metrics \
    --i-table "$FEATURE_TABLE" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA_FILE" \
    --output-dir "${OUTPUT_DIR}/core-metrics"

echo "✓ Core metrics computed"
echo ""
echo "Generated artifacts:"
ls -1 "${OUTPUT_DIR}/core-metrics/"

# -----------------------------------------------------------------------------
# Alpha Diversity Significance
# -----------------------------------------------------------------------------

echo ""
echo "=== Step 3: Alpha Diversity Significance ==="
echo ""

# Shannon diversity
echo "Testing Shannon diversity by ${GROUP_COLUMN}..."
qiime diversity alpha-group-significance \
    --i-alpha-diversity "${OUTPUT_DIR}/core-metrics/shannon_vector.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --o-visualization "${OUTPUT_DIR}/shannon-significance.qzv" \
    2>/dev/null || echo "  Warning: Could not compute Shannon significance"

# Observed features
echo "Testing observed features by ${GROUP_COLUMN}..."
qiime diversity alpha-group-significance \
    --i-alpha-diversity "${OUTPUT_DIR}/core-metrics/observed_features_vector.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --o-visualization "${OUTPUT_DIR}/observed-features-significance.qzv" \
    2>/dev/null || echo "  Warning: Could not compute observed features significance"

echo "✓ Alpha significance tests complete"

# -----------------------------------------------------------------------------
# Beta Diversity Significance (PERMANOVA)
# -----------------------------------------------------------------------------

echo ""
echo "=== Step 4: Beta Diversity Significance (PERMANOVA) ==="
echo ""

# Bray-Curtis PERMANOVA
echo "Testing Bray-Curtis by ${GROUP_COLUMN}..."
qiime diversity beta-group-significance \
    --i-distance-matrix "${OUTPUT_DIR}/core-metrics/bray_curtis_distance_matrix.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --m-metadata-column "$GROUP_COLUMN" \
    --p-method permanova \
    --p-permutations 999 \
    --o-visualization "${OUTPUT_DIR}/bray-curtis-permanova.qzv" \
    2>/dev/null || echo "  Warning: Could not compute Bray-Curtis PERMANOVA"

# Jaccard PERMANOVA
echo "Testing Jaccard by ${GROUP_COLUMN}..."
qiime diversity beta-group-significance \
    --i-distance-matrix "${OUTPUT_DIR}/core-metrics/jaccard_distance_matrix.qza" \
    --m-metadata-file "$METADATA_FILE" \
    --m-metadata-column "$GROUP_COLUMN" \
    --p-method permanova \
    --p-permutations 999 \
    --o-visualization "${OUTPUT_DIR}/jaccard-permanova.qzv" \
    2>/dev/null || echo "  Warning: Could not compute Jaccard PERMANOVA"

echo "✓ Beta significance tests complete"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

echo ""
echo "=== Diversity Analysis Complete ==="
echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Generated visualizations:"
echo "  Alpha rarefaction:      ${OUTPUT_DIR}/alpha-rarefaction.qzv"
echo "  Shannon significance:   ${OUTPUT_DIR}/shannon-significance.qzv"
echo "  Observed features sig:  ${OUTPUT_DIR}/observed-features-significance.qzv"
echo "  Bray-Curtis PERMANOVA:  ${OUTPUT_DIR}/bray-curtis-permanova.qzv"
echo "  Jaccard PERMANOVA:      ${OUTPUT_DIR}/jaccard-permanova.qzv"
echo ""
echo "Core metrics artifacts:"
echo "  ${OUTPUT_DIR}/core-metrics/"
echo ""
echo "View visualizations with:"
echo "  qiime tools view ${OUTPUT_DIR}/*.qzv"
