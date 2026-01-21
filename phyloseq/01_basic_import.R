# =============================================================================
# 01_basic_import.R
# Basic BIOM Import into Phyloseq
# =============================================================================
#
# Description:
#   Import an eDNA Explorer v2.1 BIOM file into Phyloseq. v2.1 ASV BIOM files
#   contain only the count matrix - metadata is in sidecar files:
#   - samples.tsv: Sample metadata (350+ environmental variables)
#   - {marker}-{mode}-lookup.tsv: ASV taxonomy lookup
#
# Prerequisites:
#   - phyloseq package installed
#   - Extract bundle first: zstd -d bundle.tar.zst -c | tar -xf -
#
# Input:
#   - BIOM file (e.g., "16S_Bacteria-paired-asv.biom")
#   - samples.tsv (sample metadata)
#   - lookup.tsv (optional, for taxonomy)
#
# Output:
#   - phyloseq object with otu_table, sample_data, and optionally tax_table
#
# =============================================================================

library(phyloseq)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Path to your BIOM file (v2.1 flat directory structure)
# Adjust this to match your file location
BIOM_FILE <- "16S_Bacteria-paired-asv.biom"
SAMPLES_FILE <- "samples.tsv"
LOOKUP_FILE <- "16S_Bacteria-paired-lookup.tsv"  # For taxonomy

# -----------------------------------------------------------------------------
# Import BIOM File
# -----------------------------------------------------------------------------

# Import the BIOM file (v2.1: counts only, no embedded metadata)
ps <- import_biom(BIOM_FILE)

# -----------------------------------------------------------------------------
# Load Sample Metadata from Sidecar File
# -----------------------------------------------------------------------------

if (file.exists(SAMPLES_FILE)) {
    cat("Loading sample metadata from:", SAMPLES_FILE, "\n")
    sample_meta <- read.delim(SAMPLES_FILE, row.names = 1, check.names = FALSE)

    # Filter to samples in BIOM
    common_samples <- intersect(rownames(sample_meta), sample_names(ps))
    if (length(common_samples) > 0) {
        sample_data(ps) <- sample_data(sample_meta[common_samples, , drop = FALSE])
        cat("Loaded metadata for", length(common_samples), "samples\n")
    }
} else {
    cat("Note: samples.tsv not found. Sample metadata not loaded.\n")
}

# -----------------------------------------------------------------------------
# Inspect the Result
# -----------------------------------------------------------------------------

# Print summary
print(ps)
# Example output for v2.1 bundle:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1322726 taxa and 14 samples ]
# sample_data() Sample Data:       [ 14 samples by 349 sample variables ]

# Check dimensions
cat("\nNumber of taxa (ASVs):", ntaxa(ps), "\n")
cat("Number of samples:", nsamples(ps), "\n")

# -----------------------------------------------------------------------------
# Explore Components
# -----------------------------------------------------------------------------

# View first few rows of the OTU table
cat("\n--- OTU Table (first 5 taxa, first 5 samples) ---\n")
print(otu_table(ps)[1:5, 1:5])

# View sample metadata
cat("\n--- Sample Metadata (first 5 samples) ---\n")
print(head(sample_data(ps), 5))

# View available sample variables
cat("\n--- Sample Variables ---\n")
print(sample_variables(ps))

# View taxonomy table (if present)
if (!is.null(tax_table(ps, errorIfNULL = FALSE))) {
    cat("\n--- Taxonomy Table (first 5 taxa) ---\n")
    print(head(tax_table(ps), 5))

    # View taxonomy ranks
    cat("\n--- Taxonomy Ranks ---\n")
    print(rank_names(ps))
}

# -----------------------------------------------------------------------------
# Basic Statistics
# -----------------------------------------------------------------------------

# Total reads per sample
sample_sums_vec <- sample_sums(ps)
cat("\n--- Reads per Sample ---\n")
cat("Min:", min(sample_sums_vec), "\n")
cat("Max:", max(sample_sums_vec), "\n")
cat("Mean:", round(mean(sample_sums_vec), 1), "\n")
cat("Median:", median(sample_sums_vec), "\n")

# Total reads per taxon
taxa_sums_vec <- taxa_sums(ps)
cat("\n--- Reads per Taxon ---\n")
cat("Min:", min(taxa_sums_vec), "\n")
cat("Max:", max(taxa_sums_vec), "\n")
cat("Mean:", round(mean(taxa_sums_vec), 1), "\n")

# -----------------------------------------------------------------------------
# Save for Later Use (Optional)
# -----------------------------------------------------------------------------

# Save the phyloseq object for later sessions
# saveRDS(ps, "phyloseq_object.rds")

# Load it back with:
# ps <- readRDS("phyloseq_object.rds")

cat("\n✓ BIOM file successfully imported into Phyloseq\n")
