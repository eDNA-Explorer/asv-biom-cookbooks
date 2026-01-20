# =============================================================================
# 01_basic_import.R
# Basic BIOM Import into Phyloseq
# =============================================================================
#
# Description:
#   Import an eDNA Explorer BIOM file into Phyloseq. This loads the feature
#   table (OTU/ASV counts), sample metadata, and taxonomy (if embedded).
#
# Prerequisites:
#   - phyloseq package installed
#
# Input:
#   - BIOM file (e.g., "biom/16S_Bacteria-asv.biom")
#
# Output:
#   - phyloseq object with otu_table, sample_data, and tax_table
#
# =============================================================================

library(phyloseq)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Path to your BIOM file
# Adjust this to match your file location
BIOM_FILE <- "biom/16S_Bacteria-asv.biom"

# -----------------------------------------------------------------------------
# Import BIOM File
# -----------------------------------------------------------------------------

# Import the BIOM file
# This automatically loads:
#   - Feature table (ASV/OTU counts per sample)
#   - Sample metadata (if embedded)
#   - Taxonomy table (if embedded)
ps <- import_biom(BIOM_FILE)

# -----------------------------------------------------------------------------
# Inspect the Result
# -----------------------------------------------------------------------------

# Print summary
print(ps)
# Example output:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5000 taxa and 200 samples ]
# sample_data() Sample Data:       [ 200 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 5000 taxa by 7 taxonomic ranks ]

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
