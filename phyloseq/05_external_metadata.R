# =============================================================================
# 05_external_metadata.R
# Import BIOM and Merge with External Sample Metadata
# =============================================================================
#
# Description:
#   Import an eDNA Explorer BIOM file and merge it with your own external
#   sample metadata. This is useful when you have additional variables
#   (treatments, timepoints, environmental data) not included in the
#   original eDNA Explorer export.
#
# Prerequisites:
#   - phyloseq package installed
#
# Input:
#   - BIOM file (e.g., "biom/16S_Bacteria-asv.biom")
#   - External metadata file (CSV or TSV)
#
# Output:
#   - phyloseq object with merged sample metadata
#
# =============================================================================

library(phyloseq)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Paths to your files
BIOM_FILE <- "biom/16S_Bacteria-asv.biom"

# Your external metadata file
# Should have a column that matches sample IDs in the BIOM file
METADATA_FILE <- "metadata/my_sample_metadata.csv"
# METADATA_FILE <- "metadata/my_sample_metadata.tsv"  # For TSV files

# Column in your metadata that contains sample IDs
# This should match either 'sample_id' (internal CUID) or 'sample_name'
# from the BIOM file
SAMPLE_ID_COLUMN <- "sample_name"

# -----------------------------------------------------------------------------
# Step 1: Import BIOM File
# -----------------------------------------------------------------------------

cat("=== Step 1: Importing BIOM File ===\n\n")

ps <- import_biom(BIOM_FILE)

cat("BIOM imported:\n")
cat("  Samples:", nsamples(ps), "\n")
cat("  Taxa:", ntaxa(ps), "\n")

# View existing sample metadata
cat("\nExisting sample metadata columns:\n")
print(sample_variables(ps))

cat("\nFirst few samples:\n")
print(head(sample_data(ps)[, c("sample_id", "sample_name")]))

# -----------------------------------------------------------------------------
# Step 2: Load External Metadata
# -----------------------------------------------------------------------------

cat("\n=== Step 2: Loading External Metadata ===\n\n")

# Detect file type and read accordingly
if (grepl("\\.csv$", METADATA_FILE, ignore.case = TRUE)) {
    external_meta <- read.csv(METADATA_FILE, stringsAsFactors = FALSE)
} else if (grepl("\\.tsv$|\\.txt$", METADATA_FILE, ignore.case = TRUE)) {
    external_meta <- read.delim(METADATA_FILE, stringsAsFactors = FALSE)
} else {
    stop("Unrecognized file format. Use .csv or .tsv")
}

cat("External metadata loaded:\n")
cat("  Rows:", nrow(external_meta), "\n")
cat("  Columns:", ncol(external_meta), "\n")
cat("\nColumn names:\n")
print(names(external_meta))

cat("\nFirst few rows:\n")
print(head(external_meta))

# -----------------------------------------------------------------------------
# Step 3: Match Sample IDs
# -----------------------------------------------------------------------------

cat("\n=== Step 3: Matching Sample IDs ===\n\n")

# Get sample IDs from BIOM
biom_sample_ids <- sample_names(ps)

# Get the matching column from your metadata
if (SAMPLE_ID_COLUMN == "sample_name") {
    # Match by sample_name (user-provided names)
    biom_match_ids <- sample_data(ps)$sample_name
} else if (SAMPLE_ID_COLUMN == "sample_id") {
    # Match by sample_id (internal CUIDs)
    biom_match_ids <- sample_data(ps)$sample_id
} else {
    # Use sample_names from phyloseq directly
    biom_match_ids <- biom_sample_ids
}

# Check for matching IDs
external_ids <- external_meta[[SAMPLE_ID_COLUMN]]
matching_ids <- intersect(biom_match_ids, external_ids)

cat("BIOM samples:", length(biom_match_ids), "\n")
cat("External metadata rows:", length(external_ids), "\n")
cat("Matching samples:", length(matching_ids), "\n")

if (length(matching_ids) == 0) {
    stop("No matching sample IDs found. Check your SAMPLE_ID_COLUMN setting.")
}

if (length(matching_ids) < length(biom_match_ids)) {
    missing <- setdiff(biom_match_ids, external_ids)
    cat("\nWarning:", length(missing), "BIOM samples not in external metadata:\n")
    print(head(missing, 10))
}

# -----------------------------------------------------------------------------
# Step 4: Merge Metadata
# -----------------------------------------------------------------------------

cat("\n=== Step 4: Merging Metadata ===\n\n")

# Set row names on external metadata to match sample IDs
rownames(external_meta) <- external_meta[[SAMPLE_ID_COLUMN]]

# If matching by sample_name, we need to align with phyloseq's sample_names
if (SAMPLE_ID_COLUMN == "sample_name") {
    # Create a mapping from sample_name to phyloseq sample_names
    mapping <- data.frame(
        ps_name = sample_names(ps),
        sample_name = sample_data(ps)$sample_name,
        stringsAsFactors = FALSE
    )

    # Reindex external metadata to use phyloseq sample names
    external_meta_aligned <- external_meta[mapping$sample_name, , drop = FALSE]
    rownames(external_meta_aligned) <- mapping$ps_name
    external_meta <- external_meta_aligned
}

# Convert to sample_data object
external_sample_data <- sample_data(external_meta)

# Option A: Replace all sample metadata with external
# ps_merged <- merge_phyloseq(
#     otu_table(ps),
#     tax_table(ps),
#     external_sample_data
# )

# Option B: Merge external columns into existing metadata (recommended)
# Keep existing metadata and add new columns
existing_meta <- as.data.frame(sample_data(ps))
new_columns <- setdiff(names(external_meta), names(existing_meta))

if (length(new_columns) > 0) {
    cat("Adding new columns:", paste(new_columns, collapse = ", "), "\n")

    for (col in new_columns) {
        existing_meta[[col]] <- external_meta[rownames(existing_meta), col]
    }

    sample_data(ps) <- sample_data(existing_meta)
}

# -----------------------------------------------------------------------------
# Verify Merge
# -----------------------------------------------------------------------------

cat("\n=== Verification ===\n\n")

cat("Updated sample metadata columns:\n")
print(sample_variables(ps))

cat("\nFirst few rows of merged metadata:\n")
print(head(sample_data(ps)))

# -----------------------------------------------------------------------------
# Example: Using New Metadata Variables
# -----------------------------------------------------------------------------

cat("\n=== Example Usage ===\n\n")

# If you added a 'treatment' column, you can now use it:
# plot_richness(ps, x = "treatment", measures = "Shannon")
# ordinate and color by treatment, etc.

cat("Your phyloseq object now has", length(sample_variables(ps)), "sample variables.\n")
cat("Use them in analyses like:\n")
cat("  plot_richness(ps, x = 'your_variable')\n")
cat("  plot_ordination(ps, ord, color = 'your_variable')\n")

# -----------------------------------------------------------------------------
# Save for Later Use (Optional)
# -----------------------------------------------------------------------------

# saveRDS(ps, "phyloseq_with_metadata.rds")

cat("\n✓ External metadata successfully merged\n")
