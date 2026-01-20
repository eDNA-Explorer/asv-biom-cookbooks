# =============================================================================
# 06_sample_id_mapping.R
# Map Internal Sample IDs to User-Friendly Names
# =============================================================================
#
# Description:
#   eDNA Explorer uses internal CUIDs (e.g., "clxyz123abc") as sample
#   identifiers. This script shows how to work with the user-provided
#   sample names for display and analysis.
#
# Prerequisites:
#   - phyloseq package installed
#
# Input:
#   - BIOM file with embedded sample metadata
#
# Output:
#   - phyloseq object with user-friendly sample names
#
# Key Fields:
#   - sample_id: Internal CUID (unique identifier)
#   - sample_name: User-provided name (what you entered in eDNA Explorer)
#
# =============================================================================

library(phyloseq)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

BIOM_FILE <- "biom/16S_Bacteria-asv.biom"

# -----------------------------------------------------------------------------
# Import and Explore Sample IDs
# -----------------------------------------------------------------------------

cat("=== Importing BIOM File ===\n\n")

ps <- import_biom(BIOM_FILE)

# View the current sample identifiers
cat("Current sample names (phyloseq uses these internally):\n")
print(head(sample_names(ps), 10))

cat("\n\nSample metadata fields:\n")
print(sample_variables(ps))

# View the mapping between internal IDs and user names
cat("\n\nID Mapping (first 10 samples):\n")
id_mapping <- data.frame(
    phyloseq_name = sample_names(ps),
    sample_id = sample_data(ps)$sample_id,
    sample_name = sample_data(ps)$sample_name,
    stringsAsFactors = FALSE
)
print(head(id_mapping, 10))

# -----------------------------------------------------------------------------
# Option 1: Rename Samples to User-Provided Names
# -----------------------------------------------------------------------------

cat("\n=== Option 1: Rename to User-Provided Names ===\n\n")

# Make a copy to preserve original
ps_renamed <- ps

# Replace sample names with user-provided names
# WARNING: sample_name should be unique. Check first!
user_names <- sample_data(ps)$sample_name

if (length(unique(user_names)) == length(user_names)) {
    sample_names(ps_renamed) <- user_names
    cat("✓ Sample names updated successfully\n")
    cat("\nNew sample names (first 10):\n")
    print(head(sample_names(ps_renamed), 10))
} else {
    cat("⚠ Warning: sample_name values are not unique!\n")
    cat("Duplicate names found:\n")
    dup_names <- user_names[duplicated(user_names)]
    print(unique(dup_names))
    cat("\nKeeping original sample names to avoid confusion.\n")
}

# -----------------------------------------------------------------------------
# Option 2: Create a Lookup Function
# -----------------------------------------------------------------------------

cat("\n=== Option 2: Lookup Functions ===\n\n")

#' Get user-friendly name for a sample
#' @param ps phyloseq object
#' @param internal_id Internal sample ID (CUID)
#' @return User-provided sample name
get_sample_name <- function(ps, internal_id) {
    idx <- which(sample_names(ps) == internal_id)
    if (length(idx) == 0) {
        return(NA)
    }
    sample_data(ps)$sample_name[idx]
}

#' Get internal ID for a user-provided name
#' @param ps phyloseq object
#' @param user_name User-provided sample name
#' @return Internal sample ID (first match if duplicates exist)
get_internal_id <- function(ps, user_name) {
    idx <- which(sample_data(ps)$sample_name == user_name)
    if (length(idx) == 0) {
        return(NA)
    }
    sample_names(ps)[idx[1]]  # Return first match
}

# Example usage
example_internal <- sample_names(ps)[1]
example_name <- get_sample_name(ps, example_internal)
cat("Example lookup:\n")
cat("  Internal ID:", example_internal, "\n")
cat("  User name:", example_name, "\n")

# Reverse lookup
found_id <- get_internal_id(ps, example_name)
cat("  Reverse lookup:", found_id, "\n")

# -----------------------------------------------------------------------------
# Option 3: Create Display Labels for Plots
# -----------------------------------------------------------------------------

cat("\n=== Option 3: Display Labels for Plots ===\n\n")

# Add a display_label column combining name and location
sample_data(ps)$display_label <- paste0(
    sample_data(ps)$sample_name,
    " (",
    sample_data(ps)$country,
    ")"
)

cat("Created display_label column:\n")
print(head(sample_data(ps)[, c("sample_name", "country", "display_label")]))

# Use in plots:
# plot_bar(ps, x = "display_label", fill = "Phylum")

# -----------------------------------------------------------------------------
# Option 4: Export Mapping Table
# -----------------------------------------------------------------------------

cat("\n=== Option 4: Export Mapping Table ===\n\n")

# Create comprehensive mapping table
mapping_table <- data.frame(
    internal_id = sample_names(ps),
    sample_id = sample_data(ps)$sample_id,
    sample_name = sample_data(ps)$sample_name,
    latitude = sample_data(ps)$latitude,
    longitude = sample_data(ps)$longitude,
    country = sample_data(ps)$country,
    stringsAsFactors = FALSE
)

# Export to CSV
# write.csv(mapping_table, "sample_id_mapping.csv", row.names = FALSE)

cat("Mapping table created with", nrow(mapping_table), "samples\n")
cat("Columns:", paste(names(mapping_table), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Working with Both ID Systems
# -----------------------------------------------------------------------------

cat("\n=== Tips for Working with Both ID Systems ===\n\n")

cat("1. KEEP internal IDs for data integrity:\n")
cat("   - Use internal IDs when merging datasets\n")
cat("   - Use internal IDs for programmatic operations\n\n")

cat("2. USE user names for display:\n")
cat("   - Rename samples only for final visualizations\n")
cat("   - Or use the display_label approach\n\n")

cat("3. DOCUMENT your mapping:\n")
cat("   - Export mapping table for reference\n")
cat("   - Include in supplementary materials\n")

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

# Save renamed version
# saveRDS(ps_renamed, "phyloseq_renamed.rds")

# Save original with display labels
# saveRDS(ps, "phyloseq_with_labels.rds")

cat("\n✓ Sample ID mapping complete\n")
