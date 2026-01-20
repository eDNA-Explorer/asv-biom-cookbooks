# =============================================================================
# 07_taxonomy_import.R
# Taxonomy Handling and Visualization
# =============================================================================
#
# Description:
#   Work with taxonomy data from eDNA Explorer BIOM files. This includes:
#   - Understanding the taxonomy format (Greengenes-style prefixes)
#   - Parsing and cleaning taxonomy strings
#   - Aggregating at different taxonomic levels
#   - Visualizing community composition
#
# Prerequisites:
#   - phyloseq package installed
#   - ggplot2 for visualization
#
# Input:
#   - BIOM file (ASV or Taxa version)
#
# Taxonomy Format:
#   eDNA Explorer uses Greengenes-style prefixes:
#   k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; ...
#
# =============================================================================

library(phyloseq)
library(ggplot2)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Use the taxa BIOM for pre-collapsed taxonomy
# Or the ASV BIOM for sequence-level analysis
BIOM_FILE <- "biom/16S_Bacteria-taxa.biom"
# BIOM_FILE <- "biom/16S_Bacteria-asv.biom"  # Alternative

# -----------------------------------------------------------------------------
# Import and Explore Taxonomy
# -----------------------------------------------------------------------------

cat("=== Importing BIOM File ===\n\n")

ps <- import_biom(BIOM_FILE)

print(ps)

# Check if taxonomy is present
if (is.null(tax_table(ps, errorIfNULL = FALSE))) {
    stop("No taxonomy table found in BIOM file")
}

cat("\n--- Taxonomy Table Structure ---\n")
cat("Taxa:", ntaxa(ps), "\n")
cat("Ranks:", paste(rank_names(ps), collapse = ", "), "\n")

cat("\nFirst few taxonomy entries:\n")
print(head(tax_table(ps), 5))

# -----------------------------------------------------------------------------
# Understanding the Taxonomy Format
# -----------------------------------------------------------------------------

cat("\n=== Understanding Taxonomy Format ===\n\n")

# eDNA Explorer uses Greengenes-style prefixes
# Example: "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria"

# View raw taxonomy for first taxon
first_tax <- tax_table(ps)[1, ]
cat("Example taxonomy entry:\n")
print(first_tax)

# Check for prefixes
cat("\nPrefix examples from data:\n")
for (rank in rank_names(ps)) {
    example <- tax_table(ps)[1, rank]
    cat("  ", rank, ":", example, "\n")
}

# -----------------------------------------------------------------------------
# Function: Clean Taxonomy Prefixes
# -----------------------------------------------------------------------------

#' Remove Greengenes-style prefixes from taxonomy
#' @param ps phyloseq object
#' @return phyloseq object with cleaned taxonomy
clean_taxonomy_prefixes <- function(ps) {
    # Get tax table as matrix
    tax_mat <- as(tax_table(ps), "matrix")

    # Remove prefixes like "k__", "p__", etc.
    tax_mat <- gsub("^[a-z]__", "", tax_mat)

    # Replace empty strings with NA
    tax_mat[tax_mat == ""] <- NA

    # Update tax table
    tax_table(ps) <- tax_table(tax_mat)

    return(ps)
}

# Clean the taxonomy
ps_clean <- clean_taxonomy_prefixes(ps)

cat("\n--- Cleaned Taxonomy ---\n")
print(head(tax_table(ps_clean), 5))

# -----------------------------------------------------------------------------
# Aggregate at Different Taxonomic Levels
# -----------------------------------------------------------------------------

cat("\n=== Aggregating at Taxonomic Levels ===\n\n")

# Aggregate to Phylum level
ps_phylum <- tax_glom(ps_clean, taxrank = "Rank2")  # Phylum is usually Rank2
cat("Phylum-level aggregation:", ntaxa(ps_phylum), "taxa\n")

# Aggregate to Family level (usually Rank5)
ps_family <- tax_glom(ps_clean, taxrank = "Rank5")
cat("Family-level aggregation:", ntaxa(ps_family), "taxa\n")

# Aggregate to Genus level (usually Rank6)
ps_genus <- tax_glom(ps_clean, taxrank = "Rank6")
cat("Genus-level aggregation:", ntaxa(ps_genus), "taxa\n")

# -----------------------------------------------------------------------------
# Rename Taxonomy Ranks (Optional)
# -----------------------------------------------------------------------------

cat("\n=== Renaming Taxonomy Ranks ===\n\n")

# Phyloseq imports BIOM taxonomy as "Rank1", "Rank2", etc.
# Rename to standard names for clarity

# Define standard rank names
standard_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Only rename if we have the expected number of ranks
current_ranks <- rank_names(ps_clean)
if (length(current_ranks) == length(standard_ranks)) {
    colnames(tax_table(ps_clean)) <- standard_ranks
    cat("Renamed ranks to:", paste(rank_names(ps_clean), collapse = ", "), "\n")
} else if (length(current_ranks) <= length(standard_ranks)) {
    colnames(tax_table(ps_clean)) <- standard_ranks[1:length(current_ranks)]
    cat("Renamed ranks to:", paste(rank_names(ps_clean), collapse = ", "), "\n")
} else {
    cat("Keeping original rank names (unexpected number of ranks)\n")
}

# -----------------------------------------------------------------------------
# Visualize Taxonomy
# -----------------------------------------------------------------------------

cat("\n=== Taxonomy Visualization ===\n\n")

# Transform to relative abundance
ps_rel <- transform_sample_counts(ps_clean, function(x) x / sum(x))

# --- Phylum-level bar plot ---

# Aggregate and get top phyla
ps_phylum_rel <- tax_glom(ps_rel, taxrank = "Phylum")

# Get top 10 phyla by mean abundance
phylum_means <- taxa_sums(ps_phylum_rel) / nsamples(ps_phylum_rel)
top_phyla <- names(sort(phylum_means, decreasing = TRUE))[1:10]

# Prune to top phyla
ps_top_phyla <- prune_taxa(top_phyla, ps_phylum_rel)

# Create bar plot
p1 <- plot_bar(ps_top_phyla, fill = "Phylum") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right"
    ) +
    labs(
        title = "Community Composition by Phylum",
        x = "Sample",
        y = "Relative Abundance"
    )

# Save plot
# ggsave("phylum_barplot.png", p1, width = 12, height = 6)
cat("Phylum bar plot created\n")

# --- Phylum abundance heatmap ---

p2 <- plot_heatmap(
    ps_top_phyla,
    taxa.label = "Phylum",
    method = NULL,
    low = "white",
    high = "steelblue"
) +
    labs(title = "Phylum Abundance Heatmap")

# ggsave("phylum_heatmap.png", p2, width = 10, height = 8)
cat("Phylum heatmap created\n")

# -----------------------------------------------------------------------------
# Taxonomy Summary Statistics
# -----------------------------------------------------------------------------

cat("\n=== Taxonomy Summary ===\n\n")

# Count unique taxa at each level
for (rank in rank_names(ps_clean)) {
    unique_taxa <- length(unique(tax_table(ps_clean)[, rank]))
    cat(rank, ":", unique_taxa, "unique\n")
}

# Find most abundant taxa at each level
cat("\n--- Most Abundant at Each Level ---\n")

for (rank in c("Phylum", "Family", "Genus")) {
    if (rank %in% rank_names(ps_clean)) {
        ps_agg <- tax_glom(ps_clean, taxrank = rank)

        # Get total abundance
        abundances <- taxa_sums(ps_agg)
        top_taxon_idx <- which.max(abundances)
        top_taxon_name <- tax_table(ps_agg)[top_taxon_idx, rank]
        top_taxon_pct <- round(100 * abundances[top_taxon_idx] / sum(abundances), 1)

        cat(rank, ": ", top_taxon_name, " (", top_taxon_pct, "%)\n", sep = "")
    }
}

# -----------------------------------------------------------------------------
# Export Taxonomy Table
# -----------------------------------------------------------------------------

cat("\n=== Export Options ===\n\n")

# Export taxonomy table as CSV
tax_df <- as.data.frame(tax_table(ps_clean))
tax_df$feature_id <- rownames(tax_df)

# Add abundance information
tax_df$total_reads <- taxa_sums(ps_clean)
tax_df$prevalence <- rowSums(otu_table(ps_clean) > 0)

# Reorder columns
tax_df <- tax_df[, c("feature_id", rank_names(ps_clean), "total_reads", "prevalence")]

# write.csv(tax_df, "taxonomy_summary.csv", row.names = FALSE)
cat("Taxonomy table ready for export (", nrow(tax_df), " taxa)\n", sep = "")

# Preview
cat("\nTaxonomy table preview:\n")
print(head(tax_df[order(-tax_df$total_reads), ], 10))

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

# saveRDS(ps_clean, "phyloseq_clean_taxonomy.rds")

cat("\n✓ Taxonomy import and processing complete\n")
