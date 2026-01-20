# =============================================================================
# complete_workflow.R
# Complete Phyloseq Analysis Workflow for eDNA Explorer Data
# =============================================================================
#
# Description:
#   A comprehensive end-to-end workflow for analyzing eDNA Explorer BIOM files
#   in Phyloseq. This script demonstrates:
#   - Data import with sequences
#   - Quality filtering
#   - Normalization
#   - Alpha diversity analysis
#   - Beta diversity analysis (ordination)
#   - Taxonomic composition visualization
#   - Statistical testing
#
# Prerequisites:
#   - phyloseq, Biostrings, ggplot2, vegan, dplyr
#
# Input:
#   - BIOM file (ASV or Taxa level)
#   - Forward FASTA file (optional)
#   - Sample metadata (embedded or external)
#
# =============================================================================

# Load packages
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(vegan)
library(dplyr)

# Set theme for plots
theme_set(theme_bw())

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# File paths - adjust to match your data
BIOM_FILE <- "biom/16S_Bacteria-asv.biom"
FORWARD_FASTA <- "fasta/16S_Bacteria_paired_F.fasta"  # Set to NULL if not using

# Analysis parameters
MIN_READS_PER_SAMPLE <- 1000   # Remove samples with fewer reads
MIN_READS_PER_ASV <- 10        # Remove ASVs with fewer total reads
MIN_PREVALENCE <- 2            # Remove ASVs found in fewer than N samples
RAREFACTION_DEPTH <- 5000      # Depth for rarefaction (NULL for no rarefaction)

# Grouping variable for comparisons (must exist in sample metadata)
GROUP_VARIABLE <- "country"    # Change to match your metadata

# Output directory
OUTPUT_DIR <- "results"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# =============================================================================
# PART 1: DATA IMPORT
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 1: DATA IMPORT\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Import BIOM file
if (!is.null(FORWARD_FASTA) && file.exists(FORWARD_FASTA)) {
    cat("Importing BIOM with forward sequences...\n")
    ps <- import_biom(
        BIOMfilename = BIOM_FILE,
        refseqfilename = FORWARD_FASTA,
        refseqFunction = readDNAStringSet
    )
} else {
    cat("Importing BIOM file only...\n")
    ps <- import_biom(BIOM_FILE)
}

# Print summary
cat("\n--- Initial Data Summary ---\n")
print(ps)

cat("\nSample metadata variables:\n")
print(sample_variables(ps))

# Verify grouping variable exists
if (!GROUP_VARIABLE %in% sample_variables(ps)) {
    warning(paste("GROUP_VARIABLE '", GROUP_VARIABLE, "' not found in sample data.",
                  "Using first available variable instead."))
    GROUP_VARIABLE <- sample_variables(ps)[1]
}

cat("\nGroup variable:", GROUP_VARIABLE, "\n")
cat("Group levels:", paste(unique(sample_data(ps)[[GROUP_VARIABLE]]), collapse = ", "), "\n")

# =============================================================================
# PART 2: QUALITY FILTERING
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 2: QUALITY FILTERING\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

ps_filtered <- ps

# --- Filter low-depth samples ---
sample_depths <- sample_sums(ps_filtered)
cat("Sample read depths:\n")
cat("  Min:", min(sample_depths), "\n")
cat("  Max:", max(sample_depths), "\n")
cat("  Median:", median(sample_depths), "\n")

low_depth_samples <- sample_depths < MIN_READS_PER_SAMPLE
if (sum(low_depth_samples) > 0) {
    cat("\nRemoving", sum(low_depth_samples), "samples with <", MIN_READS_PER_SAMPLE, "reads\n")
    ps_filtered <- prune_samples(!low_depth_samples, ps_filtered)
}

# --- Filter low-abundance ASVs ---
asv_totals <- taxa_sums(ps_filtered)
low_abundance_asvs <- asv_totals < MIN_READS_PER_ASV
if (sum(low_abundance_asvs) > 0) {
    cat("Removing", sum(low_abundance_asvs), "ASVs with <", MIN_READS_PER_ASV, "total reads\n")
    ps_filtered <- prune_taxa(!low_abundance_asvs, ps_filtered)
}

# --- Filter low-prevalence ASVs ---
prevalence <- rowSums(otu_table(ps_filtered) > 0)
low_prevalence_asvs <- prevalence < MIN_PREVALENCE
if (sum(low_prevalence_asvs) > 0) {
    cat("Removing", sum(low_prevalence_asvs), "ASVs found in <", MIN_PREVALENCE, "samples\n")
    ps_filtered <- prune_taxa(!low_prevalence_asvs, ps_filtered)
}

cat("\n--- After Filtering ---\n")
cat("Samples:", nsamples(ps_filtered), "(was", nsamples(ps), ")\n")
cat("Taxa:", ntaxa(ps_filtered), "(was", ntaxa(ps), ")\n")

# =============================================================================
# PART 3: NORMALIZATION
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 3: NORMALIZATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# --- Rarefaction (optional) ---
if (!is.null(RAREFACTION_DEPTH)) {
    # Check which samples can be rarefied
    can_rarefy <- sample_sums(ps_filtered) >= RAREFACTION_DEPTH
    cat("Samples that can be rarefied to", RAREFACTION_DEPTH, ":", sum(can_rarefy), "\n")

    if (sum(can_rarefy) < nsamples(ps_filtered)) {
        cat("Removing", sum(!can_rarefy), "samples below rarefaction depth\n")
        ps_filtered <- prune_samples(can_rarefy, ps_filtered)
    }

    set.seed(42)  # For reproducibility
    ps_rarefied <- rarefy_even_depth(ps_filtered, sample.size = RAREFACTION_DEPTH, verbose = FALSE)
    cat("Rarefied to", RAREFACTION_DEPTH, "reads per sample\n")
} else {
    ps_rarefied <- ps_filtered
    cat("Skipping rarefaction (RAREFACTION_DEPTH = NULL)\n")
}

# --- Relative abundance transformation ---
ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))
cat("Created relative abundance version\n")

# =============================================================================
# PART 4: ALPHA DIVERSITY
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 4: ALPHA DIVERSITY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Calculate alpha diversity metrics
alpha_div <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$sample_id <- rownames(alpha_div)
alpha_div[[GROUP_VARIABLE]] <- sample_data(ps_rarefied)[[GROUP_VARIABLE]]

cat("--- Alpha Diversity Summary by", GROUP_VARIABLE, "---\n")
alpha_summary <- alpha_div %>%
    group_by(!!sym(GROUP_VARIABLE)) %>%
    summarise(
        n = n(),
        mean_observed = round(mean(Observed), 1),
        sd_observed = round(sd(Observed), 1),
        mean_shannon = round(mean(Shannon), 2),
        sd_shannon = round(sd(Shannon), 2),
        .groups = "drop"
    )
print(as.data.frame(alpha_summary))

# Plot alpha diversity
p_alpha <- plot_richness(ps_rarefied, x = GROUP_VARIABLE, measures = c("Observed", "Shannon")) +
    geom_boxplot(alpha = 0.6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Alpha Diversity")

ggsave(file.path(OUTPUT_DIR, "alpha_diversity.png"), p_alpha, width = 10, height = 6)
cat("\nSaved: alpha_diversity.png\n")

# --- Statistical test (Kruskal-Wallis) ---
cat("\n--- Kruskal-Wallis Test for Shannon Diversity ---\n")
kw_test <- kruskal.test(Shannon ~ get(GROUP_VARIABLE), data = alpha_div)
print(kw_test)

# =============================================================================
# PART 5: BETA DIVERSITY
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 5: BETA DIVERSITY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# --- Calculate distance matrix ---
cat("Calculating Bray-Curtis distances...\n")
dist_bray <- phyloseq::distance(ps_rel, method = "bray")

# --- Ordination (PCoA) ---
cat("Running PCoA ordination...\n")
ord_pcoa <- ordinate(ps_rel, method = "PCoA", distance = dist_bray)

# Plot ordination
p_ord <- plot_ordination(ps_rel, ord_pcoa, color = GROUP_VARIABLE) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95) +
    labs(title = "PCoA Ordination (Bray-Curtis)")

ggsave(file.path(OUTPUT_DIR, "beta_diversity_pcoa.png"), p_ord, width = 8, height = 6)
cat("Saved: beta_diversity_pcoa.png\n")

# --- PERMANOVA ---
cat("\n--- PERMANOVA Test ---\n")
sample_df <- data.frame(sample_data(ps_rel))
permanova_result <- adonis2(
    dist_bray ~ get(GROUP_VARIABLE),
    data = sample_df,
    permutations = 999
)
print(permanova_result)

# --- NMDS (alternative ordination) ---
cat("\nRunning NMDS ordination...\n")
set.seed(42)
ord_nmds <- ordinate(ps_rel, method = "NMDS", distance = dist_bray)
cat("NMDS stress:", ord_nmds$stress, "\n")

p_nmds <- plot_ordination(ps_rel, ord_nmds, color = GROUP_VARIABLE) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95) +
    labs(title = paste("NMDS Ordination (stress =", round(ord_nmds$stress, 3), ")"))

ggsave(file.path(OUTPUT_DIR, "beta_diversity_nmds.png"), p_nmds, width = 8, height = 6)
cat("Saved: beta_diversity_nmds.png\n")

# =============================================================================
# PART 6: TAXONOMIC COMPOSITION
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 6: TAXONOMIC COMPOSITION\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# --- Clean taxonomy prefixes ---
clean_tax <- function(ps) {
    if (is.null(tax_table(ps, errorIfNULL = FALSE))) return(ps)
    tax_mat <- as(tax_table(ps), "matrix")
    tax_mat <- gsub("^[a-z]__", "", tax_mat)
    tax_mat[tax_mat == ""] <- NA
    tax_table(ps) <- tax_table(tax_mat)
    return(ps)
}

ps_rel_clean <- clean_tax(ps_rel)

# --- Phylum-level composition ---
if (!is.null(tax_table(ps_rel_clean, errorIfNULL = FALSE))) {
    cat("Aggregating to Phylum level...\n")

    # Get the phylum rank name (usually Rank2 or Phylum)
    rank_names_ps <- rank_names(ps_rel_clean)
    phylum_rank <- ifelse("Phylum" %in% rank_names_ps, "Phylum", rank_names_ps[2])

    ps_phylum <- tax_glom(ps_rel_clean, taxrank = phylum_rank, NArm = FALSE)

    # Get top 10 phyla
    top_phyla <- names(sort(taxa_sums(ps_phylum), decreasing = TRUE))[1:10]
    ps_top_phyla <- prune_taxa(top_phyla, ps_phylum)

    # Bar plot by group
    p_taxa <- plot_bar(ps_top_phyla, x = "Sample", fill = phylum_rank) +
        facet_wrap(as.formula(paste("~", GROUP_VARIABLE)), scales = "free_x") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        labs(
            title = "Community Composition by Phylum",
            x = "Samples",
            y = "Relative Abundance"
        )

    ggsave(file.path(OUTPUT_DIR, "taxonomic_composition.png"), p_taxa, width = 12, height = 6)
    cat("Saved: taxonomic_composition.png\n")

    # --- Mean abundance by group ---
    cat("\n--- Mean Phylum Abundance by", GROUP_VARIABLE, "---\n")

    # Melt to long format for summarization
    ps_melt <- psmelt(ps_top_phyla)

    phylum_summary <- ps_melt %>%
        group_by(!!sym(GROUP_VARIABLE), !!sym(phylum_rank)) %>%
        summarise(mean_abundance = round(mean(Abundance) * 100, 2), .groups = "drop") %>%
        arrange(!!sym(GROUP_VARIABLE), desc(mean_abundance))

    print(head(as.data.frame(phylum_summary), 20))
}

# =============================================================================
# PART 7: EXPORT RESULTS
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 7: EXPORT RESULTS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# --- Export alpha diversity ---
write.csv(alpha_div, file.path(OUTPUT_DIR, "alpha_diversity.csv"), row.names = FALSE)
cat("Saved: alpha_diversity.csv\n")

# --- Export OTU table ---
otu_df <- as.data.frame(otu_table(ps_filtered))
otu_df$feature_id <- rownames(otu_df)
write.csv(otu_df, file.path(OUTPUT_DIR, "otu_table.csv"), row.names = FALSE)
cat("Saved: otu_table.csv\n")

# --- Export taxonomy ---
if (!is.null(tax_table(ps_filtered, errorIfNULL = FALSE))) {
    tax_df <- as.data.frame(tax_table(ps_filtered))
    tax_df$feature_id <- rownames(tax_df)
    write.csv(tax_df, file.path(OUTPUT_DIR, "taxonomy.csv"), row.names = FALSE)
    cat("Saved: taxonomy.csv\n")
}

# --- Export sample metadata ---
sample_df <- as.data.frame(sample_data(ps_filtered))
sample_df$sample_id <- rownames(sample_df)
write.csv(sample_df, file.path(OUTPUT_DIR, "sample_metadata.csv"), row.names = FALSE)
cat("Saved: sample_metadata.csv\n")

# --- Save phyloseq object ---
saveRDS(ps_filtered, file.path(OUTPUT_DIR, "phyloseq_filtered.rds"))
cat("Saved: phyloseq_filtered.rds\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Final dataset:\n")
cat("  Samples:", nsamples(ps_filtered), "\n")
cat("  Taxa:", ntaxa(ps_filtered), "\n")

cat("\nOutput files saved to:", OUTPUT_DIR, "/\n")
cat("  - alpha_diversity.csv\n")
cat("  - alpha_diversity.png\n")
cat("  - beta_diversity_pcoa.png\n")
cat("  - beta_diversity_nmds.png\n")
cat("  - taxonomic_composition.png\n")
cat("  - otu_table.csv\n")
cat("  - taxonomy.csv\n")
cat("  - sample_metadata.csv\n")
cat("  - phyloseq_filtered.rds\n")

cat("\n✓ Workflow complete!\n")
