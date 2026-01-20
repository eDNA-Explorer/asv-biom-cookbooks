# =============================================================================
# 02_forward_sequences.R
# Import BIOM with Forward Reference Sequences
# =============================================================================
#
# Description:
#   Import an eDNA Explorer BIOM file along with forward FASTA sequences.
#   Forward FASTA headers match BIOM feature IDs directly - no preprocessing
#   required.
#
# Prerequisites:
#   - phyloseq package installed
#   - Biostrings package installed (for sequence handling)
#
# Input:
#   - BIOM file (e.g., "biom/16S_Bacteria-asv.biom")
#   - Forward FASTA file (e.g., "fasta/16S_Bacteria_paired_F.fasta")
#
# Output:
#   - phyloseq object with otu_table, sample_data, tax_table, and refseq
#
# =============================================================================

library(phyloseq)
library(Biostrings)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Paths to your files
# Adjust these to match your file locations
BIOM_FILE <- "biom/16S_Bacteria-asv.biom"
FORWARD_FASTA <- "fasta/16S_Bacteria_paired_F.fasta"

# -----------------------------------------------------------------------------
# Import BIOM with Forward Sequences
# -----------------------------------------------------------------------------

# Import BIOM file with forward reference sequences
# The refseqFunction parameter tells phyloseq how to read the FASTA file
ps <- import_biom(
    BIOMfilename = BIOM_FILE,
    refseqfilename = FORWARD_FASTA,
    refseqFunction = readDNAStringSet
)

# -----------------------------------------------------------------------------
# Verify Import
# -----------------------------------------------------------------------------

# Print summary
print(ps)
# Example output:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5000 taxa and 200 samples ]
# sample_data() Sample Data:       [ 200 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 5000 taxa by 7 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 5000 reference sequences ]

# Check that sequences were loaded
cat("\nReference sequences loaded:", !is.null(refseq(ps)), "\n")

# -----------------------------------------------------------------------------
# Explore Sequences
# -----------------------------------------------------------------------------

# View first few sequences
cat("\n--- Reference Sequences (first 5) ---\n")
print(refseq(ps)[1:5])

# Sequence lengths
seq_lengths <- width(refseq(ps))
cat("\n--- Sequence Length Statistics ---\n")
cat("Min length:", min(seq_lengths), "bp\n")
cat("Max length:", max(seq_lengths), "bp\n")
cat("Mean length:", round(mean(seq_lengths), 1), "bp\n")

# Verify feature IDs match
cat("\n--- Feature ID Matching ---\n")
biom_ids <- taxa_names(ps)
fasta_ids <- names(refseq(ps))
matching <- sum(biom_ids %in% fasta_ids)
cat("BIOM feature IDs:", length(biom_ids), "\n")
cat("FASTA sequence IDs:", length(fasta_ids), "\n")
cat("Matching IDs:", matching, "\n")

if (matching < length(biom_ids)) {
    cat("\nWarning: Not all BIOM features have matching sequences\n")
    cat("Missing:", length(biom_ids) - matching, "sequences\n")
}

# -----------------------------------------------------------------------------
# Example: Extract Sequence for a Specific Feature
# -----------------------------------------------------------------------------

# Get sequence for the first feature
first_feature <- taxa_names(ps)[1]
first_sequence <- as.character(refseq(ps)[first_feature])

cat("\n--- Example Feature ---\n")
cat("Feature ID:", first_feature, "\n")
cat("Sequence (first 50 bp):", substr(first_sequence, 1, 50), "...\n")
cat("Full length:", nchar(first_sequence), "bp\n")

# -----------------------------------------------------------------------------
# Example: Filter by Sequence Length
# -----------------------------------------------------------------------------

# Keep only sequences within expected length range (e.g., 200-300 bp)
min_length <- 200
max_length <- 300

seq_lengths <- width(refseq(ps))
length_filter <- seq_lengths >= min_length & seq_lengths <= max_length

cat("\n--- Length Filtering ---\n")
cat("Sequences in range", min_length, "-", max_length, "bp:",
    sum(length_filter), "of", length(length_filter), "\n")

# Apply filter (uncomment to use)
# ps_filtered <- prune_taxa(length_filter, ps)

# -----------------------------------------------------------------------------
# Save for Later Use (Optional)
# -----------------------------------------------------------------------------

# Save the phyloseq object
# saveRDS(ps, "phyloseq_with_seqs.rds")

cat("\n✓ BIOM file with forward sequences successfully imported\n")
