# =============================================================================
# 03_reverse_sequences.R
# Import BIOM with Reverse Reference Sequences (Preprocessing Required)
# =============================================================================
#
# Description:
#   Import an eDNA Explorer BIOM file with reverse FASTA sequences.
#
#   IMPORTANT: Reverse FASTA headers use "_paired_R_" while BIOM feature IDs
#   use "_paired_F_". This script transforms the headers before import.
#
# Prerequisites:
#   - phyloseq package installed
#   - Biostrings package installed
#
# Input:
#   - BIOM file (e.g., "biom/16S_Bacteria-asv.biom")
#   - Reverse FASTA file (e.g., "fasta/16S_Bacteria_paired_R.fasta")
#
# Output:
#   - phyloseq object with reverse sequences in refseq slot
#   - Transformed FASTA file (temporary or saved)
#
# =============================================================================

library(phyloseq)
library(Biostrings)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Paths to your files
BIOM_FILE <- "biom/16S_Bacteria-asv.biom"
REVERSE_FASTA <- "fasta/16S_Bacteria_paired_R.fasta"

# Output path for transformed FASTA (set to NULL to use temp file)
TRANSFORMED_FASTA <- "fasta/16S_Bacteria_paired_R_transformed.fasta"
# TRANSFORMED_FASTA <- NULL  # Use temp file instead

# -----------------------------------------------------------------------------
# Function: Preprocess Reverse FASTA
# -----------------------------------------------------------------------------

#' Preprocess reverse FASTA to match BIOM feature IDs
#'
#' Transforms FASTA headers from _paired_R_ to _paired_F_ format
#' so they match the BIOM file's feature IDs.
#'
#' @param input_fasta Path to reverse FASTA file
#' @param output_fasta Path for output FASTA (NULL for temp file)
#' @return Path to the transformed FASTA file
preprocess_reverse_fasta <- function(input_fasta, output_fasta = NULL) {

    # Create temp file if no output path specified
    if (is.null(output_fasta)) {
        output_fasta <- tempfile(fileext = ".fasta")
    }

    # Read reverse sequences
    cat("Reading reverse FASTA:", input_fasta, "\n")
    rev_seqs <- readDNAStringSet(input_fasta)

    original_count <- length(rev_seqs)
    cat("Sequences found:", original_count, "\n")

    # Show example of transformation
    if (original_count > 0) {
        cat("\nHeader transformation example:\n")
        cat("  Before:", names(rev_seqs)[1], "\n")
    }

    # Transform headers: _paired_R_ -> _paired_F_
    names(rev_seqs) <- gsub("_paired_R_", "_paired_F_", names(rev_seqs))

    if (original_count > 0) {
        cat("  After: ", names(rev_seqs)[1], "\n")
    }

    # Write transformed FASTA
    writeXStringSet(rev_seqs, output_fasta)
    cat("\nTransformed FASTA written to:", output_fasta, "\n")

    return(output_fasta)
}

# -----------------------------------------------------------------------------
# Step 1: Preprocess Reverse FASTA
# -----------------------------------------------------------------------------

cat("=== Step 1: Preprocessing Reverse FASTA ===\n\n")

transformed_fasta <- preprocess_reverse_fasta(
    input_fasta = REVERSE_FASTA,
    output_fasta = TRANSFORMED_FASTA
)

# -----------------------------------------------------------------------------
# Step 2: Import BIOM with Transformed Sequences
# -----------------------------------------------------------------------------

cat("\n=== Step 2: Importing BIOM with Reverse Sequences ===\n\n")

ps <- import_biom(
    BIOMfilename = BIOM_FILE,
    refseqfilename = transformed_fasta,
    refseqFunction = readDNAStringSet
)

# -----------------------------------------------------------------------------
# Verify Import
# -----------------------------------------------------------------------------

print(ps)

# Check sequence matching
cat("\n--- Sequence Matching ---\n")
biom_ids <- taxa_names(ps)
fasta_ids <- names(refseq(ps))
matching <- sum(biom_ids %in% fasta_ids)

cat("BIOM features:", length(biom_ids), "\n")
cat("FASTA sequences:", length(fasta_ids), "\n")
cat("Successfully matched:", matching, "\n")

if (matching == length(biom_ids)) {
    cat("\n✓ All features have matching reverse sequences\n")
} else {
    cat("\nWarning:", length(biom_ids) - matching, "features without sequences\n")
}

# -----------------------------------------------------------------------------
# View Example Reverse Sequences
# -----------------------------------------------------------------------------

cat("\n--- Example Reverse Sequences ---\n")
print(refseq(ps)[1:3])

# -----------------------------------------------------------------------------
# Cleanup (Optional)
# -----------------------------------------------------------------------------

# If using a temp file and you want to clean up:
# if (is.null(TRANSFORMED_FASTA)) {
#     unlink(transformed_fasta)
#     cat("\nTemp file cleaned up\n")
# }

cat("\n✓ BIOM file with reverse sequences successfully imported\n")
