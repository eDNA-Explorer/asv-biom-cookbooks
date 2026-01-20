# =============================================================================
# 04_paired_sequences.R
# Import BIOM with Both Forward and Reverse Sequences
# =============================================================================
#
# Description:
#   Import an eDNA Explorer BIOM file with both forward and reverse sequences.
#   Forward sequences are stored in the refseq() slot (Phyloseq's native
#   sequence storage). Reverse sequences are stored in the taxonomy table
#   as an additional column.
#
# Prerequisites:
#   - phyloseq package installed
#   - Biostrings package installed
#
# Input:
#   - BIOM file (e.g., "biom/16S_Bacteria-asv.biom")
#   - Forward FASTA file (e.g., "fasta/16S_Bacteria_paired_F.fasta")
#   - Reverse FASTA file (e.g., "fasta/16S_Bacteria_paired_R.fasta")
#
# Output:
#   - phyloseq object with:
#     - Forward sequences in refseq() slot
#     - Reverse sequences in tax_table "reverse_sequence" column
#
# =============================================================================

library(phyloseq)
library(Biostrings)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Paths to your files
BIOM_FILE <- "biom/16S_Bacteria-asv.biom"
FORWARD_FASTA <- "fasta/16S_Bacteria_paired_F.fasta"
REVERSE_FASTA <- "fasta/16S_Bacteria_paired_R.fasta"

# -----------------------------------------------------------------------------
# Function: Load Paired Phyloseq
# -----------------------------------------------------------------------------

#' Load phyloseq with both forward and reverse sequences
#'
#' @param biom_file Path to BIOM file
#' @param forward_fasta Path to forward FASTA
#' @param reverse_fasta Path to reverse FASTA
#' @return phyloseq object with forward sequences in refseq slot
#'         and reverse sequences stored in tax_table
load_paired_phyloseq <- function(biom_file, forward_fasta, reverse_fasta) {

    cat("Loading BIOM with forward sequences...\n")

    # 1. Import BIOM with forward sequences
    ps <- import_biom(
        BIOMfilename = biom_file,
        refseqfilename = forward_fasta,
        refseqFunction = readDNAStringSet
    )

    cat("  Forward sequences loaded:", length(refseq(ps)), "\n")

    # 2. Load reverse sequences
    cat("Loading reverse sequences...\n")
    rev_seqs <- readDNAStringSet(reverse_fasta)
    cat("  Reverse sequences found:", length(rev_seqs), "\n")

    # 3. Transform reverse headers to match forward IDs
    cat("Transforming reverse headers (_paired_R_ -> _paired_F_)...\n")
    names(rev_seqs) <- gsub("_paired_R_", "_paired_F_", names(rev_seqs))

    # 4. Add reverse sequences to tax_table
    if (!is.null(tax_table(ps, errorIfNULL = FALSE))) {
        # Get taxa IDs from phyloseq object
        taxa_ids <- taxa_names(ps)

        # Extract reverse sequences for each taxon (NA if not found)
        rev_seq_chars <- sapply(taxa_ids, function(id) {
            if (id %in% names(rev_seqs)) {
                as.character(rev_seqs[id])
            } else {
                NA_character_
            }
        })

        # Add to existing tax_table as new column
        tax_df <- as.data.frame(tax_table(ps))
        tax_df$reverse_sequence <- rev_seq_chars
        tax_table(ps) <- tax_table(as.matrix(tax_df))

        matched <- sum(!is.na(rev_seq_chars))
        cat("  Reverse sequences matched:", matched, "of", length(taxa_ids), "\n")
    } else {
        warning("No tax_table found. Creating one for reverse sequences.")
        taxa_ids <- taxa_names(ps)
        rev_seq_chars <- sapply(taxa_ids, function(id) {
            if (id %in% names(rev_seqs)) {
                as.character(rev_seqs[id])
            } else {
                NA_character_
            }
        })
        tax_df <- data.frame(reverse_sequence = rev_seq_chars, row.names = taxa_ids)
        tax_table(ps) <- tax_table(as.matrix(tax_df))
    }

    return(ps)
}

# -----------------------------------------------------------------------------
# Load Paired Data
# -----------------------------------------------------------------------------

cat("=== Loading Paired Sequences ===\n\n")

ps <- load_paired_phyloseq(
    biom_file = BIOM_FILE,
    forward_fasta = FORWARD_FASTA,
    reverse_fasta = REVERSE_FASTA
)

# -----------------------------------------------------------------------------
# Verify Import
# -----------------------------------------------------------------------------

cat("\n=== Verification ===\n\n")
print(ps)

# Check forward sequences
cat("\n--- Forward Sequences (refseq slot) ---\n")
print(refseq(ps)[1:3])

# Check reverse sequences
cat("\n--- Reverse Sequences (tax_table column) ---\n")
rev_seqs <- tax_table(ps)[, "reverse_sequence"]
cat("First 3 reverse sequences:\n")
for (i in 1:min(3, nrow(rev_seqs))) {
    seq <- rev_seqs[i, 1]
    if (!is.na(seq) && nchar(seq) > 50) {
        cat(rownames(rev_seqs)[i], ":", substr(seq, 1, 50), "...\n")
    } else {
        cat(rownames(rev_seqs)[i], ":", seq, "\n")
    }
}

# -----------------------------------------------------------------------------
# Helper Functions for Accessing Sequences
# -----------------------------------------------------------------------------

#' Get forward sequence for a feature
#' @param ps phyloseq object
#' @param feature_id Feature ID (e.g., "16S_Bacteria_paired_F_0")
#' @return Character string of forward sequence
get_forward_seq <- function(ps, feature_id) {
    as.character(refseq(ps)[feature_id])
}

#' Get reverse sequence for a feature
#' @param ps phyloseq object
#' @param feature_id Feature ID (e.g., "16S_Bacteria_paired_F_0")
#' @return Character string of reverse sequence
get_reverse_seq <- function(ps, feature_id) {
    tax_table(ps)[feature_id, "reverse_sequence"]
}

#' Get both sequences for a feature
#' @param ps phyloseq object
#' @param feature_id Feature ID
#' @return Named list with forward and reverse sequences
get_paired_seqs <- function(ps, feature_id) {
    list(
        forward = get_forward_seq(ps, feature_id),
        reverse = get_reverse_seq(ps, feature_id)
    )
}

# -----------------------------------------------------------------------------
# Example Usage
# -----------------------------------------------------------------------------

cat("\n=== Example: Accessing Paired Sequences ===\n\n")

# Get first feature ID
example_id <- taxa_names(ps)[1]
cat("Feature ID:", example_id, "\n\n")

# Get paired sequences
seqs <- get_paired_seqs(ps, example_id)

cat("Forward sequence (first 60 bp):\n")
cat(" ", substr(seqs$forward, 1, 60), "...\n\n")

cat("Reverse sequence (first 60 bp):\n")
cat(" ", substr(seqs$reverse, 1, 60), "...\n")

# -----------------------------------------------------------------------------
# Save for Later Use (Optional)
# -----------------------------------------------------------------------------

# saveRDS(ps, "phyloseq_paired.rds")

cat("\n✓ BIOM file with paired sequences successfully imported\n")
