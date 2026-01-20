# =============================================================================
# preprocess_reverse_fasta.R
# Transform reverse FASTA headers to match BIOM feature IDs
# =============================================================================
#
# Description:
#   eDNA Explorer BIOM files use forward read names as feature IDs:
#       16S_Bacteria_paired_F_0
#
#   Forward FASTA headers match directly:
#       >16S_Bacteria_paired_F_0
#
#   Reverse FASTA headers use _paired_R_:
#       >16S_Bacteria_paired_R_0
#
#   This script transforms reverse FASTA headers to match BIOM feature IDs
#   by replacing _paired_R_ with _paired_F_.
#
# Prerequisites:
#   - Biostrings package (from Bioconductor)
#
# Usage:
#   source("preprocess_reverse_fasta.R")
#   preprocess_reverse_fasta("input.fasta", "output.fasta")
#
# =============================================================================

# Check for Biostrings package
if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop(
        "Biostrings package required. Install with:\n",
        "  BiocManager::install('Biostrings')"
    )
}

library(Biostrings)


#' Preprocess reverse FASTA to match BIOM feature IDs
#'
#' Transforms FASTA headers from _paired_R_ to _paired_F_ format
#' so they match the BIOM file's feature IDs.
#'
#' @param input_fasta Path to input reverse FASTA file
#' @param output_fasta Path for output transformed FASTA file
#' @param verbose Print progress information (default: TRUE)
#'
#' @return Invisibly returns a list with statistics:
#'   - sequences: Total number of sequences
#'   - transformed: Number of headers transformed
#'   - unchanged: Number of headers unchanged
#'
#' @examples
#' preprocess_reverse_fasta(
#'     "fasta/16S_Bacteria_paired_R.fasta",
#'     "fasta/16S_Bacteria_paired_R_transformed.fasta"
#' )
#'
preprocess_reverse_fasta <- function(input_fasta, output_fasta, verbose = TRUE) {

    # Validate input
    if (!file.exists(input_fasta)) {
        stop("Input file not found: ", input_fasta)
    }

    # Create output directory if needed
    output_dir <- dirname(output_fasta)
    if (!dir.exists(output_dir) && output_dir != ".") {
        dir.create(output_dir, recursive = TRUE)
    }

    if (verbose) {
        cat("=" , rep("=", 58), "\n", sep = "")
        cat("  Reverse FASTA Preprocessor (R)\n")
        cat("=", rep("=", 58), "\n", sep = "")
        cat("\n")
        cat("Input: ", input_fasta, "\n")
        cat("Output:", output_fasta, "\n")
        cat("\n")
    }

    # Read sequences
    if (verbose) cat("Reading sequences...\n")
    seqs <- readDNAStringSet(input_fasta)

    total_seqs <- length(seqs)
    if (verbose) cat("  Found", total_seqs, "sequences\n")

    # Get original names for statistics
    original_names <- names(seqs)

    # Show example transformation
    if (verbose && total_seqs > 0) {
        first_name <- original_names[1]
        transformed_name <- gsub("_paired_R_", "_paired_F_", first_name)
        cat("\n")
        cat("Transformation:\n")
        cat("  Before:", first_name, "\n")
        cat("  After: ", transformed_name, "\n")
        cat("\n")
    }

    # Transform headers: _paired_R_ -> _paired_F_
    if (verbose) cat("Transforming headers...\n")
    new_names <- gsub("_paired_R_", "_paired_F_", original_names)
    names(seqs) <- new_names

    # Calculate statistics
    transformed <- sum(original_names != new_names)
    unchanged <- total_seqs - transformed

    if (verbose) {
        cat("  Transformed:", transformed, "\n")
        cat("  Unchanged:  ", unchanged, "\n")
    }

    # Write output
    if (verbose) cat("\nWriting output...\n")
    writeXStringSet(seqs, output_fasta)

    if (verbose) {
        cat("\n")
        cat("Output written to:", output_fasta, "\n")

        if (unchanged > 0) {
            cat("\n")
            cat("Note:", unchanged, "sequence(s) did not contain '_paired_R_'\n")
            cat("      These headers were written unchanged.\n")
        }
    }

    # Return statistics invisibly
    invisible(list(
        sequences = total_seqs,
        transformed = transformed,
        unchanged = unchanged
    ))
}


#' Quick transform using base R (no Biostrings dependency)
#'
#' A simpler version that uses only base R functions.
#' Works well for most files but may be slower for very large files.
#'
#' @param input_fasta Path to input FASTA file
#' @param output_fasta Path to output FASTA file
#'
preprocess_reverse_fasta_base <- function(input_fasta, output_fasta) {

    if (!file.exists(input_fasta)) {
        stop("Input file not found: ", input_fasta)
    }

    # Read all lines
    lines <- readLines(input_fasta)

    # Transform header lines only
    header_idx <- grepl("^>", lines)
    lines[header_idx] <- gsub("_paired_R_", "_paired_F_", lines[header_idx])

    # Write output
    writeLines(lines, output_fasta)

    cat("Processed", sum(header_idx), "sequences\n")
    cat("Output:", output_fasta, "\n")
}


# =============================================================================
# Command-line interface
# =============================================================================

# If run as a script (not sourced), process command line arguments
if (!interactive() && identical(environment(), globalenv())) {

    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
        cat("Usage: Rscript preprocess_reverse_fasta.R <input.fasta> <output.fasta>\n")
        cat("\n")
        cat("Example:\n")
        cat("  Rscript preprocess_reverse_fasta.R \\\n")
        cat("      fasta/16S_Bacteria_paired_R.fasta \\\n")
        cat("      fasta/16S_Bacteria_paired_R_transformed.fasta\n")
        quit(status = 0)
    }

    if (length(args) != 2) {
        cat("Error: Expected 2 arguments (input and output paths)\n")
        quit(status = 1)
    }

    input_fasta <- args[1]
    output_fasta <- args[2]

    tryCatch({
        preprocess_reverse_fasta(input_fasta, output_fasta)
    }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        quit(status = 1)
    })
}
