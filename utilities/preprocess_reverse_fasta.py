#!/usr/bin/env python3
"""
preprocess_reverse_fasta.py
Transform reverse FASTA headers to match BIOM feature IDs.

eDNA Explorer BIOM files use forward read names as feature IDs:
    16S_Bacteria_paired_F_0

Forward FASTA headers match directly:
    >16S_Bacteria_paired_F_0

Reverse FASTA headers use _paired_R_:
    >16S_Bacteria_paired_R_0

This script transforms reverse FASTA headers to match BIOM feature IDs
by replacing _paired_R_ with _paired_F_.

Usage:
    python preprocess_reverse_fasta.py <input.fasta> <output.fasta>

Example:
    python preprocess_reverse_fasta.py \\
        fasta/16S_Bacteria_paired_R.fasta \\
        fasta/16S_Bacteria_paired_R_transformed.fasta

Dependencies:
    None (Python 3 standard library only)
"""

import sys
from pathlib import Path


def preprocess_reverse_fasta(
    input_path: str, output_path: str, verbose: bool = True
) -> dict:
    """
    Transform reverse FASTA headers to match BIOM feature IDs.

    Args:
        input_path: Path to input reverse FASTA file
        output_path: Path for output transformed FASTA file
        verbose: Print progress information

    Returns:
        Dictionary with statistics:
        - sequences: Total number of sequences
        - transformed: Number of headers transformed
        - unchanged: Number of headers unchanged (no _paired_R_)

    Raises:
        FileNotFoundError: If input file doesn't exist
    """
    input_file = Path(input_path)
    output_file = Path(output_path)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)

    stats = {"sequences": 0, "transformed": 0, "unchanged": 0}

    with open(input_file) as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                stats["sequences"] += 1
                original = line

                # Transform header: _paired_R_ -> _paired_F_
                transformed_line = line.replace("_paired_R_", "_paired_F_")

                if transformed_line != original:
                    stats["transformed"] += 1
                else:
                    stats["unchanged"] += 1

                outfile.write(transformed_line)
            else:
                # Sequence line - write unchanged
                outfile.write(line)

    if verbose:
        print(f"Processed {stats['sequences']} sequences")
        print(f"  Transformed: {stats['transformed']}")
        print(f"  Unchanged: {stats['unchanged']}")

    return stats


def main():
    """Main entry point."""
    # Parse arguments
    if len(sys.argv) == 2 and sys.argv[1] in ["-h", "--help"]:
        print(__doc__)
        sys.exit(0)

    if len(sys.argv) != 3:
        print("Error: Expected 2 arguments")
        print()
        print("Usage: python preprocess_reverse_fasta.py <input.fasta> <output.fasta>")
        print()
        print("Example:")
        print("  python preprocess_reverse_fasta.py \\")
        print("      fasta/16S_Bacteria_paired_R.fasta \\")
        print("      fasta/16S_Bacteria_paired_R_transformed.fasta")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    print("=" * 60)
    print("  Reverse FASTA Preprocessor")
    print("=" * 60)
    print()
    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")
    print()

    try:
        # Show example transformation from first header
        with open(input_path) as f:
            for line in f:
                if line.startswith(">"):
                    original = line.strip()
                    transformed = original.replace("_paired_R_", "_paired_F_")
                    print("Transformation:")
                    print(f"  Before: {original}")
                    print(f"  After:  {transformed}")
                    print()
                    break

        # Process file
        stats = preprocess_reverse_fasta(input_path, output_path)

        print()
        print(f"Output written to: {output_path}")

        if stats["unchanged"] > 0:
            print()
            print(
                f"Note: {stats['unchanged']} sequence(s) did not contain '_paired_R_'"
            )
            print("      These headers were written unchanged.")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except PermissionError as e:
        print(f"Error: Permission denied - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
