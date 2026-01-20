#!/usr/bin/env python3
"""
transform_reverse_fasta.py
Preprocess reverse FASTA for QIIME2 import.

Transforms FASTA headers from _paired_R_ to _paired_F_ format
to match BIOM feature IDs.

Usage:
    python transform_reverse_fasta.py <input.fasta> <output.fasta>

Example:
    python transform_reverse_fasta.py \
        fasta/16S_Bacteria_paired_R.fasta \
        fasta/16S_Bacteria_paired_R_transformed.fasta
"""

import sys
from pathlib import Path


def transform_reverse_fasta(input_path: str, output_path: str) -> int:
    """
    Transform reverse FASTA headers to match BIOM feature IDs.

    Changes _paired_R_ to _paired_F_ in sequence headers.

    Args:
        input_path: Path to input reverse FASTA file
        output_path: Path for output transformed FASTA file

    Returns:
        Number of sequences processed
    """
    input_file = Path(input_path)
    output_file = Path(output_path)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)

    seq_count = 0
    transformed_count = 0

    with open(input_file) as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                seq_count += 1
                original = line
                # Transform header: _paired_R_ -> _paired_F_
                line = line.replace("_paired_R_", "_paired_F_")
                if line != original:
                    transformed_count += 1
            outfile.write(line)

    return seq_count


def main():
    """Main entry point."""
    if len(sys.argv) != 3:
        print(__doc__)
        print("Error: Expected 2 arguments")
        print(f"Usage: {sys.argv[0]} <input.fasta> <output.fasta>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")
    print()

    try:
        # Show example transformation
        with open(input_path) as f:
            first_line = f.readline().strip()
            if first_line.startswith(">"):
                transformed = first_line.replace("_paired_R_", "_paired_F_")
                print("Header transformation:")
                print(f"  Before: {first_line}")
                print(f"  After:  {transformed}")
                print()

        seq_count = transform_reverse_fasta(input_path, output_path)
        print(f"Transformed {seq_count} sequences")
        print(f"Output written to: {output_path}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
