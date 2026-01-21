#!/usr/bin/env python3
"""
Fast ASV lookup utility for eDNA Explorer v2.1 bundles.

Uses pyfaidx for O(1) indexed FASTA lookups (no samtools required).

Usage:
    # Single ASV lookup
    python lookup.py sequence 16S_Bacteria-paired_f.fasta 16S_Bacteria_paired_F_0

    # Multiple ASVs
    python lookup.py sequence 16S_Bacteria-paired_f.fasta ASV1 ASV2 ASV3

    # From file (one ID per line)
    python lookup.py sequence 16S_Bacteria-paired_f.fasta -f ids.txt

    # Taxonomy lookup
    python lookup.py taxonomy 16S_Bacteria-paired-lookup.tsv 16S_Bacteria_paired_F_819859

    # Filter by confidence
    python lookup.py filter 16S_Bacteria-paired-lookup.tsv --min-confidence 3

    # Full info (sequence + taxonomy)
    python lookup.py info 16S_Bacteria-paired_f.fasta 16S_Bacteria-paired-lookup.tsv 16S_Bacteria_paired_F_0
"""

import argparse
import sys
import time


def get_sequence(fasta_path: str, asv_ids: list[str], from_file: str | None = None) -> None:
    """Get sequences for ASV IDs using pyfaidx (O(1) indexed lookup)."""
    try:
        from pyfaidx import Fasta
    except ImportError:
        print("Error: pyfaidx not installed. Run: pip install pyfaidx", file=sys.stderr)
        sys.exit(1)

    # Load IDs from file if specified
    if from_file:
        with open(from_file) as f:
            asv_ids = [line.strip() for line in f if line.strip()]

    if not asv_ids:
        print("Error: No ASV IDs provided", file=sys.stderr)
        sys.exit(1)

    # Open indexed FASTA (creates index on first access)
    fasta = Fasta(fasta_path)

    for asv_id in asv_ids:
        if asv_id in fasta:
            seq = fasta[asv_id]
            print(f">{seq.name}")
            print(str(seq))
        else:
            print(f"# {asv_id} not found", file=sys.stderr)


def get_taxonomy(lookup_path: str, asv_ids: list[str]) -> None:
    """Get taxonomy for ASV IDs from lookup.tsv."""
    import pandas as pd

    lookup = pd.read_csv(lookup_path, sep="\t")
    lookup_dict = dict(zip(lookup["feature_id"], lookup["taxonomy"]))
    conf_dict = dict(zip(lookup["feature_id"], lookup["confidence"]))

    for asv_id in asv_ids:
        if asv_id in lookup_dict:
            print(f"{asv_id}\t{lookup_dict[asv_id]}\t{conf_dict[asv_id]}")
        else:
            print(f"# {asv_id} not found", file=sys.stderr)


def filter_lookup(
    lookup_path: str, min_confidence: int = 0, taxon: str | None = None, output: str | None = None
) -> None:
    """Filter lookup.tsv by confidence or taxonomy."""
    import pandas as pd

    lookup = pd.read_csv(lookup_path, sep="\t")

    if min_confidence > 0:
        lookup = lookup[lookup["confidence"] >= min_confidence]

    if taxon:
        lookup = lookup[lookup["taxonomy"].str.contains(taxon, case=False, na=False)]

    print(f"# Filtered: {len(lookup):,} ASVs", file=sys.stderr)

    if output:
        lookup.to_csv(output, sep="\t", index=False)
        print(f"# Saved to: {output}", file=sys.stderr)
    else:
        # Print just IDs to stdout for piping
        for fid in lookup["feature_id"]:
            print(fid)


def get_full_info(fasta_path: str, lookup_path: str, asv_ids: list[str]) -> None:
    """Get full info (sequence + taxonomy) for ASV IDs."""
    try:
        from pyfaidx import Fasta
    except ImportError:
        print("Error: pyfaidx not installed. Run: pip install pyfaidx", file=sys.stderr)
        sys.exit(1)

    import pandas as pd

    fasta = Fasta(fasta_path)
    lookup = pd.read_csv(lookup_path, sep="\t")
    lookup_indexed = lookup.set_index("feature_id")

    for asv_id in asv_ids:
        print(f"=== {asv_id} ===")

        # Sequence
        if asv_id in fasta:
            seq = str(fasta[asv_id])
            print(f"Sequence ({len(seq)} bp): {seq[:60]}...")
        else:
            print("Sequence: NOT FOUND")

        # Taxonomy
        if asv_id in lookup_indexed.index:
            row = lookup_indexed.loc[asv_id]
            print(f"Taxonomy: {row['taxonomy']}")
            print(f"Confidence: {row['confidence']}")
            print(f"Read set: {row['read_set']}")
            print(f"Length: {row['sequence_length']} bp")
        else:
            print("Taxonomy: NOT FOUND")

        print()


def benchmark(fasta_path: str, n: int = 100) -> None:
    """Benchmark indexed vs linear lookup performance."""
    try:
        from pyfaidx import Fasta
    except ImportError:
        print("Error: pyfaidx not installed", file=sys.stderr)
        sys.exit(1)

    import random

    # Load index
    print(f"Loading FASTA index from {fasta_path}...")
    start = time.time()
    fasta = Fasta(fasta_path)
    index_time = time.time() - start
    print(f"Index load time: {index_time:.2f}s")
    print(f"Total sequences: {len(fasta.keys()):,}")

    # Get random IDs
    all_ids = list(fasta.keys())
    test_ids = random.sample(all_ids, min(n, len(all_ids)))

    # Benchmark indexed lookups
    start = time.time()
    for asv_id in test_ids:
        _ = str(fasta[asv_id])
    indexed_time = time.time() - start

    print(f"\nIndexed lookup ({n} sequences): {indexed_time:.4f}s")
    print(f"Average per lookup: {indexed_time / n * 1000:.4f}ms")


def main():
    parser = argparse.ArgumentParser(
        description="Fast ASV lookup utility for eDNA Explorer v2.1 bundles",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # sequence command
    seq_parser = subparsers.add_parser("sequence", help="Get sequences for ASV IDs")
    seq_parser.add_argument("fasta", help="FASTA file path")
    seq_parser.add_argument("ids", nargs="*", help="ASV IDs to look up")
    seq_parser.add_argument("-f", "--file", help="File with ASV IDs (one per line)")

    # taxonomy command
    tax_parser = subparsers.add_parser("taxonomy", help="Get taxonomy for ASV IDs")
    tax_parser.add_argument("lookup", help="lookup.tsv file path")
    tax_parser.add_argument("ids", nargs="+", help="ASV IDs to look up")

    # filter command
    filter_parser = subparsers.add_parser("filter", help="Filter lookup.tsv")
    filter_parser.add_argument("lookup", help="lookup.tsv file path")
    filter_parser.add_argument("--min-confidence", "-c", type=int, default=0, help="Min confidence")
    filter_parser.add_argument("--taxon", "-t", help="Filter by taxon (substring match)")
    filter_parser.add_argument("--output", "-o", help="Output file (default: stdout)")

    # info command
    info_parser = subparsers.add_parser("info", help="Get full info for ASV IDs")
    info_parser.add_argument("fasta", help="FASTA file path")
    info_parser.add_argument("lookup", help="lookup.tsv file path")
    info_parser.add_argument("ids", nargs="+", help="ASV IDs to look up")

    # benchmark command
    bench_parser = subparsers.add_parser("benchmark", help="Benchmark lookup performance")
    bench_parser.add_argument("fasta", help="FASTA file path")
    bench_parser.add_argument("-n", type=int, default=100, help="Number of lookups")

    args = parser.parse_args()

    if args.command == "sequence":
        get_sequence(args.fasta, args.ids, args.file)
    elif args.command == "taxonomy":
        get_taxonomy(args.lookup, args.ids)
    elif args.command == "filter":
        filter_lookup(args.lookup, args.min_confidence, args.taxon, args.output)
    elif args.command == "info":
        get_full_info(args.fasta, args.lookup, args.ids)
    elif args.command == "benchmark":
        benchmark(args.fasta, args.n)


if __name__ == "__main__":
    main()
