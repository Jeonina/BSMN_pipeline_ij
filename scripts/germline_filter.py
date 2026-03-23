#!/usr/bin/env python3
"""Germline variant filter for text-format variant files.

Removes variants found in the gnomAD SNP list.
Input:  tab-separated text (chrom  pos  ref  alt  ...) — one variant per line
Output: filtered text (same format, germline variants removed)

Usage:
    python scripts/germline_filter.py --variants gnomad.txt.gz input.txt > output.txt
    cat input.txt | python scripts/germline_filter.py -V gnomad.txt.gz > output.txt
"""

import argparse
import gzip
import sys
from typing import IO, Set


def load_germline_set(variants_gz: str) -> Set[str]:
    """Load gnomAD SNPs into a lookup set.

    Expects a gzipped file where each line is:
        chrom  pos  ref  alt   (tab-separated, no 'chr' prefix)

    Returns a set of strings formatted as "chrom:pos:ref:alt".
    """
    known: Set[str] = set()
    with gzip.open(variants_gz, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                known.add(":".join(parts[:4]))
    return known


def is_germline(chrom: str, pos: str, ref: str, alt: str, known: Set[str]) -> bool:
    """Return True if the variant is a known germline variant.

    Strips 'chr' prefix from chrom for lookup (gnomAD file uses bare chrom names).
    """
    bare_chrom = chrom[3:] if chrom.startswith("chr") else chrom
    key = f"{bare_chrom}:{pos}:{ref}:{alt}"
    return key in known


def filter_txt(infile: IO[str], known_germ: Set[str], outfile: IO[str]) -> None:
    """Filter text-format variant file, removing known germline variants.

    Lines starting with '#' are passed through unchanged.
    Lines with fewer than 4 tab-separated fields are skipped.
    """
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
        if not is_germline(chrom, pos, ref, alt, known_germ):
            outfile.write(line)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove known germline variants from text-format variant file."
    )
    parser.add_argument(
        "infile",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Input text file (chrom\\tpos\\tref\\talt per line)",
    )
    parser.add_argument(
        "--variants",
        "-V",
        required=True,
        help="gzipped gnomAD SNP file (chrom\\tpos\\tref\\talt, no chr prefix)",
    )
    args = parser.parse_args()

    known_germ = load_germline_set(args.variants)
    filter_txt(args.infile, known_germ, sys.stdout)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
