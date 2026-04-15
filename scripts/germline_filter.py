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
import logging
import os
import sys
import time
from typing import IO, Set

log = logging.getLogger("germline_filter")


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


def filter_txt(infile: IO[str], known_germ: Set[str], outfile: IO[str]) -> dict:
    """Filter text-format variant file, removing known germline variants.

    Lines starting with '#' are passed through unchanged.
    Lines with fewer than 4 tab-separated fields are skipped.
    Returns a dict with filtering statistics.
    """
    stats = {"input": 0, "kept": 0, "removed": 0, "skipped": 0, "comments": 0}
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            stats["comments"] += 1
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            stats["skipped"] += 1
            continue
        stats["input"] += 1
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
        if is_germline(chrom, pos, ref, alt, known_germ):
            stats["removed"] += 1
            log.debug(
                "REMOVED germline: %s:%s %s>%s", chrom, pos, ref, alt
            )
        else:
            stats["kept"] += 1
            outfile.write(line)
    return stats


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

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    log.info("================================================================")
    log.info("START germline_filter")
    log.info("input_file=%s", args.infile.name)
    log.info("gnomad_variants=%s", args.variants)
    log.info("gnomad_file_size=%s bytes",
             os.path.getsize(args.variants) if os.path.exists(args.variants) else "N/A")
    log.info("criterion: remove variants present in gnomAD (AF > 0.001)")
    log.info("================================================================")

    t0 = time.time()
    log.info("Loading gnomAD variant set...")
    known_germ = load_germline_set(args.variants)
    t_load = time.time() - t0
    log.info("gnomAD variants loaded: %d entries in %.1f seconds", len(known_germ), t_load)

    t1 = time.time()
    stats = filter_txt(args.infile, known_germ, sys.stdout)
    sys.stdout.flush()
    t_filter = time.time() - t1

    log.info("================================================================")
    log.info("RESULTS:")
    log.info("  variants_input=%d", stats["input"])
    log.info("  variants_kept=%d", stats["kept"])
    log.info("  variants_removed=%d (germline)", stats["removed"])
    log.info("  lines_skipped=%d (malformed)", stats["skipped"])
    log.info("  comment_lines=%d", stats["comments"])
    if stats["input"] > 0:
        log.info("  pass_rate=%.1f%%", stats["kept"] / stats["input"] * 100)
    log.info("  filter_time=%.1f seconds", t_filter)
    log.info("  total_time=%.1f seconds", time.time() - t0)
    log.info("END germline_filter")
    log.info("================================================================")


if __name__ == "__main__":
    main()
