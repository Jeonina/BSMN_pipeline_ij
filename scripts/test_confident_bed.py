#!/usr/bin/env python3
"""Restrict mosaic candidates to a benchmark confident BED.

GIAB-style benchmarks only score calls inside the truth set's confident
regions. A candidate OUTSIDE the confident BED is not a false positive — it
is a no-call (no truth assertion there). This tool splits the final mosaic
candidate list into in-BED (scored) vs outside-BED (excluded) so the true FP
denominator can be established before any further filtering.

A candidate at 1-based ``pos`` is in-BED when some BED interval (0-based,
half-open ``start end``) satisfies ``start < pos <= end``.

Usage:
    python scripts/test_confident_bed.py \
        --candidates reports/benchmark_2026-06-02/HG002.final.txt \
        --bed resources/hg38/HG002_mosaic_confident.bed \
        --truth chr20:47053475
"""

from __future__ import annotations

import argparse
import gzip
from collections import defaultdict


def read_candidates(path: str) -> list[tuple[str, int, str, str]]:
    out: list[tuple[str, int, str, str]] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            out.append((p[0], int(p[1]), p[2].upper(), p[3].upper()))
    return out


def read_bed(path: str) -> dict[str, list[tuple[int, int]]]:
    """Load BED intervals per chrom as sorted (start, end), 0-based half-open."""
    opener = gzip.open if path.endswith(".gz") else open
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            by_chrom[f[0]].append((int(f[1]), int(f[2])))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def in_bed(by_chrom: dict[str, list[tuple[int, int]]], chrom: str, pos: int) -> bool:
    """True if 1-based pos falls within any 0-based half-open interval."""
    intervals = by_chrom.get(chrom)
    if not intervals:
        return False
    # Linear scan is fine for a handful of candidates.
    return any(start < pos <= end for start, end in intervals)


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--candidates", required=True, help="final mosaic list (chr pos ref alt)")
    p.add_argument("--bed", required=True, help="confident regions BED(.gz)")
    p.add_argument("--truth", default="", help="known true-positive coord, e.g. chr20:47053475")
    args = p.parse_args()

    cands = read_candidates(args.candidates)
    bed = read_bed(args.bed)
    truth_key: tuple[str, int] | None = None
    if args.truth:
        tc, tp = args.truth.split(":")
        truth_key = (tc, int(tp))

    inside: list[tuple[str, int, str, str]] = []
    outside: list[tuple[str, int, str, str]] = []
    for c, pos, ref, alt in cands:
        (inside if in_bed(bed, c, pos) else outside).append((c, pos, ref, alt))

    def tag(c: str, pos: int) -> str:
        return "  <-- TP" if (c, pos) == truth_key else ""

    print(f"Confident-BED restriction  (candidates = {len(cands)})")
    print(f"  in-BED  (scored)   : {len(inside)}")
    print(f"  outside (excluded) : {len(outside)}")

    tp_inside = truth_key in {(c, pos) for c, pos, _, _ in inside} if truth_key else None
    if truth_key is not None:
        loc = "in-BED" if tp_inside else "OUTSIDE (truth should be in-BED — check BED/coords)"
        print(f"  TP location        : {loc}")

    n_fp_inbed = sum(1 for c, pos, _, _ in inside if (c, pos) != truth_key)
    print(f"\n  => true FP denominator (in-BED, non-truth) = {n_fp_inbed}")

    print("\nIn-BED candidates (scored):")
    for c, pos, ref, alt in inside:
        print(f"  {c}:{pos} {ref}>{alt}{tag(c, pos)}")
    print("\nOutside-BED candidates (excluded from FP count):")
    for c, pos, ref, alt in outside:
        print(f"  {c}:{pos} {ref}>{alt}{tag(c, pos)}")


if __name__ == "__main__":
    main()
