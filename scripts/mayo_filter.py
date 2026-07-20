#!/usr/bin/env python3
"""Mayo filter: BSMN E-step strand-bias + repeat + multiallelic filters.

Faithful re-implementation of ``jobs/variant_filtering/E.mayo_filters.sh`` as a
pure join over the outputs of ``scripts/strand_bias_filter.py`` and
``scripts/repeat_filter.py``. A candidate is kept only when ALL hold:

  repeat:         repeat_n < 4 AND repeat_length < 10    (not inside an STR)
  multiallelic:   total == ref_n + alt_n                 (no third allele)
  both strands:   alt_fwd >= 1 AND alt_rev >= 1          (alt on both strands)
  strand balance: p_poisson >= 0.05 OR p_fisher >= 0.05  (not strand-biased)

The strand/repeat computations stay in the original BSMN scripts (which shell
out to samtools); this script performs no BAM/CRAM I/O, so its decision logic
is unit-testable without alignment data. The orchestrating Snakemake rule runs
the two upstream scripts and feeds their TSVs here.

Usage:
    python scripts/mayo_filter.py \
        --strand results/filtering/SM/SM.strand.tsv \
        --repeat results/filtering/SM/SM.repeat.tsv \
        > results/filtering/SM/SM.mayo_filtered.txt
"""

from __future__ import annotations

import argparse
import logging
import sys

log = logging.getLogger("mayo_filter")

Key = tuple[str, str, str, str]

# Columns cast to int / float when reading the upstream TSVs; everything else
# (e.g. repeat_seq, *_ratio) stays a string.
_INT_COLS = {
    "total",
    "total_fwd",
    "total_rev",
    "ref_n",
    "ref_fwd",
    "ref_rev",
    "alt_n",
    "alt_fwd",
    "alt_rev",
    "repeat_n",
    "repeat_length",
}
_FLOAT_COLS = {"p_poisson", "p_fisher"}


def read_tsv(path: str) -> dict[Key, dict[str, object]]:
    """Read a headered TSV keyed by (chr, pos, ref, alt). Header may start '#'."""
    rows: dict[Key, dict[str, object]] = {}
    with open(path) as fh:
        header_line = fh.readline()
        if not header_line:
            return rows
        cols = header_line.lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            rec: dict[str, object] = {}
            for col, val in zip(cols, vals, strict=False):
                if col in _INT_COLS:
                    rec[col] = int(val)
                elif col in _FLOAT_COLS:
                    rec[col] = float(val)
                else:
                    rec[col] = val
            key = (rec["chr"], rec["pos"], rec["ref"], rec["alt"])  # type: ignore[index]
            rows[key] = rec
    return rows


def passes_mayo(
    s: dict[str, object],
    r: dict[str, object],
    max_repeat_n: int,
    max_repeat_length: int,
    min_strand_p: float,
) -> tuple[bool, str | None]:
    """Return (kept, reject_reason). reject_reason is None when kept.

    Mirrors the four awk filters in E.mayo_filters.sh, evaluated in the same
    order so the first failing criterion is reported.
    """
    repeat_n = r["repeat_n"]
    repeat_length = r["repeat_length"]
    if not (repeat_n < max_repeat_n and repeat_length < max_repeat_length):  # type: ignore[operator]
        return False, "repeat"
    if s["total"] != s["ref_n"] + s["alt_n"]:  # type: ignore[operator]
        return False, "multiallelic"
    if not (s["alt_fwd"] >= 1 and s["alt_rev"] >= 1):  # type: ignore[operator]
        return False, "single_strand"
    if not (s["p_poisson"] >= min_strand_p or s["p_fisher"] >= min_strand_p):  # type: ignore[operator]
        return False, "strand_bias"
    return True, None


def run(args: argparse.Namespace) -> None:
    strand = read_tsv(args.strand)
    repeat = read_tsv(args.repeat)

    log.info("strand records=%d  repeat records=%d", len(strand), len(repeat))

    counts = {
        "kept": 0,
        "repeat": 0,
        "multiallelic": 0,
        "single_strand": 0,
        "strand_bias": 0,
        "missing": 0,
    }
    for key, s in strand.items():
        r = repeat.get(key)
        if r is None:
            counts["missing"] += 1
            log.warning("no repeat record for %s — dropping", ":".join(key))
            continue
        kept, reason = passes_mayo(
            s, r, args.max_repeat_n, args.max_repeat_length, args.min_strand_p
        )
        if kept:
            counts["kept"] += 1
            sys.stdout.write("\t".join(key) + "\n")
        else:
            counts[reason] += 1  # type: ignore[index]

    sys.stdout.flush()
    total_in = len(strand)
    log.info(
        "mayo: %d -> %d  (repeat=%d multiallelic=%d single_strand=%d strand_bias=%d missing=%d)",
        total_in,
        counts["kept"],
        counts["repeat"],
        counts["multiallelic"],
        counts["single_strand"],
        counts["strand_bias"],
        counts["missing"],
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="BSMN mayo filter (strand bias + repeat + multiallelic).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--strand", required=True, help="strand_bias_filter.py TSV")
    p.add_argument("--repeat", required=True, help="repeat_filter.py TSV")
    p.add_argument("--max-repeat-n", type=int, default=4, help="reject if repeat_n >= this [4]")
    p.add_argument(
        "--max-repeat-length", type=int, default=10, help="reject if repeat_length >= this [10]"
    )
    p.add_argument(
        "--min-strand-p",
        type=float,
        default=0.05,
        help="keep if p_poisson OR p_fisher >= this [0.05]",
    )
    p.set_defaults(func=run)
    return p


def main() -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(name)s] %(levelname)s: %(message)s"
    )
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
