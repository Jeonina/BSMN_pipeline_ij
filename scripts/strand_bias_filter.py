#!/usr/bin/env python3
"""Strand-bias test per SNV.

Consolidated from ``utils/strand_bias.2.py`` per SPEC-BSMN-REFACTOR-001
M2 + REQ-DUP-002. The legacy ``load_config(reference, conda_env)`` hook
has been stripped; the SAMTOOLS path is provided by the Apptainer
container on PATH. The ``--reference`` and ``--conda-env`` CLI flags are
removed.

For each SNV, emit forward/reverse strand counts for total, ref, and alt
alleles, a Poisson p-value on overall strand balance (R via rpy2), and a
Fisher exact p-value on the 2x2 fwd/rev × ref/alt table (scipy).
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from collections.abc import Generator
from multiprocessing import Pool
from typing import Any

from rpy2.robjects import r
from scipy.stats import fisher_exact

pipe_home = os.path.dirname(os.path.realpath(__file__)) + "/.."
sys.path.append(pipe_home)
from bsmn_pipeline.misc import coroutine, printer  # noqa: E402
from bsmn_pipeline.pileup import base_count  # noqa: E402


def run(args: argparse.Namespace) -> None:
    header = (
        "#chr\tpos\tref\talt\t"
        "total\ttotal_fwd\ttotal_rev\ttotal_ratio\t"
        "p_poisson\t"
        "ref_n\tref_fwd\tref_rev\tref_ratio\t"
        "alt_n\talt_fwd\talt_rev\talt_ratio\t"
        "p_fisher"
    )
    printer(header)
    if args.nproc > 1:
        with Pool(args.nproc) as p:
            for r_ in p.starmap(
                _mpileup,
                [
                    [args.bam, args.min_MQ, args.min_BQ] + snv.strip().split()[:4]
                    for snv in args.infile
                    if snv[0] != "#"
                ],
            ):
                printer(r_)
    else:
        s_info = strand_info(base_count(args.bam, args.min_MQ, args.min_BQ))
        for snv in args.infile:
            if snv[0] == "#":
                continue
            chrom, pos, ref, alt = snv.strip().split()[:4]
            printer(mpileup(s_info, chrom, pos, ref, alt))
    sys.stdout.flush()


def _mpileup(bam: str, min_MQ: int, min_BQ: int, chrom: str, pos: str, ref: str, alt: str) -> str:
    return mpileup(strand_info(base_count(bam, min_MQ, min_BQ)), chrom, pos, ref, alt)


def mpileup(
    s_info: Generator[str, tuple[str, str, str, str]],
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
) -> str:
    return f"{chrom}\t{pos}\t{ref.upper()}\t{alt.upper()}\t{s_info.send((chrom, pos, ref, alt))}"


@coroutine
def strand_info(
    target: Generator[Any, Any, Any],
) -> Generator[str, tuple[str, str, str, str]]:
    result: str | None = None
    while True:
        chrom, pos, ref, alt = yield result  # type: ignore[misc]
        base_n = target.send((chrom, pos))
        total = sum(base_n.values())
        total_fwd = sum(list(base_n.values())[:4])
        total_rev = sum(list(base_n.values())[4:8])
        try:
            total_ratio = total_fwd / total_rev
        except ZeroDivisionError:
            total_ratio = math.inf
        ref_n = base_n[ref.upper()] + base_n[ref.lower()]
        ref_fwd = base_n[ref.upper()]
        ref_rev = base_n[ref.lower()]
        try:
            ref_ratio = ref_fwd / ref_rev
        except ZeroDivisionError:
            ref_ratio = math.inf
        alt_n = base_n[alt.upper()] + base_n[alt.lower()]
        alt_fwd = base_n[alt.upper()]
        alt_rev = base_n[alt.lower()]
        try:
            alt_ratio = alt_fwd / alt_rev
        except ZeroDivisionError:
            alt_ratio = math.inf

        result = (
            f"{total}\t{total_fwd}\t{total_rev}\t{total_ratio:f}\t"
            f"{p_poisson(total_fwd, total_rev):e}\t"
            f"{ref_n}\t{ref_fwd}\t{ref_rev}\t{ref_ratio:f}\t"
            f"{alt_n}\t{alt_fwd}\t{alt_rev}\t{alt_ratio:f}\t"
            f"{p_fisher(ref_fwd, alt_fwd, ref_rev, alt_rev):e}"
        )


def p_poisson(n_fwd: int, n_rev: int) -> float:
    return float(r(f"poisson.test(c({n_fwd},{n_rev}))$p.value")[0])


def p_fisher(ref_fwd: int, alt_fwd: int, ref_rev: int, alt_rev: int) -> float:
    return float(fisher_exact([[ref_fwd, alt_fwd], [ref_rev, alt_rev]])[1])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check strand bias for SNV")
    parser.add_argument("-b", "--bam", metavar="FILE", help="bam file", required=True)
    parser.add_argument(
        "-q",
        "--min-MQ",
        metavar="INT",
        help="mapQ cutoff value [20]",
        type=int,
        default=20,
    )
    parser.add_argument(
        "-Q",
        "--min-BQ",
        metavar="INT",
        help="baseQ/BAQ cutoff value [13]",
        type=int,
        default=13,
    )
    parser.add_argument(
        "-n",
        "--nproc",
        metavar="INT",
        help="Specifies the number of processors to use [default: 1]",
        type=int,
        default=1,
    )
    parser.add_argument(
        "infile",
        metavar="snv_list.txt",
        help=(
            "SNV list. Each line format is 'chr\\tpos\\tref\\talt'. "
            "Trailing columns will be ignored. [STDIN]"
        ),
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
    )
    parser.set_defaults(func=run)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if len(vars(args)) == 0:
        parser.print_help()
    else:
        args.func(args)


if __name__ == "__main__":
    main()
