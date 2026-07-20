#!/usr/bin/env python3
"""Sum of base qualities of the alt allele at each SNV.

Consolidated from ``utils/alt_bq_sum.py`` per SPEC-BSMN-REFACTOR-001 M2 +
REQ-DUP-002 / DEF-002. The legacy ``load_config(reference, conda_env)``
hook is stripped — the SAMTOOLS path is provided by the Apptainer
container on PATH. The ``--reference`` / ``--conda-env`` CLI flags are
removed accordingly.

For each SNV (chr / pos / ref / alt), emit the count and the sum of
base qualities of reads carrying the alt allele.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Generator
from typing import Any

pipe_home = os.path.dirname(os.path.realpath(__file__)) + "/.."
sys.path.append(pipe_home)
from bsmn_pipeline.misc import coroutine, printer  # noqa: E402
from bsmn_pipeline.pileup import base_qual_tuple  # noqa: E402


def run(args: argparse.Namespace) -> None:
    alt_BQ_info = alt_BQ_sum(base_qual_tuple(args.bam, args.min_MQ, args.min_BQ))
    header = "#chr\tpos\tref\talt\talt_n\talt_BQ_sum"
    printer(header)
    for snv in args.infile:
        if snv[0] == "#":
            continue
        chrom, pos, ref, alt = snv.strip().split()[:4]
        printer(
            f"{chrom}\t{pos}\t{ref.upper()}\t{alt.upper()}\t{alt_BQ_info.send((chrom, pos, alt))}"
        )


@coroutine
def alt_BQ_sum(
    target: Generator[Any, Any, Any],
) -> Generator[str, tuple[str, str, str]]:
    result: str | None = None
    while True:
        chrom, pos, alt = yield result  # type: ignore[misc]
        alt_BQ = [q for b, q in target.send((chrom, pos)) if b == alt.upper()]
        result = f"{len(alt_BQ)}\t{sum(alt_BQ)}"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Sum of base qualities of alt allele of each SNV")
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
