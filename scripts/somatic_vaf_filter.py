#!/usr/bin/env python3
"""Somatic VAF (Variant Allele Frequency) per-SNV test.

Consolidated from ``utils/somatic_vaf.2.py`` per SPEC-BSMN-REFACTOR-001
M2 + REQ-DUP-002. The legacy ``load_config(reference, conda_env)`` hook
that switched the SAMTOOLS path via a conda env has been stripped — the
Snakemake workflow now provides SAMTOOLS through an Apptainer container
on PATH, and the resolved reference is passed via
``config/resolved_params.yaml``. The ``--reference`` and ``--conda-env``
CLI flags are removed accordingly.

For each SNV in the input list (chr / pos / ref / alt, tab-separated),
emit the VAF, depth, ref/alt allele counts, and the binomial p-value
that ``alt_n`` is observed under H0 (germline).
"""

from __future__ import annotations

import argparse
import os
import sys
from multiprocessing import Pool

from statsmodels.stats.proportion import binom_test

pipe_home = os.path.dirname(os.path.realpath(__file__)) + "/.."
sys.path.append(pipe_home)
from bsmn_pipeline.misc import coroutine, printer  # noqa: E402
from bsmn_pipeline.pileup import base_count  # noqa: E402


def run(args: argparse.Namespace) -> None:
    header = (
        "#chr\tpos\tref\talt\tvaf\t"
        "depth\tref_n\talt_n\tp_binom"
    )
    printer(header)
    if args.nproc > 1:
        with Pool(args.nproc) as p:
            for r in p.starmap(
                _mpileup,
                [
                    [args.bam, args.min_MQ, args.min_BQ]
                    + snv.strip().split()[:4]
                    for snv in args.infile
                    if snv[0] != "#"
                ],
            ):
                printer(r)
    else:
        v_info = vaf_info(base_count(args.bam, args.min_MQ, args.min_BQ))
        for snv in args.infile:
            if snv[0] == "#":
                continue
            chrom, pos, ref, alt = snv.strip().split()[:4]
            printer(mpileup(v_info, chrom, pos, ref, alt))
    sys.stdout.flush()


def _mpileup(bam, min_MQ, min_BQ, chrom, pos, ref, alt):
    return mpileup(vaf_info(base_count(bam, min_MQ, min_BQ)), chrom, pos, ref, alt)


def mpileup(v_info, chrom, pos, ref, alt):
    return "{chrom}\t{pos}\t{ref}\t{alt}\t{vaf_info}".format(
        chrom=chrom,
        pos=pos,
        ref=ref.upper(),
        alt=alt.upper(),
        vaf_info=v_info.send((chrom, pos, ref, alt)),
    )


@coroutine
def vaf_info(target):
    result = None
    while True:
        chrom, pos, ref, alt = yield result
        base_n = target.send((chrom, pos))
        depth = sum(base_n.values())
        ref_n = base_n[ref.upper()] + base_n[ref.lower()]
        alt_n = base_n[alt.upper()] + base_n[alt.lower()]
        try:
            vaf = alt_n / depth
        except ZeroDivisionError:
            vaf = 0

        result = "{vaf:f}\t{depth}\t{ref_n}\t{alt_n}\t{p_binom:e}".format(
            vaf=vaf,
            depth=depth,
            ref_n=ref_n,
            alt_n=alt_n,
            p_binom=binom_test(alt_n, depth, alternative="smaller"),
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Test whether VAF of each SNV is somatic or germline."
    )
    parser.add_argument(
        "-b", "--bam", metavar="FILE", help="bam file", required=True
    )
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
