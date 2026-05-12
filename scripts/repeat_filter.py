#!/usr/bin/env python3
"""STR (short tandem repeat) status around each SNV.

Consolidated from ``utils/repeat.2.py`` per SPEC-BSMN-REFACTOR-001 M2.
No ``load_config`` dependency to strip — this winner only shells out to
``samtools faidx`` on PATH (provided by the Apptainer container).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from multiprocessing import Pool

# Set by ``run`` before the worker pool fans out (module-level state
# is required because ``Pool.starmap`` workers fork the parent process).
ref_file: str | None = None


def ref_seq(chrom: str, pos1: int, pos2: int | None = None) -> str:
    if pos2 is None:
        site = f"{chrom}:{pos1}-{pos1}"
    else:
        site = f"{chrom}:{pos1}-{pos2}"
    if ref_file is None:
        raise RuntimeError("ref_file is not set; call run() with args.ref first")
    base = "".join(
        subprocess.run(
            ["samtools", "faidx", ref_file, site],
            stdout=subprocess.PIPE,
            encoding="utf-8",
            check=False,
        ).stdout.split("\n")[1:]
    )
    return base


def repeat(chrom: str, pos: str, alt: str) -> str:
    read_size = 100
    read = ref_seq(chrom, int(pos) - read_size, int(pos) + read_size)
    alt_p = read_size
    w_max = 5
    n_max = 0
    repeat_out = ""
    for wsize in range(1, w_max + 1):
        for i in range(wsize):
            start = alt_p - wsize + 1 + i
            end = start + wsize

            word = read[start:alt_p] + alt + read[alt_p + 1 : end]

            n = 1
            while read[start - wsize : start] == word:
                start -= wsize
                n += 1

            while read[end : end + wsize] == word:
                end += wsize
                n += 1

            repeat_seq = "{}[{}>{}]{}".format(
                read[start:alt_p], read[alt_p], alt, read[alt_p + 1 : end]
            )

            if n > n_max:
                repeat_out = f"{n}\t{end - start}\t{repeat_seq}"
                n_max = n

    return repeat_out


def run(args: argparse.Namespace) -> None:
    global ref_file
    ref_file = args.ref

    header = "#chr\tpos\tref\talt\trepeat_n\trepeat_length\trepeat_seq"
    print(header, flush=True)
    if args.nproc > 1:
        with Pool(args.nproc) as p:
            for r in p.starmap(
                faidx,
                [
                    snv.strip().split()[:4]
                    for snv in args.infile
                    if snv[0] != "#"
                ],
            ):
                print(r, flush=True)
    else:
        for snv in args.infile:
            if snv[0] == "#":
                continue
            chrom, pos, ref, alt = snv.strip().split()[:4]
            print(faidx(chrom, pos, ref, alt), flush=True)
    sys.stdout.flush()


def faidx(chrom: str, pos: str, ref: str, alt: str) -> str:
    return "{chrom}\t{pos}\t{ref}\t{alt}\t{repeat}".format(
        chrom=chrom,
        pos=pos,
        ref=ref.upper(),
        alt=alt.upper(),
        repeat=repeat(chrom, pos, alt.upper()),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="STR status around each SNV.")
    parser.add_argument(
        "-r",
        "--ref",
        metavar="FILE",
        help="reference seqeunce file",
        required=True,
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
