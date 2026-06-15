#!/usr/bin/env python3
"""Post-hoc panel-of-normals masking test for mosaic candidates.

Intersects a final mosaic candidate list (``chr pos ref alt``) against one or
more panel-of-normals VCFs (e.g. GATK ``1000g_pon.hg38.vcf.gz``) WITHOUT
re-running Mutect2. A candidate site present in a PON is a recurrent
artifact / common-variant site and would be masked if calling-time PON were
enabled (``config.yaml`` ``calling.pon``).

Reports, per PON, how many false-positive candidates it masks and whether the
known true positive survives. Use this to decide whether to turn on
calling-time PON before paying the 6 h Mutect2 re-run.

PON masking is position-level by default (Mutect2 masks a site regardless of
the specific alt); pass ``--match-alt`` to require the alt allele to match a
PON ALT as well.

Usage:
    python scripts/test_pon_mask.py \
        --candidates reports/benchmark_2026-06-02/HG002.final.txt \
        --pon resources/1000g_pon.hg38.vcf.gz \
        --truth chr20:47053475 \
        --bcftools-sif containers/bcftools_1.17.sif
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys


def run_bcftools(args_list: list[str], sif: str | None) -> str:
    if sif:
        cmd = ["apptainer", "exec", sif, "bcftools"] + args_list
    else:
        cmd = ["bcftools"] + args_list
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return str(result.stdout)


def read_candidates(path: str) -> list[tuple[str, int, str, str]]:
    out: list[tuple[str, int, str, str]] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            out.append((p[0], int(p[1]), p[2].upper(), p[3].upper()))
    return out


def pon_hits(
    pon: str, sif: str | None, cands: list[tuple[str, int, str, str]], match_alt: bool
) -> set[tuple[str, int]]:
    """Return the (chrom, pos) candidate sites present in the PON."""
    regions = ",".join(f"{c}:{pos}-{pos}" for c, pos, _, _ in cands)
    output = run_bcftools(["view", "-H", "-r", regions, pon], sif)
    cand_alt = {(c, pos): alt for c, pos, _, alt in cands}
    hits: set[tuple[str, int]] = set()
    for line in output.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        chrom, pos = f[0], int(f[1])
        key = (chrom, pos)
        if key not in cand_alt:
            continue
        if match_alt:
            pon_alts = {a.upper() for a in f[4].split(",")}
            if cand_alt[key] not in pon_alts:
                continue
        hits.add(key)
    return hits


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--candidates", required=True, help="final mosaic list (chr pos ref alt)")
    p.add_argument("--pon", required=True, nargs="+", help="one or more PON VCF(.gz) paths")
    p.add_argument("--truth", default="", help="known true-positive coord, e.g. chr20:47053475")
    p.add_argument("--match-alt", action="store_true", help="require alt to match a PON ALT")
    p.add_argument("--bcftools-sif", default=None, help="Apptainer .sif for bcftools (optional)")
    args = p.parse_args()

    cands = read_candidates(args.candidates)
    truth_key: tuple[str, int] | None = None
    if args.truth:
        tc, tp = args.truth.split(":")
        truth_key = (tc, int(tp))
    n_fp = sum(1 for c, pos, _, _ in cands if (c, pos) != truth_key)

    # Per-PON masking, plus the union across all PONs.
    per_pon: dict[str, set[tuple[str, int]]] = {}
    union: set[tuple[str, int]] = set()
    for pon in args.pon:
        if not os.path.exists(pon):
            print(f"[error] PON not found: {pon}", file=sys.stderr)
            sys.exit(1)
        hits = pon_hits(pon, args.bcftools_sif, cands, args.match_alt)
        per_pon[pon] = hits
        union |= hits

    match_mode = "position+alt" if args.match_alt else "position"
    print(f"PON masking test  (FP pool = {n_fp}, match = {match_mode})")
    print(f"{'PON':<45}{'FP masked':>12}{'TP masked?':>13}")
    print("-" * 70)
    for pon, hits in per_pon.items():
        fp_masked = sum(1 for k in hits if k != truth_key)
        tp_masked = "YES (bad)" if truth_key and truth_key in hits else (
            "no" if truth_key else "n/a"
        )
        print(f"{os.path.basename(pon):<45}{f'{fp_masked}/{n_fp}':>12}{tp_masked:>13}")
    if len(args.pon) > 1:
        fp_u = sum(1 for k in union if k != truth_key)
        tp_u = "YES (bad)" if truth_key and truth_key in union else "no"
        print(f"{'UNION (all PONs)':<45}{f'{fp_u}/{n_fp}':>12}{tp_u:>13}")

    # List the masked FP sites so they can be eyeballed.
    masked_fp = sorted(k for k in union if k != truth_key)
    if masked_fp:
        print(f"\nMasked FP sites ({len(masked_fp)}):")
        for chrom, pos in masked_fp:
            print(f"  {chrom}:{pos}")
    else:
        print("\nNo FP candidates masked by the supplied PON(s).")


if __name__ == "__main__":
    main()
