#!/usr/bin/env python3
"""Quantify the overlap between the 1KG accessibility mask and the GIAB HG002
mosaic confident region — with NO external dependencies (stdlib only).

Answers: are the two region sets independent, overlapping, or is one an
extension (superset) of the other?

  - 1KG strict mask (accessibility): population mappability mask, used as the
    production accessibility filter (Step B). BED, 0-based half-open.
  - GIAB HG002 mosaic confident BED: benchmark-scoring region for HG002 (a
    ~2.45 Gbp subset of the GIAB germline benchmark regions, autosomes only).

The two are built by different consortia from different evidence, so any
overlap is incidental, not by construction. This script measures it using a
per-chromosome interval merge + two-pointer sweep (fast, pure Python).

Usage:
    python scripts/compare_masks.py \
        [--kg  resources/hg38/1KG.20160622.strict_mask.hg38_GRCh38.bed] \
        [--giab resources/hg38/HG002_GRCh38_MosaicSNVv1.1.bed] \
        [--outdir reports/mask_overlap] \
        [--no-download]

If --giab is missing it is downloaded from the GIAB FTP (urllib, no curl).
Comparison is restricted to autosomes (chr1..chr22) since the GIAB mosaic
confident set is autosomes only; both inputs are coerced to chr-prefixed names.
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
import urllib.request

GIAB_URL = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "AshkenazimTrio/HG002_NA24385_son/mosaic_v1.10/GRCh38/SNV/"
    "HG002_GRCh38_MosaicSNVv1.1.bed"
)

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_bed(path: str) -> dict[str, list[tuple[int, int]]]:
    """Load a BED into chrom -> merged sorted [(start, end)], autosomes only.

    chrom names are coerced to chr-prefixed; intervals are 0-based half-open.
    """
    raw: dict[str, list[tuple[int, int]]] = {}
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            chrom = f[0] if f[0].startswith("chr") else f"chr{f[0]}"
            if chrom not in AUTOSOMES:
                continue
            raw.setdefault(chrom, []).append((int(f[1]), int(f[2])))
    return {c: _merge(ivs) for c, ivs in raw.items()}


def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Sort and merge overlapping/adjacent intervals."""
    if not intervals:
        return []
    intervals.sort()
    out = [intervals[0]]
    for s, e in intervals[1:]:
        ls, le = out[-1]
        if s <= le:  # overlap or touch
            out[-1] = (ls, max(le, e))
        else:
            out.append((s, e))
    return out


def total_bp(bed: dict[str, list[tuple[int, int]]]) -> int:
    return sum(e - s for ivs in bed.values() for s, e in ivs)


def intersect(
    a: dict[str, list[tuple[int, int]]], b: dict[str, list[tuple[int, int]]]
) -> dict[str, list[tuple[int, int]]]:
    """Per-chrom two-pointer intersection of two merged interval sets."""
    out: dict[str, list[tuple[int, int]]] = {}
    for chrom in a.keys() & b.keys():
        ai = bi = 0
        av, bv = a[chrom], b[chrom]
        res: list[tuple[int, int]] = []
        while ai < len(av) and bi < len(bv):
            lo = max(av[ai][0], bv[bi][0])
            hi = min(av[ai][1], bv[bi][1])
            if lo < hi:
                res.append((lo, hi))
            # advance the interval that ends first
            if av[ai][1] < bv[bi][1]:
                ai += 1
            else:
                bi += 1
        if res:
            out[chrom] = res
    return out


def maybe_download(path: str, allow: bool) -> None:
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return
    if not allow:
        sys.exit(f"ERROR: {path} not found and --no-download set")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    print(f"[info] downloading GIAB confident BED -> {path}", file=sys.stderr)
    urllib.request.urlretrieve(GIAB_URL, path)  # noqa: S310 (trusted GIAB FTP)


def write_bed(bed: dict[str, list[tuple[int, int]]], path: str) -> None:
    with open(path, "w") as fh:
        for chrom in sorted(bed, key=lambda c: int(c[3:])):
            for s, e in bed[chrom]:
                fh.write(f"{chrom}\t{s}\t{e}\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kg", default="resources/hg38/1KG.20160622.strict_mask.hg38_GRCh38.bed")
    ap.add_argument("--giab", default="resources/hg38/HG002_GRCh38_MosaicSNVv1.1.bed")
    ap.add_argument("--outdir", default="reports/mask_overlap")
    ap.add_argument("--no-download", action="store_true")
    args = ap.parse_args()

    maybe_download(args.giab, allow=not args.no_download)
    if not (os.path.exists(args.kg) and os.path.getsize(args.kg) > 0):
        sys.exit(f"ERROR: 1KG mask BED not found at {args.kg}")

    os.makedirs(args.outdir, exist_ok=True)
    print("[info] loading + merging (autosomes only)...", file=sys.stderr)
    kg = read_bed(args.kg)
    giab = read_bed(args.giab)

    kg_bp = total_bp(kg)
    giab_bp = total_bp(giab)
    inter = intersect(kg, giab)
    inter_bp = total_bp(inter)
    union_bp = kg_bp + giab_bp - inter_bp

    write_bed(inter, os.path.join(args.outdir, "intersection.bed"))

    lines = [
        "==============================================================",
        " 1KG accessibility  vs  GIAB HG002 mosaic confident  (autosomes)",
        "==============================================================",
        f"1KG mask BED        : {args.kg}",
        f"GIAB confident BED  : {args.giab}",
        "--------------------------------------------------------------",
        f"1KG accessible (A)        : {kg_bp:15,d} bp",
        f"GIAB confident (B)        : {giab_bp:15,d} bp",
        f"Intersection (A∩B)        : {inter_bp:15,d} bp",
        f"1KG-only   (A∖B)          : {kg_bp - inter_bp:15,d} bp",
        f"GIAB-only  (B∖A)          : {giab_bp - inter_bp:15,d} bp",
        f"Union      (A∪B)          : {union_bp:15,d} bp",
        "--------------------------------------------------------------",
        f"Fraction of 1KG  in GIAB  : {100 * inter_bp / kg_bp:6.2f} %   (A∩B / A)" if kg_bp else "",
        f"Fraction of GIAB in 1KG   : {100 * inter_bp / giab_bp:6.2f} %   (A∩B / B)"
        if giab_bp
        else "",
        f"Jaccard (A∩B / A∪B)       : {inter_bp / union_bp:6.4f}" if union_bp else "",
        "==============================================================",
        "Interpretation:",
        "  * both fractions high (>90%) AND neither =100%",
        "      => heavily OVERLAPPING but INDEPENDENT (neither is a subset)",
        "  * one fraction ~100% => that set is a SUBSET (other = extension)",
        "  * both fractions low => largely DISJOINT",
        "==============================================================",
    ]
    text = "\n".join(line for line in lines if line != "")
    print(text)
    summary = os.path.join(args.outdir, "mask_overlap_summary.txt")
    with open(summary, "w") as fh:
        fh.write(text + "\n")
    print(f"\n[done] summary -> {summary}", file=sys.stderr)
    print(f"[done] intersection -> {args.outdir}/intersection.bed", file=sys.stderr)


if __name__ == "__main__":
    main()
