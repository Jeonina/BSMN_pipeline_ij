#!/usr/bin/env python3
"""Cohort germline genotypes from the per-sample Mutect2 VCFs, for PLINK / KING.

There is no joint germline call set, but every sample's Mutect2 VCF
(results/calling/{s}/{s}.filtered.vcf.gz) keeps its germline SNPs with an AF.
At ~96x the AF separates het from hom-alt cleanly, and a common SNP that a
sample does not report is hom-ref.  Reading these VCFs avoids re-reading the
CRAMs (storage-bound: ~124 GB per sample over NFS).

  extract  per sample: biallelic SNPs on chr1-22/X inside --mask, written as
           CHROM POS REF ALT DOSAGE (1 = het, 2 = hom-alt, . = AF in a grey
           zone).  AF < --het-min is not kept (sequencing noise or low-level
           contamination).  One TSV per sample; finished samples are skipped.
  merge    sites non-ref in at least --min-carriers samples -> one VCF with GT
           for every sample (absent = 0/0), plus a PLINK --update-sex file.
           Streams a k-way merge of the sorted per-sample TSVs (a few GB of RAM).
           Sample names are br_region (Br5094_Caudate); use plink --double-id.

Usage (VM100, pipeline directory):
    python scripts/germline_genotypes.py extract --sample-info ~/sample_info.tsv --results results \\
        --mask resources/hg38/1KG.20160622.strict_mask.hg38_GRCh38.bed \\
        --bcftools "apptainer exec containers/bcftools_1.17.sif bcftools" \\
        --out-dir /mnt/data/bsmn/king --jobs 8
    python scripts/germline_genotypes.py merge --sample-info ~/sample_info.tsv --out-dir /mnt/data/bsmn/king \\
        --exclude Br5221_Caudate,Br5867_DLPFC,Br5296_Caudate,Br5586_DLPFC,Br6111_Caudate
"""

from __future__ import annotations

import argparse
import csv
import heapq
import os
import shlex
import subprocess
import sys
from collections.abc import Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import groupby

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CHROM_ORDER = {c: i for i, c in enumerate(CHROMS)}
Site = tuple[str, int, str, str]


def load_samples(path: str) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def af_to_dosage(af: float, het_min: float, het_max: float, hom_min: float) -> str | None:
    """1 = het, 2 = hom-alt, '.' = between bands, None = below het_min (treated as hom-ref)."""
    if af < het_min:
        return None
    if af <= het_max:
        return "1"
    if af >= hom_min:
        return "2"
    return "."


def parse_query(text: str, het_min: float, het_max: float, hom_min: float) -> list[tuple[Site, str]]:
    """bcftools query '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AF]' -> (site, dosage) for biallelic SNPs."""
    out = []
    for line in text.splitlines():
        c, p, r, a, af = line.split("\t")[:5]
        if c not in CHROM_ORDER or len(r) != 1 or len(a) != 1 or af in (".", "") or "," in af:
            continue
        d = af_to_dosage(float(af), het_min, het_max, hom_min)
        if d is not None:
            out.append(((c, int(p), r, a), d))
    return out


# --------------------------------------------------------------- extract ----

def extract_sample(s: str, args) -> str:
    vcf = os.path.join(args.results, "calling", s, f"{s}.filtered.vcf.gz")
    cmd = shlex.split(args.bcftools) + ["query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n"]
    if args.mask:
        cmd += ["-T", args.mask]                     # streams the VCF; no CRAM access
    q = subprocess.run(cmd + [vcf], capture_output=True, text=True)
    if q.returncode != 0:
        raise RuntimeError(f"{s}: bcftools query failed: {q.stderr.strip()[:300]}")
    rows = parse_query(q.stdout, args.het_min, args.het_max, args.hom_min)
    tmp = os.path.join(args.out_dir, f".{s}.gt.tsv")
    with open(tmp, "w") as fo:
        for (c, p, r, a), d in rows:
            fo.write(f"{c}\t{p}\t{r}\t{a}\t{d}\n")
    os.replace(tmp, os.path.join(args.out_dir, f"{s}.gt.tsv"))
    return f"{s}: {len(rows)} non-ref SNPs"


def cmd_extract(args) -> int:
    os.makedirs(args.out_dir, exist_ok=True)
    samples = [r["sample_id"] for r in load_samples(args.sample_info)]
    todo = [s for s in samples if not os.path.exists(os.path.join(args.out_dir, f"{s}.gt.tsv"))]
    failed = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(extract_sample, s, args): s for s in todo}
        for f in as_completed(futs):
            try:
                print(f.result(), flush=True)
            except Exception as e:            # keep going; report at the end
                failed += 1
                print(f"ERROR {e}", file=sys.stderr, flush=True)
    print(f"done {len(todo) - failed}/{len(todo)} (skipped {len(samples) - len(todo)} finished)")
    return 1 if failed else 0


# ----------------------------------------------------------------- merge ----

def merge_genotypes(per_sample: list[Iterable[tuple[Site, str]]], min_carriers: int) -> Iterator[tuple[Site, list[str]]]:
    """Sites non-ref (het/hom/grey) in >= min_carriers samples; GT per sample, absent = 0/0.

    k-way merge of per-sample rows already in genome order (extract writes them in
    VCF order), so memory stays ~ one row per sample instead of all 6e8 rows.
    """
    n = len(per_sample)
    streams = [_ordered(i, rows) for i, rows in enumerate(per_sample)]
    cur, gts = None, None
    for key, i, d in heapq.merge(*streams):
        if key != cur:
            if gts is not None and n - gts.count("0/0") >= min_carriers:
                yield _site(cur), gts
            cur, gts = key, ["0/0"] * n
        gts[i] = GT[d]
    if gts is not None and n - gts.count("0/0") >= min_carriers:
        yield _site(cur), gts


GT = {"1": "0/1", "2": "1/1", ".": "./."}


def _site(key: tuple) -> Site:
    return (CHROMS[key[0]], key[1], key[2], key[3])


def _ordered(i: int, rows: Iterable[tuple[Site, str]]) -> Iterator[tuple[tuple, int, str]]:
    """(sort key, sample index, dosage); sorts alleles within a position, refuses out-of-order input."""
    prev = (-1, -1)
    for pos, grp in groupby(rows, key=lambda r: (CHROM_ORDER[r[0][0]], r[0][1])):
        if pos <= prev:
            raise ValueError(f"sample #{i}: rows not in genome order at {CHROMS[pos[0]]}:{pos[1]}")
        prev = pos
        for (c, p, r, a), d in sorted(grp, key=lambda x: x[0][2:]):
            yield (pos[0], p, r, a), i, d


def read_gt(path: str) -> Iterator[tuple[Site, str]]:
    with open(path) as fh:
        for line in fh:
            c, p, r, a, d = line.rstrip("\n").split("\t")
            yield (c, int(p), r, a), d


def cmd_merge(args) -> int:
    excl = {x for x in args.exclude.split(",") if x}
    info = [r for r in load_samples(args.sample_info) if r["br_region"] not in excl]
    paths = [os.path.join(args.out_dir, f"{r['sample_id']}.gt.tsv") for r in info]
    for g in paths:
        if not os.path.exists(g):
            print(f"ERROR: {g} missing; run extract first", file=sys.stderr)
            return 1
    names = [r["br_region"] for r in info]
    out = os.path.join(args.out_dir, "cohort_germline.vcf")
    tmp = os.path.join(args.out_dir, ".cohort_germline.vcf")
    n_sites = n_x = 0
    with open(tmp, "w") as fo:
        fo.write("##fileformat=VCFv4.2\n")
        fo.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype from Mutect2 AF">\n')
        for c in CHROMS:
            fo.write(f"##contig=<ID={c}>\n")
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
        for (c, p, r, a), gts in merge_genotypes([read_gt(g) for g in paths], args.min_carriers):
            fo.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
            n_sites += 1
            n_x += c == "chrX"
            if n_sites % 1_000_000 == 0:
                print(f"  {n_sites:,} SNPs written (at {c}:{p})", flush=True)
    os.replace(tmp, out)
    sex_code = {"M": "1", "F": "2"}
    with open(os.path.join(args.out_dir, "sex_update.txt"), "w") as fo:
        for r in info:
            fo.write(f"{r['br_region']}\t{r['br_region']}\t{sex_code.get(r.get('sex', ''), '0')}\n")
    print(f"wrote {out}: {len(names)} samples, {n_sites} SNPs ({n_sites - n_x} autosomal, {n_x} chrX)")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("extract", "merge"):
        p = sub.add_parser(name)
        p.add_argument("--sample-info", required=True, help="TSV: sample_id, br_region, sex")
        p.add_argument("--out-dir", required=True)
    e = sub.choices["extract"]
    e.add_argument("--results", default="results")
    e.add_argument("--mask", default="", help="BED of accessible sequence (1KG strict mask)")
    e.add_argument("--bcftools", default="bcftools")
    e.add_argument("--jobs", type=int, default=8)
    e.add_argument("--het-min", type=float, default=0.25, help="AF below = not carried (hom-ref)")
    e.add_argument("--het-max", type=float, default=0.75)
    e.add_argument("--hom-min", type=float, default=0.85, help="AF at or above = hom-alt")
    m = sub.choices["merge"]
    m.add_argument("--exclude", default="", help="comma-separated br_region to leave out")
    m.add_argument("--min-carriers", type=int, default=10, help="keep SNPs non-ref in at least this many samples")
    args = ap.parse_args(argv)
    return cmd_extract(args) if args.cmd == "extract" else cmd_merge(args)


if __name__ == "__main__":
    sys.exit(main())
