#!/usr/bin/env python3
"""Write somalier .somalier files from bcftools mpileup counts (fast replacement for `somalier extract`).

`somalier extract` on ~100x CRAMs over NFS ran ~1 s per site (~5 h per sample at
17k sites; measured 2026-10-22). The cost is I/O, not the tool: each spread-out
site needs ~3 MB of CRAM and the storage tops out near 70 MB/s, so bcftools is no
faster per site. The fix is fewer sites (--auto-stride). The .somalier format is
just per-site read counts, so this script counts with bcftools and writes the same
binary; `somalier relate` and `somalier ancestry` then run unchanged against the
1000 Genomes .somalier files.

.somalier layout (somalier src/somalierpkg/common.nim, formatVersion 2), little-endian:
  uint8 version=2 | uint8 len(name) | name | uint16 n_auto | uint16 n_x | uint16 n_y
  then n_auto, n_x, n_y records of (uint32 nref, uint32 nalt, uint32 nother),
  each section in the order the sites appear in the sites VCF.

Usage (VM100, pipeline directory):
    APPTAINER_BIND=/mnt/data,/tmp python ~/somalier_counts.py \\
        --sample-info ~/sample_info.tsv --results results \\
        --sites /mnt/data/bsmn/somalier/sites.hg38.vcf.gz \\
        --ref resources/hg38/Homo_sapiens_assembly38.fasta \\
        --bcftools "apptainer exec containers/bcftools_1.17.sif bcftools" \\
        --out-dir /mnt/data/bsmn/somalier/extracted --jobs 24
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shlex
import struct
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

FORMAT_VERSION = 2
MIN_COVERED = 0.5   # fewer covered sites than this = chromosome naming / setup error


def select_sites(sites: list, auto_stride: int) -> list:
    """Every auto_stride-th autosomal site plus every X/Y site (the ones that get counted).

    Reading cost is per site (~3 MB of CRAM per site on these ~100x files), so a
    stride cuts I/O proportionally. Uncounted sites are written as 0/0/0, which
    somalier treats as no data and skips; ancestry and relatedness need far fewer
    than 17k sites.
    """
    out, k = [], 0
    for site in sites:
        if section(site[0]) == 0:
            k += 1
            if (k - 1) % auto_stride:
                continue
        out.append(site)
    return out


def load_sites(vcf: str) -> list[tuple[str, int, str, str]]:
    """Sites in VCF order: (chrom, pos, ref, alt)."""
    sites = []
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt = line.split("\t", 5)[:5]
            sites.append((chrom, int(pos), ref, alt.split(",")[0]))
    return sites


def section(chrom: str) -> int:
    c = chrom.removeprefix("chr")
    return 1 if c == "X" else 2 if c == "Y" else 0


def parse_counts(text: str, sites: list) -> dict[tuple[str, int], tuple[int, int, int]]:
    """bcftools query 'CHROM POS REF ALTS AD' -> {(chrom, pos): (nref, nalt, nother)} for the site's ALT."""
    want = {(c, p): (r, a) for c, p, r, a in sites}
    out = {}
    for line in text.splitlines():
        chrom, pos, ref, alts, ad = line.split("\t")
        key = (chrom, int(pos))
        if key not in want or want[key][0] != ref:
            continue
        counts = [int(x) if x != "." else 0 for x in ad.split(",")]
        alleles = alts.split(",")
        alt = want[key][1]
        i = 1 + alleles.index(alt) if alt in alleles else None
        nalt = counts[i] if i is not None else 0
        out[key] = (counts[0], nalt, sum(counts[1:]) - nalt)
    return out


def write_somalier(path: str, name: str, sites: list, counts: dict) -> None:
    secs: list[list[tuple[int, int, int]]] = [[], [], []]
    for c, p, _, _ in sites:
        secs[section(c)].append(counts.get((c, p), (0, 0, 0)))
    raw = name.encode()
    with open(path + ".tmp", "wb") as fo:
        fo.write(struct.pack("<BB", FORMAT_VERSION, len(raw)) + raw)
        fo.write(struct.pack("<HHH", *(len(s) for s in secs)))
        for s in secs:
            for rec in s:
                fo.write(struct.pack("<III", *rec))
    os.replace(path + ".tmp", path)


def count_sample(sid: str, args, sites: list, counted: list, regions: str) -> str:
    cram = os.path.join(args.results, "mapping", sid, f"{sid}.cram")
    bt = shlex.split(args.bcftools)
    mpileup = bt + ["mpileup", "-f", args.ref, "-R", regions, "-a", "AD",
                    "-q", "20", "-Q", "20", "-d", "10000", "-I", "-Ou", cram]
    query = bt + ["query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n", "-"]
    res = subprocess.run(["bash", "-o", "pipefail", "-c", f"{shlex.join(mpileup)} | {shlex.join(query)}"],
                         capture_output=True, text=True)
    if res.returncode != 0:
        err = res.stderr.strip()
        raise RuntimeError(f"{sid}: mpileup failed: {err[:400]}{' ... ' + err[-200:] if len(err) > 600 else ''}")
    counts = parse_counts(res.stdout, counted)
    covered = sum(1 for v in counts.values() if v[0] + v[1] > 0) / len(counted)
    if covered < MIN_COVERED:
        raise RuntimeError(f"{sid}: only {covered:.0%} of sites covered (chromosome naming?)")
    write_somalier(os.path.join(args.out_dir, f"{sid}.somalier"), sid, sites, counts)
    return f"{sid}\tcovered {covered:.1%}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-info", required=True, help="TSV with a sample_id column")
    ap.add_argument("--results", default="results")
    ap.add_argument("--sites", required=True, help="somalier sites VCF (same file the 1kg .somalier used)")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--jobs", type=int, default=24)
    ap.add_argument("--auto-stride", type=int, default=1,
                    help="count every Nth autosomal site only (X/Y always); 6 = ~2.9k sites")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    sites = load_sites(args.sites)
    counted = select_sites(sites, args.auto_stride)
    regions = os.path.join(args.out_dir, ".sites.regions")
    with open(regions, "w") as fo:
        fo.writelines(f"{c}\t{p}\n" for c, p, _, _ in counted)
    n = [0, 0, 0]
    for c, *_ in sites:
        n[section(c)] += 1
    print(f"sites : {len(sites)} (auto {n[0]}, X {n[1]}, Y {n[2]}); counted {len(counted)} "
          f"(autosomal stride {args.auto_stride})", flush=True)

    with open(args.sample_info) as fh:
        ids = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]
    todo, missing = [], []
    for sid in ids:
        if os.path.exists(os.path.join(args.out_dir, f"{sid}.somalier")):
            continue
        cram = os.path.join(args.results, "mapping", sid, f"{sid}.cram")
        (todo if os.path.exists(cram + ".crai") else missing).append(sid)
    print(f"to run: {len(todo)}  missing CRAM: {' '.join(missing) or '-'}", flush=True)

    failed = ok = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = [pool.submit(count_sample, s, args, sites, counted, regions) for s in todo]
        for fut in as_completed(futs):
            try:
                print(fut.result(), flush=True)
                ok += 1
            except Exception as exc:  # report every sample, fail at the end
                failed += 1
                print(f"FAILED\t{exc}", flush=True)
                if failed >= 3 and not ok:   # a setup error, not a bad sample: stop
                    for f in futs:
                        f.cancel()
                    print("ABORT: first 3 samples failed, nothing succeeded", flush=True)
                    return 2
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
