#!/usr/bin/env python3
"""Sample identity QC: do the two regions of each donor come from one person, and is the sex right?

Lab-meeting request (2026-10-21), and a precondition for donor-merged calling
(scripts/merge_donor_crams.py): merging a swapped sample would mix two people.

  pileup   per sample, count REF/ALT reads at common germline SNPs on chr22 and
           chrX (non-PAR) with bcftools mpileup (index jumps, cheap), and the
           mean depth of chrY (male-specific region, MQ>=30) vs a chr22 window.
           One TSV per sample under --out-dir; finished samples are skipped.
  analyze  pure Python (no numpy): pairwise Pearson correlation of chr22 alt
           allele fractions -> identity matrix; per-donor pair check; sex from
           chrX heterozygosity + chrY/chr22 depth vs the metadata sex.

Same person -> r ~ 1; unrelated people share only population frequencies, so r
is far lower.  Low-level contamination (a few %) barely moves r; it is NOT a
contamination test (see CalculateContamination / the cross-sample check).

Usage (VM100, pipeline directory):
    python scripts/sample_identity.py pileup \\
        --sample-info ~/sample_info.tsv --results results \\
        --sites-vcf resources/hg38/small_exac_common_3.hg38.vcf.gz \\
        --ref resources/hg38/Homo_sapiens_assembly38.fasta \\
        --bcftools "apptainer exec containers/bcftools_1.17.sif bcftools" \\
        --samtools "apptainer exec containers/samtools_1.17.sif samtools" \\
        --out-dir /mnt/data/bsmn/identity --jobs 16
    python scripts/sample_identity.py analyze \\
        --sample-info ~/sample_info.tsv --out-dir /mnt/data/bsmn/identity
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

# hg38 pseudo-autosomal regions are diploid in both sexes: exclude them.
X_PAR = ((10_001, 2_781_479), (155_701_383, 156_030_895))
Y_MSY = "chrY:2781480-26600000"      # euchromatic male-specific region
CHR22_WINDOW = "chr22:20000000-21000000"
MIN_DP = 20
MIN_SITES = 100                       # shared informative sites for a usable r


# ---------------------------------------------------------------- shared ----

def load_samples(path: str) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["donor"], _, r["region"] = r["br_region"].rpartition("_")
    return rows


def in_par(pos: int) -> bool:
    return any(a <= pos <= b for a, b in X_PAR)


# ---------------------------------------------------------------- pileup ----

def write_sites(vcf: str, out: str, af_min: float = 0.05) -> int:
    """Biallelic SNPs on chr22 / chrX non-PAR with af_min <= AF <= 1-af_min -> CHROM POS REF ALT."""
    n = 0
    with gzip.open(vcf, "rt") as fh, open(out, "w") as fo:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, info = line.split("\t", 8)[:8]
            if chrom not in ("chr22", "chrX") or len(ref) != 1 or len(alt) != 1:
                continue
            if chrom == "chrX" and in_par(int(pos)):
                continue
            af = next((float(f[3:]) for f in info.split(";") if f.startswith("AF=")), None)
            if af is None or not af_min <= af <= 1 - af_min:
                continue
            fo.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")
            n += 1
    return n


def parse_mpileup(text: str, sites: dict[tuple[str, int], tuple[str, str]]) -> list[tuple]:
    """bcftools query lines 'CHROM POS REF ALTS AD' -> (chrom, pos, ref_n, alt_n) for the site's ALT."""
    out = []
    for line in text.splitlines():
        chrom, pos, ref, alts, ad = line.split("\t")
        key = (chrom, int(pos))
        if key not in sites or sites[key][0] != ref:
            continue
        counts = [int(c) if c != "." else 0 for c in ad.split(",")]
        alleles = alts.split(",")
        want = sites[key][1]
        alt_n = counts[1 + alleles.index(want)] if want in alleles else 0
        out.append((chrom, key[1], counts[0], alt_n))
    return out


def pileup_sample(s: dict, args, sites: dict, regions_path: str) -> str:
    sid = s["sample_id"]
    cram = os.path.join(args.results, "mapping", sid, f"{sid}.cram")
    bt, st = shlex.split(args.bcftools), shlex.split(args.samtools)
    mpileup = bt + ["mpileup", "-f", args.ref, "-R", regions_path, "-a", "AD",
                    "-q", "20", "-Q", "20", "-d", "10000", "-I", "-Ou", cram]
    query = bt + ["query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n", "-"]
    q = subprocess.run(["bash", "-o", "pipefail", "-c", f"{shlex.join(mpileup)} | {shlex.join(query)}"],
                       capture_output=True, text=True)
    if q.returncode != 0:
        raise RuntimeError(f"{sid}: mpileup failed: {q.stderr[-300:]}")
    rows = parse_mpileup(q.stdout, sites)

    depth = {}
    for name, region in (("chrY", Y_MSY), ("chr22", CHR22_WINDOW)):
        res = subprocess.run(st + ["coverage", "--reference", args.ref, "--min-MQ", "30",
                                   "-r", region, "-H", cram], capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"{sid}: samtools coverage {region}: {res.stderr[-300:]}")
        depth[name] = float(res.stdout.split("\t")[6])

    final = os.path.join(args.out_dir, f"{sid}.tsv")
    with open(final + ".tmp", "w") as fo:
        fo.write(f"#depth\tchrY={depth['chrY']}\tchr22={depth['chr22']}\n")
        fo.write("chrom\tpos\tref_n\talt_n\n")
        fo.writelines(f"{c}\t{p}\t{r}\t{a}\n" for c, p, r, a in rows)
    os.replace(final + ".tmp", final)
    return f"{sid}\t{len(rows)} sites\tchrY/chr22={depth['chrY'] / max(depth['chr22'], 1e-9):.3f}"


def cmd_pileup(args) -> int:
    os.makedirs(args.out_dir, exist_ok=True)
    sites_path = os.path.join(args.out_dir, "sites.tsv")
    if not os.path.exists(sites_path):
        print(f"sites : {write_sites(args.sites_vcf, sites_path)} written to {sites_path}")
    sites = {}
    with open(sites_path) as fh:
        for line in fh:
            c, p, r, a = line.split()
            sites[(c, int(p))] = (r, a)
    # mpileup -R reads a 3rd column as END, so give it CHROM POS only.
    regions_path = os.path.join(args.out_dir, "sites.regions")
    with open(regions_path, "w") as fo:
        fo.writelines(f"{c}\t{p}\n" for c, p in sites)

    todo, missing = [], []
    for s in load_samples(args.sample_info):
        sid = s["sample_id"]
        if os.path.exists(os.path.join(args.out_dir, f"{sid}.tsv")):
            continue
        cram = os.path.join(args.results, "mapping", sid, f"{sid}.cram")
        (todo if os.path.exists(cram + ".crai") else missing).append(s)
    print(f"to run: {len(todo)}  missing CRAM: {' '.join(s['sample_id'] for s in missing) or '-'}")

    failed = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = [pool.submit(pileup_sample, s, args, sites, regions_path) for s in todo]
        for fut in as_completed(futs):
            try:
                print(fut.result(), flush=True)
            except Exception as exc:  # report every sample, fail at the end
                failed += 1
                print(f"FAILED\t{exc}", flush=True)
    return 1 if failed else 0


# --------------------------------------------------------------- analyze ----

def load_counts(path: str) -> tuple[dict, dict]:
    """-> ({(chrom, pos): alt_fraction} for dp >= MIN_DP, {"chrY": d, "chr22": d})."""
    af, depth = {}, {}
    with open(path) as fh:
        for kv in fh.readline().rstrip("\n").split("\t")[1:]:
            k, v = kv.split("=")
            depth[k] = float(v)
        next(fh)
        for line in fh:
            c, p, r, a = line.split("\t")
            dp = int(r) + int(a)
            if dp >= MIN_DP:
                af[(c, int(p))] = int(a) / dp
    return af, depth


def pearson(x: dict, y: dict) -> tuple[float, int]:
    keys = [k for k in x.keys() & y.keys() if k[0] == "chr22"]
    n = len(keys)
    if n < MIN_SITES:
        return float("nan"), n
    xs, ys = [x[k] for k in keys], [y[k] for k in keys]
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return (sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")), n


def x_het_rate(af: dict) -> float:
    """Fraction of informative chrX non-PAR sites that look heterozygous (males ~0)."""
    xs = [v for k, v in af.items() if k[0] == "chrX"]
    return sum(0.2 <= v <= 0.8 for v in xs) / len(xs) if xs else float("nan")


def infer_sex(het: float, y_ratio: float) -> str:
    if y_ratio >= 0.2 and het < 0.05:
        return "M"
    if y_ratio < 0.05 and het >= 0.1:
        return "F"
    return "?"


def cmd_analyze(args) -> int:
    samples = [s for s in load_samples(args.sample_info)
               if os.path.exists(os.path.join(args.out_dir, f"{s['sample_id']}.tsv"))]
    data = {s["sample_id"]: load_counts(os.path.join(args.out_dir, f"{s['sample_id']}.tsv")) for s in samples}
    ids = [s["sample_id"] for s in samples]
    by_id = {s["sample_id"]: s for s in samples}
    label = {i: by_id[i]["br_region"] for i in ids}

    r = {i: {} for i in ids}
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            v, _ = pearson(data[ids[a]][0], data[ids[b]][0])
            r[ids[a]][ids[b]] = r[ids[b]][ids[a]] = v

    with open(os.path.join(args.out_dir, "identity_matrix.tsv"), "w") as fo:
        fo.write("\t" + "\t".join(label[i] for i in ids) + "\n")
        for i in ids:
            fo.write(label[i] + "\t" + "\t".join("1" if i == j else f"{r[i][j]:.4f}" for j in ids) + "\n")

    # Per sample: expected partner (other region of the same donor) vs the best other match.
    problems = 0
    with open(os.path.join(args.out_dir, "identity_pairs.tsv"), "w") as fo:
        fo.write("sample\tbr_region\tpartner\tr_partner\tbest_match\tr_best\tnext_best_r\tn_sites\tstatus\n")
        for i in ids:
            s = by_id[i]
            partner = next((j for j in ids if j != i and by_id[j]["donor"] == s["donor"]), None)
            ranked = sorted(((v, j) for j, v in r[i].items() if not math.isnan(v)), reverse=True)
            best_r, best = ranked[0] if ranked else (float("nan"), "")
            others = [v for v, j in ranked if j != partner]
            next_r = others[0] if others else float("nan")
            rp = r[i].get(partner, float("nan")) if partner else float("nan")
            n = pearson(data[i][0], data[partner][0])[1] if partner else 0
            if partner is None:
                status = "PARTNER_NOT_RUN"
            elif math.isnan(rp):
                status = "TOO_FEW_SITES"
            elif best == partner and rp >= args.min_r and next_r < args.min_r:
                status = "OK"
            elif best != partner and best_r >= args.min_r:
                status = "SWAP?"            # matches another person's sample better
            elif rp < args.min_r:
                status = "PAIR_MISMATCH"    # partner is not the same person
            else:
                status = "DUPLICATE_OR_RELATIVE?"
            problems += status != "OK" and status != "PARTNER_NOT_RUN"
            fo.write(f"{i}\t{label[i]}\t{label.get(partner, '-')}\t{rp:.4f}\t{label.get(best, '-')}\t"
                     f"{best_r:.4f}\t{next_r:.4f}\t{n}\t{status}\n")

    sex_bad = 0
    with open(os.path.join(args.out_dir, "identity_sex.tsv"), "w") as fo:
        fo.write("sample\tbr_region\tsex_meta\tsex_inferred\tx_het_rate\tchrY_chr22_depth\tstatus\n")
        for i in ids:
            af, depth = data[i]
            het, yr = x_het_rate(af), depth["chrY"] / max(depth["chr22"], 1e-9)
            inferred = infer_sex(het, yr)
            meta = by_id[i].get("sex", "")
            status = "OK" if inferred == meta else ("UNCLEAR" if inferred == "?" else "MISMATCH")
            sex_bad += status != "OK"
            fo.write(f"{i}\t{label[i]}\t{meta}\t{inferred}\t{het:.4f}\t{yr:.4f}\t{status}\n")

    print(f"samples: {len(ids)}  pair problems: {problems}  sex not OK: {sex_bad}")
    print(f"tables : {args.out_dir}/identity_{{matrix,pairs,sex}}.tsv")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pileup")
    p.add_argument("--sample-info", required=True, help="TSV: sample_id, br_region, sex")
    p.add_argument("--results", default="results")
    p.add_argument("--sites-vcf", required=True, help="common germline SNPs with INFO/AF")
    p.add_argument("--ref", required=True)
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--jobs", type=int, default=16)
    a = sub.add_parser("analyze")
    a.add_argument("--sample-info", required=True)
    a.add_argument("--out-dir", required=True)
    a.add_argument("--min-r", type=float, default=0.8, help="r at or above = same person")
    args = ap.parse_args(argv)
    return cmd_pileup(args) if args.cmd == "pileup" else cmd_analyze(args)


if __name__ == "__main__":
    sys.exit(main())
