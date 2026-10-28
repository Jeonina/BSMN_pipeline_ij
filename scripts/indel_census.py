#!/usr/bin/env python3
"""Indel census: how many somatic indel candidates do the 200 per-sample Mutect2 VCFs hold?

The BSMN A-F cascade keeps PASS SNVs only (accessibility_filter extracts SNVs), so
every Mutect2 indel is dropped today.  Before designing an indel filter (it needs
PI approval), this counts the PASS indels per sample through a few cheap,
SMaHT #12-like stages, using only what Mutect2 already annotates:

  pass        FILTER == PASS, indel
  accessible  inside the 1000G strict mask (same BED as filter A)
  rare        INFO/POPAF > 3, i.e. population AF < 0.001 (POPAF is -log10 AF)
  low_vaf     AF < 0.3 and alt reads >= 5 (germline het/hom excluded, enough support)

At each stage it also counts indels Mutect2 flags as short tandem repeats
(INFO/STR), the main source of false indels (polymerase slippage).  This is a
census for planning, NOT a filter: nothing feeds back into the pipeline.

Usage (VM100, pipeline directory):
    python scripts/indel_census.py --sample-info ~/sample_info.tsv --results results \\
        --mask resources/hg38/1KG.20160622.strict_mask.hg38_GRCh38.bed \\
        --bcftools "apptainer exec containers/bcftools_1.17.sif bcftools" \\
        --out-dir /mnt/data/bsmn/indel_census --jobs 8
Outputs: indel_census.tsv (per sample x stage counts) and {sample}.low_vaf.tsv
(the last-stage candidates, for length / repeat-context plots on the Mac).
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

STAGES = ("pass", "accessible", "rare", "low_vaf")
FMT = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/STR\t%INFO/RU\t%INFO/RPA\t%INFO/POPAF[\t%AF\t%AD]\n"
COLS = ("chrom", "pos", "ref", "alt", "str", "ru", "rpa", "popaf", "af", "ad")


def parse_line(line: str) -> dict | None:
    """One bcftools query line -> record; None for multi-allelic or malformed lines."""
    f = line.rstrip("\n").split("\t")
    if len(f) != len(COLS) or "," in f[3]:
        return None
    r = dict(zip(COLS, f))
    ad = r["ad"].split(",")
    try:
        r["alt_reads"] = int(ad[1]) if len(ad) > 1 else 0
        r["vaf"] = float(r["af"])
        r["popaf_v"] = float(r["popaf"]) if r["popaf"] not in (".", "") else 10.0   # absent = not in gnomAD
    except ValueError:
        return None
    r["is_str"] = r["str"] not in (".", "", "0")
    r["length"] = len(r["alt"]) - len(r["ref"])          # + insertion, - deletion
    return r


def stage_of(r: dict, accessible: bool) -> list[str]:
    """Stages a record passes, in cascade order."""
    out = ["pass"]
    if not accessible:
        return out
    out.append("accessible")
    if r["popaf_v"] <= 3:
        return out
    out.append("rare")
    if r["vaf"] < 0.3 and r["alt_reads"] >= 5:
        out.append("low_vaf")
    return out


def run_query(args, vcf: str, mask: str | None) -> list[str]:
    bt = shlex.split(args.bcftools)
    view = bt + ["view", "-f", "PASS", "-v", "indels"] + (["-T", mask] if mask else []) + [vcf]
    query = bt + ["query", "-f", FMT]
    q = subprocess.run(["bash", "-o", "pipefail", "-c", f"{shlex.join(view)} | {shlex.join(query)}"],
                       capture_output=True, text=True)
    if q.returncode != 0:
        raise RuntimeError(q.stderr.strip()[:300])
    return q.stdout.splitlines()


def census_sample(s: str, args) -> dict:
    vcf = os.path.join(args.results, "calling", s, f"{s}.filtered.vcf.gz")
    keys = lambda lines: {(r["chrom"], r["pos"], r["ref"], r["alt"]): r
                          for r in map(parse_line, lines) if r}
    all_recs = keys(run_query(args, vcf, None))
    acc_keys = set(keys(run_query(args, vcf, args.mask)))
    counts = {f"{st}{suf}": 0 for st in STAGES for suf in ("", "_str")}
    last = []
    for k, r in all_recs.items():
        for st in stage_of(r, k in acc_keys):
            counts[st] += 1
            counts[f"{st}_str"] += r["is_str"]
        if "low_vaf" in stage_of(r, k in acc_keys):
            last.append(r)
    with open(os.path.join(args.out_dir, f".{s}.low_vaf.tsv"), "w") as fo:
        fo.write("\t".join(("sample", "chrom", "pos", "ref", "alt", "length", "vaf", "alt_reads", "is_str", "ru", "rpa")) + "\n")
        for r in last:
            fo.write("\t".join(map(str, (s, r["chrom"], r["pos"], r["ref"], r["alt"], r["length"], r["vaf"],
                                         r["alt_reads"], int(r["is_str"]), r["ru"], r["rpa"]))) + "\n")
    os.replace(os.path.join(args.out_dir, f".{s}.low_vaf.tsv"), os.path.join(args.out_dir, f"{s}.low_vaf.tsv"))
    return {"sample": s, **counts}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-info", required=True, help="TSV: sample_id, br_region, sex")
    ap.add_argument("--results", default="results")
    ap.add_argument("--mask", required=True, help="1000G strict mask BED (filter A)")
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    with open(args.sample_info) as fh:
        samples = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]

    rows, failed = [], 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(census_sample, s, args): s for s in samples}
        for f in as_completed(futs):
            try:
                r = f.result(); rows.append(r)
                print(f"{r['sample']}: pass {r['pass']} -> accessible {r['accessible']} -> rare {r['rare']}"
                      f" -> low_vaf {r['low_vaf']} (STR {r['low_vaf_str']})", flush=True)
            except Exception as e:            # keep going; report at the end
                failed += 1
                print(f"ERROR {futs[f]}: {e}", file=sys.stderr, flush=True)
    out = os.path.join(args.out_dir, "indel_census.tsv")
    with open(out, "w", newline="") as fo:
        w = csv.DictWriter(fo, ["sample", *[f"{st}{suf}" for st in STAGES for suf in ("", "_str")]], delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["sample"]))
    print(f"wrote {out} ({len(rows)} samples, {failed} failed)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
