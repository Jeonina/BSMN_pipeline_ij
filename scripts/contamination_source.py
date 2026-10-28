#!/usr/bin/env python3
"""Cross-sample contamination check for every sample: which donor's germline do its calls match?

The GATK CalculateContamination estimate says how much foreign DNA a sample
carries, not whose, and the report's 2% cut-off was ours, not GATK's.  This
extends the cross-sample check (until now run only on samples above 2%) to all
samples, so the cut-off can be read off the data: at what estimate do foreign
germline variants start to appear among the final calls?

  query    per sample, look up the union of all final-call sites in its Mutect2
           VCF (results/calling/{s}/{s}.filtered.vcf.gz) and keep the alleles it
           carries as germline (AF >= --germline-af).  One TSV per sample under
           --out-dir; finished samples are skipped.  Index jumps only (~6k sites).
  analyze  pure Python: for each sample, count its final calls that are germline
           in another donor (the donor's own other region is excluded: shared
           somatic calls are real).  Reports the top and second donor.

Usage (VM100, pipeline directory):
    python scripts/contamination_source.py query \\
        --sample-info ~/sample_info.tsv --results results \\
        --bcftools "apptainer exec containers/bcftools_1.17.sif bcftools" \\
        --out-dir /mnt/data/bsmn/contam_source --jobs 8
    python scripts/contamination_source.py analyze \\
        --sample-info ~/sample_info.tsv --results results --out-dir /mnt/data/bsmn/contam_source
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

Site = tuple[str, int, str, str]      # chrom, pos, ref, alt


def load_samples(path: str) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["donor"] = r["br_region"].rpartition("_")[0]
    return rows


def load_final(path: str) -> set[Site]:
    """final.txt: CHROM POS REF ALT per line."""
    out = set()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                out.add((f[0], int(f[1]), f[2], f[3]))
    return out


def parse_query(text: str, germline_af: float) -> set[Site]:
    """bcftools query '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AF]' -> alleles with AF >= germline_af."""
    out = set()
    for line in text.splitlines():
        chrom, pos, ref, alts, afs = line.split("\t")[:5]
        for alt, af in zip(alts.split(","), afs.split(",")):
            if af not in (".", "") and float(af) >= germline_af:
                out.add((chrom, int(pos), ref, alt))
    return out


def source_table(finals: dict[str, set[Site]], germ: dict[str, set[Site]],
                 donor_of: dict[str, str]) -> list[dict]:
    """Per sample: final calls found as germline in each other donor -> top and second donor."""
    germ_by_donor: dict[str, set[Site]] = {}
    for s, sites in germ.items():
        germ_by_donor.setdefault(donor_of[s], set()).update(sites)
    rows = []
    for s, calls in finals.items():
        hits = Counter({d: len(calls & g) for d, g in germ_by_donor.items() if d != donor_of[s]})
        ranked = hits.most_common(2) + [("", 0)] * 2
        (top, n1), (second, n2) = ranked[0], ranked[1]
        rows.append({
            "sample": s, "donor": donor_of[s], "n_final": len(calls),
            "top_source": top if n1 else "", "top_hits": n1,
            "top_frac": round(n1 / len(calls), 4) if calls else 0.0,
            "second_source": second if n2 else "", "second_hits": n2,
            "any_foreign_germline": len(calls & set().union(
                *(g for d, g in germ_by_donor.items() if d != donor_of[s]))),
        })
    return rows


# ----------------------------------------------------------------- query ----

def query_sample(s: str, args, regions: str) -> str:
    vcf = os.path.join(args.results, "calling", s, f"{s}.filtered.vcf.gz")
    cmd = shlex.split(args.bcftools) + ["query", "-R", regions,
                                        "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n", vcf]
    q = subprocess.run(cmd, capture_output=True, text=True)
    if q.returncode != 0:
        raise RuntimeError(f"{s}: bcftools query failed: {q.stderr.strip()[:300]}")
    sites = parse_query(q.stdout, args.germline_af)
    tmp = os.path.join(args.out_dir, f".{s}.germline.tsv")
    with open(tmp, "w") as fo:
        for c, p, r, a in sorted(sites):
            fo.write(f"{c}\t{p}\t{r}\t{a}\n")
    os.replace(tmp, os.path.join(args.out_dir, f"{s}.germline.tsv"))
    return f"{s}: {len(sites)} germline alleles at final-call sites"


def cmd_query(args) -> int:
    os.makedirs(args.out_dir, exist_ok=True)
    samples = [r["sample_id"] for r in load_samples(args.sample_info)]
    union: set[tuple[str, int]] = set()
    for s in samples:
        union |= {(c, p) for c, p, _, _ in load_final(os.path.join(args.results, "filtering", s, f"{s}.final.txt"))}
    regions = os.path.join(args.out_dir, "final_sites.regions.tsv")
    with open(regions, "w") as fo:
        for c, p in sorted(union):
            fo.write(f"{c}\t{p}\n")
    print(f"{len(union)} final-call sites from {len(samples)} samples", flush=True)

    todo = [s for s in samples if not os.path.exists(os.path.join(args.out_dir, f"{s}.germline.tsv"))]
    failed = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(query_sample, s, args, regions): s for s in todo}
        for f in as_completed(futs):
            try:
                print(f.result(), flush=True)
            except Exception as e:            # keep going; report at the end
                failed += 1
                print(f"ERROR {e}", file=sys.stderr, flush=True)
    print(f"done {len(todo) - failed}/{len(todo)} (skipped {len(samples) - len(todo)} finished)")
    return 1 if failed else 0


# --------------------------------------------------------------- analyze ----

def cmd_analyze(args) -> int:
    info = load_samples(args.sample_info)
    donor_of = {r["sample_id"]: r["donor"] for r in info}
    finals, germ = {}, {}
    for s in donor_of:
        finals[s] = load_final(os.path.join(args.results, "filtering", s, f"{s}.final.txt"))
        g = os.path.join(args.out_dir, f"{s}.germline.tsv")
        if not os.path.exists(g):
            print(f"ERROR: {g} missing; run query first", file=sys.stderr)
            return 1
        germ[s] = load_final(g)
    rows = source_table(finals, germ, donor_of)
    region = {r["sample_id"]: r["br_region"] for r in info}
    out = os.path.join(args.out_dir, "contamination_source.tsv")
    with open(out, "w", newline="") as fo:
        w = csv.DictWriter(fo, ["sample", "br_region", "donor", "n_final", "top_source", "top_hits", "top_frac",
                                "second_source", "second_hits", "any_foreign_germline"], delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda r: -r["top_frac"]):
            w.writerow({**r, "br_region": region[r["sample"]]})
    print(f"wrote {out}")
    for r in sorted(rows, key=lambda r: -r["top_frac"])[:12]:
        print(f"  {r['sample']} {region[r['sample']]:16} n={r['n_final']:4} top={r['top_source'] or '-':7} "
              f"{r['top_hits']:4} ({r['top_frac']:.0%})  2nd={r['second_hits']}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("query", "analyze"):
        p = sub.add_parser(name)
        p.add_argument("--sample-info", required=True, help="TSV: sample_id, br_region, sex")
        p.add_argument("--results", default="results")
        p.add_argument("--out-dir", required=True)
    q = sub.choices["query"]
    q.add_argument("--bcftools", default="bcftools")
    q.add_argument("--jobs", type=int, default=8)
    q.add_argument("--germline-af", type=float, default=0.3, help="AF at or above = germline in that sample")
    args = ap.parse_args(argv)
    return cmd_query(args) if args.cmd == "query" else cmd_analyze(args)


if __name__ == "__main__":
    sys.exit(main())
