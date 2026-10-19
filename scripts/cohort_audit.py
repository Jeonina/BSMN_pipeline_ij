#!/usr/bin/env python3
"""Cohort audit: every expected sample mapped, called and filtered exactly once.

The cohort was processed on two hosts (VM100, VM102) that write into ONE shared
NFS ``results/`` tree but keep separate ``logs/`` directories.  This script
checks, per expected sample (from ``barcode_map.tsv``):

  mapping    CRAM + CRAI + flagstat present; CRAM @RG set == expected readgroups
             (PU = <fc>_<udb>, SM = sample) — catches merged/missing/foreign rgs
  calling    filtered.vcf.gz + .tbi present; VCF sample column == sample
  filtering  final.txt present and pon_mask_filter log reached END
  once-only  each stage's log exists on exactly ONE host (two = ran twice)
  staleness  outputs are newer than their inputs (cram < vcf < final)

and collects the per-step filtering cascade from the logs (the intermediate
``.txt`` files are temp() and already deleted).

Usage (on VM100, after copying VM102's logs onto the shared NFS):
    python scripts/cohort_audit.py \\
        --barcode-map barcode_map.tsv --results results \\
        --log-dir vm100=logs --log-dir vm102=/mnt/data/bsmn/logs_vm102 \\
        --samtools "apptainer exec containers/samtools_1.17.sif samtools" \\
        --out reports/cohort_audit.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shlex
import statistics
import subprocess
import sys
from collections import defaultdict

# (column, log file under logs/<stage>/<sample>/, regex) — first match wins.
CASCADE = [
    ("pass_snv", "filtering", "accessibility_filter.log", r"variants_input=(\d+)"),
    ("accessible", "filtering", "accessibility_filter.log", r"variants_kept=(\d+)"),
    ("germline", "filtering", "germline_filter.log", r"variants_kept=(\d+)"),
    ("vaf", "filtering", "vaf_filter.log", r"variants_kept=(\d+)"),
    ("cnvnator", "filtering", "cnvnator_filter.log", r"cnvnator: \d+ -> (\d+)"),
    ("mosaicforecast", "filtering", "mosaicforecast_filter.log", r"mosaicforecast: \d+ -> (\d+)"),
    ("pon", "filtering", "pon_mask_filter.log", r"variants_kept=(\d+)"),
]
MF_MODEL_RE = r"sample=\S+  model=(\S+)"
CONTAM_RE = r"contamination_value=([0-9.eE+-]+)"
# One log per stage that is written exactly once per sample — the host-attribution key.
STAGE_LOGS = {
    "calling": ("calling", "filter_mutect_calls.log"),
    "filtering": ("filtering", "pon_mask_filter.log"),
}


def load_expected(path: str) -> dict[str, set[str]]:
    """barcode_map.tsv (sample_id, fc, udb) -> {sample: {"<fc>_<udb>", ...}}."""
    expected: dict[str, set[str]] = defaultdict(set)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            expected[row["sample_id"]].add(f"{row['fc']}_{row['udb']}")
    return dict(expected)


def parse_flagstat(path: str) -> dict[str, float]:
    out: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"(\d+) \+ \d+ (.+?)(?: \(|$)", line.strip())
            if m:
                out.setdefault(m.group(2).strip(), int(m.group(1)))
    return out


def cram_readgroups(cram: str, samtools: list[str]) -> list[tuple[str, str]]:
    """[(PU, SM), ...] from the CRAM header."""
    res = subprocess.run(samtools + ["view", "-H", cram], capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(res.stderr.strip()[:200])
    rgs = []
    for line in res.stdout.splitlines():
        if line.startswith("@RG"):
            tags = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            rgs.append((tags.get("PU", ""), tags.get("SM", "")))
    return rgs


def vcf_samples(vcf: str) -> list[str]:
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    return []


def grep(path: str, pattern: str) -> str | None:
    try:
        with open(path, errors="replace") as fh:
            m = re.search(pattern, fh.read())
    except OSError:
        return None
    return m.group(1) if m else None


def count_lines(path: str) -> int:
    with open(path) as fh:
        return sum(1 for line in fh if line.strip() and not line.startswith("#"))


def mtime(path: str) -> float | None:
    try:
        return os.path.getmtime(path)
    except OSError:
        return None


def audit_sample(s: str, rgs_expected: set[str], args, log_dirs) -> dict:
    row: dict = {"sample": s, "n_rg_expected": len(rgs_expected)}
    problems: list[str] = []
    mdir = os.path.join(args.results, "mapping", s)
    cram, crai = os.path.join(mdir, f"{s}.cram"), os.path.join(mdir, f"{s}.cram.crai")
    flagstat = os.path.join(mdir, "flagstat.txt")
    vcf = os.path.join(args.results, "calling", s, f"{s}.filtered.vcf.gz")
    final = os.path.join(args.results, "filtering", s, f"{s}.final.txt")

    # --- mapping ---
    row["mapped_ok"] = all(os.path.exists(p) for p in (cram, crai, flagstat))
    if not row["mapped_ok"]:
        problems.append("mapping_incomplete")
    else:
        fs = parse_flagstat(flagstat)
        primary = fs.get("primary") or fs.get("in total") or 0
        row["primary_reads"] = primary
        if primary:
            row["mapped_pct"] = round(
                100 * fs.get("primary mapped", fs.get("mapped", 0)) / primary, 3
            )
            row["dup_pct"] = round(
                100 * fs.get("primary duplicates", fs.get("duplicates", 0)) / primary, 3
            )
            row["est_depth"] = round(
                (fs.get("primary mapped", 0) - fs.get("primary duplicates", 0))
                * args.read_length
                / args.genome_size,
                1,
            )
        if args.samtools:
            try:
                rgs = cram_readgroups(cram, args.samtools)
                pus = {pu for pu, _ in rgs}
                sms = {sm for _, sm in rgs}
                row["n_rg_cram"] = len(rgs)
                if pus != rgs_expected:
                    problems.append(
                        f"rg_mismatch(missing={sorted(rgs_expected - pus)},"
                        f"extra={sorted(pus - rgs_expected)})"
                    )
                if sms != {s}:
                    problems.append(f"rg_sm={','.join(sorted(sms))}")
            except RuntimeError as exc:
                problems.append(f"cram_header_error:{exc}")

    # --- calling ---
    row["called_ok"] = os.path.exists(vcf) and os.path.exists(vcf + ".tbi")
    if not row["called_ok"]:
        problems.append("calling_incomplete")
    else:
        names = vcf_samples(vcf)
        if names != [s]:
            problems.append(f"vcf_sample={','.join(names) or 'none'}")

    # --- filtering ---
    row["filtered_ok"] = os.path.exists(final)
    if row["filtered_ok"]:
        row["final_snv"] = count_lines(final)
    else:
        problems.append("filtering_incomplete")

    # --- staleness: each output must be newer than the input it was built from ---
    for older, newer, tag in (
        (cram, vcf, "vcf_older_than_cram"),
        (vcf, final, "final_older_than_vcf"),
    ):
        a, b = mtime(older), mtime(newer)
        if a and b and b < a:
            problems.append(tag)

    # --- host attribution + once-only ---
    for stage, (sub, fname) in STAGE_LOGS.items():
        hosts = [lbl for lbl, d in log_dirs if os.path.exists(os.path.join(d, sub, s, fname))]
        row[f"{stage}_host"] = "+".join(hosts) or "none"
        if len(hosts) > 1:
            problems.append(f"{stage}_ran_on_{'+'.join(hosts)}")
    fhost = row["filtering_host"].split("+")[0]
    chost = row["calling_host"].split("+")[0]
    fdir = dict(log_dirs).get(fhost)
    cdir = dict(log_dirs).get(chost)

    # --- cascade: filtering host's logs first, then any other host (a sample whose
    # early steps ran in a separate pilot clone has them only there) ---
    search = [fdir] + [d for _, d in log_dirs if d != fdir] if fdir else [d for _, d in log_dirs]
    for col, sub, fname, pat in CASCADE:
        v = next((x for d in search if (x := grep(os.path.join(d, sub, s, fname), pat))), None)
        row[col] = int(v) if v is not None else ""
    if fdir:
        row["mf_model"] = (
            grep(os.path.join(fdir, "filtering", s, "mosaicforecast_filter.log"), MF_MODEL_RE) or ""
        )
        row["mf_model"] = os.path.basename(row["mf_model"])
        if not grep(
            os.path.join(fdir, "filtering", s, "pon_mask_filter.log"), r"(END) pon_mask_filter"
        ):
            problems.append("pon_log_no_END")
    if cdir:
        row["contamination"] = (
            grep(os.path.join(cdir, "calling", s, "calculate_contamination.log"), CONTAM_RE) or ""
        )
    if (
        row.get("pon") != ""
        and row.get("final_snv") is not None
        and row.get("pon") != row["final_snv"]
    ):
        problems.append(f"final_lines({row['final_snv']})!=pon_log({row['pon']})")

    row["problems"] = ";".join(problems)
    return row


def stray_dirs(results: str, expected: set[str]) -> dict[str, list[str]]:
    out = {}
    for stage in ("mapping", "calling", "filtering"):
        d = os.path.join(results, stage)
        names = set(os.listdir(d)) if os.path.isdir(d) else set()
        out[stage] = sorted(n for n in names - expected if not n.startswith((".", "_")))
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--barcode-map", required=True)
    p.add_argument("--results", default="results")
    p.add_argument(
        "--log-dir",
        action="append",
        required=True,
        metavar="LABEL=DIR",
        help="host log dir, repeatable (e.g. vm100=logs)",
    )
    p.add_argument("--samtools", default="", help="samtools command for @RG check (omit to skip)")
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--genome-size", type=int, default=3_099_734_149)
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    args.samtools = shlex.split(args.samtools)
    log_dirs = [tuple(spec.split("=", 1)) for spec in args.log_dir]

    expected = load_expected(args.barcode_map)
    rows = []
    for i, s in enumerate(sorted(expected), 1):
        rows.append(audit_sample(s, expected[s], args, log_dirs))
        print(f"\r[{i}/{len(expected)}] {s}", end="", file=sys.stderr)
    print(file=sys.stderr)

    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # --- console summary ---
    n = len(rows)
    print(f"expected samples      : {n} (readgroups {sum(len(v) for v in expected.values())})")
    for key, label in (
        ("mapped_ok", "mapping"),
        ("called_ok", "calling"),
        ("filtered_ok", "filtering"),
    ):
        print(f"{label:<22}: {sum(bool(r[key]) for r in rows)} / {n}")
    for stage in STAGE_LOGS:
        tally = defaultdict(int)
        for r in rows:
            tally[r[f"{stage}_host"]] += 1
        print(f"{stage} host split   : " + ", ".join(f"{k}={v}" for k, v in sorted(tally.items())))
    bad = [r for r in rows if r["problems"]]
    print(f"samples with problems : {len(bad)}")
    for r in bad[:40]:
        print(f"  {r['sample']}: {r['problems']}")
    for stage, names in stray_dirs(args.results, set(expected)).items():
        if names:
            print(f"unexpected dirs in results/{stage}: {len(names)} -> {' '.join(names[:15])}")
    print("cascade (median [min-max]):")
    for col, *_ in CASCADE:
        vals = [r[col] for r in rows if r.get(col) != ""]
        if vals:
            print(
                f"  {col:<15} n={len(vals):<4} {statistics.median(vals):>10,.0f} "
                f"[{min(vals):,}-{max(vals):,}]  total={sum(vals):,}"
            )
    print(f"wrote {args.out}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
