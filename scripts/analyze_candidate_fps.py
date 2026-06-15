#!/usr/bin/env python3
"""Characterize final mosaic candidates against Mutect2 annotations.

Given the final mosaic candidate list (``chr pos ref alt``) and the Mutect2
``*.filtered.vcf.gz`` they were derived from, pull the per-site annotations
(strand table, read-position, mapping/base quality, VAF, depth) and report
which false-positive heuristics would remove each candidate while preserving
a known true-positive coordinate.

This is a *diagnostic* tool, not a pipeline filter. It answers the question
"why did these candidates survive, and what extra signal separates the
artifacts from the real mosaic?" before any new filter is wired into the
cascade (``workflow/rules/filtering.smk``).

Strand and read-position signals come straight from the Mutect2 INFO/FORMAT
fields, so no BAM access is required. For the BAM-based confirmation, run
``scripts/strand_bias_filter.py`` and ``scripts/repeat_filter.py`` on the
same candidate list.

Usage:
    python scripts/analyze_candidate_fps.py \
        --vcf results/calling/HG002/HG002.filtered.vcf.gz \
        --candidates reports/benchmark_2026-06-02/HG002.final.txt \
        --truth chr20:47053475 \
        --bcftools-sif containers/bcftools_1.17.sif
"""

from __future__ import annotations

import argparse
import subprocess
import sys


def fisher_strand_p(rf: int, rr: int, af: int, ar: int) -> float | None:
    """Fisher exact p on the ref/alt x fwd/rev table. None if scipy absent.

    scipy ships in the Apptainer container (see strand_bias_filter.py) but is
    optional locally; the rest of the diagnostic works without it.
    """
    try:
        from scipy.stats import fisher_exact
    except ImportError:
        return None
    return float(fisher_exact([[rf, af], [rr, ar]])[1])


def run_bcftools(args_list: list[str], sif: str | None) -> str:
    """Run bcftools, optionally inside an Apptainer container (repo convention)."""
    if sif:
        cmd = ["apptainer", "exec", sif, "bcftools"] + args_list
    else:
        cmd = ["bcftools"] + args_list
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return str(result.stdout)


def read_candidates(path: str) -> list[tuple[str, int, str, str]]:
    """Parse a ``chr pos ref alt`` candidate list (extra columns ignored)."""
    out: list[tuple[str, int, str, str]] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
            out.append((chrom, pos, ref.upper(), alt.upper()))
    return out


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in info.split(";"):
        if "=" in field:
            key, val = field.split("=", 1)
            out[key] = val
        else:
            out[field] = ""
    return out


def alt_value(field: str | None) -> float | None:
    """Mutect2 MBQ/MMQ are 'ref,alt'; MPOS is alt-only. Take the last value."""
    if field is None or field == "" or field == ".":
        return None
    try:
        return float(field.split(",")[-1])
    except ValueError:
        return None


def strand_counts(info: dict[str, str], fmt: dict[str, str]) -> tuple[int, int, int, int] | None:
    """Return (ref_fwd, ref_rev, alt_fwd, alt_rev) from AS_SB_TABLE or FORMAT/SB."""
    table = info.get("AS_SB_TABLE")
    if table and "|" in table:
        ref_part, alt_part = table.split("|")[:2]
        rf, rr = (int(x) for x in ref_part.split(","))
        af, ar = (int(x) for x in alt_part.split(","))
        return rf, rr, af, ar
    sb = fmt.get("SB")
    if sb:
        vals = [int(x) for x in sb.split(",")]
        if len(vals) == 4:
            return vals[0], vals[1], vals[2], vals[3]
    return None


def vaf_and_depth(fmt: dict[str, str]) -> tuple[float, int, int]:
    """Return (vaf, alt_count, depth) from FORMAT AF/AD/DP."""
    alt_n = 0
    depth = 0
    if "AD" in fmt:
        ad = [int(x) for x in fmt["AD"].split(",")]
        alt_n = ad[1] if len(ad) > 1 else 0
        depth = sum(ad)
    if "DP" in fmt:
        depth = int(fmt["DP"]) or depth
    if "AF" in fmt:
        try:
            vaf = float(fmt["AF"].split(",")[0])
        except ValueError:
            vaf = alt_n / depth if depth else 0.0
    else:
        vaf = alt_n / depth if depth else 0.0
    return vaf, alt_n, depth


def collect(vcf: str, sif: str | None, want: set[tuple[str, int]]) -> dict[tuple[str, int], dict]:
    """Pull annotations for the wanted (chrom, pos) sites from the VCF."""
    output = run_bcftools(["view", "-H", vcf], sif)
    rows: dict[tuple[str, int], dict] = {}
    for line in output.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 10:
            continue
        chrom, pos = f[0], int(f[1])
        if (chrom, pos) not in want:
            continue
        info = parse_info(f[7])
        fmt = dict(zip(f[8].split(":"), f[9].split(":"), strict=False))
        vaf, alt_n, depth = vaf_and_depth(fmt)
        sc = strand_counts(info, fmt)
        p_strand = None
        alt_fwd = alt_rev = None
        if sc is not None:
            rf, rr, af, ar = sc
            alt_fwd, alt_rev = af, ar
            p_strand = fisher_strand_p(rf, rr, af, ar)
        rows[(chrom, pos)] = {
            "filter": f[6],
            "vaf": vaf,
            "alt_n": alt_n,
            "depth": depth,
            "mpos": alt_value(info.get("MPOS")),
            "mmq_alt": alt_value(info.get("MMQ")),
            "mbq_alt": alt_value(info.get("MBQ")),
            "popaf": alt_value(info.get("POPAF")),
            "tlod": alt_value(info.get("TLOD")),
            "alt_fwd": alt_fwd,
            "alt_rev": alt_rev,
            "p_strand": p_strand,
        }
    return rows


# Heuristic thresholds — tunable; defaults from BSMN-classic practice.
HEURISTICS = {
    "strand_bias (Fisher p<0.01)": lambda r: r["p_strand"] is not None and r["p_strand"] < 0.01,
    "single-strand alt": lambda r: r["alt_fwd"] == 0 or r["alt_rev"] == 0,
    "low MPOS (<10)": lambda r: r["mpos"] is not None and r["mpos"] < 10,
    "low alt MMQ (<40)": lambda r: r["mmq_alt"] is not None and r["mmq_alt"] < 40,
    "low alt MBQ (<25)": lambda r: r["mbq_alt"] is not None and r["mbq_alt"] < 25,
    "min_alt raise (alt<8)": lambda r: r["alt_n"] < 8,
}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="Mutect2 *.filtered.vcf.gz")
    p.add_argument("--candidates", required=True, help="final mosaic list (chr pos ref alt)")
    p.add_argument("--truth", default="", help="known true-positive coord, e.g. chr20:47053475")
    p.add_argument("--bcftools-sif", default=None, help="Apptainer .sif for bcftools (optional)")
    args = p.parse_args()

    cands = read_candidates(args.candidates)
    want = {(c, pos) for c, pos, _, _ in cands}
    rows = collect(args.vcf, args.bcftools_sif, want)

    truth_key: tuple[str, int] | None = None
    if args.truth:
        tc, tp = args.truth.split(":")
        truth_key = (tc, int(tp))

    header = (
        f"{'site':<22}{'TP':<4}{'VAF':>7}{'alt/DP':>10}"
        f"{'MPOS':>6}{'MMQ':>6}{'MBQ':>6}{'altF/R':>9}{'pStr':>9}"
    )
    print(header)
    print("-" * len(header))
    missing: list[str] = []
    for chrom, pos, ref, alt in cands:
        r = rows.get((chrom, pos))
        site = f"{chrom}:{pos} {ref}>{alt}"
        if r is None:
            missing.append(site)
            continue
        is_tp = "TP" if (chrom, pos) == truth_key else ""
        ps = f"{r['p_strand']:.1e}" if r["p_strand"] is not None else "."
        afr = f"{r['alt_fwd']}/{r['alt_rev']}" if r["alt_fwd"] is not None else "."
        mpos = f"{r['mpos']:.0f}" if r["mpos"] is not None else "."
        mmq = f"{r['mmq_alt']:.0f}" if r["mmq_alt"] is not None else "."
        mbq = f"{r['mbq_alt']:.0f}" if r["mbq_alt"] is not None else "."
        alt_dp = f"{r['alt_n']}/{r['depth']}"
        print(
            f"{site:<22}{is_tp:<4}{r['vaf']:>7.3f}{alt_dp:>10}"
            f"{mpos:>6}{mmq:>6}{mbq:>6}{afr:>9}{ps:>9}"
        )

    if missing:
        print(f"\n[warn] {len(missing)} candidate(s) not found in VCF:", file=sys.stderr)
        for s in missing:
            print(f"  {s}", file=sys.stderr)

    # Heuristic effectiveness: how many FPs removed vs TP preserved.
    fp_rows = [
        rows[(c, pos)] for c, pos, _, _ in cands
        if (c, pos) in rows and (c, pos) != truth_key
    ]
    tp_row = rows.get(truth_key) if truth_key else None
    n_fp = len(fp_rows)
    print(f"\nHeuristic effectiveness  (FP pool = {n_fp}" + (", TP tracked)" if tp_row else ")"))
    print(f"{'heuristic':<30}{'FP removed':>12}{'TP removed?':>14}")
    print("-" * 56)
    for name, fn in HEURISTICS.items():
        removed = sum(1 for r in fp_rows if fn(r))
        tp_hit = "YES (bad)" if tp_row and fn(tp_row) else ("no" if tp_row else "n/a")
        print(f"{name:<30}{f'{removed}/{n_fp}':>12}{tp_hit:>14}")


if __name__ == "__main__":
    main()
