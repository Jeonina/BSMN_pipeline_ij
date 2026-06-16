#!/usr/bin/env python3
"""CNVnator genotype filter — BSMN D.CNVnator_genotype_filter.sh step.

Faithful re-implementation: genotype a +/-window region around each candidate
with CNVnator and drop candidates whose estimated copy number is >= a threshold
(2.5), i.e. in a duplicated/amplified region where low apparent VAF is a
copy-number / paralog artifact rather than true mosaicism.

The sample ROOT file (read-depth histogram) is built upstream by the
cnvnator_root rule. This script only runs `cnvnator -genotype` (interactive,
regions on stdin) via the CNVnator Apptainer image and parses the output.

BSMN parsed `$9 < 2.5` from `paste CAND(5col) genotype`; that maps to the 4th
whitespace field of the genotype line (0-based index 3) = the CN estimate. The
field index and threshold are configurable, and raw genotype lines are logged
so the column can be verified on first run.

Genotype-line parsing is a pure function (unit-tested); the container call is
not exercised without a ROOT file.

Usage:
    python scripts/cnvnator_filter.py \
        --candidates results/filtering/SM/SM.vaf_filtered.txt \
        --root results/filtering/SM/cnvnator/SM.root \
        --cnvnator-sif containers/cnvnator.sif \
        --binsize 100 --cn-threshold 2.5 \
        > results/filtering/SM/SM.cnvnator_filtered.txt
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys

log = logging.getLogger("cnvnator_filter")

# Genotype-line field holding the CN estimate (0-based), matching BSMN's $9
# over a 5-column candidate paste.
CN_FIELD = 3


def read_candidates(path: str) -> list[tuple[str, int, str, str]]:
    out: list[tuple[str, int, str, str]] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            out.append((p[0], int(p[1]), p[2].upper(), p[3].upper()))
    return out


def region_str(chrom: str, pos: int, window: int) -> str:
    """+/-window region around a 1-based position (start clamped to 1)."""
    return f"{chrom}:{max(1, pos - window)}-{pos + window}"


def parse_genotype(text: str, cn_field: int = CN_FIELD) -> dict[str, float]:
    """Map region -> CN estimate from `cnvnator -genotype` output.

    Each emitted line looks like: ``Genotype <region> <...> <CN> ...``.
    Lines without a parseable CN at cn_field are skipped.
    """
    out: dict[str, float] = {}
    for line in text.splitlines():
        if "Genotype" not in line:
            continue
        f = line.split()
        # locate the 'Genotype' keyword; region is the next field
        try:
            gi = f.index("Genotype")
        except ValueError:
            continue
        if len(f) <= gi + cn_field:
            continue
        region = f[gi + 1]
        try:
            cn = float(f[gi + cn_field])
        except ValueError:
            continue
        out[region] = cn
    return out


def run(args: argparse.Namespace) -> None:
    cands = read_candidates(args.candidates)
    log.info("candidates in=%d  root=%s  binsize=%d", len(cands), args.root, args.binsize)
    if not cands:
        return

    regions = [region_str(c, p, args.window) for c, p, _, _ in cands]
    stdin = "\n".join(regions) + "\nexit\n"

    cmd = [
        "apptainer",
        "exec",
        args.cnvnator_sif,
        "cnvnator",
        "-root",
        args.root,
        "-genotype",
        str(args.binsize),
    ]
    proc = subprocess.run(cmd, input=stdin, capture_output=True, text=True)
    if proc.returncode != 0:
        log.error(
            "cnvnator -genotype failed (rc=%d): %s", proc.returncode, (proc.stderr or "")[-500:]
        )
        raise SystemExit(1)

    cn = parse_genotype(proc.stdout, args.cn_field)
    # Log a few raw lines so the CN column can be verified.
    sample_lines = [ln for ln in proc.stdout.splitlines() if "Genotype" in ln][:3]
    for ln in sample_lines:
        log.info("raw genotype: %s", ln.strip())

    kept = removed = missing = 0
    for c, p, ref, alt in cands:
        reg = region_str(c, p, args.window)
        val = cn.get(reg)
        if val is None:
            # No genotype for this region — keep (sensitivity-first) and warn.
            missing += 1
            log.warning("no genotype for %s — keeping", reg)
            sys.stdout.write(f"{c}\t{p}\t{ref}\t{alt}\n")
            kept += 1
        elif val < args.cn_threshold:
            sys.stdout.write(f"{c}\t{p}\t{ref}\t{alt}\n")
            kept += 1
        else:
            removed += 1
    sys.stdout.flush()
    log.info(
        "cnvnator: %d -> %d  (removed CN>=%.1f: %d, missing/kept: %d)",
        len(cands),
        kept,
        args.cn_threshold,
        removed,
        missing,
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="CNVnator genotype filter (BSMN D-step).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--candidates", required=True, help="input candidates (chr pos ref alt)")
    p.add_argument("--root", required=True, help="sample CNVnator ROOT file")
    p.add_argument("--cnvnator-sif", required=True, help="CNVnator Apptainer .sif")
    p.add_argument("--binsize", type=int, default=100, help="CNVnator bin size [100]")
    p.add_argument("--cn-threshold", type=float, default=2.5, help="drop if CN >= this [2.5]")
    p.add_argument("--window", type=int, default=1000, help="+/-bp window per candidate [1000]")
    p.add_argument("--cn-field", type=int, default=CN_FIELD, help=f"CN field index [{CN_FIELD}]")
    p.set_defaults(func=run)
    return p


def main() -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(name)s] %(levelname)s: %(message)s"
    )
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
