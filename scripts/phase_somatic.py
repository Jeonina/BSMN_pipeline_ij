#!/usr/bin/env python3
"""Read-backed phasing of the final somatic SNVs (MosaicForecast Phase.py).

Bae et al. 2022 (Science) checked hypermutable brains by phasing each call to
nearby germline het SNPs: a real mosaic SNV sits on ONE parental haplotype,
so reads show three haplotypes (ref-A, ref-B, alt-B; "hap=3"); a germline het
shows two ("hap=2"); alt reads on both parental haplotypes ("hap>3") point to
an artifact (sequencing/mapping error, CNV).  Phasing only tests calls with an
informative het SNP within read-pair reach, so most calls stay unphased.

Phasing does NOT separate contamination from mosaicism: a contaminant's
germline allele also rides on one of its own haplotypes, so it can look like
hap=3.  The contaminated (QC tier A) samples are kept in the run as that
reference.  This is a validation readout, not a filter.

Runs Phase.py once over every sample's final.txt, inside the same
yanmei/mosaicforecast image the E-step uses (scripts + k24 bigwig in
/usr/local/bin), with htslib pointed at the local REF_CACHE.  Phase.py wants
one directory holding {sample}.cram, so this builds a symlink directory.

Usage (VM100, pipeline directory):
    python ~/phase_somatic.py --sample-info ~/sample_info.tsv --results results \\
        --ref resources/hg38/Homo_sapiens_assembly38.fasta \\
        --ref-cache resources/hg38/ref_cache \\
        --mf-sif containers/mosaicforecast.sif --bind /mnt/data \\
        --out-dir /mnt/data/bsmn/phasing --threads 8
    add --dry-run to write the inputs and print the command only.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys

PHASE_SCRIPT = "/usr/local/bin/Phase.py"
UMAP_BW = "resources/hg38/k24.Umap.MultiTrackMappability.bw"   # image's own track is hg19 (output column only)


def read_final(path: str) -> list[tuple[str, int, str, str]]:
    """final.txt: CHROM POS REF ALT (extra columns ignored); SNVs only."""
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[1].isdigit() and len(f[2]) == 1 and len(f[3]) == 1:
                out.append((f[0], int(f[1]), f[2], f[3]))
    return out


def position_lines(calls: dict[str, list[tuple[str, int, str, str]]]) -> list[str]:
    """Phase.py input: chrom  pos-1  pos  ref  alt  sample, sorted by sample then position."""
    lines = []
    for s in sorted(calls):
        for c, p, r, a in sorted(calls[s], key=lambda x: (x[0], x[1])):
            lines.append(f"{c}\t{p - 1}\t{p}\t{r}\t{a}\t{s}")
    return lines


def mf_env(ref_cache: str | None) -> dict[str, str]:
    """Same as mosaicforecast_filter.mf_env: local REF_CACHE so CRAM decoding never hits EBI."""
    env = dict(os.environ)
    if ref_cache:
        tmpl = os.path.join(os.path.realpath(ref_cache), "%2s", "%2s", "%s")
        env["APPTAINERENV_REF_PATH"] = tmpl
        env["APPTAINERENV_REF_CACHE"] = tmpl
    return env


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-info", required=True, help="TSV with a sample_id column")
    ap.add_argument("--results", default="results")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--ref-cache", default=None)
    ap.add_argument("--mf-sif", required=True)
    ap.add_argument("--bind", default=None, help="extra apptainer bind, e.g. /mnt/data")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-dp-informative", type=int, default=20)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    with open(args.sample_info) as fh:
        samples = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]
    bam_dir = os.path.join(args.out_dir, "crams")
    phase_dir = os.path.join(args.out_dir, "phase")
    os.makedirs(bam_dir, exist_ok=True)
    os.makedirs(phase_dir, exist_ok=True)

    calls, missing = {}, []
    for s in samples:
        final = os.path.join(args.results, "filtering", s, f"{s}.final.txt")
        cram = os.path.realpath(os.path.join(args.results, "mapping", s, f"{s}.cram"))
        if not (os.path.exists(final) and os.path.exists(cram) and os.path.exists(cram + ".crai")):
            missing.append(s)
            continue
        calls[s] = read_final(final)
        # Phase.py's last step opens "{sample}.bam" even in cram mode (mode "rc"), so the
        # CRAM is linked under both names; htslib finds the index as "<name>.crai".
        for src, dst in ((cram, f"{s}.cram"), (cram + ".crai", f"{s}.cram.crai"),
                         (cram, f"{s}.bam"), (cram + ".crai", f"{s}.bam.crai")):
            link = os.path.join(bam_dir, dst)
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(src, link)
    positions = os.path.join(args.out_dir, "input_positions.txt")
    lines = position_lines(calls)
    with open(positions, "w") as fo:
        fo.write("\n".join(lines) + "\n")
    print(f"{len(calls)} samples, {len(lines)} SNVs -> {positions}; missing: {', '.join(missing) or 'none'}", flush=True)

    cmd = ["apptainer", "exec"] + (["-B", args.bind] if args.bind else []) + [
        args.mf_sif, "python3", PHASE_SCRIPT, bam_dir, phase_dir, os.path.realpath(args.ref), positions,
        str(args.min_dp_informative), UMAP_BW, str(args.threads), "cram"]
    print("command:", " ".join(cmd), flush=True)
    if args.dry_run:
        return 0
    rc = subprocess.run(cmd, env=mf_env(args.ref_cache)).returncode
    print(f"Phase.py exit {rc}; outputs: {', '.join(sorted(os.listdir(phase_dir)))}", flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
