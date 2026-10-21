#!/usr/bin/env python3
"""Merge each donor's DLPFC + Caudate CRAMs into one donor CRAM for joint-depth calling.

Lab-meeting decision (2026-10-21): call somatic SNVs per donor on the pooled
reads of both brain regions (~2x depth), then go back to the per-region CRAMs
to label each call shared / DLPFC-only / Caudate-only.  This script makes the
donor CRAMs and a bam-schema sample sheet so the pipeline runs on them with
ZERO rule changes (rules/ingest.smk links a CRAM whose @RG already carries
SM == sample_id).

Per donor (idempotent; a finished donor is skipped):
  1. samtools merge  <Br>_DLPFC.cram + <Br>_Caudate.cram  -> tmp CRAM
     (read groups stay distinct: mapping wrote ID:<D-number>.<fc>_<udb>)
  2. rewrite every @RG SM:<D-number> to SM:<Br>  (samtools reheader, no re-encode)
  3. index + quickcheck, then rename into place

Inputs come from sample_info.tsv (sample_id, br_region = "<Br>_<Region>",
exported from the SampleInfo sheet of the barcode workbook); each donor must
have exactly one DLPFC and one Caudate sample.

Usage (on VM100, from the pipeline directory):
    python scripts/merge_donor_crams.py \\
        --sample-info sample_info.tsv --results results \\
        --out-dir /mnt/data/bsmn/donor_merge/cram \\
        --ref resources/hg38/Homo_sapiens_assembly38.fasta \\
        --samtools "apptainer exec containers/samtools_1.17.sif samtools" \\
        --jobs 4 --threads 8 \\
        --sheet-out config/samples.donor_merge.tsv
Add --dry-run first: it checks inputs and disk space and merges nothing.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

REGIONS = ("DLPFC", "Caudate")


def load_donors(path: str) -> dict[str, dict[str, str]]:
    """sample_info.tsv -> {Br: {"DLPFC": D-number, "Caudate": D-number}}."""
    donors: dict[str, dict[str, str]] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            donor, _, region = row["br_region"].rpartition("_")
            if region not in REGIONS or not donor:
                raise ValueError(f"{row['sample_id']}: unexpected br_region {row['br_region']!r}")
            slot = donors.setdefault(donor, {})
            if region in slot:
                raise ValueError(f"{donor}: two {region} samples ({slot[region]}, {row['sample_id']})")
            slot[region] = row["sample_id"]
    incomplete = sorted(d for d, s in donors.items() if len(s) != len(REGIONS))
    if incomplete:
        raise ValueError(f"donors without both regions: {', '.join(incomplete)}")
    return donors


def rewrite_sm(header: str, donor: str, expected: set[str]) -> str:
    """Set SM:<donor> on every @RG; refuse a header whose SMs are not exactly `expected`."""
    seen: set[str] = set()
    out = []
    for line in header.splitlines():
        if line.startswith("@RG"):
            fields = line.split("\t")
            for i, f in enumerate(fields):
                if f.startswith("SM:"):
                    seen.add(f[3:])
                    fields[i] = f"SM:{donor}"
            line = "\t".join(fields)
        out.append(line)
    if seen != expected:
        raise ValueError(f"{donor}: @RG SM set {sorted(seen)} != expected {sorted(expected)}")
    return "\n".join(out) + "\n"


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    res = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if res.returncode != 0:
        raise RuntimeError(f"{shlex.join(cmd)}\n{res.stderr.strip()[-500:]}")
    return res


def merge_donor(donor: str, inputs: list[str], args) -> str:
    st = shlex.split(args.samtools)
    ddir = os.path.join(args.out_dir, donor)
    final = os.path.join(ddir, f"{donor}.cram")
    os.makedirs(ddir, exist_ok=True)
    merged, reheaded = (os.path.join(ddir, f".{donor}.{t}.cram") for t in ("merged", "reheaded"))
    t0 = time.time()

    run(st + ["merge", "-f", "-@", str(args.threads), "--reference", args.ref,
              "-O", "CRAM", "-o", merged] + inputs)
    header = run(st + ["view", "-H", merged]).stdout
    expected = {os.path.basename(p).removesuffix(".cram") for p in inputs}
    hdr_path = os.path.join(ddir, f".{donor}.header.sam")
    with open(hdr_path, "w") as fh:
        fh.write(rewrite_sm(header, donor, expected))
    with open(reheaded, "wb") as fh:
        res = subprocess.run(st + ["reheader", hdr_path, merged], stdout=fh, stderr=subprocess.PIPE)
    if res.returncode != 0:
        raise RuntimeError(f"reheader {donor}: {res.stderr.decode()[-500:]}")
    run(st + ["quickcheck", reheaded])
    run(st + ["index", "-@", str(args.threads), reheaded])

    os.replace(reheaded + ".crai", final + ".crai")
    os.replace(reheaded, final)   # CRAM last: its presence + .crai marks the donor done
    for p in (merged, hdr_path):
        os.remove(p)
    return f"{donor}\tmerged\t{(time.time() - t0) / 60:.0f} min"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-info", required=True, help="TSV: sample_id, br_region")
    ap.add_argument("--results", default="results", help="pipeline results/ (mapping/<D>/<D>.cram)")
    ap.add_argument("--out-dir", required=True, help="donor CRAMs go to <out-dir>/<Br>/<Br>.cram")
    ap.add_argument("--ref", required=True, help="reference FASTA the CRAMs were written against")
    ap.add_argument("--samtools", default="samtools", help='e.g. "apptainer exec samtools.sif samtools"')
    ap.add_argument("--jobs", type=int, default=4, help="donors merged at once")
    ap.add_argument("--threads", type=int, default=8, help="samtools threads per donor")
    ap.add_argument("--donors", nargs="*", help="only these Br IDs (default: all)")
    ap.add_argument("--sheet-out", default="config/samples.donor_merge.tsv")
    ap.add_argument("--dry-run", action="store_true", help="check inputs and disk space only")
    args = ap.parse_args(argv)

    donors = load_donors(args.sample_info)
    if args.donors:
        unknown = set(args.donors) - donors.keys()
        if unknown:
            ap.error(f"unknown donors: {', '.join(sorted(unknown))}")
        donors = {d: donors[d] for d in args.donors}

    todo: dict[str, list[str]] = {}
    done, missing = [], []
    for donor, slot in sorted(donors.items()):
        final = os.path.join(args.out_dir, donor, f"{donor}.cram")
        if os.path.exists(final) and os.path.exists(final + ".crai"):
            done.append(donor)
            continue
        inputs = [os.path.join(args.results, "mapping", slot[r], f"{slot[r]}.cram") for r in REGIONS]
        absent = [p for p in inputs if not (os.path.exists(p) and os.path.exists(p + ".crai"))]
        if absent:
            missing.append(f"{donor}: {', '.join(absent)}")
        else:
            todo[donor] = inputs

    # Each donor CRAM ~ sum of its inputs; the merged and reheaded copies coexist
    # briefly, so in-flight donors need 2x.
    need = sum(os.path.getsize(p) for ins in todo.values() for p in ins)
    peak = sorted((sum(map(os.path.getsize, ins)) for ins in todo.values()), reverse=True)[: args.jobs]
    need += sum(peak)
    os.makedirs(args.out_dir, exist_ok=True)
    free = shutil.disk_usage(args.out_dir).free
    print(f"donors: {len(donors)}  done: {len(done)}  to merge: {len(todo)}  missing input: {len(missing)}")
    print(f"space : need ~{need / 1e12:.2f} TB, free {free / 1e12:.2f} TB on {args.out_dir}")
    for m in missing:
        print(f"  MISSING {m}")
    if need > free:
        print("ERROR: not enough free space", file=sys.stderr)
        return 1
    if args.dry_run:
        return 0

    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(merge_donor, d, ins, args): d for d, ins in todo.items()}
        for fut in as_completed(futs):
            try:
                print(fut.result(), flush=True)
                done.append(futs[fut])
            except Exception as exc:  # report every donor, then fail at the end
                failed.append(futs[fut])
                print(f"{futs[fut]}\tFAILED\t{exc}", flush=True)

    with open(args.sheet_out, "w") as fh:
        fh.write("sample_id\tbam\n")
        for donor in sorted(done):
            fh.write(f"{donor}\t{os.path.abspath(os.path.join(args.out_dir, donor, donor + '.cram'))}\n")
    print(f"sheet : {args.sheet_out} ({len(done)} donors)")
    if failed or missing:
        print(f"incomplete: {len(failed)} failed, {len(missing)} missing input", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
