#!/usr/bin/env python3
"""Build config/samples.tsv from a barcode map + the FASTQs on disk.

The Hypermutability cohort multiplexes each biological sample (a D-number) with
TWO barcodes per flowcell, and REUSES barcodes across flowcells — so a (FC, UDB)
pair, not the UDB alone, identifies a readgroup, and the true sample is the
D-number. Grouping by UDB would merge different samples.

Input: a TSV mapping each (sample_id, fc, udb) triple — one row per readgroup:

    sample_id   fc            udb
    D19508      E250080627    UDB-385
    D19508      E250080627    UDB-386
    D19508      E250084677    UDB-497
    ...

For each triple this locates <data-root>/**/<fc>/<fc>_L01_<udb>_{1,2}.fq.gz and
emits one sample sheet row (readgroup = "<fc>_<udb>"). All readgroups of a
sample_id are merged into one CRAM downstream.

Usage:
    python scripts/build_cohort_samples.py --map barcode_map.tsv \
        --data-root /mnt/data -o config/samples.tsv
"""

import argparse
import os
import re
import subprocess
import sys


def build_fastq_index(data_root: str) -> dict[tuple[str, str], str]:
    """Index every real UDB R1 FASTQ under data_root by (fc, udb).

    One `find` pass (metadata only). Empty/placeholder dirs contribute nothing.
    The UV_sensitive tree (different cohort/naming) is excluded.
    """
    res = subprocess.run(
        ["find", data_root, "-name", "*_L01_UDB-*_1.fq.gz", "!", "-path", "*UV_sensitive*"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    index: dict[tuple[str, str], str] = {}
    for p in res.stdout.split("\n"):
        if not p:
            continue
        try:
            if os.path.getsize(p) == 0:
                continue
        except OSError:
            continue
        # re.match (not search) anchors at the start of the basename, so macOS
        # AppleDouble junk like "._E250...UDB-438_1.fq.gz" (which starts with
        # ".") is excluded — it must not poison a readgroup with a 4 KB stub.
        m = re.match(r"(E\d+)_L01_(UDB-\d+)_1\.fq\.gz$", os.path.basename(p))
        if m:
            index[(m.group(1), m.group(2))] = p
    return index


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map", default="barcode_map.tsv", help="barcode map TSV (sample_id, fc, udb)")
    ap.add_argument("--data-root", default="/mnt/data", help="root of the mounted FASTQ storage")
    ap.add_argument("-o", "--output", default="config/samples.tsv", help="output sample sheet")
    args = ap.parse_args()

    print(f"[build_cohort_samples] indexing FASTQs under {args.data_root} ...")
    index = build_fastq_index(args.data_root)
    print(f"[build_cohort_samples] found {len(index)} (fc, udb) FASTQ pairs on disk")

    rows: list[tuple[str, str, str, str]] = []
    missing: list[tuple[str, str, str]] = []
    with open(args.map) as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or not parts[0]:
                continue
            sid, fc, udb = parts[0], parts[1], parts[2]
            r1 = index.get((fc, udb))
            if r1 is None:
                missing.append((sid, fc, udb))
                continue
            r2 = re.sub(r"_1\.fq\.gz$", "_2.fq.gz", r1)
            if not os.path.exists(r2):
                missing.append((sid, fc, udb))
                continue
            rows.append((sid, f"{fc}_{udb}", os.path.abspath(r1), os.path.abspath(r2)))

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write("sample_id\treadgroup\tfq1\tfq2\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    nsamp = len({r[0] for r in rows})
    print(f"[build_cohort_samples] wrote {len(rows)} readgroups / {nsamp} samples -> {args.output}")
    if missing:
        print(f"[build_cohort_samples] WARNING: {len(missing)} (fc, udb) from the map not found on disk:",
              file=sys.stderr)
        for m in missing[:15]:
            print(f"    {m[0]}  {m[1]}  {m[2]}", file=sys.stderr)


if __name__ == "__main__":
    main()
