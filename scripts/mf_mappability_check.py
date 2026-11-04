#!/usr/bin/env python3
"""Re-score the E-step (MosaicForecast) with an hg38 mappability track.

The yanmei/mosaicforecast image ships k24.umap.wg.bw for hg19, so the RF feature
`mappability` (column 3 of every results/filtering/{s}/mf/{s}.features) was read
at hg38 positions from an hg19 track.  The other features come from the reads, so
no CRAM is re-read: this swaps column 3 for the value MosaicForecast would have
got from an hg38 track, re-runs Prediction.R with each sample's model, and lists
what changes.  The lookup is MF's own (ReadLevel_Features_extraction.py): BED
`chrom pos pos+1 id`, score = column 6 of bigWigAverageOverBed, one call per
variant as in the E-step.  It is first run against the image's track and must
reproduce the stored column 3 exactly.

Outputs (--out-dir; the originals are not touched):
  features/{s}.features, predictions/{s}.predictions   re-scored copies
  variants.tsv  every candidate: mappability old/new, prediction old/new, in final.txt
  samples.tsv   per sample: candidates, mappability changed, mosaic old/new,
                lost from final (old mosaic in final.txt, new not mosaic), gained

Usage (VM100, pipeline directory, APPTAINER_BIND=/mnt/data,/tmp):
    python scripts/mf_mappability_check.py --results results --mf-sif containers/mosaicforecast.sif \\
        --bw /mnt/data/bsmn/resources/umap_hg38/k24.Umap.MultiTrackMappability.bw \\
        --out-dir /mnt/data/bsmn/mf_umap_hg38 --jobs 16
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

IMAGE_BW = "/usr/local/bin/k24.umap.wg.bw"
PRED_SCRIPT = "/usr/local/bin/Prediction.R"
MAP_COL = 2              # 0-based column of `mappability` in {s}.features
PRED_COL, PROB_COL = 34, 36   # same prediction columns as mosaicforecast_filter.py


def bed_line(vid: str) -> str:
    """MF's mappability BED row for a feature id 'sample~chrom~pos~ref~alt'."""
    _, chrom, pos = vid.split("~")[:3]
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return f"{chrom}\t{pos}\t{int(pos) + 1}\t{vid}"


def parse_bw_average(text: str) -> dict[str, float]:
    """bigWigAverageOverBed output -> {id: mean over covered bases} (column 6, as MF reads it)."""
    out = {}
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) >= 6:
            out[f[0]] = float(f[5])
    return out


def parse_predictions(text: str) -> dict[str, tuple[str, str]]:
    """Prediction.R output -> {id: (label, probability)}; header and short rows skipped."""
    out = {}
    for line in text.splitlines():
        f = [x.strip('"') for x in line.split("\t")]
        if len(f) > PRED_COL and "~" in f[0]:
            out[f[0]] = (f[PRED_COL], f[PROB_COL] if len(f) > PROB_COL else "")
    return out


def read_final(path: str, sample: str) -> set[str]:
    if not os.path.exists(path):
        return set()
    with open(path) as fh:
        return {f"{sample}~" + "~".join(line.split("\t")[:4]) for line in fh
                if line.strip() and not line.startswith("#")}


def swap_mappability(lines: list[str], new: dict[str, float]) -> list[str]:
    """Feature lines (header first) with column 3 replaced from `new`; unknown ids unchanged."""
    out = [lines[0]]
    for line in lines[1:]:
        f = line.rstrip("\n").split("\t")
        if f[0] in new:
            f[MAP_COL] = repr(new[f[0]])
        out.append("\t".join(f) + "\n")
    return out


def run(cmd: list[str]) -> str:
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd[:4])} ... failed: {p.stderr.strip()[-300:]}")
    return p.stdout


def lookup(args, bw: str, bed: str, tag: str) -> dict[str, float]:
    """One bigWigAverageOverBed call per variant, as the E-step ran it.

    A single call over all variants is not equivalent: for a position past the
    track's chromosome end (hg38 chr18:79.9 Mb in the hg19 track) the batch call
    reports 1 where the per-variant call reported 0 (24 of 10,576 on VM100).
    """
    with open(bed) as fh:
        lines = fh.readlines()
    loop = ('while IFS= read -r l; do printf "%s\\n" "$l" > "$1.one"; '
            'bigWigAverageOverBed "$2" "$1.one" "$1.out" && cat "$1.out"; done < "$1"; rm -f "$1.one" "$1.out"')
    parts = []
    for i in range(args.jobs):
        part = os.path.join(args.out_dir, f"umap_{tag}.part{i}")
        with open(part, "w") as fo:
            fo.writelines(lines[i::args.jobs])
        parts.append(part)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        outs = list(pool.map(lambda p: run(["apptainer", "exec", "--cleanenv", args.mf_sif,
                                            "bash", "-c", loop, "_", p, bw]), parts))
    for p in parts:
        os.remove(p)
    with open(os.path.join(args.out_dir, f"umap_{tag}.tab"), "w") as fo:
        fo.write("".join(outs))
    return parse_bw_average("".join(outs))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results")
    ap.add_argument("--mf-sif", required=True)
    ap.add_argument("--bw", required=True, help="hg38 k24 Umap bigWig")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--mode", default="Refine")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args(argv)
    for d in ("features", "predictions"):
        os.makedirs(os.path.join(args.out_dir, d), exist_ok=True)

    feats = {}
    for path in sorted(glob.glob(os.path.join(args.results, "filtering", "*", "mf", "*.features"))):
        s = os.path.basename(path)[: -len(".features")]
        with open(path) as fh:
            feats[s] = fh.readlines()
    old_map = {line.split("\t")[0]: float(line.split("\t")[MAP_COL])
               for lines in feats.values() for line in lines[1:]}
    print(f"{len(feats)} samples, {len(old_map)} candidates", flush=True)

    bed = os.path.join(args.out_dir, "candidates.bed")
    with open(bed, "w") as fo:
        fo.writelines(bed_line(v) + "\n" for v in sorted(old_map))
    image = lookup(args, IMAGE_BW, bed, "image")
    bad = [v for v in old_map if abs(image.get(v, -1) - old_map[v]) > 1e-6]
    print(f"image track reproduces stored mappability: {len(old_map) - len(bad)}/{len(old_map)}", flush=True)
    if bad:
        print(f"ERROR: lookup rule does not match MF, e.g. {bad[:3]}", file=sys.stderr)
        return 1
    new_map = lookup(args, args.bw, bed, "hg38")

    def rescore(s: str) -> tuple[str, dict, dict]:
        f_new = os.path.join(args.out_dir, "features", f"{s}.features")
        p_new = os.path.join(args.out_dir, "predictions", f"{s}.predictions")
        with open(f_new, "w") as fo:
            fo.writelines(swap_mappability(feats[s], new_map))
        with open(os.path.join(args.results, "filtering", s, f"{s}.mf_model.tsv")) as fh:
            model = fh.read().splitlines()[1].split("\t")[1]
        run(["apptainer", "exec", "--cleanenv", args.mf_sif, "Rscript", PRED_SCRIPT, f_new, model, args.mode, p_new])
        with open(os.path.join(args.results, "filtering", s, "mf", f"{s}.predictions")) as fh:
            old = parse_predictions(fh.read())
        with open(p_new) as fh:
            new = parse_predictions(fh.read())
        return s, old, new

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        scored = list(pool.map(rescore, sorted(feats)))

    unexplained = flips = total = 0
    with open(os.path.join(args.out_dir, "variants.tsv"), "w") as fv, \
         open(os.path.join(args.out_dir, "samples.tsv"), "w") as fs:
        fv.write("sample\tid\tmap_old\tmap_hg38\tpred_old\tprob_old\tpred_new\tprob_new\tin_final\n")
        fs.write("sample\tcandidates\tmap_changed\tmosaic_old\tmosaic_new\tlost_from_final\tgained\n")
        for s, old, new in scored:
            final = read_final(os.path.join(args.results, "filtering", s, f"{s}.final.txt"), s)
            n = changed = m_old = m_new = lost = gained = 0
            for v, (lab_o, pr_o) in old.items():
                lab_n, pr_n = new.get(v, ("NA", ""))
                mo, mn = lab_o.startswith("mosaic"), lab_n.startswith("mosaic")
                moved = abs(old_map.get(v, 0) - new_map.get(v, 0)) > 1e-6
                n += 1
                changed += moved
                m_old += mo
                m_new += mn
                lost += mo and not mn and v in final
                gained += mn and not mo
                unexplained += (mo != mn) and not moved
                flips += mo != mn
                fv.write(f"{s}\t{v}\t{old_map.get(v, 'NA')}\t{new_map.get(v, 'NA')}\t{lab_o}\t{pr_o}\t{lab_n}\t{pr_n}\t{int(v in final)}\n")
            fs.write(f"{s}\t{n}\t{changed}\t{m_old}\t{m_new}\t{lost}\t{gained}\n")
            total += n
    print(f"done: {total} predictions, {flips} change mosaic call; {unexplained} changed without a mappability change "
          f"(should be 0); see {args.out_dir}/samples.tsv", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
