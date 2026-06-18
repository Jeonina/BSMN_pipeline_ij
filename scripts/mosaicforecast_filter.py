#!/usr/bin/env python3
"""MosaicForecast filter — BSMN E.MosaicForecast.sh step via Apptainer.

Faithful re-implementation of ``jobs/variant_filtering/E.MosaicForecast.sh``.
Runs MosaicForecast (read-level feature extraction + trained-RF prediction) on
the surviving candidates and keeps only those predicted ``mosaic``.

The MosaicForecast scripts and the k24 mappability bigwig live inside the
``yanmei/mosaicforecast`` image (``/usr/local/bin``); the trained RF model is
provided on the host (cloned MosaicForecast repo). Reference, CRAM, candidate
BED, model and temp files are reached via Apptainer's automatic bind of the
project working directory.

Pipeline per BSMN:
  1. candidates (chr pos ref alt) -> BED  (chr  pos-1  pos  ref  alt  sample)
  2. ReadLevel_Features_extraction.py per variant (timeout + retry; MF can hang)
  3. Prediction.R  <features>  <model>  Refine  <prediction>
  4. keep rows whose prediction starts with 'mosaic' (optional prob cutoff)

The BED build and prediction parsing are pure functions (unit-tested); the
container calls cannot run without the image + alignment data.

Usage:
    python scripts/mosaicforecast_filter.py \
        --candidates results/filtering/SM/SM.cnvnator_filtered.txt \
        --sample SM --bam-dir results/mapping/SM --fmt cram \
        --ref resources/hg38/Homo_sapiens_assembly38.fasta \
        --model resources/MosaicForecast/models_trained/250xRFmodel_addRMSK_Refine.rds \
        --mf-sif containers/mosaicforecast.sif \
        --workdir results/filtering/SM/mf \
        > results/filtering/SM/SM.mosaicforecast_filtered.txt
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys

log = logging.getLogger("mosaicforecast_filter")

# In-image paths (yanmei/mosaicforecast:0.0.1).
RLF_SCRIPT = "/usr/local/bin/ReadLevel_Features_extraction.py"
PRED_SCRIPT = "/usr/local/bin/Prediction.R"
UMAP_BW = "/usr/local/bin/k24.umap.wg.bw"

# Prediction-file columns, matching the BSMN awk ($35 prediction, $37 prob).
PRED_COL = 34  # 0-based
PROB_COL = 36  # 0-based


def read_candidates(path: str) -> list[tuple[str, str, str, str]]:
    out: list[tuple[str, str, str, str]] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            out.append((p[0], p[1], p[2].upper(), p[3].upper()))
    return out


def to_bed_line(chrom: str, pos: str, ref: str, alt: str, sample: str) -> str:
    """1-based SNV -> MosaicForecast BED row (0-based start, sample in col6)."""
    return f"{chrom}\t{int(pos) - 1}\t{pos}\t{ref}\t{alt}\t{sample}"


def parse_predictions(text: str, min_prob: float = 0.0) -> list[tuple[str, str, str, str]]:
    """Return (chrom, pos, ref, alt) for rows predicted mosaic.

    Mirrors ``awk '$35~/^mosaic/ {print $1}' | cut -f2- -d~ | tr '~' '\\t'``.
    Column 1 is ``sample~chr~pos~ref~alt``; the header row is skipped naturally
    because its prediction column does not start with 'mosaic'.
    """
    out: list[tuple[str, str, str, str]] = []
    for line in text.splitlines():
        f = line.split()
        if len(f) <= PRED_COL:
            continue
        if not f[PRED_COL].startswith("mosaic"):
            continue
        if min_prob > 0.0:
            if len(f) <= PROB_COL:
                continue
            try:
                if float(f[PROB_COL]) < min_prob:
                    continue
            except ValueError:
                continue
        parts = f[0].split("~")
        if len(parts) >= 5:
            out.append((parts[1], parts[2], parts[3], parts[4]))
    return out


def _apptainer(sif: str, *cmd: str) -> list[str]:
    return ["apptainer", "exec", sif, *cmd]


def extract_features(args: argparse.Namespace, beds: list[str]) -> str | None:
    """Per-variant feature extraction with timeout + retry. Returns features path."""
    os.makedirs(args.workdir, exist_ok=True)
    features = os.path.join(args.workdir, f"{args.sample}.features")
    tin = os.path.join(args.workdir, f"{args.sample}.mf.in")
    tout = os.path.join(args.workdir, f"{args.sample}.mf.out")

    wrote_header = False
    with open(features, "w") as feat:
        for i, bed in enumerate(beds, 1):
            with open(tin, "w") as fh:
                fh.write(bed + "\n")
            log.info(">> [%d/%d] %s", i, len(beds), bed)
            # Image MF signature (no read_length): input output bam_dir ref umap nthreads fmt.
            # capture_output keeps MF's stdout/stderr OUT of our stdout (the result file).
            proc = None
            for attempt in range(args.retries):
                if os.path.exists(tout):
                    os.remove(tout)
                cmd = _apptainer(
                    args.mf_sif,
                    "python3",
                    RLF_SCRIPT,
                    tin,
                    tout,
                    args.bam_dir,
                    args.ref,
                    UMAP_BW,
                    str(args.threads),
                    args.fmt,
                )
                try:
                    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    log.warning("timeout (attempt %d/%d) — retrying", attempt + 1, args.retries)
                    continue
                break  # completed (success measured by TOUT below); only timeouts retry
            if os.path.exists(tout) and os.path.getsize(tout) > 0:
                with open(tout) as fh:
                    lines = fh.read().splitlines()
                if lines:
                    if not wrote_header:
                        feat.write(lines[0] + "\n")
                        wrote_header = True
                    feat.write(lines[-1] + "\n")  # last line = this variant's feature row
            elif proc is not None:
                log.warning(
                    "no features (rc=%s): %s",
                    proc.returncode,
                    (proc.stderr or proc.stdout or "").strip()[-400:],
                )
            else:
                log.warning("no features: all %d attempts timed out", args.retries)
    for tmp in (tin, tout):
        if os.path.exists(tmp):
            os.remove(tmp)
    if not wrote_header:
        log.warning("no features extracted for any candidate")
        return None
    return features


def predict(args: argparse.Namespace, features: str) -> str:
    prediction = os.path.join(args.workdir, f"{args.sample}.predictions")
    cmd = _apptainer(
        args.mf_sif, "Rscript", PRED_SCRIPT, features, args.model, args.mode, prediction
    )
    # capture_output so Prediction.R chatter never reaches our stdout (the result file)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log.error("Prediction.R failed (rc=%d): %s", proc.returncode, (proc.stderr or "")[-500:])
        raise SystemExit(1)
    return prediction


def run(args: argparse.Namespace) -> None:
    cands = read_candidates(args.candidates)
    log.info("candidates in=%d  sample=%s  model=%s", len(cands), args.sample, args.model)
    if not cands:
        log.info("no candidates — empty output")
        return

    beds = sorted(
        (to_bed_line(c, p, r, a, args.sample) for c, p, r, a in cands),
        key=lambda b: (b.split("\t")[0], int(b.split("\t")[2])),
    )
    features = extract_features(args, beds)
    if features is None:
        return

    prediction = predict(args, features)
    if not os.path.exists(prediction) or os.path.getsize(prediction) == 0:
        log.warning("no predictions produced — empty output")
        return

    with open(prediction) as fh:
        mosaics = parse_predictions(fh.read(), args.min_prob)
    mosaics = sorted(set(mosaics), key=lambda v: (v[0], int(v[1])))
    for chrom, pos, ref, alt in mosaics:
        sys.stdout.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")
    sys.stdout.flush()
    log.info("mosaicforecast: %d -> %d (min_prob=%g)", len(cands), len(mosaics), args.min_prob)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="MosaicForecast filter (BSMN E-step) via Apptainer.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--candidates", required=True, help="input candidates (chr pos ref alt)")
    p.add_argument("--sample", required=True, help="sample name (must match the CRAM basename)")
    p.add_argument(
        "--bam-dir", required=True, help="dir holding {sample}.{fmt} (e.g. results/mapping/SM)"
    )
    p.add_argument("--fmt", default="cram", choices=["cram", "bam"], help="alignment format [cram]")
    p.add_argument("--ref", required=True, help="reference fasta")
    p.add_argument("--model", required=True, help="trained RF model .rds")
    p.add_argument("--mf-sif", required=True, help="MosaicForecast Apptainer .sif")
    p.add_argument("--workdir", required=True, help="scratch dir for BED/features/predictions")
    p.add_argument("--mode", default="Refine", help="Prediction.R mode [Refine]")
    p.add_argument("--threads", type=int, default=4, help="threads for feature extraction [4]")
    p.add_argument("--timeout", type=int, default=300, help="per-variant timeout seconds [300]")
    p.add_argument("--retries", type=int, default=5, help="per-variant retries on timeout [5]")
    p.add_argument(
        "--min-prob",
        type=float,
        default=0.0,
        help="keep mosaic rows with prob >= this (0 = all mosaic; 0.6 = high-conf) [0.0]",
    )
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
