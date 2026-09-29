#!/usr/bin/env python3
"""Per-sample MosaicForecast model selection by measured sequencing depth.

The MosaicForecast RF models are depth-specific (50x/100x/150x/... in
models_trained/). Scoring a sample with a model trained at a different depth
silently changes every mosaic call, and a cohort rarely has one depth: CRAM
sizes in a single cohort easily span ~2x. So the model is chosen PER SAMPLE:

  1. measure mean depth over one region with `samtools coverage` (default
     chr20 - reads only that region through the CRAM index, and samtools'
     default flag filter already drops unmapped/secondary/QC-fail/duplicate
     reads);
  2. pick the model whose depth label is nearest (a tie goes to the lower
     label, so the choice is deterministic);
  3. write a one-row TSV recording the depth and the model, so every sample's
     calls can be traced back to the model that produced them.

With --model PATH the given model is used as-is (a deliberate pin), but the
depth is still measured and a warning is logged when a different model would
have been nearer - the mismatch is otherwise invisible.

Usage:
    python scripts/select_mf_model.py \
        --cram results/mapping/SM/SM.cram --ref ref.fasta \
        --samtools-sif containers/samtools.sif \
        --models-dir resources/MosaicForecast/models_trained --mode Refine \
        --model auto --region chr20 \
        > results/filtering/SM/SM.mf_model.tsv
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import subprocess
import sys

log = logging.getLogger("select_mf_model")

HEADER = "meandepth\tmodel\tmodel_depth\tsource"


def model_pattern(mode: str) -> re.Pattern[str]:
    """Trained-model file names carry the depth: 100xRFmodel_addRMSK_Refine.rds."""
    return re.compile(r"^(\d+)x[A-Za-z_]*_" + re.escape(mode) + r"\.rds$")


def list_models(models_dir: str, mode: str) -> dict[int, str]:
    """{depth_label: path} for every model in models_dir matching the mode."""
    pat = model_pattern(mode)
    out: dict[int, str] = {}
    for name in sorted(os.listdir(models_dir)):
        m = pat.match(name)
        if m:
            out[int(m.group(1))] = os.path.join(models_dir, name)
    return dict(sorted(out.items()))


def nearest_depth(depth: float, labels: list[int]) -> int:
    """Nearest depth label; a tie goes to the lower label."""
    if not labels:
        raise ValueError("no model depth labels to choose from")
    return min(sorted(labels), key=lambda d: abs(d - depth))


def parse_meandepth(text: str) -> float:
    """meandepth from `samtools coverage` output (single-region run)."""
    header: list[str] | None = None
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if line.startswith("#"):
            header = [f.lstrip("#") for f in fields]
            continue
        if header is None:
            raise ValueError("samtools coverage output has no header line")
        return float(fields[header.index("meandepth")])
    raise ValueError("samtools coverage output has no data row")


def measure_depth(cram: str, ref: str, region: str, samtools_sif: str) -> float:
    cmd = [
        "apptainer",
        "exec",
        samtools_sif,
        "samtools",
        "coverage",
        "--reference",
        ref,
        "-r",
        region,
        cram,
    ]
    log.info("measuring depth: %s", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        log.error("samtools coverage failed (exit %d):\n%s", res.returncode, res.stderr.strip())
        if "No such file" in res.stderr and os.path.exists(cram):
            # Visible on the host but not in the container: a path outside the
            # working directory (e.g. results/ symlinked to NFS) that is not bound.
            log.error(
                "%s exists on the host but not inside the container - "
                "export APPTAINER_BIND to include its real location (%s)",
                cram,
                os.path.realpath(cram),
            )
        raise SystemExit(1)
    return parse_meandepth(res.stdout)


def choose(depth: float, models: dict[int, str], pinned: str) -> tuple[str, int | None, str]:
    """(model_path, model_depth_label, source) for the sample.

    pinned == "auto" selects by depth; any other value is a user-pinned path.
    """
    if pinned != "auto":
        label_m = re.match(r"^(\d+)x", os.path.basename(pinned))
        label = int(label_m.group(1)) if label_m else None
        if models:
            best = nearest_depth(depth, list(models))
            if label is not None and best != label:
                log.warning(
                    "pinned model %s (%dx) but measured depth %.1fx is nearest %dx (%s)",
                    pinned,
                    label,
                    depth,
                    best,
                    models[best],
                )
        return pinned, label, "pinned"
    if not models:
        raise ValueError("no trained models to choose from")
    best = nearest_depth(depth, list(models))
    return models[best], best, "auto"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cram", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--samtools-sif", required=True)
    ap.add_argument("--models-dir", required=True)
    ap.add_argument("--mode", default="Refine")
    ap.add_argument("--model", default="auto", help="'auto' or a model path to pin")
    ap.add_argument("--region", default="chr20")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    models = list_models(args.models_dir, args.mode)
    log.info("models (%s): %s", args.mode, ", ".join(f"{d}x" for d in models) or "none")
    depth = measure_depth(args.cram, args.ref, args.region, args.samtools_sif)
    model, label, source = choose(depth, models, args.model)
    log.info("meandepth=%.2fx over %s -> %s (%s)", depth, args.region, model, source)
    print(HEADER)
    print(f"{depth:.2f}\t{model}\t{label if label is not None else 'NA'}\t{source}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
