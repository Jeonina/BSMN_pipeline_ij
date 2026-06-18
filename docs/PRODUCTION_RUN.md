# Production Run Guide — Brain WGS Cohort / 본 분석 풀런 가이드

Runbook for the full cohort run (e.g. 200 brain WGS samples, ~200X, hg38) on a
SLURM cluster. Assumes setup + validation are already done:
- `scripts/server_setup.sh` resource gate passed (all resources present)
- `docs/PRODUCTION_VALIDATION.md` runtime contracts confirmed on chr20

> A 200-sample WGS cohort MUST run on SLURM — never local `--cores`.

---

## 0. Pre-flight (once)

```bash
cd ~/BSMN_pipeline/BSMN_pipeline_ij
git pull origin refactor/snakemake     # ensure latest (parallel MosaicForecast)
git log --oneline -1
df -h .                                 # disk: 200 x WGS(200X) is large — ensure TB-scale free
```

Resources were already verified by `server_setup.sh` (14 files ✓); no need to
re-download.

---

## 1. Sample sheet

Two paths depending on whether the brain data is raw **FASTQ** or pre-aligned
**BAM/CRAM**.

### (A) FASTQ — `run.py` auto-builds `config/samples.tsv`

```bash
python run.py /data/brain_fastq/ --cluster slurm --stage all --dry-run
# if sample names are mis-parsed, use recursion / a custom regex:
#   --recursive  --pattern '<regex>'
```

### (B) Pre-aligned BAM/CRAM — write a bam-schema `config/samples.tsv` by hand

```
sample_id	bam
BR001	/data/brain/BR001.cram
BR002	/data/brain/BR002.cram
...   (200 rows)
```

Always `--dry-run` first and confirm all 200 samples appear and the DAG is sane.

---

## 2. Production config (`config/config.yaml`)

- `calling.chromosomes`: already chr1–22, X, Y (whole genome) — leave as is.
- `filtering.mosaicforecast.workers`: **raise to the node core count** — this is
  the pipeline bottleneck. Example for a 16-core node:

  ```yaml
  filtering:
    mosaicforecast:
      workers: 16
  ```

  Per-variant RLF cost is ~4.5 min (inherent locus I/O in deep WGS), so
  throughput = workers / 4.5 min. Memory auto-scales at 2 GB/worker
  (16 workers → ~32 GB requested per MF job).

- SLURM partition / account (override the profile placeholders without editing
  the file):

  ```bash
  export SBATCH_PARTITION=<your-partition>
  export SBATCH_ACCOUNT=<your-account>
  ```

See `docs/SLURM_USAGE.md` for more on overrides, monitoring, and tuning.

---

## 3. Pilot first (1–2 samples) — strongly recommended

Whole-genome WGS produces far more MosaicForecast candidates than the chr20
validation, so measure real MF runtime on a couple of samples before launching
all 200.

```bash
# put only 1–2 samples in samples.tsv, then:
python run.py /data/brain_fastq_pilot/ --cluster slurm --stage all
python scripts/run_summary.py --sample <id>   # per-step counts + timing
```

⚠️ If an MF job hits the SLURM `TIME LIMIT` (rule `runtime` is 720 min / 12 h),
there are too many candidates for that window — raise `runtime` in
`workflow/rules/filtering.smk` (`mosaicforecast_filter`) and/or raise `workers`.

---

## 4. Full run (200 samples)

### (A) FASTQ

```bash
python run.py /data/brain_fastq/ --cluster slurm --stage all
```

### (B) BAM/CRAM (manual samples.tsv)

```bash
snakemake --snakefile workflow/Snakefile \
  --profile workflow/profiles/slurm \
  --config stage=filtering
```

Tune max concurrent submitted jobs via `workflow/profiles/slurm/config.yaml`
`jobs:` (default 50) to match your fair-share policy.

---

## 5. Monitoring

```bash
squeue -u $USER                                            # your job queue
tail -f .snakemake/log/$(ls -t .snakemake/log/ | head -1)  # snakemake driver log
ls -lt .snakemake/slurm_logs/ | head                       # per-job SLURM output
```

---

## 6. Outputs

- Final mosaic SNV list per sample (`chr pos ref alt`):
  `results/filtering/{sample}/{sample}.final.txt`
- Per-step cascade summary for all samples:

  ```bash
  for s in $(cut -f1 config/samples.tsv | tail -n +2); do
      python scripts/run_summary.py --sample "$s"
  done
  ```

Default filtering cascade (see `workflow/rules/filtering.smk`):
`accessibility → germline → vaf → cnvnator → mosaicforecast → pon_mask`.

---

## 7. Interrupt / resume

```bash
snakemake ... --unlock                                     # if a kill left a lock
python run.py ... --snakemake-args "--rerun-incomplete"    # resume from where it stopped
```

The SLURM profile sets `restart-times: 2`, so transient node failures retry
automatically.

---

## Key reminders

1. **MosaicForecast is the dominant cost** — set `workers` = node cores; size
   `runtime` from the pilot.
2. **Pilot before the cohort** — gauge WGS candidate counts and MF time on 1–2
   samples first.
3. **Disk** — mapping CRAM + calling VCF intermediates are large. Intermediates
   are `temp()` (cleaned as the cascade advances), but watch peak usage with 200
   samples in flight.
