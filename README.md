# Somatic / Mosaic Variant Calling Pipeline (BSMN refactor)

A tumor-only somatic / mosaic SNV calling pipeline, refactored from the original
[BSMN pipeline](https://github.com/bsmn/bsmn-pipeline) for TJBaeLab.

**This is a reimplementation** — it does NOT match the original BSMN at the tool
or scheduler level:

| | Original BSMN | This pipeline |
|---|---|---|
| Orchestration | shell jobs + SGE | **Snakemake** (local or SLURM) |
| Software | conda (`bp`/`bp_frozen`) | **Apptainer containers** (pinned in `config/containers.yaml`) |
| Caller | GATK3/GATK4 HaplotypeCaller (multi-ploidy) | **GATK4 Mutect2** (tumor-only) |

> The legacy shell pipeline still lives under `jobs/` for reference, but the
> supported entry points are `run.py` and the Snakemake workflow below.

---

## Pipeline stages

```
mapping     FASTQ -> CRAM   (bwa-mem -> sambamba sort -> markdup -> BQSR)
calling     CRAM  -> VCF    (Mutect2 tumor-only, scatter per chromosome -> FilterMutectCalls)
filtering   VCF   -> final  (6-step somatic/mosaic cascade)
```

### Filtering cascade (per `workflow/rules/filtering.smk`)

```
1 accessibility   1KG strict mask (keep accessible positions)
2 germline        remove gnomAD AF > 0.001 known germline
3 vaf             binomial test p < 1e-6 AND alt_count >= 5
4 cnvnator        BSMN D-step: drop CNV-region calls (CN >= 2.5)
5 mosaicforecast  BSMN E-step: MosaicForecast trained-RF mosaic prediction
6 pon_mask        IUPAC panel-of-normals FASTA mask
                                  -> results/filtering/{sample}/{sample}.final.txt
```

The BSMN E-step is **mayo OR MosaicForecast** (two alternatives, not a chain).
This pipeline defaults to **MosaicForecast only**; the `mayo_filter` rule is
retained but unwired (re-point `mosaicforecast_filter`'s input at
`{sample}.mayo_filtered.txt` to enable it).

---

## Requirements

- **Apptainer** (or Singularity) on every compute node — all tools run in
  containers pulled from `config/containers.yaml`.
- **Python 3.9+ with `venv`** on the submit host (deps from `requirements.txt`; no conda).
- `git`, `wget`, `lftp`, `samtools` for resource download/assembly.
- ~60 GB free disk for `resources/hg38/`.
- A SLURM cluster is optional (see [docs/SLURM_USAGE.md](docs/SLURM_USAGE.md)).

---

## 1. Setup (once, on the server)

```bash
git clone <repo-url> BSMN_pipeline_ij && cd BSMN_pipeline_ij
bash scripts/server_setup.sh
```

`server_setup.sh` performs all of:
1. create the Python venv (`.venv`) and install `requirements.txt`
2. pull Apptainer containers (`scripts/prepare_containers.py`)
3. download + index public references (`scripts/download_and_index_hg38.sh`):
   reference FASTA, dbSNP, Mills, 1000G SNPs, contamination resource,
   af-only-gnomAD (Mutect2 germline resource) + its AF>0.001 SNP lookup,
   and the 1KG strict mask BED
4. assemble the PON mask FASTA and clone the MosaicForecast RF model
5. verify every required resource is present (fails loudly if any is missing)

---

## 2. Configure

Edit `config/config.yaml` (production) — key fields: `ref`, `known_sites`,
`calling.chromosomes`, `calling.germline_resource`, and the `filtering` block.
Resource and memory/thread values are auto-resolved at runtime by
`scripts/auto_params.py`.

### Sample sheet (`config/samples.tsv`)

Two schemas, auto-detected by column names:

```
# fastq mode — runs full mapping -> calling -> filtering
sample_id    readgroup    fq1                 fq2
HG002        RG1          /data/HG002.R1.fq.gz /data/HG002.R2.fq.gz

# bam mode — pre-aligned BAM/CRAM, skips mapping (see config/samples.bam.example.tsv)
sample_id    bam
AN02255      /data/AN02255.cram
```

`run.py` can generate `samples.tsv` for you from a FASTQ directory. For cohorts
laid out as **one directory per sample** (each dir holding one or more readgroup
FASTQ pairs), add `--sample-per-dir`: the parent directory becomes `sample_id`,
the files within become readgroups, and demux leftovers (`undecoded` /
`Undetermined`) are excluded automatically. Example:

```bash
python run.py /data/cohort/ --sample-per-dir --stage mapping --cores 150
```

---

## 3. Run

### Local

```bash
# one-shot from a FASTQ directory (auto-builds samples.tsv, runs all stages)
python run.py /data/fastq/ --cores 40 --stage all

# stage by stage
python run.py /data/fastq/ --stage mapping
python run.py /data/fastq/ --stage calling
python run.py /data/fastq/ --stage filtering

# plan only
python run.py /data/fastq/ --dry-run
```

Or drive Snakemake directly (e.g. for bam-mode sample sheets):

```bash
snakemake --snakefile workflow/Snakefile --config stage=filtering --cores 16
```

### SLURM cluster

```bash
export SBATCH_PARTITION=cpu-long SBATCH_ACCOUNT=mygroup
python run.py /data/cohort/ --cluster slurm --stage all
```

See [docs/SLURM_USAGE.md](docs/SLURM_USAGE.md) for partition/account overrides,
monitoring, and tuning, and [docs/PRODUCTION_RUN.md](docs/PRODUCTION_RUN.md) for
the full cohort runbook (sample sheet, MosaicForecast worker tuning, pilot-first,
outputs).

---

## 4. Output

Final somatic/mosaic SNV list per sample (4-column `chrom pos ref alt`):

```
results/filtering/{sample}/{sample}.final.txt
```

`python scripts/run_summary.py --sample {sample}` prints a per-step cascade
summary (input -> kept counts) from the run logs.

---

## Validation / examples

- `config/config.validation.yaml` — chr22 validation run.
- `config/config.chr20.yaml` + `config/samples.bam.example.tsv` — end-to-end
  chr20 example from a pre-aligned BAM/CRAM (HG002).
- [docs/PRODUCTION_VALIDATION.md](docs/PRODUCTION_VALIDATION.md) — pre-production
  runtime checks (external-tool output contracts that must be confirmed on real
  data before the first full run).
- [PIPELINE.md](PIPELINE.md) — architecture and design notes.

---

## Contributing

The `main` branch is protected. Fork, branch with a descriptive name, open a PR,
and request a review.
