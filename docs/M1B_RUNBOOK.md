# M1b Cluster Runbook — Golden Snapshot Generation

SPEC: `SPEC-BSMN-REFACTOR-001` v1.2.0
Milestone: **M1b** (Cluster-side golden generation)
Audience: pipeline maintainer (`inajeon`) and any future on-call running this
on the lab SLURM cluster.

---

## 1. Purpose

### English

M1b drives the **legacy** `jobs/*.py` pipeline against the canonical FASTQ
fixtures (`tests/data/sample_R{1,2}.fastq.gz`) on the SLURM cluster, capturing
its BAM/CRAM/VCF outputs as **golden snapshots** under
`tests/data/golden/{mapping,calling,filtering,auto_params}/`. These goldens
are the byte-equivalence baselines that M4 (mapping), M5 (calling), and M6
(filtering) characterize the new Snakemake rules against. Without M1b, those
milestones cannot exit. M1b is **not** executable on a developer laptop: it
requires SLURM, Apptainer ≥ 1.2.5, the bio-tool conda env (`bp`), and the
~30 GB indexed hg38_no_alt reference.

### 한국어

M1b는 클러스터의 SLURM 환경에서 **레거시** `jobs/*.py` 파이프라인을
`tests/data/sample_R{1,2}.fastq.gz` 픽스처에 대해 실행하여, 그
BAM/CRAM/VCF 산출물을 `tests/data/golden/{mapping,calling,filtering,
auto_params}/` 아래에 **골든 스냅샷**으로 저장합니다. 이 골든 산출물은
M4(매핑)/M5(변이호출)/M6(필터링) 마일스톤에서 새로운 Snakemake 규칙의
출력이 레거시와 바이트 단위로 동일한지 확인하는 기준이 됩니다. M1b가
완료되기 전에는 M4/M5/M6의 종료 조건을 충족시킬 수 없습니다. M1b는
랩톱에서는 실행 불가능하며, SLURM·Apptainer 1.2.5 이상·`bp` 콘다 환경·
약 30 GB의 인덱싱된 hg38_no_alt 참조 게놈을 필요로 합니다.

---

## 2. Prerequisites

Mandatory on the cluster login node before running this script:

- [ ] SLURM access (`sbatch`, `squeue`, `sacct` on PATH).
- [ ] Apptainer ≥ **1.2.5** (or SingularityCE ≥ 3.11.0) — required by
      `REQ-SNK-008`.
- [ ] Conda environment named `bp` (or `bp_frozen` for the version-pinned
      variant) containing: `bwa`, `picard`, `gatk` (GATK4), `samtools`,
      `sambamba`, `bcftools`, `bgzip`, `tabix`.
- [ ] hg38_no_alt reference downloaded and indexed via
      `scripts/download_and_index_hg38.sh`. Verify by checking
      `config.hg38_no_alt.ini` `[RESOURCES] REF=…` path exists.
- [ ] gnomAD AF VCF, BSMN PON mask FASTA, 1KG strict mask FASTA downloaded
      via `download_resources.sh`. All four paths in `[RESOURCES]` of
      `config.hg38_no_alt.ini` must resolve to non-empty files.
- [ ] Repository cloned to a cluster path with ≥ 50 GB free scratch space.
- [ ] Branch `refactor/snakemake` checked out (see §3 Step 1 for SHA
      selection guidance).

---

## 3. Step-by-Step Instructions

### Step 1 — SSH to the cluster and check out the right SHA

```bash
ssh you@cluster
cd /path/to/somatic-variant-pipeline   # your existing checkout
git fetch origin
git checkout refactor/snakemake
git pull --ff-only
```

**Strict byte-equivalence path** (recommended per SPEC §M1b constraint):
M1b should run against an *unmodified* `library/` and `jobs/` tree. Since
M2 has already consolidated `library/{config,parser,...}.py` into
`bsmn_pipeline/`, the safest course is to roll back to the M1a tip:

```bash
git log --oneline | grep -iE 'M1a|M2' | head -20
# Pick the last commit BEFORE the first M2 commit
git checkout <pre-M2-sha>
```

If that is impractical, the current HEAD will still execute because
`library/job_queue.py` was deliberately retained (deletion deferred to M7).
The legacy `jobs/run_*.py` still find every symbol they need:
`bsmn_pipeline.config`, `bsmn_pipeline.parser`, `library.job_queue.GridEngineQueue`.
Document your SHA choice in the commit message that lands the goldens.

### Step 2 — Edit the configuration block of the script

```bash
nano scripts/m1b_generate_goldens.sh
# Edit the top "1. Configuration" section:
#   QUEUE       <-- your SLURM partition (e.g. "bigmem", "all", "shared")
#   CONDA_ENV   <-- "bp" or "bp_frozen"
#   WORKDIR     <-- scratch path with >= 50 GB free
# Leave REFERENCE=hg38_no_alt, SAMPLE_NAME=TEST001, and FQ1/FQ2 at defaults.
```

### Step 3 — Dry run

```bash
bash scripts/m1b_generate_goldens.sh --dry-run
```

This executes every pre-flight check (SLURM presence, Apptainer version,
conda env existence, reference + resource paths, FASTQ presence, disk space)
but submits **no** sbatch jobs. Fix every error reported before continuing.

### Step 4 — Submit the run

```bash
bash scripts/m1b_generate_goldens.sh 2>&1 | tee m1b.log
```

Optionally run inside `screen` or `tmux` since the script blocks while
polling `squeue` between stages.

### Step 5 — Monitor SLURM activity

In another shell on the cluster:

```bash
squeue -u $USER
# or, more detail:
sacct -u $USER --starttime now-2hours --format=JobID,JobName,State,Elapsed
```

### Step 6 — Verify outputs

When the script prints `M1b golden generation COMPLETE.`:

```bash
ls -la tests/data/golden/
find tests/data/golden/ -type f | sort
```

You should see at minimum:

```
tests/data/golden/mapping/TEST001.recal.cram
tests/data/golden/mapping/TEST001.recal.cram.crai
tests/data/golden/mapping/TEST001.flagstat.txt
tests/data/golden/calling/TEST001.ploidy_2.vcf.gz
tests/data/golden/calling/TEST001.ploidy_2.vcf.gz.tbi
tests/data/golden/filtering/TEST001.final.vcf.gz
tests/data/golden/filtering/TEST001.final.vcf.gz.tbi
tests/data/golden/auto_params/resolved_params.yaml
tests/data/golden/MANIFEST.yaml
```

### Step 7 — Inspect the MANIFEST

```bash
cat tests/data/golden/MANIFEST.yaml
```

Confirm fields: `spec_version: 1.2.0`, `milestone: M1b`, non-empty
`git_sha`, `reference.fasta_sha256`, all `tools.*` lines, all
`resources.*_sha256`, `sample.r1_sha256`, `sample.r2_sha256`, `cluster.host`,
`cluster.slurm_version`.

### Step 8 — Commit and push

On the cluster (or pull to laptop then commit):

```bash
git checkout refactor/snakemake          # if you ran from a pre-M2 SHA
git add tests/data/golden/
git status                                # review additions
git commit -m "chore(M1b): populate golden snapshots (SPEC-BSMN-REFACTOR-001)

Generated on $(hostname -s) at $(date -u +%FT%TZ).
Manifest: tests/data/golden/MANIFEST.yaml.
"
git push origin refactor/snakemake
```

---

## 4. Expected Runtime

Order-of-magnitude estimates for the ~9 MB test FASTQ pair on a typical
academic SLURM partition. Real numbers depend on partition wait time.

| Stage     | Wall time (compute) | Notes                                  |
| --------- | ------------------- | -------------------------------------- |
| Mapping   | ~10–30 min          | bwa + markdup + indel realign + BQSR   |
| Calling   | ~5–15 min           | GATK HaplotypeCaller (ploidy 2)        |
| Filtering | ~5 min              | accessibility + germline + VAF + PON   |
| Total     | < 1 hour            | excludes SLURM queue wait              |

If SLURM queue wait dominates, total wall-clock may extend to several hours.

---

## 5. Expected Output Sizes

| File                                                | Approx. size |
| --------------------------------------------------- | ------------ |
| `mapping/TEST001.recal.cram`                        | 5–15 MB      |
| `mapping/TEST001.recal.cram.crai`                   | < 100 KB     |
| `mapping/TEST001.flagstat.txt`                      | < 2 KB       |
| `calling/TEST001.ploidy_2.vcf.gz`                   | 50–500 KB    |
| `calling/TEST001.ploidy_2.vcf.gz.tbi`               | < 10 KB      |
| `filtering/TEST001.final.vcf.gz`                    | 5–100 KB     |
| `filtering/TEST001.final.vcf.gz.tbi`                | < 10 KB      |
| `auto_params/resolved_params.yaml`                  | < 2 KB       |
| `MANIFEST.yaml`                                     | ~3 KB        |
| **Total**                                           | **< 50 MB**  |

If the total exceeds 50 MB, surface to the user before committing — `git lfs`
adoption is **out of scope** for this SPEC (see Exclusion 7).

---

## 6. Troubleshooting

| Symptom                                                                   | Likely cause                                | Fix                                                                                                              |
| ------------------------------------------------------------------------- | ------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `ERROR: sbatch not found`                                                 | Running on laptop, not cluster              | SSH to the cluster login node and retry.                                                                         |
| `ERROR: conda env 'bp' not found`                                         | Env not created                             | `conda env create -f environment.yml` (or `environment_frozen.yml` for `bp_frozen`).                             |
| `ERROR: Reference fasta not found at …`                                   | hg38 not downloaded                         | `bash scripts/download_and_index_hg38.sh` then `bash download_resources.sh`.                                     |
| `apptainer version 1.1.x is below required floor 1.2.5`                   | Old Apptainer module                        | `module load apptainer/1.2.5` (cluster-specific) or ask cluster admins to upgrade.                               |
| SLURM job stuck in `PD` (pending) for > 30 min                            | Wrong `--partition` or partition full       | Edit `QUEUE` in the script. Use `sinfo` to list available partitions.                                            |
| `WARNING: no file matched *.recal.cram` in Stage 5                        | Mapping stage failed silently               | Inspect `${WORKDIR}/${SAMPLE_NAME}/logs/` for sbatch stderr.                                                      |
| VCF empty or contains zero records                                        | Expected — the 9 MB test FASTQ is too small | Not a failure. Char-tests will still assert header structure + record count parity.                              |
| `INSUFFICIENT_FREE_GB` failure                                            | Scratch quota exhausted                     | Set `WORKDIR=/some/other/scratch/path` env var and re-run.                                                       |
| Script aborts mid-run                                                     | Any of the above                            | `$WORKDIR` is **preserved on failure**. Inspect `${WORKDIR}/${SAMPLE_NAME}/run_jid` + logs, fix, then re-run.    |

---

## 7. Resuming After Partial Failure

The legacy `jobs/run_*.py` use **per-sample `run_jid` files** to track SLURM
job state. If a stage fails partway:

```bash
# Inspect what jobs got submitted
cat ${WORKDIR}/${SAMPLE_NAME}/run_jid
sacct -j $(paste -sd, ${WORKDIR}/${SAMPLE_NAME}/run_jid) \
      --format=JobID,JobName,State,ExitCode
```

To reset and re-run from scratch (safest):

```bash
rm -rf ${WORKDIR}/${SAMPLE_NAME}
bash scripts/m1b_generate_goldens.sh    # picks up clean state
```

To re-run only the failed stage, manually invoke that single `jobs/run_*.py`
with the same arguments shown in the script's Stage 2 / Stage 3 blocks, then
re-execute the script with `--force` to refresh `tests/data/golden/`.

If `tests/data/golden/` is already populated from a previous successful run,
the script exits 0 with a notice. Pass `--force` to overwrite.

---

## 8. What Happens Next

Once `tests/data/golden/` is committed and pushed to `refactor/snakemake`:

1. **On any machine** (cluster or laptop):

   ```bash
   git checkout refactor/snakemake
   git pull --ff-only
   ls tests/data/golden/MANIFEST.yaml   # should exist
   ```

2. **In a fresh Claude Code session** at the project root:

   ```
   /moai run SPEC-BSMN-REFACTOR-001 M4
   ```

   This unblocks the mapping characterization tests. M5 and M6 follow once
   M4 lands.

3. **Manifest auditability:** any future contributor can trace any golden
   back to the cluster invocation that produced it via
   `tests/data/golden/MANIFEST.yaml` — git SHA, host, SLURM version, tool
   versions, and resource SHA-256s are all recorded.

---

## 9. References

- SPEC: `.moai/specs/SPEC-BSMN-REFACTOR-001/spec.md` §5 (M1b), §3.5
  (REQ-TST-002/003), §4 (AC-22), §6 (R-NEW-1, R-NEW-2).
- Script: `scripts/m1b_generate_goldens.sh`.
- Legacy entry points (do **not** edit): `jobs/run_genome_mapping.py`,
  `jobs/run_variant_calling.py`, `jobs/run_variant_filtering.py`.
- SLURM helper: `library/job_queue.py` (held until M7 by design).
- Reference config: `config.hg38_no_alt.ini`.
