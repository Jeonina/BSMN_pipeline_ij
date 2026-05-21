# SLURM Cluster Usage Guide / SLURM 클러스터 사용 가이드

This guide explains how to run the BSMN pipeline on a multi-node SLURM cluster.
다중 노드 SLURM 클러스터에서 BSMN 파이프라인을 실행하는 방법입니다.

## 1. Prerequisites / 사전 준비

- A SLURM-managed cluster with `sbatch` and `squeue` on `PATH`.
- Apptainer (or Singularity) on every compute node — containers in
  `workflow/envs/` provide the bioconda toolchain.
- The conda `bp` environment (or compatible Python 3.13+) on the submit host.
- Reference data downloaded once (see `download_resources.sh`).

한국어: `sbatch`/`squeue`가 PATH에 있는 SLURM 클러스터, 모든 계산 노드에 Apptainer
(또는 Singularity), 제출 호스트에 conda `bp` 환경, 한 번의 레퍼런스 데이터 다운로드.

## 2. First-time setup / 최초 설정

```bash
git clone <repo-url> bsmn_pipeline && cd bsmn_pipeline
conda env create -f environment.yml
conda activate bp
bash download_resources.sh           # downloads hg38 references + indexes
```

## 3. Per-cluster customization / 클러스터별 커스터마이즈

The shipped profile uses placeholder partition `default` and an empty account.
Override these for your site without editing the profile file:

배포된 프로파일은 placeholder 값(`default` 파티션, 빈 account)을 사용합니다.
프로파일 파일을 수정하지 말고 다음 중 한 가지 방법으로 override 하세요.

```bash
# Option A — environment variables (recommended; the slurm plugin reads them)
export SBATCH_PARTITION=cpu-long
export SBATCH_ACCOUNT=mygroup

# Option B — pass through to snakemake via --snakemake-args
python run.py /data/cohort/ --cluster slurm \
    --snakemake-args "--default-resources slurm_partition=cpu-long slurm_account=mygroup"
```

## 4. Quick start / 빠른 실행

```bash
# Run the full pipeline on a cohort directory via SLURM
python run.py /data/cohort/ --cluster slurm --stage all

# Stage-by-stage
python run.py /data/cohort/ --cluster slurm --stage mapping
python run.py /data/cohort/ --cluster slurm --stage calling
python run.py /data/cohort/ --cluster slurm --stage filtering

# Dry-run (plan only, no job submission)
python run.py /data/cohort/ --cluster slurm --dry-run
```

## 5. Monitoring / 모니터링

```bash
# Job queue status (only your jobs)
squeue -u $USER

# Live snakemake driver log
tail -f .snakemake/log/$(ls -t .snakemake/log/ | head -1)

# Per-job SLURM output (the plugin writes them under .snakemake/slurm_logs/)
ls -lt .snakemake/slurm_logs/ | head -20
```

## 6. Common errors / 흔한 에러

| Symptom | Cause | Fix |
|---|---|---|
| `sbatch: error: Invalid partition specified: default` | Profile placeholder not overridden | Set `SBATCH_PARTITION` (see Section 3) |
| `slurmstepd: error: TIME LIMIT` | `runtime` resource too low for that rule | Edit the rule's `resources.runtime` or pass `--default-resources runtime=...` |
| `slurmstepd: error: Out of memory` | `mem_mb` too low | Increase the rule's `resources.mem_mb` |
| `FATAL: container creation failed` | Apptainer not on compute node, or container path wrong | Confirm Apptainer module is loaded in the SLURM prolog |
| `sbatch가 PATH에서 검출되지 않았습니다` | `run.py` failed fast — not on a SLURM submit host | Run on the submit node, not your laptop |

## 7. Performance tuning / 성능 튜닝

- **Concurrent jobs**: edit `workflow/profiles/slurm/config.yaml` `jobs:` (default
  50) to respect your fair-share policy.
- **Per-rule overrides**: bump `mem_mb` or `runtime` in `workflow/rules/*.smk`
  for outlier rules (e.g., variant calling on deep WGS samples).
- **Restart-on-failure**: the profile sets `restart-times: 2`. Increase if your
  cluster has flaky nodes; decrease to fail fast.
- **Resume partial runs**: `python run.py ... --cluster slurm \
  --snakemake-args "--rerun-incomplete"` continues after manual cancellation.

한국어: 동시 작업 수는 `jobs:`로 조절. 특정 규칙의 자원은 `workflow/rules/*.smk`의
`resources:` 블록을 수정. `--rerun-incomplete`로 중단된 실행 재개 가능.

## 8. References / 참고

- Snakemake SLURM executor plugin:
  https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html
- Snakemake profiles guide:
  https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
- SPEC: `.moai/specs/SPEC-BSMN-REFACTOR-001/spec.md` (Phase 3 — SLURM cluster profile)
