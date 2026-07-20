#!/usr/bin/env bash
# =============================================================================
# M1b - Generate golden snapshots for SPEC-BSMN-REFACTOR-001 v1.2.0
#
# Run this on the SLURM cluster, NOT on a laptop.
#
# Purpose
# -------
# Drive the *legacy* jobs/*.py pipeline against the canonical test fixtures at
# tests/data/sample_R{1,2}.fastq.gz to produce byte-equivalent golden BAM/VCF
# outputs under tests/data/golden/. These goldens become the characterization
# baselines for the new Snakemake rules (M4/M5/M6).
#
# This script assumes:
#   - You are on a Linux SLURM cluster with sbatch, squeue, apptainer.
#   - You have a conda env (default: bp) with the BSMN bio-tool chain.
#   - You have downloaded hg38_no_alt reference and resource files (see
#     scripts/download_and_index_hg38.sh and download_resources.sh).
#   - You have checked out branch refactor/snakemake (or a pre-M2 commit
#     for strict byte-equivalence; see docs/M1B_RUNBOOK.md).
#
# Idempotency
# -----------
# Re-running this script after tests/data/golden/ is populated exits 0 with a
# message, unless --force is supplied.
#
# Flags
# -----
#   --dry-run   Run all pre-flight checks but do not submit any sbatch jobs.
#   --force     Overwrite any existing tests/data/golden/ contents.
#   --help      Print this header and exit.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# 1. Configuration - REVIEW AND EDIT BEFORE RUNNING
# -----------------------------------------------------------------------------

# Absolute path to the cloned repo on the cluster (must contain pyproject.toml
# and jobs/run_genome_mapping.py at its root).
PIPELINE_HOME="${PIPELINE_HOME:-$(pwd)}"

# SLURM partition / queue name. CHANGE THIS to match your cluster.
QUEUE="${QUEUE:-bigmem}"

# Conda env carrying bwa, picard, GATK4, samtools, sambamba, bcftools.
# Default 'bp' per existing README. Use 'bp_frozen' for the version-pinned env.
CONDA_ENV="${CONDA_ENV:-bp}"

# Reference selector passed to jobs/*.py. The hg38 config file used is
# config.hg38_no_alt.ini regardless (REQ-SNK-006). The -r flag value MUST be
# the basename of the .ini file without the 'config.' prefix and '.ini' suffix.
REFERENCE="${REFERENCE:-hg38_no_alt}"

# Working directory where intermediate sample folders, run_jid files, and
# per-sample outputs accumulate before being copied to tests/data/golden/.
# Use a scratch path with at least 50 GB free space.
WORKDIR="${WORKDIR:-${PIPELINE_HOME}/_m1b_work}"

# Logical sample identifier. Used as the per-sample subdirectory name. Keep
# stable across stages.
SAMPLE_NAME="${SAMPLE_NAME:-TEST001}"

# Input FASTQ paths. Default to the canonical fixtures in tests/data/.
FQ1="${FQ1:-${PIPELINE_HOME}/tests/data/sample_R1.fastq.gz}"
FQ2="${FQ2:-${PIPELINE_HOME}/tests/data/sample_R2.fastq.gz}"

# Polling interval (seconds) while waiting for SLURM jobs to drain.
POLL_INTERVAL="${POLL_INTERVAL:-60}"

# Minimum required free disk space (GB) under WORKDIR.
MIN_FREE_GB="${MIN_FREE_GB:-50}"

# Apptainer / SingularityCE version floor (REQ-SNK-008).
APPTAINER_FLOOR="1.2.5"

# Git SHA at start - captured for MANIFEST.yaml provenance.
GIT_SHA_AT_START=""

# Flags (set by argument parser below).
DRY_RUN=0
FORCE=0

# Final golden output directory inside the repo.
GOLDEN_DIR="${PIPELINE_HOME}/tests/data/golden"

# -----------------------------------------------------------------------------
# 2. Argument parsing
# -----------------------------------------------------------------------------

usage() {
    sed -n '2,40p' "$0"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=1; shift ;;
        --force)   FORCE=1;   shift ;;
        --help|-h) usage ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

# -----------------------------------------------------------------------------
# 3. Helpers
# -----------------------------------------------------------------------------

log()  { printf '[m1b] %s\n' "$*" >&2; }
die()  { printf '[m1b] ERROR: %s\n' "$*" >&2; exit 1; }
note() { printf '[m1b]        %s\n' "$*" >&2; }

# Version comparator: returns 0 (true) if $1 >= $2 (dot-separated semver-ish).
version_ge() {
    [ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n1)" = "$2" ]
}

# Strip the leading 'config.' and trailing '.ini' from a config filename to
# get the -r/--reference selector value the legacy jobs/*.py expects.
config_ini_path() {
    printf '%s/config.%s.ini' "$PIPELINE_HOME" "$REFERENCE"
}

# Parse a key from the [RESOURCES] section of an .ini file. Substitutes the
# {PIPEHOME} placeholder with $PIPELINE_HOME.
ini_get() {
    local ini="$1" key="$2"
    awk -F'=' -v k="$key" '
        /^\[/ { section=$0; next }
        section ~ /\[RESOURCES\]/ && $1 ~ "^[[:space:]]*"k"[[:space:]]*$" {
            sub(/^[[:space:]]+/, "", $2); sub(/[[:space:]]+$/, "", $2); print $2; exit
        }
    ' "$ini" | sed "s|{PIPEHOME}|${PIPELINE_HOME}|g"
}

sha256_of() {
    # Cross-platform sha256.
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | cut -d' ' -f1
    else
        shasum -a 256 "$1" | cut -d' ' -f1
    fi
}

# -----------------------------------------------------------------------------
# 4. Idempotency check
# -----------------------------------------------------------------------------

if [[ -d "$GOLDEN_DIR" ]] && find "$GOLDEN_DIR" -mindepth 1 -type f 2>/dev/null | read -r _; then
    if [[ $FORCE -eq 0 ]]; then
        log "tests/data/golden/ already populated. Re-run with --force to overwrite."
        log "Existing layout:"
        ls -la "$GOLDEN_DIR" >&2 || true
        exit 0
    else
        log "--force supplied; will overwrite $GOLDEN_DIR after pre-flight checks pass."
    fi
fi

# -----------------------------------------------------------------------------
# 5. Pre-flight checks
# -----------------------------------------------------------------------------

log "===== M1b pre-flight checks ====="

# 5.1 cwd / repo identity
[[ -f "${PIPELINE_HOME}/pyproject.toml" ]] || die "pyproject.toml not found at $PIPELINE_HOME"
[[ -f "${PIPELINE_HOME}/jobs/run_genome_mapping.py" ]]    || die "jobs/run_genome_mapping.py missing - not a BSMN pipeline checkout"
[[ -f "${PIPELINE_HOME}/jobs/run_variant_calling.py" ]]   || die "jobs/run_variant_calling.py missing"
[[ -f "${PIPELINE_HOME}/jobs/run_variant_filtering.py" ]] || die "jobs/run_variant_filtering.py missing"
[[ -f "${PIPELINE_HOME}/library/job_queue.py" ]]          || die "library/job_queue.py missing - cluster runbook requires the legacy GridEngineQueue"
note "repo root verified at $PIPELINE_HOME"

# 5.2 git SHA capture
if command -v git >/dev/null 2>&1 && git -C "$PIPELINE_HOME" rev-parse HEAD >/dev/null 2>&1; then
    GIT_SHA_AT_START="$(git -C "$PIPELINE_HOME" rev-parse HEAD)"
    note "git SHA: $GIT_SHA_AT_START"
else
    GIT_SHA_AT_START="unknown"
    note "WARNING: git SHA could not be captured (not a git checkout?)"
fi

# 5.3 SLURM
command -v sbatch >/dev/null 2>&1 || die "sbatch not found - this script must run on a SLURM cluster"
command -v squeue >/dev/null 2>&1 || die "squeue not found"
note "SLURM tools present"

# 5.4 Apptainer / SingularityCE
APPTAINER_BIN=""
if command -v apptainer >/dev/null 2>&1; then
    APPTAINER_BIN="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    APPTAINER_BIN="singularity"
else
    die "neither apptainer nor singularity found in PATH (REQ-SNK-008 requires apptainer >= ${APPTAINER_FLOOR})"
fi
APPTAINER_VER="$("$APPTAINER_BIN" --version 2>&1 | awk '{print $NF}' | head -n1)"
version_ge "$APPTAINER_VER" "$APPTAINER_FLOOR" \
    || die "$APPTAINER_BIN version $APPTAINER_VER is below required floor $APPTAINER_FLOOR (REQ-SNK-008)"
note "$APPTAINER_BIN $APPTAINER_VER OK (>= $APPTAINER_FLOOR)"

# 5.5 Conda env
command -v conda >/dev/null 2>&1 || die "conda not found in PATH - source your miniconda/anaconda profile first"
conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV" \
    || die "conda env '$CONDA_ENV' not found. Create via: conda env create -f environment.yml (or environment_frozen.yml)"
note "conda env '$CONDA_ENV' exists"

# 5.6 hg38_no_alt reference and resources
CONFIG_INI="$(config_ini_path)"
[[ -f "$CONFIG_INI" ]] || die "config file not found at $CONFIG_INI (REFERENCE='$REFERENCE' is wrong?)"
note "config: $CONFIG_INI"

REFERENCE_FASTA="$(ini_get "$CONFIG_INI" REF)"
[[ -n "$REFERENCE_FASTA" ]] || die "REF key not parsed from $CONFIG_INI"
[[ -f "$REFERENCE_FASTA" ]] || die "Reference fasta not found at $REFERENCE_FASTA - run scripts/download_and_index_hg38.sh"
[[ -f "${REFERENCE_FASTA}.fai" ]] || die "Reference fai index missing at ${REFERENCE_FASTA}.fai"
note "reference fasta: $REFERENCE_FASTA"

PONFA="$(ini_get "$CONFIG_INI" PONFA)"
[[ -f "$PONFA" ]] || die "PON mask fasta missing at $PONFA - run download_resources.sh"
note "PON mask: $PONFA"

GNOMAD="$(ini_get "$CONFIG_INI" GNOMAD_SNP)"
[[ -f "$GNOMAD" ]] || die "gnomAD SNP file missing at $GNOMAD - run download_resources.sh"
note "gnomAD: $GNOMAD"

MASK1KG="$(ini_get "$CONFIG_INI" MASK1KG)"
[[ -f "$MASK1KG" ]] || die "1KG strict mask missing at $MASK1KG - run download_resources.sh"
note "1KG mask: $MASK1KG"

# 5.7 Input FASTQs
[[ -s "$FQ1" ]] || die "FQ1 not found or empty: $FQ1"
[[ -s "$FQ2" ]] || die "FQ2 not found or empty: $FQ2"
note "input: $FQ1"
note "input: $FQ2"

# 5.8 Output writability
mkdir -p "$GOLDEN_DIR"
[[ -w "$GOLDEN_DIR" ]] || die "$GOLDEN_DIR is not writable"
note "$GOLDEN_DIR is writable"

# 5.9 WORKDIR + disk space
mkdir -p "$WORKDIR"
[[ -w "$WORKDIR" ]] || die "$WORKDIR is not writable"
FREE_GB="$(df -BG --output=avail "$WORKDIR" 2>/dev/null | tail -n1 | tr -dc '0-9' || echo 0)"
if [[ -z "$FREE_GB" || "$FREE_GB" -eq 0 ]]; then
    # macOS / BSD df fallback (cluster is Linux but be defensive).
    FREE_GB="$(df -g "$WORKDIR" 2>/dev/null | awk 'NR==2 {print $4}' || echo 0)"
fi
[[ "${FREE_GB:-0}" -ge "$MIN_FREE_GB" ]] \
    || die "Insufficient free space at $WORKDIR: ${FREE_GB} GB available, need >= ${MIN_FREE_GB} GB"
note "free space at $WORKDIR: ${FREE_GB} GB"

log "Pre-flight checks PASSED."

if [[ $DRY_RUN -eq 1 ]]; then
    log "--dry-run supplied; exiting before any sbatch submission."
    exit 0
fi

# -----------------------------------------------------------------------------
# 6. Failure trap (preserve $WORKDIR on failure, remove on success)
# -----------------------------------------------------------------------------

CLEAN_ON_EXIT=0
on_exit() {
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        log "FAILED (exit $rc). Preserving $WORKDIR for debugging."
    elif [[ $CLEAN_ON_EXIT -eq 1 ]]; then
        log "Success - $WORKDIR preserved (manual cleanup recommended once goldens are committed)."
    fi
}
trap on_exit EXIT

# -----------------------------------------------------------------------------
# 7. sample_list.txt generation (legacy format from library/parser.py)
# -----------------------------------------------------------------------------

SAMPLE_LIST="${WORKDIR}/sample_list.txt"
log "===== Writing sample_list.txt ====="
{
    printf '#sample_id\tfile_name\tlocation\n'
    printf '%s\t%s\t%s\n' "$SAMPLE_NAME" "$(basename "$FQ1")" "$FQ1"
    printf '%s\t%s\t%s\n' "$SAMPLE_NAME" "$(basename "$FQ2")" "$FQ2"
} > "$SAMPLE_LIST"
note "wrote $SAMPLE_LIST"
cat "$SAMPLE_LIST" >&2

# -----------------------------------------------------------------------------
# 8. SLURM job-wait helper
# -----------------------------------------------------------------------------

# Collect all job IDs from a sample's run_jid file and squeue-wait until empty.
wait_for_sample_jobs() {
    local sample_dir="$1"
    local run_jid_file="${sample_dir}/run_jid"
    if [[ ! -f "$run_jid_file" ]]; then
        log "WARNING: no run_jid file at $run_jid_file - nothing to wait for"
        return 0
    fi
    local jids
    jids="$(paste -sd, "$run_jid_file" | tr -d '[:space:]')"
    if [[ -z "$jids" ]]; then
        log "WARNING: run_jid file is empty"
        return 0
    fi
    log "Waiting on SLURM jobs: $jids"
    local elapsed=0
    while :; do
        local remaining
        remaining="$(squeue -h --jobs "$jids" 2>/dev/null | wc -l | tr -d '[:space:]')"
        if [[ "${remaining:-0}" -eq 0 ]]; then
            log "All jobs drained after ${elapsed}s"
            break
        fi
        log "  ${remaining} job(s) still running/pending (elapsed ${elapsed}s)"
        sleep "$POLL_INTERVAL"
        elapsed=$((elapsed + POLL_INTERVAL))
    done

    # Post-mortem: any FAILED / CANCELLED / TIMEOUT?
    local bad
    bad="$(sacct -j "$jids" --format=JobID,State --noheader --parsable2 2>/dev/null \
            | awk -F'|' '$2 !~ /COMPLETED|RUNNING|PENDING/ && $2 != "" {print $1"="$2}' || true)"
    if [[ -n "$bad" ]]; then
        die "One or more jobs did not complete cleanly: $bad"
    fi
}

# -----------------------------------------------------------------------------
# 9. Stage 1: genome mapping
# -----------------------------------------------------------------------------

log "===== Stage 1: genome mapping ====="
(
    cd "$WORKDIR"
    python3 "${PIPELINE_HOME}/jobs/run_genome_mapping.py" \
        --queue       "$QUEUE" \
        --conda-env   "$CONDA_ENV" \
        --reference   "$REFERENCE" \
        --align-fmt   "cram" \
        --sample-list "$SAMPLE_LIST"
)
wait_for_sample_jobs "${WORKDIR}/${SAMPLE_NAME}"

# -----------------------------------------------------------------------------
# 10. Stage 2: variant calling
# -----------------------------------------------------------------------------

log "===== Stage 2: variant calling ====="
(
    cd "$WORKDIR"
    python3 "${PIPELINE_HOME}/jobs/run_variant_calling.py" \
        --queue       "$QUEUE" \
        --conda-env   "$CONDA_ENV" \
        --reference   "$REFERENCE" \
        --align-fmt   "cram" \
        --run-gatk-hc 2 \
        --sample-list "$SAMPLE_LIST"
)
wait_for_sample_jobs "${WORKDIR}/${SAMPLE_NAME}"

# -----------------------------------------------------------------------------
# 11. Stage 3: variant filtering
# -----------------------------------------------------------------------------

log "===== Stage 3: variant filtering ====="
(
    cd "$WORKDIR"
    python3 "${PIPELINE_HOME}/jobs/run_variant_filtering.py" \
        --queue       "$QUEUE" \
        --conda-env   "$CONDA_ENV" \
        --reference   "$REFERENCE" \
        --align-fmt   "cram" \
        --run-gatk-hc 2 \
        --sample-list "$SAMPLE_LIST"
)
wait_for_sample_jobs "${WORKDIR}/${SAMPLE_NAME}"

# -----------------------------------------------------------------------------
# 12. Stage 4: VCF header normalization (R-NEW-2 mitigation)
#     Strip absolute cluster paths from VCF headers so committed goldens are
#     environment-portable.
# -----------------------------------------------------------------------------

log "===== Stage 4: VCF header normalization ====="

# Activate conda env so bcftools/bgzip/tabix are on PATH.
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

shopt -s globstar nullglob
NORMALIZED_COUNT=0
for vcf in "$WORKDIR"/**/*.vcf.gz; do
    [[ -f "$vcf" ]] || continue
    log "  normalizing $vcf"
    tmpfile="${vcf%.vcf.gz}.norm.vcf.gz"
    bcftools view "$vcf" \
        | sed -e "s|${PIPELINE_HOME}|/PIPELINE|g" \
              -e "s|${WORKDIR}|/PIPELINE/_work|g" \
              -e "s|/scratch/[^[:space:]\"',;]*|/PIPELINE/scratch|g" \
              -e "s|/home/[^[:space:]\"',;/]*|/PIPELINE/home|g" \
        | bgzip > "$tmpfile"
    mv "$tmpfile" "$vcf"
    tabix -f -p vcf "$vcf"
    NORMALIZED_COUNT=$((NORMALIZED_COUNT + 1))
done
shopt -u globstar nullglob
log "Normalized $NORMALIZED_COUNT VCF file(s)."

# -----------------------------------------------------------------------------
# 13. Stage 5: copy results into tests/data/golden/
# -----------------------------------------------------------------------------

log "===== Stage 5: assembling tests/data/golden/ ====="

if [[ $FORCE -eq 1 ]]; then
    rm -rf "${GOLDEN_DIR}/mapping" "${GOLDEN_DIR}/calling" "${GOLDEN_DIR}/filtering" "${GOLDEN_DIR}/auto_params"
fi

mkdir -p \
    "${GOLDEN_DIR}/mapping" \
    "${GOLDEN_DIR}/calling" \
    "${GOLDEN_DIR}/filtering" \
    "${GOLDEN_DIR}/auto_params"

SAMPLE_OUTDIR="${WORKDIR}/${SAMPLE_NAME}"

# Locate produced artifacts (paths follow the legacy jobs/genome_mapping/ and
# jobs/variant_filtering/ conventions). The find-first-then-copy approach is
# defensive: if the small fixture produces unexpected layouts we capture
# everything matching the artifact type.

copy_first() {
    local glob="$1" dest="$2"
    local found
    # shellcheck disable=SC2086
    found="$(find "$SAMPLE_OUTDIR" -type f -name "$glob" 2>/dev/null | head -n1)"
    if [[ -n "$found" ]]; then
        cp -v "$found" "$dest" >&2
        return 0
    fi
    log "  WARNING: no file matched $glob under $SAMPLE_OUTDIR"
    return 1
}

# Mapping artifacts
copy_first "*.recal.cram"      "${GOLDEN_DIR}/mapping/${SAMPLE_NAME}.recal.cram"      || true
copy_first "*.recal.cram.crai" "${GOLDEN_DIR}/mapping/${SAMPLE_NAME}.recal.cram.crai" || true
copy_first "*.recal.bam"       "${GOLDEN_DIR}/mapping/${SAMPLE_NAME}.recal.bam"      || true
copy_first "*.recal.bam.bai"   "${GOLDEN_DIR}/mapping/${SAMPLE_NAME}.recal.bam.bai"  || true
copy_first "*flagstat*"        "${GOLDEN_DIR}/mapping/${SAMPLE_NAME}.flagstat.txt"   || true

# Calling artifacts (Mutect2 / HaplotypeCaller ploidy_2)
copy_first "*ploidy_2*.vcf.gz"     "${GOLDEN_DIR}/calling/${SAMPLE_NAME}.ploidy_2.vcf.gz"     || true
copy_first "*ploidy_2*.vcf.gz.tbi" "${GOLDEN_DIR}/calling/${SAMPLE_NAME}.ploidy_2.vcf.gz.tbi" || true

# Filtering artifacts (final filtered VCF + per-filter outputs)
copy_first "*final*.vcf.gz"     "${GOLDEN_DIR}/filtering/${SAMPLE_NAME}.final.vcf.gz"     || true
copy_first "*final*.vcf.gz.tbi" "${GOLDEN_DIR}/filtering/${SAMPLE_NAME}.final.vcf.gz.tbi" || true
# Per-filter snapshots (accessibility, germline, vaf, pon_mask). Best-effort.
for sub in accessibility germline vaf pon_mask; do
    mkdir -p "${GOLDEN_DIR}/filtering/${sub}"
    while IFS= read -r f; do
        [[ -n "$f" ]] && cp -v "$f" "${GOLDEN_DIR}/filtering/${sub}/" >&2
    done < <(find "$SAMPLE_OUTDIR" -type f -iname "*${sub}*" 2>/dev/null || true)
done

# auto_params snapshot
log "===== Generating resolved_params.yaml golden ====="
(
    cd "$WORKDIR"
    python3 "${PIPELINE_HOME}/scripts/auto_params.py" \
        --r1 "$FQ1" --r2 "$FQ2" \
        --output "${GOLDEN_DIR}/auto_params/resolved_params.yaml" \
        2>&1 | tee "${WORKDIR}/auto_params.log" \
        || log "  WARNING: auto_params.py invocation failed; check args"
)

# -----------------------------------------------------------------------------
# 14. Stage 6: MANIFEST.yaml generation
# -----------------------------------------------------------------------------

log "===== Stage 6: writing MANIFEST.yaml ====="

# Capture tool versions inside their respective containers/conda env. These
# calls are best-effort - any one that fails contributes an empty string.
get_version() {
    # Run a command and return its first line of stdout/stderr (untrimmed).
    "$@" 2>&1 | head -n1 || true
}

BWA_VER="$(get_version bwa 2>/dev/null | grep -i 'version' || true)"
SAMTOOLS_VER="$(get_version samtools --version 2>/dev/null || true)"
GATK_VER="$(get_version gatk --version 2>/dev/null || true)"
PICARD_VER="$(get_version picard MarkDuplicates --version 2>/dev/null || true)"
SAMBAMBA_VER="$(get_version sambamba --version 2>/dev/null || true)"
BCFTOOLS_VER="$(get_version bcftools --version 2>/dev/null || true)"

REF_SHA="$(sha256_of "$REFERENCE_FASTA")"
PONFA_SHA="$(sha256_of "$PONFA")"
GNOMAD_SHA="$(sha256_of "$GNOMAD")"
MASK1KG_SHA="$(sha256_of "$MASK1KG")"
FQ1_SHA="$(sha256_of "$FQ1")"
FQ2_SHA="$(sha256_of "$FQ2")"
SLURM_VER="$(sbatch --version 2>/dev/null | head -n1)"
HOSTNAME_FQDN="$(hostname -f 2>/dev/null || hostname)"
GEN_TIMESTAMP="$(date -u +%FT%TZ)"

# Use python3 to write YAML safely (avoids shell quoting hazards).
MANIFEST_PATH="${GOLDEN_DIR}/MANIFEST.yaml"
python3 - "$MANIFEST_PATH" <<'PYEOF'
import os, sys, datetime, json

manifest_path = sys.argv[1]
env = os.environ

data = {
    "spec": "SPEC-BSMN-REFACTOR-001",
    "spec_version": "1.2.0",
    "milestone": "M1b",
    "generated_at": env.get("GEN_TIMESTAMP", ""),
    "git_sha": env.get("GIT_SHA_AT_START", "unknown"),
    "reference": {
        "build":         "hg38_no_alt",
        "fasta_path":    env.get("REFERENCE_FASTA", ""),
        "fasta_sha256":  env.get("REF_SHA", ""),
    },
    "tools": {
        "bwa":          env.get("BWA_VER", ""),
        "samtools":     env.get("SAMTOOLS_VER", ""),
        "gatk":         env.get("GATK_VER", ""),
        "picard":       env.get("PICARD_VER", ""),
        "sambamba":     env.get("SAMBAMBA_VER", ""),
        "bcftools":     env.get("BCFTOOLS_VER", ""),
        "apptainer":    env.get("APPTAINER_VER", ""),
    },
    "resources": {
        "gnomad_af_vcf":           env.get("GNOMAD", ""),
        "gnomad_af_vcf_sha256":    env.get("GNOMAD_SHA", ""),
        "pon_mask_fasta":          env.get("PONFA", ""),
        "pon_mask_fasta_sha256":   env.get("PONFA_SHA", ""),
        "kg_strict_mask":          env.get("MASK1KG", ""),
        "kg_strict_mask_sha256":   env.get("MASK1KG_SHA", ""),
    },
    "sample": {
        "id":         env.get("SAMPLE_NAME", ""),
        "r1":         "tests/data/sample_R1.fastq.gz",
        "r2":         "tests/data/sample_R2.fastq.gz",
        "r1_sha256":  env.get("FQ1_SHA", ""),
        "r2_sha256":  env.get("FQ2_SHA", ""),
    },
    "cluster": {
        "host":           env.get("HOSTNAME_FQDN", ""),
        "slurm_version":  env.get("SLURM_VER", ""),
        "partition":      env.get("QUEUE", ""),
        "conda_env":      env.get("CONDA_ENV", ""),
    },
    "runtime_seconds": int(env.get("SECONDS", "0")),
    "commands": {
        "stage1_mapping":   "python3 jobs/run_genome_mapping.py --queue $QUEUE --conda-env $CONDA_ENV --reference $REFERENCE --align-fmt cram --sample-list $SAMPLE_LIST",
        "stage2_calling":   "python3 jobs/run_variant_calling.py  --queue $QUEUE --conda-env $CONDA_ENV --reference $REFERENCE --align-fmt cram --run-gatk-hc 2 --sample-list $SAMPLE_LIST",
        "stage3_filtering": "python3 jobs/run_variant_filtering.py --queue $QUEUE --conda-env $CONDA_ENV --reference $REFERENCE --align-fmt cram --run-gatk-hc 2 --sample-list $SAMPLE_LIST",
    },
}

# Minimal hand-written YAML serializer to avoid PyYAML dependency. Keys are
# stable strings, values are str/int only -> safe to emit as JSON which is a
# valid YAML 1.2 subset.
with open(manifest_path, "w") as fh:
    fh.write("# SPEC-BSMN-REFACTOR-001 v1.2.0 M1b golden snapshot manifest\n")
    fh.write("# This file is auto-generated by scripts/m1b_generate_goldens.sh\n")
    fh.write(json.dumps(data, indent=2, sort_keys=False, default=str))
    fh.write("\n")

print(f"wrote {manifest_path}")
PYEOF

export GEN_TIMESTAMP GIT_SHA_AT_START REFERENCE_FASTA REF_SHA \
       BWA_VER SAMTOOLS_VER GATK_VER PICARD_VER SAMBAMBA_VER BCFTOOLS_VER APPTAINER_VER \
       GNOMAD GNOMAD_SHA PONFA PONFA_SHA MASK1KG MASK1KG_SHA \
       SAMPLE_NAME FQ1_SHA FQ2_SHA \
       HOSTNAME_FQDN SLURM_VER QUEUE CONDA_ENV SECONDS

# (Note: the python heredoc above reads from env; the export above ensures the
# variables propagate. If using set -u, all referenced env vars must exist - the
# defaults via .get("KEY","") cover any miss.)

# -----------------------------------------------------------------------------
# 15. Stage 7: summary report
# -----------------------------------------------------------------------------

log "===== Stage 7: summary ====="
log "Golden directory: $GOLDEN_DIR"
find "$GOLDEN_DIR" -type f -printf '  %p  (%s bytes)\n' 2>/dev/null \
    || find "$GOLDEN_DIR" -type f -exec ls -la {} \;
TOTAL_BYTES="$(du -sb "$GOLDEN_DIR" 2>/dev/null | cut -f1 || du -sk "$GOLDEN_DIR" | awk '{print $1*1024}')"
log "Total golden size: $TOTAL_BYTES bytes"
log "Total runtime:     ${SECONDS}s"

log "===================================================="
log "M1b golden generation COMPLETE."
log "Next steps:"
log "  1. cd $PIPELINE_HOME"
log "  2. git status   # review tests/data/golden/ additions"
log "  3. git add tests/data/golden/"
log "  4. git commit -m 'chore(M1b): populate golden snapshots (SPEC-BSMN-REFACTOR-001)'"
log "  5. git push origin refactor/snakemake"
log "  6. In a new Claude session: /moai run SPEC-BSMN-REFACTOR-001 M4"
log "===================================================="

CLEAN_ON_EXIT=1
exit 0
