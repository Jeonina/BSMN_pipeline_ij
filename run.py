#!/usr/bin/env python3
"""
BSMN Somatic Variant Calling Pipeline — one-command entry point.

Usage
-----
  # FASTQ 디렉토리 (가장 일반적)
  python run.py /data/fastq/

  # R1 파일 하나 (R2는 파일명 규칙으로 자동 탐색)
  python run.py /data/WES_R1.fastq.gz

  # R1 + R2 명시
  python run.py /data/WES_R1.fastq.gz /data/WES_R2.fastq.gz

  # 옵션 예시
  python run.py /data/fastq/ --cores 40 --stage mapping
  python run.py /data/fastq/ --dry-run

  # 하위 디렉토리까지 재귀 탐색
  python run.py /data/fastq/ --recursive
"""

import argparse
import shutil
import subprocess
import sys
import os
from pathlib import Path

# scripts/ 디렉토리를 경로에 추가 (make_samples_tsv 직접 임포트)
_SCRIPTS = Path(__file__).resolve().parent / "scripts"
sys.path.insert(0, str(_SCRIPTS))
from make_samples_tsv import build_table, write_tsv, find_r2, R_MARKERS

# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────

def _row_for_pair(r1: Path, r2: Path, pattern: str | None) -> dict:
    """Build a single samples.tsv row from an explicit R1/R2 pair.

    @MX:NOTE: bypasses build_table's directory scan to prevent the
    Bug 2 leak (sibling FASTQ pairs in the same directory contaminating
    samples.tsv). Reuses extract_sample_rg for the naming convention.
    """
    from make_samples_tsv import extract_sample_rg

    sample_id, rg = extract_sample_rg(
        r1,
        custom_pattern=pattern,
        split_fields=None,
        split_delim="_",
    )
    return {
        "sample_id": sample_id,
        "readgroup": rg,
        "fq1": str(r1.resolve()),
        "fq2": str(r2.resolve()),
    }


def _resolve_inputs(inputs: list[str], recursive: bool, pattern: str | None) -> list[dict]:
    """입력 형태를 자동 판별하고 샘플 row 리스트를 반환."""
    p0 = Path(inputs[0])

    # Case 1: 디렉토리 — scan all R1 files in the directory
    if len(inputs) == 1 and p0.is_dir():
        return build_table(
            directory=str(p0),
            recursive=recursive,
            custom_pattern=pattern,
            split_fields=None,
            split_delim="_",
        )

    # Case 2: R1 파일 하나 → R2 자동 탐색 (single pair only, no directory scan)
    if len(inputs) == 1 and p0.is_file():
        r2 = find_r2(p0)
        if r2 is None:
            _die(f"R2 파일을 찾을 수 없습니다: {p0}")
        assert r2 is not None  # mypy: _die() never returns
        return [_row_for_pair(p0, r2, pattern)]

    # Case 3: R1 + R2 명시 (single pair only, no directory scan)
    if len(inputs) == 2:
        r1, r2 = Path(inputs[0]), Path(inputs[1])
        if not r1.is_file():
            _die(f"파일 없음: {r1}")
        if not r2.is_file():
            _die(f"파일 없음: {r2}")
        return [_row_for_pair(r1, r2, pattern)]

    _die("입력을 인식할 수 없습니다. --help 참고.")
    return []  # unreachable; satisfies type checker


def _die(msg: str) -> None:
    print(f"[run.py] ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


# @MX:ANCHOR: runtime_minutes_to_hms — used by SLURM profile docs and tests
# @MX:REASON: stable contract between Snakemake `resources.runtime` (minutes,
#             integer) and SLURM `--time=HH:MM:SS`. Centralized to avoid
#             format drift between docs, tests, and any future helper scripts.
def runtime_minutes_to_hms(minutes: int) -> str:
    """Convert integer minutes to SLURM-compatible ``HH:MM:SS`` string.

    Examples:
        >>> runtime_minutes_to_hms(30)
        '00:30:00'
        >>> runtime_minutes_to_hms(1440)
        '24:00:00'
    """
    if minutes < 0:
        raise ValueError(f"minutes must be non-negative, got {minutes}")
    hours, mins = divmod(int(minutes), 60)
    return f"{hours:02d}:{mins:02d}:00"


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="BSMN 파이프라인 — 한 줄 실행",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # ── 필수: 입력 경로 ──────────────────────────────────────────────────────
    parser.add_argument(
        "input",
        nargs="+",
        metavar="INPUT",
        help=(
            "FASTQ 디렉토리, R1 파일, 또는 R1+R2 파일 쌍.\n"
            "예) /data/fastq/  |  sample_R1.fastq.gz  |  R1.fq.gz R2.fq.gz"
        ),
    )

    # ── 실행 옵션 ────────────────────────────────────────────────────────────
    parser.add_argument(
        "--stage",
        choices=["mapping", "calling", "filtering", "all"],
        default="all",
        help="실행할 파이프라인 단계 (기본값: all)",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=80,
        help="사용할 CPU 코어 수 (기본값: 80, --cluster slurm 사용 시 무시됨)",
    )
    parser.add_argument(
        "--cluster",
        choices=["none", "slurm"],
        default="none",
        help=(
            "실행 백엔드 선택 (기본값: none = 로컬 --cores). "
            "'slurm' 사용 시 workflow/profiles/slurm 프로파일이 적용되며 "
            "sbatch가 PATH에 있어야 합니다."
        ),
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        dest="dry_run",
        help="Snakemake dry-run (실제 실행하지 않고 계획만 출력)",
    )

    # ── 샘플 탐색 옵션 ───────────────────────────────────────────────────────
    parser.add_argument(
        "--recursive", "-r",
        action="store_true",
        help="FASTQ 디렉토리를 하위 폴더까지 재귀 탐색",
    )
    parser.add_argument(
        "--pattern",
        default=None,
        metavar="REGEX",
        help="샘플명 추출용 커스텀 정규식 (make_samples_tsv --pattern 과 동일)",
    )

    # ── Snakemake 추가 인자 passthrough ─────────────────────────────────────
    parser.add_argument(
        "--snakemake-args",
        default="",
        metavar="ARGS",
        help='추가 Snakemake 인자 (예: "--rerun-incomplete --keep-going")',
    )

    args = parser.parse_args()

    # ── 1. samples.tsv 생성 ─────────────────────────────────────────────────
    print("\n[run.py] ① 샘플 목록 생성 중...")

    # Defensive: remove any stale samples.tsv from a previous run before
    # building the new one. write_tsv() already opens in "w" mode, but
    # explicit deletion makes the contract obvious and protects against
    # an interrupted run leaving partial content.
    samples_tsv = Path("config/samples.tsv")
    if samples_tsv.exists():
        samples_tsv.unlink()

    rows = _resolve_inputs(args.input, args.recursive, args.pattern)

    if not rows:
        _die("유효한 FASTQ 쌍을 찾지 못했습니다.")

    write_tsv(rows, "config/samples.tsv")

    # ── 2. 이전 resolved_params.yaml 삭제 (샘플 바뀌면 재생성 필요) ─────────
    resolved = Path("config/resolved_params.yaml")
    if resolved.exists():
        resolved.unlink()
        print("[run.py]    기존 resolved_params.yaml 삭제 (Snakemake가 재생성)")

    # ── 3. Snakemake 실행 ───────────────────────────────────────────────────
    cmd = [
        "snakemake",
        "--snakefile", "workflow/Snakefile",
        "--config", f"stage={args.stage}",
    ]

    if args.cluster == "slurm":
        # Validate sbatch is available — fail fast with a clear message rather
        # than letting snakemake hit a cryptic plugin error mid-pipeline.
        if shutil.which("sbatch") is None:
            _die(
                "sbatch가 PATH에서 검출되지 않았습니다. SLURM 클러스터 노드에서 "
                "실행하거나 SLURM 클라이언트 도구를 설치하세요."
            )
        cmd.extend(["--profile", "workflow/profiles/slurm"])
        print("[run.py]    Running on SLURM cluster (profile: workflow/profiles/slurm)")
    else:
        cmd.extend(["--cores", str(args.cores)])

    if args.dry_run:
        cmd.append("--dry-run")
    if args.snakemake_args:
        cmd.extend(args.snakemake_args.split())

    print(f"\n[run.py] ② Snakemake 실행:")
    print(f"   {' '.join(cmd)}\n")

    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
