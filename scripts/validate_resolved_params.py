#!/usr/bin/env python3
"""Validate resolved_params.yaml against the REQ-SNK-009 schema.

Required top-level keys (PIPELINE.md §169-177):

  - ``resolved_at``                       — str (ISO-8601 timestamp)
  - ``input_fastq``                       — str (absolute path)
  - ``sequencer``                         — str
  - ``bwa_threads``                       — int (>= 1)
  - ``optical_duplicate_pixel_distance``  — int (>= 0)
  - ``bqsr_memory_gb``                    — int (>= 1)

Additional keys MAY be present; the validator does not reject them.

On any violation the script exits with status 2 and writes a clear
diagnostic to stderr.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class ResolvedParams:
    """Validated resolved_params payload (only the required keys)."""

    resolved_at: str
    input_fastq: str
    sequencer: str
    bwa_threads: int
    optical_duplicate_pixel_distance: int
    bqsr_memory_gb: int


def _die(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(2)


def _require_str(payload: dict[str, Any], key: str) -> str:
    if key not in payload:
        _die(f"resolved_params.yaml: missing required key '{key}'")
    value = payload[key]
    if not isinstance(value, str) or isinstance(value, bool):
        _die(f"resolved_params.yaml: key '{key}' must be a string, got {type(value).__name__}")
    if not value.strip():
        _die(f"resolved_params.yaml: key '{key}' must be a non-empty string")
    return str(value)


def _require_int(payload: dict[str, Any], key: str, *, minimum: int) -> int:
    if key not in payload:
        _die(f"resolved_params.yaml: missing required key '{key}'")
    value = payload[key]
    # In Python, bool is a subclass of int — reject explicitly.
    if isinstance(value, bool) or not isinstance(value, int):
        _die(f"resolved_params.yaml: key '{key}' must be an integer, got {type(value).__name__}")
    if value < minimum:
        _die(f"resolved_params.yaml: key '{key}' must be >= {minimum}, got {value}")
    return int(value)


def _load_yaml(path: Path | str) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        _die(f"resolved_params.yaml not found: {p}")
    try:
        with p.open("r", encoding="utf-8") as fh:
            payload = yaml.safe_load(fh)
    except yaml.YAMLError as exc:
        _die(f"resolved_params.yaml: malformed YAML: {exc}")
    if not isinstance(payload, dict):
        _die(
            f"resolved_params.yaml: top-level value must be a mapping, got {type(payload).__name__}"
        )
    return dict(payload)


def validate(path: Path | str) -> ResolvedParams:
    """Validate ``path`` and return a :class:`ResolvedParams` instance.

    Exits with status 2 on any violation.
    """
    payload = _load_yaml(path)
    return ResolvedParams(
        resolved_at=_require_str(payload, "resolved_at"),
        input_fastq=_require_str(payload, "input_fastq"),
        sequencer=_require_str(payload, "sequencer"),
        bwa_threads=_require_int(payload, "bwa_threads", minimum=1),
        optical_duplicate_pixel_distance=_require_int(
            payload, "optical_duplicate_pixel_distance", minimum=0
        ),
        bqsr_memory_gb=_require_int(payload, "bqsr_memory_gb", minimum=1),
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="validate_resolved_params",
        description=("Validate resolved_params.yaml against the REQ-SNK-009 schema."),
    )
    parser.add_argument("path", help="Path to resolved_params.yaml")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    params = validate(args.path)
    print(
        f"OK: resolved_params.yaml valid (sequencer={params.sequencer})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
