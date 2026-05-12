#!/usr/bin/env python3
"""Reference identity contract check (REQ-SNK-006 / AC-15).

For all CI golden snapshots and characterization tests, the reference
genome **shall** be ``hg38_no_alt``. The aliases ``hg38_decoy``,
``hg38_v0``, and bare ``hg38`` are explicitly rejected here.

Configuration sources, checked in this order:

  1. ``BSMN_REFERENCE`` environment variable (highest precedence)
  2. ``BSMN_REFERENCE_FILE`` environment variable → path to a one-line
     marker file (e.g. ``config/ci_reference.txt``)

If neither is set, the check fails fast.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

EXPECTED_REFERENCE: str = "hg38_no_alt"
ENV_VAR: str = "BSMN_REFERENCE"
ENV_FILE_VAR: str = "BSMN_REFERENCE_FILE"


def _die(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(2)


def _read_marker_file(path: str) -> str:
    p = Path(path)
    if not p.exists():
        _die(f"{ENV_FILE_VAR} points to missing file: {p}")
    return p.read_text(encoding="utf-8").strip()


def _resolve_source() -> tuple[str, str]:
    """Return ``(source_description, value)`` for the configured marker."""
    env_value = os.environ.get(ENV_VAR)
    if env_value is not None:
        return (f"env:{ENV_VAR}", env_value.strip())

    file_path = os.environ.get(ENV_FILE_VAR)
    if file_path is not None:
        return (f"file:{file_path}", _read_marker_file(file_path))

    _die(
        "reference identity not configured: set either "
        f"${ENV_VAR}={EXPECTED_REFERENCE!r} or "
        f"${ENV_FILE_VAR}=/path/to/ci_reference.txt"
    )
    return ("", "")  # pragma: no cover


def check() -> str:
    """Verify the configured reference equals :data:`EXPECTED_REFERENCE`.

    Returns the validated reference string.
    Exits with status 2 on mismatch.
    """
    source, value = _resolve_source()
    if value != EXPECTED_REFERENCE:
        _die(
            f"reference identity mismatch from {source}: "
            f"expected {EXPECTED_REFERENCE!r}, got {value!r}. "
            "Aliases such as 'hg38_decoy', 'hg38_v0', or bare 'hg38' are "
            "not permitted in CI (REQ-SNK-006)."
        )
    return value


def _build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        prog="check_reference_identity",
        description=(
            f"Verify the CI reference genome marker equals {EXPECTED_REFERENCE!r} (REQ-SNK-006)."
        ),
    )


def main(argv: list[str] | None = None) -> int:
    _build_parser().parse_args(argv)
    value = check()
    print(f"OK: reference identity = {value}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
