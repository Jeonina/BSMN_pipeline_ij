#!/usr/bin/env python3
"""Apptainer / SingularityCE runtime pre-flight check (REQ-SNK-008).

Required floor:

  - Apptainer       >= 1.2.5
  - SingularityCE   >= 3.11.0   (fallback when ``apptainer`` is absent)

The check fails fast with status 2 and a clear diagnostic naming both
expected version floors when neither runtime is available, or when an
available runtime is below its floor.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass

APPTAINER_FLOOR: tuple[int, int, int] = (1, 2, 5)
SINGULARITY_FLOOR: tuple[int, int, int] = (3, 11, 0)

_VERSION_RE = re.compile(r"(\d+)\.(\d+)\.(\d+)")


@dataclass(frozen=True)
class PreflightResult:
    """Outcome of a successful pre-flight check."""

    tool: str  # "apptainer" or "singularity"
    version: tuple[int, int, int]


def _die(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(2)


def parse_version(raw: str) -> tuple[int, int, int]:
    """Extract the first ``X.Y.Z`` triple from ``raw``.

    Raises :class:`ValueError` if no triple is found.
    """
    match = _VERSION_RE.search(raw)
    if not match:
        raise ValueError(f"cannot parse version from: {raw!r}")
    return (int(match.group(1)), int(match.group(2)), int(match.group(3)))


def _query(tool: str) -> tuple[int, int, int] | None:
    """Return the parsed version of ``tool``, or ``None`` if not in PATH."""
    try:
        proc = subprocess.run(
            [tool, "--version"],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
    except FileNotFoundError:
        return None
    except subprocess.TimeoutExpired:
        _die(f"{tool} --version timed out")
        return None  # pragma: no cover
    raw = (proc.stdout or proc.stderr or "").strip()
    try:
        return parse_version(raw)
    except ValueError as exc:
        _die(f"{tool}: {exc}")
        return None  # pragma: no cover


def preflight() -> PreflightResult:
    """Verify Apptainer or SingularityCE meets the version floor.

    Exits with status 2 if neither runtime is present, or if the available
    runtime is below its floor.
    """
    apptainer = _query("apptainer")
    if apptainer is not None:
        if apptainer < APPTAINER_FLOOR:
            _die(
                "Apptainer "
                f"{apptainer[0]}.{apptainer[1]}.{apptainer[2]} is below the "
                f"required floor {APPTAINER_FLOOR[0]}.{APPTAINER_FLOOR[1]}."
                f"{APPTAINER_FLOOR[2]} (or SingularityCE >= "
                f"{SINGULARITY_FLOOR[0]}.{SINGULARITY_FLOOR[1]}."
                f"{SINGULARITY_FLOOR[2]})"
            )
        return PreflightResult(tool="apptainer", version=apptainer)

    singularity = _query("singularity")
    if singularity is not None:
        if singularity < SINGULARITY_FLOOR:
            _die(
                "SingularityCE "
                f"{singularity[0]}.{singularity[1]}.{singularity[2]} is "
                f"below the required floor "
                f"{SINGULARITY_FLOOR[0]}.{SINGULARITY_FLOOR[1]}."
                f"{SINGULARITY_FLOOR[2]} (or Apptainer >= "
                f"{APPTAINER_FLOOR[0]}.{APPTAINER_FLOOR[1]}."
                f"{APPTAINER_FLOOR[2]})"
            )
        return PreflightResult(tool="singularity", version=singularity)

    _die(
        "neither 'apptainer' nor 'singularity' was found in PATH. "
        f"Install Apptainer >= {APPTAINER_FLOOR[0]}.{APPTAINER_FLOOR[1]}."
        f"{APPTAINER_FLOOR[2]} or SingularityCE >= "
        f"{SINGULARITY_FLOOR[0]}.{SINGULARITY_FLOOR[1]}."
        f"{SINGULARITY_FLOOR[2]} (see https://apptainer.org/docs/)."
    )
    # Unreachable; _die exits, but mypy needs a return path.
    return PreflightResult(tool="", version=(0, 0, 0))  # pragma: no cover


def _build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        prog="apptainer_preflight",
        description=(
            "Verify Apptainer >= "
            f"{APPTAINER_FLOOR[0]}.{APPTAINER_FLOOR[1]}.{APPTAINER_FLOOR[2]} "
            "or SingularityCE >= "
            f"{SINGULARITY_FLOOR[0]}.{SINGULARITY_FLOOR[1]}."
            f"{SINGULARITY_FLOOR[2]} is installed."
        ),
    )


def main(argv: list[str] | None = None) -> int:
    _build_parser().parse_args(argv)
    result = preflight()
    major, minor, patch = result.version
    print(f"OK: {result.tool} {major}.{minor}.{patch}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
