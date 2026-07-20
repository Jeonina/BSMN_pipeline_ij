#!/usr/bin/env python3
"""PyPI package-name uniqueness check (REQ-DUP-005 / AC-20).

Hits ``https://pypi.org/pypi/<name>/json`` and asserts the response is
``404 Not Found`` — i.e. no existing distribution shadows the chosen
internal package name.

Offline mode (``--offline``) is a no-op success, useful in air-gapped
dev environments and during M1a scaffolding when network access is not
required. Full enforcement is the M2 / AC-20 gate.
"""

from __future__ import annotations

import argparse
import sys
import urllib.error
from urllib.request import urlopen

DEFAULT_NAME: str = "bsmn_pipeline"
PYPI_URL_TEMPLATE: str = "https://pypi.org/pypi/{name}/json"
_TIMEOUT_SECONDS: float = 10.0


def _die(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(2)


def _check_online(name: str) -> bool:
    """Return ``True`` if ``name`` is free on PyPI (HTTP 404), else die."""
    url = PYPI_URL_TEMPLATE.format(name=name)
    try:
        with urlopen(url, timeout=_TIMEOUT_SECONDS) as resp:
            status = getattr(resp, "status", 200)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return True
        _die(
            f"PyPI lookup for {name!r} returned HTTP {exc.code}: {exc.reason}. "
            "Cannot determine package-name uniqueness; re-run with --offline "
            "to skip this check, or investigate."
        )
        return False  # pragma: no cover
    except urllib.error.URLError as exc:
        _die(
            f"PyPI lookup for {name!r} failed: {exc.reason}. "
            "Re-run with --offline if the dev environment lacks network "
            "access."
        )
        return False  # pragma: no cover

    if status == 404:
        return True

    _die(
        f"PyPI lookup for {name!r} returned HTTP {status}; an existing "
        "distribution may collide with this internal package name. "
        "Rename to a lab-prefixed alternative (e.g. tjbaelab_bsmn_pipeline) "
        "or rerun with --offline if this is a known sandbox response."
    )
    return False  # pragma: no cover


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="check_pkg_name",
        description=(
            "Verify that the internal package name is free on PyPI "
            "(REQ-DUP-005). Use --offline in air-gapped environments."
        ),
    )
    parser.add_argument(
        "--name",
        default=DEFAULT_NAME,
        help=f"Package name to check (default: {DEFAULT_NAME!r})",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Skip the PyPI request and return 0 (logs 'skipped (offline)').",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    name: str = args.name

    if args.offline:
        print(f"OK: skipped (offline) — name={name}")
        return 0

    if _check_online(name):
        print(f"OK: {name} is not on PyPI (free)")
        return 0

    return 2  # pragma: no cover (unreachable; _check_online dies on failure)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
