"""Tests for scripts.check_pkg_name (REQ-DUP-005 / AC-20).

Verifies that the package name ``bsmn_pipeline`` does not collide with
an existing PyPI distribution. Online mode hits PyPI; offline mode is a
no-op success used in air-gapped dev environments.
"""

from __future__ import annotations

import http
import urllib.error
from typing import Any

import pytest
from scripts import check_pkg_name as c

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class _FakeResp:
    def __init__(self, status: int) -> None:
        self.status = status

    def __enter__(self) -> _FakeResp:
        return self

    def __exit__(self, *exc: object) -> None:
        return None

    def read(self) -> bytes:
        return b"{}"


def _make_urlopen_returning(status: int):
    def _fake(_url: str, timeout: float = 0) -> _FakeResp:
        return _FakeResp(status)

    return _fake


def _make_urlopen_raising(exc: BaseException):
    def _fake(_url: str, timeout: float = 0) -> Any:
        raise exc

    return _fake


# ---------------------------------------------------------------------------
# Offline mode
# ---------------------------------------------------------------------------


def test_offline_mode_returns_skipped(
    capsys: pytest.CaptureFixture[str],
) -> None:
    rc = c.main(["--offline"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "skipped" in out.lower()
    assert "offline" in out.lower()


# ---------------------------------------------------------------------------
# Online mode — name is free
# ---------------------------------------------------------------------------


def test_online_404_means_name_is_free(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(
        c,
        "urlopen",
        _make_urlopen_raising(
            urllib.error.HTTPError(
                url="https://pypi.org/pypi/bsmn_pipeline/json",
                code=404,
                msg="Not Found",
                hdrs=None,  # type: ignore[arg-type]
                fp=None,
            )
        ),
    )
    rc = c.main([])
    captured = capsys.readouterr()
    assert rc == 0
    assert "free" in captured.out.lower() or "not on pypi" in captured.out.lower()


# ---------------------------------------------------------------------------
# Online mode — collision
# ---------------------------------------------------------------------------


def test_online_200_means_collision(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(c, "urlopen", _make_urlopen_returning(200))
    with pytest.raises(SystemExit) as exc:
        c.main([])
    assert exc.value.code == 2


# ---------------------------------------------------------------------------
# Online mode — network / unexpected error
# ---------------------------------------------------------------------------


def test_online_url_error_exits_2(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        c,
        "urlopen",
        _make_urlopen_raising(urllib.error.URLError("dns failure")),
    )
    with pytest.raises(SystemExit) as exc:
        c.main([])
    assert exc.value.code == 2


def test_online_unexpected_http_status_exits_2(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        c,
        "urlopen",
        _make_urlopen_raising(
            urllib.error.HTTPError(
                url="https://pypi.org/pypi/bsmn_pipeline/json",
                code=int(http.HTTPStatus.INTERNAL_SERVER_ERROR),
                msg="Server Error",
                hdrs=None,  # type: ignore[arg-type]
                fp=None,
            )
        ),
    )
    with pytest.raises(SystemExit) as exc:
        c.main([])
    assert exc.value.code == 2


# ---------------------------------------------------------------------------
# Custom package name
# ---------------------------------------------------------------------------


def test_custom_name_arg(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(
        c,
        "urlopen",
        _make_urlopen_raising(
            urllib.error.HTTPError(
                url="https://pypi.org/pypi/tjbaelab_bsmn_pipeline/json",
                code=404,
                msg="Not Found",
                hdrs=None,  # type: ignore[arg-type]
                fp=None,
            )
        ),
    )
    rc = c.main(["--name", "tjbaelab_bsmn_pipeline"])
    captured = capsys.readouterr()
    assert rc == 0
    assert "tjbaelab_bsmn_pipeline" in captured.out
