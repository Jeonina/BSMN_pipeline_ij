"""Unit tests for scripts/build_ref_cache.py (htslib REF_CACHE layout + M5)."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
_spec = importlib.util.spec_from_file_location(
    "build_ref_cache", SCRIPTS_DIR / "build_ref_cache.py"
)
assert _spec and _spec.loader
brc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(brc)


def _fasta(tmp_path: Path) -> Path:
    fa = tmp_path / "ref.fa"
    # Mixed case, wrapped lines, CRLF and a header description: all must be
    # normalised away exactly as the SAM spec's @SQ M5 does.
    fa.write_bytes(b">chrA desc here\nacgtN\r\nACG\n>chrB\nTTTT\n\n")
    return fa


def test_md5_matches_sam_spec_normalisation(tmp_path):
    entries = brc.build(str(_fasta(tmp_path)), str(tmp_path / "cache"))
    assert entries == [
        ("chrA", 8, hashlib.md5(b"ACGTNACG").hexdigest()),
        ("chrB", 4, hashlib.md5(b"TTTT").hexdigest()),
    ]


def test_cache_layout_and_content(tmp_path):
    out = tmp_path / "cache"
    brc.build(str(_fasta(tmp_path)), str(out))
    md5 = hashlib.md5(b"ACGTNACG").hexdigest()
    f = out / md5[:2] / md5[2:4] / md5[4:]
    assert f.read_bytes() == b"ACGTNACG"
    assert not list(out.glob(".part.*")), "temporary files must not be left behind"


def test_rerun_is_idempotent(tmp_path):
    fa, out = _fasta(tmp_path), str(tmp_path / "cache")
    assert brc.build(str(fa), out) == brc.build(str(fa), out)


def test_ref_path_template(tmp_path):
    assert brc.ref_path_template("/x/cache") == "/x/cache/%2s/%2s/%s"
