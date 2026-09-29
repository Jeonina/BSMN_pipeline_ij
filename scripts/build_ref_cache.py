#!/usr/bin/env python3
"""Build an htslib REF_CACHE (MD5-keyed reference sequences) from a FASTA.

Why: a CRAM stores each contig's MD5 plus the reference path it was written
with (@SQ UR:). Any tool that decodes the CRAM WITHOUT being handed the
reference looks up UR: - usually a path on whichever host did the mapping -
and, when that is missing, falls back to downloading the sequence from the EBI
reference server. On a host without outbound access that download does not
fail, it hangs until the TCP timeout. MosaicForecast's feature extraction runs
exactly such a call (`samtools view <cram> | head`, no -T) once per candidate,
which cost ~10-15 min of pure waiting per variant.

With REF_PATH/REF_CACHE pointing at this cache, htslib resolves every contig
locally and never touches the network. The layout is htslib's own
(seq_cache_populate.pl): <out>/<md5[0:2]>/<md5[2:4]>/<md5[4:]>, each file the
contig's sequence uppercased with all whitespace removed - the same bytes the
SAM spec hashes for @SQ M5. Contigs are hashed while streaming, so memory stays
flat regardless of contig size. Existing entries are kept (safe to re-run).

Writes <out>/manifest.tsv (name, length, md5) LAST, as the completion marker.

Usage:
    python scripts/build_ref_cache.py --ref resources/hg38/ref.fasta \
        --out resources/hg38/ref_cache
    export REF_PATH=resources/hg38/ref_cache/%2s/%2s/%s REF_CACHE=$REF_PATH
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import os
import sys
import tempfile

log = logging.getLogger("build_ref_cache")

# SAM spec M5: uppercase, only printable non-space bytes (33-126) are hashed.
_DROP = bytes(range(0, 33)) + bytes(range(127, 256))


def cache_path(out_dir: str, md5: str) -> str:
    return os.path.join(out_dir, md5[:2], md5[2:4], md5[4:])


def normalise(line: bytes) -> bytes:
    """One FASTA sequence line as htslib hashes it."""
    return line.translate(None, _DROP).upper()


def ref_path_template(out_dir: str) -> str:
    """REF_PATH / REF_CACHE value for this cache."""
    return os.path.join(out_dir, "%2s", "%2s", "%s")


def _finish(out_dir: str, name: str, tmp, md5, length: int) -> tuple[str, int, str]:
    tmp.close()
    digest = md5.hexdigest()
    dest = cache_path(out_dir, digest)
    if os.path.exists(dest):
        os.remove(tmp.name)
    else:
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        os.replace(tmp.name, dest)
    log.info("%s\t%d\t%s", name, length, digest)
    return name, length, digest


def build(fasta: str, out_dir: str) -> list[tuple[str, int, str]]:
    os.makedirs(out_dir, exist_ok=True)
    entries: list[tuple[str, int, str]] = []
    name, tmp, md5, length = None, None, None, 0
    with open(fasta, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    entries.append(_finish(out_dir, name, tmp, md5, length))
                name = line[1:].split()[0].decode()
                tmp = tempfile.NamedTemporaryFile(dir=out_dir, prefix=".part.", delete=False)
                md5, length = hashlib.md5(), 0
                continue
            if name is None:
                continue
            seq = normalise(line)
            md5.update(seq)
            tmp.write(seq)
            length += len(seq)
    if name is not None:
        entries.append(_finish(out_dir, name, tmp, md5, length))
    return entries


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", required=True, help="reference FASTA (uncompressed)")
    ap.add_argument("--out", required=True, help="cache directory")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    entries = build(args.ref, args.out)
    manifest = os.path.join(args.out, "manifest.tsv")
    with open(manifest + ".tmp", "w") as fh:
        fh.write("name\tlength\tmd5\n")
        for name, length, digest in entries:
            fh.write(f"{name}\t{length}\t{digest}\n")
    os.replace(manifest + ".tmp", manifest)
    log.info("%d contigs -> %s (REF_PATH=%s)", len(entries), args.out, ref_path_template(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
