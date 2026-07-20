import os
import re
import shutil
import subprocess
import sys
from collections.abc import Generator
from typing import Any

from .config import read_config
from .misc import coroutine

# By default, use a samtools in PATH.
SAMTOOLS: str | None = shutil.which("samtools")


def load_config(reference: str, conda_env: str) -> None:
    global SAMTOOLS
    config = read_config(reference, conda_env)
    SAMTOOLS = config["TOOLS"]["SAMTOOLS"]
    if SAMTOOLS is None or not os.path.isfile(SAMTOOLS) or not os.access(SAMTOOLS, os.X_OK):
        SAMTOOLS = shutil.which("samtools")


def base_count(bam: str, min_MQ: int, min_BQ: int) -> Generator[Any, Any, Any]:
    return pileup(bam, min_MQ, min_BQ, base_n())


def base_qual_tuple(bam: str, min_MQ: int, min_BQ: int) -> Generator[Any, Any, Any]:
    return pileup(bam, min_MQ, min_BQ, base_qual())


@coroutine
def pileup(
    bam: str, min_MQ: int, min_BQ: int, target: Generator[Any, Any, Any]
) -> Generator[Any, tuple[str, int]]:
    result = None
    while True:
        chrom, pos = yield result
        cmd = [
            str(SAMTOOLS),
            "mpileup",
            "-d",
            "8000",
            "-q",
            str(min_MQ),
            "-Q",
            str(min_BQ),
            "-r",
            f"{chrom}:{pos}-{pos}",
        ] + bam.split()

        max_retries = 5
        n_retries = 0
        while True:
            cmd_out = subprocess.run(
                cmd,
                text=True,
                capture_output=True,
                check=False,
            )
            try:
                cmd_out.check_returncode()
                break
            except subprocess.CalledProcessError:
                if n_retries == 0:
                    sys.stderr.write("\n")
                sys.stderr.write(
                    f"{cmd_out.stderr}source bam: {bam}\npileup location: {chrom}:{pos}-{pos}\n"
                )
                if n_retries < max_retries:
                    n_retries += 1
                    sys.stderr.write("Retry pileup...\n")
                else:
                    sys.exit("Failed in pileup.")
        try:
            pileup_outs = cmd_out.stdout.split()
            bases = ""
            quals = ""
            for i in range(4, len(pileup_outs), 3):
                _bases, _quals = pileup_outs[i : i + 2]
                bases += bases_clean(_bases)
                quals += _quals
        except ValueError:
            bases, quals = ("", "")
        result = target.send((bases, quals))


def bases_clean(bases: str) -> str:
    bases = re.sub(r"\^.", "", bases)
    bases = re.sub(r"\$", "", bases)
    for n in set(re.findall(r"-(\d+)", bases)):
        bases = re.sub(rf"-{n}[ACGTNacgtn]{{{n}}}", "", bases)
    for n in set(re.findall(r"\+(\d+)", bases)):
        bases = re.sub(rf"\+{n}[ACGTNacgtn]{{{n}}}", "", bases)
    return bases


@coroutine
def base_n() -> Generator[dict[str, int], tuple[str, str]]:
    result: dict[str, int] | None = None
    while True:
        bases, quals = yield result  # type: ignore[misc]
        base_n_dict: dict[str, int] = {}
        base_n_dict["A"] = bases.count("A")
        base_n_dict["C"] = bases.count("C")
        base_n_dict["G"] = bases.count("G")
        base_n_dict["T"] = bases.count("T")
        base_n_dict["a"] = bases.count("a")
        base_n_dict["c"] = bases.count("c")
        base_n_dict["g"] = bases.count("g")
        base_n_dict["t"] = bases.count("t")
        base_n_dict["dels"] = bases.count("*")
        result = base_n_dict


@coroutine
def base_qual() -> Generator[list[tuple[str, int]], tuple[str, str]]:
    result: list[tuple[str, int]] | None = None
    while True:
        bases, quals = yield result  # type: ignore[misc]
        bases = re.sub(r"\*", "", bases)
        result = list(map(lambda b, q: (b.upper(), ord(q) - 33), bases, quals))
