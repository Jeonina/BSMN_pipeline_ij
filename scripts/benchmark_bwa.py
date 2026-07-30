#!/usr/bin/env python3
"""
benchmark_bwa.py — Empirically find the best bwa-mem per-job thread width for a host.

``bwa_threads`` is the one [HEURISTIC] knob in the auto-tuner with no published
optimum (see scripts/auto_params.py :func:`derive_bwa_threads`). This harness
measures it on a fixed subsample, in two modes:

    single      one bwa|sort job, vary -t              -> per-job scaling curve
                                                           (validates the "~16 knee")
    throughput  at fixed --cores, run cores//-t jobs    -> host throughput per split
                in parallel, vary -t                       (the decision-relevant metric)

It reports the split that maximises reads/second and compares it to what
``derive_bwa_threads`` would pick. Runs standalone (not via Snakemake) using the
same apptainer bwa/sambamba SIFs + reference as the pipeline, subsampling ONCE and
reusing it, so each config is minutes — not a full WGS map.

Usage:
    python scripts/benchmark_bwa.py \
        --fq1 /data/.../R1.fq.gz --fq2 /data/.../R2.fq.gz \
        --ref resources/hg38/Homo_sapiens_assembly38.fasta \
        --containers config/containers.yaml \
        --cores 24 --subsample 15000000 \
        --threads-list 6,8,12,16,24 --mode both --reps 2 \
        --scratch /tmp/bwabench \
        --out benchmarks/results/bwa_host.tsv

    # preview the run matrix + rough time, WITHOUT running anything:
    python scripts/benchmark_bwa.py ... --dry-run
"""

import argparse
import os
import shlex
import statistics
import subprocess
import threading
import time
from typing import Any

import yaml

_READS_PER_PAIR = 2  # a "read pair" = 2 reads; flagstat/throughput count reads


# ---------------------------------------------------------------------------
# Pure planning / summary helpers (unit-tested; no host or bwa dependency)
# ---------------------------------------------------------------------------


def plan_runs(threads_list: list[int], cores: int, modes: list[str], reps: int) -> list[dict[str, Any]]:
    """Expand the sweep into individual runs. In throughput mode a config uses
    ``concurrency = max(1, cores // threads)`` parallel jobs; single mode uses 1."""
    runs: list[dict[str, Any]] = []
    for mode in modes:
        for threads in threads_list:
            concurrency = 1 if mode == "single" else max(1, cores // threads)
            for rep in range(1, reps + 1):
                runs.append(
                    {"mode": mode, "threads": threads, "concurrency": concurrency, "rep": rep}
                )
    return runs


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Median reads/s per (mode, threads); pick the throughput-mode argmax."""
    by_key: dict[tuple[str, int], list[float]] = {}
    for r in rows:
        by_key.setdefault((r["mode"], r["threads"]), []).append(r["throughput_reads_s"])
    medians = {k: statistics.median(v) for k, v in by_key.items()}
    tput = {threads: mt for (mode, threads), mt in medians.items() if mode == "throughput"}
    best_threads = max(tput, key=tput.get) if tput else None
    return {"medians": medians, "throughput_by_threads": tput, "best_threads": best_threads}


# ---------------------------------------------------------------------------
# Container / IO helpers
# ---------------------------------------------------------------------------


def load_sifs(containers_path: str) -> tuple[str, str]:
    with open(containers_path) as fh:
        c = yaml.safe_load(fh)
    return c["bwa"]["sif"], c["sambamba"]["sif"]


def bind_flags(paths: list[str]) -> list[str]:
    """--bind for each unique existing parent dir, so apptainer can see inputs on
    NFS mounts (e.g. /data) and the scratch dir."""
    dirs = sorted({os.path.dirname(os.path.abspath(p)) for p in paths if p})
    flags: list[str] = []
    for d in dirs:
        if os.path.isdir(d):
            flags += ["--bind", d]
    return flags


def subsample(fq1: str, fq2: str, n_pairs: int, scratch: str) -> tuple[str, str]:
    """First ``n_pairs`` pairs written UNCOMPRESSED to scratch (removes the
    decompression variable and per-config subsample cost). n_pairs<=0 -> originals.
    Deterministic (head), which is fine for timing."""
    if n_pairs <= 0:
        return fq1, fq2
    os.makedirs(scratch, exist_ok=True)
    # Key the file by pair count so a different --subsample never silently reuses a
    # stale (wrong-size) subsample left in the same scratch dir.
    out1 = os.path.join(scratch, f"sub_{n_pairs}_R1.fastq")
    out2 = os.path.join(scratch, f"sub_{n_pairs}_R2.fastq")
    n_lines = n_pairs * 4
    for src, dst in ((fq1, out1), (fq2, out2)):
        if os.path.exists(dst) and os.path.getsize(dst) > 0:
            print(f"[bench] reuse existing subsample {dst}")
            continue
        print(f"[bench] subsampling {n_pairs:,} pairs -> {dst}")
        cat = "zcat" if src.endswith(".gz") else "cat"
        with open(dst, "wb") as out:
            p1 = subprocess.Popen([cat, src], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["head", "-n", str(n_lines)], stdin=p1.stdout, stdout=out)
            p1.stdout.close()
            p2.communicate()
            p1.wait()
    return out1, out2


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------


def _peak_bwa_rss_gb(stop: threading.Event) -> list[float]:
    """Poll summed RSS of all live `bwa` processes; return [peak_gb]. Best-effort;
    works across containers because bwa runs as a host-visible process."""
    peak = 0.0
    while not stop.is_set():
        try:
            pids = subprocess.run(
                ["pgrep", "-x", "bwa"], capture_output=True, text=True
            ).stdout.split()
            total_kb = 0
            for pid in pids:
                try:
                    with open(f"/proc/{pid}/status") as fh:
                        for line in fh:
                            if line.startswith("VmRSS:"):
                                total_kb += int(line.split()[1])
                                break
                except OSError:
                    pass
            peak = max(peak, total_kb / (1024 * 1024))
        except Exception:
            pass
        time.sleep(2)
    return [peak]


def build_cmd(
    threads: int, sub1: str, sub2: str, ref: str, bwa_sif: str, samb_sif: str,
    out_bam: str, fix_k: int, pipe: str, sort_mem: str, binds: list[str],
) -> str:
    sort_threads = max(2, threads // 4)
    k_flag = f"-K {fix_k} " if fix_k else ""
    rg = r"@RG\tID:bench\tSM:bench\tPL:illumina\tLB:bench\tPU:bench"
    b = " ".join(binds)
    bwa = (
        f"apptainer exec {b} {bwa_sif} bwa mem -M -t {threads} {k_flag}"
        f"-R '{rg}' {ref} {sub1} {sub2}"
    )
    if pipe == "bwa-only":
        return f"{bwa} > /dev/null"
    return (
        f"{bwa} "
        f"| apptainer exec {b} {samb_sif} sambamba view -S -f bam -l 0 /dev/stdin "
        f"| apptainer exec {b} {samb_sif} sambamba sort -m {sort_mem} "
        f"-t {sort_threads} -o {out_bam} /dev/stdin"
    )


def run_config(cmd_builder, concurrency: int, scratch: str) -> tuple[float, float]:
    """Launch ``concurrency`` identical jobs (distinct outputs), wait for all.
    Returns (wall_s, peak_rss_gb). wall = slowest job."""
    stop = threading.Event()
    peak_holder: list[float] = [0.0]

    def poll():
        peak_holder[0] = _peak_bwa_rss_gb(stop)[0]

    poller = threading.Thread(target=poll, daemon=True)
    poller.start()

    procs = []
    t0 = time.monotonic()
    for i in range(concurrency):
        out_bam = os.path.join(scratch, f"bench_out_{i}.bam")
        cmd = cmd_builder(out_bam)
        procs.append(subprocess.Popen(["bash", "-c", f"set -o pipefail; {cmd}"]))
    rc = [p.wait() for p in procs]
    wall = time.monotonic() - t0
    stop.set()
    poller.join(timeout=5)

    for i in range(concurrency):
        for suf in (".bam", ".bam.bai"):
            f = os.path.join(scratch, f"bench_out_{i}{suf}")
            if os.path.exists(f):
                os.remove(f)
    if any(c != 0 for c in rc):
        raise RuntimeError(f"a bwa job exited non-zero: {rc}")
    return wall, peak_holder[0]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fq1", required=True)
    ap.add_argument("--fq2", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--containers", default="config/containers.yaml")
    ap.add_argument("--cores", type=int, required=True, help="Host core budget (for throughput concurrency)")
    ap.add_argument("--subsample", type=int, default=15_000_000, help="Read pairs to use (0 = full input)")
    ap.add_argument("--threads-list", default="6,8,12,16,24", help="CSV of bwa -t values to test")
    ap.add_argument("--mode", choices=["single", "throughput", "both"], default="both")
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--pipe", choices=["full", "bwa-only"], default="full", help="full = bwa|sambamba sort (like the rule)")
    ap.add_argument("--fix-K", type=int, default=100_000_000, help="Fix bwa -K for cross-thread reproducibility (0 to disable)")
    ap.add_argument("--sort-mem", default="4GB")
    ap.add_argument("--scratch", default="/tmp/bwabench")
    ap.add_argument("--warmup", action="store_true", default=True)
    ap.add_argument("--no-warmup", dest="warmup", action="store_false")
    ap.add_argument("--out", default="benchmarks/results/bwa_bench.tsv")
    ap.add_argument("--dry-run", action="store_true", help="Print the run matrix + rough time and exit")
    args = ap.parse_args()

    threads_list = [int(t) for t in args.threads_list.split(",") if t.strip()]
    modes = ["single", "throughput"] if args.mode == "both" else [args.mode]
    runs = plan_runs(threads_list, args.cores, modes, args.reps)

    print(f"[bench] host cores={args.cores}  subsample={args.subsample:,} pairs  pipe={args.pipe}  fix_K={args.fix_K}")
    print(f"[bench] plan: {len(runs)} runs ({', '.join(modes)} x {len(threads_list)} threads x {args.reps} reps)")
    for m in modes:
        splits = ", ".join(f"{t}t/{max(1, args.cores // t) if m == 'throughput' else 1}job" for t in threads_list)
        print(f"[bench]   {m:10s}: {splits}")
    print("[bench] rough time ~= (runs / reps) x (one subsample-map wall). Each map is a few min on a subsample.")
    if args.dry_run:
        print("[bench] --dry-run: nothing executed.")
        return

    if not os.path.exists("/usr/bin/time"):
        print("[bench] note: /usr/bin/time not found; wall clock measured in-process (fine).")

    bwa_sif, samb_sif = load_sifs(args.containers)
    os.makedirs(args.scratch, exist_ok=True)
    sub1, sub2 = subsample(args.fq1, args.fq2, args.subsample, args.scratch)
    binds = bind_flags([args.ref, sub1, sub2, args.scratch])
    sub_pairs = args.subsample if args.subsample > 0 else _count_pairs(sub1)

    def make_builder(threads: int):
        return lambda out_bam: build_cmd(
            threads, sub1, sub2, args.ref, bwa_sif, samb_sif, out_bam,
            args.fix_K, args.pipe, args.sort_mem, binds,
        )

    if args.warmup:
        print("[bench] warm-up (load index into cache)...")
        try:
            run_config(make_builder(max(threads_list)), 1, args.scratch)
        except Exception as e:
            print(f"[bench] warm-up failed ({e}); continuing.")

    rows: list[dict[str, Any]] = []
    for i, r in enumerate(runs, 1):
        threads, conc = r["threads"], r["concurrency"]
        print(f"[bench] ({i}/{len(runs)}) mode={r['mode']} -t{threads} x{conc} rep{r['rep']} ...", flush=True)
        try:
            wall, rss = run_config(make_builder(threads), conc, args.scratch)
        except Exception as e:
            print(f"[bench]   FAILED: {e}")
            continue
        reads = conc * sub_pairs * _READS_PER_PAIR
        row = {
            **r,
            "wall_s": round(wall, 1),
            "reads": reads,
            "throughput_reads_s": round(reads / wall) if wall > 0 else 0,
            "throughput_per_core": round(reads / wall / args.cores) if wall > 0 else 0,
            "peak_rss_gb": round(rss, 1),
        }
        rows.append(row)
        print(f"[bench]   wall={wall:.0f}s  {row['throughput_reads_s']:,} reads/s  peak_rss={rss:.1f}GB")

    _write_tsv(rows, args.out)
    _report(rows, args.cores)


def _count_pairs(path: str) -> int:
    n = subprocess.run(["wc", "-l", path], capture_output=True, text=True).stdout.split()[0]
    return int(n) // 4


def _write_tsv(rows: list[dict[str, Any]], out: str) -> None:
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    cols = ["mode", "threads", "concurrency", "rep", "wall_s", "reads",
            "throughput_reads_s", "throughput_per_core", "peak_rss_gb"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"[bench] wrote {len(rows)} rows -> {out}")


def _report(rows: list[dict[str, Any]], cores: int) -> None:
    if not rows:
        print("[bench] no successful runs.")
        return
    s = summarize(rows)
    print("\n=== throughput by per-job threads (median reads/s) ===")
    for threads in sorted(s["throughput_by_threads"]):
        conc = max(1, cores // threads)
        print(f"  -t {threads:>3}  x{conc:<2} jobs : {s['throughput_by_threads'][threads]:>14,.0f} reads/s")
    best = s["best_threads"]
    if best is not None:
        try:
            from scripts.auto_params import derive_bwa_threads, get_total_memory_gb
            auto = derive_bwa_threads(cores, get_total_memory_gb())
            verdict = "matches auto" if auto == best else f"auto picked {auto} -> revisit"
            print(f"\n[bench] host-best bwa_threads = {best}  ({verdict})")
            print(f"[bench] to pin: set  mapping.bwa_threads: {best}  in your host config")
        except Exception:
            print(f"\n[bench] host-best bwa_threads = {best}")


if __name__ == "__main__":
    main()
