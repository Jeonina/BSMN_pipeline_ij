# bwa thread-width benchmark

`bwa_threads` is the one auto-tuner knob with no published optimum
(`scripts/auto_params.py:derive_bwa_threads`, labelled `[HEURISTIC]`). This harness
measures the best per-job width **on a given host** so you can confirm the auto value
or pin an override.

## Run

Always preview first (`--dry-run` executes nothing):

```bash
python scripts/benchmark_bwa.py \
  --fq1 /data/.../E250....R1.fq.gz --fq2 /data/.../E250....R2.fq.gz \
  --ref resources/hg38/Homo_sapiens_assembly38.fasta \
  --containers config/containers.yaml \
  --cores $(nproc) --subsample 15000000 \
  --threads-list 6,8,12,16,24 --mode both --reps 2 \
  --scratch /tmp/bwabench \
  --out benchmarks/results/bwa_$(hostname).tsv \
  --dry-run          # <- drop this to actually run
```

Run the real sweep inside `tmux` (it is minutes-per-config, ~1.5–2.5 h total for the
default matrix; the subsample is taken once and reused).

## Modes

- `throughput` (the decision): at fixed `--cores`, run `cores // -t` jobs in parallel
  and vary `-t`. The split with the highest total reads/s wins → that is `bwa_threads`.
- `single` (diagnostic): one job, vary `-t` → the per-job scaling curve; shows where
  bwa's speed-up flattens on this host (validates the "~16 knee" claim on-metal).

## Reading the output

The TSV has one row per run; the console prints median reads/s per split and the
`host-best bwa_threads`, compared against what `derive_bwa_threads` would pick:

- **matches auto** → the heuristic is right for this host-class; nothing to change.
- **auto picked N → revisit** → pin `mapping.bwa_threads: <best>` in a host config,
  or, if several host-classes agree on a different value, adjust `_BWA_TARGET_THREADS`.

## Caveats

- **Per-host.** A 24-core result characterises the 24-core class only; a 168-core /
  RAM-bound host is a different regime — benchmark it separately when it is idle.
- **Scratch I/O confound.** If `--scratch` (and the subsample) sit on NFS/HDD, absolute
  throughput is depressed, but the *ranking* of splits still holds (the cost is equal
  across configs). Use a fast local `--scratch` when available for cleaner absolute numbers.
- `--fix-K 100000000` pins bwa's batch size so timing is comparable across thread counts
  (bwa couples batch size to `-t` by default). This differs from the pipeline default (no
  `-K`); it is a benchmarking control, not a pipeline change.
