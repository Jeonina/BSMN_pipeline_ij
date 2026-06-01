# ERR194146 — WGS Production Run Statistics

First end-to-end production run of the refactored Snakemake BSMN somatic
variant-calling pipeline. SPEC-BSMN-REFACTOR-001.

## Run overview

| Item | Value |
|------|-------|
| Sample | ERR194146 (1000 Genomes HG00096, LCL) |
| Data | 190X WGS, hg38 (`Homo_sapiens_assembly38.fasta`) |
| Input | R1 + R2 = 813,180,578 reads each; ~60 GB × 2 compressed (~120 GB raw) |
| Environment | 80 core / 503 GB RAM workstation, Apptainer 1.4.5 |
| Run dates | 2026-05-18 ~ 05-19 (first end-to-end), final verification 05-21 |
| Result | 66 / 66 steps succeeded (100%), ~30 h total |

## Per-stage wall-clock

| # | Rule | Time |
|---|------|------|
| 1 | validate_fastq_pair | 25 min |
| 2 | bwa_mem_sort | 7 h (largest) |
| 3 | mark_duplicates | 2 h |
| 4 | base_recalibrator (24 chr) | 3 h |
| 5 | apply_bqsr | 4 h |
| 6 | Mutect2 (24 chr scatter) | 6 h |
| 7 | filter_mutect_calls | 0.5 h |
| 8 | 4 filters | 0.5 h |
| | **Total** | **~30 h** |

## Mapping QC (samtools flagstat)

| Metric | Value |
|--------|-------|
| Total reads | 1,628,690,072 |
| Mapped | 99.80% |
| Properly paired | 98.69% |
| Duplicates | 1.22% |
| Singletons | 0.14% |

Mapping, duplicate marking, and BQSR all within normal ranges.

## Output artifacts

| File | Size |
|------|------|
| `ERR194146.cram` (final) | 76 GB |
| `ERR194146.markduped.bam` (intermediate) | 112 GB |
| `ERR194146.filtered.vcf.gz` (Mutect2) | 312 MB |

## Filter cascade (post-fix, M-FIX-004 applied)

gnomAD lookup corrected from hg19 to hg38 coordinates (M-FIX-004) before
this cascade was captured.

```
Mutect2 raw      4,573,788
  → PASS            82,014
accessibility   47,508 → 21,695   (45.7% pass)
germline        kept = 21,437     (258 removed — the fix)
vaf             21,437 → 0        (0.0% pass)
                low_alt=1,020 · high_p=20,385 · no_cov=32
pon_mask        kept = 0
──────────────────────────────────────────────
Final mosaic                0
```

Upstream stages (Mutect2 raw → PASS) are unaffected by the gnomAD fix,
which only touches the germline stage.

### vaf_filter rejection breakdown (why all 21,437 dropped to 0)

| Reason | Count | Meaning |
|--------|-------|---------|
| high_p | 20,385 | VAF too high (germline-like, not mosaic) — Poisson/Fisher p-value above threshold |
| low_alt | 1,020 | alt allele support below minimum |
| no_cov | 32 | no coverage at site |
| **total rejected** | **21,437** | → final mosaic = 0 |

### Before / after M-FIX-004

| Stage | Pre-fix | Post-fix |
|-------|---------|----------|
| germline_filter removed | 26 (0.12%) | 258 (1.2%) |
| vaf_filter input | 21,669 | 21,437 |
| vaf_filter output | 0 | 0 |
| Final mosaic | 0 | 0 |

## Interpretation

The gnomAD coordinate fix raised direct germline removal 10x (26 → 258),
but the dominant filter remains `vaf_filter`, which rejected 20,385 of the
survivors as `high_p` — VAF at germline levels (~0.5 / 1.0) rather than
mosaic (< ~0.35). Both paths converge on **mosaic = 0**, the biologically
correct result for a clonally-expanded lymphoblastoid cell line (LCL): a
single progenitor clone leaves essentially no somatic mosaicism.

A spot check of 5 vaf-rejected variants (chr1:1,651,137–1,654,143, a 3 kb
window, all VAF = 1.0) confirmed they are absent from gnomAD and consistent
with a paralogous-gene-cluster mapping artifact.

Brain tissue (multiple cell lineages) is expected to yield 10–100 mosaic
candidates per sample; ERR194146 (LCL) validates the pipeline mechanics
without expecting positive calls.

## Source

Compiled from the 2026-05-21 production diagnostic run. Per-stage filter
counts were recovered from the run session record; headline figures also
appear in `.moai/specs/SPEC-BSMN-REFACTOR-001/progress.md` (M-FIX-004) and
`PRESENTATION_2026-05-21.md`.
