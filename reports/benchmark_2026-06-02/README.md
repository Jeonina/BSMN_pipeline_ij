# Benchmark validation — 2026-06-02

Pipeline: refactored Snakemake BSMN somatic mosaic caller (commit ac0c340).
Reference: GRCh38 (Homo_sapiens_assembly38). Mode: tumor-only Mutect2 + 4-filter cascade.

## VAF-filter bug found & fixed (ac0c340)
`vaf_filter.py` computed pileup depth from mismatch bases only, excluding
reference-matching reads ('.'/','). Every variant scored as VAF~1.0 → one-sided
binomial test p~1.0 → all candidates rejected (final mosaic = 0) on BOTH samples.
Fix: include '.'/',' in depth. Regression test added.

## Filter cascade comparison

| Stage | ERR194146 pre-fix (buggy) | ERR194146 post-fix | HG002 chr20 post-fix |
|-------|---------------------------|--------------------|----------------------|
| accessibility | 47,508 -> 21,695 | -> 21,695 | 3,564 -> 666 |
| germline kept | 21,437 | 21,418 | 638 |
| vaf pass | 21,437 -> 0 | 21,418 -> 300 | 638 -> 31 |
| final mosaic | 0 | 291 | 30 |
| nature | (bug artifact) | LCL whole-genome FP baseline | known mosaic, sensitivity 1/1 |

## HG002 chr20 smoke test
- Truth mosaic: chr20:47,053,475 G>T (the only chr20 variant among 85 genome-wide).
- Detected at calling (FILTER=PASS, VAF 0.115, 31/279 reads) AND survives filtering.
- Sensitivity = 1/1. Final candidates on chr20 = 30 (1 truth + 29 others).

## Notes / caveats
- chr20 has only 1 truth mosaic -> functional smoke test, not a quantitative benchmark.
- HG002 = external NHGRI 300X novoalign BAM (calling+filtering only, no pipeline mapping).
- ERR194146 (HG00096 LCL) is a near-clonal germline sample; the 291 post-fix
  candidates are largely the genome-wide false-positive baseline (specificity).
- Prior ERR194146 "mosaic = 0, biologically correct" conclusion was a bug artifact
  (see reports/ERR194146_production_run_stats.md — needs correction).
- Current thresholds: Mutect2 PON OFF; vaf alt>=5, binom p<1e-6; gnomAD AF>0.001;
  FilterMutectCalls --min-reads-per-strand 1. FP-reduction levers tracked separately.
