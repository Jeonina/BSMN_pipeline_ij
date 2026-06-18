# Production Validation Checklist / 실전 투입 전 검증 체크리스트

The pipeline's wiring, DAG, container manifest, and config keys are statically
verified by `tests/` + `ruff` + the workflow audit. Three things, however,
depend on **external-tool output contracts** that cannot be checked without real
data / the actual container image — and several of them **fail silently** (the
filter becomes a no-op rather than crashing). Confirm these once on the server
before the first full run.

정적 검증(테스트/ruff/감사)으로 잡히지 않는, 외부 도구의 출력 형식에 의존하는
3가지 항목입니다. 일부는 실패해도 에러 없이 "필터가 그냥 통과"가 되므로(silent
no-op) 첫 풀런 전에 실제 데이터로 1회 확인하세요.

Replace `{SAMPLE}` and SIF paths (from `config/containers.yaml`) accordingly.

---

## 0. Setup gate (must pass first)

```bash
bash scripts/server_setup.sh        # ends with the resource verify gate
```

Confirm the final `[5/5] Verifying resources` step lists every file with `✓`,
including `af-only-gnomad.hg38.vcf.gz`, `gnomAD.hg38.AFover0.001.snps.txt.gz`,
the `.bed` 1KG mask, and
`resources/MosaicForecast/models_trained/250xRFmodel_addRMSK_Refine.rds`.

Then dry-run the DAG:

```bash
snakemake --snakefile workflow/Snakefile --config stage=filtering -n
```

---

## 1. CNVnator `-genotype` output parsing  ⚠ fails OPEN

`scripts/cnvnator_filter.py` parses `cnvnator -genotype` stdout assuming the
copy-number is at a fixed column. If the column layout differs, the parser keeps
**every** candidate (the D-step silently becomes a no-op) instead of crashing.

**Check** — after a real filtering run, inspect the log:

```bash
grep -E "raw|kept|dropped|CN" logs/filtering/{SAMPLE}/cnvnator_filter.log | head
```

Expected: the 3 logged raw genotype lines parse to plausible CN values, and the
kept/dropped counts make sense (you should see *some* drops in a real WGS sample,
unless none of the candidates fall in CN≥2.5 regions).

**If suspicious**, reproduce one genotype call directly in the container and read
off the actual columns (CN is what `cnvnator_filter.py` must pick up):

```bash
echo "chr20:47053475-47053475" | \
  apptainer exec <cnvnator.sif> cnvnator -root results/calling/{SAMPLE}/cnvnator/{SAMPLE}.root -genotype 100
```

Confirm `cnvnator_filter.py`'s `cn_field` index matches the CN position in that
output. Source of truth for the format: CNVnator 0.4.1.

---

## 2. gnomAD lookup chrom encoding  ⚠ would fail OPEN

`germline_filter.py` strips the `chr` prefix from each candidate before lookup
(`is_germline`, line 49), so the gnomAD lookup table **must be bare-chrom**.

By construction this is consistent: `extract_hg38_gnomad_snps.py` strips `chr`
by default, and `download_and_index_hg38.sh` Step 6 calls it without
`--no-strip-chr`. Just confirm it:

```bash
zcat resources/hg38/gnomAD.hg38.AFover0.001.snps.txt.gz | head -1
# expect bare chrom, e.g.  "1   12345   A   G"   (NOT "chr1 ...")
```

**Runtime confirmation** (the real proof it isn't a no-op):

```bash
grep -E "variants_kept|variants_removed|pass_rate" logs/filtering/{SAMPLE}/germline_filter.log
# removed must be > 0 on a real sample
```

---

## 3. MosaicForecast image contract  ⚠ default E-step

`scripts/mosaicforecast_filter.py` hardcodes in-image script paths, the
`ReadLevel_Features_extraction.py` positional arg order, and the prediction
column indices (`PRED_COL=34`, `PROB_COL=36`, i.e. BSMN 1-based `$35`/`$37`).
These are pinned to `yanmei/mosaicforecast:0.0.1`; a previous commit (`b639496`)
already had to correct the RLF signature once.

**Check the image has the expected entry points:**

```bash
apptainer exec <mosaicforecast.sif> ls -la \
  /usr/local/bin/ReadLevel_Features_extraction.py \
  /usr/local/bin/Prediction.R \
  /usr/local/bin/k24.umap.wg.bw
```

**Smoke-test on 1–2 candidates** and confirm the prediction TSV has the label in
column 35 and the probability in column 37 (1-based), and that the filter emits
sensible output:

```bash
# run the filtering cascade on a tiny region (e.g. config.chr20.yaml), then:
grep -E "input|kept|prediction|mosaic" logs/filtering/{SAMPLE}/mosaicforecast_filter.log | head
awk '{print NF; exit}' results/filtering/{SAMPLE}/mf/*.pred*  # column count >= 37
```

If the column count or labels differ, update `PRED_COL`/`PROB_COL` and the RLF
arg order in `scripts/mosaicforecast_filter.py`.

---

## Summary

| # | Item | Failure mode | Quick check |
|---|------|--------------|-------------|
| 1 | CNVnator `-genotype` parse | silent no-op (keeps all) | cnvnator_filter.log: drops > 0 |
| 2 | gnomAD lookup chrom | silent no-op (germline leak) | lookup head = bare chrom; removed > 0 |
| 3 | MosaicForecast image contract | wrong/empty mosaic set | image scripts exist; pred TSV ≥ 37 cols |

All three are also the items flagged WARNING/CORRECTNESS in the pre-production
code audit. Once confirmed, the pipeline is cleared for the full cohort run.
