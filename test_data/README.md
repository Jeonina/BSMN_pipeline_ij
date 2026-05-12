# `test_data/` — manual smoke-test fixtures only

This directory contains fixtures intended for **manual** developer
smoke-tests. **No CI job or characterization test depends on the files
here.**

The canonical fixture path used by CI, unit tests, and golden-snapshot
generation is `tests/data/` (see SPEC-BSMN-REFACTOR-001 §1.5).

## Contents

| File | Status | Purpose |
|---|---|---|
| `TEST001_1.fastq.gz` | retained | manual smoke-test (paired-end R1) |
| `TEST001_2.fastq.gz` | retained | manual smoke-test (paired-end R2) |

## Removed in M1a (2026-05-12)

The zero-byte placeholders `SAMPLE001_R1.fastq.gz` and
`SAMPLE001_R2.fastq.gz` were deleted in M1a of SPEC-BSMN-REFACTOR-001;
they carried no data and misled anyone inspecting the tree.

## Why this directory exists separately from `tests/data/`

`test_data/TEST001_*.fastq.gz` predates the Snakemake migration and is
kept untouched to support manual reruns of legacy `jobs/*.py` wrappers
that some developers still invoke locally. After M7 retires the `jobs/`
directory, this fixture set may be archived or removed in a follow-up
SPEC.
