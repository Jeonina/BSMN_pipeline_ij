# BSMN Snakemake Mapping Pipeline - Test Report

## 개요

BSMN somatic variant calling pipeline의 mapping 모듈을 Snakemake 기반으로 리팩토링한 후,
실제 WES 데이터를 사용하여 전체 매핑 파이프라인의 정상 작동 여부를 검증하였다.

## 테스트 환경

| 항목 | 내용 |
|------|------|
| OS | Linux (WSL2) |
| Snakemake | conda env `snakemake` |
| Reference | hg38 (chr22 subset) |
| 입력 데이터 | WES paired-end FASTQ (1 sample, 1 readgroup) |
| Cores | 2 |

## 파이프라인 단계 및 결과

| 단계 | Rule | 상태 | 비고 |
|------|------|------|------|
| 1 | `bwa_mem_sort` | 완료 | BWA-MEM → sambamba view → sambamba sort |
| 2 | `merge_bams` | 완료 | Single RG → rename |
| 3 | `mark_duplicates` | 완료 | Picard MarkDuplicates |
| 4 | `base_recalibrator` | 완료 | GATK4 BaseRecalibrator |
| 5 | `apply_bqsr` | 완료 | GATK4 ApplyBQSR → CRAM 변환 |
| 6 | `samtools_flagstat` | 완료 | QC 통계 생성 |

## Flagstat QC 결과

| 항목 | 수치 |
|------|------|
| Total reads | 83,850,746 |
| Primary reads | 83,551,700 |
| Mapped reads | 10,902,654 (13.00%) |
| Properly paired | 7,375,692 (8.83%) |
| Duplicates | 1,914,587 |

> **Note**: Mapped 비율이 낮은 것은 chr22 레퍼런스만 사용했기 때문이며, 정상 동작 확인 목적으로는 문제 없음.

## 수정 사항 (버그 수정)

### `apply_bqsr` rule — `-O /dev/stdout` 옵션 추가

- **파일**: `workflow/rules/mapping.smk` (rule `apply_bqsr`)
- **증상**: GATK ApplyBQSR 출력이 samtools로 전달되지 않아 `fail to read the header from "-"` 에러 발생
- **원인**: ApplyBQSR의 기본 출력이 stdout이 아니므로, pipe로 연결된 samtools가 입력을 받지 못함
- **수정**: `-O /dev/stdout` 옵션을 명시적으로 추가

```diff
  gatk ... ApplyBQSR \
      -R {params.ref} \
      --bqsr-recal-file {input.table} \
      -I {input.bam} \
+     -O /dev/stdout \
  | samtools view -@ {threads} -C -T {params.ref} \
      -o {output.cram}
```

## 결론

Snakemake 기반 mapping 파이프라인이 전체 단계에서 정상 작동함을 확인하였다.
`apply_bqsr` rule의 버그 수정 1건을 반영하였으며, 전체 게놈 레퍼런스 적용 시
실제 분석에 사용할 수 있는 상태이다.

---

- **테스트 일자**: 2026-03-20
- **작성자**: Jeonina
