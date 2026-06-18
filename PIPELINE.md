# BSMN Somatic Variant Calling Pipeline — Agent Instructions

## 프로젝트 목적

200명 bulk brain WGS (200X, hg38) 데이터에서 **brain somatic mutation (초돌연변이)** 을 발굴하는 파이프라인이다.
matched normal 없는 **tumor-only** 모드로 동작한다.

핵심 목표는 **데이터 경로만 입력하면 파라미터 결정부터 최종 결과까지 자동으로 돌아가는 one-command 플랫폼**을 만드는 것이다.
사용자는 샘플 이름과 FASTQ 경로만 입력하면 된다. 시퀀서 종류, 스레드 수, 메모리 등 모든 파라미터는 파이프라인이 데이터를 분석해서 자동으로 결정한다.
논문 제출 후 reviewer가 재현성을 요구할 때 자동 생성된 `resolved_params.yaml`을 공유하면 동일한 결과를 낼 수 있어야 한다.

---

## 전체 아키텍처 (지향점)

```
사용자 입력: 샘플 이름 + FASTQ 경로만 (samples.tsv)
        ↓
auto_params.py: FASTQ 분석 → 파라미터 자동 결정 → resolved_params.yaml 저장
        ↓
prepare_containers.py: 필요한 SIF 없으면 자동 pull
        ↓
Snakemake (execution engine)
        ↓
Apptainer (tool isolation)
        ↓
결과 + run manifest (재현성 기록)
```

- **samples.tsv**: 사용자가 입력하는 유일한 것. 샘플 이름 + FASTQ 경로만 기입.
- **auto_params.py**: FASTQ를 분석해서 시퀀서 종류, 스레드 수, 메모리 등 모든 파라미터를 자동으로 결정.
- **resolved_params.yaml**: 자동 결정된 파라미터 저장. 재현성 보장용.
- **Snakemake**: 파이프라인의 핵심 execution engine. DAG 관리, 병렬 실행, SLURM 연동.
- **Apptainer**: 각 툴을 container로 격리 실행. conda 환경 대신 사용한다.
- **run manifest**: 실행에 사용된 툴 버전, 파라미터, 입출력 경로를 자동으로 기록 (추후 추가 예정).

---

## 현재 구현 단계 (Phase 1)

**지금은 매핑 파이프라인을 Apptainer + TDD 방식으로 재구축하는 것이 목표다.**

### 지금 만들어야 할 것
1. `config/containers.yaml` — 툴별 container 정보
2. `scripts/auto_params.py` — 파라미터 자동 최적화
3. `tests/test_mapping.py` — TDD 테스트 케이스
4. `workflow/rules/mapping.smk` — 매핑 rule (Apptainer 기반)

### 나중에 붙일 것 (Phase 2~3)
- `workflow/rules/calling.smk` — Mutect2 tumor-only
- `workflow/rules/filtering.smk` — gnomAD, VAF, PON mask
- run manifest 자동 생성
- CLI wrapper (`bio-map`, `bio-call`, `bio-filter`)
- SLURM profile

---

## 파이프라인 단계별 상세

### Step 1: Mapping (현재 구현 중)
BSMN 원본 aln_1~5 기반. aln_4 (IndelRealigner)는 GATK4에서 불필요하므로 제외.

```
FASTQ (per readgroup)
  → bwa_mem_sort      : BWA-MEM + sambamba sort (per RG)
  → merge_bams        : sambamba merge (per sample)
  → mark_duplicates   : Picard MarkDuplicates (OPTICAL_DUPLICATE_PIXEL_DISTANCE 자동 결정)
  → base_recalibrator : GATK4 BaseRecalibrator
  → apply_bqsr        : GATK4 ApplyBQSR → CRAM
  → samtools_flagstat : QC
```

### Step 2: Calling (추후 구현)
- Mutect2 tumor-only mode (matched normal 없음)
- `--panel-of-normals`: BSMN PON (별도 파일)
- `--germline-resource`: gnomAD
- chromosome-level scatter → merge → FilterMutectCalls

### Step 3: Filtering (구현 완료)
BSMN 원본 필터링 로직 이식. 현재 기본 cascade는 6단계이며, 정본은
`workflow/rules/filtering.smk` (아래는 요약):

1. **Accessibility filter**: 1KG strict mask (mappability)
2. **Germline filter**: gnomAD AF > 0.001 제거
3. **VAF filter**: binomial test (p < 1e-6) AND alt_count >= 5 (BAM 직접 pileup 방식 유지)
4. **CNVnator filter**: BSMN D-step, CN >= 2.5 영역 제거
5. **MosaicForecast filter**: BSMN E-step, 학습된 RF mosaic 예측 (mayo는 미연결 대안)
6. **PON mask**: IUPAC FASTA 방식 (GATK Mutect2 `--panel-of-normals`와 별개로 추가 적용)

---

## Apptainer 사용 규칙

### container 정보 관리
모든 container 정보는 `config/containers.yaml`에서 관리한다:

```yaml
bwa:
  name: bwa
  version: v0.7.17
  uri: docker://biocontainers/bwa:v0.7.17
  sif: containers/bwa_0.7.17.sif

sambamba:
  name: sambamba
  version: v0.8.2
  uri: docker://quay.io/biocontainers/sambamba:0.8.2--h98b6b92_2
  sif: containers/sambamba_0.8.2.sif

picard:
  name: picard
  version: v3.1.1
  uri: docker://broadinstitute/picard:3.1.1
  sif: containers/picard_3.1.1.sif

gatk:
  name: gatk
  version: v4.4.0.0
  uri: docker://broadinstitute/gatk:4.4.0.0
  sif: containers/gatk_4.4.0.0.sif

samtools:
  name: samtools
  version: v1.17
  uri: docker://biocontainers/samtools:v1.17
  sif: containers/samtools_1.17.sif
```

### container 자동 pull
SIF 파일이 없으면 실행 전에 자동으로 pull한다:

```python
# scripts/prepare_containers.py
import subprocess, os, yaml

def prepare_containers(containers_yaml):
    with open(containers_yaml) as f:
        containers = yaml.safe_load(f)
    for tool, spec in containers.items():
        sif = spec['sif']
        if not os.path.exists(sif):
            os.makedirs(os.path.dirname(sif), exist_ok=True)
            subprocess.run(['apptainer', 'pull', sif, spec['uri']], check=True)
```

### Snakemake rule에서 호출
```python
rule bwa_mem_sort:
    params:
        sif = config["containers"]["bwa"]["sif"]
    shell:
        "apptainer exec {params.sif} bwa mem ..."
```

---

## 파라미터 자동 최적화 규칙

사용자는 데이터 경로만 입력한다. 아래 파라미터는 파이프라인이 FASTQ를 분석해서 자동으로 결정한다.

### 자동 결정 항목
| 파라미터 | 결정 방법 | 저장 위치 |
|---|---|---|
| BWA threads | 가용 CPU 수 (`os.cpu_count()`) | config/resolved_params.yaml |
| OPTICAL_DUPLICATE_PIXEL_DISTANCE | read name 패턴으로 기기 판별 | config/resolved_params.yaml |
| BQSR memory | 가용 메모리 기반 (`psutil`) | config/resolved_params.yaml |

### 기기 판별 기준
- read name에 `:` 구분자가 7개 이상 → NovaSeq → ODPD=2500
- 그 외 → HiSeq → ODPD=100

### resolved_params.yaml 예시
```yaml
resolved_at: "2026-03-18T12:00:00"
input_fastq: "/data/sample.R1.fastq.gz"
sequencer: "NovaSeq"
bwa_threads: 16
optical_duplicate_pixel_distance: 2500
bqsr_memory_gb: 32
```

---

## TDD 개발 규칙

### 테스트 구조
```
tests/
  test_mapping.py     # pytest 기반 전체 테스트
  data/
    tiny.R1.fastq.gz  # 소형 테스트 데이터 (직접 생성)
    tiny.R2.fastq.gz
```

### 각 rule 테스트 항목
1. output 파일이 생성됐는가
2. output 파일이 비어있지 않은가
3. 파일 포맷이 올바른가 (BAM header 확인, CRAM 등)
4. 예상 wildcard가 올바르게 치환됐는가

### 개발 순서
테스트 작성 → 실패 확인 → rule 구현 → 테스트 통과 → 다음 rule

---

## 절대 지켜야 할 규칙

- **사용자 입력은 samples.tsv 하나뿐**: 샘플 이름 + FASTQ 경로만. 파라미터는 자동 결정.
- **경로 하드코딩 금지**: 모든 경로는 config에서 읽는다
- **파라미터 하드코딩 금지**: 모든 파라미터는 auto_params.py가 결정하고 resolved_params.yaml에 저장
- **conda 사용 금지**: 모든 툴은 Apptainer container로 실행한다
- **중간 파일**: `temp()`로 마킹해서 자동 삭제
- **로그**: 모든 rule에 `log:` 섹션 포함
- **IndelRealigner 사용 금지**: GATK4에서 불필요
- **VQSR 사용 금지**: Mutect2 FilterMutectCalls로 대체

---

## 참고 파이프라인

| 파이프라인 | 역할 |
|---|---|
| `jobs/genome_mapping/` | 원본 BSMN 쉘 스크립트 — 로직 참고 |
| `tjbencomo/ngs-pipeline` | Snakemake 구조 참고 |
| `brian-arnold/SomaticVarCall` | Mutect2 rule 참고 |
