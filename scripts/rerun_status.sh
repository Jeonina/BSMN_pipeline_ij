#!/usr/bin/env bash
# Per-sample progress for a run that goes mapping -> calling -> filtering (stage=all).
# Usage: bash scripts/rerun_status.sh [sample_sheet] [log_glob]
#   defaults: config/samples.rerun3.tsv  logs/rerun3_*.log
# Stage markers are checked newest-first; temp() files exist while the NEXT step runs.
for r in "$REPO" ~/IJ/BSMN_pipeline_ij ~/BSMN_pipeline_ij; do [ -n "$r" ] && [ -d "$r/config" ] && { cd "$r"; break; }; done
SHEET=${1:-config/samples.rerun3.tsv}
LOG=$(ls -t ${2:-logs/rerun3_*.log} 2>/dev/null | head -1)
# On NFS, stat() of a file another process is writing waits for that writer's
# dirty pages to flush (nfs_getattr) -- minutes on a 50 GB BAM. So never stat
# outputs here: test existence by reading the directory listing only.
ex(){ \ls -f "${1%/*}" 2>/dev/null | grep -qxF -- "${1##*/}"; }
cnt(){ \ls -f "${1%/*}" 2>/dev/null | grep -c -- "$2"; }
NCHR=24
PS=$(ps -eo comm,args)
echo "작업폴더 : $PWD"; echo "시트     : $SHEET"; echo "로그     : $LOG"

for s in $(tail -n +2 "$SHEET" | cut -f1 | sort -u); do
  m=results/mapping/$s c=results/calling/$s f=results/filtering/$s
  R=$(awk -F'\t' -v s="$s" '$1==s' "$SHEET" | wc -l | tr -d ' ')
  if   ex "$f/$s.final.txt";                    then st="F6 FINAL ($(wc -l < $f/$s.final.txt) SNVs)"
  elif ex "$f/$s.mosaicforecast_filtered.txt";  then st="F5 MF 완료 → PON"
  elif ex "$f/$s.cnvnator_filtered.txt";        then st="F4 CNV 필터 완료 → MosaicForecast"
  elif ex "$f/cnvnator/$s.root";                then st="F3 CNVnator ROOT 완료"
  elif ex "$f/$s.vaf_filtered.txt";             then st="F2 VAF 완료 (cnvnator_root 대기/진행)"
  elif { ex "$f/$s.accessible.txt" || ex "$f/$s.germline_filtered.txt"; }; then st="F1 필터 A/B"
  elif ex "$c/$s.filtered.vcf.gz.tbi";          then st="C4 콜링 완료 → 필터링 대기"
  elif ex "$c/$s.merged.vcf.gz";                then st="C3 gather (GetPileupSummaries/FilterMutectCalls)"
  elif ex "$c/scatter";                         then st="C2 Mutect2 $(cnt "$c/scatter/x" "\.vcf\.gz\.stats$")/$NCHR chrom"
  elif ex "$m/flagstat.txt";                    then st="M6 매핑 완료 → 콜링 대기"
  elif ex "$m/$s.cram";                         then st="M5 ApplyBQSR (CRAM 작성 중)"
  elif ex "$m/recal_data.table";                then st="M4 BQSR gather 완료 → ApplyBQSR"
  elif ex "$m/bqsr_scatter";                    then st="M4 BaseRecalibrator $(cnt "$m/bqsr_scatter/x" "\.recal\.table$")/$NCHR"
  elif ex "$m/$s.markduped.bam";                then st="M3 MarkDuplicates"
  elif ex "$m/$s.merged.bam";                   then st="M2 merge 완료/진행"
  elif ex "$m/align";                           then
    # Judge each readgroup by the REAL running binary (ps comm), not by files or
    # by the job's bash -c line (that line names both bwa and the output BAM):
    #   bwa reading this rg's fq1 -> aligning; sambamba writing its BAM -> merging.
    done_n=0 merge_n=0 bwa_n=0
    while IFS=$'\t' read -r sid rg fq1 _; do
      [ "$sid" = "$s" ] || continue
      b="$m/align/$s.$rg.sorted.bam"
      if   echo "$PS" | awk -v f="$fq1" '$1=="bwa" && index($0, f)' | grep -q .; then bwa_n=$((bwa_n+1))
      elif ex "$b" && echo "$PS" | awk -v f="$b" '$1=="sambamba" && index($0, f)' | grep -q .; then merge_n=$((merge_n+1))
      elif ex "$b"; then done_n=$((done_n+1)); fi
    done < <(tail -n +2 "$SHEET")
    st="M1 bwa 완료 $done_n/$R  (bwa 중 $bwa_n, 정렬 병합 중 $merge_n)"
  elif ex "$m/validation";                      then st="M0 FASTQ 검증 $(cnt "$m/validation/x" "\.ok$")/$R"
  else st="_QUEUED"; fi
  printf "  %-8s %s\n" "$s" "$st"
done

echo "진행     : $(grep -oE '[0-9]+ of [0-9]+ steps \([0-9]+%\)' "$LOG" | tail -1)"
echo "에러     : $(grep -c 'Error in rule' "$LOG")"
grep 'Error in rule' "$LOG" | sort | uniq -c | head -5
# Count top-level drivers only: Snakemake re-invokes itself (python -m snakemake
# --target-jobs ...) to run some rules, and those children are not extra drivers.
echo "드라이버 : $(ps -eo args | grep '[s]nakemake --snakefile' | grep -v -e '--target-jobs' -e '^/bin/sh -c' | wc -l)  (작업 실행용 하위 프로세스 제외)"
# Real JVMs only (wrappers also carry the tool name) -> count comm=="java".
ps -eo comm,args | awk '
  $1=="bwa"{b++} $1=="cnvnator"{cn++}
  $1=="java"{for(t in T) if (index($0,t)) n[t]++}
  BEGIN{split("MarkDuplicates BaseRecalibrator ApplyBQSR Mutect2 GetPileupSummaries LearnReadOrientationModel FilterMutectCalls",a," "); for(i in a) T[a[i]]=1}
  END{printf "실행 중  : bwa %d", b; for(t in n) printf "  %s %d", t, n[t]; printf "  cnvnator %d\n", cn}'
echo "가용메모리: $(free -g | awk '/Mem:/{print $7}') GB"
