#!/usr/bin/env bash
# Per-sample progress for a run that goes mapping -> calling -> filtering (stage=all).
# Usage: bash scripts/rerun_status.sh [sample_sheet] [log_glob]
#   defaults: config/samples.rerun3.tsv  logs/rerun3_*.log
# Stage markers are checked newest-first; temp() files exist while the NEXT step runs.
for r in "$REPO" ~/IJ/BSMN_pipeline_ij ~/BSMN_pipeline_ij; do [ -d "$r/config" ] && { cd "$r"; break; }; done
SHEET=${1:-config/samples.rerun3.tsv}
LOG=$(ls -t ${2:-logs/rerun3_*.log} 2>/dev/null | head -1)
NCHR=24
echo "작업폴더 : $PWD"; echo "시트     : $SHEET"; echo "로그     : $LOG"

for s in $(tail -n +2 "$SHEET" | cut -f1 | sort -u); do
  m=results/mapping/$s c=results/calling/$s f=results/filtering/$s
  R=$(awk -F'\t' -v s="$s" '$1==s' "$SHEET" | wc -l | tr -d ' ')
  cnt(){ ls $1 2>/dev/null | wc -l | tr -d " "; }
  if   [ -s "$f/$s.final.txt" ];                    then st="F6 FINAL ($(wc -l < $f/$s.final.txt) SNVs)"
  elif [ -e "$f/$s.mosaicforecast_filtered.txt" ];  then st="F5 MF 완료 → PON"
  elif [ -e "$f/$s.cnvnator_filtered.txt" ];        then st="F4 CNV 필터 완료 → MosaicForecast"
  elif [ -s "$f/cnvnator/$s.root" ];                then st="F3 CNVnator ROOT 완료"
  elif [ -e "$f/$s.vaf_filtered.txt" ];             then st="F2 VAF 완료 (cnvnator_root 대기/진행)"
  elif [ -e "$f/$s.accessible.txt" -o -e "$f/$s.germline_filtered.txt" ]; then st="F1 필터 A/B"
  elif [ -e "$c/$s.filtered.vcf.gz.tbi" ];          then st="C4 콜링 완료 → 필터링 대기"
  elif [ -e "$c/$s.merged.vcf.gz" ];                then st="C3 gather (GetPileupSummaries/FilterMutectCalls)"
  elif [ -d "$c/scatter" ];                         then st="C2 Mutect2 $(cnt "$c/scatter/*.vcf.gz.stats")/$NCHR chrom"
  elif [ -s "$m/flagstat.txt" ];                    then st="M6 매핑 완료 → 콜링 대기"
  elif [ -e "$m/$s.cram" ];                         then st="M5 ApplyBQSR (CRAM 작성 중)"
  elif [ -e "$m/recal_data.table" ];                then st="M4 BQSR gather 완료 → ApplyBQSR"
  elif [ -d "$m/bqsr_scatter" ];                    then st="M4 BaseRecalibrator $(cnt "$m/bqsr_scatter/*.recal.table")/$NCHR"
  elif [ -e "$m/$s.markduped.bam" ];                then st="M3 MarkDuplicates"
  elif [ -e "$m/$s.merged.bam" ];                   then st="M2 merge 완료/진행"
  elif [ -d "$m/align" ];                           then
    # A readgroup is done only when its bwa|sambamba job has exited: sambamba
    # writes the BAM and its .bai while still merging, so files alone mislead.
    done_n=0 merge_n=0 bwa_n=0
    for b in "$m"/align/*.sorted.bam; do [ -e "$b" ] || continue
      if   ps -eo args | grep -F -- "$b" | grep -q '[b]wa mem'; then bwa_n=$((bwa_n+1))
      elif ps -eo args | grep -F -- "$b" | grep -q '[s]ambamba sort'; then merge_n=$((merge_n+1))
      else done_n=$((done_n+1)); fi
    done
    st="M1 bwa 완료 $done_n/$R readgroup (정렬 병합 중 $merge_n)"
  elif [ -d "$m/validation" ];                      then st="M0 FASTQ 검증 $(cnt "$m/validation/*.ok")/$R"
  else st="_QUEUED"; fi
  printf "  %-8s %s\n" "$s" "$st"
done

echo "진행     : $(grep -oE '[0-9]+ of [0-9]+ steps \([0-9]+%\)' "$LOG" | tail -1)"
echo "에러     : $(grep -c 'Error in rule' "$LOG")"
grep 'Error in rule' "$LOG" | sort | uniq -c | head -5
echo "드라이버 : $(pgrep -fc 'snakemake --snakefile')"
# Real JVMs only (wrappers also carry the tool name) -> count comm=="java".
ps -eo comm,args | awk '
  $1=="bwa"{b++} $1=="cnvnator"{cn++}
  $1=="java"{for(t in T) if (index($0,t)) n[t]++}
  BEGIN{split("MarkDuplicates BaseRecalibrator ApplyBQSR Mutect2 GetPileupSummaries LearnReadOrientationModel FilterMutectCalls",a," "); for(i in a) T[a[i]]=1}
  END{printf "실행 중  : bwa %d", b; for(t in n) printf "  %s %d", t, n[t]; printf "  cnvnator %d\n", cn}'
echo "가용메모리: $(free -g | awk '/Mem:/{print $7}') GB"
