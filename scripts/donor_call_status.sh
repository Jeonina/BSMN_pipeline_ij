#!/usr/bin/env bash
# Progress of donor-merged calling on this server (VM100 or VM102).
# Both servers write into the same NFS results/ tree, so every count is restricted
# to this server's donor sheet. Mutect2 progress counts scatter .tbi files (written
# only when a chromosome finishes; the .vcf.gz grows while Mutect2 runs) plus 24 per
# donor whose merged filtered VCF already exists (scatter outputs are temp).
#
# Usage: bash ~/donor_call_status.sh            (auto-detects VM100 / VM102)
set -uo pipefail

if [ -d ~/IJ/BSMN_pipeline_ij ]; then
  host=VM100; cd ~/IJ/BSMN_pipeline_ij; sheet=config/samples.donor_merge.vm100.tsv; logs='logs/donor_call_[0-9]*.log'
else
  host=VM102; cd ~/BSMN_pipeline_ij; sheet=config/samples.donor_merge.vm102.tsv; logs='logs/donor_call_vm102_*.log'
fi
donors=$(tail -n +2 "$sheet" | cut -f1)
n_donor=$(echo "$donors" | wc -l)
total=$((n_donor * 24))
log=$(ls -t $logs 2>/dev/null | head -1)

m2_done=0; m2_day=0; called=0; final=0
for d in $donors; do
  if [ -s results/calling/$d/$d.filtered.vcf.gz ]; then
    m2_done=$((m2_done + 24)); called=$((called + 1))
  else
    m2_done=$((m2_done + $(ls results/calling/$d/scatter/*.vcf.gz.tbi 2>/dev/null | wc -l)))
    m2_day=$((m2_day + $(find results/calling/$d/scatter -name '*.vcf.gz.tbi' -mmin -1440 2>/dev/null | wc -l)))
  fi
  [ -s results/filtering/$d/$d.final.txt ] && final=$((final + 1))
done

echo "== $host  $(date '+%m/%d %H:%M')  donors $n_donor  log ${log:-none}"
pgrep -f '\.venv/bin/snakemake --cores' >/dev/null && echo "driver   : running" || echo "driver   : NOT RUNNING"
free -g | awk 'NR==2 {print "memory   : available " $7 " / " $2 " GB"} NR==3 {print "swap     : used " $3 " / " $2 " GB"}'
echo "Mutect2 java processes: $(ps -eo args | grep -c '^java .*Mutect2')"
if [ -n "$log" ]; then   # jobs started in this log and not yet finished, by rule
  python3 - "$log" <<'PY'
import collections, re, sys
started, finished, rule = {}, set(), None
for line in open(sys.argv[1]):
    m = re.match(r"\s*(?:local)?rule (\w+):", line)
    if m: rule = m.group(1)
    m = re.match(r"\s+jobid: (\d+)", line)
    if m and rule: started[m.group(1)] = rule; rule = None
    m = re.search(r"Finished job(?:id)?:? (\d+)", line)
    if m: finished.add(m.group(1))
run = collections.Counter(r for j, r in started.items() if j not in finished)
print("running  :", ", ".join(f"{r} {n}" for r, n in run.most_common()) or "none")
print("finished :", len(finished), "jobs in this log")
PY
fi
echo "Mutect2  : $m2_done / $total chromosomes  ($((100 * m2_done / total))%), last 24 h: $m2_day"
echo "donors   : calling done $called / $n_donor, final SNVs $final / $n_donor"
if [ -n "$log" ]; then
  echo "errors   : $(grep -c 'Error in rule' "$log") in current log"
  grep -oP 'Error in rule \K\w+' "$log" | sort | uniq -c | sed 's/^/           /'
fi
if [ "$m2_day" -gt 0 ]; then
  echo "estimate : ~$(( (total - m2_done + m2_day - 1) / m2_day )) more days of Mutect2 at the last-24h rate (rough)"
fi
