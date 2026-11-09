#!/bin/bash
# Restart-proof batch runner. usage: resume.sh [--dry] MTB37 MTB38 ...
#  * decides what is DONE from disk (verified report + check record with FAILURES: 0), never from a log line
#  * runs only pending samples, strictly one at a time (a lock prevents two copies running)
#  * a sample that does not end verified is logged as FAILED, retried once, then reported honestly
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
IN=/mnt/user-data/uploads/Downloads/FASTQ_55_samples
L=$P/logs/queue.log
DRY=0; [ "${1:-}" = "--dry" ] && { DRY=1; shift; }
done_ok(){ S=$1
  [ -f "$P/out/$S/${S}_complete_analysis_report_2.pdf" ] && [ -f "$P/work/$S/results.json" ] && [ -f "$P/work/$S/.asm.ok" ] \
  && grep -q "FAILURES: 0" "$P/logs/$S.check" 2>/dev/null \
  && [ "$P/out/$S/${S}_complete_analysis_report_2.pdf" -nt "$P/work/$S/results.json" ] \
  && [ -f "/mnt/user-data/outputs/${S}_report2/${S}_complete_analysis_report_2.pdf" ]; }
if [ $DRY -eq 1 ]; then for S in "$@"; do done_ok $S && echo "$S DONE" || echo "$S pending"; done; exit 0; fi
bash $P/code/ensure_swap.sh
exec 9>$P/logs/.resume.lock; flock -n 9 || { echo "another batch is already running"; exit 1; }
bash $P/code/check_frozen.sh >> $L 2>&1 || { echo "[resume] ABORT: pipeline changed" >> $L; exit 1; }
echo "[resume] start: $* at $(date -u +%H:%M:%S)" >> $L
for S in "$@"; do
  done_ok $S && { echo "[resume] $S already verified, skipping" >> $L; continue; }
  for try in 1 2; do
    until [ "$(ls $IN/${S}_S*_R[12]*.gz 2>/dev/null | wc -l)" = 2 ] && [ -f "$IN/.verified_$S" ]; do sleep 30; done   # start only after both reads are md5-verified against the originals
    echo "[resume] $S attempt $try at $(date -u +%H:%M:%S)" >> $L
    bash $P/code/run_sample.sh $S
    if done_ok $S; then
      echo "[resume] $S VERIFIED $(grep -h 'checks passed' $P/logs/$S.check) at $(date -u +%H:%M:%S)" >> $L
      # free disk only AFTER verification: staged reads and alignment files are not used by any report
      rm -f "${IN:?}"/"${S:?}"_S*_R[12]*.gz "${IN:?}"/.verified_"${S:?}"
      rm -f "${P:?}"/work/"${S:?}"/"${S:?}".bam "${P:?}"/work/"${S:?}"/"${S:?}".bam.bai "${P:?}"/work/"${S:?}"/raw.bam "${P:?}"/work/"${S:?}"/raw.bam.bai "${P:?}"/work/"${S:?}"/t1.fq.gz "${P:?}"/work/"${S:?}"/t2.fq.gz
      break
    fi
    echo "[resume] $S FAILED attempt $try at $(date -u +%H:%M:%S)" >> $L
  done
  done_ok $S || echo "[resume] $S STILL FAILED after 2 attempts - needs attention" >> $L
done
echo "[resume] ALL DONE at $(date -u +%H:%M:%S)" >> $L
