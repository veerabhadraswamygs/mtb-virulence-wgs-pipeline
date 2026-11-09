#!/bin/bash
# Rebuild report + checks + deliverables for already-analysed samples under the CURRENT frozen pipeline.
# Uses cached results.json (no alignment, no assembly). Low priority so it never slows a running batch.
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
bash $P/code/check_frozen.sh >> $P/logs/queue.log 2>&1 || { echo "[regen] ABORT: pipeline changed" >> $P/logs/queue.log; exit 1; }
for S in "$@"; do
  echo "[regen] $S start $(date -u +%H:%M:%S)" >> $P/logs/queue.log
  nice -n 10 python3 -I -c "
import sys; sys.path.insert(0,'$P/code')
import report; report.build('$S','$P/out/$S')" >> $P/logs/${S}.regen.log 2>&1
  nice -n 10 python3 -I $P/code/check.py $S $P/out/$S > $P/logs/$S.check.new 2>&1
  if grep -q "FAILURES: 0" $P/logs/$S.check.new; then
    mv $P/logs/$S.check.new $P/logs/$S.check
    bash $P/code/package_inplace.sh $S >> $P/logs/${S}.regen.log 2>&1
    echo "[regen] $S OK $(grep 'checks passed' $P/logs/$S.check) at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
  else
    echo "[regen] $S CHECK FAILED - old deliverables left untouched: $(tail -3 $P/logs/$S.check.new | tr '\n' ' ')" >> $P/logs/queue.log
  fi
done
echo "[regen] ALL DONE $(date -u +%H:%M:%S)" >> $P/logs/queue.log
