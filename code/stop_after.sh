#!/bin/bash
# Halt the sample queue once TARGET has a FRESH report (built after this watcher armed).
# Runs as its own process: editing the live queue.sh while bash is still reading it
# is what broke the queue earlier today.
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
T=${1:-MTB25}
R="$P/out/$T/${T}_complete_analysis_report_2.pdf"
A="$P/work/$T/.asm.ok"
B="$P/logs/.stop_baseline"
echo "[stop] armed for $T at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
while true; do
  if [ -f "$R" ] && [ -f "$A" ] && [ "$R" -nt "$A" ] && [ "$R" -nt "$B" ]; then
      sleep 25
      pkill -f "code/queue.sh"      >/dev/null 2>&1
      pkill -f "code/run_sample.sh" >/dev/null 2>&1
      pkill -f "spades"             >/dev/null 2>&1
      echo "[stop] $T finished - queue halted at $(date -u +%H:%M:%S) as requested" >> $P/logs/queue.log
      exit 0
  fi
  sleep 60
done
