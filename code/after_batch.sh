#!/bin/bash
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
while pgrep -f "code/batch.sh" >/dev/null; do sleep 20; done
echo "[batch] retrying MTB36 (analysis was OOM-killed by overlap with another job) at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
bash $P/code/run_sample.sh MTB36
echo "[batch] MTB36 retry finished $(date -u +%H:%M:%S) - $(grep -h 'checks passed' $P/logs/MTB36.check 2>/dev/null)" >> $P/logs/queue.log
