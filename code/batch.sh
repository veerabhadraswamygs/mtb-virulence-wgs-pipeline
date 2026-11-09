#!/bin/bash
# Run an explicit list of samples, in order, then EXIT. Never rolls on to other samples.
# usage: batch.sh MTB29 MTB30 ...
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
IN=/mnt/user-data/uploads/Downloads/FASTQ_55_samples
bash $P/code/check_frozen.sh >> $P/logs/queue.log 2>&1 || { echo "[batch] ABORT: pipeline changed" >> $P/logs/queue.log; exit 1; }
echo "[batch] start $# samples: $* at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
for S in "$@"; do
  # wait (max 3h) for both read files to be staged and complete
  w=0
  until [ "$(ls $IN/${S}_S*_R[12]*.gz 2>/dev/null | wc -l)" = 2 ] && [ ! -e "$IN/.staging_$S" ]; do
    sleep 30; w=$((w+30)); [ $w -gt 10800 ] && { echo "[batch] $S reads never arrived - skipping" >> $P/logs/queue.log; continue 2; }
  done
  echo "[batch] starting $S at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
  bash $P/code/run_sample.sh $S
  echo "[batch] finished $S at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
done
echo "[batch] ALL DONE - exiting at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
