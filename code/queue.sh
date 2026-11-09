#!/bin/bash
# Sequential queue: processes any sample whose reads are staged and which has no
# finished report yet. Safe to leave running; picks up newly staged samples.
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
IN=/mnt/user-data/uploads/Downloads/FASTQ_55_samples
ORDER="MTB3 MTB4 MTB9 MTB14 MTB18 MTB21 MTB22 MTB24 MTB25 MTB28 MTB29 MTB30 MTB31 MTB32 MTB35 MTB36 MTB37 MTB38 MTB39 MTB40 MTB41 MTB42 MTB43 MTB44 MTB45 MTB46 MTB48 MTB52 MTB53 MTB54 MTB56 MTB57 MTB58 MTB59 MTB60 MTB61 MTB62 MTB63 MTB64 MTB65 MTB66 MTB67 MTB68 MTB69 MTB71 MTB72 MTB77 MTB78 MTB79 MTB80 MTB81 MTB82 MTB83 MTB84 MTB131"
while true; do
  did=0
  for S in $ORDER; do
    [ -f "$P/out/$S/${S}_complete_analysis_report_2.pdf" ] && [ -f "$P/work/$S/.asm.ok" ] && continue
    ls $IN/${S}_S*_R1*.gz >/dev/null 2>&1 || continue
    echo "[queue] starting $S at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
    bash $P/code/run_sample.sh $S
    echo "[queue] finished $S at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
    did=1
  done
  [ $did -eq 0 ] && { echo "[queue] nothing to do, sleeping" >> $P/logs/queue.log; sleep 120; }
done
