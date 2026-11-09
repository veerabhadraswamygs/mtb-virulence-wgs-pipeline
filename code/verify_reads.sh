#!/bin/bash
# usage: verify_reads.sh MTB65 ... -> a sample gets .verified_<S> only if BOTH reads
#   (1) match the original's md5, (2) pass gzip integrity, and (3) contain the same number of reads
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
IN=/mnt/user-data/uploads/Downloads/FASTQ_55_samples
for S in "$@"; do
  exp=$(grep " ${S}_S" $P/logs/batch16.md5)
  [ "$(echo "$exp" | wc -l)" = 2 ] || { echo "$S: no expected checksums"; continue; }
  ( cd $IN && echo "$exp" | md5sum -c --quiet - ) 2>/dev/null || { echo "$S NOT verified: checksum mismatch or file incomplete"; continue; }
  R1=$(ls $IN/${S}_S*_R1*.gz); R2=$(ls $IN/${S}_S*_R2*.gz)
  gzip -t "$R1" 2>/dev/null && gzip -t "$R2" 2>/dev/null || { echo "$S NOT verified: gzip integrity failed (truncated file)"; echo "[stage] $S failed gzip integrity at $(date -u +%H:%M:%S)" >> $P/logs/queue.log; continue; }
  n1=$(zcat "$R1" | awk 'END{print NR}'); n2=$(zcat "$R2" | awk 'END{print NR}')
  [ "$n1" = "$n2" ] || { echo "$S NOT verified: read counts differ (R1 $((n1/4)) vs R2 $((n2/4)))"; continue; }
  touch $IN/.verified_$S; echo "$S verified ($((n1/4)) read pairs)"
  echo "[stage] $S reads verified: md5 + gzip + $((n1/4)) pairs at $(date -u +%H:%M:%S)" >> $P/logs/queue.log
done
