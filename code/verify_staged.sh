#!/bin/bash
# verify_staged.sh SAMPLE : compare staged reads to the md5 taken on the user's computer; clear the staging marker only if both match
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
IN=/mnt/user-data/uploads/Downloads/FASTQ_55_samples
S=$1; ok=0
for f in $(ls "${IN:?}"/"${S:?}"_S*_R[12]*.gz 2>/dev/null); do
  b=$(basename "$f"); want=$(grep "  $b\$" "$P"/logs/device_md5.txt | tail -1 | cut -d' ' -f1); got=$(md5sum < "$f" | cut -d' ' -f1)
  [ -n "$want" ] && [ "$want" = "$got" ] && ok=$((ok+1)) || echo "MISMATCH $b want=$want got=$got"
done
if [ $ok -eq 2 ]; then rm -f "${IN:?}"/.staging_"${S:?}"; echo "[stage] $S reads verified (md5 match) at $(date -u +%H:%M:%S)" >> "$P"/logs/queue.log; echo "$S verified"; else echo "$S NOT verified ($ok/2)"; fi
