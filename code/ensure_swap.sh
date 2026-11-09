#!/bin/bash
# The analysis step's memory grows with sequencing depth (about 6 GB at ~190x). Keep a 12 GB swap file active.
[ -n "$(swapon --show --noheadings 2>/dev/null)" ] && exit 0
[ -f /swapfile ] || { fallocate -l 12G /swapfile && chmod 600 /swapfile && mkswap /swapfile >/dev/null; }
swapon /swapfile 2>/dev/null; sysctl -w vm.swappiness=10 >/dev/null 2>&1
echo "[swap] enabled at $(date -u +%H:%M:%S)" >> /tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/logs/queue.log
