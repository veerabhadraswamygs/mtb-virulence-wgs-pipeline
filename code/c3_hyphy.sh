#!/bin/bash
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2; C=$P/cohort/c3; g=$1; d=$C/hyphy/$g; mkdir -p $d; cd $d
[ -s fel.json ] && [ -s meme.json ] && exit 0
timeout 600 $P/bin/hyphy LIBPATH=$P/share/hyphy CPU=1 fel --alignment $C/aln_c/$g.fas --tree $C/aln_c/$g.nwk --pvalue 0.1 --output fel.json > fel.log 2>&1 </dev/null
timeout 900 $P/bin/hyphy LIBPATH=$P/share/hyphy CPU=1 meme --alignment $C/aln_c/$g.fas --tree $C/aln_c/$g.nwk --pvalue 0.1 --output meme.json > meme.log 2>&1 </dev/null
