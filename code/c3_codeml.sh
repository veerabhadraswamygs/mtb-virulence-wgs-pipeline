#!/bin/bash
# usage: c3_codeml.sh Rvxxxx  -> runs M0 (ω), M1a, M2a in one codeml run
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2; C=$P/cohort/c3
g=$1; NS=${2:-"0 1 2"}; d=$C/cml/$g; mkdir -p $d; cd $d
[ -s $d/out.txt ] && grep -q "lnL(ntime" out.txt && exit 0
cat > ctl <<CTL
seqfile = $C/aln_c/$g.phy
treefile = $C/aln_c/$g.nwk
outfile = out.txt
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
model = 0
NSsites = $NS
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
getSE = 0
RateAncestor = 0
Small_Diff = 1e-6
cleandata = 0
fix_blength = 0
method = 1
CTL
timeout 900 $P/bin/codeml ctl > log.txt 2>&1 </dev/null
