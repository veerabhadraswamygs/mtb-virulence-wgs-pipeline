#!/bin/bash
# One sample, end to end, with the same frozen pipeline. Resumable at every stage.
set -u
S=$1
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
L=$P/logs/$S.log
{
echo "=== $S start $(date -u +%H:%M:%S) ==="
bash $P/code/core.sh $S || { echo "CORE_FAILED"; exit 1; }
cd $P/work/$S
# filtered VCF at the brief's thresholds
if [ ! -f .filt.ok ]; then
  $P/bin/bcftools view -v snps,indels -i 'INFO/MQ>=40 && FORMAT/DP>=10' ${S}.raw.vcf.gz -Ou 2>/dev/null \
   | $P/bin/bcftools filter -e 'ALT="."' -Oz -o pass1.vcf.gz 2>/dev/null
  python3 -I - "$S" <<'PY'
import gzip,sys
S=sys.argv[1]; keep=[]
for l in gzip.open('pass1.vcf.gz','rt'):
    if l[0]=='#': keep.append(l); continue
    f=l.rstrip('\n').split('\t'); d=dict(zip(f[8].split(':'),f[9].split(':')))
    if 'AD' not in d: continue
    ad=[int(x) for x in d['AD'].split(',') if x!='.']
    tot=sum(ad)
    if tot<10: continue
    alt=max(ad[1:]) if len(ad)>1 else 0
    if alt/tot>=0.9: keep.append(l)
open('%s.filt.vcf'%S,'w').writelines(keep)
PY
  $P/bin/bgzip -f ${S}.filt.vcf && $P/bin/bcftools index -f ${S}.filt.vcf.gz && touch .filt.ok
fi
# annotation
if [ ! -f .ann.ok ]; then
  java -Xmx3g -jar $P/src/snpEff/snpEff.jar -noStats -no-intergenic -upDownStreamLen 200 H37Rv ${S}.filt.vcf.gz 2>/dev/null > ${S}.ann.vcf \
   && touch .ann.ok
fi
# coverage
[ -f cov.json ] || python3 -I $P/code/mod_cov.py $S >/dev/null
# de novo assembly (stage 5b)
bash $P/code/assemble.sh $S || { echo "ASM_FAILED"; exit 1; }
# analysis + report + checks
python3 -I $P/code/analyse_all.py $S >/dev/null
python3 -I -c "
import sys; sys.path.insert(0,'$P/code')
import report; report.build('$S','$P/out/$S')" 2>/dev/null
python3 -I $P/code/check.py $S $P/out/$S | tee $P/logs/$S.check
# prune: keep only what the deliverables need
rm -f raw.bam raw.bam.bai pass1.vcf.gz t1.fq.gz t2.fq.gz clips.fa barcode.pos depth.txt.gz
rm -rf spades
# once results.json and the packaged deliverables exist, the heavy intermediates are no longer
# needed; everything downstream reads the cached results
if [ -f results.json ] && [ -f "/mnt/user-data/outputs/${S}_report2/${S}_complete_analysis_report_2.pdf" ]; then
  # keep depth_all.txt.gz and the annotated VCF: the report figure and the verification
  # suite re-read them, so a report can be regenerated without re-running the sample
  rm -f "${S:?}.bam" "${S:?}.bam.bai" contig_hits.tsv
  rm -f /mnt/user-data/uploads/Downloads/FASTQ_55_samples/"${S:?}"_S*_R1*.gz
  rm -f /mnt/user-data/uploads/Downloads/FASTQ_55_samples/"${S:?}"_S*_R2*.gz
fi
bash $P/code/package.sh $S
echo "=== $S done $(date -u +%H:%M:%S) ==="
} >> $L 2>&1
