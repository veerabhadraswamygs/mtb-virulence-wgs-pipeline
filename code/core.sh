#!/bin/bash
# Pipeline 2 - core read-level stages. Resumable: each stage skipped if its .ok flag exists.
set -u
S=$1                       # sample id
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
export PATH=$P/bin:$PATH
IN=/mnt/user-data/uploads/Downloads/FASTQ_55_samples
W=$P/work/$S; mkdir -p $W; cd $W
R1=$(ls $IN/${S}_S*_R1*.gz 2>/dev/null | head -1)
R2=$(ls $IN/${S}_S*_R2*.gz 2>/dev/null | head -1)
[ -z "$R1" ] && { echo "FATAL no R1 for $S"; exit 1; }
echo "R1=$(basename "$R1")  R2=$(basename "$R2")"
ts(){ date +%s; }
T0=$(ts)
# --- S02 fastp ---
if [ ! -f .fastp.ok ]; then
  fastp -i "$R1" -I "$R2" -o t1.fq.gz -O t2.fq.gz --qualified_quality_phred 20 \
        --detect_adapter_for_pe --thread 2 -j ${S}_fastp.json -h /dev/null >/dev/null 2>fastp.err && touch .fastp.ok
fi
echo "fastp done $(( $(ts)-T0 ))s"
# --- S03 align ---
T1=$(ts)
if [ ! -f .bam.ok ]; then
  bwa mem -t 2 -R "@RG\tID:$S\tSM:$S\tPL:ILLUMINA" $P/ref/H37Rv.fa t1.fq.gz t2.fq.gz 2>bwa.err \
   | samtools sort -@1 -m 1G -o raw.bam - 2>>bwa.err && samtools index raw.bam && touch .bam.ok
fi
echo "align done $(( $(ts)-T1 ))s"
# --- S03b mark duplicates ---
T2=$(ts)
if [ ! -f .dedup.ok ]; then
  samtools collate -@1 -O -u raw.bam 2>/dev/null | samtools fixmate -m -u - - 2>/dev/null \
   | samtools sort -@1 -m 1G -u - 2>/dev/null | samtools markdup -@1 - ${S}.bam 2>markdup.err \
   && samtools index ${S}.bam && touch .dedup.ok
fi
echo "dedup done $(( $(ts)-T2 ))s"
# --- S04 depth ---
T3=$(ts)
if [ ! -f .depth.ok ]; then
  samtools depth -a -Q 20 -q 20 ${S}.bam | gzip -1 > depth.txt.gz
  samtools depth -a            ${S}.bam | gzip -1 > depth_all.txt.gz && touch .depth.ok
fi
echo "depth done $(( $(ts)-T3 ))s"
# --- S06 variant calling ---
T4=$(ts)
if [ ! -f .vcf.ok ]; then
  bcftools mpileup -f $P/ref/H37Rv.fa -a AD,DP,ADF,ADR -q 20 -Q 20 -d 1000 -Ou ${S}.bam 2>mpileup.err \
   | bcftools call --ploidy 1 -m -Ou 2>>mpileup.err \
   | bcftools norm -f $P/ref/H37Rv.fa -Oz -o ${S}.raw.vcf.gz 2>>mpileup.err \
   && bcftools index ${S}.raw.vcf.gz && touch .vcf.ok
fi
echo "call done $(( $(ts)-T4 ))s"
# --- S07 filter to the brief's thresholds: depth>=10, MQ>=40, AF>=0.9 ---
T5=$(ts)
if [ ! -f .filt.ok ]; then
  bcftools view -v snps,indels -i 'INFO/MQ>=40 && FORMAT/DP>=10' ${S}.raw.vcf.gz -Ou 2>/dev/null \
   | bcftools filter -e 'ALT="."' -Oz -o pass1.vcf.gz 2>/dev/null
  python3 -I - <<'PYF'
import gzip
keep=[]
for l in gzip.open('pass1.vcf.gz','rt'):
    if l[0]=='#': keep.append(l); continue
    f=l.rstrip('\n').split('\t'); d=dict(zip(f[8].split(':'),f[9].split(':')))
    if 'AD' not in d: continue
    ad=[int(x) for x in d['AD'].split(',') if x!='.']
    tot=sum(ad)
    if tot<10 or len(ad)<2: continue
    if max(ad[1:])/tot>=0.9: keep.append(l)
open('tmp.filt.vcf','w').writelines(keep)
PYF
  mv tmp.filt.vcf ${S}.filt.vcf && bgzip -f ${S}.filt.vcf && bcftools index -f ${S}.filt.vcf.gz && touch .filt.ok
fi
echo "filter done $(( $(ts)-T5 ))s"
echo "TOTAL $(( $(ts)-T0 ))s"
samtools flagstat ${S}.bam | head -5
