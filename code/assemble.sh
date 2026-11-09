#!/bin/bash
# Pipeline 2 stage 5b - de novo assembly. Identical treatment for every isolate.
#   assemble ALL quality-trimmed read pairs with SPAdes --isolate
#   then apply one contig filter to every sample:
#     length >= 500 bp, k-mer coverage >= 5, and a BLASTn hit to the MTBC reference
#     panel (H37Rv + M. canettii + 12 public MTBC genomes) of >=80% identity over >=200 bp.
#   The panel includes ancestral genomes, so sequence absent from H37Rv - TbD1 above all -
#   is retained; contigs from a non-mycobacterial organism are not.
set -u
S=$1
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
export PATH=$P/bin:$P/src/SPAdes-4.3.0-Linux/bin:$PATH
W=$P/work/$S; cd $W
[ -f .asm.ok ] && { echo "assembly already done"; exit 0; }
T0=$(date +%s)
# resumable: a container restart or an interrupted run costs only the current k-mer stage
if [ -d spades ] && [ -f spades/params.txt ]; then
  spades.py --continue -o spades >>spades.log 2>&1 || true
else
  spades.py --isolate -t 2 -m 6 -1 t1.fq.gz -2 t2.fq.gz -o spades >spades.log 2>&1 || true
fi
[ -s spades/contigs.fasta ] || { echo "SPADES_INCOMPLETE - rerun to continue"; exit 2; }
blastn -query spades/contigs.fasta -db $P/db/mtbc_panel -outfmt '6 qseqid pident length' \
       -evalue 1e-20 -num_threads 2 -max_target_seqs 1 > contig_hits.tsv 2>/dev/null
python3 -I - "$S" <<'PY'
import re,sys,collections
S=sys.argv[1]
ok=set()
for l in open('contig_hits.tsv'):
    q,pid,ln=l.split('\t')[:3]
    if float(pid)>=80 and int(ln)>=200: ok.add(q)
keep=[];drop=collections.Counter();name=None;seq=[]
def flush():
    global name,seq
    if name is None: return
    cid=name.split()[0]
    m=re.search(r'length_(\d+)_cov_([0-9.]+)',cid)
    L=int(m.group(1)) if m else len(''.join(seq)); cov=float(m.group(2)) if m else 0
    s=''.join(seq)
    if L<500: drop['short']+=1
    elif cov<5: drop['low_coverage']+=1
    elif cid not in ok: drop['no_MTBC_hit']+=1
    else: keep.append((cid,s))
for l in open('spades/contigs.fasta'):
    if l[0]=='>': flush(); name=l[1:].strip(); seq=[]
    else: seq.append(l.strip())
flush()
with open('%s_denovo_p2.fasta'%S,'w') as fh:
    for i,(n,s) in enumerate(keep,1):
        fh.write('>%s_contig%d %s\n'%(S,i,n))
        for j in range(0,len(s),80): fh.write(s[j:j+80]+'\n')
tot=sum(len(s) for _,s in keep)
print('contigs kept %d (%d bp) | dropped: %s'%(len(keep),tot,dict(drop)))
open('asm_filter.json','w').write('{"kept":%d,"total_bp":%d,"dropped_short":%d,"dropped_low_cov":%d,"dropped_no_MTBC_hit":%d}'
    %(len(keep),tot,drop['short'],drop['low_coverage'],drop['no_MTBC_hit']))
PY
rm -rf spades/K* spades/tmp spades/misc spades/corrected
touch .asm.ok
echo "assembly done $(( $(date +%s)-T0 ))s"
