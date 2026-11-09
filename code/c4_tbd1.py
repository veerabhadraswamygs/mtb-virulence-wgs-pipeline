"""C4 analysis 28: TbD1 locus (mmpS6/mmpL6) across lineage 1. Locus recovered from de novo assemblies (local) / public assemblies by BLASTn against the 2,139 bp M. canettii TbD1 segment (no H37Rv mapping).
Segments are tiled in reference order; unrecovered stretches are N. MAFFT --auto with the canettii reference as sequence 1."""
import json,os,re,subprocess,gzip,collections,csv,glob
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c4'; TL=2139
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
ref=''.join(l.strip() for l in open(P+'/ref/TbD1.fa') if l[0]!='>').upper(); assert len(ref)==TL
comp=str.maketrans('ACGTN','TGCAN')
def readfa(f):
    op=gzip.open if f.endswith('.gz') else open; d={};cur=None
    for l in op(f,'rt'):
        if l[0]=='>': cur=l[1:].split()[0]; d[cur]=[]
        else: d[cur].append(l.strip())
    return {k:''.join(v).upper() for k,v in d.items()}
def locus(fa_path,tag):
    tmp='%s/tmp_%s.fa'%(O,tag)
    if fa_path.endswith('.gz'):
        with gzip.open(fa_path,'rt') as i,open(tmp,'w') as o: o.write(i.read())
    else: tmp=fa_path
    out=subprocess.run("%s/bin/blastn -query %s -db %s/db/tbd1 -outfmt '6 qseqid qstart qend sstart send pident length' -evalue 1e-20 -num_threads 2"%(P,tmp,P),shell=True,capture_output=True,text=True).stdout
    hits=[]
    for l in out.splitlines():
        q,qs,qe,ss,se,pid,ln=l.split('\t'); qs,qe,ss,se,ln=map(int,(qs,qe,ss,se,ln))
        if float(pid)>=90 and ln>=200: hits.append((q,qs,qe,ss,se,float(pid),ln))
    if not hits:
        if tmp!=fa_path: os.remove(tmp)
        return None
    seqs=readfa(tmp)
    # greedy tiling by length
    acc=[]
    for h in sorted(hits,key=lambda x:-x[6]):
        a,b=sorted((h[3],h[4]))
        if all(min(b,y)-max(a,x)<0.3*(b-a+1) for (x,y,_) in acc): acc.append((a,b,h))
    acc.sort(key=lambda x:x[0])
    # build consensus in reference coordinates: segment seq oriented to ref strand, trimmed of reference overlap with previous
    pieces=[];cur=1;covered=0
    for a,b,h in acc:
        q,qs,qe,ss,se,pid,ln=h
        s=seqs[q][min(qs,qe)-1:max(qs,qe)]
        if ss>se: s=s.translate(comp)[::-1]
        if a>cur: pieces.append('N'*(a-cur))
        if a<cur: # overlapping previous: trim
            trim=cur-a; s=s[trim:] ; a=cur
        pieces.append(s); cur=b+1; covered+=b-a+1
    if cur<=TL: pieces.append('N'*(TL-cur+1))
    if tmp!=fa_path: os.remove(tmp)
    return ''.join(pieces),len(acc),covered/TL,sorted({h[0] for _,_,h in acc})
res={}
loc={}
for s in samples:
    r=json.load(open('%s/work/%s/results.json'%(P,s)))
    res[s]={'lineage':r['lineage']['lineage'],'tbd1_call':r['tbd1'].get('call'),'frac':r['tbd1'].get('fraction_of_TbD1_recovered'),'flank':r['tbd1'].get('flank_fraction_covered'),'source':'local'}
    if r['tbd1'].get('call')=='present':
        x=locus('%s/work/%s/%s_denovo_p2.fasta'%(P,s,s),s)
        if x: loc[s]=x
pub=[]
for f in sorted(glob.glob(P+'/pub/calls/*.json')):
    j=json.load(open(f)); lin=j.get('lineage'); acc=os.path.basename(f)[:-5]
    if lin and str(lin).startswith('lineage1'):
        x=locus('%s/pub/fna/%s.fna.gz'%(P,acc),acc)
        res[acc]={'lineage':lin,'tbd1_call':'present' if x and x[2]>=0.8 else ('partial' if x else 'absent'),'frac':round(x[2],3) if x else 0,'flank':'','source':'public'}
        if x: loc[acc]=x
# fasta: reference + L1 isolates with >=80% recovered
sel=[k for k,v in loc.items() if v[2]>=0.8 and res[k]['lineage'].startswith('lineage1')]
sel_local=[k for k in sel if k.startswith('MTB')]; sel_pub=[k for k in sel if not k.startswith('MTB')]
def write_aln(names,tag):
    with open('%s/tbd1_%s_unaligned.fa'%(O,tag),'w') as o:
        o.write('>canettii_ref\n%s\n'%ref)
        for k in names: o.write('>%s\n%s\n'%(k,loc[k][0]))
    out=subprocess.run('%s/bin/mafft --auto --thread 2 %s/tbd1_%s_unaligned.fa'%(P,O,tag),shell=True,capture_output=True,text=True,env={**os.environ,'MAFFT_BINARIES':P+'/libexec/mafft','PATH':P+'/bin:'+os.environ['PATH']})
    open('%s/tbd1_%s_aligned.fa'%(O,tag),'w').write(out.stdout)
    return out.stderr[-300:]
e1=write_aln(sel_local,'L1_local'); e2=write_aln(sel_local+sel_pub,'L1_local_plus_public')
json.dump({'res':res,'n_local_L1_with_locus':len(sel_local),'n_public_L1_with_locus':len(sel_pub),'n_pieces':{k:loc[k][1] for k in loc},'frac':{k:round(loc[k][2],3) for k in loc},'contigs':{k:loc[k][3] for k in loc}},open(O+'/tbd1_extract.json','w'),indent=1)
print(len(sel_local),len(sel_pub),e1,e2)
