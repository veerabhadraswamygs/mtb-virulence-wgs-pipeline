"""C3 step 1: per-gene codon alignments (SNP-based, 54 isolates, MTB25 excluded) + pruned unrooted trees. 
Codons containing any masked position are removed for all taxa; codons with any non-callable position (depth<10 etc.) or an N at a variable site are set to '---' for that isolate;
premature stop codons are set to '---' for that isolate (counted). Reference terminal stop removed."""
import numpy as np,json,os,re,bisect,collections,sys
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; C1=P+'/cohort/c1'; O=P+'/cohort/c3'; L=4411532
os.makedirs(O+'/aln',exist_ok=True)
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
use=[s for s in samples if s!='MTB25']
ref=''.join(l.strip() for l in open(P+'/ref/H37Rv.fa') if l[0]!='>').upper(); assert len(ref)==L
G=json.load(open(P+'/ref/genes.json')); cds={rv:g for rv,g in G.items() if g.get('type')=='CDS'}
mask=np.unpackbits(np.load(C1+'/mask.npy'))[:L+1].astype(bool)
cal={s:np.unpackbits(np.load('%s/cal_%s.npy'%(C1,s)))[:L+1].astype(bool)&~mask for s in use}
seq={};cur=None
for l in open(C1+'/aln55.var.fasta'):
    if l[0]=='>': cur=l[1:].strip(); seq[cur]=[]
    else: seq[cur].append(l.strip())
seq={k:''.join(v) for k,v in seq.items()}
sj=json.load(open(C1+'/aln55.sites.json')); sites=sj['sites']
assert all(s in seq for s in use) and len(seq[use[0]])==len(sites)
comp=str.maketrans('ACGT','TGCA')
code={}
b='TCAG'; aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i=0
for x in b:
    for y in b:
        for z in b: code[x+y+z]=aa[i]; i+=1
stops={'TAA','TAG','TGA'}
# ---- tree utils
def parse(s):
    s=s.strip().rstrip(';'); pos=0
    def node():
        nonlocal pos
        if s[pos]=='(':
            pos+=1; ch=[]
            while True:
                ch.append(node())
                if s[pos]==',': pos+=1; continue
                if s[pos]==')': pos+=1; break
            while pos<len(s) and s[pos] not in ':,);': pos+=1   # support label
            n={'ch':ch,'name':None,'len':0.0}
        else:
            st=pos
            while s[pos] not in ':,);': pos+=1
            n={'ch':[],'name':s[st:pos],'len':0.0}
        if pos<len(s) and s[pos]==':':
            pos+=1; st=pos
            while s[pos] not in ',);': pos+=1
            n['len']=float(s[st:pos])
        return n
    return node()
def prune(n,keep):
    if not n['ch']: return n if n['name'] in keep else None
    ch=[c for c in (prune(c,keep) for c in n['ch']) if c is not None]
    if not ch: return None
    if len(ch)==1: ch[0]['len']+=n['len']; return ch[0]
    n['ch']=ch; return n
def unroot(n):
    while len(n['ch'])==2:
        a,b=n['ch']
        if a['ch']: big,oth=a,b
        elif b['ch']: big,oth=b,a
        else: break
        oth['len']+=big['len']; n['ch']=big['ch']+[oth]
    return n
def newick(n):
    if not n['ch']: return '%s:%.8f'%(n['name'],max(n['len'],1e-7))
    return '('+','.join(newick(c) for c in n['ch'])+')'+(':%.8f'%max(n['len'],1e-7) if n.get('len') else '')
import copy
tree0=parse(open(C1+'/tree55.treefile').read())
rows=[]
col_of={s:np.frombuffer(seq[s].encode(),dtype=np.uint8) for s in use}
sites_a=np.array(sites)
for rv,g in sorted(cds.items()):
    a,bb,strand=g['start'],g['end'],g['strand']; ln=bb-a+1
    if ln%3: rows.append({'rv':rv,'status':'length not multiple of 3'}); continue
    unm=(~mask[a:bb+1]); n_unm=int(unm.sum())
    if n_unm<0.9*ln: rows.append({'rv':rv,'len':ln,'status':'<90% unmasked'}); continue
    usable=[s for s in use if cal[s][a:bb+1].sum()>=0.9*ln]
    if len(usable)<10: rows.append({'rv':rv,'len':ln,'n_usable':len(usable),'status':'<10 usable isolates'}); continue
    j0=bisect.bisect_left(sites,a); j1=bisect.bisect_right(sites,bb)
    nc=ln//3
    # codon positions are in transcript order; build forward-strand arrays then revcomp whole sequence
    fw=np.array(list(ref[a-1:bb]))
    keepc=np.ones(nc,bool)
    # masked codons: on forward coordinates codon k in transcript order covers forward block
    # work in forward coords: block index i//3 (for - strand reversed later)
    blk=np.arange(ln)//3
    mk=np.zeros(nc,bool); np.logical_or.at(mk,blk,mask[a:bb+1])
    seqs={}; nstop=0; nmiss=0
    for s in usable:
        arr=fw.copy()
        bad=np.zeros(nc,bool); np.logical_or.at(bad,blk,~cal[s][a:bb+1])
        for j in range(j0,j1):
            c=chr(col_of[s][j])
            p=sites[j]-a
            if c in 'ACGT': arr[p]=c
            else: bad[p//3]=True
        st=''.join(arr)
        if strand=='-': # reverse complement; codon k (transcript) = forward block nc-1-k
            st=st[::-1].translate(comp); bad=bad[::-1]; mkk=mk[::-1]
        else: mkk=mk
        cod=[st[3*k:3*k+3] for k in range(nc)]
        for k in range(nc):
            if bad[k]: cod[k]='---'; nmiss+=1
            elif cod[k] in stops and k<nc-1: cod[k]='---'; nstop+=1
        seqs[s]=cod
    if strand=='-': mkk=mk[::-1]
    else: mkk=mk
    last_stop=any(seqs[s][-1] in stops for s in usable) or ''.join(fw)[-3:] in stops or True
    keep=[k for k in range(nc-1) if not mkk[k]]   # drop terminal codon (stop)
    if len(keep)<0.5*(nc-1): rows.append({'rv':rv,'len':ln,'n_usable':len(usable),'status':'<50% codons retained'}); continue
    final={s:''.join(seqs[s][k] for k in keep) for s in usable}
    # variable codons (among non-missing)
    var=0; nsyn_cod=0
    ncod=len(keep)
    kv=[]
    for i,k in enumerate(keep):
        cs={seqs[s][k] for s in usable if seqs[s][k]!='---'}
        if len(cs)>1: kv.append(i)
    # classify variable codons by pairwise changes relative to reference codon: count syn/nonsyn alleles
    ncs=0; nns=0
    for i in kv:
        k=keep[i]; cs=collections.Counter(seqs[s][k] for s in usable if seqs[s][k]!='---')
        refc=max(cs,key=cs.get)
        for c in cs:
            if c==refc: continue
            d=[x for x in range(3) if c[x]!=refc[x]]
            if code.get(c)==code.get(refc): ncs+=1
            else: nns+=1
    tn=prune(copy.deepcopy(tree0),set(usable)); tn=unroot(tn)
    row={'rv':rv,'gene':g.get('name',''),'len':ln,'codons_kept':ncod,'n_usable':len(usable),'variable_codons':len(kv),'alt_codon_alleles_syn':ncs,'alt_codon_alleles_nonsyn':nns,'premature_stop_codons_set_missing':nstop,'status':'ok' if kv else 'invariant'}
    rows.append(row)
    if kv:
        with open('%s/aln/%s.phy'%(O,rv),'w') as o:
            o.write(' %d %d\n'%(len(usable),3*ncod))
            for s in usable: o.write('%s  %s\n'%(s,final[s]))
        with open('%s/aln/%s.fas'%(O,rv),'w') as o:
            for s in usable: o.write('>%s\n%s\n'%(s,final[s]))
        open('%s/aln/%s.nwk'%(O,rv),'w').write(newick(tn)+';\n')
import csv
keys=['rv','gene','len','codons_kept','n_usable','variable_codons','alt_codon_alleles_syn','alt_codon_alleles_nonsyn','premature_stop_codons_set_missing','status']
with open(O+'/c3_gene_alignment_summary.tsv','w') as o:
    w=csv.DictWriter(o,fieldnames=keys,delimiter='\t',restval=''); w.writeheader(); w.writerows(rows)
print(collections.Counter(r['status'] for r in rows))
