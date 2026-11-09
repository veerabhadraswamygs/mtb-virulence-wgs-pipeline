"""C1 step 2: masked SNP alignment. Soft-core = unmasked positions callable (depth>=10) in >=51/55 isolates.
Missing (N) = low coverage, or any raw variant call that failed the depth/MQ/AF filters (sub-threshold variant), or inside an indel.
Taxa: 55 isolates + M. canettii (outgroup) [+ 64 public lineage-1 assemblies for the expanded alignment]."""
import numpy as np,subprocess,os,re,json,gzip,collections
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
O=P+'/cohort/c1'; L=4411532; M=P+'/src/minimap2-2.31_x64-linux'
REF=''.join(l.strip() for l in open(P+'/ref/H37Rv.fa') if not l.startswith('>')).upper(); assert len(REF)==L
mask=np.unpackbits(np.load(O+'/mask.npy'))[:L+1].astype(bool)
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
cnt=np.zeros(L+1,np.int16); cal={}
for s in samples:
    c=np.unpackbits(np.load('%s/cal_%s.npy'%(O,s)))[:L+1].astype(bool); cal[s]=c; cnt+=c
SOFT=(cnt>=51)&~mask; SOFT[0]=False
print('soft-core positions',int(SOFT.sum()),'strict core (all 55)',int(((cnt==55)&~mask).sum()))
# local isolate calls
snp={}; miss={}
for s in samples:
    sn={}; ms=set()
    filt={}
    for l in gzip.open('%s/work/%s/%s.filt.vcf.gz'%(P,s,s),'rt'):
        if l[0]=='#': continue
        f=l.split('\t',6); p=int(f[1]); filt[p]=(f[3],f[4])
        if len(f[3])==1 and len(f[4])==1: sn[p]=f[4]
        else:
            for q in range(p,p+len(f[3])): ms.add(q)
    rawv=subprocess.run("zcat %s/work/%s/%s.raw.vcf.gz | awk '!/^#/ && $5!=\".\" {print $2}'"%(P,s,s),shell=True,capture_output=True,text=True).stdout.split()
    for x in rawv:
        if int(x) not in filt: ms.add(int(x))
    snp[s]=sn; miss[s]=ms
# canettii + public
def aln_asm(fa,tag):
    j='%s/asm_%s.json'%(O,tag)
    if os.path.exists(j): return json.load(open(j))
    paf=subprocess.run("%s/minimap2 -cx asm20 --cs %s/ref/H37Rv.fa %s 2>/dev/null | sort -k6,6 -k8,8n"%(M,P,fa),shell=True,capture_output=True,text=True).stdout
    vcf=subprocess.run("%s/k8 %s/paftools.js call -f %s/ref/H37Rv.fa -l 1000 -L 1000 -"%(M,M,P),shell=True,input=paf,capture_output=True,text=True).stdout
    cov=sorted((int(f[7]),int(f[8])) for f in (l.split('\t') for l in paf.splitlines()) if int(f[11])>=20 and int(f[10])>=1000)
    mg=[]
    for s,e in cov:
        if mg and s<=mg[-1][1]: mg[-1][1]=max(mg[-1][1],e)
        else: mg.append([s,e])
    sn={};ms=[]
    for l in vcf.splitlines():
        if l[0]=='#': continue
        f=l.split('\t',6); p=int(f[1])
        if len(f[3])==1 and len(f[4])==1 and f[4]!='.': sn[p]=f[4]
        else: ms.append((p,p+len(f[3])))
    d={'snp':{str(k):v for k,v in sn.items()},'miss':ms,'cov':mg}; json.dump(d,open(j,'w')); return d
def to_arrays(d):
    c=np.zeros(L+1,bool)
    for s,e in d['cov']: c[s+1:e+1]=True
    for s,e in d['miss']: c[s:e]=False
    return c,{int(k):v for k,v in d['snp'].items()}
can=aln_asm(P+'/ref/canettii.fa','canettii'); cc,cs=to_arrays(can)
print('canettii covered soft-core',int((cc&SOFT).sum()),'SNPs in soft-core',sum(1 for p in cs if SOFT[p]))
pubL1=[]
for f in sorted(os.listdir(P+'/pub/calls')):
    if f.endswith('.json'):
        d=json.load(open(P+'/pub/calls/'+f))
        if d['lineage'].startswith('lineage1'): pubL1.append(d['acc'])
print('public L1',len(pubL1))
pub={}
for a in pubL1:
    d=json.load(open(P+'/pub/calls/%s.json'%a)); c=np.zeros(L+1,bool)
    for s,e in d['intervals']: c[s+1:e+1]=True
    sn={}
    for l in gzip.open(P+'/pub/calls/%s.vcf.gz'%a,'rt'):
        if l[0]=='#': continue
        f=l.split('\t',6); p=int(f[1])
        if len(f[3])==1 and len(f[4])==1 and f[4]!='.': sn[p]=f[4]
        else:
            for q in range(p,p+len(f[3])): c[q]=False
    pub[a]=(c,sn)
def build(taxa,name,vary_from):
    # taxa: list of (label, callable_array, snp_dict)
    var=set()
    for lab,c,sn in taxa:
        if lab in vary_from:
            for p,a in sn.items():
                if SOFT[p] and c[p]: var.add(p)
    sites=sorted(var)
    seqs={}
    for lab,c,sn in taxa:
        out=[]
        for p in sites:
            if not c[p]: out.append('N')
            else: out.append(sn.get(p,REF[p-1]))
        seqs[lab]=''.join(out)
    # constant-site counts: soft-core positions not variable
    constmask=SOFT.copy()
    for p in sites: constmask[p]=False
    idx=np.nonzero(constmask)[0]
    cn=collections.Counter(REF[i-1] for i in idx)
    with open('%s/%s.var.fasta'%(O,name),'w') as o:
        for lab,sq in seqs.items(): o.write('>%s\n%s\n'%(lab,sq))
    json.dump({'sites':sites,'ref':[REF[p-1] for p in sites],'const':{k:cn[k] for k in 'ACGT'},'n_sites':len(sites),
               'taxa':list(seqs),'missing_frac':{k:round(v.count('N')/max(1,len(v)),4) for k,v in seqs.items()}},open('%s/%s.sites.json'%(O,name),'w'))
    print(name,len(seqs),'taxa',len(sites),'variable sites; const',dict(cn))
loc=[(s,cal[s]&~np.zeros(L+1,bool),{p:a for p,a in snp[s].items()}) for s in samples]
# local callable: remove raw-unreliable positions and indel spans
loc2=[]
for s,c,sn in loc:
    c=c.copy()
    for p in miss[s]:
        if p<=L: c[p]=False
    loc2.append((s,c,sn))
can_t=('Mcanettii',cc,cs)
# MTB25: 77.6% of reads are foreign (ANI 74.9%); its ~47k private SNPs are contamination artefacts, so sites are defined WITHOUT MTB25; MTB25 is still genotyped (flagged) at those sites
build(loc2+[can_t],'aln55',(set(samples)-{'MTB25'})|{'Mcanettii'})
pubt=[('PUB_'+a,pub[a][0],pub[a][1]) for a in pubL1]
build(loc2+pubt+[can_t],'alnX',(set(samples)-{'MTB25'})|{'PUB_'+a for a in pubL1}|{'Mcanettii'})
json.dump({'soft_core_positions':int(SOFT.sum()),'strict_core_positions':int(((cnt==55)&~mask).sum()),'masked_bp':int(mask.sum()),'public_L1':pubL1},open(O+'/c1_core_summary.json','w'))
