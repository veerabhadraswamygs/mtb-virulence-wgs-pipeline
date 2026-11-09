"""Analysis 24: pairwise SNP distances (pairwise-complete: positions missing in either isolate are not compared) and 5- / 12-SNP clusters."""
import numpy as np,json,os,re,csv,itertools
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c1'; L=4411532
seq={};cur=None
for l in open(O+'/aln55.var.fasta'):
    if l[0]=='>': cur=l[1:].strip(); seq[cur]=[]
    else: seq[cur].append(l.strip())
seq={k:''.join(v) for k,v in seq.items()}
iso=[k for k in seq if k!='Mcanettii']; iso.sort(key=lambda x:int(x[3:]))
can=np.frombuffer(seq['Mcanettii'].encode(),dtype=np.uint8)
ref=np.array(json.load(open(O+'/aln55.sites.json'))['ref'])
refb=np.frombuffer(''.join(ref).encode(),dtype=np.uint8)
# ingroup-variable sites only (exclude sites where only canettii differs)
A=np.array([np.frombuffer(seq[k].encode(),dtype=np.uint8) for k in iso])
N=ord('N'); known=A!=N
invar=((A!=refb)&known).any(axis=0)
A=A[:,invar]; known=known[:,invar]; print('ingroup variable sites',A.shape[1])
n=len(iso); D=np.zeros((n,n),int); S=np.zeros((n,n),int); U=np.zeros((n,n),int)
altm=(A!=refb[invar])&known
for i in range(n):
    for j in range(i+1,n):
        both=known[i]&known[j]; U[i,j]=U[j,i]=int(((altm[i]&~known[j])|(altm[j]&~known[i])).sum()); D[i,j]=D[j,i]=int(((A[i]!=A[j])&both).sum()); S[i,j]=S[j,i]=int(both.sum())
# comparable callable positions (soft-core) per pair
mask=np.unpackbits(np.load(O+'/mask.npy'))[:L+1].astype(bool)
cal=np.array([np.unpackbits(np.load('%s/cal_%s.npy'%(O,s)))[:L+1].astype(bool)&~mask for s in iso])
cs=json.load(open(O+'/callable_stats.json'))
comp=np.zeros((n,n),int)
for i in range(n):
    for j in range(i+1,n): comp[i,j]=comp[j,i]=int((cal[i]&cal[j]).sum())
with open(O+'/a24_snp_distance_matrix.tsv','w') as o:
    o.write('\t'+'\t'.join(iso)+'\n')
    for i in range(n): o.write(iso[i]+'\t'+'\t'.join(map(str,D[i]))+'\n')
with open(O+'/a24_comparable_positions.tsv','w') as o:
    o.write('\t'+'\t'.join(iso)+'\n')
    for i in range(n): o.write(iso[i]+'\t'+'\t'.join(map(str,comp[i]))+'\n')
LOW=[s for s in iso if cs[s]['callable_frac_unmasked']<0.5]   # MTB4, MTB25, MTB32
def clusters(th,M=None):
    M=D if M is None else M
    par=list(range(n))
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i in range(n):
        for j in range(i+1,n):
            if iso[i] in LOW or iso[j] in LOW: continue
            if M[i,j]<=th: par[f(i)]=f(j)
    g={}
    for i in range(n): g.setdefault(f(i),[]).append(iso[i])
    return [sorted(v,key=lambda x:int(x[3:])) for v in g.values() if len(v)>1]
res={'n_isolates':n,'ingroup_variable_sites':int(A.shape[1]),'not_assessable_for_clustering_callable_lt_50pct':LOW,
     'clusters_5':clusters(5),'clusters_12':clusters(12),'clusters_5_upper_bound':clusters(5,D+U),'clusters_12_upper_bound':clusters(12,D+U),
     'close_pairs':[{'a':iso[i],'b':iso[j],'snps':int(D[i,j]),'unresolved_positions':int(U[i,j]),'upper_bound':int(D[i,j]+U[i,j]),'comparable_bp':int(comp[i,j])} for i in range(n) for j in range(i+1,n) if D[i,j]<=12 and iso[i] not in LOW and iso[j] not in LOW]}
# closest partner of the low-coverage isolates (lower bounds)
res['low_coverage_closest']={}
for s in LOW:
    i=iso.index(s); k=[j for j in range(n) if j!=i]; jmin=min(k,key=lambda j:D[i,j])
    res['low_coverage_closest'][s]={'closest':iso[jmin],'snps':int(D[i,jmin]),'comparable_sites':int(S[i,jmin]),'comparable_positions':int(comp[i,jmin])}
tri=[D[i,j] for i in range(n) for j in range(i+1,n) if iso[i] not in LOW and iso[j] not in LOW]
res['median_pairwise']=float(np.median(tri)); res['min_pairwise']=int(min(tri))
for g in (5,12): print(g,res['clusters_%d'%g]); print(g,'upper',res['clusters_%d_upper_bound'%g])
print(sorted(set(c['unresolved_positions'] for c in res['close_pairs'])))
json.dump(res,open(O+'/a24_summary.json','w'),indent=1)
for a,b in [('MTB69','MTB78'),('MTB69','MTB131'),('MTB78','MTB131')]:
    i,j=iso.index(a),iso.index(b); print(a,b,'SNPs',D[i,j],'comparable sites',S[i,j],'comparable bp',comp[i,j])
print(res['low_coverage_closest']); print('median',res['median_pairwise'],'min',res['min_pairwise'])
