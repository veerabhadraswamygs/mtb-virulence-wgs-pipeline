"""Recombination test (replaces Gubbins): Phi-style permutation test of the four-gamete incompatibility of NEIGHBOURING informative sites.
If recombination has occurred, neighbouring sites are more compatible than random pairs; clonal evolution (no recombination) gives no spatial pattern."""
import numpy as np,json
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c1'
seq={};cur=None
for l in open(O+'/aln55.var.fasta'):
    if l[0]=='>': cur=l[1:].strip(); seq[cur]=[]
    else: seq[cur].append(l.strip())
iso=[k for k in seq if k not in ('Mcanettii','MTB25')]
sj=json.load(open(O+'/aln55.sites.json')); sites=np.array(sj['sites'])
A=np.array([list(''.join(seq[k])) for k in iso])
keep=[];enc=[]
for i in range(A.shape[1]):
    col=A[:,i]; st=sorted(set(col)-{'N'})
    if len(st)!=2: continue
    cnt=[(col==s).sum() for s in st]
    if min(cnt)<2: continue
    e=np.full(len(iso),-1,np.int8); e[col==st[0]]=0; e[col==st[1]]=1
    keep.append(i); enc.append(e)
E_=np.array(enc); pos=sites[keep]; print('informative biallelic sites',len(pos))
def incompat(x,y):
    m=(x>=0)&(y>=0)
    if m.sum()<4: return 0
    g={(a,b) for a,b in zip(x[m],y[m])}
    return int(len(g)==4)
def stat(order):
    return np.mean([incompat(E_[order[i]],E_[order[i+1]]) for i in range(len(order)-1)])
rng=np.random.default_rng(20261103)
obs=stat(np.arange(len(pos)))
perm=[stat(rng.permutation(len(pos))) for _ in range(999)]
p=(1+sum(1 for x in perm if x<=obs))/1000      # recombination => neighbours MORE compatible (lower incompatibility)
# genome-wide incompatible fraction of random pairs
idx=rng.integers(0,len(pos),(20000,2)); rnd=np.mean([incompat(E_[a],E_[b]) for a,b in idx if a!=b])
res={'informative_sites':int(len(pos)),'adjacent_incompat_fraction':round(float(obs),4),'permuted_mean':round(float(np.mean(perm)),4),
     'permuted_5th_percentile':round(float(np.percentile(perm,5)),4),'p_value_recombination_signal':p,'random_pair_incompat_fraction':round(float(rnd),4),
     'interpretation':'no significant spatial clustering of compatibility (p>=0.05): consistent with clonal evolution, so no recombination masking was applied' if p>=0.05 else 'significant signal of recombination (p<0.05)'}
json.dump(res,open(O+'/c1_recombination_test.json','w'),indent=1); print(res)
