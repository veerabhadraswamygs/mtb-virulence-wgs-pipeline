"""E12: T-cell epitope density vs dN/dS. (a) gene level: epitope density (fraction of codons in mapped epitopes) vs M0 omega (Spearman, genes with >=1 epitope & omega defined);
(b) within-gene pooled pN/pS in epitope vs non-epitope codons (alt codon alleles / NG86 sites), permutation test shuffling epitope labels among codons within each gene."""
import os,json,csv,collections,numpy as np
from scipy import stats
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c3'; C1=P+'/cohort/c1'; L=4411532
import importlib.util
G=json.load(open(P+'/ref/genes.json')); cds={rv:g for rv,g in G.items() if g.get('type')=='CDS'}
ref=''.join(l.strip() for l in open(P+'/ref/H37Rv.fa') if l[0]!='>').upper()
mask=np.unpackbits(np.load(C1+'/mask.npy'))[:L+1].astype(bool)
comp=str.maketrans('ACGT','TGCA'); b='TCAG'; aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; code={};i=0
for x in b:
    for y in b:
        for z in b: code[x+y+z]=aa[i]; i+=1
def sites(c):
    S=0.0
    for p in range(3):
        for n in 'ACGT':
            if n==c[p]: continue
            m=c[:p]+n+c[p+1:]
            if code[m]=='*': continue
            if code[m]==code[c]: S+=1/3
    return 9-S-0 if False else (S, 3-S)   # (syn sites, nonsyn sites) NG86 (stop changes excluded roughly)
ep=collections.defaultdict(list)
for e in json.load(open(P+'/ref/epitopes_mapped.json')): ep[e['rv']].append((e['aa_start'],e['aa_end']))
rows=[];perm_data=[]
for rv in sorted(ep):
    f='%s/aln/%s.fas'%(O,rv)
    if rv not in cds or not os.path.exists(f): continue
    g=cds[rv]; a,bb=g['start'],g['end']; nc=(bb-a+1)//3
    m=np.zeros(nc,bool); np.logical_or.at(m,np.arange(bb-a+1)//3,mask[a:bb+1])
    fw=ref[a-1:bb]
    if g['strand']=='-': fw=fw.translate(comp)[::-1]; m=m[::-1]
    keep=[k for k in range(nc-1) if not m[k]]
    seqs=[]; cur=None
    for l in open(f):
        if l[0]!='>': seqs.append(l.strip())
    isep=np.zeros(len(keep),bool)
    for i,k in enumerate(keep):
        aap=k+1
        isep[i]=any(s<=aap<=e for s,e in ep[rv])
    N=np.zeros(len(keep));S=np.zeros(len(keep));nn=np.zeros(len(keep));ns=np.zeros(len(keep))
    for i,k in enumerate(keep):
        rc=fw[3*k:3*k+3]
        if any(c not in 'ACGT' for c in rc) or rc not in code: continue
        S[i],N[i]=sites(rc)
        alts={s[3*i:3*i+3] for s in seqs}-{'---',rc}
        for c in alts:
            if c in code and '-' not in c:
                if code[c]==code[rc]: ns[i]+=1
                else: nn[i]+=1
    n_ep=int(isep.sum())
    if n_ep==0: continue
    def pp(sel):
        return (nn[sel].sum(),ns[sel].sum(),N[sel].sum(),S[sel].sum())
    rows.append({'rv':rv,'gene':g.get('name',''),'codons':len(keep),'epitope_codons':n_ep,'epitope_density':round(n_ep/len(keep),4),'ep_nonsyn':nn[isep].sum(),'ep_syn':ns[isep].sum(),'non_nonsyn':nn[~isep].sum(),'non_syn':ns[~isep].sum()})
    perm_data.append((nn,ns,N,S,isep))
def pNpS(a,b,c,d): return (a/c)/(b/d) if b>0 and c>0 and d>0 else float('nan')
def pooled(sel_epitope):
    A=B=Cc=D=0.0
    for nn,ns,N,S,isep in perm_data:
        s=isep if sel_epitope else ~isep
        A+=nn[s].sum();B+=ns[s].sum();Cc+=N[s].sum();D+=S[s].sum()
    return A,B,Cc,D
Ae=pooled(True); An=pooled(False)
pe,pn=pNpS(*Ae),pNpS(*An)
rng=np.random.default_rng(20261103); cnt=0; null=[]
for _ in range(2000):
    AE=[0,0,0,0]; AN=[0,0,0,0]
    for nn,ns,N,S,isep in perm_data:
        s=rng.permutation(isep)
        for arr,sel in ((AE,s),(AN,~s)):
            arr[0]+=nn[sel].sum();arr[1]+=ns[sel].sum();arr[2]+=N[sel].sum();arr[3]+=S[sel].sum()
    d=pNpS(*AE)-pNpS(*AN); null.append(d)
obs=pe-pn; null=np.array(null); p=(np.sum(np.abs(null[~np.isnan(null)])>=abs(obs))+1)/(np.sum(~np.isnan(null))+1)
# gene-level
om={}
for r in csv.DictReader(open(O+'/e15_genomewide_dNdS_per_gene.tsv'),delimiter='\t'):
    if r['omega_M0'] not in ('','NA (dS=0)') and float(r['omega_M0'])<50: om[r['rv']]=float(r['omega_M0'])
x=[r['epitope_density'] for r in rows if r['rv'] in om]; y=[om[r['rv']] for r in rows if r['rv'] in om]
rho,pr=stats.spearmanr(x,y) if len(x)>5 else (float('nan'),float('nan'))
for r in rows: r['omega_M0']=om.get(r['rv'],'')
with open(O+'/e12_epitope_density_vs_dNdS_per_gene.tsv','w') as o:
    w=csv.DictWriter(o,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
res={'genes_with_epitopes_and_alignment':len(rows),'epitope_codons_total':int(sum(r['epitope_codons'] for r in rows)),'epitope':{'nonsyn_alleles':Ae[0],'syn_alleles':Ae[1],'N_sites':round(Ae[2],1),'S_sites':round(Ae[3],1),'pN_pS':round(pe,4)},
     'non_epitope':{'nonsyn_alleles':An[0],'syn_alleles':An[1],'N_sites':round(An[2],1),'S_sites':round(An[3],1),'pN_pS':round(pn,4)},
     'difference_epitope_minus_non':round(obs,4),'permutation_p_within_gene_label_shuffle_2000':float(p),
     'gene_level_spearman':{'n_genes':len(x),'rho':round(float(rho),4),'p':float(pr)}}
json.dump(res,open(O+'/e12_summary.json','w'),indent=1); print(json.dumps(res,indent=1))
