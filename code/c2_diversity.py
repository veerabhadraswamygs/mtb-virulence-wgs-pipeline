"""Analysis 18 (pi, Tajima's D per gene; panel vs genome background) and E4 (COG-category burden of unique non-synonymous polymorphic sites)."""
import numpy as np,json,os,re,csv,collections,math
from scipy import stats
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c2'; C1=P+'/cohort/c1'; L=4411532
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
R={s:json.load(open('%s/work/%s/results.json'%(P,s))) for s in samples}
def coarse(l):
    m=re.match(r'lineage(\d)',l); return 'L'+m.group(1) if m else 'unassigned'
LIN={s:coarse(R[s]['lineage']['lineage']) for s in samples}
EXC={'MTB25'}; use=[s for s in samples if s not in EXC]
G=json.load(open(P+'/ref/genes.json')); cds={rv:g for rv,g in G.items() if g.get('type')=='CDS'}
panel=set(); 
import glob
for f in glob.glob(P+'/ref/panel/*.bed'):
    for l in open(f): panel.add(l.split('\t')[3].split('_')[0])
mask=np.unpackbits(np.load(C1+'/mask.npy'))[:L+1].astype(bool)
cum_um=np.cumsum(~mask)   # unmasked prefix
def seg(a,b): return cum_um[b]-cum_um[a-1]
cc={}
for s in use:
    c=np.unpackbits(np.load('%s/cal_%s.npy'%(C1,s)))[:L+1].astype(bool)&~mask
    cc[s]=np.cumsum(c)
genes=sorted(cds)
glen={rv:int(seg(cds[rv]['start'],cds[rv]['end'])) for rv in genes}
usable={rv:[s for s in use if glen[rv]>0 and (cc[s][cds[rv]['end']]-cc[s][cds[rv]['start']-1])/glen[rv]>=0.9] for rv in genes}
# alignment
seq={};cur=None
for l in open(C1+'/aln55.var.fasta'):
    if l[0]=='>': cur=l[1:].strip(); seq[cur]=[]
    else: seq[cur].append(l.strip())
seq={k:''.join(v) for k,v in seq.items()}
sj=json.load(open(C1+'/aln55.sites.json')); sites=sj['sites']
col={s:np.frombuffer(seq[s].encode(),dtype=np.uint8) for s in use}
idx={s:i for i,s in enumerate(use)}
M=np.array([col[s] for s in use])
gid=np.full(L+2,-1,np.int32); gi={rv:i for i,rv in enumerate(genes)}
for rv in genes: gid[cds[rv]['start']:cds[rv]['end']+1]=gi[rv]
sites_by_gene=collections.defaultdict(list)
for j,p in enumerate(sites):
    g=gid[p]
    if g>=0: sites_by_gene[genes[g]].append(j)
def tajima(S,pi_sum,n):
    if n<4 or S==0: return None,None
    a1=sum(1/i for i in range(1,n)); a2=sum(1/i**2 for i in range(1,n))
    b1=(n+1)/(3*(n-1)); b2=2*(n*n+n+3)/(9*n*(n-1))
    c1=b1-1/a1; c2=b2-(n+2)/(a1*n)+a2/a1**2
    e1=c1/a1; e2=c2/(a1**2+a2)
    D=(pi_sum-S/a1)/math.sqrt(e1*S+e2*S*(S-1))
    Dmin=(2/n-1/a1)/math.sqrt(e2) if False else (2/n-1/a1)/math.sqrt(e2)  # Tajima 1989 bounds
    Dmax=((n+1)/(2*n)-1/a1)/math.sqrt(e2)
    # beta approximation
    try:
        a_=-(1+Dmin*Dmax)*Dmax/(Dmax-Dmin); b_=(1+Dmin*Dmax)*Dmin/(Dmax-Dmin)
        x=(D-Dmin)/(Dmax-Dmin); p=stats.beta.cdf(x,b_,a_)   # Tajima 1989: x ~ Beta(b_, a_) so that E[D]=0
        p2=2*min(p,1-p)
    except Exception: p2=None
    return D,p2
def compute(group,tag):
    rows=[]
    gset=set(group)
    for rv in genes:
        us=[s for s in usable[rv] if s in gset]; n=len(us)
        base={'rv':rv,'gene':cds[rv].get('name',''),'panel':rv in panel,'unmasked_bp':glen[rv],'n_isolates_usable':n}
        if n<10 or glen[rv]<150: rows.append({**base,'pi_per_kb':'','S':'','tajimas_D':'','D_p':'','status':'not estimable'}); continue
        ii=[idx[s] for s in us]; pis=0.0; S=0
        for j in sites_by_gene.get(rv,[]):
            c=M[ii,j]; c=c[c!=ord('N')]
            if len(c)<max(10,int(0.8*n)): continue
            vals,cnt=np.unique(c,return_counts=True)
            if len(vals)<2: continue
            m=len(c); pr=cnt/m; pis+=(m/(m-1))*(1-(pr**2).sum()); S+=1
        D,p=tajima(S,pis,n)
        rows.append({**base,'pi_per_kb':round(1000*pis/glen[rv],4),'S':S,'tajimas_D':round(D,3) if D is not None else '','D_p':round(p,4) if p is not None else '','status':'ok'})
    with open('%s/a18_pi_tajimasD_%s.tsv'%(O,tag),'w') as o:
        w=csv.DictWriter(o,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
    ok=[r for r in rows if r['status']=='ok']
    pa=np.array([r['pi_per_kb'] for r in ok if r['panel']]); bg=np.array([r['pi_per_kb'] for r in ok if not r['panel']])
    res={'group':tag,'n_isolates':len(gset),'genes_estimable_panel':len(pa),'genes_estimable_background':len(bg),
         'median_pi_per_kb_panel':round(float(np.median(pa)),4),'median_pi_per_kb_background':round(float(np.median(bg)),4),
         'mean_pi_per_kb_panel':round(float(pa.mean()),4),'mean_pi_per_kb_background':round(float(bg.mean()),4)}
    u,p=stats.mannwhitneyu(pa,bg,alternative='two-sided'); res['mannwhitney_p']=float(p); res['rank_biserial']=round(float(2*u/(len(pa)*len(bg))-1),3)
    dp=np.array([r['tajimas_D'] for r in ok if r['panel'] and r['tajimas_D']!='']); db=np.array([r['tajimas_D'] for r in ok if not r['panel'] and r['tajimas_D']!=''])
    if len(dp)>5 and len(db)>5: res['median_D_panel']=round(float(np.median(dp)),3); res['median_D_background']=round(float(np.median(db)),3); res['genes_D_significant_panel']=int(sum(1 for r in ok if r['panel'] and r['D_p']!='' and r['D_p']<0.05))
    return res
out={'a18':[]}
for tag,grp in [('pooled',use)]+[(l,[s for s in use if LIN[s]==l]) for l in ('L1','L3','L4')]:
    r=compute(grp,tag); out['a18'].append(r); print(r)
# ---- E4
cog=json.load(open(P+'/ref/cog_h37rv.json'))
nonsyn=set(); sitegene=collections.defaultdict(set)
for s in use:
    for v in R[s]['variants']:
        if v['type']=='snp' and not v.get('in_mask') and v.get('rv') and any(k in v['effect'] for k in ('missense','stop_gained','start_lost','stop_lost')):
            sitegene[v['rv']].add((v['pos'],v['alt']))
tot_len=sum(glen[g] for g in genes if glen[g]>0); tot_ns=sum(len(sitegene[g]) for g in genes)
bg_rate=tot_ns/tot_len
def rr_test(k,len_,K,Ltot):
    # exact conditional binomial: k events in genes of total length len_, versus rest (K-k) in Ltot-len_
    p0=len_/Ltot; pv=stats.binomtest(k,K,p0).pvalue
    rest_k=K-k; rest_len=Ltot-len_
    rr=(k/len_)/(rest_k/rest_len) if k>0 and rest_k>0 else None
    ci=None
    if rr:
        se=math.sqrt(1/k+1/rest_k); ci=(round(rr*math.exp(-1.96*se),3),round(rr*math.exp(1.96*se),3))
    return rr,ci,pv
e4=[]; groups={'virulence panel (946 genes)':[g for g in genes if g in panel],'non-panel genes':[g for g in genes if g not in panel]}
for name,gl in groups.items():
    k=sum(len(sitegene[g]) for g in gl); ln=sum(glen[g] for g in gl)
    rr,ci,pv=rr_test(k,ln,tot_ns,tot_len)
    e4.append({'group':name,'genes':len(gl),'unmasked_bp':ln,'unique_nonsyn_sites':k,'sites_per_kb':round(1000*k/ln,3),'rate_ratio_vs_rest':round(rr,3) if rr else '','CI95':ci,'p':pv,'BH_q':''})
catrows=[]
cats=cog['category_names']
for c in cats:
    gl=[g for g in genes if c in cog['genes'].get(g,{}).get('categories',[])]
    if len(gl)<20: continue
    k=sum(len(sitegene[g]) for g in gl); ln=sum(glen[g] for g in gl)
    rr,ci,pv=rr_test(k,ln,tot_ns,tot_len)
    pan=sum(1 for g in gl if g in panel)
    catrows.append({'COG':c,'name':cats[c],'genes':len(gl),'panel_genes':pan,'unique_nonsyn_sites':k,'sites_per_kb':round(1000*k/ln,3),'rate_ratio_vs_rest':round(rr,3) if rr else '','CI95':ci,'p':pv})
# no-COG genes
gl=[g for g in genes if g not in cog['genes']]; k=sum(len(sitegene[g]) for g in gl); ln=sum(glen[g] for g in gl); rr,ci,pv=rr_test(k,ln,tot_ns,tot_len)
catrows.append({'COG':'none','name':'no COG assignment','genes':len(gl),'panel_genes':sum(1 for g in gl if g in panel),'unique_nonsyn_sites':k,'sites_per_kb':round(1000*k/ln,3),'rate_ratio_vs_rest':round(rr,3) if rr else '','CI95':ci,'p':pv})
ps=np.array([r['p'] for r in catrows]); o=np.argsort(ps); q=np.empty(len(ps)); prev=1
for rank,i in list(enumerate(o))[::-1]:
    prev=min(prev,ps[i]*len(ps)/(rank+1)); q[i]=prev
for r,qq in zip(catrows,q): r['BH_q']=round(float(qq),4)
# panel within COG (panel genes vs non-panel genes of the same category: Mantel-Haenszel-style summary via stratified rate ratio)
mh_num=mh_den=0.0
for c in cats:
    gp=[g for g in genes if c in cog['genes'].get(g,{}).get('categories',[]) and g in panel]; gb=[g for g in genes if c in cog['genes'].get(g,{}).get('categories',[]) and g not in panel]
    if not gp or not gb: continue
    a=sum(len(sitegene[g]) for g in gp); la=sum(glen[g] for g in gp); b=sum(len(sitegene[g]) for g in gb); lb=sum(glen[g] for g in gb); T=la+lb
    mh_num+=a*lb/T; mh_den+=b*la/T
out['e4']={'unique_nonsyn_sites_total':tot_ns,'genome_rate_per_kb':round(1000*bg_rate,3),'panel_vs_background':e4,'categories':catrows,'MH_rate_ratio_panel_vs_background_within_COG':round(mh_num/mh_den,3),
           'excluded_isolate':sorted(EXC),'note':'unique polymorphic sites counted once; shared-by-descent variants are one site, so lineage structure does not multiply the count'}
for fn,rows in (('e4_panel_vs_background',e4),('e4_cog_category_burden',catrows)):
    with open('%s/%s.tsv'%(O,fn),'w') as o_:
        w=csv.DictWriter(o_,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
json.dump(out,open(O+'/c2_diversity_summary.json','w'),indent=1,default=str)
print(json.dumps(out['e4']['panel_vs_background'],default=str)); print('MH',out['e4']['MH_rate_ratio_panel_vs_background_within_COG'])
print([ (r['COG'],r['rate_ratio_vs_rest'],r['BH_q']) for r in catrows if r['BH_q']<0.05])
