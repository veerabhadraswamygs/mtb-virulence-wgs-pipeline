"""C3 step 3: parse codeml output -> analysis 19 (M1a vs M2a LRT, BH within tier) and E15 (genome-wide dN/dS)."""
import os,re,json,csv,glob,math,collections
import numpy as np
from scipy import stats
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; C1=P+'/cohort/c1'; O=P+'/cohort/c3'; L=4411532
G=json.load(open(P+'/ref/genes.json')); cds={rv:g for rv,g in G.items() if g.get('type')=='CDS'}
tier={}
names={'tier1a_core_virulence':'1a','tier1b_toxin_antitoxin':'1b','tier2_pe_ppe':'2','tier3_exploratory':'3'}
for f in glob.glob(P+'/ref/panel/*.bed'):
    t=[v for k,v in names.items() if k in f][0]
    for l in open(f): tier[l.split('\t')[3].split('_')[0]]=t
mask=np.unpackbits(np.load(C1+'/mask.npy'))[:L+1].astype(bool)
summ={r['rv']:r for r in csv.DictReader(open(O+'/c3_gene_alignment_summary.tsv'),delimiter='\t')}
def keeplist(rv):
    g=cds[rv]; a,b=g['start'],g['end']; nc=(b-a+1)//3
    m=np.zeros(nc,bool); np.logical_or.at(m,np.arange(b-a+1)//3,mask[a:b+1])
    if g['strand']=='-': m=m[::-1]
    return [k for k in range(nc-1) if not m[k]]
def parse(rv):
    f='%s/cml/%s/out.txt'%(O,rv)
    if not os.path.exists(f): return None
    t=open(f).read()
    if 'lnL(ntime' not in t: return None
    secs=re.split(r'\nNSsites Model (\d):',t)
    res={}
    ln=re.findall(r'lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\):\s+(-?[\d.]+)',t)
    secd={}
    for i in range(1,len(secs),2): secd[int(secs[i])]=secs[i+1]
    if not secd: secd[0]=t
    out={}
    for k,s in secd.items():
        m=re.search(r'lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\):\s+(-?[\d.]+)',t[t.find('Model %d:'%k) if False else 0:])
    # lnL lines appear before each section header: take them in order
    lnls=[(int(a),int(b),float(c)) for a,b,c in ln]
    ks=sorted(secd)
    # find lnL positions preceding each section
    pos=[m.start() for m in re.finditer(r'lnL\(ntime',t)]
    for i,k in enumerate(ks): out[k]={'lnL':lnls[i][2],'np':lnls[i][1]}
    s0=secd[0]
    m=re.search(r'omega \(dN/dS\) =\s+([\d.]+)',s0); out[0]['omega']=float(m.group(1))
    m=re.search(r'kappa \(ts/tv\) =\s+([\d.]+)',s0); out[0]['kappa']=float(m.group(1)) if m else None
    m=re.search(r'tree length for dN:\s+([\d.]+)',s0); out[0]['tl_dN']=float(m.group(1)) if m else None
    m=re.search(r'tree length for dS:\s+([\d.]+)',s0); out[0]['tl_dS']=float(m.group(1)) if m else None
    nn=ss=0.0; Ns=[];Ss=[]
    tab=s0.split('dN & dS for each branch')[1].split('\n')
    for l in tab:
        mm=re.match(r'\s+\d+\.\.\d+\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)',l)
        if mm:
            Ns.append(float(mm.group(2))); Ss.append(float(mm.group(3))); nn+=float(mm.group(7)); ss+=float(mm.group(8))
        elif l.strip().startswith('tree length') or l.strip().startswith('w ratios'): break
    out[0].update({'N_sites':Ns[0] if Ns else None,'S_sites':Ss[0] if Ss else None,'N_subs':nn,'S_subs':ss})
    if 2 in secd:
        s2=secd[2]
        m=re.search(r'\np:\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\nw:\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)',s2)
        if m: out[2]['p2']=float(m.group(3)); out[2]['w2']=float(m.group(6))
        beb=[]
        if 'Bayes Empirical Bayes' in s2:
            blk=s2.split('Bayes Empirical Bayes')[1]
            for l in blk.split('\n'):
                mm=re.match(r'\s+(\d+)\s+([A-Z\*])\s+([\d.]+)(\**)',l)
                if mm and float(mm.group(3))>0.95: beb.append((int(mm.group(1)),mm.group(2),float(mm.group(3))))
        out[2]['beb']=beb
    return out
def bh(p):
    p=np.array(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); prev=1.0
    for r,i in enumerate(o[::-1]):
        k=n-r; prev=min(prev,p[i]*n/k); q[i]=prev
    return q
rows=[];e15=[]
for rv,s in summ.items():
    if s['status']!='ok' and s['status']!='invariant': continue
    base={'rv':rv,'gene':s['gene'],'tier':tier.get(rv,'non-panel'),'codons':s['codons_kept'],'n_isolates':s['n_usable'],'variable_codons':s['variable_codons'],'syn_alleles':s['alt_codon_alleles_syn'],'nonsyn_alleles':s['alt_codon_alleles_nonsyn']}
    if s['status']=='invariant': e15.append({**base,'status':'invariant (no polymorphic codon)','omega_M0':'','N_subs':0,'S_subs':0,'N_sites':'','S_sites':''}); continue
    o=parse(rv)
    if o is None: 
        sk=os.path.exists('%s/aln_c/%s.skip'%(O,rv))
        e15.append({**base,'status':'<4 distinct haplotypes' if sk else 'not run/failed','omega_M0':'','N_subs':'','S_subs':'','N_sites':'','S_sites':''}); 
        if rv in tier: rows.append({**base,'status':'untestable (<4 distinct haplotypes)' if sk else 'not run/failed'})
        continue
    m0=o[0]; om=m0['omega']
    e15.append({**base,'status':'ok','omega_M0':round(om,4) if m0['tl_dS'] and m0['tl_dS']>0.0 else 'NA (dS=0)','tree_len_dN':m0['tl_dN'],'tree_len_dS':m0['tl_dS'],'N_subs':round(m0['N_subs'],2),'S_subs':round(m0['S_subs'],2),'N_sites':m0['N_sites'],'S_sites':m0['S_sites']})
    if rv in tier and 2 in o:
        l1=o[1]['lnL'];l2=o[2]['lnL']; lrt=max(0.0,2*(l2-l1)); p=float(stats.chi2.sf(lrt,2))
        kp=keeplist(rv); st=cds[rv]['strand']
        beb=[(kp[i-1]+1,a,pr) for i,a,pr in o[2].get('beb',[])]
        rows.append({**base,'status':'tested','lnL_M1a':l1,'lnL_M2a':l2,'LRT_2dlnL':round(lrt,4),'p':p,'M2a_p_positive_class':o[2].get('p2',''),'M2a_omega_positive':o[2].get('w2',''),'omega_M0':round(om,4) if m0['tl_dS'] and m0['tl_dS']>0 else 'NA (dS=0)','BEB_sites_P>0.95':';'.join('%s%d(%.3f)'%(a,i,pr) for i,a,pr in beb)})
# BH within tier
for t in ('1a','1b','2','3'):
    r=[x for x in rows if x['tier']==t and x['status']=='tested']
    q=bh([x['p'] for x in r]) if r else []
    for x,qq in zip(r,q): x['BH_q_within_tier']=float(qq)
keys=['rv','gene','tier','codons','n_isolates','variable_codons','syn_alleles','nonsyn_alleles','status','omega_M0','lnL_M1a','lnL_M2a','LRT_2dlnL','p','BH_q_within_tier','M2a_p_positive_class','M2a_omega_positive','BEB_sites_P>0.95']
with open(O+'/a19_codeml_M1a_vs_M2a_panel.tsv','w') as o:
    w=csv.DictWriter(o,fieldnames=keys,delimiter='\t',restval=''); w.writeheader(); w.writerows(sorted(rows,key=lambda x:(x['tier'],x.get('p',2) if x['status']=='tested' else 2)))
k2=['rv','gene','tier','codons','n_isolates','variable_codons','syn_alleles','nonsyn_alleles','status','omega_M0','tree_len_dN','tree_len_dS','N_subs','S_subs','N_sites','S_sites']
with open(O+'/e15_genomewide_dNdS_per_gene.tsv','w') as o:
    w=csv.DictWriter(o,fieldnames=k2,delimiter='\t',restval=''); w.writeheader(); w.writerows(e15)
# summaries
S={}
for t in ('1a','1b','2','3'):
    r=[x for x in rows if x['tier']==t]; te=[x for x in r if x['status']=='tested']
    S['tier'+t]={'panel_genes_total_in_tier':sum(1 for g,v in tier.items() if v==t),'genes_with_alignment_rows':len(r),'tested':len(te),'p<0.05_uncorrected':sum(1 for x in te if x['p']<0.05),'BH_q<0.05':sum(1 for x in te if x['BH_q_within_tier']<0.05),'BH_q<0.10':sum(1 for x in te if x['BH_q_within_tier']<0.10),
      'significant_genes':[(x['rv'],x['gene'],round(x['p'],6),round(x['BH_q_within_tier'],4),x['BEB_sites_P>0.95']) for x in te if x['BH_q_within_tier']<0.05]}
# E15 pooled
rng=np.random.default_rng(20261103)
def pooled(sel):
    N=np.array([float(x['N_sites']) for x in sel]); Sx=np.array([float(x['S_sites']) for x in sel]); nn=np.array([float(x['N_subs']) for x in sel]); ss=np.array([float(x['S_subs']) for x in sel])
    return N,Sx,nn,ss
def om(N,Sx,nn,ss): 
    return (nn.sum()/N.sum())/(ss.sum()/Sx.sum()) if ss.sum()>0 else float('nan')
ok=[x for x in e15 if x['status']=='ok']
inv=[x for x in e15 if x['status'].startswith('invariant')]
groups={'panel':[x for x in ok if x['tier']!='non-panel'],'non-panel':[x for x in ok if x['tier']=='non-panel']}
for t in ('1a','1b','2','3'): groups['tier'+t]=[x for x in ok if x['tier']==t]
E={}
for k,sel in groups.items():
    if len(sel)<5: continue
    N,Sx,nn,ss=pooled(sel); est=om(N,Sx,nn,ss); bs=[]
    for _ in range(2000):
        i=rng.integers(0,len(sel),len(sel)); bs.append(om(N[i],Sx[i],nn[i],ss[i]))
    E[k]={'genes_polymorphic':len(sel),'N_subs':float(nn.sum()),'S_subs':float(ss.sum()),'pooled_dNdS':round(est,4),'CI95_gene_bootstrap':(round(float(np.nanpercentile(bs,2.5)),4),round(float(np.nanpercentile(bs,97.5)),4))}
# genes with omega defined
def wvals(sel): return np.array([float(x['omega_M0']) for x in sel if x['omega_M0'] not in ('','NA (dS=0)') and float(x['omega_M0'])<50])
wp,wb=wvals(groups['panel']),wvals(groups['non-panel'])
u,p=stats.mannwhitneyu(wp,wb,alternative='two-sided')
E['per_gene_omega_compare']={'n_panel':len(wp),'n_nonpanel':len(wb),'median_panel':round(float(np.median(wp)),4),'median_nonpanel':round(float(np.median(wb)),4),'mannwhitney_p':float(p),'rank_biserial':round(float(2*u/(len(wp)*len(wb))-1),3),'note':'genes with dS=0 or omega>=50 excluded (undefined/saturated)'}
# panel vs non-panel difference in pooled omega: permutation of gene labels
sel=groups['panel']+groups['non-panel']; N,Sx,nn,ss=pooled(sel); lab=np.array([1]*len(groups['panel'])+[0]*len(groups['non-panel']))
obs=om(N[lab==1],Sx[lab==1],nn[lab==1],ss[lab==1])-om(N[lab==0],Sx[lab==0],nn[lab==0],ss[lab==0]); cnt=0
for _ in range(5000):
    pl=rng.permutation(lab); d=om(N[pl==1],Sx[pl==1],nn[pl==1],ss[pl==1])-om(N[pl==0],Sx[pl==0],nn[pl==0],ss[pl==0]); cnt+= abs(d)>=abs(obs)
E['panel_minus_nonpanel_pooled_omega']={'difference':round(obs,4),'permutation_p_5000':(cnt+1)/5001}
S['E15']=E; S['E15_counts']={'genes_total_cds':len(cds),'with_alignment_status':collections.Counter(x['status'] for x in e15)}
json.dump(S,open(O+'/c3_codeml_summary.json','w'),indent=1,default=str)
print(json.dumps(S,indent=1,default=str)[:6000])
