"""C3 step 4: collect HyPhy FEL (pervasive) and MEME (episodic) site results for panel genes. Reports sites at p<=0.1 (HyPhy default), p<=0.05, and BH q<0.1 over all tested sites per method."""
import os,json,csv,glob,numpy as np,collections
exec(open(P:='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/code/c3_parse.py').read().split('def parse(rv):')[0]) if False else None
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c3'; C1=P+'/cohort/c1'; L=4411532
G=json.load(open(P+'/ref/genes.json')); cds={rv:g for rv,g in G.items() if g.get('type')=='CDS'}
names={'tier1a_core_virulence':'1a','tier1b_toxin_antitoxin':'1b','tier2_pe_ppe':'2','tier3_exploratory':'3'}
tier={}
for f in glob.glob(P+'/ref/panel/*.bed'):
    t=[v for k,v in names.items() if k in f][0]
    for l in open(f): tier[l.split('\t')[3].split('_')[0]]=t
mask=np.unpackbits(np.load(C1+'/mask.npy'))[:L+1].astype(bool)
def keeplist(rv):
    g=cds[rv]; a,b=g['start'],g['end']; nc=(b-a+1)//3
    m=np.zeros(nc,bool); np.logical_or.at(m,np.arange(b-a+1)//3,mask[a:b+1])
    if g['strand']=='-': m=m[::-1]
    return [k for k in range(nc-1) if not m[k]]
def bh(p):
    p=np.array(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); prev=1.0
    for r,i in enumerate(o[::-1]):
        k=n-r; prev=min(prev,p[i]*n/k); q[i]=prev
    return q
rows=[]; ngenes=0; nsites={'fel':0,'meme':0}
for rv in sorted(tier):
    d='%s/hyphy/%s'%(O,rv)
    if not os.path.exists(d+'/fel.json'): continue
    ngenes+=1; kp=keeplist(rv)
    for m,pi in (('fel',4),('meme',6)):
        if not os.path.exists('%s/%s.json'%(d,m)): continue
        j=json.load(open('%s/%s.json'%(d,m))); cont=j['MLE']['content']['0']
        for i,r in enumerate(cont):
            nsites[m]+=1
            rows.append({'rv':rv,'gene':cds[rv].get('name',''),'tier':tier[rv],'method':m.upper(),'aa_pos':kp[i]+1,'alpha':r[0],'beta':r[1],'p':r[pi],'LRT':r[pi-1],'direction':('positive (beta>alpha)' if r[1]>r[0] else 'negative (beta<alpha)') if m=='fel' else 'episodic'})
for m in ('FEL','MEME'):
    sub=[r for r in rows if r['method']==m and (r['alpha']>0 or r['beta']>0)]; q=bh([r['p'] for r in sub])
    for r,qq in zip(sub,q): r['BH_q_all_panel_sites']=float(qq)
    for r in rows:
        if r['method']==m and 'BH_q_all_panel_sites' not in r: r['BH_q_all_panel_sites']=1.0
sig=[r for r in rows if r['p']<=0.1]
keys=['rv','gene','tier','method','aa_pos','alpha','beta','LRT','p','BH_q_all_panel_sites','direction']
with open(O+'/a19_hyphy_FEL_MEME_sites_p0.1.tsv','w') as o:
    w=csv.DictWriter(o,fieldnames=keys,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(sorted(sig,key=lambda r:(r['method'],r['p'])))
S={'genes_run':ngenes,'sites_tested':nsites}
for m in ('FEL','MEME'):
    allm=[r for r in rows if r['method']==m]; sub=allm
    S[m]={'sites':len(sub),'variable_sites_tested_for_BH':sum(1 for r in allm if r['alpha']>0 or r['beta']>0),'p<=0.1':sum(r['p']<=0.1 for r in sub),'p<=0.05':sum(r['p']<=0.05 for r in sub),'BH_q<0.1':sum(r['BH_q_all_panel_sites']<0.1 for r in sub)}
    if m=='FEL': S[m]['positive_p<=0.1']=sum(r['p']<=0.1 and r['beta']>r['alpha'] for r in sub); S[m]['negative_p<=0.1']=sum(r['p']<=0.1 and r['beta']<r['alpha'] for r in sub)
    S[m]['sites_p<=0.1_list']=[(r['rv'],r['gene'],r['aa_pos'],round(r['p'],5),r.get('direction')) for r in sub if r['p']<=0.1 and (m=='MEME' or r['beta']>r['alpha'])]
    S[m]['BH_q<0.1_sites']=[(r['rv'],r['gene'],r['aa_pos'],round(r['p'],6),round(r['BH_q_all_panel_sites'],4),r.get('direction')) for r in sub if r['BH_q_all_panel_sites']<0.1]
json.dump(S,open(O+'/c3_hyphy_summary.json','w'),indent=1); print(json.dumps(S,indent=1)[:3500])
