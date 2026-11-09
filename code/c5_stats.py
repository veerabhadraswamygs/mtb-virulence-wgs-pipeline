"""C5: statistics 31-36 (diabetes status and other clinical variables).
Pre-specified: Tier 1 only; BH within Tier 1a and within Tier 1b separately (and within analysis family); variants/genes tested only if present in >=3 callable isolates;
Mann-Whitney U with rank-biserial effect and 95% bootstrap CI; Fisher exact with OR (Haldane-Anscombe) and 95% Woolf CI; lineage-stratified analyses by van Elteren / exact stratified (conditional) test."""
import pickle,json,collections,numpy as np,pandas as pd,sys
from scipy import stats
from scipy.stats import hypergeom
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c5'
D=pickle.load(open(O+'/c5_data.pkl','rb')); MISS=pickle.load(open(O+'/c5_missing.pkl','rb'))
S=D['samples']; LIN=D['LIN']; C=D['clin']; V=D['V']; B=D['B']; PRES=D['pres']; INFO=D['info']; LOW=set(D['lowcov'])
rng=np.random.default_rng(20261105)
NS={'missense_variant','stop_gained','start_lost','stop_lost','stop_lost&splice_region_variant','frameshift_variant','disruptive_inframe_deletion','disruptive_inframe_insertion','inframe_deletion','inframe_insertion'}
SNV_NS={'missense_variant','stop_gained','start_lost','stop_lost','stop_lost&splice_region_variant'}
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); prev=1.0
    for i in range(n-1,-1,-1):
        j=o[i]; prev=min(prev,p[j]*n/(i+1)); q[j]=prev
    return q
# ------------------------------------------------------------------ effect sizes
def mwu(x1,x0,nboot=5000):
    x1=np.asarray(x1,float); x0=np.asarray(x0,float)
    u=stats.mannwhitneyu(x1,x0,alternative='two-sided',method='auto')
    n1,n0=len(x1),len(x0); rb=2*u.statistic/(n1*n0)-1
    diffs=(x1[:,None]-x0[None,:]).ravel(); hl=float(np.median(diffs))
    rbs=[];hls=[]
    for _ in range(nboot):
        a=rng.choice(x1,n1); b=rng.choice(x0,n0)
        uu=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic').statistic
        rbs.append(2*uu/(n1*n0)-1); hls.append(np.median((a[:,None]-b[None,:]).ravel()))
    return dict(U=float(u.statistic),p=float(u.pvalue),rank_biserial=rb,rb_lo=float(np.percentile(rbs,2.5)),rb_hi=float(np.percentile(rbs,97.5)),
                prob_superiority=(rb+1)/2,HL_shift=hl,HL_lo=float(np.percentile(hls,2.5)),HL_hi=float(np.percentile(hls,97.5)))
def strat_vanelteren(vals,grp,strata,nperm=10000,nboot=3000):
    """van Elteren stratified rank-sum test, permutation p; effect = weighted mean stratum rank-biserial."""
    vals=np.asarray(vals,float); grp=np.asarray(grp,bool); strata=np.asarray(strata)
    ks=[k for k in np.unique(strata) if grp[strata==k].any() and (~grp[strata==k]).any()]
    if not ks: return None
    idx={k:np.where(strata==k)[0] for k in ks}
    ranks={k:stats.rankdata(vals[idx[k]]) for k in ks}
    def stat(g):
        t=0.0
        for k in ks:
            r=ranks[k]; gg=g[idx[k]]; t+=r[gg].sum()/(len(r)+1)
        return t
    obs=stat(grp); E=sum(grp[idx[k]].sum()*(len(idx[k])+1)/2/(len(idx[k])+1) for k in ks)
    cnt=0
    for _ in range(nperm):
        g=grp.copy()
        for k in ks: g[idx[k]]=rng.permutation(grp[idx[k]])
        if abs(stat(g)-E)>=abs(obs-E)-1e-12: cnt+=1
    p=(cnt+1)/(nperm+1)
    def eff(vv,gg,ss):
        num=den=0.0
        for k in ks:
            m=ss==k; a=vv[m&gg]; b=vv[m&~gg]
            if len(a)==0 or len(b)==0: continue
            u=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic').statistic
            w=len(a)*len(b)/(m.sum()+1); num+=w*(2*u/(len(a)*len(b))-1); den+=w
        return num/den
    est=eff(vals,grp,strata); bs=[]
    for _ in range(nboot):
        sel=np.concatenate([rng.choice(idx[k],len(idx[k])) for k in ks])
        if grp[sel].all() or (~grp[sel]).all(): continue
        try: bs.append(eff(vals[sel],grp[sel],strata[sel]))
        except Exception: pass
    return dict(p=float(p),rank_biserial=float(est),rb_lo=float(np.percentile(bs,2.5)),rb_hi=float(np.percentile(bs,97.5)),strata_used=','.join(ks),n_used=int(sum(len(idx[k]) for k in ks)))
def fisher_row(a,b,c,d):
    """a=exposed carriers, b=exposed non-carriers, c=unexposed carriers, d=unexposed non-carriers"""
    p=stats.fisher_exact([[a,b],[c,d]])[1]
    A,Bq,Cq,Dq=a+.5,b+.5,c+.5,d+.5; orr=A*Dq/(Bq*Cq); se=np.sqrt(1/A+1/Bq+1/Cq+1/Dq)
    return p,orr,orr*np.exp(-1.96*se),orr*np.exp(1.96*se)
def exact_strat(tabs):
    """tabs: list of (a,b,c,d) per stratum. Exact conditional test of common association; returns p, MH OR, RGB CI."""
    pm=np.array([1.0]); obs=0; mh_n=mh_d=0.0; Pq=Rq=Sq=Qq=0.0; used=[]
    for (a,b,c,d) in tabs:
        n=a+b+c+d; n1=a+b; m=a+c
        if n<2 or n1==0 or n1==n or m==0 or m==n: continue
        used.append((a,b,c,d)); obs+=a
        lo=max(0,m-(n-n1)); hi=min(m,n1); pk=hypergeom.pmf(np.arange(lo,hi+1),n,m,n1)
        new=np.zeros(len(pm)+len(pk)-1);
        for i,x in enumerate(pk): new[i:i+len(pm)]+=pm*x
        # shift handled by offset
        pm=new
    if not used: return None
    # recompute with offsets
    offs=0; pm=np.array([1.0])
    for (a,b,c,d) in used:
        n=a+b+c+d; n1=a+b; m=a+c; lo=max(0,m-(n-n1)); hi=min(m,n1); pk=hypergeom.pmf(np.arange(lo,hi+1),n,m,n1)
        new=np.zeros(len(pm)+len(pk)-1)
        for i,x in enumerate(pk): new[i:i+len(pm)]+=pm*x
        pm=new; offs+=lo
    po=pm[obs-offs]; p=float(pm[pm<=po*(1+1e-7)].sum())
    # Mantel-Haenszel OR with RGB variance; 0.5 continuity if any stratum has a zero in a/b/c/d
    zero=any(min(t)==0 for t in used)
    R=Sx=0.0; Pn=Qn=Rn=Sn=0.0
    for (a,b,c,d) in used:
        if zero: a,b,c,d=a+.5,b+.5,c+.5,d+.5
        n=a+b+c+d; r=a*d/n; s=b*c/n; Pk=(a+d)/n; Qk=(b+c)/n
        Rn+=r; Sn+=s; Pn+=Pk*r; Qn+=Qk*s; Sx+=Pk*s+Qk*r
    orr=Rn/Sn if Sn>0 else np.nan
    var=Pn/(2*Rn**2)+Sx/(2*Rn*Sn)+Qn/(2*Sn**2) if Rn>0 and Sn>0 else np.nan
    se=np.sqrt(var); lo=orr*np.exp(-1.96*se); hi=orr*np.exp(1.96*se)
    return dict(p=p,mh_or=orr,mh_lo=lo,mh_hi=hi,n_strata=len(used))
# ------------------------------------------------------------------ groups
def norm_group(s,var,pos): return C[s][var]==pos
def make_burden(sel,tier,defn):
    return np.array([B[(s,tier,defn)]['n']/(B[(s,tier,defn)]['callable_bp']/1000) for s in sel])
out={}
# ================================================================== 31 and 33 for a binary exposure
def burden_tests(var,pos,label,sel=None,sens=''):
    sel=sel or S; rows=[]
    g=np.array([C[s][var]==pos for s in sel]); lin=np.array([LIN[s] for s in sel])
    for tier in('1','1a','1b'):
        for defn in('snv','all'):
            x=make_burden(sel,tier,defn); n=np.array([B[(s,tier,defn)]['n'] for s in sel])
            r=mwu(x[g],x[~g]); st=strat_vanelteren(x,g,lin)
            row=dict(analysis='31/33',variable=label,exposed=pos,sensitivity=sens or 'all 55 isolates',tier='Tier '+tier,variant_definition='nonsynonymous SNV' if defn=='snv' else 'all protein-altering (SNV+indel)',
                     n_exposed=int(g.sum()),n_unexposed=int((~g).sum()),
                     median_burden_per_kb_exposed=float(np.median(x[g])),iqr_exposed='%.2f-%.2f'%tuple(np.percentile(x[g],[25,75])),
                     median_burden_per_kb_unexposed=float(np.median(x[~g])),iqr_unexposed='%.2f-%.2f'%tuple(np.percentile(x[~g],[25,75])),
                     median_count_exposed=float(np.median(n[g])),median_count_unexposed=float(np.median(n[~g])),
                     U=r['U'],p_unadjusted=r['p'],rank_biserial=r['rank_biserial'],rb_ci_lo=r['rb_lo'],rb_ci_hi=r['rb_hi'],prob_superiority=r['prob_superiority'],
                     HL_shift_per_kb=r['HL_shift'],HL_ci_lo=r['HL_lo'],HL_ci_hi=r['HL_hi'])
            if st: row.update(p_lineage_stratified=st['p'],rank_biserial_stratified=st['rank_biserial'],rb_strat_ci_lo=st['rb_lo'],rb_strat_ci_hi=st['rb_hi'],strata_used=st['strata_used'],n_in_informative_strata=st['n_used'])
            rows.append(row)
    return rows
# within-lineage Mann-Whitney (33)
def within_lineage(var,pos,label):
    rows=[]
    for L in sorted(set(LIN.values())):
        sel=[s for s in S if LIN[s]==L]; g=np.array([C[s][var]==pos for s in sel])
        if g.sum()<1 or (~g).sum()<1:
            rows.append(dict(analysis='33',variable=label,lineage=L,n=len(sel),n_exposed=int(g.sum()),n_unexposed=int((~g).sum()),note='not testable (one group absent)')); continue
        x=make_burden(sel,'1','snv');
        u=stats.mannwhitneyu(x[g],x[~g],alternative='two-sided',method='auto')
        rows.append(dict(analysis='33',variable=label,lineage=L,n=len(sel),n_exposed=int(g.sum()),n_unexposed=int((~g).sum()),median_exposed=float(np.median(x[g])),median_unexposed=float(np.median(x[~g])),
                         p=float(u.pvalue),rank_biserial=2*u.statistic/(g.sum()*(~g).sum())-1,note='Tier 1, nonsynonymous SNV burden per kb; n<10 per lineage = very low power' if len(sel)<20 else 'Tier 1, nonsynonymous SNV burden per kb'))
    return rows
# ================================================================== 32 variant / gene tests
def gene_callable(s,rv):
    x=INFO[rv]  # length from per-isolate record below
    return True
GL={}  # gene callable fraction per isolate from catalogue_presence
import json as _j
for s in S:
    R=_j.load(open('%s/work/%s/results.json'%(P,s)))
    GL[s]={x['rv']:(x['callable_bp'],x['length_bp'],x['breadth10']) for x in R['catalogue_presence']}
def gene_pres_call(s,rv):
    c=PRES[s][rv]
    if c=='absent' and s in LOW: return 'not determined'      # single gene-absence rule (v5)
    return c
def assoc_tables(var,pos,label,sel=None,sens=''):
    sel=sel or S; selset=set(sel); g={s:(C[s][var]==pos) for s in sel}
    rows=[]
    # ---- variants
    for k,d in V.items():
        miss=set(MISS[k])
        usable=[s for s in sel if s not in miss]
        carr=[s for s in usable if s in d['isolates']]
        if len(carr)<3: continue
        fam='A protein-altering variant' if d['effect'] in NS else 'D other variant (synonymous / non-coding)'
        a=sum(g[s] for s in carr); c_=len(carr)-a
        n1=sum(g[s] for s in usable); n0=len(usable)-n1
        p,orr,lo,hi=fisher_row(a,n1-a,c_,n0-c_)
        tabs=[]
        for L in sorted(set(LIN[s] for s in usable)):
            u=[s for s in usable if LIN[s]==L]; aa=sum(g[s] for s in u if s in carr); nn1=sum(g[s] for s in u)
            cc=sum(1 for s in u if s in carr and not g[s]); nn0=len(u)-nn1
            tabs.append((aa,nn1-aa,cc,nn0-cc))
        ex=exact_strat(tabs)
        carr_lin=collections.Counter(LIN[s] for s in carr)
        rows.append(dict(analysis='32/33',level='variant',family=fam,variable=label,exposed=pos,sensitivity=sens or 'all 55 isolates',tier='Tier '+d['tier'],rv=d['rv'],gene=d['gene'],
                         pos=d['pos'],ref=d['ref'],alt=d['alt'],effect=d['effect'],hgvs_p=d['hgvs_p'],hgvs_c=d['hgvs_c'],lineage_defining=d['lineage_defining'],
                         n_callable=len(usable),n_exposed_callable=n1,n_unexposed_callable=n0,carriers_exposed=int(a),carriers_unexposed=int(c_),
                         freq_exposed=a/n1 if n1 else np.nan,freq_unexposed=c_/n0 if n0 else np.nan,risk_difference=(a/n1-c_/n0) if n1 and n0 else np.nan,
                         OR=orr,OR_ci_lo=lo,OR_ci_hi=hi,p_fisher=p,
                         carriers_by_lineage=';'.join('%s:%d'%(l,n) for l,n in sorted(carr_lin.items())),lineage_confounded='yes' if len(carr_lin)==1 else ('lineage-defining' if d['lineage_defining'] else 'no'),
                         p_lineage_stratified=ex['p'] if ex else np.nan,OR_MH=ex['mh_or'] if ex else np.nan,OR_MH_ci_lo=ex['mh_lo'] if ex else np.nan,OR_MH_ci_hi=ex['mh_hi'] if ex else np.nan,n_informative_strata=ex['n_strata'] if ex else 0))
    # ---- gene-level carriage of protein-altering variants
    gens=collections.defaultdict(dict)
    for rv in D['t1']:
        ok=[s for s in sel if GL[s][rv][0]/GL[s][rv][1]>=0.9 and gene_pres_call(s,rv)!='absent']
        carr=[s for s in ok if any((s in d['isolates']) and d['rv']==rv and d['effect'] in NS for d in [])]
        gens[rv]['ok']=ok
    gcarr=collections.defaultdict(set)
    for k,d in V.items():
        if d['effect'] in NS:
            for s in d['isolates']: gcarr[d['rv']].add(s)
    for rv,dd in gens.items():
        ok=dd['ok']; carr=[s for s in ok if s in gcarr[rv]]
        if len(carr)<3: continue
        a=sum(g[s] for s in carr); c_=len(carr)-a; n1=sum(g[s] for s in ok); n0=len(ok)-n1
        p,orr,lo,hi=fisher_row(a,n1-a,c_,n0-c_)
        tabs=[]
        for L in sorted(set(LIN[s] for s in ok)):
            u=[s for s in ok if LIN[s]==L]; aa=sum(g[s] for s in u if s in carr); nn1=sum(g[s] for s in u)
            cc=sum(1 for s in u if s in carr and not g[s]); tabs.append((aa,nn1-aa,cc,len(u)-nn1-cc))
        ex=exact_strat(tabs); cl=collections.Counter(LIN[s] for s in carr)
        rows.append(dict(analysis='32/33',level='gene (any protein-altering variant)',family='B gene-level carriage',variable=label,exposed=pos,sensitivity=sens or 'all 55 isolates',tier='Tier '+D['t1'][rv],rv=rv,gene=INFO[rv]['gene'],
                         n_callable=len(ok),n_exposed_callable=n1,n_unexposed_callable=n0,carriers_exposed=int(a),carriers_unexposed=int(c_),freq_exposed=a/n1,freq_unexposed=c_/n0,risk_difference=a/n1-c_/n0,
                         OR=orr,OR_ci_lo=lo,OR_ci_hi=hi,p_fisher=p,carriers_by_lineage=';'.join('%s:%d'%(l,n) for l,n in sorted(cl.items())),lineage_confounded='yes' if len(cl)==1 else 'no',
                         p_lineage_stratified=ex['p'] if ex else np.nan,OR_MH=ex['mh_or'] if ex else np.nan,OR_MH_ci_lo=ex['mh_lo'] if ex else np.nan,OR_MH_ci_hi=ex['mh_hi'] if ex else np.nan,n_informative_strata=ex['n_strata'] if ex else 0))
    # ---- presence/absence (>=3 isolates absent among determined)
    for rv in D['t1']:
        det=[s for s in sel if gene_pres_call(s,rv)!='not determined']; ab=[s for s in det if gene_pres_call(s,rv)=='absent']
        if len(ab)<3: continue
        a=sum(g[s] for s in ab); c_=len(ab)-a; n1=sum(g[s] for s in det); n0=len(det)-n1
        p,orr,lo,hi=fisher_row(a,n1-a,c_,n0-c_)
        tabs=[]
        for L in sorted(set(LIN[s] for s in det)):
            u=[s for s in det if LIN[s]==L]; aa=sum(g[s] for s in u if s in ab); nn1=sum(g[s] for s in u); cc=sum(1 for s in u if s in ab and not g[s]); tabs.append((aa,nn1-aa,cc,len(u)-nn1-cc))
        ex=exact_strat(tabs); cl=collections.Counter(LIN[s] for s in ab)
        rows.append(dict(analysis='32/33',level='gene absence',family='C gene presence/absence',variable=label,exposed=pos,sensitivity=sens or 'all 55 isolates',tier='Tier '+D['t1'][rv],rv=rv,gene=INFO[rv]['gene'],
                         n_callable=len(det),n_exposed_callable=n1,n_unexposed_callable=n0,carriers_exposed=int(a),carriers_unexposed=int(c_),freq_exposed=a/n1,freq_unexposed=c_/n0,risk_difference=a/n1-c_/n0,
                         OR=orr,OR_ci_lo=lo,OR_ci_hi=hi,p_fisher=p,carriers_by_lineage=';'.join('%s:%d'%(l,n) for l,n in sorted(cl.items())),lineage_confounded='yes' if len(cl)==1 else 'no',
                         p_lineage_stratified=ex['p'] if ex else np.nan,OR_MH=ex['mh_or'] if ex else np.nan,OR_MH_ci_lo=ex['mh_lo'] if ex else np.nan,OR_MH_ci_hi=ex['mh_hi'] if ex else np.nan,n_informative_strata=ex['n_strata'] if ex else 0))
    df=pd.DataFrame(rows)
    if len(df):
        df['q_fisher_BH']=np.nan; df['q_lineage_stratified_BH']=np.nan; df['n_tests_in_stratum']=0
        for (tier,fam),idx in df.groupby(['tier','family']).groups.items():
            df.loc[idx,'q_fisher_BH']=bh(df.loc[idx,'p_fisher'].values); df.loc[idx,'n_tests_in_stratum']=len(idx)
            pp=df.loc[idx,'p_lineage_stratified']; m=pp.notna()
            if m.any(): df.loc[pp[m].index,'q_lineage_stratified_BH']=bh(pp[m].values)
    return df
# ================================================================== run diabetes (31-33)
res={}
burden=burden_tests('diabetes','Yes','Diabetes mellitus')
burden_sens=burden_tests('diabetes','Yes','Diabetes mellitus',sel=[s for s in S if s not in LOW],sens='excluding 4 low-coverage isolates (MTB4, MTB25, MTB32, MTB35)')
wl=within_lineage('diabetes','Yes','Diabetes mellitus')
print('31/33 done'); sys.stdout.flush()
A=assoc_tables('diabetes','Yes','Diabetes mellitus')
As=assoc_tables('diabetes','Yes','Diabetes mellitus',sel=[s for s in S if s not in LOW],sens='excluding 4 low-coverage isolates (MTB4, MTB25, MTB32, MTB35)')
print('32 done',len(A),len(As)); sys.stdout.flush()
pd.DataFrame(burden+burden_sens).to_csv(O+'/a31_a33_burden_diabetes.tsv',sep='\t',index=False)
pd.DataFrame(wl).to_csv(O+'/a33_burden_within_lineage_diabetes.tsv',sep='\t',index=False)
A.to_csv(O+'/a32_a33_variant_gene_tests_diabetes.tsv',sep='\t',index=False)
As.to_csv(O+'/a32_sensitivity_excl_lowcov_diabetes.tsv',sep='\t',index=False)
# lineage x diabetes
ct=pd.crosstab(pd.Series([LIN[s] for s in S],name='lineage'),pd.Series([C[s]['diabetes'] for s in S],name='diabetes'))
ct.to_csv(O+'/a33_lineage_by_diabetes_crosstab.tsv',sep='\t')
# permutation Fisher-Freeman-Halton style test lineage vs diabetes
lin=np.array([LIN[s] for s in S]); dm=np.array([C[s]['diabetes']=='Yes' for s in S])
def chi(l,d):
    t=pd.crosstab(l,d).values.astype(float); e=t.sum(1)[:,None]*t.sum(0)[None,:]/t.sum(); return ((t-e)**2/e).sum()
o=chi(lin,dm); cnt=sum(chi(lin,rng.permutation(dm))>=o-1e-9 for _ in range(10000)); lin_p=(cnt+1)/10001
json.dump({'lineage_vs_diabetes_chi2':o,'permutation_p':lin_p},open(O+'/a33_lineage_vs_diabetes_test.json','w'))
print('lineage~diabetes perm p',lin_p)
# ================================================================== 34 epitope partition
rows=[]
def a34(sel,sens):
    out=[]
    vals=[];grp=[];ln=[]
    for s in sel:
        r=D['epi'][s]['ratio_epitope_to_non_epitope_pi']
        if r is None: continue
        vals.append(r);grp.append(C[s]['diabetes']=='Yes');ln.append(LIN[s])
    vals=np.array(vals,float); grp=np.array(grp); ln=np.array(ln)
    def run(stat_name,f):
        obs=f(vals[grp])-f(vals[~grp]); cnt=0
        for _ in range(10000):
            g=grp.copy()
            for L in np.unique(ln):
                ix=np.where(ln==L)[0]; g[ix]=rng.permutation(grp[ix])
            if abs(f(vals[g])-f(vals[~g]))>=abs(obs)-1e-12: cnt+=1
        # bootstrap CI for difference
        bs=[f(rng.choice(vals[grp],grp.sum()))-f(rng.choice(vals[~grp],(~grp).sum())) for _ in range(5000)]
        return dict(analysis='34',sensitivity=sens,statistic=stat_name,n_exposed=int(grp.sum()),n_unexposed=int((~grp).sum()),value_exposed=float(f(vals[grp])),value_unexposed=float(f(vals[~grp])),
                    difference=float(obs),diff_ci_lo=float(np.percentile(bs,2.5)),diff_ci_hi=float(np.percentile(bs,97.5)),permutation_p_within_lineage=(cnt+1)/10001,permutations=10000)
    out.append(run('mean epitope:non-epitope pi ratio',np.mean)); out.append(run('median epitope:non-epitope pi ratio',np.median))
    st=strat_vanelteren(vals,grp,ln); r=mwu(vals[grp],vals[~grp])
    out.append(dict(analysis='34',sensitivity=sens,statistic='rank-based: Mann-Whitney (unadjusted) / van Elteren (lineage-stratified)',n_exposed=int(grp.sum()),n_unexposed=int((~grp).sum()),
                    difference=r['rank_biserial'],diff_ci_lo=r['rb_lo'],diff_ci_hi=r['rb_hi'],permutation_p_within_lineage=st['p'] if st else np.nan,mwu_p_unadjusted=r['p'],value_exposed=np.nan,value_unexposed=np.nan,permutations=10000))
    return out
rows+=a34(S,'all 55 isolates with a defined ratio')
rows+=a34([s for s in S if s not in LOW],'excluding 4 low-coverage isolates')
# secondary per-isolate quantities
for nm,key,sub in [('epitope pi per kb','pi_per_kb','epitope'),('non-epitope pi per kb','pi_per_kb','non_epitope')]:
    x=np.array([D['epi'][s][sub][key] for s in S],float); g=np.array([C[s]['diabetes']=='Yes' for s in S]); r=mwu(x[g],x[~g]); st=strat_vanelteren(x,g,np.array([LIN[s] for s in S]))
    rows.append(dict(analysis='34',sensitivity='all 55 isolates (secondary)',statistic=nm,n_exposed=int(g.sum()),n_unexposed=int((~g).sum()),value_exposed=float(np.median(x[g])),value_unexposed=float(np.median(x[~g])),
                     difference=r['rank_biserial'],diff_ci_lo=r['rb_lo'],diff_ci_hi=r['rb_hi'],mwu_p_unadjusted=r['p'],permutation_p_within_lineage=st['p'] if st else np.nan,permutations=10000))
pd.DataFrame(rows).to_csv(O+'/a34_epitope_partition_diabetes.tsv',sep='\t',index=False)
print('34 done'); sys.stdout.flush()
# ================================================================== 35 other clinical variables
audit=[]
CAND=[('sex','Male','Sex (male vs female)'),('smoking','Yes','Smoking (ADDITIVE HABITS includes Smoker)'),('alcohol','Yes','Alcohol (ADDITIVE HABITS includes Alcoholic)'),('ses','BPL','Socioeconomic status (BPL vs APL)'),
      ('smear_3plus','Yes','Sputum smear 3+ vs lower'),('prev_tb','Yes','Previous TB'),('died','Yes','Clinical outcome: dead'),('cough','Yes','Persistent cough'),('fever','Yes','Fever'),
      ('weight_loss','Yes','Weight loss'),('chest_pain','Yes','Chest pain'),('breathlessness','Yes','Breathlessness'),('haemoptysis','Yes','Haemoptysis')]
tested=[]
for var,pos,label in CAND:
    n1=sum(C[s][var]==pos for s in S); n0=len(S)-n1
    ok=min(n1,n0)>=5
    audit.append(dict(variable=label,column=var,n_exposed=n1,n_unexposed=n0,tested='yes' if ok else 'no',reason='' if ok else 'minority group < 5 isolates (not testable)'))
    if ok: tested.append((var,pos,label))
age=np.array([C[s]['age'] for s in S]); audit.append(dict(variable='Age (continuous)',column='age',n_exposed=len(S),n_unexposed='',tested='yes (Spearman with burden)',reason=''))
pd.DataFrame(audit).to_csv(O+'/a35_metadata_audit_variables.tsv',sep='\t',index=False)
b35=[];w35=[];v35=[]
for var,pos,label in tested:
    b35+=burden_tests(var,pos,label); w35+=within_lineage(var,pos,label)
    df=assoc_tables(var,pos,label); v35.append(df)
V35=pd.concat(v35,ignore_index=True)
pd.DataFrame(b35).to_csv(O+'/a35_burden_other_clinical.tsv',sep='\t',index=False)
pd.DataFrame(w35).to_csv(O+'/a35_burden_within_lineage_other_clinical.tsv',sep='\t',index=False)
V35.to_csv(O+'/a35_variant_gene_tests_other_clinical.tsv',sep='\t',index=False)
x=make_burden(S,'1','snv'); rho=stats.spearmanr(age,x)
# bootstrap CI for rho
bs=[]
for _ in range(3000):
    ix=rng.integers(0,len(S),len(S)); bs.append(stats.spearmanr(age[ix],x[ix])[0])
pd.DataFrame([dict(analysis='35',variable='Age (continuous)',outcome='Tier 1 nonsynonymous SNV burden per kb',n=len(S),spearman_rho=rho[0],ci_lo=np.percentile(bs,2.5),ci_hi=np.percentile(bs,97.5),p=rho[1])]).to_csv(O+'/a35_age_burden_spearman.tsv',sep='\t',index=False)
age_dm=mwu(age[dm],age[~dm])
print('35 done'); sys.stdout.flush()
# ================================================================== 36 multivariable
events=int(dm.sum()); nonev=int((~dm).sum())
lin_levels=sorted(set(LIN.values())); npar=1+1+(len(lin_levels)-1)+1   # age, sex, lineage dummies, smoking
epv=min(events,nonev)/npar
desc=[]
for nm,f in [('age (median, IQR)',lambda g:'%.0f (%.0f-%.0f)'%(np.median(age[g]),*np.percentile(age[g],[25,75]))),
             ('male n (%)',lambda g:'%d (%.0f%%)'%(sum(C[s]['sex']=='Male' for s,gg in zip(S,g) if gg),100*sum(C[s]['sex']=='Male' for s,gg in zip(S,g) if gg)/g.sum())),
             ('smoker n (%)',lambda g:'%d (%.0f%%)'%(sum(C[s]['smoking']=='Yes' for s,gg in zip(S,g) if gg),100*sum(C[s]['smoking']=='Yes' for s,gg in zip(S,g) if gg)/g.sum()))]+\
            [('lineage %s n (%%)'%L,(lambda L:lambda g:'%d (%.0f%%)'%(sum(LIN[s]==L for s,gg in zip(S,g) if gg),100*sum(LIN[s]==L for s,gg in zip(S,g) if gg)/g.sum()))(L)) for L in lin_levels]:
    desc.append(dict(characteristic=nm,diabetic=f(dm),non_diabetic=f(~dm)))
sx=[C[s]['sex']=='Male' for s in S]; sm=[C[s]['smoking']=='Yes' for s in S]
fs=fisher_row(sum(a and b for a,b in zip(dm,sx)),sum((not a) and False for a,b in zip(dm,sx)) if False else int(dm.sum())-sum(a and b for a,b in zip(dm,sx)),sum((not a) and b for a,b in zip(dm,sx)),int((~dm).sum())-sum((not a) and b for a,b in zip(dm,sx)))
fm=fisher_row(sum(a and b for a,b in zip(dm,sm)),int(dm.sum())-sum(a and b for a,b in zip(dm,sm)),sum((not a) and b for a,b in zip(dm,sm)),int((~dm).sum())-sum((not a) and b for a,b in zip(dm,sm)))
res36=dict(analysis=36,outcome='diabetes mellitus',events=events,non_events=nonev,candidate_parameters=npar,parameters='age, sex, lineage (%s as %d levels -> %d dummies), smoking'%('/'.join(lin_levels),len(lin_levels),len(lin_levels)-1),
           events_per_variable=round(epv,2),threshold=10,condition_met=bool(epv>=10),decision='Condition not met (EPV %.1f < 10): logistic regression NOT fitted; descriptive results only'%epv if epv<10 else 'fitted',
           age_MWU_p=age_dm['p'],age_rank_biserial=age_dm['rank_biserial'],age_rb_ci=[age_dm['rb_lo'],age_dm['rb_hi']],sex_male_OR=fs[1],sex_OR_ci=[fs[2],fs[3]],sex_fisher_p=fs[0],smoking_OR=fm[1],smoking_OR_ci=[fm[2],fm[3]],smoking_fisher_p=fm[0],
           lineage_vs_diabetes_permutation_p=lin_p,descriptives=desc)
# smallest-EPV statement for other binary outcomes
res36['other_outcomes_EPV']={lab:round(min(sum(C[s][v]==p for s in S),sum(C[s][v]!=p for s in S))/npar,2) for v,p,lab in tested}
json.dump(res36,open(O+'/a36_multivariable_decision.json','w'),indent=1,default=float)
print(json.dumps({k:res36[k] for k in('events','non_events','candidate_parameters','events_per_variable','decision')}))
# ---------- console summary
for nm,df in[('A (diabetes)',A),('A sens',As)]:
    print(nm,df.groupby(['tier','family']).size().to_dict())
    print(' min q_fisher',df.groupby(['tier','family'])['q_fisher_BH'].min().round(3).to_dict())
    print(' min q_strat',df.groupby(['tier','family'])['q_lineage_stratified_BH'].min().round(3).to_dict())
print(pd.DataFrame(burden)[['tier','variant_definition','p_unadjusted','rank_biserial','rb_ci_lo','rb_ci_hi','p_lineage_stratified','rank_biserial_stratified']].round(3).to_string())
