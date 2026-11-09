"""C5 deliverable 7: publication figures + source data (TbD1 figure was made in C4)."""
import json,os,re,collections,pickle,numpy as np,pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,Rectangle
from matplotlib.lines import Line2D
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; C=P+'/cohort'
OUT='/mnt/user-data/outputs/FINAL_DELIVERABLES/7_figures'; os.makedirs(OUT,exist_ok=True)
D=pickle.load(open(C+'/c5/c5_data.pkl','rb')); S=D['samples']; LIN=D['LIN']; CL=D['clin']; LOW=set(D['lowcov'])
INK='#0b0b0b'; INK2='#52514e'; MUTED='#898781'; GRID='#e1e0d9'
LC={'L1':'#2a78d6','L2':'#eb6834','L3':'#1baf7a','L4':'#6250d6','unassigned':'#898781'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.edgecolor':'#c3c2b7','axes.labelcolor':INK2,'xtick.color':INK2,'ytick.color':INK2,'text.color':INK,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
    fig.savefig('%s/%s.png'%(OUT,name),dpi=300,bbox_inches='tight',facecolor='white'); fig.savefig('%s/%s.svg'%(OUT,name),bbox_inches='tight',facecolor='white'); plt.close(fig)
# ---------------- newick
def parse_newick(t):
    t=t.strip().rstrip(';'); pos=0
    def node():
        nonlocal pos
        if t[pos]=='(':
            pos+=1; ch=[]
            while True:
                ch.append(node())
                if t[pos]==',': pos+=1; continue
                if t[pos]==')': pos+=1; break
            m=re.match(r'([^:,();]*)(?::([0-9.eE+-]+))?',t[pos:]); pos+=m.end()
            return dict(children=ch,support=float(m.group(1)) if m.group(1) else None,bl=float(m.group(2) or 0),name=None)
        m=re.match(r'([^:,();]+)(?::([0-9.eE+-]+))?',t[pos:]); pos+=m.end()
        return dict(children=[],name=m.group(1),bl=float(m.group(2) or 0),support=None)
    return node()
def drop(n,name):
    n['children']=[drop(c,name) for c in n['children'] if c['name']!=name]
    return n
def collapse_unary(n):
    n['children']=[collapse_unary(c) for c in n['children']]
    if len(n['children'])==1 and n['name'] is None:
        c=n['children'][0]; c['bl']+=n['bl']; return c
    return n
def nleaves(n): return 1 if not n['children'] else sum(nleaves(c) for c in n['children'])
def ladderize(n):
    n['children']=sorted([ladderize(c) for c in n['children']],key=nleaves); return n
def layout(n,x=0,ctr=[0]):
    n['x']=x+n['bl']
    if not n['children']: n['y']=ctr[0]; ctr[0]+=1
    else:
        for c in n['children']: layout(c,n['x'],ctr)
        n['y']=np.mean([c['y'] for c in n['children']])
    return n
def leaves(n): return [n] if not n['children'] else [l for c in n['children'] for l in leaves(c)]
def draw(n,ax,rows):
    for c in n['children']:
        ax.plot([n['x'],n['x']],[n['y'],c['y']],color=INK2,lw=0.6,solid_capstyle='butt'); ax.plot([n['x'],c['x']],[c['y'],c['y']],color=INK2,lw=0.6,solid_capstyle='butt')
        draw(c,ax,rows)
    if n['children'] and n.get('support') is not None:
        rows.append((n['x'],n['y'],n['support']))
def load_tree(path,dropname='Mcanettii'):
    t=parse_newick(open(path).read()); t=collapse_unary(drop(t,dropname)); t=ladderize(t); layout(t,0,[0]); return t
def load_tree2(path):
    t=parse_newick(open(path).read()); t=collapse_unary(drop(t,'Mcanettii')); 
    big=max(t['children'],key=nleaves); big['bl']=0.00012     # long stem to the ingroup is drawn shortened (axis break)
    t=ladderize(t); layout(t,0,[0]); return t,big
tree,bigclade=load_tree2(C+'/c1/tree55.treefile'); tips=leaves(tree); order=[l['name'] for l in sorted(tips,key=lambda l:l['y'])]
# =================================================================== F2 annotated tree
res_dm={}
for s in S:
    R=json.load(open('%s/work/%s/results.json'%(P,s))); res_dm[s]=sorted({r['drug'] for r in R['resistance'] if r.get('confidence') in('Assoc w R','Assoc w R - interim')})
fig=plt.figure(figsize=(7.2,9.2)); ax=fig.add_axes([0.04,0.04,0.5,0.93]); rows=[]; draw(tree,ax,rows)
xmax=max(l['x'] for l in tips)
for l in tips: ax.plot([l['x'],xmax],[l['y'],l['y']],color=GRID,lw=0.4,ls=(0,(1,2)),zorder=0)
sup=[(x,y) for x,y,v in rows if v>=95]; ax.scatter(*zip(*sup),s=6,color='#0b0b0b',zorder=3)
ax.plot([0.00005,0.00007],[bigclade['y']-0.9,bigclade['y']+0.9],color='white',lw=3); ax.text(0.00006,bigclade['y']-1.2,'//',ha='center',fontsize=9); ax.set_ylim(len(tips)-0.3,-3.0); ax.set_xlim(0,xmax*1.02); ax.set_yticks([]); ax.spines['left'].set_visible(False); ax.set_xlabel('substitutions per site (GTR+F+I, 1,000 ultrafast bootstraps)')
ax.set_xticks([0,0.0001,0.0002,0.0003,0.0004]); ax.set_xticklabels(['0','1e-4','2e-4','3e-4','4e-4'])
ax2=fig.add_axes([0.545,0.04,0.42,0.93]); ax2.set_xlim(0,10); ax2.set_ylim(len(tips)-0.3,-3.0); ax2.axis('off')
src=[]
dcol={'Yes':'#0b0b0b','No':'#e1e0d9'}
for l in tips:
    s=l['name']; y=l['y']; c=CL[s]
    ax2.add_patch(Rectangle((0.0,y-0.4),0.9,0.8,color=LC[LIN[s]],lw=0)); ax2.add_patch(Rectangle((1.1,y-0.4),0.9,0.8,color=dcol[c['diabetes']],lw=0))
    mdr=any(d in res_dm[s] for d in['isoniazid']) and any(d in res_dm[s] for d in['rifampicin'])
    ax2.add_patch(Rectangle((2.2,y-0.4),0.9,0.8,color='#d03b3b' if res_dm[s] else '#e1e0d9',lw=0))
    ax2.text(3.4,y,s,va='center',fontsize=6.5,color=INK)
    ax2.text(5.0,y,('flagged: %s'%('contaminated' if s=='MTB25' else 'low coverage')) if s in LOW else '',va='center',fontsize=6,color=INK2)
    src.append(dict(isolate=s,lineage=LIN[s],diabetes=c['diabetes'],predicted_resistance=';'.join(res_dm[s]),low_quality_flag='yes' if s in LOW else '',tree_y=int(round(y))))
for x,t in[(0.45,'Lineage'),(1.55,'DM'),(2.65,'Drug-R')]: ax2.text(x,-0.9,t,rotation=90,ha='center',va='bottom',fontsize=7)
fig.legend(handles=[Patch(color=LC[k],label=k) for k in('L1','L2','L3','L4','unassigned')]+[Patch(color='#0b0b0b',label='Diabetes: yes'),Patch(color='#e1e0d9',label='Diabetes: no'),Patch(color='#d03b3b',label='Any WHO-catalogue drug-resistance mutation'),Line2D([0],[0],marker='o',color='w',markerfacecolor='#0b0b0b',markersize=4,label='Node bootstrap >= 95')],
           loc='lower center',ncol=4,frameon=False,fontsize=7,bbox_to_anchor=(0.5,-0.07))
fig.suptitle('Maximum-likelihood phylogeny of 55 isolates (rooted on M. canettii; outgroup and root stem not to scale)',fontsize=9,y=1.0)
fig.suptitle('Maximum-likelihood phylogeny of 55 isolates (rooted on M. canettii; outgroup and root stem not to scale)',fontsize=9,y=0.995)
save(fig,'fig2_annotated_phylogeny'); pd.DataFrame(src).to_csv(OUT+'/fig2_annotated_phylogeny_source_data.tsv',sep='\t',index=False)
# =================================================================== F1 presence/absence heatmap
pm=pd.read_csv('/mnt/user-data/outputs/cohort_C1_C2/C2/a16_presence_matrix.tsv',sep='\t')
pm=pm[pm.tier.astype(str).isin(['1a','1b'])]
good=[s for s in S if s not in LOW]
var=pm[(pm[good]=='absent').any(axis=1)].copy(); var['n_absent']=(var[S]=='absent').sum(axis=1); var=var.sort_values(['n_absent'],ascending=False)
sm=np.zeros((len(var),len(S)),int)
for j,s in enumerate(order):
    sm[:,j]=var[s].map({'present':0,'absent':1,'not determined':2}).values
fig=plt.figure(figsize=(9,0.17*len(var)+3)); ax=fig.add_axes([0.14,0.17,0.78,0.71])
from matplotlib.colors import ListedColormap
ax.imshow(sm,aspect='auto',cmap=ListedColormap(['#cde2fb','#0d366b','#c3c2b7']),interpolation='nearest',vmin=0,vmax=2)
for i in range(len(var)):
    for j in range(len(order)):
        if sm[i,j]==2: ax.add_patch(Rectangle((j-0.5,i-0.5),1,1,fill=False,hatch='////',edgecolor='white',lw=0))
ax.set_yticks(range(len(var))); ax.set_yticklabels(['%s (%s, T%s)'%(g,r,t) for g,r,t in zip(var.gene,var.rv,var.tier)],fontsize=6.5); ax.set_xticks(range(len(order))); ax.set_xticklabels(order,rotation=90,fontsize=6)
ax.spines[:].set_visible(False); ax.tick_params(length=0)
a=fig.add_axes([0.14,0.895,0.78,0.02]); a.imshow([[list(LC).index(LIN[s]) for s in order]],aspect='auto',cmap=ListedColormap(list(LC.values())),vmin=0,vmax=4); a.axis('off'); a.text(-0.01,0.5,'Lineage',transform=a.transAxes,ha='right',va='center',fontsize=7)
b=fig.add_axes([0.14,0.92,0.78,0.02]); b.imshow([[1 if CL[s]['diabetes']=='Yes' else 0 for s in order]],aspect='auto',cmap=ListedColormap(['#e1e0d9','#0b0b0b']),vmin=0,vmax=1); b.axis('off'); b.text(-0.01,0.5,'Diabetes',transform=b.transAxes,ha='right',va='center',fontsize=7)
fig.legend(handles=[Patch(color='#cde2fb',label='present'),Patch(color='#0d366b',label='absent (>=90% genome at 10x and gene <10% covered)'),Patch(facecolor='#c3c2b7',hatch='////',edgecolor='white',label='not determined')]+[Patch(color=LC[k],label=k) for k in LC]+[Patch(color='#0b0b0b',label='diabetes yes'),Patch(color='#e1e0d9',label='diabetes no')],
           loc='upper center',ncol=4,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.085))
fig.text(0.5,-0.035,'esxR/esxS (Rv3019c/Rv3020c) are absent in about two thirds of isolates by both reads and assemblies; these two calls sit in a repeated ESX gene family and are not interpreted biologically.',ha='center',fontsize=6.5,color=INK2)
fig.suptitle('Tier 1 virulence genes with at least one absence call (isolates ordered as in the phylogeny)',fontsize=9,y=0.985)
save(fig,'fig1_virulence_presence_absence_heatmap')
out=var[['rv','gene','tier','system','n_absent']].copy()
for s in order: out[s]=var[s].values
out['lineage_row']=''; hdr=pd.DataFrame([dict(rv='lineage',**{s:LIN[s] for s in order}),dict(rv='diabetes',**{s:CL[s]['diabetes'] for s in order})])
pd.concat([hdr,out],ignore_index=True).to_csv(OUT+'/fig1_virulence_presence_absence_heatmap_source_data.tsv',sep='\t',index=False)
# =================================================================== F3 dN/dS by system
cat=json.load(open(P+'/ref/catalogue.json'))
e=pd.read_csv('/mnt/user-data/outputs/cohort_C3/e15_genomewide_dNdS_per_gene.tsv',sep='\t')
e['rv']=e['rv'].astype(str); e['system']=e['rv'].map(lambda r:cat[r]['system'] if r in cat else 'Non-panel genes')
e=e[e.status=='ok'].copy()
for c in['N_subs','S_subs','N_sites','S_sites','omega_M0']: e[c]=pd.to_numeric(e[c],errors='coerce')
rng=np.random.default_rng(20261106)
def pooled(d):
    ns=d.N_subs.sum(); ss=d.S_subs.sum();
    return (ns/d.N_sites.sum())/(ss/d.S_sites.sum()) if ss>0 else np.nan
rows=[]
for sysn,d in e.groupby('system'):
    d=d.dropna(subset=['N_sites','S_sites'])
    if len(d)<8: continue
    pt=pooled(d); bs=[pooled(d.iloc[rng.integers(0,len(d),len(d))]) for _ in range(1000)]; bs=[b for b in bs if np.isfinite(b)]
    om=d.omega_M0[(d.omega_M0<50)&(d.omega_M0>0)]
    rows.append(dict(system=sysn,n_genes=len(d),N_subs=d.N_subs.sum(),S_subs=d.S_subs.sum(),pooled_dN_dS=pt,ci_lo=np.percentile(bs,2.5),ci_hi=np.percentile(bs,97.5),median_gene_omega=om.median() if len(om) else np.nan,n_genes_omega=len(om)))
sd=pd.DataFrame(rows).sort_values('pooled_dN_dS'); sd.to_csv(OUT+'/fig3_dNdS_by_system_source_data.tsv',sep='\t',index=False)
e[['rv','gene','tier','system','status','omega_M0','N_subs','S_subs','N_sites','S_sites']].to_csv(OUT+'/fig3_dNdS_per_gene_source_data.tsv',sep='\t',index=False)
fig,(a1,a2)=plt.subplots(1,2,figsize=(9,5.2),gridspec_kw=dict(width_ratios=[1,1.1]),sharey=True)
ys=np.arange(len(sd)); a1.errorbar(sd.pooled_dN_dS,ys,xerr=[sd.pooled_dN_dS-sd.ci_lo,sd.ci_hi-sd.pooled_dN_dS],fmt='o',color='#2a78d6',ecolor='#2a78d6',ms=4,lw=1.2,capsize=0)
a1.axvline(1,color='#d03b3b',lw=1,ls='--'); a1.text(1.02,len(sd)-0.3,'neutral (dN/dS = 1)',color='#d03b3b',fontsize=7,va='top')
a1.set_yticks(ys); a1.set_yticklabels(['%s (n=%d)'%(s,n) for s,n in zip(sd.system,sd.n_genes)],fontsize=7); a1.set_xlabel('pooled dN/dS (95% gene-bootstrap CI)'); a1.set_xlim(0,max(1.3,sd.ci_hi.max()*1.05)); a1.grid(axis='x',color=GRID,lw=0.5)
for i,s in enumerate(sd.system):
    d=e[(e.system==s)&(e.omega_M0>0)&(e.omega_M0<50)];
    a2.scatter(d.omega_M0,i+rng.uniform(-0.28,0.28,len(d)),s=5,color='#6250d6',alpha=0.45,lw=0)
    a2.plot([d.omega_M0.median()]*2,[i-0.35,i+0.35],color=INK,lw=1.4)
a2.text(1.2e-4,len(sd)-0.2,'points at 1e-4 = no nonsynonymous change (omega at lower bound)',fontsize=6,color=INK2,va='bottom'); a2.set_xscale('log'); a2.axvline(1,color='#d03b3b',lw=1,ls='--'); a2.set_xlabel('per-gene M0 omega (dN/dS), log scale; bar = median'); a2.grid(axis='x',color=GRID,lw=0.5)
fig.suptitle('Selection pressure by functional system: pooled and per-gene dN/dS (codeml M0, 54 isolates)',fontsize=9)
save(fig,'fig3_dNdS_by_gene_and_system')
# =================================================================== F4 variant plots
doms=json.load(open(P+'/ref/domains.json')); epi=json.load(open(P+'/ref/epitopes_mapped.json')); genes=json.load(open(P+'/ref/genes.json'))
GOI=[('Rv1908c','katG'),('Rv3805c','aftB'),('Rv3617','ephA'),('Rv3826','fadD23'),('Rv0836c','vapC36'),('Rv0422c','mycP2')]
# resolve rv by gene name from catalogue / genes.json
name2rv={}
for rv,g in genes.items():
    if g.get('name'): name2rv.setdefault(g['name'],rv)
GOI=[(name2rv.get(n,r),n) for r,n in GOI]
vrows=[]
AAL={}
for s in S:
    R=json.load(open('%s/work/%s/results.json'%(P,s)))
    for v in R['variants']:
        if v.get('in_mask') or v['type']!='snp' and v['effect']=='' : continue
        for rv,n in GOI:
            if v.get('rv')==rv and v.get('codon'):
                vrows.append(dict(rv=rv,gene=n,isolate=s,codon=v['codon'],effect=v['effect'],hgvs_p=v['hgvs_p'],lineage=LIN[s]))
vd=pd.DataFrame(vrows); vd.to_csv(OUT+'/fig4_variant_distribution_source_data_per_isolate.tsv',sep='\t',index=False)
fig,axes=plt.subplots(len(GOI),1,figsize=(8,2.0*len(GOI)),gridspec_kw=dict(hspace=0.9))
EC={'missense_variant':'#e34948','synonymous_variant':'#898781'}
dsrc=[]
for ax,(rv,n) in zip(axes,GOI):
    L=(genes[rv]['end']-genes[rv]['start']+1)//3; ax.set_xlim(0,L+1); ax.set_ylim(-1.7,4.2)
    ax.add_patch(Rectangle((0,-0.1),L,0.2,color='#c3c2b7',lw=0))
    ds=[d for d in doms.get(rv,[]) if d['type']=='DOMAIN']
    for d in ds:
        ax.add_patch(Rectangle((d['aa_start'],-0.55),d['aa_end']-d['aa_start'],0.45,color='#9ec5f4',lw=0)); ax.text((d['aa_start']+d['aa_end'])/2,-0.33,d['note'][:30],ha='center',va='center',fontsize=5.5,clip_on=True); dsrc.append(dict(gene=n,rv=rv,track='domain',start=d['aa_start'],end=d['aa_end'],label=d['note']))
    ep=[x for x in epi if x['rv']==rv]
    for x in ep:
        ax.add_patch(Rectangle((x['aa_start'],-1.35),x['aa_end']-x['aa_start']+1,0.5,color='#eb6834',alpha=0.55,lw=0)); dsrc.append(dict(gene=n,rv=rv,track='T-cell epitope',start=x['aa_start'],end=x['aa_end'],label=x['pep']))
    if not ds: ax.text(L*0.98,-0.33,'no UniProt domain annotation',ha='right',fontsize=6,color=INK2)
    if not ep: ax.text(L*0.98,-1.1,'no mapped T-cell epitope',ha='right',fontsize=6,color=INK2)
    g=vd[vd.rv==rv].groupby(['codon','effect','hgvs_p']).size().reset_index(name='n_isolates') if len(vd) else pd.DataFrame()
    for _,r in g.iterrows():
        col=EC.get(r.effect,'#6250d6'); ax.plot([r.codon,r.codon],[0.1,0.1+min(r.n_isolates,40)/40*3.4],color=col,lw=0.8); ax.scatter([r.codon],[0.1+min(r.n_isolates,40)/40*3.4],s=14,color=col,zorder=3)
        if r.n_isolates>=10: ax.text(r.codon,0.1+min(r.n_isolates,40)/40*3.4+0.25,r.hgvs_p.replace('p.',''),ha='center',fontsize=6)
        dsrc.append(dict(gene=n,rv=rv,track='variant',start=r.codon,end=r.codon,label=r.hgvs_p,effect=r.effect,n_isolates=r.n_isolates))
    ax.set_yticks([]); ax.spines['left'].set_visible(False); ax.set_title('%s (%s), %d aa  -  %d variant sites'%(n,rv,L,len(g)),fontsize=8,loc='left'); ax.set_xlabel('amino-acid position',fontsize=7)
fig.legend(handles=[Line2D([0],[0],color='#e34948',marker='o',ms=4,label='missense'),Line2D([0],[0],color='#898781',marker='o',ms=4,label='synonymous'),Line2D([0],[0],color='#6250d6',marker='o',ms=4,label='frameshift / stop / other'),Patch(color='#9ec5f4',label='UniProt domain'),Patch(color='#eb6834',alpha=0.55,label='mapped T-cell epitope (IEDB)')],loc='lower center',ncol=5,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.05))
fig.suptitle('Variant positions with protein domains and T-cell epitopes on one axis (stem height = number of isolates, capped at 40)',fontsize=9,y=0.93)
save(fig,'fig4_variant_distribution_domains_epitopes'); pd.DataFrame(dsrc).to_csv(OUT+'/fig4_variant_distribution_source_data.tsv',sep='\t',index=False)
# =================================================================== F5 lineage-stratified allele frequencies
a20=pd.read_csv('/mnt/user-data/outputs/cohort_C1_C2/C2/a20_lineage_allele_profiles.tsv',sep='\t')
a20=a20[a20.tier.astype(str).isin(['1a','1b'])&(a20.n_isolates>=3)&(~a20.only_in_MTB25_flagged)&a20.effect.isin(['missense_variant','stop_gained','start_lost','frameshift_variant','disruptive_inframe_deletion','disruptive_inframe_insertion'])].copy()
lin_cols=['freq_L1','freq_L3','freq_L4']
nL={'L1':sum(LIN[s]=='L1' for s in S),'L3':sum(LIN[s]=='L3' for s in S),'L4':sum(LIN[s]=='L4' for s in S)}
a20['lab']=a20.gene.astype(str)+' '+a20.hgvs_p.fillna('').astype(str)
a20['dominant']=a20[lin_cols].values.argmax(1); a20['maxf']=a20[lin_cols].max(axis=1)
a20=a20.sort_values(['dominant','maxf','n_isolates'],ascending=[True,False,False]).reset_index(drop=True)
half=int(np.ceil(len(a20)/2)); fig,axs=plt.subplots(1,2,figsize=(9,0.115*half+1.6),gridspec_kw=dict(wspace=0.8))
import matplotlib.colors as mc
cm=mc.LinearSegmentedColormap.from_list('seq',['#f6f9fe','#b7d3f6','#3987e5','#184f95','#0d366b'])
for ax,part in zip(axs,[a20.iloc[:half],a20.iloc[half:]]):
    im=ax.imshow(part[lin_cols].values,aspect='auto',cmap=cm,vmin=0,vmax=1,interpolation='nearest')
    ax.set_yticks(range(len(part))); ax.set_yticklabels(part.lab,fontsize=4.8); ax.set_xticks(range(3)); ax.set_xticklabels(['L1\n(n=%d)'%nL['L1'],'L3\n(n=%d)'%nL['L3'],'L4\n(n=%d)'%nL['L4']],fontsize=7); ax.xaxis.tick_top(); ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
cb=fig.colorbar(im,ax=axs,orientation='horizontal',fraction=0.02,pad=0.02,shrink=0.4); cb.set_label('allele frequency within lineage',fontsize=7)
fig.suptitle('Tier 1 protein-altering variants carried by >=3 isolates: allele frequency within each lineage\n(L2 has 1 isolate and 3 isolates are unassigned; both are in the source table, not shown)',fontsize=8,y=1.0)
save(fig,'fig5_lineage_stratified_allele_frequencies'); a20.drop(columns=['dominant','maxf']).to_csv(OUT+'/fig5_lineage_stratified_allele_frequencies_source_data.tsv',sep='\t',index=False)
print('figures done',os.listdir(OUT))
print(sd[['system','n_genes','pooled_dN_dS','ci_lo','ci_hi']].round(3).to_string())
print(len(var),'genes in fig1;',len(a20),'variants in fig5; GOI',GOI)
