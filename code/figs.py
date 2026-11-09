"""Publication figures. Each figure is written as its own PNG at 300 dpi so it can be
screenshotted or placed individually; nothing is combined into multi-panel images."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, gzip, collections, json, os
import palette as C
plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','axes.linewidth':0.8,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'figure.dpi':300,'savefig.dpi':300,'savefig.bbox':'tight'})
def _save(fig,path):
    fig.savefig(path,facecolor='white'); plt.close(fig); return path

def fig_coverage(depth_all_gz,out,window=5000):
    d=[]
    with gzip.open(depth_all_gz,'rt') as fh:
        for l in fh: d.append(int(l.split('\t')[2]))
    a=np.array(d); n=len(a)//window*window
    m=a[:n].reshape(-1,window).mean(axis=1)
    x=np.arange(len(m))*window/1e6
    fig,ax=plt.subplots(figsize=(7.2,2.4))
    ax.fill_between(x,m,color=C.BLUE,alpha=0.85,linewidth=0)
    ax.axhline(np.median(a),color=C.ORANGE,lw=1,ls='--',label='median %.0f×'%np.median(a))
    ax.set_xlabel('H37Rv position (Mb)'); ax.set_ylabel('depth (×)')
    ax.set_xlim(0,x[-1]); ax.set_ylim(0,np.percentile(m,99.5)*1.15)
    ax.legend(frameon=False,loc='upper right')
    return _save(fig,out)

def fig_effects(rows,out):
    c=collections.Counter(r['effect'] for r in rows if r['effect'])
    items=c.most_common(10)[::-1]
    if not items: return None
    fig,ax=plt.subplots(figsize=(5.6,max(2.2,0.32*len(items)+0.8)))
    y=np.arange(len(items))
    ax.barh(y,[v for _,v in items],color=C.BLUE,height=0.65)
    ax.set_yticks(y); ax.set_yticklabels([k.replace('_',' ') for k,_ in items])
    ax.set_xlabel('variants')
    for i,(k,v) in enumerate(items): ax.text(v+max(1,0.01*items[-1][1]),i,str(v),va='center',fontsize=8)
    ax.set_xlim(0,items[-1][1]*1.12)
    return _save(fig,out)

def fig_systems(rows,out):
    c=collections.Counter(r['functional_system'] for r in rows
                          if r['in_virulence_catalogue'] and not r['in_mask'] and r['functional_system'])
    items=c.most_common(12)[::-1]
    if not items: return None
    fig,ax=plt.subplots(figsize=(6.4,max(2.2,0.34*len(items)+0.8)))
    y=np.arange(len(items))
    ax.barh(y,[v for _,v in items],color=C.GREEN,height=0.65)
    ax.set_yticks(y); ax.set_yticklabels([k[:42] for k,_ in items],fontsize=8)
    ax.set_xlabel('variants (unmasked, virulence catalogue)')
    for i,(k,v) in enumerate(items): ax.text(v+0.3,i,str(v),va='center',fontsize=8)
    ax.set_xlim(0,items[-1][1]*1.14)
    return _save(fig,out)

def fig_af(rows,subclonal,out):
    fig,ax=plt.subplots(figsize=(5.2,2.4))
    af=[r['allele_fraction'] for r in rows]
    sc=[s['allele_fraction'] for s in subclonal]
    ax.hist(af,bins=np.arange(0,1.025,0.025),color=C.BLUE,label='fixed variants (n=%d)'%len(af))
    if sc: ax.hist(sc,bins=np.arange(0,1.025,0.025),color=C.ORANGE,label='sub-clonal (n=%d)'%len(sc))
    ax.axvline(0.9,color=C.GREY,ls='--',lw=0.9)
    ax.set_xlabel('alternate allele fraction'); ax.set_ylabel('variants'); ax.set_yscale('log')
    ax.legend(frameon=False,fontsize=8)
    return _save(fig,out)

def fig_genome_map(dels,junctions,out):
    GL=4411532
    fig,ax=plt.subplots(figsize=(7.2,1.7))
    ax.add_patch(plt.Rectangle((0,0.42),GL/1e6,0.16,color='#E8E8E8',lw=0))
    for d in dels:
        ax.add_patch(plt.Rectangle((d['start']/1e6,0.42),max(d['length']/1e6,0.004),0.16,color=C.ORANGE,lw=0))
    for j in junctions:
        ax.plot([j['start']/1e6],[0.72],marker='v',color=C.BLUE,ms=4)
    ax.set_xlim(0,GL/1e6); ax.set_ylim(0.2,0.95); ax.set_yticks([])
    ax.set_xlabel('H37Rv position (Mb)')
    ax.plot([],[],marker='v',color=C.BLUE,ls='',label='IS6110 junction (n=%d)'%len(junctions))
    ax.add_patch(plt.Rectangle((0,0),0,0,color=C.ORANGE,label='deletion ≥1 kb (n=%d)'%len(dels)))
    ax.legend(frameon=False,fontsize=8,ncol=2,loc='upper center',bbox_to_anchor=(0.5,1.35))
    ax.spines['left'].set_visible(False)
    return _save(fig,out)

def fig_epitope(pp,out):
    fig,ax=plt.subplots(figsize=(4.0,2.6))
    ks=['epitope','non_epitope']
    v=[pp[k]['pi_per_kb'] or 0 for k in ks]
    ax.bar([0,1],v,color=[C.ORANGE,C.BLUE],width=0.55)
    for i,x in enumerate(v): ax.text(i,x,'%.3f'%x,ha='center',va='bottom',fontsize=8)
    ax.set_xticks([0,1]); ax.set_xticklabels(['epitope','non-epitope'])
    ax.set_ylabel('variants per kb of callable sequence')
    ax.set_ylim(0,max(v)*1.25 if max(v)>0 else 1)
    return _save(fig,out)
