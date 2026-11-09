"""Coverage: per-gene breadth, callable mask (brief 5), exclusion rule,
large deletions / RD (29), TbD1 locus (28a)."""
import gzip,json,collections,sys
import numpy as np
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
GL=4411532
def load_depth(wd,fn='depth.txt.gz'):
    d=np.zeros(GL+1,dtype=np.int32)
    with gzip.open(wd+'/'+fn,'rt') as fh:
        for line in fh:
            f=line.split('\t')
            d[int(f[1])]=int(f[2])
    return d
def gene_breadth(d,genes,mind=10):
    out={}
    for rv,g in genes.items():
        if not g['start']: continue
        seg=d[g['start']:g['end']+1]
        L=len(seg)
        if L==0: continue
        out[rv]={'len':L,'mean':float(seg.mean()),'breadth10':float((seg>=mind).mean()),
                 'breadth1':float((seg>0).mean()),'median':float(np.median(seg))}
    return out
def callable_bed(d,mind=10):
    ok=d[1:]>=mind
    bed=[];start=None
    idx=np.flatnonzero(np.diff(np.concatenate(([0],(~ok).view(np.int8),[0]))))
    # build low-coverage intervals
    low=~ok; runs=[]
    i=0; n=len(low)
    changes=np.flatnonzero(np.diff(low.astype(np.int8)))
    pos=0; cur=low[0]
    for c in changes:
        if cur: runs.append((pos+1,c+1))
        pos=c+1; cur=low[pos]
    if cur: runs.append((pos+1,n))
    return runs
def deletions(d,genes,minlen=1000):
    """contiguous runs with zero coverage >= minlen"""
    z=d[1:]==0
    runs=[];pos=0;cur=z[0];n=len(z)
    for c in np.flatnonzero(np.diff(z.astype(np.int8))):
        if cur: runs.append((pos+1,c+1))
        pos=c+1; cur=z[pos]
    if cur: runs.append((pos+1,n))
    out=[]
    for s,e in runs:
        if e-s+1<minlen: continue
        hit=[rv for rv,g in genes.items() if g['start'] and not(g['end']<s or g['start']>e)]
        out.append({'start':int(s),'end':int(e),'length':int(e-s+1),'genes':sorted(hit)})
    return out
if __name__=='__main__':
    s=sys.argv[1]; wd=P+'/work/'+s
    genes=json.load(open(P+'/ref/genes.json'))
    d=load_depth(wd)
    tot=d[1:]
    print('genome mean %.1fx | >=10x %.2f%% | >=1x %.2f%%'%(tot.mean(),100*(tot>=10).mean(),100*(tot>0).mean()))
    gb=gene_breadth(d,genes)
    low=[rv for rv,v in gb.items() if v['breadth10']<0.95]
    print('genes with <95%% at >=10x: %d of %d'%(len(low),len(gb)))
    cat=json.load(open(P+'/ref/catalogue.json'))
    lowc=[rv for rv in cat if rv in gb and gb[rv]['breadth10']<0.95]
    print('catalogue genes not fully callable: %d of %d'%(len(lowc),len(cat)))
    dall=load_depth(wd,'depth_all.txt.gz')
    dels=deletions(dall,genes)   # MQ-UNFILTERED: repeats such as IS6110 map with MQ 0 and
                                 # would otherwise appear as spurious deletions
    print('zero-coverage deletions >=1kb:',len(dels))
    for x in dels[:8]: print('   %d-%d (%d bp) %s'%(x['start'],x['end'],x['length'],','.join(g for g in x['genes'][:5])))
    # TbD1
    for rv in ['Rv2073c','Rv2074']:
        v=gb.get(rv); print('TbD1 %s: mean %.1fx breadth10 %.3f'%(rv,v['mean'],v['breadth10']))
    low_runs=callable_bed(d)
    bp=sum(e-s+1 for s,e in low_runs)
    print('non-callable (<10x) regions: %d intervals, %d bp (%.2f%%)'%(len(low_runs),bp,100*bp/GL))
    print('genome mean (MQ-unfiltered) %.1fx | >=1x %.2f%%'%(dall[1:].mean(),100*(dall[1:]>0).mean()))
    json.dump({'gene_breadth':gb,'deletions':dels,'mean_depth_mqfiltered':float(tot.mean()),
               'mean_depth_all':float(dall[1:].mean()),
               'noncallable':[[int(a),int(b)] for a,b in low_runs]},open(wd+'/cov.json','w'))
