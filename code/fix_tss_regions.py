"""Define promoter and 5'UTR intervals per gene, circular-genome safe.
promoter = 100 bp immediately upstream of the TSS, through the TSS
5'UTR    = TSS to the start codon, only when the UTR is <= 500 nt; longer ones are
           flagged distal (likely operon-level) and no UTR interval is emitted.
"""
import json
R='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/ref'
GL=4411532
genes=json.load(open(R+'/genes.json'))
tss=json.load(open(R+'/tss.json'))
def norm(p):
    p=((p-1)%GL)+1
    return p
def span(a,b):
    """inclusive circular interval from a to b going forward; returns 1 or 2 plain intervals"""
    a,b=norm(a),norm(b)
    return [[a,b]] if a<=b else [[a,GL],[1,b]]
UTR_CAP=500
out={}
for rv,t in tss.items():
    g=genes[rv]; s=t['strand']; p=t['tss']
    start=g['start'] if s=='+' else g['end']
    utr=(start-p)%GL if s=='+' else (p-start)%GL
    if utr>GL//2: utr=utr-GL
    prom = span(p-100,p) if s=='+' else span(p,p+100)
    rec={'tss':p,'strand':s,'class':t['class'],'utr':int(utr),
         'distal':bool(utr>UTR_CAP or utr<0),'promoter':prom}
    if 0<utr<=UTR_CAP:
        rec['utr_region']= span(p,start-1) if s=='+' else span(start+1,p)
    else:
        rec['utr_region']=[]
    out[rv]=rec
json.dump(out,open(R+'/tss_regions.json','w'))
w=[sum(b-a+1 for a,b in v['promoter']) for v in out.values()]
u=[sum(b-a+1 for a,b in v['utr_region']) for v in out.values() if v['utr_region']]
print('genes:',len(out))
print('promoter width: min %d max %d (all should be 101)'%(min(w),max(w)))
print('UTR regions emitted:',len(u),'| max width',max(u) if u else 0)
print('flagged distal (UTR>500 or negative):',sum(1 for v in out.values() if v['distal']))
for rv in ['Rv0001','Rv3924c','Rv0757','Rv2031c']:
    print(' ',rv,{k:out[rv][k] for k in ('tss','utr','class','distal','promoter','utr_region')})
# post-check: every emitted UTR interval must have exactly the stated length
bad=0
for rv,v in out.items():
    w=sum(b-a+1 for a,b in v['utr_region'])
    if v['utr_region'] and w!=v['utr']: bad+=1; print('LENGTH MISMATCH',rv,v['utr'],w)
print('UTR length check failures:',bad)
