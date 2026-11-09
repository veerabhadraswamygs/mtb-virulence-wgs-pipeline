"""Regulatory reference: TSS / 5'UTR (Shell 2015 + Cortes 2013) and ChIP-seq peaks (Minch 2015).
Columns are read BY POSITION after locating the header row, because these workbooks
contain duplicated column names (a known trap in Cortes S1/S6).
"""
import json,re
from openpyxl import load_workbook
U='/mnt/user-data/uploads/Downloads/limitation fix files'
R='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/ref'
genes=json.load(open(R+'/genes.json'))

def header_row(ws,needle,maxscan=12):
    for i,row in enumerate(ws.iter_rows(values_only=True)):
        vals=[str(x).strip() if x is not None else '' for x in row]
        if any(needle.lower() in v.lower() for v in vals): return i,vals
        if i>maxscan: break
    return None,None

# ---- Shell 2015 TSS ----
wb=load_workbook(U+'/supplementary_tables/pgen.1005641.s010.xlsx',read_only=True,data_only=True)
ws=wb['All TSSs']
rows=list(ws.iter_rows(values_only=True))
hdr=[str(x).strip() if x is not None else '' for x in rows[0]]
iCoord=[i for i,h in enumerate(hdr) if h.lower().startswith('coordinate')][0]
iStrand=[i for i,h in enumerate(hdr) if h.lower()=='strand'][0]
iDown=[i for i,h in enumerate(hdr) if h.lower().startswith('downstream same-strand')][0]
iRel=[i for i,h in enumerate(hdr) if 'relative to d' in h.lower()][0]
tss=[]
for r in rows[1:]:
    if r[iCoord] is None: continue
    try: c=int(r[iCoord])
    except: continue
    tss.append({'pos':c,'strand':str(r[iStrand]).strip(),'gene':(str(r[iDown]).strip() if r[iDown] else ''),
                'rel':(r[iRel] if r[iRel] not in (None,'') else None),'src':'Shell2015'})
wb.close()
print('Shell TSS:',len(tss),'| with downstream gene:',sum(1 for t in tss if t['gene']))

# map TSS to genes, pick the one with smallest |offset| to the start codon
by_gene={}
for t in tss:
    rv=t['gene']
    if rv not in genes: continue
    g=genes[rv]
    if g['strand']!=t['strand']: continue
    start=g['start'] if g['strand']=='+' else g['end']
    off=(start-t['pos']) if g['strand']=='+' else (t['pos']-start)
    if off<0 or off>500: continue          # TSS must be upstream, within 500 nt
    cur=by_gene.get(rv)
    if cur is None or abs(off)<abs(cur['utr']):
        by_gene[rv]={'tss':t['pos'],'utr':off,'strand':g['strand'],'src':'Shell2015'}
print('genes with mapped TSS:',len(by_gene))
ll=sum(1 for v in by_gene.values() if v['utr']<=5)
print('leaderless (5\'UTR <= 5 nt, derived from coordinates):',ll)

# promoter window: 5'UTR plus 100 bp upstream of the TSS
for rv,v in by_gene.items():
    g=genes[rv]
    if v['strand']=='+': v['prom_start'],v['prom_end']=v['tss']-100,g['start']-1
    else:                v['prom_start'],v['prom_end']=g['end']+1,v['tss']+100
    if v['prom_start']>v['prom_end']: v['prom_start'],v['prom_end']=v['prom_end'],v['prom_start']
json.dump(by_gene,open(R+'/tss.json','w'))

# ---- Minch 2015 ChIP-seq peaks ----
wb=load_workbook(U+'/ncomms6829-s2.xlsx',read_only=True,data_only=True)
ws=wb['Sheet1']; rows=list(ws.iter_rows(values_only=True))
hi,hdr=None,None
for i,r in enumerate(rows[:12]):
    vals=[str(x).strip() if x is not None else '' for x in r]
    if 'Regulator' in vals: hi,hdr=i,vals; break
iReg=hdr.index('Regulator'); iGene=hdr.index('Gene')
iS=[i for i,h in enumerate(hdr) if h.lower().startswith('peak start')][0]
iE=[i for i,h in enumerate(hdr) if h.lower().startswith('peak stop')][0]
peaks=[]
for r in rows[hi+1:]:
    try: s,e=int(r[iS]),int(r[iE])
    except: continue
    peaks.append({'reg':str(r[iReg]).strip(),'gene':str(r[iGene]).strip(),'start':min(s,e),'end':max(s,e)})
wb.close()
print('Minch peaks:',len(peaks),'| regulators:',len({p['reg'] for p in peaks}))
cov=sum(p['end']-p['start']+1 for p in peaks)
json.dump(peaks,open(R+'/chipseq.json','w'))
print('peak bp (unmerged): %d'%cov)
