"""TSS / 5'UTR reference from Shell 2015, using the authors' own Leadered / Leaderless
classification sheets. Columns read by position after locating them in the header."""
import warnings,json; warnings.filterwarnings('ignore')
from openpyxl import load_workbook
U='/mnt/user-data/uploads/Downloads/limitation fix files'
R='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/ref'
genes=json.load(open(R+'/genes.json'))
wb=load_workbook(U+'/supplementary_tables/pgen.1005641.s010.xlsx',read_only=True,data_only=True)
def cols(hdr):
    f=lambda s:[i for i,h in enumerate(hdr) if h.lower().startswith(s)][0]
    return (f('coordinate'),[i for i,h in enumerate(hdr) if h.lower()=='strand'][0],
            f('downstream same-strand gene'),[i for i,h in enumerate(hdr) if 'relative to d' in h.lower()][0],
            f('overlapping same-strand gene'))
out={}; stat={'leadered':0,'leaderless':0,'skipped_no_gene':0,'skipped_strand':0,'skipped_range':0}
for sheet,kind in [('Leadered','leadered'),('Leaderless','leaderless')]:
    ws=wb[sheet]; rows=list(ws.iter_rows(values_only=True))
    hdr=[str(x).strip() if x is not None else '' for x in rows[0]]
    iC,iS,iD,iR,iO=cols(hdr)
    for r in rows[1:]:
        try: pos=int(r[iC])
        except: continue
        strand=str(r[iS]).strip()
        rv=str(r[iD]).strip() if r[iD] else ''
        if kind=='leaderless':
            ov=str(r[iO]).strip() if r[iO] else ''
            if ov and ov.lower()!='none' and ov in genes: rv=ov
        if rv not in genes: stat['skipped_no_gene']+=1; continue
        g=genes[rv]
        if g['strand']!=strand: stat['skipped_strand']+=1; continue
        start=g['start'] if strand=='+' else g['end']
        utr=(start-pos) if strand=='+' else (pos-start)
        if kind=='leadered' and not (0<=utr<=5000): stat['skipped_range']+=1; continue
        if kind=='leaderless': utr=max(0,utr)
        cur=out.get(rv)
        if cur is None or abs(utr)<abs(cur['utr']):
            out[rv]={'tss':pos,'utr':int(utr),'strand':strand,'class':kind,'src':'Shell2015','long_utr':bool(utr>300)}
        stat[kind]+=1
wb.close()
for rv,v in out.items():
    g=genes[rv]
    if v['strand']=='+': a,b=v['tss']-100,g['start']-1
    else:                a,b=g['end']+1,v['tss']+100
    v['prom_start'],v['prom_end']=min(a,b),max(a,b)
print('TSS rows used:',stat)
print('genes with TSS:',len(out))
print('   leaderless genes:',sum(1 for v in out.values() if v['class']=='leaderless'),
      '| leadered genes:',sum(1 for v in out.values() if v['class']=='leadered'))
u=[v['utr'] for v in out.values() if v['class']=='leadered']
u.sort(); print('leadered 5\'UTR length: median %d nt, range %d-%d'%(u[len(u)//2],u[0],u[-1]))
json.dump(out,open(R+'/tss.json','w'))
for rv in ['Rv3875','Rv3874','Rv0757','Rv1908c','Rv2031c']:
    print('  ',rv,out.get(rv))
