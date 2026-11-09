"""C4 analysis 25b: FoldX 5.1 RepairPDB (once per gene) + BuildModel (3 runs) on AlphaFold DB models for pre-specified candidates:
Tier 1a/1b missense, >=3 non-MTB25 isolates, pLDDT>=70, buried or <=8 A from an annotated functional site, protein <=700 aa (compute limit); plus katG S315T as a known-resistance reference."""
import csv,json,os,sys,collections,concurrent.futures as cf
sys.path.insert(0,'/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/code')
import mod_struct as ms
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c4'; os.makedirs(O+'/foldx',exist_ok=True)
rs=json.load(open(P+'/ref/rv2struct.json'))
r=list(csv.DictReader(open(O+'/a25_variant_structure_descriptors.tsv'),delimiter='\t'))
def fl(v,d=999):
    try: return float(v)
    except: return d
cand=[x for x in r if x['status']=='ok' and int(x['n_isolates'])>=3 and fl(x['plddt'])>=70 and (x['burial']=='buried' or fl(x['nearest_functional_site_A'])<=8) and rs[x['rv']]['struct_len']<=700]
kat=[x for x in r if x['rv']=='Rv1908c' and x['hgvs_p']=='p.Ser315Thr']
byg=collections.defaultdict(list)
for x in cand+kat: byg[x['rv']].append(x)
def work(rv):
    out='%s/foldx/%s.json'%(O,rv)
    if os.path.exists(out): return rv
    mp=rs[rv]; model=ms.load_model(mp['acc'])
    muts=[]; meta=[]
    for x in byg[rv]:
        mres=int(x['codon'])-(mp['offset']-1)
        muts.append((x['wt'],mres,x['mt'])); meta.append({'hgvs_p':x['hgvs_p'],'model_residue':mres,'key':'%s%s%d%s'%(x['wt'],'A',mres,x['mt'])})
    wd='/tmp/fx_%s'%rv
    res=ms.foldx_ddg(model,muts,wd)
    json.dump({'rv':rv,'meta':meta,'foldx':res},open(out,'w'),indent=1)
    return rv
genes=sorted(byg,key=lambda g:rs[g]['struct_len'])
print(len(genes),'genes',sum(len(v) for v in byg.values()),'variants',flush=True)
with cf.ThreadPoolExecutor(2) as ex:
    for g in ex.map(work,genes): print('done',g,flush=True)
