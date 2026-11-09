"""C4 analysis 26: benchmarking against public genomes. (i) accession list with country/lineage/in-global-tree flag; (ii) for every Tier 1a/1b missense variant seen locally (>=1 non-MTB25 isolate): public frequency overall, in L1, and in Indian public genomes (assemblies), plus FoldX ddG where run."""
import json,gzip,os,csv,glob,collections,re
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c4'
meta={}
import itertools
for l in itertools.chain(open(P+'/pub/ncbi_dataset/data/assembly_data_report.jsonl'),open(P+'/pub2/meta2/ncbi_dataset/data/assembly_data_report.jsonl')):
    j=json.loads(l); at={a.get('name'):a.get('value','') for a in j.get('assemblyInfo',{}).get('biosample',{}).get('attributes',[])}
    acc=j['accession']; meta[acc]={'assembly_name':j['assemblyInfo'].get('assemblyName',''),'level':j['assemblyInfo'].get('assemblyLevel',''),'biosample':j['assemblyInfo'].get('biosample',{}).get('accession',''),'country':(at.get('geo_loc_name','') or '').split(':')[0].strip(),'collection_date':at.get('collection_date',''),'host':at.get('host','')}
sel=set(json.load(open(O+'/c4_global_selection.json'))['selected_public'])
sel300=set(json.load(open(P+'/pub/l1_selection.json'))['selected_300'])
calls={}
rows=[]
for acc in sorted(meta):
    f=P+'/pub/calls/%s.json'%acc
    if not os.path.exists(f): continue
    j=json.load(open(f)); calls[acc]=j
    rows.append({'accession':acc,**meta[acc],'lineage':j.get('lineage'),'covered_bp_vs_H37Rv':j.get('covered_bp'),'in_global_all_lineage_tree':acc in sel,'in_L1_context_tree':acc in sel300,'in_L1_tbd1_analysis':str(j.get('lineage','')).startswith('lineage1')})
with open(O+'/a26_public_accession_list.tsv','w') as o:
    w=csv.DictWriter(o,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
ind=[a for a in calls if meta[a]['country'].lower()=='india']
print('public total',len(calls),'India',len(ind),collections.Counter(re.sub(r'(lineage\d).*',r'\1',str(calls[a]['lineage'])) for a in ind))
# India variants
V=list(csv.DictReader(open(O+'/a25_variant_structure_descriptors.tsv'),delimiter='\t'))
need={(int(v['pos']),) for v in V}
altat=collections.defaultdict(dict); cov={}
import numpy as np
for a in ind:
    snp={}
    for l in gzip.open(P+'/pub/calls/%s.vcf.gz'%a,'rt'):
        if l[0]=='#': continue
        f=l.split('\t',6)
        if len(f[3])==1 and len(f[4])==1: snp[int(f[1])]=f[4]
    altat[a]=snp
    cov[a]=calls[a]['intervals']
def covered(a,p):
    import bisect
    iv=cov[a]; i=bisect.bisect_right([x[0] for x in iv],p-1)-1
    return i>=0 and iv[i][0]<p<=iv[i][1]
e5={(x['pos'],x['alt']):x for x in json.load(open(P+'/cohort/e5e6/e5e6.json'))['e5']}
fx={}
for f in glob.glob(O+'/foldx/*.json'):
    j=json.load(open(f))
    for m in j['meta']:
        r=j['foldx'].get(m['key'])
        if isinstance(r,dict) and 'ddG_kcal_mol' in r: fx[(j['rv'],m['hgvs_p'])]=r
out=[]
for v in V:
    if v['status']!='ok' and v['status'] is None: continue
    pos=int(v['pos']); e=e5.get((pos,v['mt'] and None)) 
    # find e5 record by pos and nucleotide alt: e5 keyed by nucleotide alt; recover from effects (match pos + gene + aa)
    cand=[x for (p,a),x in e5.items() if p==pos and x['rv']==v['rv']]
    e=cand[0] if cand else None
    ic=[a for a in ind if covered(a,pos)]
    ia=[a for a in ic if e and altat[a].get(pos)==e['alt']]
    out.append({'rv':v['rv'],'gene':v['gene'],'tier':v['tier'],'hgvs_p':v['hgvs_p'],'pos':pos,'n_local_isolates':v['n_isolates'],'lineage_defining':v['lineage_defining'],
        'public_covered':e['public_covered'] if e else '','public_with_alt':e['public_with_alt'] if e else '','public_freq':e['public_freq'] if e else '','public_L1_freq':e['per_lineage']['L1']['freq'] if e else '','public_class':e['class'] if e else '',
        'India_covered':len(ic),'India_with_alt':len(ia),'India_freq':round(len(ia)/len(ic),3) if ic else '','pLDDT':v['plddt'],'burial':v['burial'],
        'FoldX_ddG_kcal_mol':fx.get((v['rv'],v['hgvs_p']),{}).get('ddG_kcal_mol',''),'FoldX_sd':fx.get((v['rv'],v['hgvs_p']),{}).get('sd','')})
with open(O+'/a26_variant_public_frequency.tsv','w') as o:
    w=csv.DictWriter(o,fieldnames=list(out[0]),delimiter='\t'); w.writeheader(); w.writerows(out)
print(len(out),collections.Counter(x['public_class'] for x in out))
