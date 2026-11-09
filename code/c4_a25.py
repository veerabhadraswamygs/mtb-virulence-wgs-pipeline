"""C4 analysis 25: predicted functional impact descriptors (substitution score BLOSUM62 + AlphaFold pLDDT, secondary structure, RSA, distance to functional sites) for all Tier 1a/1b missense variants seen in >=1 non-MTB25 isolate. No FoldX here."""
import json,os,re,collections,sys,csv
sys.path.insert(0,'/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/code')
import mod_struct as ms
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c4'
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
LIN={}
V={}
for s in samples:
    r=json.load(open('%s/work/%s/results.json'%(P,s))); LIN[s]=r['lineage']['lineage']
    for v in r['variants']:
        if v['type']=='snp' and not v.get('in_mask') and v.get('in_virulence_catalogue') and 'missense' in v['effect'] and v['panel_tier'] in ('1a','1b'):
            k=(v['rv'],v['pos'],v['alt'])
            d=V.setdefault(k,{'iso':set(),'v':v}); d['iso'].add(s)
rows=[]
for k,d in sorted(V.items(),key=lambda x:(x[0][0],x[0][1])):
    iso=sorted(d['iso']-{'MTB25'})
    if not iso: continue
    v=d['v']; wt,mt=v['aa_ref'],v['aa_alt']
    a=ms.analyse_variant(v['rv'],wt,int(v['codon']),mt,run_foldx=False)
    rs=a.get('solvent_accessibility') or {}
    fs=a.get('functional_sites') or {}
    near=min([x['distance_A'] for x in fs.values()],default=None) if isinstance(fs,dict) else None
    rows.append({'rv':v['rv'],'gene':v['gene'],'tier':v['panel_tier'],'pos':v['pos'],'hgvs_p':v['hgvs_p'],'wt':wt,'mt':mt,'codon':v['codon'],'n_isolates':len(iso),'isolates':','.join(iso),
        'lineage_defining':v['lineage_defining'],'blosum62':v['blosum62'],'in_tcell_epitope':v['in_tcell_epitope'],'domains':v['domains'],'status':a.get('status'),'plddt':a.get('plddt',''),'secondary_structure':a.get('secondary_structure',''),
        'rsa_pct':rs.get('rsa',''),'burial':rs.get('burial',''),'nearest_functional_site_A':near if near is not None else '','transmembrane':a.get('transmembrane_residue',''),'uniprot':a.get('uniprot','')})
keys=list(rows[0]); 
with open(O+'/a25_variant_structure_descriptors.tsv','w') as o:
    w=csv.DictWriter(o,fieldnames=keys,delimiter='\t'); w.writeheader(); w.writerows(rows)
print(len(rows),collections.Counter(r['status'] for r in rows))
