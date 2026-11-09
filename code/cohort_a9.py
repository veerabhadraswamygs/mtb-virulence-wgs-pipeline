"""Analysis 9: whole-genome pan-genome. CD-HIT protein clustering (>=90% identity, >=80% length both ways) over the proteomes of all
isolates with a gene-content-assessable assembly plus the H37Rv reference proteome. Isolates that are not assessable
(see a8_summary.json) are listed, not used to call gene absence."""
import json,subprocess,os,collections
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
O=P+'/cohort/a9'; os.makedirs(O,exist_ok=True)
S=json.load(open(P+'/cohort/a8/a8_summary.json'))
ok=[x['sample'] for x in S if x['assessable_for_gene_content']]
with open(O+'/all.faa','w') as o:
    for l in open(P+'/cohort/a8/H37Rv_proteome.faa'): o.write(l.replace('>','>REF|',1) if l.startswith('>') else l)
    for s in ok:
        for l in open('%s/work/%s/%s_proteome_p2.faa'%(P,s,s)):
            if l.startswith('>'): o.write('>%s|%s\n'%(s,l[1:].split(' ')[0].split('_',1)[1] if False else l[1:].split(' ')[0]))
            elif l.strip(): o.write(l.rstrip('*\n').replace('*','')+'\n')
subprocess.run('%s/bin/cd-hit -i %s/all.faa -o %s/clusters -c 0.9 -aS 0.8 -aL 0.8 -n 5 -M 3000 -T 2 -d 0 > %s/cdhit.log 2>&1'%(P,O,O,O),shell=True,check=True)
cl=collections.defaultdict(list); cur=None
for l in open(O+'/clusters.clstr'):
    if l.startswith('>Cluster'): cur=int(l.split()[1])
    else:
        name=l.split('>')[1].split('...')[0]; cl[cur].append(name)
n=len(ok); rows=[]
G=json.load(open(P+'/ref/genes.json'))
for c,m in cl.items():
    iso=sorted({x.split('|')[0] for x in m if not x.startswith('REF|')}); ref=[x[4:] for x in m if x.startswith('REF|')]
    f=len(iso)/n
    cat='core' if f>=0.99 else 'soft-core' if f>=0.95 else 'shell' if f>=0.15 else 'cloud'
    rows.append({'cluster':c,'H37Rv_gene':ref[0] if ref else '','name':G.get(ref[0],{}).get('name','') if ref else '','n_isolates':len(iso),'frac':round(f,3),'class':cat,'isolates':','.join(iso)})
json.dump({'n_isolates_used':n,'isolates_used':ok,'not_assessable':[x['sample'] for x in S if not x['assessable_for_gene_content']],
           'n_clusters':len(rows),'by_class':collections.Counter(r['class'] for r in rows),'rows':rows},open(O+'/a9_pangenome.json','w'))
print(n,'isolates;',len(rows),'clusters;',collections.Counter(r['class'] for r in rows))
