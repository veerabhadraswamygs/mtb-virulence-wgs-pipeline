"""C6 analysis 25 (PROVEAN): PROVEAN 1.1.5 (patched for non-gi sequence IDs) against a mycobacterial UniRef90 database (UniProt taxonomy 1763, 309,177 clusters),
for all Tier 1a/1b missense variants seen in >=1 non-MTB25 isolate (analysis 25 list). One run per gene, all that gene's variants."""
import csv,json,os,subprocess,collections,re,sys
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
BL=P+'/src/ncbi-blast-2.17.0+/bin'
env=dict(os.environ,PATH='/tmp/claude-0/prov/inst/bin:%s:%s/bin:%s'%(BL,P,os.environ['PATH']))
W='/mnt/user-data/outputs/_provean/run'; os.makedirs(W,exist_ok=True)
h=json.load(open(P+'/ref/h37rv_prot.json'))
r=list(csv.DictReader(open(P+'/cohort/c4/a25_variant_structure_descriptors.tsv'),delimiter='\t'))
by=collections.defaultdict(list)
for x in r:
    s=h.get(x['rv'])
    if not s: continue
    c=int(x['codon'])
    if 1<c<=len(s) and s[c-1]==x['wt'] and x['mt'] not in ('*',''): by[x['rv']].append(x)
out=P+'/cohort/c4/a25_provean_scores.tsv'; res={}
for i,(rv,xs) in enumerate(sorted(by.items())):
    d='%s/%s'%(W,rv); os.makedirs(d,exist_ok=True)
    open(d+'/q.fa','w').write('>%s\n%s\n'%(rv,h[rv]))
    vs=sorted({'%s%s%s'%(x['wt'],x['codon'],x['mt']) for x in xs}); open(d+'/q.var','w').write('\n'.join(vs)+'\n')
    p=subprocess.run(['provean.sh','-q',d+'/q.fa','-v',d+'/q.var','--num_threads','2','--tmp_dir',d],capture_output=True,text=True,env=env,cwd=d)
    txt=p.stdout; open(d+'/out.txt','w').write(txt+p.stderr)
    m=re.search(r'Number of supporting sequences used:\s*(\d+)',txt); n=int(m.group(1)) if m else 0
    cl=re.search(r'Number of clusters:\s*(\d+)',txt)
    for l in txt.split('\n'):
        mm=re.match(r'^([A-Z]\d+[A-Z])\t(-?[\d.]+)$',l.strip())
        if mm: res[(rv,mm.group(1))]=(float(mm.group(2)),n,int(cl.group(1)) if cl else 0)
    print(i,rv,len(vs),n,flush=True)
rows=[]
for x in r:
    k=(x['rv'],'%s%s%s'%(x['wt'],x['codon'],x['mt']))
    v=res.get(k)
    rows.append({'rv':x['rv'],'gene':x['gene'],'tier':x['tier'],'hgvs_p':x['hgvs_p'],'n_isolates':x['n_isolates'],
      'PROVEAN_score':v[0] if v else '','supporting_sequences':v[1] if v else '','supporting_clusters':v[2] if v else '',
      'PROVEAN_call_threshold_-2.5':('deleterious' if v and v[0]<=-2.5 else 'neutral') if v else 'not scored'})
w=csv.DictWriter(open(out,'w'),fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
print('scored',len(res),'of',len(r))
