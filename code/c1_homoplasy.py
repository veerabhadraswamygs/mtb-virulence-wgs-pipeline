"""Analysis 23 (+E7): per-site homoplasy by Fitch parsimony on the rooted ML tree (root = M. canettii).
usage: c1_homoplasy.py <alignment prefix e.g. aln55> <tree prefix e.g. tree55>"""
import numpy as np,json,sys,re,glob,collections
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c1'
aln,tre=sys.argv[1],sys.argv[2]
EXCL={'MTB25'}     # 77.6% foreign reads: genotypes at shared sites not trusted for convergent-evolution inference
seq={};cur=None
for l in open('%s/%s.var.fasta'%(O,aln)):
    if l[0]=='>': cur=l[1:].strip(); seq[cur]=[]
    else: seq[cur].append(l.strip())
seq={k:''.join(v) for k,v in seq.items()}
sj=json.load(open('%s/%s.sites.json'%(O,aln))); sites=sj['sites']; ref=sj['ref']
code={'A':1,'C':2,'G':4,'T':8}
tip={k:np.array([code.get(c,15) for c in v],dtype=np.uint8) for k,v in seq.items()}
for k in EXCL:
    if k in tip: tip[k][:]=15
# newick
txt=open('%s/%s.treefile'%(O,tre)).read().strip()
pos=0
def parse():
    global pos
    if txt[pos]=='(':
        pos+=1; ch=[]
        while True:
            ch.append(parse())
            if txt[pos]==',': pos+=1; continue
            if txt[pos]==')': pos+=1; break
        m=re.match(r'[^:,;)]*',txt[pos:]); lab=m.group(0); pos+=len(lab)
        if txt[pos]==':':
            m=re.match(r':[0-9eE.+-]+',txt[pos:]); pos+=len(m.group(0))
        return {'ch':ch,'name':None,'support':lab}
    m=re.match(r'[^:,;)(]+',txt[pos:]); nm=m.group(0); pos+=len(nm)
    if txt[pos]==':':
        m=re.match(r':[0-9eE.+-]+',txt[pos:]); pos+=len(m.group(0))
    return {'ch':[],'name':nm}
root=parse()
def leaves(n): return [n['name']] if not n['ch'] else [x for c in n['ch'] for x in leaves(c)]
# root on the M. canettii edge
if 'Mcanettii' not in tip: raise SystemExit('no outgroup')
def find_path(n,t):
    if not n['ch']: return [n] if n['name']==t else None
    for c in n['ch']:
        p=find_path(c,t)
        if p: return [n]+p
if root['ch'] and any(c['name']=='Mcanettii' for c in root['ch']):
    rest=[c for c in root['ch'] if c['name']!='Mcanettii']
    inner={'ch':rest,'name':None,'support':''} if len(rest)>1 else rest[0]
    can=[c for c in root['ch'] if c['name']=='Mcanettii'][0]
    root={'ch':[can,inner],'name':None}
else: raise SystemExit('outgroup not at root; reroot needed')
n=len(sites); steps=np.zeros(n,np.int32)
def fitch(nd):
    if not nd['ch']: return tip[nd['name']]
    sets=[fitch(c) for c in nd['ch']]
    cur=sets[0]
    for s in sets[1:]:
        inter=cur&s
        empty=inter==0
        steps[:]+=empty
        cur=np.where(empty,cur|s,inter)
    return cur
fitch(root)
# minimum steps = distinct observed states - 1 (excluding missing)
nstates=np.zeros(n,np.int32)
M=np.array([tip[k] for k in tip])
for b in (1,2,4,8):
    nstates+=((M&b)!=0)&(M!=15)==True if False else (((M==b).any(axis=0)))
minsteps=np.maximum(nstates-1,0); homo=steps>minsteps
# genes / panel
G=json.load(open(P+'/ref/genes.json')); L=4411532
gid=np.full(L+2,'',dtype=object)
for rv,g in G.items():
    if g.get('type')=='CDS':
        gid[g['start']:g['end']+1]=rv
tier={}
for f in glob.glob(P+'/ref/panel/*.bed'):
    t=re.search(r'tier(\w+?)_',f).group(1)
    for l in open(f):
        x=l.split('\t'); tier[x[3].split('_')[0]]=t
rows=[]
for i,p in enumerate(sites):
    if nstates[i]<2: continue
    al=''.join(sorted({c for k in tip for c in 'ACGT' if tip[k][i]==code[c]}))
    rv=gid[p]
    rows.append({'pos':p,'ref':ref[i],'alleles':al,'steps':int(steps[i]),'min_steps':int(minsteps[i]),'homoplastic':bool(homo[i]),
                 'CI':round(float(minsteps[i]/steps[i]),3) if steps[i] else 1.0,'rv':rv,'gene':G.get(rv,{}).get('name','') if rv else '','panel_tier':tier.get(rv,'') if rv else ''})
import csv
with open('%s/a23_homoplasy_%s.tsv'%(O,tre),'w') as o:
    w=csv.DictWriter(o,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
hs=[r for r in rows if r['homoplastic']]
tot_steps=sum(r['steps'] for r in rows); tot_min=sum(r['min_steps'] for r in rows)
pan=[r for r in hs if r['panel_tier']]
summ={'tree':tre,'sites_analysed':len(rows),'homoplastic_sites':len(hs),'excluded_taxa_set_to_missing':sorted(EXCL),
      'consistency_index_overall':round(tot_min/tot_steps,4),
      'homoplastic_in_panel':len(pan),'by_tier':dict(collections.Counter(r['panel_tier'] for r in pan)),
      'top_panel':[(r['pos'],r['rv'],r['gene'],r['panel_tier'],r['steps'],r['alleles']) for r in sorted(pan,key=lambda r:-r['steps'])[:15]]}
json.dump(summ,open('%s/a23_summary_%s.json'%(O,tre),'w'),indent=1)
print(json.dumps(summ)[:1500])
