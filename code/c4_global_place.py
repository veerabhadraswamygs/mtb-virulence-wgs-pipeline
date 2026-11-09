import re,json,csv,collections,os,sys
sys.setrecursionlimit(10000)
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c4'
s=open(O+'/tree/treeG.treefile').read().strip().rstrip(';')
pos=0
def node():
    global pos
    if s[pos]=='(':
        pos+=1; ch=[]
        while True:
            ch.append(node())
            if s[pos]==',': pos+=1; continue
            if s[pos]==')': pos+=1; break
        while pos<len(s) and s[pos] not in ':,);': pos+=1
        n={'ch':ch,'name':None,'len':0.0}
    else:
        st=pos
        while s[pos] not in ':,);': pos+=1
        n={'ch':[],'name':s[st:pos],'len':0.0}
    if pos<len(s) and s[pos]==':':
        pos+=1; st=pos
        while pos<len(s) and s[pos] not in ',);': pos+=1
        n['len']=float(s[st:pos])
    return n
t=node()
tips={}
def walk(n,par,d):
    n['par']=par; n['d']=d+n['len']
    if not n['ch']: tips[n['name']]=n
    for c in n['ch']: walk(c,n,n['d'])
walk(t,None,0)
def anc(n):
    r=[]
    while n: r.append(n); n=n['par']
    return r
def dist(a,b):
    A=anc(a); sa={id(x) for x in A}
    for x in anc(b):
        if id(x) in sa: return a['d']+b['d']-2*x['d']
local=[k for k in tips if k.startswith('MTB')]; pub=[k for k in tips if k.startswith('PUB_')]
meta={r['accession']:r for r in csv.DictReader(open(O+'/a26_public_accession_list.tsv'),delimiter='\t')}
LIN={}
for k in local: LIN[k]=json.load(open('%s/work/%s/results.json'%(P,k)))['lineage']['lineage']
rows=[]
for k in local:
    best=min(pub,key=lambda p:dist(tips[k],tips[p])); acc=best[4:]
    # nearest local too
    bl=min([x for x in local if x!=k],key=lambda p:dist(tips[k],tips[p]))
    rows.append({'isolate':k,'local_lineage':LIN[k],'nearest_public':acc,'public_lineage':meta[acc]['lineage'],'public_country':meta[acc]['country'],'patristic_distance_subs_per_site':round(dist(tips[k],tips[best]),6),'approx_SNV_difference':round(dist(tips[k],tips[best])*3863041),'nearest_local':bl,'nearest_local_distance_SNV':round(dist(tips[k],tips[bl])*3863041)})
w=csv.DictWriter(open(O+'/a26_local_isolates_nearest_public_genome.tsv','w'),fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
# lineage concordance with nearest public (coarse)
co=lambda l: re.match(r'lineage(\d)',str(l)).group(1) if re.match(r'lineage(\d)',str(l)) else 'x'
agree=sum(1 for r in rows if co(r['local_lineage'])==co(r['public_lineage'])); known=sum(1 for r in rows if co(r['local_lineage'])!='x')
print('concordant coarse lineage with nearest public:',agree,'of',known)
print(collections.Counter(r['public_country'] for r in rows).most_common(6))
print(sorted(r['approx_SNV_difference'] for r in rows))
# monophyly of public lineages in tree
def mono(names):
    nodes=[tips[x] for x in names]
    # find MRCA
    A=anc(nodes[0]); 
    for x in nodes[1:]:
        sa={id(y) for y in anc(x)}; A=[y for y in A if id(y) in sa]
    m=A[0]
    def leaves(n): return [n['name']] if not n['ch'] else [l for c in n['ch'] for l in leaves(c)]
    L=set(leaves(m)); return L<=set(names)
sets=collections.defaultdict(list)
for k in tips:
    if k.startswith('PUB_'): sets['pub_'+co(meta[k[4:]]['lineage'])].append(k)
    elif k in LIN: sets['local_'+co(LIN[k])].append(k)
for lab in ('1','2','3','4'):
    names=[k for k in tips if (k.startswith('PUB_') and co(meta[k[4:]]['lineage'])==lab) or (k in LIN and co(LIN[k])==lab)]
    print('L'+lab,len(names),'monophyletic (local+public):',mono(names))
for lab in ('1','4'):
    names=[k for k in tips if (k.startswith('PUB_') and co(meta[k[4:]]['lineage'])==lab) or (k in LIN and co(LIN[k])==lab)]
    nodes=[tips[x] for x in names]; A=anc(nodes[0])
    for x in nodes[1:]:
        sa={id(y) for y in anc(x)}; A=[y for y in A if id(y) in sa]
    def leaves(n): return [n['name']] if not n['ch'] else [l for c in n['ch'] for l in leaves(c)]
    intr=[l for l in leaves(A[0]) if l not in names]
    print('L'+lab,'intruders in MRCA clade:',[(l, LIN.get(l) or meta[l[4:]]['lineage']) for l in intr])
