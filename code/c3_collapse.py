"""C3 step 2: collapse tips whose clade has a single identical codon sequence (exactly likelihood-preserving for reversible models; LRTs unchanged). Writes aln_c/<rv>.{phy,fas,nwk}"""
import os,sys,re,json
sys.path.insert(0,os.path.dirname(__file__))
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c3'; os.makedirs(O+'/aln_c',exist_ok=True)
def parse(s):
    s=s.strip().rstrip(';'); pos=0
    def node():
        nonlocal pos
        if s[pos]=='(':
            pos+=1; ch=[]
            while True:
                ch.append(node())
                if s[pos]==',': pos+=1; continue
                if s[pos]==')': pos+=1; break
            n={'ch':ch,'name':None,'len':0.0}
        else:
            st=pos
            while pos<len(s) and s[pos] not in ':,);': pos+=1
            n={'ch':[],'name':s[st:pos],'len':0.0}
        if pos<len(s) and s[pos]==':':
            pos+=1; st=pos
            while pos<len(s) and s[pos] not in ',);': pos+=1
            n['len']=float(s[st:pos])
        return n
    return node()
def newick(n,top=True):
    if not n['ch']: return '%s:%.8f'%(n['name'],max(n['len'],1e-7))
    return '('+','.join(newick(c,False) for c in n['ch'])+')'+('' if top else ':%.8f'%max(n['len'],1e-7))
def collapse(n,seq):
    """return (node, sig) ; sig = sequence string if all tips below identical else None"""
    if not n['ch']: return n,seq[n['name']]
    res=[collapse(c,seq) for c in n['ch']]
    n['ch']=[r[0] for r in res]
    sigs={r[1] for r in res}
    if len(sigs)==1 and None not in sigs:
        keep=res[0][0]; return {'ch':[],'name':keep['name'],'len':n['len']},res[0][1]
    return n,None
nt=0
for f in sorted(os.listdir(O+'/aln')):
    if not f.endswith('.nwk'): continue
    rv=f[:-4]
    seq={};cur=None
    for l in open('%s/aln/%s.fas'%(O,rv)):
        if l[0]=='>': cur=l[1:].strip(); seq[cur]=''
        else: seq[cur]+=l.strip()
    t=parse(open('%s/aln/%s'%(O,f)).read())
    t,_=collapse(t,seq)
    # root may have become single tip? handle
    tips=[]
    def lt(n):
        if not n['ch']: tips.append(n['name'])
        for c in n['ch']: lt(c)
    lt(t)
    if len(tips)<4:  # too few distinct haplotypes for a codeml tree (needs >=3 tips; unrooted trifurcation)
        open('%s/aln_c/%s.skip'%(O,rv),'w').write('%d tips after collapse\n'%len(tips)); continue
    # root must be multifurcating >=3 children; if bifurcating, unroot
    while len(t['ch'])==2:
        a,b=t['ch']
        if a['ch']: big,oth=a,b
        elif b['ch']: big,oth=b,a
        else: break
        oth['len']+=big['len']; t['ch']=big['ch']+[oth]
    L=len(next(iter(seq.values())))
    with open('%s/aln_c/%s.phy'%(O,rv),'w') as o:
        o.write(' %d %d\n'%(len(tips),L))
        for x in tips: o.write('%s  %s\n'%(x,seq[x]))
    with open('%s/aln_c/%s.fas'%(O,rv),'w') as o:
        for x in tips: o.write('>%s\n%s\n'%(x,seq[x]))
    open('%s/aln_c/%s.nwk'%(O,rv),'w').write(newick(t)+';\n'); nt+=1
print('collapsed alignments',nt, 'skipped',len([x for x in os.listdir(O+'/aln_c') if x.endswith('.skip')]))
