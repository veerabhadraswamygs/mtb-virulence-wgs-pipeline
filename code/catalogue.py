"""Freeze the virulence gene catalogue (brief section 4).
Sources: VFDB (setA+setB), Mycobrowser functional category, BV-BRC (both spellings),
and the brief's appendix curated locus list. Deduplicated by Rv number.
"""
import json,csv,re,collections,sys
R='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/ref'
U='/mnt/user-data/uploads/Downloads/limitation fix files'
genes=json.load(open(R+'/genes.json'))
vfdb=json.load(open(R+'/vfdb_rv.json'))

name2rv={}; lname={}
for rv,g in genes.items():
    n=(g['name'] or '').lower()
    if n: name2rv.setdefault(n,rv); lname[rv]=n
    else: lname[rv]=rv.lower()
def by_name(n):
    return name2rv.get(n.lower())
def by_prefix(pfx):
    pfx=pfx.lower()
    return [rv for rv,n in lname.items() if n.startswith(pfx)]
def rv_range(a,b):
    """Rv0169..Rv0174 style inclusive numeric range, same suffix handling."""
    m1=re.match(r'Rv(\d+)([c]?)$',a); m2=re.match(r'Rv(\d+)([c]?)$',b)
    lo,hi=int(m1.group(1)),int(m2.group(1)); suf=m1.group(2)
    out=[]
    for rv in genes:
        m=re.match(r'Rv(\d+)([c]?)$',rv)
        if m and lo<=int(m.group(1))<=hi: out.append(rv)
    return out

# ---------- appendix curated list (brief appendix, verbatim groups) ----------
GROUPS=[
 ('TbD1 region',                          dict(exact=['Rv2073c','Rv2074'])),
 ('ESX / type VII secretion',             dict(names=['esxA','esxB','espA','espC','espD','espR','whiB6'],
                                               exact=['Rv3849','Rv3862c'],
                                               prefixes=['esx','ecc','esp','mycP'])),
 ('PE / PPE / PE_PGRS',                   dict(names=['PE_PGRS33','PPE38','PPE68','PE35'],
                                               exact=['Rv1818c','Rv2352c','Rv3873'],
                                               category=['PE/PPE'])),
 ('Cell envelope and lipid biosynthesis', dict(names=['pks15','pks1','mas','fadD26','fadD28','drrA','drrB','drrC',
                                                      'mmpL3','mmpL7','mmpL11','papA5','tesA','pcaA','fbpA'],
                                               exact=['Rv2962c'],
                                               prefixes=['pks','pps','mmpl','mmps','drr'])),
 ('mce operons',                          dict(prefixes=['mce','yrbe'])),
 ('Phagosome maturation arrest',          dict(names=['ptpA','ptpB','ndkA','sapM','eis','pknG','ndk'])),
 ('Oxidative and nitrosative stress',     dict(names=['katG','sodA','sodC','ahpC','ahpD','mel2','Rv2466c'],
                                               exact=['Rv2466c'])),
 ('Regulatory systems',                   dict(names=['phoP','phoR','sigA','sigC','sigE','sigF','sigH','sigK',
                                                      'whiB3','mprA','mprB','dosR','dosS','dosT','devR','devS'],
                                               exact=['Rv0757','Rv0758'])),
 ('Iron acquisition',                     dict(names=['irtA','irtB'], prefixes=['mbt','irt'])),
 ('Metabolic adaptation and persistence', dict(names=['icl1','icl2','icl','aceA','pckA','relA','hspX'])),
 ('Toxin-antitoxin systems',              dict(prefixes=['vap','maz','relb','rele','rela','higa','higb','parde','pare','pard'])),
 ('Antigens and epitope carriers',        dict(names=['fbpA','fbpB','fbpC','esxH','mpt64','mpt63','lprG'],
                                               exact=['Rv3804c','Rv1886c','Rv0129c','Rv0309','Rv2945c'])),
]
appendix={}
detail=collections.OrderedDict()
for gname,rule in GROUPS:
    s=set()
    for n in rule.get('names',[]):
        rv=by_name(n)
        if rv: s.add(rv)
    for rv in rule.get('exact',[]):
        if rv in genes: s.add(rv)
    for p in rule.get('prefixes',[]):
        s.update(by_prefix(p))
    for c in rule.get('category',[]):
        s.update(rv for rv,g in genes.items() if g['category']==c)
    s={rv for rv in s if genes.get(rv,{}).get('type')=='CDS'}
    detail[gname]=sorted(s)
    for rv in s: appendix.setdefault(rv,gname)
print('appendix loci:',len(appendix))
for k,v in detail.items(): print('   %-40s %4d'%(k,len(v)))

# ---------- BV-BRC ----------
bvbrc={}
for r in csv.DictReader(open(U+'/BVBRC_sp_gene.csv',encoding='utf8',errors='replace')):
    if r['Property'] not in ('Virulence Factor','Virulance factor'): continue
    rv=(r['RefSeq Locus Tag'] or '').strip()
    if not rv: rv=by_name((r['Gene'] or '').strip()) or ''
    if not rv or rv not in genes: continue
    e=bvbrc.setdefault(rv,{'sources':set(),'evidence':set()})
    e['sources'].add(r['Source']); e['evidence'].add(r['Evidence'])
print('BV-BRC virulence loci mapped:',len(bvbrc))

# ---------- Mycobrowser category 1 ----------
myco={rv for rv,g in genes.items() if g['category']=='virulence, detoxification, adaptation' and g['type']=='CDS'}
print('Mycobrowser cat-1 loci:',len(myco))
print('VFDB loci:',len(vfdb))

# ---------- merge ----------
cat={}
for rv in set(list(vfdb)+list(bvbrc)+list(myco)+list(appendix)):
    g=genes.get(rv)
    if not g or g['type']!='CDS': continue
    src=[]
    if rv in vfdb: src.append('VFDB')
    if rv in myco: src.append('Mycobrowser')
    if rv in bvbrc: src.append('BV-BRC')
    if rv in appendix: src.append('Appendix')
    cat[rv]={'rv':rv,'name':g['name'] or rv,'product':g['product'] or '','start':g['start'],'end':g['end'],
             'strand':g['strand'],'mycobrowser_category':g['category'] or '',
             'system':appendix.get(rv,''),'sources':src,
             'vfdb_vf':vfdb.get(rv,{}).get('vf',''),'vfdb_vfid':vfdb.get(rv,{}).get('vfid',''),
             'vfdb_category':vfdb.get(rv,{}).get('vfcat',''),
             'vfdb_sets':','.join(vfdb.get(rv,{}).get('sets',[])),
             'bvbrc_sources':','.join(sorted(bvbrc.get(rv,{}).get('sources',[]))) if rv in bvbrc else ''}
# assign a functional system to loci that only came from VFDB/BV-BRC/Mycobrowser
VFC2SYS={'Effector delivery system':'ESX / type VII secretion','Immune modulation':'Cell envelope and lipid biosynthesis',
         'Nutritional/Metabolic factor':'Metabolic adaptation and persistence','Regulation':'Regulatory systems',
         'Stress survival':'Oxidative and nitrosative stress','Adherence':'Cell envelope and lipid biosynthesis',
         'Exotoxin':'Toxin-antitoxin systems','Enzyme':'Other / unclassified','Others':'Other / unclassified',
         'Antimicrobial activity/Competitive advantage':'Other / unclassified'}
for rv,e in cat.items():
    if not e['system']:
        e['system']=VFC2SYS.get(e['vfdb_category'],'') or 'Other / unclassified'
print('\nCATALOGUE TOTAL:',len(cat))
print('by source combination:')
for k,v in collections.Counter(','.join(e['sources']) for e in cat.values()).most_common(): print('   %4d  %s'%(v,k))
print('by functional system:')
for k,v in collections.Counter(e['system'] for e in cat.values()).most_common(): print('   %4d  %s'%(v,k))
bp=sum(e['end']-e['start']+1 for e in cat.values())
print('total catalogue length: %d bp (%.2f%% of genome)'%(bp,100*bp/4411532))
json.dump(cat,open(R+'/catalogue.json','w'))
