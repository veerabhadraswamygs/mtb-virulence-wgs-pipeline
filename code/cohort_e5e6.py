"""E5: global allele frequency (per lineage) of every local virulence-panel SNV, from all public assemblies (1,675: RefSeq and GenBank).
E6: local vs global panel-gene presence (read-depth presence in 55 isolates vs H37Rv-region coverage in public genomes).
Public assemblies are aligned to H37Rv; a position counts as 'covered' only inside an aligned block (>=1 kb, MAPQ>=20)."""
import json,glob,gzip,os,collections,bisect,re
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
O=P+'/cohort/e5e6'
def coarse(l):
    if re.match(r'lineage\d',l): return 'L'+re.match(r'lineage(\d)',l).group(1)
    if l.startswith('La'): return 'animal-adapted'
    return 'unassigned'
pub={}
for f in glob.glob(P+'/pub/calls/*.json'):
    d=json.load(open(f)); pub[d['acc']]=d
accs=sorted(pub); lin={a:coarse(pub[a]['lineage']) for a in accs}
starts={a:[s for s,e in pub[a]['intervals']] for a in accs}
def covered(a,p):
    i=bisect.bisect_left(starts[a],p)-1
    return i>=0 and pub[a]['intervals'][i][0]<p<=pub[a]['intervals'][i][1]
# local variants
loc=collections.defaultdict(lambda:{'iso':[],'info':None}); samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
loclin={}
for s in samples:
    r=json.load(open('%s/work/%s/results.json'%(P,s))); loclin[s]=coarse(r['lineage']['lineage'])
    for v in r['variants']:
        if v['type']=='snp' and v.get('in_virulence_catalogue') and not v.get('in_mask'):
            k=(v['pos'],v['ref'],v['alt']); loc[k]['iso'].append(s)
            loc[k]['info']=(v['rv'],v['gene'],v.get('hgvs_p',''),v.get('effect',''))
# public SNVs
snv=collections.defaultdict(list)
want={k[0] for k in loc}
for a in accs:
    for l in gzip.open(P+'/pub/calls/%s.vcf.gz'%a,'rt'):
        if l[0]=='#': continue
        f=l.split('\t',6)
        p=int(f[1])
        if p in want and len(f[3])==1 and len(f[4])==1 and f[4]!='.': snv[(p,f[3],f[4])].append(a)
LIN=sorted(set(lin.values()))
rows=[]
for k,v in sorted(loc.items()):
    p=k[0]; has=set(snv.get(k,[])); cov={a for a in accs if covered(a,p)}
    per={}
    for L in LIN:
        c=[a for a in cov if lin[a]==L]; n=len([a for a in has if lin[a]==L])
        per[L]={'covered':len(c),'with_alt':n,'freq':round(n/len(c),4) if c else None}
    tot=len([a for a in has if a in cov]); 
    ll=collections.Counter(loclin[s] for s in v['iso'])
    fr=tot/len(cov) if cov else None
    L1=per.get('L1',{}); 
    cls='not seen in public genomes' if tot==0 else 'rare globally (<1%)' if fr<0.01 else 'common (>=1%)'
    rows.append({'pos':p,'ref':k[1],'alt':k[2],'rv':v['info'][0],'gene':v['info'][1],'aa_change':v['info'][2],'effect':v['info'][3],
                 'n_local_isolates':len(v['iso']),'local_isolates':','.join(v['iso']),'local_lineages':dict(ll),
                 'public_covered':len(cov),'public_with_alt':tot,'public_freq':round(fr,4) if fr is not None else None,
                 'per_lineage':per,'class':cls,'only_in_MTB25_flagged_mixed_sample':v['iso']==['MTB25'],'lineage1_public_freq':L1.get('freq')})
# E6
panel={}
for b in glob.glob(P+'/ref/panel/*.bed'):
    for l in open(b):
        f=l.rstrip('\n').split('\t'); panel[f[3].split('_')[0]]=(int(f[1]),int(f[2]),f[3])
def frac_cov(a,s,e):
    tot=0
    for x,y in pub[a]['intervals']:
        lo=max(x,s); hi=min(y,e)
        if hi>lo: tot+=hi-lo
    return tot/(e-s)
glob_pres={}
for rv,(s,e,name) in panel.items():
    glob_pres[rv]=sum(frac_cov(a,s,e)>=0.9 for a in accs)/len(accs)
loc_pres=collections.defaultdict(int); nloc=0
for sm in samples:
    r=json.load(open('%s/work/%s/results.json'%(P,sm))); nloc+=1
    for g in r['catalogue_presence']:
        if g['call']=='present': loc_pres[g['rv']]+=1
e6=[]
for rv,(s,e,name) in sorted(panel.items()):
    lp=loc_pres.get(rv,0)/nloc; gp=glob_pres[rv]
    e6.append({'rv':rv,'gene':name,'local_present_frac':round(lp,3),'global_present_frac':round(gp,3),
               'flag':('local lower than global' if gp-lp>=0.2 else 'local higher than global' if lp-gp>=0.2 else '')})
json.dump({'n_public':len(accs),'public_lineage_counts':collections.Counter(lin.values()),'local_lineage_counts':collections.Counter(loclin.values()),'e5':rows,'e6':e6},open(O+'/e5e6.json','w'))
print(len(accs),'public;',collections.Counter(lin.values()));print(len(rows),'local virulence SNVs;',collections.Counter(r['class'] for r in rows));nm=[r for r in rows if not r['only_in_MTB25_flagged_mixed_sample']];print(len(nm),'excluding variants private to MTB25;',collections.Counter(r['class'] for r in nm));print(sum(1 for x in e6 if x['flag']),'E6 flagged genes')
