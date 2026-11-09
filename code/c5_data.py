"""C5 step 1: analysis dataset (isolate x variant, isolate x gene, clinical variables). No patient names / phones / addresses are read into any output."""
import json,os,re,collections,pickle,numpy as np,openpyxl
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c5'; os.makedirs(O,exist_ok=True)
XL='/mnt/user-data/uploads/Downloads/Dr. Ranjitha TB-Pt final list.xlsx'
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
R={s:json.load(open('%s/work/%s/results.json'%(P,s))) for s in samples}
def coarse(l):
    m=re.match(r'lineage(\d)',l); return 'L'+m.group(1) if m else 'unassigned'
LIN={s:coarse(R[s]['lineage']['lineage']) for s in samples}
# ---- clinical (only the columns used; identifiers dropped)
ws=openpyxl.load_workbook(XL,data_only=True)['Sheet1']; rows=list(ws.iter_rows(values_only=True)); H=[str(x).strip() for x in rows[0]]
norm=lambda s:re.sub(r'\s+','',str(s)).upper()
M={norm(r[0]):dict(zip(H,r)) for r in rows[1:] if r[0]}
def g(s,c):
    v=M[norm(s)][c]; return None if v is None else str(v).strip()
clin={}
for s in samples:
    hab=g(s,'ADDITIVE HABITS'); smear=g(s,'SPUTUM SMEAR GRADING')
    clin[s]={'diabetes':g(s,'DIABETES MELLITUS'),'age':float(M[norm(s)]['AGE']),'sex':g(s,'GENDER'),
      'smoking':'Yes' if hab and 'Smoker' in hab else 'No','alcohol':'Yes' if hab and 'Alcoholic' in hab else 'No','additive_habits_raw':hab,
      'ses':g(s,'SOCIOECONOMIC STATUS'),'smear':smear,'smear_3plus':'Yes' if smear=='3+' else 'No',
      'prev_tb':g(s,'PREVIOUS  HISTORY OF TB'),'outcome':g(s,'CLINICAL OUTCOME'),'died':'Yes' if g(s,'CLINICAL OUTCOME')=='Dead' else 'No',
      'cough':g(s,'PERSISTENT COUGH WITH EXPECTORATION'),'fever':g(s,'FEVER'),'weight_loss':g(s,'WEIGHT LOSS'),'chest_pain':g(s,'CHEST PAIN'),
      'breathlessness':g(s,'BREATHLESSNESS'),'haemoptysis':g(s,'HAEMOPTYSIS'),'cbnaat':g(s,'CBNAAT RESULTS')}
    clin[s]['sex']=clin[s]['sex'].strip()
G=json.load(open(P+'/ref/genes.json')); cat=json.load(open(P+'/ref/catalogue.json'))
# ---- callable arrays
cal={s:np.unpackbits(np.load('%s/cohort/c1/cal_%s.npy'%(P,s)))[:4411533].astype(bool) for s in samples}
mask=np.zeros(4411533,bool)
for l in open(P+'/ref/mask_static.bed'):
    f=l.split('\t'); mask[int(f[1])+1:int(f[2])+1]=True
NS={'missense_variant','stop_gained','start_lost','stop_lost','stop_lost&splice_region_variant','frameshift_variant','disruptive_inframe_deletion','disruptive_inframe_insertion','inframe_deletion','inframe_insertion'}
SNV_NS={'missense_variant','stop_gained','start_lost','stop_lost','stop_lost&splice_region_variant'}
# Tier 1 gene list from per-isolate catalogue_presence
info={x['rv']:x for x in R[samples[0]]['catalogue_presence']}
t1={rv:x['tier'] for rv,x in info.items() if x['tier'] in('1a','1b')}
print('tier1 genes',collections.Counter(t1.values()))
# ---- variants
V={}   # key -> dict
for s in samples:
    for v in R[s]['variants']:
        if not v.get('in_virulence_catalogue') or v.get('in_mask'): continue
        rv=v['rv']
        if rv not in t1: continue
        if not (v['effect'] in NS or v['effect']=='synonymous_variant' or 'upstream' in v['effect']): pass
        k=(v['pos'],v['ref'],v['alt'])
        d=V.setdefault(k,{'pos':v['pos'],'ref':v['ref'],'alt':v['alt'],'type':v['type'],'rv':rv,'gene':v['gene'],'tier':t1[rv],'effect':v['effect'],'hgvs_p':v['hgvs_p'],'hgvs_c':v['hgvs_c'],
                          'lineage_defining':v.get('lineage_defining',''),'in_tcell_epitope':v.get('in_tcell_epitope',False),'isolates':[]})
        d['isolates'].append(s)
print('tier1 variant records',len(V))
# per-isolate nonsynonymous burden
tc={s:R[s]['tier_callable'] for s in samples}
def burden(s,tier,defn):
    n=0
    for k,d in V.items():
        if s not in d['isolates']: continue
        if tier!='1' and d['tier']!=tier: continue
        if d['effect'] in (SNV_NS if defn=='snv' else NS): n+=1
    return n
# variant-sharing note: counts are per isolate (all alleles vs reference), including lineage-defining variants
B={}
for s in samples:
    for tier in('1','1a','1b'):
        cbp=(tc[s]['1a']['callable_bp'] if tier=='1a' else tc[s]['1b']['callable_bp'] if tier=='1b' else tc[s]['1a']['callable_bp']+tc[s]['1b']['callable_bp'])
        for defn in('snv','all'):
            B[(s,tier,defn)]={'n':burden(s,tier,defn),'callable_bp':cbp}
pickle.dump({'samples':samples,'LIN':LIN,'clin':clin,'V':V,'B':B,'t1':t1,'info':info,
             'pres':{s:{x['rv']:x['call'] for x in R[s]['catalogue_presence']} for s in samples},
             'lowcov':[s for s in samples if R[s]['exclusion']['genome_breadth_10x_percent']<90],
             'breadth':{s:R[s]['exclusion']['genome_breadth_10x_percent'] for s in samples},
             'epi':{s:R[s]['epitope_partition'] for s in samples}},open(O+'/c5_data.pkl','wb'))
# callable status of each variant position in each isolate (missing if position not callable)
miss={}
for k,d in V.items():
    p=d['pos']; miss[k]=[s for s in samples if not cal[s][p]]
pickle.dump(miss,open(O+'/c5_missing.pkl','wb'))
print('variants with any missing isolate',sum(1 for k in miss if miss[k]),'of',len(miss))
print(collections.Counter(len(x) for x in miss.values()).most_common(8))
print('lowcov',[s for s in samples if R[s]['exclusion']['genome_breadth_10x_percent']<90])
print(collections.Counter((clin[s]['diabetes'],LIN[s]) for s in samples))
