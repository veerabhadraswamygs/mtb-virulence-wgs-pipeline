"""C2 cohort tables from the 55 frozen per-isolate results (analyses 16, 17, 20, 21 cohort part, E11, 27, 29, 30)."""
import json,os,re,csv,collections,numpy as np
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c2'
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
R={s:json.load(open('%s/work/%s/results.json'%(P,s))) for s in samples}
def coarse(l):
    m=re.match(r'lineage(\d)',l); return 'L'+m.group(1) if m else 'unassigned'
LIN={s:coarse(R[s]['lineage']['lineage']) for s in samples}
cohort_status={s:R[s]['exclusion'].get('cohort_status','') for s in samples}
def tsv(name,rows,cols):
    with open('%s/%s.tsv'%(O,name),'w') as o:
        w=csv.writer(o,delimiter='\t'); w.writerow(cols)
        for r in rows: w.writerow([r.get(c,'') if isinstance(r,dict) else r for c in cols] if isinstance(r,dict) else r)
G=json.load(open(P+'/ref/genes.json'))
out={}
# ---- 16 presence/absence matrix
genes=[g['rv'] for g in R[samples[0]]['catalogue_presence']]
info={g['rv']:g for g in R[samples[0]]['catalogue_presence']}
with open(O+'/a16_presence_matrix.tsv','w') as o:
    w=csv.writer(o,delimiter='\t'); w.writerow(['rv','gene','tier','system']+samples)
    for rv in genes:
        row=[]
        for s in samples:
            d={g['rv']:g for g in R[s]['catalogue_presence']}[rv] if False else None
        w.writerow([rv,info[rv]['gene'],info[rv]['tier'],info[rv]['system']])
pres={s:{g['rv']:g['call'] for g in R[s]['catalogue_presence']} for s in samples}
LOWCOV=[s for s in samples if R[s]['exclusion']['genome_breadth_10x_percent']<90]
for s_ in LOWCOV:
    for rv in pres[s_]:
        if pres[s_][rv]=='absent': pres[s_][rv]='not determined'   # cohort rule: no gene is called absent in an isolate with <90% of the genome at >=10x
out['a16_low_coverage_isolates_absent_calls_set_to_not_determined']=LOWCOV
with open(O+'/a16_presence_matrix.tsv','w') as o:
    w=csv.writer(o,delimiter='\t'); w.writerow(['rv','gene','tier','system']+samples)
    for rv in genes: w.writerow([rv,info[rv]['gene'],info[rv]['tier'],info[rv]['system']]+[pres[s][rv] for s in samples])
# ---- 17 pan-virulome
a8=json.load(open(P+'/cohort/a8/a8_summary.json')); ass={x['sample'] for x in a8 if x['assessable_for_gene_content']}
found={}
for s in ass: found[s]=set(json.load(open('%s/cohort/a8/%s.found.json'%(P,s)))['found'])
rows=[]
for rv in genes:
    calls=[pres[s][rv] for s in samples]; det=[c for c in calls if c!='not determined']
    n_det=len(det); n_pres=sum(c=='present' for c in det)
    f=n_pres/n_det if n_det else None
    cls='not classifiable' if n_det<30 else ('core' if f>=0.99 else 'soft-core' if f>=0.95 else 'shell' if f>=0.15 else 'cloud')
    af=sum(rv in found[s] for s in ass)/len(ass)
    rows.append({'rv':rv,'gene':info[rv]['gene'],'tier':info[rv]['tier'],'system':info[rv]['system'],'n_determined':n_det,'n_present':n_pres,
                 'frac_present_reads':round(f,3) if f is not None else '','class_reads':cls,'frac_present_assembly_49':round(af,3),
                 'absent_isolates':','.join(s for s in samples if pres[s][rv]=='absent')})
tsv('a17_pan_virulome',rows,list(rows[0]))
out['a17']=dict(collections.Counter(r['class_reads'] for r in rows)); out['a17_non_core']=[(r['rv'],r['gene'],r['class_reads'],r['n_present'],r['n_determined']) for r in rows if r['class_reads'] not in('core','not classifiable')][:60]
# ---- 20 lineage-stratified allele profiles (panel variants)
vv=collections.defaultdict(lambda:{'iso':[],'info':None})
for s in samples:
    for v in R[s]['variants']:
        if v['type']=='snp' and v.get('in_virulence_catalogue') and not v.get('in_mask'):
            k=(v['pos'],v['ref'],v['alt']); vv[k]['iso'].append(s); vv[k]['info']=v
nl=collections.Counter(LIN.values()); rows=[]
for k,d in sorted(vv.items()):
    v=d['info']; c=collections.Counter(LIN[s] for s in d['iso'])
    big=[l for l in c if c[l]>=1]
    restricted=''
    if len(c)==1 and len(d['iso'])>=3: restricted='restricted to %s'%list(c)[0]
    only25=(d['iso']==['MTB25'])
    if only25: restricted=''
    row={'only_in_MTB25_flagged':only25,'pos':k[0],'ref':k[1],'alt':k[2],'rv':v['rv'],'gene':v['gene'],'tier':v['panel_tier'],'effect':v['effect'],'hgvs_p':v['hgvs_p'],'n_isolates':len(d['iso']),'lineage_restricted':restricted,'in_tcell_epitope':v['in_tcell_epitope'],'lineage_defining':v['lineage_defining']}
    for l in sorted(nl): row['freq_'+l]=round(c[l]/nl[l],3)
    row['isolates']=','.join(d['iso']); rows.append(row)
tsv('a20_lineage_allele_profiles',rows,list(rows[0]))
out['a20']={'cohort_lineages':dict(nl),'variants':len(rows),'restricted':sum(1 for r in rows if r['lineage_restricted']),'variants_only_in_MTB25':sum(1 for r in rows if r['only_in_MTB25_flagged']),'nonsyn_panel_variants':sum(1 for r in rows if 'missense' in r['effect'] or 'stop' in r['effect'] or 'frameshift' in r['effect'])}
# ---- 21 epitopes: variants in epitopes (cohort) + E11 conservation ranking + per-isolate partition table
rows=[]; epv=collections.defaultdict(lambda:collections.defaultdict(list))
for s in samples:
    for v in R[s]['variants']:
        if v['type']=='snp' and v.get('in_tcell_epitope') and not v.get('in_mask') and v.get('in_virulence_catalogue'):
            for e in str(v['epitope_ids']).split(';'):
                if e: epv[e][(v['rv'],v['gene'],v['hgvs_p'],v['effect'])].append(s)
tsv('a21_variants_in_epitopes',[{'epitope_id':e,'rv':k[0],'gene':k[1],'hgvs_p':k[2],'effect':k[3],'n_isolates':len(v),'isolates':','.join(v)} for e,d in epv.items() for k,v in d.items()],['epitope_id','rv','gene','hgvs_p','effect','n_isolates','isolates'])
rows=[]
for s in samples:
    ep=R[s]['epitope_partition']; 
    rows.append({'sample':s,'lineage':LIN[s],'epitope_variants':ep['epitope']['variants'],'epitope_pi_per_kb':ep['epitope']['pi_per_kb'],'epitope_pN_pS':ep['epitope']['pN_pS'],'nonepitope_variants':ep['non_epitope']['variants'],'nonepitope_pi_per_kb':ep['non_epitope']['pi_per_kb'],'nonepitope_pN_pS':ep['non_epitope']['pN_pS'],'ratio_pi_epitope_to_nonepitope':ep['ratio_epitope_to_non_epitope_pi']})
tsv('a21_per_isolate_epitope_partition',rows,list(rows[0]))
# ---- E11 epitope conservation ranking
ep=json.load(open(P+'/ref/epitopes_mapped.json'))
callf={s:{g['rv']:(g['callable_bp']/g['length_bp'] if g['length_bp'] else 0) for g in R[s]['catalogue_presence']} for s in samples}
vbygene=collections.defaultdict(list)
for s in samples:
    for v in R[s]['variants']:
        if v.get('rv') and v.get('in_virulence_catalogue') and not v.get('in_mask') and v.get('codon') and 'synonymous' not in v['effect'].replace('non_synonymous','x') and v['effect']!='synonymous_variant':
            vbygene[v['rv']].append((v['codon'],s,v['hgvs_p'] or v['effect']))
rows=[]
for e in ep:
    rv=e['rv']
    if rv not in callf[samples[0]]: continue
    ok=[s for s in samples if callf[s][rv]>=0.9]
    hits=[(c,s,h) for c,s,h in vbygene.get(rv,[]) if e['aa_start']<=c<=e['aa_end'] and s in ok]
    iso=sorted({s for c,s,h in hits}); sites=sorted({(c,h) for c,s,h in hits})
    L=e['aa_end']-e['aa_start']+1
    rows.append({'rv':rv,'gene':G.get(rv,{}).get('name',''),'peptide':e['pep'],'aa_start':e['aa_start'],'aa_end':e['aa_end'],'assessable_isolates':len(ok),
                 'isolates_with_change':len(iso),'frac_isolates_with_change':round(len(iso)/len(ok),3) if ok else '',
                 'distinct_changes':len(sites),'changes_per_residue':round(len(sites)/L,3),'conservation':('not assessable' if len(ok)<30 else 'invariant in all assessable isolates' if not iso else 'variable'),
                 'changes':';'.join('%s'%h for c,h in sites)})
rows.sort(key=lambda r:(r['conservation']!='variable',-(r['frac_isolates_with_change'] or 0) if r['conservation']=='variable' else 0))
tsv('e11_epitope_conservation_ranking',rows,list(rows[0]))
out['e11']=dict(collections.Counter(r['conservation'] for r in rows)); out['e11_most_variable']=[(r['rv'],r['gene'],r['peptide'],r['frac_isolates_with_change']) for r in rows if r['conservation']=='variable'][:8]
# ---- 27 mixed infection / within-host heterogeneity
rows=[]
for s in samples:
    r=R[s]; sc=r['subclonal']; pu=r['purity']; ex=r['exclusion']
    afs=[x['allele_fraction'] for x in sc]
    rows.append({'sample':s,'lineage':LIN[s],'foreign_read_percent':r['decontamination']['foreign_read_percent'],'mean_depth':r['coverage']['mean_depth_mqfiltered'],
                 'subclonal_sites_0.10-0.75_or_<0.9':len(sc),'subclonal_in_panel':sum(1 for x in sc if x.get('in_virulence_catalogue')),
                 'median_subclonal_AF':round(float(np.median(afs)),3) if afs else '', 'cohort_status':ex.get('cohort_status',''),
                 'mixed_infection_assessment':('not assessable (depth<10x or foreign DNA)' if (r['coverage']['mean_depth_mqfiltered'] or 0)<10 else 'possible mixed (>=30 subclonal sites)' if len(sc)>=30 else 'no evidence of mixed infection')})
tsv('a27_mixed_infection',rows,list(rows[0]))
out['a27']=dict(collections.Counter(r['mixed_infection_assessment'] for r in rows))
# ---- 29 deletions
rows=[]; dl=collections.defaultdict(list)
for s in samples:
    for d in R[s]['coverage']['deletions']:
        rows.append({'sample':s,'lineage':LIN[s],'start':d['start'],'end':d['end'],'length':d['length'],'genes':','.join(d['genes'])}); dl[(d['start']//500,d['genes'][0] if d['genes'] else '')].append(s)
tsv('a29_deletions',rows,['sample','lineage','start','end','length','genes'])
out['a29']={'deletion_rows':len(rows),'isolates_with_any':len({r['sample'] for r in rows}),'recurrent_regions':[(k[1],len(set(v))) for k,v in dl.items() if len(set(v))>=3][:20]}
# ---- 30 IS6110
rows=[]; ins=collections.defaultdict(set)
for s in samples:
    r=R[s]['is6110']; cn=R[s]['is6110_copy_number']
    cands=[j for j in r['junctions'] if not j['reference_copy']]
    for j in cands: ins[(j['start']//200)].add(s)
    rows.append({'sample':s,'lineage':LIN[s],'est_copy_number':cn['estimated_copy_number'],'junctions_supported':r['junctions_supported'],'reference_copies_supported':r['junctions_at_reference_copies'],'novel_insertions':len(cands),'novel_insertion_positions':';'.join('%d%s'%(j['start'],('('+j['gene']+')') if j['gene'] else '') for j in cands)})
tsv('a30_is6110',rows,list(rows[0]))
out['a30']={'recurrent_novel_sites':sum(1 for v in ins.values() if len(v)>=3),'median_copy_number':float(np.median([r['est_copy_number'] for r in rows]))}
json.dump(out,open(O+'/c2_tables_summary.json','w'),indent=1)
print(json.dumps(out,default=str)[:2500])
