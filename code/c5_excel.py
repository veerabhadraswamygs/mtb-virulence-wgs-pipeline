"""C5 deliverable 2: Excel workbook (one sheet per analysis) + TSV copies of every table."""
import json,os,re,glob,collections,pickle,numpy as np,pandas as pd
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; C=P+'/cohort'
F='/mnt/user-data/outputs/FINAL_DELIVERABLES'; TS=F+'/2_tables/tsv'; os.makedirs(TS,exist_ok=True)
O='/mnt/user-data/outputs'
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
R={s:json.load(open('%s/work/%s/results.json'%(P,s))) for s in samples}
SET=set(samples)
D=pickle.load(open(C+'/c5/c5_data.pkl','rb')); LIN=D['LIN']
def flat(d,pre=''):
    o={}
    for k,v in d.items():
        kk=(pre+'_'+k) if pre else k
        if isinstance(v,dict): o.update(flat(v,kk))
        elif isinstance(v,list): o[kk]=v if (len(v)>0 and not isinstance(v[0],(dict,list))) and len(v)<12 else ('[%d items]'%len(v) if v else '')
        else: o[kk]=v
    return o
def listify(df):
    for c in df.columns:
        df[c]=df[c].map(lambda x:';'.join(map(str,x)) if isinstance(x,list) else x)
    return df
T={}   # name -> df
def per_iso(key,sub=None,drop=()):
    rows=[]
    for s in samples:
        d=R[s][key] if sub is None else R[s][key][sub]
        f=flat(d); [f.pop(k,None) for k in list(f) if any(k.startswith(x) for x in drop)]
        rows.append(dict(isolate=s,lineage=LIN[s],**f))
    return listify(pd.DataFrame(rows))
T['A01_assembly_stats']=per_iso('assembly')
T['A02_completeness']=per_iso('completeness')
ident=[]; rdrows=[]
for s in samples:
    I=R[s]['identification']; r=dict(isolate=s,**{('ANI_'+k if not k.startswith('ANI') else k):v for k,v in (I.get('ani') or {}).items()})
    for m in (I.get('markers') or []): r['marker_%s_identity'%m['target'].split()[0]]=m['identity_percent']
    r.update({'decontam_'+k:v for k,v in R[s]['decontamination'].items() if k!='method'}); r.update({'purity_'+k:v for k,v in R[s]['purity'].items() if k!='method'})
    ident.append(r)
    for x in ((I.get('rd') or {}).get('regions') or []): rdrows.append(dict(isolate=s,**x))
T['A03_species_purity']=pd.DataFrame(ident); T['A03b_RD_regions']=pd.DataFrame(rdrows)
T['A04_read_QC']=per_iso('fastp')
T['A04b_read_filtering']=per_iso('fastp_filtering')
rows=[]
for s in samples:
    for p in R[s]['catalogue_presence']: rows.append(dict(isolate=s,**p))
T['A05_gene_callable_breadth']=pd.DataFrame(rows)
rows=[]
for s in samples:
    for d in R[s]['coverage']['deletions']: rows.append(dict(isolate=s,start=d['start'],end=d['end'],length=d['length'],genes=';'.join(d['genes'])))
T['A05b_coverage_gaps_deletions']=pd.DataFrame(rows)
T['A05c_denovo_assembly_filter']=per_iso('assembly',drop=('source',))
rows=[]
for s in samples:
    L=R[s]['lineage']; sp=R[s]['spoligotype']
    rows.append(dict(isolate=s,lineage_call=L['lineage'],major_lineage=L['major_lineage'],barcode_positions_called='%s/%s'%(L['positions_called'],L['barcode_positions']),support=json.dumps(L['all_supported']),
                     spoligotype_binary=sp['binary'],spoligotype_octal=sp['octal'],spacers_present=sp['spacers_present'],SIT=sp['sit'],family=sp['family']))
T['A06_lineage_spoligotype']=pd.DataFrame(rows)
rows=[]
for s in samples:
    for r in R[s]['resistance']+R[s]['resistance_uncertain']: rows.append(dict(isolate=s,**r))
T['A07_resistance']=pd.DataFrame(rows)
rows=[]
for s in samples:
    for h in R[s]['amrfinder'].get('detail',[]): rows.append(dict(isolate=s,**h))
T['A07b_AMRFinder_crosscheck']=pd.DataFrame(rows)
T['A08_reannotation']=pd.DataFrame(json.load(open(C+'/a8/a8_summary.json')))
a9=json.load(open(C+'/a9/a9_pangenome.json'))
T['A09_pangenome_summary']=pd.DataFrame([dict(item=k,value=json.dumps(v) if not isinstance(v,(int,float,str)) else v) for k,v in a9.items() if k!='rows'])
T['A09_pangenome_clusters']=pd.DataFrame(a9['rows'])
rows=[];eff=[];dom=[];ld=[]
for s in samples:
    V=R[s]['variants']; snp=[v for v in V if v['type']=='snp']
    rows.append(dict(isolate=s,lineage=LIN[s],mean_depth=R[s]['coverage']['mean_depth_mqfiltered'],breadth_10x_pct=R[s]['exclusion']['genome_breadth_10x_percent'],cohort_status=R[s]['exclusion']['cohort_status'],
                     snps=len(snp),indels=len(V)-len(snp),snps_outside_mask=sum(1 for v in snp if not v['in_mask']),subclonal_sites=len(R[s]['subclonal'])))
    c=collections.Counter(v['effect'] for v in V if not v['in_mask']); eff.append(dict(isolate=s,**c))
    for v in V:
        if v.get('in_mask'): continue
        if v['effect'] in('missense_variant','stop_gained','start_lost','stop_lost','frameshift_variant') and v.get('rv'):
            dom.append(dict(isolate=s,rv=v['rv'],gene=v['gene'],hgvs_p=v['hgvs_p'],codon=v['codon'],effect=v['effect'],domains=v['domains'],in_tcell_epitope=v['in_tcell_epitope'],epitope_ids=v['epitope_ids'],blosum62=v['blosum62'],panel_tier=v['panel_tier']))
        if v.get('lineage_defining'): ld.append(dict(isolate=s,pos=v['pos'],ref=v['ref'],alt=v['alt'],gene=v['gene'],lineage_defining=v['lineage_defining']))
T['A10_A11_variant_calling_summary']=pd.DataFrame(rows)
T['A12_functional_annotation_counts']=pd.DataFrame(eff).fillna(0)
T['A13_nonsyn_codon_protein_domain']=pd.DataFrame(dom)
T['A14_lineage_defining_variants']=pd.DataFrame(ld)
allv=[]
for s in samples:
    for v in R[s]['variants']: allv.append(dict(isolate=s,**v))
T['A11_variants_all']=pd.DataFrame(allv)
rows=[]
for s in samples:
    if R[s]['pe_ppe'].get('rows'):
        for r in R[s]['pe_ppe']['rows']: rows.append(dict(isolate=s,**r))
T['A15_pe_ppe_reassembly']=pd.DataFrame(rows)
# existing cohort tables
def ld_(path,name):
    T[name]=pd.read_csv(path,sep='\t',low_memory=False).rename(columns={'Unnamed: 0':'isolate'})
C2=O+'/cohort_C1_C2/C2'; C1=O+'/cohort_C1_C2/C1'; C3=O+'/cohort_C3'; C4=O+'/cohort_C4'; C5=C+'/c5'
ld_(C2+'/a16_presence_matrix.tsv','A16_virulence_presence_matrix'); ld_(C2+'/a17_pan_virulome.tsv','A17_pan_virulome')
for k in['pooled','L1','L3','L4']: ld_(C2+'/a18_pi_tajimasD_%s.tsv'%k,'A18_pi_TajimaD_%s'%k)
ld_(C2+'/a20_lineage_allele_profiles.tsv','A20_lineage_allele_profiles')
ld_(C3+'/a19_codeml_M1a_vs_M2a_panel.tsv','A19_codeml_M1a_vs_M2a'); ld_(C3+'/a19_hyphy_FEL_MEME_sites_p0.1.tsv','A19b_HyPhy_FEL_MEME')
ld_(C2+'/a21_per_isolate_epitope_partition.tsv','A21b_epitope_partition_isolate'); ld_(C2+'/a21_variants_in_epitopes.tsv','A21a_variants_in_epitopes')
ld_(C1+'/a24_snp_distance_matrix.tsv','A24_snp_distance_matrix'); ld_(C1+'/a24_comparable_positions.tsv','A24b_comparable_positions')
ld_(C1+'/a23_homoplasy_tree55.tsv','A23_homoplasy_55_tree'); ld_(C1+'/a23_homoplasy_treeX.tsv','A23b_homoplasy_L1context_tree')
ld_(C4+'/a25_variant_structure_descriptors.tsv','A25_structure_descriptors'); ld_(C4+'/a25b_foldx_ddG_candidates.tsv','A25b_FoldX_ddG'); ld_(C4+'/a25b_foldx_ddG_large_proteins.tsv','A25b_FoldX_large_proteins'); ld_(C4+'/a25_provean_scores.tsv','A25c_PROVEAN_scores'); ld_(C4+'/a25_sift4g_scores.tsv','A25c_SIFT4G_scores'); ld_(C4+'/a25d_foldx_alphafold_complexes.tsv','A25d_FoldX_AF_complexes')
ld_(C4+'/a26_variant_public_frequency.tsv','A26_variant_public_frequency'); ld_(C4+'/a26_local_isolates_nearest_public_genome.tsv','A26b_nearest_public_all_lineage'); ld_(C4+'/a26_L1_local_nearest_public_L1_genome.tsv','A26c_nearest_public_L1'); ld_(C4+'/a26_public_accession_list.tsv','A26d_public_accessions')
ld_(C2+'/a27_mixed_infection.tsv','A27_mixed_infection'); ld_(C4+'/a28_tbd1_variants_L1_local_plus_public.tsv','A28_TbD1_variants'); ld_(C2+'/a29_deletions.tsv','A29_deletions_CNV'); ld_(C2+'/a30_is6110.tsv','A30_IS6110')
ld_(C5+'/a31_a33_burden_diabetes.tsv','A31_A33_burden_diabetes'); ld_(C5+'/a33_burden_within_lineage_diabetes.tsv','A33_burden_within_lineage')
ld_(C5+'/a32_a33_variant_gene_tests_diabetes.tsv','A32_A33_variant_gene_tests'); ld_(C5+'/a34_epitope_partition_diabetes.tsv','A34_epitope_partition_perm')
T['A33_lineage_by_diabetes']=pd.read_csv(C5+'/a33_lineage_by_diabetes_crosstab.tsv',sep='\t')
lt=json.load(open(C5+'/a33_lineage_vs_diabetes_test.json')); T['A33_lineage_by_diabetes']['permutation_p_lineage_vs_diabetes']=lt['permutation_p']
ld_(C5+'/a35_metadata_audit_variables.tsv','A35a_metadata_audit'); ld_(C5+'/a35_burden_other_clinical.tsv','A35b_burden_other_clinical'); ld_(C5+'/a35_burden_within_lineage_other_clinical.tsv','A35c_within_lineage_other')
ld_(C5+'/a35_variant_gene_tests_other_clinical.tsv','A35d_variant_gene_tests_other'); ld_(C5+'/a35_age_burden_spearman.tsv','A35e_age_spearman')
for suf,nm in[('cluster_collapsed','S_cluster_collapsed'),('strict_sensitivity','S_strict')]:
    ld_(C5+'/a31_a33_burden_diabetes_%s.tsv'%suf,'Sens_burden_%s'%nm); ld_(C5+'/a32_a33_variant_gene_tests_diabetes_%s.tsv'%suf,'Sens_variant_tests_%s'%nm)
ld_(C5+'/a32_sensitivity_excl_lowcov_diabetes.tsv','Sens_variant_tests_excl_lowcov')
for suf,nm in[('cluster_collapsed','cluster_collapsed'),('strict_sensitivity','strict')]:
    ld_(C5+'/a35_burden_other_clinical_%s.tsv'%suf,'Sens_A35_burden_%s'%nm); ld_(C5+'/a35_variant_gene_tests_other_clinical_%s.tsv'%suf,'Sens_A35_tests_%s'%nm)
ld_(C5+'/a35_burden_other_clinical_excl_lowcov.tsv','Sens_A35_burden_excl_lowcov'); ld_(C5+'/c5_snp5_clusters.tsv','Sens_snp5_clusters')
a36=json.load(open(C5+'/a36_multivariable_decision.json'))
T['A36_multivariable']=pd.DataFrame([dict(item=k,value=json.dumps(v,default=float) if isinstance(v,(list,dict)) else v) for k,v in a36.items()])
# trees summary
c1=json.load(open(C1+'/c1_core_summary.json')); T['A22_tree_summary']=pd.DataFrame([dict(item=k,value=json.dumps(v) if not isinstance(v,(int,float,str)) else v) for k,v in c1.items()])
ld_(C3+'/e15_genomewide_dNdS_per_gene.tsv','E15_genomewide_dNdS'); ld_(C3+'/e12_epitope_density_vs_dNdS_per_gene.tsv','E12_epitope_density_dNdS')
ld_(C2+'/e4_cog_category_burden.tsv','E4_cog_category_burden'); ld_(C2+'/e4_panel_vs_background.tsv','E4_panel_vs_background'); ld_(C2+'/e11_epitope_conservation_ranking.tsv','E11_epitope_conservation')
ld_(C4+'/e5_global_allele_frequency_per_lineage.tsv','E5_global_allele_freq'); ld_(C4+'/e6_local_vs_global_gene_presence_flagged.tsv','E6_local_vs_global_presence')
# ID map
import openpyxl
ws=openpyxl.load_workbook('/mnt/user-data/uploads/Downloads/Dr. Ranjitha TB-Pt final list.xlsx',data_only=True)['Sheet1']
meta={re.sub(r'\s+','',str(r[0])).upper():str(r[0]).strip() for r in ws.iter_rows(min_row=2,values_only=True) if r[0]}
T['ID_map']=pd.DataFrame([dict(pipeline_id=s,metadata_sheet_id=meta[s.upper()],lineage=LIN[s],diabetes=D['clin'][s]['diabetes'],low_quality_flag='MTB25 contaminated (77.6% foreign reads)' if s=='MTB25' else ('low coverage' if s in D['lowcov'] else '')) for s in samples])
# ---------------- write TSV (pipeline IDs) and Excel (metadata-sheet spelling)
pat=re.compile(r'\bMTB(\d+)\b')
def xid(x):
    if isinstance(x,str): return pat.sub(lambda m:('MTB %s'%m.group(1)) if 'MTB'+m.group(1) in SET else m.group(0),x)
    return x
def to_xl(df):
    d=df.copy(); d.columns=[xid(str(c)) for c in d.columns]
    for c in d.columns:
        if not pd.api.types.is_numeric_dtype(d[c]) and not pd.api.types.is_bool_dtype(d[c]): d[c]=d[c].map(xid)
    return d
index=[]
with pd.ExcelWriter(F+'/2_tables/MTB_virulence_all_analyses_tables.xlsx',engine='xlsxwriter',engine_kwargs={'options':{'constant_memory':False}}) as xw:
    readme=pd.DataFrame([dict(sheet='README',note='One sheet per analysis (A##/E##). Sheet names start with the brief analysis number. Isolate IDs in this workbook are written exactly as in the metadata sheet (for example "MTB 3"); the TSV copies and the trees/VCFs use the file-safe form "MTB3". Sheet ID_map links the two.'),
                         dict(sheet='',note='Statistics 31-36: all comparisons performed are listed, including null results. Tests restricted to Tier 1; BH within Tier 1a and Tier 1b separately (and within analysis family); variants and genes tested only if present in >=3 isolates.'),
                         dict(sheet='',note='Large sheets are split only by analysis; no sheet was filtered to significant results.')])
    readme.to_excel(xw,sheet_name='README',index=False); xw.sheets['README'].set_column(0,0,18); xw.sheets['README'].set_column(1,1,150)
    for name in sorted(T,key=lambda n:(n[0]!='A' and n[0]!='E', n)):
        df=T[name]; df.to_csv('%s/%s.tsv'%(TS,name),sep='\t',index=False)
        x=to_xl(df); nm=name[:31]
        x.to_excel(xw,sheet_name=nm,index=False,freeze_panes=(1,0))
        ws_=xw.sheets[nm]; ws_.set_column(0,min(len(x.columns)-1,40),14); ws_.autofilter(0,0,len(x)-0 if len(x) else 0,max(len(x.columns)-1,0))
        index.append(dict(sheet=nm,rows=len(df),columns=len(df.columns)))
pd.DataFrame(index).to_csv(F+'/2_tables/sheet_index.tsv',sep='\t',index=False)
print(len(index),'sheets'); print(pd.DataFrame(index).to_string())
