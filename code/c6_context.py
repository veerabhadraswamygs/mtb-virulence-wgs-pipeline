# Build ref/cohort_context.json: per-isolate cohort-level facts for the "Cohort context" report section (v5).
import csv,json,statistics
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
def rd(f): return list(csv.DictReader(open(f),delimiter='\t'))
C=P+'/cohort/'
dist={};cmp_={}
for fn,store in (('c1/a24_snp_distance_matrix.tsv',dist),('c1/a24_comparable_positions.tsv',cmp_)):
    L=[l.rstrip('\n').split('\t') for l in open(C+fn)]
    h=L[0][1:]
    for r in L[1:]:
        store[r[0]]={h[i]:int(v) for i,v in enumerate(r[1:]) if v!=''}
a24=json.load(open(C+'c1/a24_summary.json'))
c5={r['cluster_rep']:r for r in rd(C+'c5/c5_snp5_clusters.tsv')}
members={}
for r in c5.values():
    for m in r['members'].split(','): members[m.strip()]=r['members'].replace(',',', ')
mix={r['sample']:r for r in rd(C+'c2/a27_mixed_infection.tsv')}
is6={r['sample']:r for r in rd(C+'c2/a30_is6110.tsv')}
gl={r['isolate']:r for r in rd(C+'c4/a26_local_isolates_nearest_public_genome.tsv')}
l1={r['isolate']:r for r in rd(C+'c4/a26_L1_local_nearest_public_L1_genome.tsv')}
tb=json.load(open(C+'c4/a28_tbd1_summary.json'))
tbv=rd(C+'c4/a28_tbd1_variants_L1_local.tsv')
p328=[x for x in tbv if x['polymorphic_within_L1']=='True']
minor={m for x in p328 for m in x['minor_allele_isolates'].split(',') if m}
ps=[json.loads(l) if False else None for l in []]
iso=sorted(dist,key=lambda s:int(s[3:]))
cn=[float(r['est_copy_number']) for r in is6.values() if r['est_copy_number']]
med_cn=statistics.median(cn)
low_ass=set(a24['not_assessable_for_clustering_callable_lt_50pct'])
c12=[c for c in a24['clusters_12']]
out={'_meta':{'n_isolates':len(iso),'is6110_cohort_median':round(med_cn,1),'ingroup_variable_sites':a24['ingroup_variable_sites'],
  'n_clusters_5':len(a24['clusters_5']),'n_clusters_12':len(c12)}}
for s in iso:
    d={k:v for k,v in dist[s].items() if k!=s}
    nn=sorted(d.items(),key=lambda kv:(kv[1],kv[0]))[:3]
    g=gl.get(s,{}); l=l1.get(s,{})
    cl5=[c for c in a24['clusters_5'] if s in c]; cl12=[c for c in c12 if s in c]
    out[s]={'cluster5':cl5[0] if cl5 else None,'cluster12':cl12[0] if cl12 else None,
      'not_assessable_for_clustering':s in low_ass,
      'nearest_local':[{'isolate':k,'snps':v,'comparable_positions':cmp_[s].get(k)} for k,v in nn],
      'nearest_public_all_lineage_tree':{'accession':g.get('nearest_public'),'lineage':g.get('public_lineage'),'country':g.get('public_country'),'snv':g.get('approx_SNV_difference')} if g else None,
      'nearest_public_L1_context_tree':{'accession':l.get('nearest_public'),'lineage':l.get('public_lineage'),'country':l.get('public_country'),'snv':l.get('approx_SNV_difference')} if l and (l.get('local_lineage','').startswith('lineage1') or s=='MTB4') else None,
      'mixed_infection':mix[s]['mixed_infection_assessment'] if s in mix else None,
      'subclonal_sites':mix[s]['subclonal_sites_0.10-0.75_or_<0.9'] if s in mix else None,
      'cohort_status':mix[s]['cohort_status'] if s in mix else None,
      'is6110_copies':is6[s]['est_copy_number'] if s in is6 else None,
      'tbd1_L1_pos328_variant':s in minor}
json.dump(out,open(P+'/ref/cohort_context.json','w'),indent=1)
print(len(out)-1,out['MTB29'],out['MTB64']['tbd1_L1_pos328_variant'])
