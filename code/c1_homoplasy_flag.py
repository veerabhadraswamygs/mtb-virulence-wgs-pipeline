"""Flag homoplastic sites that sit within 10 bp of another homoplastic site (probable alignment/assembly artefacts) and summarise."""
import csv,json,collections,sys
O='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/cohort/c1'
for tre in sys.argv[1:]:
    rows=list(csv.DictReader(open('%s/a23_homoplasy_%s.tsv'%(O,tre),newline=''),delimiter='\t'))
    hp=sorted(int(r['pos']) for r in rows if r['homoplastic']=='True')
    for r in rows:
        p=int(r['pos']); r['clustered_with_other_homoplastic_site_within_10bp']=bool(r['homoplastic']=='True' and any(q!=p and abs(q-p)<=10 for q in hp))
    with open('%s/a23_homoplasy_%s.tsv'%(O,tre),'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
    h=[r for r in rows if r['homoplastic']=='True']; cl=[r for r in h if r['clustered_with_other_homoplastic_site_within_10bp']]
    pan=[r for r in h if r['panel_tier'] and not r['clustered_with_other_homoplastic_site_within_10bp']]
    s=json.load(open('%s/a23_summary_%s.json'%(O,tre))); s['homoplastic_clustered_within_10bp']=len(cl); s['homoplastic_in_panel_not_clustered']=len(pan)
    s['by_tier_not_clustered']=dict(collections.Counter(r['panel_tier'] for r in pan))
    s['top_panel_not_clustered']=[(r['pos'],r['rv'],r['gene'],r['panel_tier'],int(r['steps']),r['alleles']) for r in sorted(pan,key=lambda r:-int(r['steps']))[:15]]
    json.dump(s,open('%s/a23_summary_%s.json'%(O,tre),'w'),indent=1)
    print(tre,'homoplastic',len(h),'clustered',len(cl),'CI',s['consistency_index_overall'],'panel not clustered',len(pan),s['by_tier_not_clustered'],s['top_panel_not_clustered'][:4])
