"""C6 analysis 25 (SIFT4G, CPU build) on the same Tier 1a/1b missense variants as PROVEAN; mycobacterial UniRef90 as database."""
import csv,json,os,collections
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
W='/mnt/user-data/outputs/_sift'
h=json.load(open(P+'/ref/h37rv_prot.json'))
r=list(csv.DictReader(open(P+'/cohort/c4/a25_variant_structure_descriptors.tsv'),delimiter='\t'))
by=collections.defaultdict(set)
for x in r:
    s=h.get(x['rv']);
    if not s: continue
    c=int(x['codon'])
    if 1<c<=len(s) and s[c-1]==x['wt'] and x['mt'] not in ('*',''): by[x['rv']].add('%s%s%s'%(x['wt'],x['codon'],x['mt']))
with open('/tmp/claude-0/sift/q.fa','w') as q:
    for rv in sorted(by):
        q.write('>%s\n%s\n'%(rv,h[rv])); open('%s/subst/%s.subst'%(W,rv),'w').write('\n'.join(sorted(by[rv]))+'\n')
print(len(by),'genes',sum(len(v) for v in by.values()),'variants')
# ---- parse
import glob,re
res={}
for f in glob.glob(W+'/out/*.SIFTprediction'):
    rv=os.path.basename(f).split('.')[0]
    for l in open(f):
        t=l.rstrip('\n').split('\t')
        if len(t)>=6 and re.match(r'^[A-Z]\d+[A-Z]$',t[0]):
            res[(rv,t[0])]=(t[1],t[2],t[3],t[4],t[5])
rows=[]
for x in r:
    v=res.get((x['rv'],'%s%s%s'%(x['wt'],x['codon'],x['mt'])))
    rows.append({'rv':x['rv'],'gene':x['gene'],'tier':x['tier'],'hgvs_p':x['hgvs_p'],'n_isolates':x['n_isolates'],
     'SIFT4G_score':v[1] if v else '','SIFT4G_call_(<=0.05_deleterious)':v[0].lower() if v else 'not scored',
     'median_sequence_info':v[2] if v else '','sequences_at_position':v[3] if v else '','total_sequences':v[4] if v else '',
     'low_support_(<10_sequences_at_position)':('yes' if v and int(v[3])<10 else 'no') if v else ''})
w=csv.DictWriter(open(P+'/cohort/c4/a25_sift4g_scores.tsv','w'),fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
print('scored',len(res),'of',len(r))
