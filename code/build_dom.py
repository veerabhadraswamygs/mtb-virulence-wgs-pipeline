"""UniProt feature reference: domains, regions, motifs, sites, TM, signal peptide (brief 13)."""
import gzip,csv,re,json
U='/mnt/user-data/uploads/Downloads/limitation fix files'
R='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/ref'
genes=json.load(open(R+'/genes.json'))
FT={'Domain [FT]':'DOMAIN','Region':'REGION','Motif':'MOTIF','Active site':'ACT_SITE',
    'Binding site':'BINDING','Transmembrane':'TRANSMEM','Signal peptide':'SIGNAL'}
pat=re.compile(r'([A-Z_]+)\s+(\d+)(?:\.\.(\d+))?')
note=re.compile(r'/note="([^"]*)"')
out={}
with gzip.open(U+'/uniprotkb_proteome_UP000001584_2026_10_27.tsv.gz','rt',errors='replace') as fh:
    for r in csv.DictReader(fh,delimiter='\t'):
        rvs=[x for x in re.findall(r'Rv\d+[A-Za-z]?c?',r['Gene Names'] or '') if x in genes]
        if not rvs: continue
        rv=rvs[0]
        feats=[]
        for col,kind in FT.items():
            v=(r.get(col) or '').strip()
            if not v: continue
            # split on the feature keyword so each record is parsed whole
            parts=re.split(r'(?=%s\s+\d)'%kind,v)
            for p in parts:
                m=pat.match(p.strip())
                if not m: continue
                s=int(m.group(2)); e=int(m.group(3)) if m.group(3) else s
                nm=note.search(p)
                feats.append({'type':kind,'aa_start':s,'aa_end':e,'note':nm.group(1) if nm else ''})
        if feats: out[rv]=feats
json.dump(out,open(R+'/domains.json','w'))
import collections
c=collections.Counter(f['type'] for v in out.values() for f in v)
print('genes with features:',len(out),'| features:',sum(len(v) for v in out.values()))
for k,v in c.most_common(): print('   %5d  %s'%(v,k))
for rv in ['Rv1908c','Rv0757','Rv3875']:
    print(rv,[ (f['type'],f['aa_start'],f['aa_end'],f['note'][:30]) for f in out.get(rv,[])][:4])
