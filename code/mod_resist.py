"""Resistance calling (brief 7) with WHO v2 confidence grading (E1) and
variants of uncertain significance in resistance genes (E2)."""
import csv,collections,re,json
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
T=P+'/src/tbdb/tbdb-master'
GRADE_ORDER=['Assoc w R','Assoc w R - Interim','Uncertain significance',
             'Not assoc w R - Interim','Not assoc w R']
def load_catalogue():
    cat=collections.defaultdict(list)
    with open(T+'/mutations.csv',encoding='utf8',errors='replace') as fh:
        for r in csv.DictReader(fh):
            cat[(r['Gene'].strip(),r['Mutation'].strip())].append(r)
    genes={g for g,_ in cat}
    return cat,genes
def parse_ann(info):
    out=[]
    for blk in info.split(';'):
        if not blk.startswith('ANN='): continue
        for a in blk[4:].split(','):
            f=a.split('|')
            if len(f)>10: out.append({'effect':f[1],'impact':f[2],'gene':f[3],'hgvs_c':f[9],'hgvs_p':f[10]})
    return out
def call(vcf_annotated):
    cat,rgenes=load_catalogue()
    hits=[]; vus=[]
    for line in open(vcf_annotated):
        if line[0]=='#': continue
        f=line.rstrip('\n').split('\t')
        pos=int(f[1])
        for a in parse_ann(f[7]):
            if a['gene'] not in rgenes: continue
            for mut in [a['hgvs_p'],a['hgvs_c']]:
                if not mut: continue
                recs=cat.get((a['gene'],mut))
                if not recs: continue
                for r in recs:
                    row={'pos':pos,'ref':f[3],'alt':f[4],'gene':a['gene'],'mutation':mut,
                         'effect':a['effect'],'drug':r['drug'],'confidence':r['confidence'] or 'unspecified',
                         'source':r['source'],'type':r['type']}
                    # an 'Uncertain significance' grade is never a resistance call, whichever source lists it
                    if row['confidence']=='Uncertain significance':
                        vus.append(row)
                    elif row['confidence'] in ('Assoc w R','Assoc w R - Interim') or r['type']=='drug_resistance':
                        hits.append(row)
                break
    def dedup(rows):
        seen=set(); out=[]
        for x in rows:
            k=(x['gene'],x['mutation'],x['drug'])
            if k not in seen: seen.add(k); out.append(x)
        return out
    return dedup(hits),dedup(vus),sorted(rgenes)
if __name__=='__main__':
    import sys
    h,v,g=call(P+'/work/%s/%s.ann.vcf'%(sys.argv[1],sys.argv[1]))
    print('resistance genes in catalogue:',len(g))
    print('RESISTANCE-ASSOCIATED hits:',len(h))
    for r in h: print('   %-8s %-10s %-14s %-22s %s'%(r['gene'],r['mutation'],r['drug'],r['confidence'],r['pos']))
    print('UNCERTAIN-SIGNIFICANCE variants in resistance genes:',len(v))
    for r in v[:8]: print('   %-8s %-12s %-14s %s'%(r['gene'],r['mutation'],r['drug'],r['confidence']))
