"""Map IEDB epitopes onto H37Rv genome coordinates.
Positions in the IEDB export mix precursor and mature-protein numbering and include
strain polymorphisms, so peptides are LOCATED BY SEQUENCE SEARCH in the H37Rv protein,
not by the reported coordinates. Exact match first, then a single-mismatch search.
"""
import gzip,csv,json,re,zipfile,io,collections,sys
U='/mnt/user-data/uploads/Downloads/limitation fix files'
R='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/ref'
genes=json.load(open(R+'/genes.json'))
up2rv={}
with gzip.open(U+'/uniprotkb_proteome_UP000001584_2026_09_29.tsv.gz','rt',errors='replace') as fh:
    for r in csv.DictReader(fh,delimiter='\t'):
        for rv in re.findall(r'Rv\d+[A-Za-z]?c?',r['Gene Names'] or ''):
            if rv in genes: up2rv.setdefault(r['Entry'],rv)
seq=[l.strip() for l in open(R+'/H37Rv.fa') if l[0]!='>']
G=''.join(seq)
B='TCAG'; AAs='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON={}; i=0
for a in B:
    for b in B:
        for c in B: CODON[a+b+c]=AAs[i]; i+=1
comp=str.maketrans('ACGT','TGCA')
PROT={}
def protein(rv):
    if rv in PROT: return PROT[rv]
    g=genes[rv]; s=G[g['start']-1:g['end']]
    if g['strand']=='-': s=s.translate(comp)[::-1]
    p=''.join(CODON.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
    PROT[rv]=p; return p
def find(pep,P):
    i=P.find(pep)
    if i>=0: return i+1,0
    L=len(pep); best=None
    for i in range(len(P)-L+1):
        d=0
        for a,b in zip(pep,P[i:i+L]):
            if a!=b:
                d+=1
                if d>1: break
        if d<=1:
            if best is None or d<best[1]: best=(i+1,d)
            if d==0: break
    return best if best else (None,None)
z=zipfile.ZipFile(U+'/epitope_table_export_1793116069.csv.zip')
rows=list(csv.reader(io.StringIO(z.read(z.namelist()[0]).decode('utf8','replace'))))
h1,h2,data=rows[0],rows[1],rows[2:]
def idx(g,n):
    for i,(a,b) in enumerate(zip(h1,h2)):
        if a==g and b==n: return i
iID=idx('Epitope ID','IEDB IRI') if idx('Epitope ID','IEDB IRI') is not None else 0
iName,iPar=idx('Epitope','Name'),idx('Epitope','Molecule Parent IRI')
st=collections.Counter(); out=[]
for r in data:
    pep=r[iName].strip().upper()
    if not pep or not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY]+',pep): st['not_plain_peptide']+=1; continue
    m=re.search(r'uniprot/([A-Z0-9]+)',r[iPar] or '')
    rv=up2rv.get(m.group(1)) if m else None
    if not rv: st['no_Rv_mapping']+=1; continue
    P=protein(rv)
    pos,mm=find(pep,P)
    if pos is None: st['peptide_not_found']+=1; continue
    st['exact' if mm==0 else 'one_mismatch']+=1
    g=genes[rv]; L=len(pep)
    aa_s,aa_e=pos,pos+L-1
    if g['strand']=='+':
        gs=g['start']+(aa_s-1)*3; ge=g['start']+aa_e*3-1
    else:
        ge=g['end']-(aa_s-1)*3; gs=g['end']-aa_e*3+1
    out.append({'iedb':r[iID].rsplit('/',1)[-1],'rv':rv,'gene':g['name'] or rv,'pep':pep,
                'aa_start':aa_s,'aa_end':aa_e,'g_start':gs,'g_end':ge,'mismatch':mm})
for k,v in st.most_common(): print('%-20s %5d'%(k,v))
print('mapped epitopes:',len(out),'| distinct genes:',len({o['rv'] for o in out}))
json.dump(out,open(R+'/epitopes.json','w'))
# verify genome coordinates round-trip for a sample
import random
random.seed(1); bad=0
for o in random.sample(out,300):
    g=genes[o['rv']]; s=G[o['g_start']-1:o['g_end']]
    if g['strand']=='-': s=s.translate(comp)[::-1]
    p=''.join(CODON.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
    d=sum(1 for a,b in zip(p,o['pep']) if a!=b)
    if d>o['mismatch']: bad+=1
print('genome round-trip check on 300 random epitopes: %d failures'%bad)
