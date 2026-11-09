"""Analysis 26 prep: align each public assembly to H37Rv (minimap2 asm5), call SNPs/indels
(paftools), record covered intervals, Coll 62-SNP barcode lineage (same tbdb barcode as the 55 isolates)."""
import subprocess,os,sys,gzip,json,collections
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
M=P+'/src/minimap2-2.31_x64-linux'; T=P+'/src/tbdb/tbdb-master'
OUT=P+'/pub/calls'; os.makedirs(OUT,exist_ok=True)
bed=[]
for l in open(T+'/barcode.bed'):
    f=l.rstrip('\n').split('\t')
    if len(f)>=5: bed.append((int(f[2]),f[3],f[4],f[5] if len(f)>5 else ''))
def process(acc):
    j=OUT+'/%s.json'%acc
    if os.path.exists(j) and 'intervals' in json.load(open(j)): return
    fa='/tmp/pub_%s.fna'%acc
    subprocess.run("zcat %s/pub/fna/%s.fna.gz > %s"%(P,acc,fa),shell=True,check=True)
    paf=subprocess.run("%s/minimap2 -cx asm5 --cs %s/ref/H37Rv.fa %s 2>/dev/null | sort -k6,6 -k8,8n"%(M,P,fa),shell=True,capture_output=True,text=True).stdout
    vcf=subprocess.run("%s/k8 %s/paftools.js call -f %s/ref/H37Rv.fa -l 1000 -L 1000 -"%(M,M,P),shell=True,input=paf,capture_output=True,text=True).stdout
    os.remove(fa)
    cov=[]
    for l in paf.splitlines():
        f=l.split('\t')
        if int(f[11])>=20 and int(f[10])>=1000: cov.append((int(f[7]),int(f[8])))
    cov.sort(); merged=[]
    for s,e in cov:
        if merged and s<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],e)
        else: merged.append([s,e])
    covered=sum(e-s for s,e in merged)
    def isc(p):
        for s,e in merged:
            if s<p<=e: return True
        return False
    snps={}
    with gzip.open(OUT+'/%s.vcf.gz'%acc,'wt') as o:
        for l in vcf.splitlines():
            o.write(l+'\n')
            if l.startswith('#'): continue
            f=l.split('\t')
            if len(f[3])==1 and len(f[4])==1 and f[4]!='.': snps[int(f[1])]=f[4]
    hits=collections.Counter(); tot=collections.Counter(); called=0
    for pos,lin,allele,major in bed:
        tot[lin]+=1
        if isc(pos):
            called+=1
            if snps.get(pos, None) is None:
                # covered and no variant => reference allele; need ref base
                if allele==REF[pos-1]: hits[lin]+=1
            elif snps[pos]==allele: hits[lin]+=1
    good=[l for l in tot if hits[l]>0 and hits[l]/tot[l]>=0.9]
    good.sort(key=lambda l:(l.count('.'),len(l)))
    major=''
    for pos,lin,allele,mj in bed:
        if good and lin==good[-1]: major=mj
    json.dump({'acc':acc,'covered_bp':covered,'covered_frac':round(covered/4411532,4),'n_snps':len(snps),
               'barcode_called':called,'lineage':good[-1] if good else 'unassigned','major_lineage':major,
               'candidates':good,'intervals':merged},open(j,'w'))
REF=''.join(l.strip() for l in open(P+'/ref/H37Rv.fa') if not l.startswith('>'))
if __name__=='__main__':
    for a in sys.argv[1:]: process(a)
