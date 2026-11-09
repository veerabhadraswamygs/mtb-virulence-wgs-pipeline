"""Analysis 28 - TbD1 locus (mmpS6 / mmpL6).

TbD1 is deleted from H37Rv, which is a lineage 4 (modern) strain, so the region has no
coordinates in the reference and NO mapping-based test against H37Rv can detect it.
The locus is therefore recovered from the de novo assembly.

Reference sequence: the 2,139 bp segment present in M. canettii CIPT 140010059 and absent
from H37Rv, located by aligning the H37Rv region spanning Rv1556-Rv1559 (the deletion
junction lies at mmpL6, Rv1557) against M. canettii. Its two open reading frames match
MmpS-family (Rv0403c 58%, Rv0451c 55%) and MmpL-family (Rv0450c 60%, Rv0676c 59%)
proteins, confirming mmpS6 and mmpL6. Presence across a 13-genome MTBC panel follows the
expected ancestral-present / modern-absent pattern.

Rv1557 in the catalogue is annotated mmpL6 but is a DIFFERENT gene; it is not used here.
"""
import subprocess,os,json,sys,re
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
TBD1_LEN=2139
# An 'absent' call is only meaningful if the assembly actually contains the region around the deletion
# junction; otherwise absence of TbD1 sequence cannot be told apart from absence of that part of the
# assembly. FLANK is the H37Rv region (NC_000962.3:1,757,997-1,767,717) that spans the junction.
FLANK=P+'/ref/h37rv_tbd1_flank.fa'
FLANK_LEN=9721
MIN_FLANK=0.80
def run(sample):
    fa=P+'/work/%s/%s_denovo_p2.fasta'%(sample,sample)
    if not os.path.exists(fa):
        return {'available':False,'reason':'no de novo assembly'}
    out=subprocess.run("%s/bin/blastn -query %s -db %s/db/tbd1 -outfmt "
        "'6 qseqid qstart qend sstart send pident length evalue' -evalue 1e-20 -num_threads 2"%(P,fa,P),
        shell=True,capture_output=True,text=True).stdout
    hits=[]
    for l in out.splitlines():
        q,qs,qe,ss,se,pid,ln,ev=l.split('\t')
        if float(pid)>=90 and int(ln)>=200:
            hits.append({'contig':q,'q_start':int(qs),'q_end':int(qe),'ref_start':int(ss),
                         'ref_end':int(se),'identity':float(pid),'length':int(ln)})
    covered=set()
    for h in hits:
        a,b=sorted((h['ref_start'],h['ref_end']))
        covered.update(range(a,b+1))
    frac=len(covered)/TBD1_LEN
    call=('present' if frac>=0.80 else ('partial' if frac>=0.20 else 'absent'))
    # fraction of the junction-spanning H37Rv region found on a SINGLE assembly contig
    flank_frac=0.0
    if os.path.getsize(fa)>0:
        fo=subprocess.run("%s/bin/blastn -query %s -subject %s -outfmt '6 sseqid qstart qend pident length' "
            "-evalue 1e-20 -num_threads 2"%(P,FLANK,fa),shell=True,capture_output=True,text=True).stdout
        per={}
        for l in fo.splitlines():
            sid,qs,qe,pid,ln=l.split('\t')
            if float(pid)>=90 and int(ln)>=200:
                a,b=sorted((int(qs),int(qe))); per.setdefault(sid,set()).update(range(a,b+1))
        flank_frac=max((len(v) for v in per.values()),default=0)/FLANK_LEN
    if call=='absent' and flank_frac<MIN_FLANK:
        call='not_determined'
    res={'available':True,'call':call,'fraction_of_TbD1_recovered':round(frac,3),
         'aligned_bp':len(covered),'reference_length':TBD1_LEN,'hits':hits,
         'reference':'M. canettii CIPT 140010059, 2,139 bp segment absent from H37Rv',
         'flank_fraction_covered':round(flank_frac,3),'flank_length':FLANK_LEN,
         'flank_reference':'H37Rv NC_000962.3:1,757,997-1,767,717, spans the TbD1 deletion junction',
         'method':'BLASTn of the de novo assembly against the TbD1 reference; no H37Rv mapping is used'}
    if call=='not_determined':
        res['reason']=('the assembly does not contain the region around the TbD1 deletion junction '
                       '(%.1f%% of it on a single contig, below the %d%% needed), so absence of TbD1 sequence '
                       'cannot be distinguished from absence of that part of the assembly'%(100*flank_frac,int(100*MIN_FLANK)))
    if call=='present' and hits:
        best=max(hits,key=lambda h:h['length'])
        seq={}
        name=None
        for line in open(fa):
            if line[0]=='>': name=line[1:].split()[0]; seq[name]=[]
            else: seq[name].append(line.strip())
        s=''.join(seq.get(best['contig'],[]))
        a,b=sorted((best['q_start'],best['q_end']))
        pad=200
        locus=s[max(0,a-1-pad):min(len(s),b+pad)]
        res['locus_sequence_length']=len(locus)
        res['locus_contig']=best['contig']
        open(P+'/work/%s/%s_TbD1.fasta'%(sample,sample),'w').write(
            '>%s_TbD1 contig=%s %d-%d identity=%.2f\n%s\n'%(sample,best['contig'],a,b,best['identity'],locus))
    return res
if __name__=='__main__':
    r=run(sys.argv[1])
    print(json.dumps({k:v for k,v in r.items() if k!='hits'},indent=1))
    for h in r.get('hits',[])[:5]: print('  ',h)
