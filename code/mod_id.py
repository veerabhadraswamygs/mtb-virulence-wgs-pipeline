"""Initial identification and validation, performed before the virulence analysis is interpreted.

  A. Species and complex-member identification
       ANI against H37Rv; identity of the standard mycobacterial identification targets
       (16S rRNA, rpoB, hsp65, gyrB); IS6110 as an MTBC marker; and the region-of-difference
       profile, which distinguishes members of the complex (Brosch et al. deletion scheme).
  B. Sequencing validation
       insert size, duplicate rate, coverage evenness, heterozygosity.
  C. Laboratory cross-check
       the identification and drug-susceptibility tests already performed on this isolate
       (smear, CBNAAT, first- and second-line LPA) set against the genome-based calls.
"""
import subprocess,json,os,re,sys,gzip
import numpy as np
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
R=P+'/ref'
MARKERS=[('rrs (16S rRNA)','MTB000019'),('rpoB','Rv0667'),('hsp65 (groEL2)','Rv0440'),('gyrB','Rv0005')]
RDS=[('RD1','Rv3871','Rv3879c','present in M. tuberculosis and M. bovis; deleted in BCG'),
     ('RD4','Rv1505c','Rv1516c','deleted in M. bovis'),
     ('RD7','Rv1964','Rv1977','deleted in M. africanum subtype II, M. microti, M. bovis'),
     ('RD8','Rv0309','Rv0314c','deleted in M. africanum subtype II, M. microti, M. bovis'),
     ('RD9','Rv2072c','Rv2075c','present in M. tuberculosis sensu stricto; deleted in M. africanum, M. bovis, M. microti'),
     ('RD10','Rv0221','Rv0223c','deleted in M. africanum subtype II, M. microti, M. bovis'),
     ('RD12','Rv3120','Rv3124','deleted in M. bovis and M. caprae')]

def _depth(sample):
    d=np.zeros(4411533,dtype=np.int32)
    with gzip.open(P+'/work/%s/depth_all.txt.gz'%sample,'rt') as fh:
        for l in fh:
            f=l.split('\t'); d[int(f[1])]=int(f[2])
    return d

def ani(sample):
    """ANIb-style: fragment the assembly into 1,020 nt pieces, align to H37Rv, average the
    identity of alignments passing 30% identity and 70% fragment coverage."""
    fa=P+'/work/%s/%s_denovo_p2.fasta'%(sample,sample)
    if not os.path.exists(fa): return None
    seq=[];name=None;frag=[]
    for l in open(fa):
        if l[0]=='>':
            if seq: frag.append(''.join(seq))
            seq=[]
        else: seq.append(l.strip())
    if seq: frag.append(''.join(seq))
    tmp=P+'/work/%s/_ani_frag.fa'%sample
    n=0
    with open(tmp,'w') as fh:
        for c in frag:
            for i in range(0,len(c)-1019,1020):
                n+=1; fh.write('>f%d\n%s\n'%(n,c[i:i+1020]))
    out=subprocess.run("%s/bin/blastn -query %s -db %s/db/h37rv -outfmt '6 qseqid pident length' "
                       "-evalue 1e-15 -max_target_seqs 1 -num_threads 2 -task blastn"%(P,tmp,P),
                       shell=True,capture_output=True,text=True).stdout
    best={}
    for l in out.splitlines():
        q,pid,ln=l.split('\t')
        if int(ln)>=714 and float(pid)>=30:
            if float(pid)>best.get(q,0): best[q]=float(pid)
    os.remove(tmp)
    if not best: return None
    return {'ANI_percent':round(sum(best.values())/len(best),3),'fragments_total':n,
            'fragments_aligned':len(best),'aligned_fraction':round(len(best)/n,3),
            'method':'ANIb-style: 1,020 nt assembly fragments aligned to H37Rv by BLASTn'}

def markers(sample):
    fa=P+'/work/%s/%s_denovo_p2.fasta'%(sample,sample)
    if not os.path.exists(fa): return []
    genes=json.load(open(R+'/genes.json'))
    s=''.join(l.strip() for l in open(R+'/H37Rv.fa') if l[0]!='>')
    q=P+'/work/%s/_markers.fa'%sample
    with open(q,'w') as fh:
        for nm,rv in MARKERS:
            g=genes.get(rv)
            if g: fh.write('>%s\n%s\n'%(rv,s[g['start']-1:g['end']]))   # Rv number only: BLAST truncates IDs at whitespace
    db=P+'/work/%s/_asm'%sample
    subprocess.run("%s/bin/makeblastdb -in %s -dbtype nucl -out %s >/dev/null 2>&1"%(P,fa,db),shell=True)
    out=subprocess.run("%s/bin/blastn -query %s -db %s -outfmt '6 qseqid pident length qlen' "
                       "-evalue 1e-20 -max_target_seqs 1 -num_threads 2"%(P,q,db),
                       shell=True,capture_output=True,text=True).stdout
    best={}
    for l in out.splitlines():
        qq,pid,ln,qlen=l.split('\t')
        cov=100*int(ln)/int(qlen)
        if (float(pid),cov)>best.get(qq,(0,0)): best[qq]=(float(pid),cov)
    for f in (q,): os.path.exists(f) and os.remove(f)
    for ext in ('.ndb','.nhr','.nin','.not','.nsq','.ntf','.nto'):
        p2=db+ext
        os.path.exists(p2) and os.remove(p2)
    return [{'target':nm,'locus':rv,
             'identity_percent':round(best[rv][0],2) if rv in best else None,
             'coverage_percent':round(best[rv][1],1) if rv in best else None} for nm,rv in MARKERS]

def rd_profile(sample,d=None):
    genes=json.load(open(R+'/genes.json'))
    if d is None: d=_depth(sample)
    rows=[]
    for rd,a,b,note in RDS:
        ga,gb=genes.get(a),genes.get(b)
        if not ga or not gb: continue
        s,e=min(ga['start'],gb['start']),max(ga['end'],gb['end'])
        seg=d[s:e+1]
        frac=float((seg>=5).mean())
        call='present' if frac>=0.90 else ('partial' if frac>=0.20 else 'deleted')
        rows.append({'region':rd,'span':'%d-%d'%(s,e),'length_bp':int(e-s+1),
                     'fraction_covered_5x':round(frac,3),'mean_depth':round(float(seg.mean()),1),
                     'call':call,'interpretation':note})
    by={r['region']:r['call'] for r in rows}
    if by.get('RD9')=='present':
        ident='Mycobacterium tuberculosis sensu stricto (RD9 intact)'
    elif by.get('RD4')=='deleted' and by.get('RD12')=='deleted':
        ident='Mycobacterium bovis (RD9, RD4 and RD12 deleted)'
    elif by.get('RD9')=='deleted':
        ident='Mycobacterium africanum / M. microti group (RD9 deleted, RD4 intact)'
    else:
        ident='MTBC member not resolved by this RD set'
    return {'regions':rows,'complex_member':ident,
            'note':'Region-of-difference scheme of Brosch et al.; coverage-based, computed on '
                   'mapping-quality-unfiltered depth so repeats are not mistaken for deletions.'}

def seq_validation(sample,d=None):
    if d is None: d=_depth(sample)
    g=d[1:]
    st=subprocess.run("%s/bin/samtools stats %s/work/%s/%s.bam"%(P,P,sample,sample),
                      shell=True,capture_output=True,text=True).stdout
    def sn(key):
        m=re.search(r'^SN\s+%s\s+([0-9.eE+-]+)'%re.escape(key),st,re.M)
        return float(m.group(1)) if m else None
    dup=sn('reads duplicated:'); tot=sn('sequences:')
    het=None
    rj=P+'/work/%s/results.json'%sample
    if os.path.exists(rj):
        R_=json.load(open(rj))
        nc=sum(b-a+1 for a,b in R_['coverage']['noncallable'])
        callable_bp=4411532-nc
        het=round(1e6*len(R_['subclonal'])/callable_bp,2) if callable_bp else None
    return {'insert_size_mean':sn('insert size average:'),'insert_size_sd':sn('insert size standard deviation:'),
            'duplicate_rate_percent':round(100*dup/tot,2) if (dup and tot) else None,
            'error_rate':sn('error rate:'),
            'mean_depth':round(float(g.mean()),1),'median_depth':float(np.median(g)),
            'depth_cv':round(float(g.std()/g.mean()),3) if g.mean() else None,
            'percent_within_2fold_of_median':round(100*float(((g>=np.median(g)/2)&(g<=np.median(g)*2)).mean()),2),
            'heterozygous_sites_per_Mb':het}

LAB=None
def lab_tests(sample):
    """Identification and susceptibility tests already performed on this isolate in the
    laboratory, read from the study metadata. Clinical fields only; no patient identifiers."""
    global LAB
    pth='/mnt/user-data/uploads/Downloads/Dr. Ranjitha TB-Pt final list.xlsx'
    if not os.path.exists(pth): return None
    if LAB is None:
        import warnings; warnings.filterwarnings('ignore')
        from openpyxl import load_workbook
        wb=load_workbook(pth,read_only=True,data_only=True)
        ws=wb['Sheet1']; it=ws.iter_rows(values_only=True); hdr=[str(x).strip() if x else '' for x in next(it)]
        LAB={}
        for r in it:
            if not r[0]: continue
            LAB[re.sub(r'\s+','',str(r[0])).upper()]=dict(zip(hdr,r))
        wb.close()
    row=LAB.get(sample.upper())
    if not row: return None
    keep=['SPUTUM SMEAR GRADING','CBNAAT RESULTS','FIRST LINE LPA RESULTS','SECOND LINE LPA RESULTS',
          'PREVIOUS  HISTORY OF TB','CLINICAL OUTCOME']
    return {k:(str(row.get(k)).strip() if row.get(k) is not None else '') for k in keep}

def run(sample):
    d=_depth(sample)
    return {'ani':ani(sample),'markers':markers(sample),'rd':rd_profile(sample,d),
            'validation':seq_validation(sample,d),'lab':lab_tests(sample)}
if __name__=='__main__':
    r=run(sys.argv[1])
    print('ANI:',r['ani'])
    print('markers:',r['markers'])
    print('complex member:',r['rd']['complex_member'])
    for x in r['rd']['regions']: print('   %-5s %-8s %s'%(x['region'],x['call'],x['fraction_covered_5x']))
    print('validation:',r['validation'])
    print('lab:',r['lab'])
