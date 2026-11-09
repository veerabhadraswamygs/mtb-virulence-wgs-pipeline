"""Analysis 15 - PE/PPE/PE_PGRS secondary analysis (Tier 2, masked from the primary analysis).

Short-read mapping is unreliable across these repetitive, GC-rich loci, so the brief asks for
assembly rather than whole-genome mapping. Each locus is resolved from the isolate's own de novo
assembly: a locus is reported ONLY when a single contig aligns unambiguously over at least 80%
of its length at >=90% identity AND read depth over the locus is >=20x. Depth and mapping
quality are given for every call. Loci failing either test are reported as not resolved, never
as reference-identical.
"""
import subprocess,os,json,sys,collections,gzip
import numpy as np
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
MIN_DEPTH=20; MIN_ID=90.0; MIN_COV=0.80
def locus_depth_mq(sample,genes,rvs):
    d=np.zeros(4411533,dtype=np.int32)
    with gzip.open(P+'/work/%s/depth_all.txt.gz'%sample,'rt') as fh:
        for l in fh:
            f=l.split('\t'); d[int(f[1])]=int(f[2])
    out={}
    for rv in rvs:
        g=genes[rv]; seg=d[g['start']:g['end']+1]
        mq=subprocess.run("%s/bin/samtools view -q 0 %s/work/%s/%s.bam NC_000962.3:%d-%d | awk '{s+=$5; n++} END{if(n)printf \"%%.1f\",s/n; else print 0}'"
                          %(P,P,sample,sample,g['start'],g['end']),shell=True,capture_output=True,text=True).stdout.strip()
        out[rv]={'mean_depth':round(float(seg.mean()),1),'min_depth':int(seg.min()),
                 'mean_mapq':float(mq or 0)}
    return out
def run(sample):
    fa=P+'/work/%s/%s_denovo_p2.fasta'%(sample,sample)
    if not os.path.exists(fa): return {'available':False,'reason':'no de novo assembly'}
    genes=json.load(open(P+'/ref/genes.json')); cat=json.load(open(P+'/ref/catalogue.json'))
    rvs=sorted([rv for rv,e in cat.items() if e.get('tier')=='2'])
    out=subprocess.run("%s/bin/blastn -query %s -db %s/db/pe_ppe -outfmt "
        "'6 sseqid qseqid pident length slen qstart qend sstart send' -evalue 1e-20 -num_threads 2"%(P,fa,P),
        shell=True,capture_output=True,text=True).stdout
    byloc=collections.defaultdict(list)
    for l in out.splitlines():
        s_,q,pid,ln,slen,qs,qe,ss,se=l.split('\t')
        byloc[s_].append({'contig':q,'identity':float(pid),'length':int(ln),'slen':int(slen),
                          'q':(int(qs),int(qe)),'s':(int(ss),int(se))})
    dm=locus_depth_mq(sample,genes,rvs)
    rows=[]; res=collections.Counter()
    for rv in rvs:
        hits=[h for h in byloc.get(rv,[]) if h['identity']>=MIN_ID]
        g=genes[rv]; L=g['end']-g['start']+1
        cov=set()
        for h in hits:
            a,b=sorted(h['s']); cov.update(range(a,b+1))
        frac=len(cov)/L if L else 0
        contigs={h['contig'] for h in hits}
        d=dm[rv]
        if not hits or frac<MIN_COV:           call='not resolved (insufficient assembly coverage)'
        elif len(contigs)>1:                   call='not resolved (ambiguous: multiple contigs)'
        elif d['mean_depth']<MIN_DEPTH:        call='not resolved (read depth below 20x)'
        else:                                  call='resolved'
        res[call.split(' (')[0]]+=1
        best=max(hits,key=lambda h:h['length']) if hits else None
        rows.append({'rv':rv,'gene':(genes[rv]['name'] or rv),'length_bp':L,'call':call,
                     'assembly_coverage':round(frac,3),'contigs_aligned':len(contigs),
                     'best_identity':round(best['identity'],2) if best else None,
                     'mean_depth':d['mean_depth'],'min_depth':d['min_depth'],'mean_mapq':d['mean_mapq'],
                     'percent_divergence_from_H37Rv':round(100-best['identity'],2) if best else None})
    # pe_pgrs copy number from read mass over the family, as for IS6110
    pgrs=[rv for rv in rvs if 'PE_PGRS' in (genes[rv]['name'] or '')]
    tot=sum(dm[rv]['mean_depth']*(genes[rv]['end']-genes[rv]['start']+1) for rv in pgrs)
    allbp=sum(genes[rv]['end']-genes[rv]['start']+1 for rv in pgrs)
    gm=json.load(open(P+'/work/%s/results.json'%sample))['is6110_copy_number']['genome_mean_depth_unfiltered'] \
       if os.path.exists(P+'/work/%s/results.json'%sample) else None
    return {'available':True,'loci':len(rvs),'summary':dict(res),'rows':rows,
            'pe_pgrs_genes':len(pgrs),
            'pe_pgrs_mean_depth':round(tot/allbp,1) if allbp else None,
            'pe_pgrs_copy_equivalents':round((tot/allbp)/gm,2) if (allbp and gm) else None,
            'thresholds':{'min_identity':MIN_ID,'min_assembly_coverage':MIN_COV,'min_read_depth':MIN_DEPTH}}
if __name__=='__main__':
    r=run(sys.argv[1])
    print('loci %d | %s'%(r['loci'],r['summary']))
    print('PE_PGRS genes %d | mean depth %s | copy equivalents %s'%(r['pe_pgrs_genes'],r['pe_pgrs_mean_depth'],r['pe_pgrs_copy_equivalents']))
    for x in r['rows'][:5]: print('  ',x)
