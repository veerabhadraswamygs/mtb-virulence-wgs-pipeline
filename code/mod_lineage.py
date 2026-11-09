"""Lineage (Coll 62-SNP barcode) + spoligotype from mapped reads."""
import subprocess,collections,json,sys,gzip,os
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
T=P+'/src/tbdb/tbdb-master'
def run(cmd): return subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout

def lineage(sample,wd):
    bed=[]
    for l in open(T+'/barcode.bed'):
        f=l.rstrip('\n').split('\t')
        if len(f)<5: continue
        bed.append((int(f[2]),f[3],f[4],f[5] if len(f)>5 else '',f[6] if len(f)>6 else ''))
    regions=wd+'/barcode.pos'
    with open(regions,'w') as fh:
        for pos,*_ in bed: fh.write('NC_000962.3\t%d\n'%pos)
    out=run("%s/bin/bcftools mpileup -f %s/ref/H37Rv.fa -R %s -a AD -q 20 -Q 20 -d 1000 %s/%s.bam 2>/dev/null | "
            "%s/bin/bcftools call --ploidy 1 -m 2>/dev/null | grep -v '^#'"%(P,P,regions,wd,sample,P))
    gt={}
    for l in out.splitlines():
        f=l.split('\t')
        if len(f)<10: continue
        pos=int(f[1]); ref=f[3]; alts=f[4].split(',')
        d=dict(zip(f[8].split(':'),f[9].split(':')))
        if 'AD' not in d: continue
        ad=[int(x) for x in d['AD'].split(',') if x!='.']
        if sum(ad)<5: continue
        alleles=[ref]+[a for a in alts if a!='<*>']
        best=max(range(len(ad[:len(alleles)])),key=lambda i:ad[i]) if alleles else None
        if best is None: continue
        gt[pos]=(alleles[best],ad[best],sum(ad))
    hits=collections.Counter(); tot=collections.Counter(); called=0
    for pos,lin,allele,major,spol in bed:
        tot[lin]+=1
        if pos in gt:
            called+=1
            if gt[pos][0]==allele: hits[lin]+=1
    support={l:(hits[l],tot[l]) for l in tot if hits[l]>0}
    # the called lineage = most specific barcode with >=90% of its markers present
    good=[l for l,(h,t) in support.items() if t>0 and h/t>=0.90]
    good.sort(key=lambda l:(l.count('.'),len(l)))
    major=''
    for pos,lin,allele,maj,spol in bed:
        if good and lin==good[-1]: major=maj
    return {'barcode_positions':len(bed),'positions_called':called,
            'lineage':good[-1] if good else 'unassigned','major_lineage':major,
            'all_supported':{l:'%d/%d'%(h,t) for l,(h,t) in sorted(support.items())},
            'candidates':good}

DR_REGION='NC_000962.3:3115000-3132000'   # direct repeat locus in H37Rv
def spoligotype(sample,wd):
    """In-silico 43-spacer spoligotyping. Reads are taken from the DR region of the BAM
    plus unmapped reads, which is far faster than scanning the whole FASTQ and loses
    nothing: spacer-bearing reads either map to the DR locus or fail to map."""
    spacers=[]
    for l in open(T+'/spoligotype_spacers.txt'):
        f=l.rstrip('\n').split('\t')
        if len(f)>=2: spacers.append((f[1],f[0]))      # (name, sequence)
    comp=str.maketrans('ACGT','TGCA')
    pats=[(n,sq,sq.translate(comp)[::-1]) for n,sq in spacers]
    counts=collections.Counter()
    cmds=["%s/bin/samtools view %s/%s.bam %s"%(P,wd,sample,DR_REGION),
          "%s/bin/samtools view -f 4 %s/%s.bam"%(P,wd,sample)]
    seen=0
    for cmd in cmds:
        out=subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout
        for line in out.splitlines():
            f=line.split('\t',11)
            if len(f)<10: continue
            seq=f[9]; seen+=1
            for n,fwd,rev in pats:
                if fwd in seq or rev in seq: counts[n]+=1
    binary=''.join('1' if counts[n]>=5 else '0' for n,_ in spacers)
    oct_=''
    for i in range(0,len(binary)-1,3):
        chunk=binary[i:i+3]
        if len(chunk)==3: oct_+=str(int(chunk,2))
    oct_+=binary[-1]
    sit={}
    import csv as _csv
    with open(T+'/spoligotype_list.csv',encoding='utf8',errors='replace') as fh:
        for r in _csv.DictReader(fh):
            b=(r.get('Spoligo Binary') or '').strip()
            if not b: continue
            sit.setdefault(b.replace('n','1').replace('o','0'),(r.get('SIT','').strip(),r.get('Lineage (SITVIT2)','').strip()))
    hit=sit.get(binary)
    return {'binary':binary,'octal':oct_,'spacers_present':binary.count('1'),'reads_scanned':seen,
            'counts':dict(counts),'sit':hit[0] if hit else None,'family':hit[1] if hit else None,
            'n_patterns_in_db':len(sit)}

if __name__=='__main__':
    s=sys.argv[1]; wd=P+'/work/'+s
    L=lineage(s,wd); print(json.dumps(L,indent=1)[:900])
    S=spoligotype(s,wd); print('spoligo binary:',S['binary']); print('octal:',S['octal'],'| spacers:',S['spacers_present'],'| SIT:',S['sit'],S['family'])
