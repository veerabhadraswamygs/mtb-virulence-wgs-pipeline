"""Pipeline 2 - per-sample analysis. Reads the stage outputs in work/<sample>/ and
writes tables + a JSON summary. Every module records what it could not determine."""
import json,gzip,os,sys,re,collections,subprocess,csv
import numpy as np
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
R=P+'/ref'; GL=4411532
genes=json.load(open(R+'/genes.json'))
cat=json.load(open(R+'/catalogue.json'))
tss=json.load(open(R+'/tss_regions.json'))
dom=json.load(open(R+'/domains.json'))
epi=json.load(open(R+'/epitopes.json'))
chip=json.load(open(R+'/chipseq.json'))
MASK=[]
for l in open(R+'/mask_static.bed'):
    f=l.split('\t'); MASK.append((int(f[1])+1,int(f[2]),f[3].strip()))
MASK.sort()
mask_arr=np.zeros(GL+2,dtype=bool)
for s,e,_ in MASK: mask_arr[s:e+1]=True
epi_by_gene=collections.defaultdict(list)
for e in epi: epi_by_gene[e['rv']].append(e)
chip_by_pos=sorted((p['start'],p['end'],p['reg'],p['gene']) for p in chip)
chip_arr=np.zeros(GL+2,dtype=bool)
for s,e,_,_ in chip_by_pos:
    if 0<s<=GL: chip_arr[max(1,s):min(GL,e)+1]=True

def gene_at(pos):
    for rv,g in genes.items():
        if g['start'] and g['start']<=pos<=g['end']: return rv
    return None

def parse_ann(info):
    for blk in info.split(';'):
        if blk.startswith('ANN='):
            out=[]
            for a in blk[4:].split(','):
                f=a.split('|')
                if len(f)>10:
                    out.append({'effect':f[1],'impact':f[2],'gene':f[3],'gene_id':f[4],
                                'hgvs_c':f[9],'hgvs_p':f[10]})
            return out
    return []

AA3={'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H',
     'Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W',
     'Tyr':'Y','Val':'V','Ter':'*'}
def parse_hgvsp(h):
    m=re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)',h or '')
    if not m: return None,None,None
    return AA3.get(m.group(1),'X'),int(m.group(2)),AA3.get(m.group(3),m.group(3))

# BLOSUM62 for the impact substitute (brief 25)
B62_RAW="""A R N D C Q E G H I L K M F P S T W Y V
A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4"""
_l=B62_RAW.strip().splitlines(); _h=_l[0].split()
B62={}
for row in _l[1:]:
    p=row.split(); 
    for j,v in enumerate(p[1:]): B62[(p[0],_h[j])]=int(v)

def in_mask(pos): return bool(mask_arr[pos]) if 0<pos<=GL else False
def feats_at(rv,aa):
    out=[]
    for f in dom.get(rv,[]):
        if f['aa_start']<=aa<=f['aa_end']: out.append(f)
    return out
def epitopes_at(rv,aa):
    return [e for e in epi_by_gene.get(rv,[]) if e['aa_start']<=aa<=e['aa_end']]

# ---------------- variant table (brief 11,12,13,14,21a; extras E8,E9,E10) ----------------
def load_barcode():
    bc={}
    for l in open(P+'/src/tbdb/tbdb-master/barcode.bed'):
        f=l.rstrip('\n').split('\t')
        if len(f)>=5: bc[int(f[2])]=(f[3],f[4])
    return bc
BARCODE=load_barcode()

# position -> list of (rv, region kind); built once, circular-safe
_PROM=collections.defaultdict(list)
for _rv,_t in tss.items():
    for _a,_b in _t['promoter']:
        _PROM[('P',_rv)].append((_a,_b))
    for _a,_b in _t['utr_region']:
        _PROM[('U',_rv)].append((_a,_b))
_IDX=[]
for (_k,_rv),_ivs in _PROM.items():
    for _a,_b in _ivs: _IDX.append((_a,_b,_rv,_k))
_IDX.sort()
def promoter_context(pos):
    """Which gene's mapped promoter (P) or 5'UTR (U) does this position fall in?"""
    hits=[]
    for a,b,rv,k in _IDX:
        if a>pos: break
        if a<=pos<=b:
            t=tss[rv]; g=genes[rv]
            hits.append({'rv':rv,'gene':g['name'] or rv,'region':'promoter' if k=='P' else "5'UTR",
                         'class':t['class'],'utr':t['utr'],'tss':t['tss'],'distal':t['distal']})
    return hits

def chip_context(pos):
    return [{'regulator':r,'gene':g,'start':s,'end':e}
            for s,e,r,g in chip_by_pos if s<=pos<=e]

def variant_table(sample,wd):
    rows=[]
    for line in open(wd+'/%s.ann.vcf'%sample):
        if line[0]=='#': continue
        f=line.rstrip('\n').split('\t')
        pos=int(f[1]); ref=f[3]; alt=f[4]
        d=dict(zip(f[8].split(':'),f[9].split(':')))
        ad=[int(x) for x in d.get('AD','0').split(',') if x not in ('','.')]
        dp=sum(ad) if ad else int(d.get('DP',0) or 0)
        af=(max(ad[1:])/dp) if (len(ad)>1 and dp) else 0.0
        mq=re.search(r'MQ=(\d+)',f[7]); mq=int(mq.group(1)) if mq else None
        anns=parse_ann(f[7])
        a=anns[0] if anns else {}
        rv=a.get('gene_id') or a.get('gene') or ''
        if rv not in genes:
            g2=gene_at(pos); rv=g2 or ''
        wt,aa,mut=parse_hgvsp(a.get('hgvs_p',''))
        vtype=('indel' if len(ref)!=len(alt) else 'snp')
        row={'pos':pos,'ref':ref,'alt':alt,'type':vtype,'depth':dp,'allele_fraction':round(af,4),
             'mapping_quality':mq,'effect':a.get('effect',''),'impact':a.get('impact',''),
             'rv':rv,'gene':(genes[rv]['name'] if rv in genes and genes[rv]['name'] else rv),
             'hgvs_c':a.get('hgvs_c',''),'hgvs_p':a.get('hgvs_p',''),
             'codon':aa,'aa_ref':wt,'aa_alt':mut,
             'in_mask':in_mask(pos),'in_virulence_catalogue':rv in cat,
             'functional_system':cat[rv]['system'] if rv in cat else '',
             'panel_tier':cat[rv].get('tier','') if rv in cat else '',
             'fdr_stratum':cat[rv].get('fdr_stratum','') if rv in cat else '',
             'tested_panel':bool(cat[rv].get('tested')) if rv in cat else False,
             'lineage_defining':BARCODE.get(pos,('',''))[0],
             'in_chipseq_peak':bool(chip_arr[pos]) if 0<pos<=GL else False}
        # domain / site context
        if rv and aa:
            fts=feats_at(rv,aa)
            row['domains']=';'.join('%s:%s(%d-%d)'%(x['type'],x['note'] or '-',x['aa_start'],x['aa_end']) for x in fts)
            eps=epitopes_at(rv,aa)
            row['in_tcell_epitope']=bool(eps)
            row['epitope_ids']=','.join(e['iedb'] for e in eps[:6])
            row['blosum62']=B62.get((wt,mut)) if (wt and mut and wt!='*' and mut!='*') else None
        else:
            row['domains']=''; row['in_tcell_epitope']=False; row['epitope_ids']=''; row['blosum62']=None
        pc=promoter_context(pos)
        row['promoter_of']=';'.join('%s:%s(%s,5UTR=%d)'%(h['gene'],h['region'],h['class'],h['utr']) for h in pc)
        cc=chip_context(pos)
        row['chipseq_regulators']=','.join(sorted({c['regulator'] for c in cc}))
        rows.append(row)
    return rows

def subclonal(sample,wd,lo=0.10,hi=0.90):
    """Extra E14 + brief 27: minor alleles between lo and hi."""
    out=[]
    for line in gzip.open(wd+'/%s.raw.vcf.gz'%sample,'rt'):
        if line[0]=='#': continue
        f=line.rstrip('\n').split('\t')
        if f[4] in ('.','<*>'): continue
        d=dict(zip(f[8].split(':'),f[9].split(':')))
        ad=[int(x) for x in d.get('AD','0').split(',') if x not in ('','.')]
        if len(ad)<2: continue
        dp=sum(ad)
        if dp<20: continue
        mq=re.search(r'MQ=(\d+)',f[7]); mq=int(mq.group(1)) if mq else 0
        if mq<40: continue
        pos=int(f[1])
        if in_mask(pos): continue
        af=max(ad[1:])/dp
        if lo<=af<hi:
            rv=gene_at(pos) or ''
            out.append({'pos':pos,'ref':f[3],'alt':f[4].split(',')[0],'depth':dp,
                        'allele_fraction':round(af,4),'rv':rv,
                        'gene':(genes[rv]['name'] if rv in genes and genes[rv]['name'] else rv),
                        'in_virulence_catalogue':rv in cat})
    return out

# ---------------- per-isolate diversity in epitope vs non-epitope partitions (brief 21b) --------
def partition_pi(sample,rows,covjson):
    """Per-isolate pi and pN/pS computed separately for epitope and non-epitope partitions
    of the SAME genes (genes that carry at least one mapped epitope)."""
    gb=covjson['gene_breadth']
    epi_genes=set(epi_by_gene)
    callable_bp={'epitope':0,'non_epitope':0}
    for rv in epi_genes:
        g=genes.get(rv); b=gb.get(rv)
        if not g or not b: continue
        frac=b['breadth10']
        elen=0
        seen=set()
        for e in epi_by_gene[rv]:
            for a in range(e['aa_start'],e['aa_end']+1): seen.add(a)
        elen=len(seen)*3
        glen=g['end']-g['start']+1
        callable_bp['epitope']+=elen*frac
        callable_bp['non_epitope']+=max(0,glen-elen)*frac
    cnt={'epitope':{'n':0,'syn':0,'nonsyn':0},'non_epitope':{'n':0,'syn':0,'nonsyn':0}}
    for r in rows:
        if r['rv'] not in epi_genes or r['in_mask']: continue
        part='epitope' if r['in_tcell_epitope'] else 'non_epitope'
        cnt[part]['n']+=1
        if r['effect']=='synonymous_variant': cnt[part]['syn']+=1
        elif r['effect'] in ('missense_variant','stop_gained','stop_lost'): cnt[part]['nonsyn']+=1
    out={}
    for part in ('epitope','non_epitope'):
        bp=callable_bp[part]
        out[part]={'variants':cnt[part]['n'],'callable_bp':round(bp,1),
                   'pi_per_kb':round(1000*cnt[part]['n']/bp,4) if bp>0 else None,
                   'syn':cnt[part]['syn'],'nonsyn':cnt[part]['nonsyn'],
                   'pN_pS':round(cnt[part]['nonsyn']/cnt[part]['syn'],3) if cnt[part]['syn'] else None}
    out['genes_considered']=len([rv for rv in epi_genes if rv in gb])
    out['ratio_epitope_to_non_epitope_pi']=(round(out['epitope']['pi_per_kb']/out['non_epitope']['pi_per_kb'],3)
        if out['epitope']['pi_per_kb'] and out['non_epitope']['pi_per_kb'] else None)
    return out

# ---------------- IS6110 insertion sites (brief 30) ----------------
IS_SITES=json.load(open(R+'/is6110_sites.json'))['copies_in_H37Rv']
def is6110(sample,wd,minsup=5,win=50):
    """Cluster soft-clipped reads, then KEEP ONLY clusters whose clipped sequence
    matches the IS6110 element by BLAST. Clusters within `win` of a reference IS6110
    copy are reported as reference copies, the rest as candidate novel insertions."""
    out=subprocess.run("%s/bin/samtools view -F 0x904 %s/%s.bam"%(P,wd,sample),shell=True,
                       capture_output=True,text=True).stdout
    clip=collections.defaultdict(list)
    for line in out.splitlines():
        f=line.split('\t',11)
        if len(f)<11: continue
        cig,pos,seq=f[5],int(f[3]),f[9]
        m=re.match(r'^(\d+)S',cig)
        if m and int(m.group(1))>=20: clip[pos].append(seq[:int(m.group(1))])
        n=re.search(r'(\d+)S$',cig)
        if n and int(n.group(1))>=20:
            end=pos+sum(int(x) for x,o in re.findall(r'(\d+)([MDN])',cig))
            clip[end].append(seq[-int(n.group(1)):])
    # merge nearby positions
    pts=sorted(clip)
    groups=[]
    for pt in pts:
        if groups and pt-groups[-1][1]<=20: groups[-1][1]=pt; groups[-1][2]+=clip[pt]
        else: groups.append([pt,pt,list(clip[pt])])
    groups=[g for g in groups if len(g[2])>=minsup]
    # blast the clipped sequences
    fa=wd+'/clips.fa'
    with open(fa,'w') as fh:
        for gi,(a,b,seqs) in enumerate(groups):
            for si,sq in enumerate(seqs[:30]): fh.write('>g%d_%d\n%s\n'%(gi,si,sq))
    hits=collections.Counter()
    if groups:
        bl=subprocess.run("%s/bin/blastn -query %s -db %s/db/is6110 -outfmt '6 qseqid pident length' "
                          "-evalue 1e-5 -max_target_seqs 1 -num_threads 2"%(P,fa,P),
                          shell=True,capture_output=True,text=True).stdout
        for l in bl.splitlines():
            q,pid,ln=l.split('\t')[:3]
            if float(pid)>=95 and int(ln)>=20: hits[int(q.split('_')[0][1:])]+=1
    res=[]
    for gi,(a,b,seqs) in enumerate(groups):
        n=min(len(seqs),30)
        if hits.get(gi,0) < max(2,0.5*n): continue
        ref=any(abs(a-x)<=win or abs(b-x)<=win for x in IS_SITES)
        rv=gene_at(a) or ''
        res.append({'start':a,'end':b,'clipped_reads':len(seqs),'is6110_matching':hits[gi],
                    'reference_copy':ref,'rv':rv,
                    'gene':(genes[rv]['name'] if rv in genes and genes[rv]['name'] else rv)})
    novel=[r for r in res if not r['reference_copy']]
    return {'H37Rv_reference_copies':len(IS_SITES),'junctions_supported':len(res),
            'junctions_at_reference_copies':len(res)-len(novel),
            'candidate_insertions_absent_from_H37Rv':len(novel),
            'junctions':sorted(res,key=lambda x:-x['clipped_reads'])}

def is6110_copy_number(wd):
    """Copy number from read mass over the 16 reference IS6110 elements relative to
    genome mean depth, using MQ-UNFILTERED depth (repeat reads have MQ 0)."""
    d=np.zeros(GL+1,dtype=np.int32)
    with gzip.open(wd+'/depth_all.txt.gz','rt') as fh:
        for l in fh:
            f=l.split('\t'); d[int(f[1])]=int(f[2])
    gmean=float(d[1:].mean())
    per=[float(d[s:s+1355].mean()) for s in IS_SITES]
    total=sum(per)
    return {'genome_mean_depth_unfiltered':round(gmean,2),
            'mean_depth_per_reference_element':round(float(np.mean(per)),2),
            'estimated_copy_number':round(total/gmean,1) if gmean else None,
            'method':'sum of mean depth over the 16 H37Rv IS6110 elements divided by genome mean depth'}

def _unused_is6110_reference_retention(covjson):
    """Which of the 16 H37Rv IS6110 copies are retained in this isolate, judged by
    read depth across the element (absent copies show as zero/low coverage)."""
    nc=[(a,b) for a,b in covjson['noncallable']]
    out=[]
    for site in IS_SITES:
        a,b=site,site+1354
        low=sum(max(0,min(b,y)-max(a,x)+1) for x,y in nc if not(y<a or x>b))
        frac=low/1355.0
        out.append({'h37rv_site':site,'fraction_below_10x':round(frac,3),
                    'call':'absent' if frac>0.8 else ('partial' if frac>0.2 else 'retained')})
    ret=sum(1 for x in out if x['call']=='retained')
    return {'retained':ret,'partial':sum(1 for x in out if x['call']=='partial'),
            'absent':sum(1 for x in out if x['call']=='absent'),'detail':out}

# ---------------- assembly stats (brief 1), completeness (2), purity (3) ----------------
def assembly_path(sample):
    return P+'/work/%s/%s_denovo_p2.fasta'%(sample,sample)
def proteome_path(sample):
    return P+'/work/%s/%s_proteome_p2.faa'%(sample,sample)
def assembly_stats(sample):
    fa=assembly_path(sample)
    if not os.path.exists(fa): return {'available':False,'reason':'de novo assembly not produced'}
    lens=[];cur=0;gc=0;n=0;tot=0
    for line in open(fa):
        if line[0]=='>':
            if cur: lens.append(cur)
            cur=0
        else:
            t=line.strip().upper(); cur+=len(t); gc+=t.count('G')+t.count('C'); n+=t.count('N')
    if cur: lens.append(cur)
    lens.sort(reverse=True); tot=sum(lens); acc=0; n50=0; l50=0
    for i,L in enumerate(lens,1):
        acc+=L
        if acc>=tot/2: n50=L; l50=i; break
    return {'available':True,'contigs':len(lens),'total_length':tot,'largest_contig':lens[0] if lens else 0,
            'N50':n50,'L50':l50,'GC_percent':round(100*gc/tot,2) if tot else None,'ambiguous_bases':n,
            'source':'SPAdes 4.3.0 --isolate on this isolate\u2019s quality-trimmed reads (analysis 5b); contigs '
                     'filtered to \u2265500 bp, k-mer coverage \u22655, and a BLASTn hit to the MTBC reference panel',
            'filter':json.load(open(P+'/work/%s/asm_filter.json'%sample)) if os.path.exists(P+'/work/%s/asm_filter.json'%sample) else None}

def completeness(covjson):
    """Substitute for CheckM2: breadth of coverage across single-copy reference CDS."""
    gb=covjson['gene_breadth']
    sc=[rv for rv,g in genes.items() if g['type']=='CDS' and g['category'] not in ('PE/PPE','insertion seqs and phages')]
    vals=[gb[rv]['breadth10'] for rv in sc if rv in gb]
    present=sum(1 for v in vals if v>=0.95)
    return {'single_copy_cds_assessed':len(vals),'cds_with_>=95%_at_10x':present,
            'completeness_percent':round(100*present/len(vals),2) if vals else None,
            'method':'breadth of coverage >=10x across single-copy H37Rv CDS (CheckM2 unavailable offline)'}

def purity(sample,wd,fastp_json):
    j=json.load(open(fastp_json))
    tot=j['summary']['after_filtering']['total_reads']
    fl=subprocess.run("%s/bin/samtools flagstat %s/%s.bam"%(P,wd,sample),shell=True,
                      capture_output=True,text=True).stdout
    mapped=int(re.search(r'(\d+) \+ \d+ mapped',fl).group(1))
    prim=int(re.search(r'(\d+) \+ \d+ primary',fl).group(1))
    pm=int(re.search(r'(\d+) \+ \d+ primary mapped',fl).group(1)) if 'primary mapped' in fl else mapped
    rate=100*pm/prim if prim else 0
    return {'reads_after_qc':tot,'primary_reads':prim,'primary_mapped':pm,
            'percent_mapped_to_H37Rv':round(rate,2),'percent_unmapped':round(100-rate,2),
            'method':'mapping rate against H37Rv (Kraken2 standard DB does not fit the available disk)'}

# ---------------- decontamination (applied identically to every isolate) ----------------
def decontamination(covjson,pur):
    """Reads that do not belong to M. tuberculosis are excluded from every analysis:
    variant calling uses only reads aligned to H37Rv, and the assembly keeps only contigs
    with a BLASTn hit to the MTBC reference panel. This records how much was removed and
    what the isolate's own data looks like once it has been.
    Note that coverage is measured from aligned reads only, so the depth figures were never
    inflated by the foreign fraction - decontamination removes reads, it cannot add any."""
    nc=sum(b-a+1 for a,b in covjson['noncallable'])
    breadth=100*(GL-nc)/GL
    return {'foreign_read_percent':pur['percent_unmapped'],
            'retained_read_percent':pur['percent_mapped_to_H37Rv'],
            'depth_after_decontamination':round(covjson['mean_depth_all'],1),
            'breadth_10x_after_decontamination':round(breadth,2),
            'method':'read level: only reads aligned to H37Rv enter variant calling and coverage; '
                     'contig level: only contigs with a BLASTn hit to the MTBC reference panel are kept',
            'sufficient_after_decontamination':bool(breadth>=90)}

# ---------------- exclusion rule (brief section 2) ----------------
def exclusion(covjson,pur):
    d=covjson
    br=None
    # genome-wide breadth at >=10x from the non-callable intervals
    nc=sum(b-a+1 for a,b in d['noncallable'])
    br=100*(GL-nc)/GL
    reasons=[]
    if br<90: reasons.append('less than 90%% of the genome covered at >=10x (%.2f%%)'%br)
    if pur['percent_unmapped']>5: reasons.append('more than 5%% of reads not mapping to H37Rv (%.2f%%)'%pur['percent_unmapped'])
    # the brief's rule is applied to the RAW library. Where the only failure is the foreign-read
    # share and the isolate's own data is sufficient once decontaminated, that is recorded
    # separately so the decision is visible rather than silent.
    only_foreign = bool(reasons) and all('not mapping' in r for r in reasons)
    return {'genome_breadth_10x_percent':round(br,2),'percent_unmapped':pur['percent_unmapped'],
            'excluded':bool(reasons),'reasons':reasons,
            'fails_only_on_foreign_reads':only_foreign,
            'usable_after_decontamination':bool(br>=90),
            # study decision: every isolate enters the cohort analyses; weak ones are flagged, never dropped
            'cohort_eligible':True,
            'cohort_status':('full' if not reasons else
                             'retained_after_decontamination' if (only_foreign and br>=90) else 'low_depth'),
            'rule':'report (never silently drop) any isolate with <90% of H37Rv covered at >=10x or >5% non-M. tuberculosis reads'}

# ---------------- virulence gene presence / absence (brief 16) ----------------
def vfdb_presence(sample,minid=90.0,mincov=60.0):
    """BLASTn of the isolate assembly against VFDB (M. tuberculosis entries).
    Genes whose H37Rv locus is not callable are reported 'not determined', never 'absent'."""
    fa=assembly_path(sample)
    if not os.path.exists(fa): return None
    out=subprocess.run("%s/bin/blastn -query %s -db %s/db/vfdb_setB_mtb -outfmt "
        "'6 sseqid pident length slen qseqid' -evalue 1e-20 -num_threads 2 -max_target_seqs 5"%(P,fa,P),
        shell=True,capture_output=True,text=True).stdout
    best=collections.defaultdict(float)
    for l in out.splitlines():
        s,pid,ln,slen,q=l.split('\t')
        if float(pid)<minid: continue
        cov=100*int(ln)/int(slen)
        if cov>best[s]: best[s]=cov
    return best

def catalogue_presence(sample,covjson,minbreadth=0.95):
    """Per-isolate presence call for every catalogue gene, using coverage of the H37Rv locus.
    Low-coverage loci are 'not determined' (brief: never call them absent)."""
    gb=covjson['gene_breadth']
    rows=[]
    for rv,e in cat.items():
        b=gb.get(rv)
        if not b: rows.append({'rv':rv,'gene':e['name'],'system':e['system'],'tier':e.get('tier',''),
                               'call':'not determined','breadth10':None,'mean_depth':None,
                               'length_bp':e.get('length_bp'),'callable_bp':0}); continue
        if b['breadth10']>=minbreadth: call='present'
        elif b['breadth1']<0.10:       call='absent'
        else:                          call='not determined'
        rows.append({'rv':rv,'gene':e['name'],'system':e['system'],'tier':e.get('tier',''),
                     'call':call,'breadth10':round(b['breadth10'],3),'mean_depth':round(b['mean'],1),
                     'length_bp':e.get('length_bp'),
                     'callable_bp':int(round((b['breadth10'] or 0)*(e.get('length_bp') or 0)))})
    return rows

# ---------------- AMRFinder concordance (extra E3) ----------------
def amrfinder(sample):
    fa=proteome_path(sample)
    if not os.path.exists(fa):
        src=assembly_path(sample)
        if not os.path.exists(src): return {'available':False,'reason':'no assembly available'}
        subprocess.run('%s/bin/prodigal -i %s -a %s -p single -q >/dev/null 2>&1'%(P,src,fa),shell=True)
    if not os.path.exists(fa): return {'available':False,'reason':'proteome prediction failed'}
    out=subprocess.run("%s/bin/blastp -query %s -db %s/db/amrprot -outfmt "
        "'6 qseqid sseqid pident length slen evalue' -evalue 1e-30 -num_threads 2 -max_target_seqs 1"%(P,fa,P),
        shell=True,capture_output=True,text=True).stdout
    hits=[]
    for l in out.splitlines():
        q,s,pid,ln,slen,ev=l.split('\t')
        if float(pid)>=90 and 100*int(ln)/int(slen)>=80: hits.append({'query':q,'amr_ref':s,'identity':float(pid)})
    return {'available':True,'hits':len(hits),'detail':hits[:200]}


def tier_callable(pres):
    """Per-isolate callable length per FDR stratum. This is the denominator the brief
    specifies for analysis 31 (variant burden normalised by callable sequence length)."""
    out={}
    for t in ('1a','1b','2','3'):
        sub=[p for p in pres if p['tier']==t]
        nominal=sum(p['length_bp'] or 0 for p in sub)
        callable_=sum(p['callable_bp'] or 0 for p in sub)
        out[t]={'genes':len(sub),'nominal_bp':nominal,'callable_bp':callable_,
                'callable_fraction':round(callable_/nominal,4) if nominal else None}
    return out
