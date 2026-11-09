"""Verification suite. Every number printed in a report is re-derived independently from
the source data and compared. A report is not released unless every check passes."""
import json,os,sys,re,subprocess,collections,csv
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import analyse as A
P=A.P
class Checker:
    def __init__(s,sample): s.sample=sample; s.ok=0; s.fail=[]
    def eq(s,name,a,b,tol=0):
        good=(abs(a-b)<=tol) if isinstance(a,(int,float)) and isinstance(b,(int,float)) else (a==b)
        if good: s.ok+=1
        else: s.fail.append('%s: got %r expected %r'%(name,a,b))
    def true(s,name,cond,detail=''):
        if cond: s.ok+=1
        else: s.fail.append('%s%s'%(name,(' — '+detail) if detail else ''))


def run_panel_checks(c,R,txt,outdir,sample):
    """Checks specific to the supplied frozen panel."""
    import hashlib,os
    EXPECTED={'panel_tier1a_core_virulence.bed':'a73ead619009c734f7d2a25cc64c04ed',
              'panel_tier1b_toxin_antitoxin.bed':'dff8d6ae0d20f7b4c4dddf5061608f03',
              'panel_tier2_pe_ppe.bed':'75917c35aa6f89b09070547b1e7bf09f',
              'panel_tier3_exploratory.bed':'c30cdbf19adfe295f17f627e6214b4e6'}
    for f,h in EXPECTED.items():
        p=os.path.join(P,'ref','panel',f)
        c.true('panel BED present: %s'%f,os.path.exists(p))
        if os.path.exists(p):
            c.eq('panel BED md5: %s'%f,hashlib.md5(open(p,'rb').read()).hexdigest(),h)
    import collections
    t=collections.Counter(e['tier'] for e in A.cat.values())
    for k,n in [('1a',211),('1b',120),('2',168),('3',447)]: c.eq('tier %s gene count'%k,t[k],n)
    for k,bp in [('1a',338667),('1b',38262),('2',281073),('3',483184)]:
        c.eq('tier %s nominal bp'%k,sum(e['length_bp'] for e in A.cat.values() if e['tier']==k),bp)
    c.eq('total panel bp',sum(e['length_bp'] for e in A.cat.values()),1141186)
    c.eq('tested genes (1a+1b)',sum(1 for e in A.cat.values() if e['tested']),331)
    c.true('tier 2 and 3 are not tested',all(not e['tested'] for e in A.cat.values() if e['tier'] in ('2','3')))
    c.true('length_bp equals end-start+1',all(e['length_bp']==e['end']-e['start']+1 for e in A.cat.values()))
    tc=R.get('tier_callable',{})
    for k in ('1a','1b','2','3'):
        c.true('tier %s callable <= nominal'%k,tc[k]['callable_bp']<=tc[k]['nominal_bp'])
        c.eq('tier %s nominal bp in results'%k,tc[k]['nominal_bp'],
             {'1a':338667,'1b':38262,'2':281073,'3':483184}[k])
    for s_ in ['338,667','38,262','281,073','483,184','1,141,186','946']:
        c.true('panel figure in text: %s'%s_,s_ in txt)
    c.true('tiering explained in text','Benjamini' in txt and 'separate' in txt.lower())
    c.true('Tier 3 stated as never tested','never formally tested' in txt)
    c.true('md5 verification stated','md5' in txt)
    pr=R['catalogue_presence']
    c.eq('presence rows carry a tier',sum(1 for p in pr if p['tier'] in ('1a','1b','2','3')),len(A.cat))
    return c


def run(sample,outdir):
    c=Checker(sample)
    W=P+'/work/'+sample
    R=json.load(open(W+'/results.json'))
    cov=json.load(open(W+'/cov.json'))
    pdf=os.path.join(outdir,'%s_complete_analysis_report_2.pdf'%sample)
    htm=os.path.join(outdir,'%s_report2.html'%sample)
    raw=subprocess.run(['pdftotext','-layout',pdf,'-'],capture_output=True,text=True).stdout
    # pdftotext wraps lines, so phrase searches run against a whitespace-normalised copy
    txt=re.sub(r'\s+',' ',raw)
    H=open(htm).read()
    rows=R['variants']

    # ---- structural integrity of the document
    for sec in ['1. Organism identification','2. Read quality control','3. Species purity',
                '4. Coverage, callable regions','5. De novo assembly','6. Lineage, sublineage',
                '7. Drug resistance','8. Variant calling','9. Variants in the frozen virulence panel',
                '10. T-cell epitope variation','11. Promoter, 5','12. Sub-clonal variants',
                '13. Structural variation','14. TbD1 locus','15. PE / PPE','16. Virulence gene presence',
                '17. Cohort context','18. Limitations','19. Methods']:
        c.true('section present: %s'%sec, sec.replace('5′','5').replace('5,','5') in txt.replace('′',''))
    c.true('PDF has pages', len(txt)>5000, 'extracted %d chars'%len(txt))
    # ---- input data integrity and cross-report comparison statements (pipeline v4)
    _fj=json.load(open(W+'/%s_fastp.json'%sample))
    _n1=_fj['read1_before_filtering']['total_reads']; _n2=_fj['read2_before_filtering']['total_reads']
    c.true('R1 and R2 read counts equal in the analysed data',_n1==_n2,'%d vs %d'%(_n1,_n2))
    c.true('read pairing row shown with both counts','{:,} / {:,}'.format(_n1,_n2) in txt)
    c.true('input data integrity statement present','Input data integrity.' in txt)
    if os.path.exists(W+'/input_integrity.json'):
        _ii=json.load(open(W+'/input_integrity.json'))
        c.true('truncated input is stated plainly, with pairs analysed',
               'was truncated' in txt and '{:,} read pairs complete in both files'.format(_ii['pairs_analysed']) in txt)
    else:
        c.true('complete input stated, and no truncation claimed','Both read files were complete' in txt and 'was truncated' not in txt)
    c.true('methods explain that count differences between reports are not genetic distance',
           'Comparing reports between isolates.' in txt and 'not a measure of genetic distance' in txt)
    for bad in ['None','nan','%s','{}','TODO','PLACEHOLDER','FIXME','N/A%']:
        c.true('no placeholder %r in document'%bad, ('>'+bad+'<') not in H)
    c.true('no unformatted python dict in text', '{\'' not in txt)
    # the brief REQUIRES tool and database release dates, so those are expected; what must not
    # appear is a run/analysis date, which would differ between otherwise identical reports
    c.true('no analysis/run date printed', 'Analysis date' not in txt and 'Report date' not in txt
           and 'Generated on' not in txt)
    import datetime as _dt
    _today=_dt.date.today()
    for _fmt in ('%d %B %Y','%Y-%m-%d','%d/%m/%Y','%B %d, %Y'):
        c.true('today\'s date not stamped (%s)'%_fmt,_today.strftime(_fmt) not in txt)
    _dates=set(re.findall(r'\b20[0-9]{2}-[01][0-9]-[0-3][0-9]\b',txt))
    _allowed=set(re.findall(r'\b20[0-9]{2}-[01][0-9]-[0-3][0-9]\b',json.dumps(json.load(open(P+'/ref/versions.json')))))
    c.true('any ISO date present is a declared tool/database release date',
           _dates<=_allowed,'unexpected: %s'%sorted(_dates-_allowed))
    # figures
    figs=re.findall(r'<img src="([^"]+)"',H)
    _no_effects=not any(r['effect'] for r in rows)
    _no_systems=not any(r['in_virulence_catalogue'] and not r['in_mask'] and r['functional_system'] for r in rows)
    _undrawn={n for n,miss in ((2,_no_effects),(4,_no_systems)) if miss}
    c.eq('figure count',len(figs),6-len(_undrawn))
    for f in figs: c.true('figure file exists: %s'%os.path.basename(f), os.path.exists(f))
    c.eq('figure blocks',H.count('class="figblock"'),len(figs))
    for n in range(1,7):
        if n in _undrawn:
            # no data to plot: the absence must be stated in the report, not left silent
            c.true('figure %d not drawn because there is no data, and the report says so'%n,
                   'Figure not drawn' in H)
        else:
            c.true('figure %d captioned'%n,'<b>Figure %d.</b>'%n in H)
    c.eq('undrawn figures each announced',H.count('Figure not drawn'),len(_undrawn))

    # ---- numbers: variants
    snp=sum(1 for r in rows if r['type']=='snp'); ind=sum(1 for r in rows if r['type']=='indel')
    c.eq('snp+indel == total',snp+ind,len(rows))
    unm=[r for r in rows if not r['in_mask']]; msk=[r for r in rows if r['in_mask']]
    c.eq('masked+unmasked == total',len(unm)+len(msk),len(rows))
    c.true('variant total appears in text','{:,}'.format(len(rows)) in txt or str(len(rows)) in txt)
    c.true('unmasked count appears in text','{:,}'.format(len(unm)) in txt or str(len(unm)) in txt)
    # recount variants straight from the annotated VCF
    nvcf=sum(1 for l in open(W+'/%s.ann.vcf'%sample) if l[0]!='#')
    c.eq('variant table length == annotated VCF lines',len(rows),nvcf)
    # every variant must satisfy the brief's filters
    c.true('all variants depth>=10',all(r['depth']>=10 for r in rows))
    c.true('all variants AF>=0.9',all(r['allele_fraction']>=0.9 for r in rows))
    c.true('all variants MQ>=40',all((r['mapping_quality'] or 99)>=40 for r in rows))
    # sub-clonal band must not overlap the fixed set
    c.true('subclonal AF strictly below 0.9',all(0.10<=s['allele_fraction']<0.90 for s in R['subclonal']))
    c.true('subclonal depth>=20',all(s['depth']>=20 for s in R['subclonal']))
    c.true('no position in both fixed and subclonal sets',
           not (set(r['pos'] for r in rows) & set(s['pos'] for s in R['subclonal'])))

    # ---- mask consistency
    mb=0
    for l in open(P+'/ref/mask_static.bed'):
        f=l.split('\t'); mb+=int(f[2])-int(f[1])
    c.eq('mask size',mb,490058)
    c.true('mask size stated in text','490,058' in txt)
    c.true('in_mask flag agrees with mask file',all(r['in_mask']==A.in_mask(r['pos']) for r in rows[:400]))

    # ---- catalogue
    c.eq('catalogue size',len(A.cat),946)
    c.true('catalogue size stated','946' in txt)
    cpres=collections.Counter(p['call'] for p in R['catalogue_presence'])
    c.eq('presence calls sum to catalogue',sum(cpres.values()),len(A.cat))
    c.true('no absent call with coverage',all(not(p['call']=='absent' and (p['breadth10'] or 0)>=0.95) for p in R['catalogue_presence']))

    # ---- coverage / exclusion
    nc=sum(b-a+1 for a,b in cov['noncallable'])
    br=100*(A.GL-nc)/A.GL
    c.eq('breadth recomputed',round(br,2),round(R['exclusion']['genome_breadth_10x_percent'],2),0.01)
    c.true('breadth stated in text','%.2f%%'%R['exclusion']['genome_breadth_10x_percent'] in txt)
    excl_expected=(br<90) or (R['purity']['percent_unmapped']>5)
    c.eq('exclusion decision',R['exclusion']['excluded'],excl_expected)
    c.true('mapped+unmapped == 100%',abs(R['purity']['percent_mapped_to_H37Rv']+R['purity']['percent_unmapped']-100)<0.011)

    # ---- lineage / spoligotype
    L=R['lineage']; S=R['spoligotype']
    c.true('lineage stated in text',L['lineage'] in txt)
    c.eq('spoligotype binary length',len(S['binary']),43)
    c.eq('spacers present matches binary',S['spacers_present'],S['binary'].count('1'))
    c.true('spoligotype binary in text',S['binary'] in txt.replace(' ',''))
    c.eq('SITVIT2 patterns searched',S['n_patterns_in_db'],9658)
    LF={'lineage1':'Indo-Oceanic','lineage2':'East-Asian','lineage3':'East-African-Indian',
        'lineage4':'Euro-American','lineage5':'West-Africa 1','lineage6':'West-Africa 2'}
    base=re.match(r'(lineage\d)',L['lineage'])
    if base and base.group(1) in LF:
        # tbdb qualifies some names, e.g. 'East-Asian (Beijing)'; the stem must match
        c.true('major lineage name matches barcode lineage',
               (L['major_lineage'] or '').startswith(LF[base.group(1)]),
               'got %r expected it to start with %r'%(L['major_lineage'],LF[base.group(1)]))
    FAM={'lineage1':('EAI',),'lineage2':('Beijing',),'lineage3':('CAS',),
         'lineage4':('LAM','T','S','X','H','Haarlem','Cameroon','Ural','TUR','URAL')}
    if base and S['family'] and base.group(1) in FAM:
        c.true('spoligotype family consistent with lineage',
               any(S['family'].upper().startswith(f.upper()) for f in FAM[base.group(1)]),
               'lineage %s but family %s'%(L['lineage'],S['family']))

    # ---- resistance
    for r in R['resistance']:
        c.true('resistance grade is an association grade: %s %s'%(r['gene'],r['mutation']),
               r['confidence'] in ('Assoc w R','Assoc w R - Interim') or r['type']=='drug_resistance')
    for r in R['resistance_uncertain']:
        c.true('uncertain grade correct: %s %s'%(r['gene'],r['mutation']),r['confidence']=='Uncertain significance')
    c.true('no resistance mutation duplicated into the uncertain table',
           not ({(r['gene'],r['mutation'],r['drug']) for r in R['resistance']} &
                {(r['gene'],r['mutation'],r['drug']) for r in R['resistance_uncertain']}))
    for r in R['resistance']: c.true('resistance row in text: %s'%r['mutation'],r['mutation'] in txt)

    # ---- epitopes
    pp=R['epitope_partition']
    for k in ('epitope','non_epitope'):
        c.true('%s partition has callable bp'%k,(pp[k]['callable_bp'] or 0)>0)
        if pp[k]['pi_per_kb'] is not None:
            c.eq('%s pi recomputed'%k,round(1000*pp[k]['variants']/pp[k]['callable_bp'],4),pp[k]['pi_per_kb'],0.001)
    ev=[r for r in rows if r['in_tcell_epitope']]
    c.true('every epitope variant has an IEDB id',all(r['epitope_ids'] for r in ev))
    c.true('every epitope variant has a codon',all(r['codon'] for r in ev))

    # ---- IS6110 / structural
    iz=R['is6110']
    c.eq('H37Rv IS6110 copies',iz['H37Rv_reference_copies'],16)
    c.true('all reported junctions are sequence-confirmed',all(j['is6110_matching']>=2 for j in iz['junctions']))
    c.eq('junction total consistent',iz['junctions_supported'],len(iz['junctions']))
    c.true('deletions are at least 1 kb',all(d['length']>=1000 for d in cov['deletions']))
    c.true('deletion count stated',str(len(cov['deletions'])) in txt)
    c.true('copy number positive',(R['is6110_copy_number']['estimated_copy_number'] or 0)>0)

    # ---- tables on disk match the report
    T=os.path.join(outdir,'tables')
    def nrows(f):
        p=os.path.join(T,'%s_%s.tsv'%(sample,f))
        return sum(1 for _ in open(p))-1 if os.path.exists(p) else -1
    c.eq('variants TSV rows',nrows('variants'),len(rows))
    c.eq('subclonal TSV rows',nrows('subclonal'),len(R['subclonal']))
    c.eq('gene_presence TSV rows',nrows('gene_presence'),len(A.cat))
    c.eq('deletions TSV rows',nrows('deletions'),len(cov['deletions']))
    c.eq('is6110 TSV rows',nrows('is6110_junctions'),len(iz['junctions']))
    c.eq('resistance TSV rows',nrows('resistance'),len(R['resistance'])+len(R['resistance_uncertain']))

    # ---- provenance statements that must be present
    for s_ in ['NC_000962.3','depth ≥ 10','0.9','Coll','WHO mutation catalogue version 2','VFDB','BV-BRC','IEDB',
               'Shell et al. 2015','Minch et al. 2015','frozen']:
        c.true('provenance mentioned: %s'%s_, s_ in txt)
    c.true('substitutions section lists Kraken2','Kraken2' in txt)
    c.true('substitutions section lists CheckM2','CheckM2' in txt)
    c.true('not-completed section lists TbD1','TbD1' in txt)
    c.true('independence from earlier analysis stated','Independent of any earlier analysis' in txt)
    c.true('assembly produced by this pipeline',R['assembly'].get('available') is True)
    c.true('assembly provenance states SPAdes','SPAdes 4.3.0' in R['assembly'].get('source',''))
    c.true('no assembly carried over','carried over' not in R['assembly'].get('source',''))
    c.eq('no ambiguous bases in assembly',R['assembly']['ambiguous_bases'],0)
    _a=R['assembly']
    if _a.get('available') and _a['contigs']==0:
        c.true('empty assembly is stated plainly in the report','No contigs passed the assembly filters' in txt)
    elif _a.get('available') and (_a['contigs']>500 or _a['N50']<20000):
        c.true('fragmented assembly is flagged in the report','Fragmented assembly' in txt)
    EXCLUDED=R['exclusion']['excluded']
    if not EXCLUDED:
        c.true('assembly length plausible for M. tuberculosis',
               4.0e6 < R['assembly']['total_length'] < 4.8e6,'%d bp'%R['assembly']['total_length'])
        c.true('assembly GC plausible',64.5 < R['assembly']['GC_percent'] < 66.5,'%.2f%%'%R['assembly']['GC_percent'])
    else:
        # an excluded isolate must say so plainly and must not be presented as a valid result
        c.true('report states the isolate is excluded','EXCLUDED' in txt)
        c.true('exclusion reasons given in the report',
               all(any(w in txt for w in rs.split()[:4]) for rs in R['exclusion']['reasons']))
        c.true('exclusion flagged on the cover','Meets inclusion criteria' in txt and 'NO' in txt)
    # cohort inclusion: every isolate is retained; weak ones must be flagged, never dropped or left unstated
    _ex=R['exclusion']
    _low=_ex['excluded'] and not (_ex.get('fails_only_on_foreign_reads') and _ex.get('usable_after_decontamination'))
    c.true('cover states the isolate is included in the cohort analyses',
           re.search(r'Included in cohort analyses\s+yes',txt) is not None)
    if _low:
        c.true('low-depth isolate is flagged LOW-DEPTH on the cover',
               re.search(r'Included in cohort analyses\s+yes\s+[—-]\s+flagged LOW-DEPTH',txt) is not None)
        c.true('low-depth inclusion is explained in section 4',
               'is included in the cohort analyses' in txt and 'LOW-DEPTH' in txt)
    else:
        c.true('isolate with adequate data is not flagged low-depth','LOW-DEPTH' not in txt)
    c.true('stored record agrees: every isolate is cohort-eligible',_ex.get('cohort_eligible') is True)
    # identification and validation
    I=R['identification']
    # isolate-specific notes (ref/sample_notes.json) must each appear in the report
    _nf=P+'/ref/sample_notes.json'
    for _n in (json.load(open(_nf)).get(sample,[]) if os.path.exists(_nf) else []):
        c.true('isolate note present: %s'%_n['title'],_n['title'] in txt and ' '.join(_n['text'].split()[:12]) in ' '.join(txt.split()))
    if _a.get('available') and _a['contigs']==0:
        c.true('ANI not computable (empty assembly) and the report states this',
               I['ani'] is None and 'Average nucleotide identity could not be computed' in txt)
    else:
        c.true('ANI computed',I['ani'] is not None)
    c.eq('identification targets assessed',len(I['markers']),4)
    if I['ani']: c.true('ANI is not impossibly high',I['ani']['ANI_percent']<=100)
    if not EXCLUDED:
        if I['ani']:
            c.true('ANI confirms M. tuberculosis (>=95%)',I['ani']['ANI_percent']>=95,
                   '%.3f%%'%I['ani']['ANI_percent'])
        c.true('all identification targets recovered',
               all(m['identity_percent'] is not None for m in I['markers']))
        c.true('identification targets match M. tuberculosis',
               all((m['identity_percent'] or 0)>=97 for m in I['markers']))
    else:
        # for an excluded isolate the identification evidence must be shown, not asserted away
        c.true('identification evidence reported for the excluded isolate',
               ('average nucleotide identity' in txt or 'Average nucleotide identity could not be computed' in txt)
               and 'Identification:' in txt)
    c.eq('RD regions assessed',len(I['rd']['regions']),7)
    c.true('complex member assigned',bool(I['rd']['complex_member']))
    c.true('RD9 call drives the species statement',
           ('sensu stricto' in I['rd']['complex_member']) ==
           (next(r['call'] for r in I['rd']['regions'] if r['region']=='RD9')=='present'))
    v=I['validation']
    c.true('insert size measured',(v['insert_size_mean'] or 0)>0)
    c.true('duplicate rate plausible',0<=(v['duplicate_rate_percent'] or 0)<=60)
    if not EXCLUDED:
        c.true('error rate plausible',0<(v['error_rate'] or 1)<0.1,'%s'%v['error_rate'])
    else:
        c.true('error rate reported for the excluded isolate',v['error_rate'] is not None)
    c.true('coverage evenness computed',(v['depth_cv'] or 0)>0)
    c.true('no patient identifier in the report',
           not any(k in txt for k in ['PATIENT NAME','NIKSHAY','PHONE NUMBER','ADDRESS']))
    if I.get('lab'):
        c.true('laboratory tests section present','Laboratory tests already performed' in txt)
    # TbD1 (analysis 28)
    t=R['tbd1']
    c.true('TbD1 assessed',t.get('available') is True)
    c.true('TbD1 call is one of present/partial/absent/not_determined',
           t.get('call') in ('present','partial','absent','not_determined'))
    c.true('TbD1 not derived from H37Rv mapping','no H37Rv mapping' in t.get('method',''))
    c.true('TbD1 reference is M. canettii','canettii' in t.get('reference',''))
    c.true('TbD1 fraction consistent with call',
           (t['call']=='present' and t['fraction_of_TbD1_recovered']>=0.80) or
           (t['call']=='partial' and 0.20<=t['fraction_of_TbD1_recovered']<0.80) or
           (t['call']=='absent' and t['fraction_of_TbD1_recovered']<0.20 and t.get('flank_fraction_covered',0)>=0.80) or
           (t['call']=='not_determined' and t['fraction_of_TbD1_recovered']<0.20 and t.get('flank_fraction_covered',1)<0.80))
    c.true('an absent TbD1 call is supported by the junction region being observed',
           t['call']!='absent' or t.get('flank_fraction_covered',0)>=0.80)
    if t['call']=='not_determined':
        c.true('TbD1 not-determined is stated plainly in the report','could not be determined for this isolate' in txt)
    c.true('TbD1 junction-region evidence is shown','junction region found on a single contig' in txt)
    _lin=R['lineage']['lineage']
    if (not EXCLUDED) and _lin.startswith(('lineage2','lineage3','lineage4')):
        c.true('modern lineage: TbD1 absent, or honestly not determined',t['call'] in ('absent','not_determined'))
    c.true('TbD1 section states Rv1557 is not a proxy','not used as a proxy' in txt)
    # PE/PPE (analysis 15)
    pe=R['pe_ppe']
    c.true('PE/PPE assessed',pe.get('available') is True)
    c.eq('PE/PPE loci count',pe.get('loci'),168)
    c.eq('PE/PPE rows match loci',len(pe.get('rows',[])),168)
    c.eq('PE/PPE outcomes sum to loci',sum(pe['summary'].values()),168)
    c.true('PE/PPE depth threshold is 20x',pe['thresholds']['min_read_depth']==20)
    c.true('every resolved PE/PPE locus meets depth >=20',
           all(r['mean_depth']>=20 for r in pe['rows'] if r['call']=='resolved'))
    c.true('every resolved PE/PPE locus is a single contig',
           all(r['contigs_aligned']==1 for r in pe['rows'] if r['call']=='resolved'))
    c.true('PE/PPE labelled secondary','secondary' in txt)
    run_panel_checks(c,R,txt,outdir,sample)
    return c

if __name__=='__main__':
    c=run(sys.argv[1],sys.argv[2])
    print('checks passed: %d'%c.ok)
    print('FAILURES: %d'%len(c.fail))
    for f in c.fail: print('   FAIL',f)
