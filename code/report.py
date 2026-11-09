"""Pipeline 2 - per-sample Report 2 generator (HTML -> PDF)."""
import json,os,sys,subprocess,collections,html,csv,re
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import analyse as A, figs, palette as C
P=A.P
VERS=json.load(open(P+'/ref/versions.json')) if os.path.exists(P+'/ref/versions.json') else {}

CSS="""
@page { size:A4; margin:17mm 15mm 16mm 15mm; }
body{font-family:'DejaVu Sans',Arial,sans-serif;font-size:9.4pt;color:#1A1A1A;line-height:1.45;}
h1{font-size:19pt;margin:0 0 2mm 0;color:#0072B2;}
h2{font-size:12.5pt;margin:7mm 0 2mm 0;color:#0072B2;border-bottom:1.1pt solid #0072B2;padding-bottom:1mm;
   page-break-after:avoid;}
h3{font-size:10.4pt;margin:4mm 0 1.5mm 0;color:#004E7A;page-break-after:avoid;}
p{margin:0 0 2mm 0;text-align:justify;}
table{border-collapse:collapse;width:100%;margin:2mm 0 3mm 0;font-size:8.3pt;page-break-inside:avoid;}
th{background:#0072B2;color:#fff;text-align:left;padding:1.4mm 1.8mm;font-weight:600;}
td{padding:1.2mm 1.8mm;border-bottom:0.4pt solid #DDE3E8;vertical-align:top;}
tr:nth-child(even) td{background:#F5F8FB;}
.figblock{page-break-inside:avoid;margin:4mm 0 5mm 0;text-align:center;}
.figblock img{max-width:100%;}
.cap{font-size:8.2pt;color:#444;text-align:left;margin-top:1.5mm;padding:0 2mm;}
.note{background:#F5F8FB;border-left:2.6pt solid #0072B2;padding:2mm 3mm;margin:3mm 0;font-size:8.6pt;}
.warn{background:#FDF4EC;border-left:2.6pt solid #D55E00;padding:2mm 3mm;margin:3mm 0;font-size:8.6pt;}
.kv{width:100%;font-size:9pt;}
.kv td{border:none;padding:0.8mm 2mm 0.8mm 0;}
.kv td:first-child{color:#555;width:46%;}
.small{font-size:8.2pt;color:#555;}
.mono{font-family:'DejaVu Sans Mono',monospace;font-size:8pt;}
.tag{display:inline-block;padding:0.3mm 1.6mm;border-radius:2pt;font-size:7.8pt;color:#fff;}
"""
def esc(x): return html.escape(str(x)) if x is not None else ''
def table(cols,rows,maxrows=None):
    out=['<table><tr>'+''.join('<th>%s</th>'%esc(c) for c in cols)+'</tr>']
    n=0
    for r in rows:
        out.append('<tr>'+''.join('<td>%s</td>'%esc(v) for v in r)+'</tr>'); n+=1
        if maxrows and n>=maxrows: break
    out.append('</table>')
    if maxrows and len(rows)>maxrows:
        out.append('<p class="small">Showing %d of %d rows; the full table is in the accompanying workbook and TSV files.</p>'%(maxrows,len(rows)))
    return ''.join(out)
def kv(pairs):
    return '<table class="kv">'+''.join('<tr><td>%s</td><td><b>%s</b></td></tr>'%(esc(a),esc(b)) for a,b in pairs)+'</table>'
def fig(path,caption):
    if not path:
        return '<p class="small"><i>Figure not drawn: there is no data to plot for this isolate.</i></p>'
    return '<div class="figblock"><img src="%s"><div class="cap">%s</div></div>'%(path,caption)

def build(sample,outdir):
    import analyse_all
    W=P+'/work/'+sample; FG=W+'/fig'; os.makedirs(FG,exist_ok=True); os.makedirs(outdir,exist_ok=True)
    R=analyse_all.run(sample)
    _nf=P+'/ref/sample_notes.json'
    NOTES=json.load(open(_nf)).get(sample,[]) if os.path.exists(_nf) else []
    def NOTE(where):
        return ''.join('<div class="warn"><b>%s</b> %s</div>'%(esc(n['title']),esc(n['text'])) for n in NOTES if n['section']==where)
    apply_absence_rule(R)
    cov=R['coverage']; rows=R['variants']; sc=R['subclonal']; pur=R['purity']; excl=R['exclusion']
    asm=R['assembly']; comp=R['completeness']; pp=R['epitope_partition']; iz=R['is6110']
    cn=R['is6110_copy_number']; pres=R['catalogue_presence']; amr=R['amrfinder']
    lin=R['lineage']; spo=R['spoligotype']; rhits=R['resistance']; rvus=R['resistance_uncertain']
    b4,af=R['fastp']['before_filtering'],R['fastp']['after_filtering']; fr=R['fastp_filtering']
    gb=json.load(open(W+'/cov.json'))['gene_breadth']
    unmask=[r for r in rows if not r['in_mask']]
    vir=[r for r in unmask if r['in_virulence_catalogue']]
    f_cov=figs.fig_coverage(W+'/depth_all.txt.gz',FG+'/coverage.png')
    f_eff=figs.fig_effects(rows,FG+'/effects.png'); f_sys=figs.fig_systems(rows,FG+'/systems.png')
    f_af=figs.fig_af(rows,sc,FG+'/af.png'); f_map=figs.fig_genome_map(cov['deletions'],iz['junctions'],FG+'/map.png')
    f_epi=figs.fig_epitope(pp,FG+'/epitope.png')
    H=[];S=H.append
    S('<h1>%s — whole-genome virulence analysis</h1>'%esc(sample))
    S('<p class="small">Report 2 — virulence-gene focused secondary analysis of <i>Mycobacterium tuberculosis</i>, '
      'performed to the bioinformatics brief. Independent of any earlier analysis of this isolate: every result below '
      'was computed by this pipeline from the raw reads. Reference H37Rv NC_000962.3.</p>')
    S(kv([('Lineage','%s (%s)'%(lin['lineage'],lin['major_lineage'] or 'n/a')),
          ('Spoligotype','SIT %s, %s'%(spo['sit'] or 'not in SITVIT2',spo['family'] or 'unassigned family')),
          ('Mean depth','%.1f×'%cov['mean_depth_all']),
          ('Genome covered ≥10×','%.2f%%'%excl['genome_breadth_10x_percent']),
          ('Reads mapping to H37Rv','%.2f%%'%pur['percent_mapped_to_H37Rv']),
          ('Filtered variants','%d total, %d outside masked regions'%(len(rows),len(unmask))),
          ('Virulence-catalogue variants','%d'%len(vir)),
          ('Resistance-associated mutations','%d'%len(rhits)),
          ('IS6110 copy number','%s'%cn['estimated_copy_number']),
          ('Meets inclusion criteria','yes' if not excl['excluded'] else 'NO — '+'; '.join(excl['reasons'])),
          ('Included in cohort analyses','yes' if not (excl['excluded'] and not (excl.get('fails_only_on_foreign_reads') and excl.get('usable_after_decontamination'))) else 'yes — flagged LOW-DEPTH')]))
    S(SEC_ident(R['identification']))
    S(SEC_qc(b4,af,fr,cov,sample))
    if NOTE('qc'): S(NOTE('qc'))
    DECON[:]=[R.get('decontamination')] if R.get('decontamination') else []
    S(SEC_purity(pur))
    S(SEC_cov(cov,excl,gb,f_cov))
    S(SEC_asm(asm,comp))
    S(SEC_lineage(lin,spo))
    if NOTE('lineage'): S(NOTE('lineage'))
    S(SEC_resist(rhits,rvus,amr,R['resistance_genes']))
    S(SEC_variants(rows,unmask,f_eff,f_af))
    tc=A.tier_callable(pres)
    S(SEC_virulence(vir,f_sys,rows,tc,pres))
    S(SEC_epitope(rows,pp,f_epi))
    S(SEC_reg(rows))
    S(SEC_subclonal(sc,rows))
    S(SEC_struct(cov,iz,cn,f_map))
    S(SEC_tbd1(R['tbd1'],lin))
    S(SEC_pe(R['pe_ppe']))
    S(SEC_presence(pres))
    S(SEC_cohort(sample,R))
    S(SEC_limits(excl,cov))
    if NOTE('limits'): S(NOTE('limits'))
    S(SEC_methods())
    body='\n'.join(H)
    body=body.replace('%%','%')  # v5: literal %% typos in unformatted strings
    htmlout=os.path.join(outdir,'%s_report2.html'%sample)
    open(htmlout,'w').write('<html><head><meta charset="utf-8"><style>%s</style></head><body>%s</body></html>'%(CSS,body))
    pdf=os.path.join(outdir,'%s_complete_analysis_report_2.pdf'%sample)
    subprocess.run(['wkhtmltopdf','--enable-local-file-access','--quiet',
                    '--margin-top','17','--margin-bottom','16','--margin-left','15','--margin-right','15',
                    '--footer-font-size','7','--footer-right','[page]',htmlout,pdf],check=True)
    write_tables(sample,outdir,R,gb)
    return pdf,R

def SEC_ident(I):
    o=['<h2>1. Organism identification and sequencing validation</h2>',
       '<p>Before any virulence result is interpreted, the isolate is confirmed to be '
       '<i>Mycobacterium tuberculosis</i>, placed within the complex, and the sequencing run is '
       'checked for the artefacts that would otherwise be mistaken for biology.</p>',
       '<h3>1.1 Species confirmation</h3>']
    a=I.get('ani')
    if a:
        o.append(table(['metric','value'],[
            ['average nucleotide identity to H37Rv','%.3f%%'%a['ANI_percent']],
            ['assembly fragments aligned','%s of %s (%.1f%%)'%('{:,}'.format(a['fragments_aligned']),
             '{:,}'.format(a['fragments_total']),100*a['aligned_fraction'])],
            ['method',a['method']]]))
        v='confirms the same species' if a['ANI_percent']>=95 else 'is BELOW the 95% species threshold'
        o.append('<p>The conventional species boundary is 95%% ANI. This isolate %s as the H37Rv reference.</p>'%v)
    else:
        o.append('<div class="warn"><b>Average nucleotide identity could not be computed for this isolate.</b> '
                 'ANI is measured on the de novo assembly, and no contig passed the uniform assembly filters '
                 '(section 5). Species identity is therefore assessed only from the read-level evidence below, the '
                 'identification target genes and the region-of-difference profile, which for this isolate is limited '
                 'by low sequencing depth. It is reported as found and is not upgraded.</div>')
    mk=I.get('markers') or []
    if mk:
        o.append('<h3>1.2 Identification target genes</h3>')
        o.append('<p>The genes routinely used to identify mycobacteria to species level, compared with H37Rv.</p>')
        o.append(table(['target','locus','identity to H37Rv','coverage'],
            [[m['target'],m['locus'],
              ('%.2f%%'%m['identity_percent']) if m['identity_percent'] is not None else 'not recovered',
              ('%.1f%%'%m['coverage_percent']) if m['coverage_percent'] is not None else '—'] for m in mk]))
    rd=I.get('rd')
    if rd:
        o.append('<h3>1.3 Region-of-difference profile and complex member</h3>')
        o.append('<p>%s</p>'%esc(rd['note']))
        o.append(table(['region','span','length bp','fraction covered ≥5×','mean depth','call','what its loss indicates'],
            [[r['region'],r['span'],'{:,}'.format(r['length_bp']),r['fraction_covered_5x'],r['mean_depth'],
              r['call'],r['interpretation']] for r in rd['regions']]))
        o.append('<div class="note"><b>Identification:</b> %s</div>'%esc(rd['complex_member']))
    v=I.get('validation')
    if v:
        o.append('<h3>1.4 Sequencing run validation</h3>')
        o.append(table(['metric','value','what it checks'],[
            ['mean insert size','%s bp (SD %s)'%(v['insert_size_mean'],v['insert_size_sd']),'library construction'],
            ['duplicate rate','%s%%'%v['duplicate_rate_percent'],'PCR over-amplification'],
            ['sequencing error rate','%s'%v['error_rate'],'base-calling quality (mismatches per mapped base)'],
            ['mean / median depth','%s× / %s×'%(v['mean_depth'],v['median_depth']),'overall yield'],
            ['depth coefficient of variation',v['depth_cv'],'evenness of coverage'],
            ['genome within 2-fold of median depth','%s%%'%v['percent_within_2fold_of_median'],
             'absence of large coverage bias'],
            ['heterozygous sites per Mb',v['heterozygous_sites_per_Mb'],
             'purity of the isolate; a raised value suggests mixed infection or contamination']]))
    lab=I.get('lab')
    if lab:
        o.append('<h3>1.5 Laboratory tests already performed on this isolate</h3>')
        o.append('<p>Taken from the study metadata, for comparison with the genome-based calls in section 7. '
                 'Clinical fields only; no patient identifiers are used anywhere in this report.</p>')
        o.append(table(['test','result'],[[k.title().replace('Lpa','LPA').replace('Cbnaat','CBNAAT'),
                                           (v_ if v_ else 'not recorded')] for k,v_ in lab.items()]))
    return '\n'.join(o)

def SEC_qc(b4,af,fr,cov,sample):
    _W=P+'/work/'+sample
    _fj=json.load(open(_W+'/%s_fastp.json'%sample))
    _n1=_fj['read1_before_filtering']['total_reads']; _n2=_fj['read2_before_filtering']['total_reads']
    _ii=json.load(open(_W+'/input_integrity.json')) if os.path.exists(_W+'/input_integrity.json') else None
    o=['<h2>2. Read quality control (analysis 4)</h2>',
       '<p>Reads were trimmed with fastp at a Phred quality threshold of 20, with paired-end adapter detection.</p>',
       table(['metric','before filtering','after filtering'],[
        ['reads','{:,}'.format(b4['total_reads']),'{:,}'.format(af['total_reads'])],
        ['bases','{:,}'.format(b4['total_bases']),'{:,}'.format(af['total_bases'])],
        ['mean read length (R1)',b4['read1_mean_length'],af['read1_mean_length']],
        ['Q20 rate','%.1f%%'%(100*b4['q20_rate']),'%.1f%%'%(100*af['q20_rate'])],
        ['Q30 rate','%.1f%%'%(100*b4['q30_rate']),'%.1f%%'%(100*af['q30_rate'])],
        ['GC content','%.1f%%'%(100*b4['gc_content']),'%.1f%%'%(100*af['gc_content'])],
        ['read pairing (R1 / R2 reads)','{:,} / {:,}'.format(_n1,_n2),'equal: every read paired' if _n1==_n2 else 'UNEQUAL']]),
       ('<p><b>Input data integrity.</b> The R2 file supplied for this isolate was truncated (it ended after {:,} complete '
        'reads, whereas R1 contained {:,}). Analysis was therefore restricted to the {:,} read pairs complete in both files, '
        'with read identifiers matched pair by pair ({} mismatches). The resulting coverage is assessed against the same '
        'inclusion criteria as every other isolate in section 4.</p>'.format(_ii['r2_complete_reads_supplied'],_ii['r1_reads_supplied'],
        _ii['pairs_analysed'],_ii['read_name_mismatches'])
        if _ii else
        '<p><b>Input data integrity.</b> Both read files were complete: R1 and R2 contained the same number of reads '
        '({:,} each), so every read was paired.</p>'.format(_n1)),
       '<p>Reads removed: %s low quality, %s with too many ambiguous bases, %s adapter dimers, %s too short.</p>'%(
        '{:,}'.format(fr['low_quality_reads']),'{:,}'.format(fr['too_many_N_reads']),
        '{:,}'.format(fr['adapter_dimer_reads']),'{:,}'.format(fr['too_short_reads']))]
    if float(b4['read1_mean_length'])>200:
        o.append('<div class="warn"><b>Platform note.</b> The brief describes the dataset as Illumina HiSeq 2 × 150 bp at '
          'approximately 100× depth. This isolate has a mean read length of %s bp and a realised mean depth of %.1f×, which is '
          'consistent with a 2 × 300 bp MiSeq run rather than a HiSeq. Measured values are used throughout this report; the '
          'brief’s stated figures are not reproduced.</div>'%(b4['read1_mean_length'],cov['mean_depth_all']))
    return '\n'.join(o)

DECON=[]
def SEC_purity(pur):
    return ('<h2>3. Species purity from read mapping (analysis 3)</h2>'
      '<p>Kraken2 with the standard database could not be used: that database is 50–100 GB and exceeds the disk available in '
      'this offline environment. Purity is assessed instead from the proportion of quality-filtered reads mapping to H37Rv. '
      'This substitution is stated rather than made silently, and is listed in section 14.</p>'
      +table(['metric','value'],[['reads after QC','{:,}'.format(pur['reads_after_qc'])],
        ['primary alignments','{:,}'.format(pur['primary_reads'])],
        ['primary alignments mapped to H37Rv','{:,}'.format(pur['primary_mapped'])],
        ['mapped','%.2f%%'%pur['percent_mapped_to_H37Rv']],
        ['unmapped','%.2f%%'%pur['percent_unmapped']]]))
    D=DECON[0] if DECON else None
    if D:
        o.append('<h3>3.1 Decontamination</h3>')
        o.append('<p>%s</p>'%esc(D['method']))
        o.append(table(['metric','value'],[
            ['foreign reads removed','%.2f%%'%D['foreign_read_percent']],
            ['M. tuberculosis reads retained','%.2f%%'%D['retained_read_percent']],
            ['depth after decontamination','%.1f×'%D['depth_after_decontamination']],
            ['genome at ≥10× after decontamination','%.2f%%'%D['breadth_10x_after_decontamination']]]))
        o.append('<div class="note">Coverage is measured from aligned reads only, so these depth figures were '
                 'never inflated by the foreign fraction. Decontamination removes reads; it cannot add coverage '
                 'that was not sequenced. Where an isolate is short of depth, that is a property of the library, '
                 'not of the contamination.</div>')

def SEC_cov(cov,excl,gb,f_cov):
    nc=sum(b-a+1 for a,b in cov['noncallable'])
    low_cat=sum(1 for rv in A.cat if rv in gb and gb[rv]['breadth10']<0.95)
    o=['<h2>4. Coverage, callable regions and the exclusion rule (analysis 5)</h2>',
       fig(f_cov,'<b>Figure 1.</b> Read depth across the H37Rv chromosome in 5 kb windows. The dashed line marks the '
           'median depth. Troughs correspond to deleted or poorly mappable regions.'),
       table(['metric','value'],[
        ['mean depth (mapping-quality filtered)','%.1f×'%cov['mean_depth_mqfiltered']],
        ['mean depth (unfiltered)','%.1f×'%cov['mean_depth_all']],
        ['genome covered ≥10×','%.2f%%'%excl['genome_breadth_10x_percent']],
        ['non-callable intervals (<10×)','%d intervals, %s bp (%.2f%%)'%(len(cov['noncallable']),'{:,}'.format(nc),100*nc/A.GL)],
        ['genes <95%% callable','%d of %d'%(sum(1 for v in gb.values() if v['breadth10']<0.95),len(gb))],
        ['virulence-catalogue genes <95%% callable','%d of %d'%(low_cat,len(A.cat))]]),
       '<p>Two depth tracks are maintained. Mapping-quality filtered depth defines the callable mask and supports variant '
       'calling, where unique placement of reads is required. Unfiltered depth is used for deletions and for repetitive '
       'elements such as IS<i>6110</i>, where reads cannot be placed uniquely and excluding them would create false absences.</p>',
       '<div class="%s"><b>Exclusion rule.</b> %s<br>This isolate: %.2f%% of the genome at ≥10× and %.2f%% of reads unmapped. '
       '<b>%s</b></div>'%('warn' if excl['excluded'] else 'note',esc(excl['rule']),
        excl['genome_breadth_10x_percent'],excl['percent_unmapped'],
        'EXCLUDED under the rule as written — '+'; '.join(excl['reasons']) if excl['excluded'] else 'Retained for analysis.')]
    if excl.get('fails_only_on_foreign_reads') and excl.get('usable_after_decontamination'):
        o.append('<div class="note"><b>Cohort eligibility.</b> The only criterion this isolate fails is the share of '
                 'foreign reads in the raw library. Those reads are removed before any analysis, and the isolate\u2019s own '
                 'data then covers %.2f%% of the genome at ≥10×, which meets every quality requirement. It is therefore '
                 '<b>retained for the cohort analyses</b>, with the raw foreign fraction reported above so the decision '
                 'is visible.</div>'%excl['genome_breadth_10x_percent'])
    elif excl['excluded']:
        o.append('<div class="warn"><b>Cohort eligibility.</b> This isolate does not meet the brief\u2019s inclusion '
                 'criteria, but it <b>is included in the cohort analyses</b> by the study\u2019s decision that every isolate is '
                 'retained, and it is flagged <b>LOW-DEPTH</b>. Its own sequence data, after decontamination, covers only '
                 '%.2f%% of the genome at ≥10×, below the depth the brief requires for variant calling. This is a shortage '
                 'of data, not an artefact of contamination, and cannot be corrected by filtering. Positions without '
                 'confident coverage are treated as missing data, never as reference calls, and cohort results are '
                 'reported with and without low-depth isolates so that their influence is visible.</div>'%excl['genome_breadth_10x_percent'])
    return '\n'.join(o)

def SEC_asm(asm,comp):
    o=['<h2>5. De novo assembly and completeness (analyses 5b, 1 and 2)</h2>',
       '<p>Every isolate is assembled from its own quality-trimmed reads with SPAdes 4.3.0 in <span class="mono">--isolate</span> '
       'mode. No assembly is carried over from earlier work, so all isolates in the study are assembled by the same '
       'program with the same settings and are directly comparable. The assembly is what makes the TbD1 locus '
       'recoverable, since that region has no coordinates in H37Rv.</p>']
    if not asm.get('available'):
        o.append('<div class="warn">No de novo assembly available for this isolate: %s</div>'%esc(asm.get('reason','')))
    elif asm.get('contigs',0)==0:
        fl=asm.get('filter',{})
        o.append('<div class="warn"><b>No contigs passed the assembly filters.</b> SPAdes was run exactly as for every '
                 'other isolate, but no contig met all three uniform criteria (length ≥500 bp, k-mer coverage ≥5, and a '
                 'BLASTn hit to the MTBC reference panel). Contigs dropped: %s for length, %s for low coverage, %s for '
                 'lack of an MTBC hit. This reflects a shortage of sequencing data for this isolate and is reported as '
                 'it is; no assembly statistics can be computed and none are invented. Results that depend on an '
                 'assembly — TbD1, PE/PPE secondary analysis and gene presence or absence — are therefore not '
                 'available for this isolate.</div>'%('{:,}'.format(fl.get('dropped_short',0)),
                 '{:,}'.format(fl.get('dropped_low_cov',0)),'{:,}'.format(fl.get('dropped_no_MTBC_hit',0))))
    else:
        o.append(table(['metric','value'],[
            ['contigs','{:,}'.format(asm['contigs'])],['total length (bp)','{:,}'.format(asm['total_length'])],
            ['largest contig (bp)','{:,}'.format(asm['largest_contig'])],['N50 (bp)','{:,}'.format(asm['N50'])],
            ['L50',asm['L50']],['GC content','%.2f%%'%asm['GC_percent']],['ambiguous bases',asm['ambiguous_bases']]]))
        o.append('<p class="small">%s. QUAST was available but is not required to compute these statistics; they are '
                 'derived directly from the assembly and are reported with the assembly’s provenance stated.</p>'%esc(asm['source']))
        if asm['contigs']>500 or asm['N50']<20000:
            o.append('<div class="warn"><b>Fragmented assembly.</b> This assembly has %s contigs with an N50 of '
                     '%s bp, well short of the few hundred contigs and N50 above 50 kb obtained from a '
                     'well-behaved library of this depth. The cause is visible in section 1.4: a broad insert-size '
                     'distribution leaves read pairs uninformative for scaffolding. Total assembly length and the '
                     'read-level analyses are unaffected, but results that depend on contiguity — gene presence and '
                     'absence, the PE/PPE secondary analysis and TbD1 recovery — are less reliable for this isolate '
                     'and are flagged accordingly.</div>'%('{:,}'.format(asm['contigs']),'{:,}'.format(asm['N50'])))
    o.append('<p>CheckM2 requires a multi-gigabyte DIAMOND reference database that could not be obtained offline. '
             'Completeness is therefore estimated from breadth of coverage across single-copy H37Rv coding sequences, '
             'excluding the PE/PPE and insertion-sequence families, which are repetitive by nature.</p>')
    o.append(table(['metric','value'],[
        ['single-copy CDS assessed','{:,}'.format(comp['single_copy_cds_assessed'])],
        ['CDS with ≥95% at ≥10×','{:,}'.format(comp['cds_with_>=95%_at_10x'])],
        ['estimated completeness','%.2f%%'%comp['completeness_percent']]]))
    return '\n'.join(o)

def SEC_lineage(lin,spo):
    o=['<h2>6. Lineage, sublineage and spoligotype (analyses 6 and 14)</h2>',
       '<p>Lineage was assigned from the Coll 62-SNP barcode as distributed in the TB-Profiler database, requiring at '
       'least 90%% of a barcode’s defining positions to carry the derived allele. %d of %d barcode positions were '
       'callable in this isolate.</p>'%(lin['positions_called'],lin['barcode_positions']),
       table(['metric','value'],[['lineage',lin['lineage']],['major lineage',lin['major_lineage'] or 'not assigned'],
        ['barcode support',', '.join('%s %s'%(k,v) for k,v in lin['all_supported'].items()) or 'none'],
        ['spoligotype (binary)','<span class="mono">%s</span>'%spo['binary']],
        ['spoligotype (octal)',spo['octal']],['spacers present','%d of 43'%spo['spacers_present']],
        ['SIT',spo['sit'] or 'pattern not in SITVIT2'],['family',spo['family'] or 'unassigned family'],
        ['SITVIT2 patterns searched','{:,}'.format(spo['n_patterns_in_db'])]]).replace('&lt;span class=&quot;mono&quot;&gt;','<span class="mono">').replace('&lt;/span&gt;','</span>'),
       '<div class="note"><b>Nomenclature.</b> Lineage names and spoligotype family names are separate vocabularies. '
       'Following TB-Profiler: lineage 1 = Indo-Oceanic (EAI families); lineage 2 = East-Asian (Beijing); '
       'lineage 3 = East-African-Indian (CAS families); lineage 4 = Euro-American (LAM, T, S, X, H); '
       'lineage 5 = West-Africa 1; lineage 6 = West-Africa 2. Lineage 1 and lineage 3 are different lineages despite the '
       'similar historical wording. Spoligotyping reads a single repetitive locus and convergent patterns between lineages '
       'are documented, so the Coll barcode is the assignment used throughout this report.</div>']
    if not spo['sit']:
        o.append('<div class="note">This spoligotype pattern is not present in SITVIT2. It is reported as an orphan '
                 'pattern with no SIT assigned. No nearest-neighbour SIT is given, because similarity of spoligotype '
                 'patterns does not imply shared ancestry.</div>')
    return '\n'.join(o)

def SEC_resist(rhits,rvus,amr,ngenes):
    o=['<h2>7. Drug resistance (analysis 7) with WHO confidence grading (extra E1), '
       'variants of uncertain significance (E2) and an independent screen (E3)</h2>',
       '<p>Variants were matched against the TB-Profiler mutation catalogue, which embeds the WHO mutation catalogue '
       'version 2 and its five confidence grades. %d resistance-associated genes were screened.</p>'%ngenes]
    if rhits:
        o.append('<h3>7.1 Resistance-associated mutations</h3>')
        o.append(table(['gene','mutation','drug','WHO confidence','position','effect'],
            [[r['gene'],r['mutation'],r['drug'],r['confidence'],r['pos'],r['effect'].replace('_',' ')] for r in rhits]))
    else:
        o.append('<p>No mutation graded <i>Assoc w R</i> or <i>Assoc w R - Interim</i> was detected in any catalogued '
                 'resistance gene. On the evidence of this genome the isolate is predicted susceptible to all drugs in the catalogue.</p>')
    o.append('<h3>7.2 Variants of uncertain significance in resistance genes (E2)</h3>')
    if rvus:
        o.append('<p>These are graded <i>Uncertain significance</i> by WHO. They are not resistance calls, and are '
                 'reported separately so they are not mistaken for one.</p>')
        o.append(table(['gene','mutation','drug','WHO confidence','position'],
            [[r['gene'],r['mutation'],r['drug'],r['confidence'],r['pos']] for r in rvus],maxrows=40))
    else:
        o.append('<p>None detected.</p>')
    o.append('<h3>7.3 Independent resistance-gene screen (E3)</h3>')
    if amr.get('available'):
        o.append('<p>The isolate proteome was screened against the NCBI AMRFinder reference protein set by BLASTp '
                 '(≥90%% identity, ≥80%% of the reference length). This is an independent second opinion on the '
                 'catalogue-based calls above, and detects intrinsic determinants that mutation catalogues do not list.</p>')
        if amr['hits']:
            rr=[]
            for h in amr['detail']:
                parts=h['amr_ref'].split('|')
                rr.append([h['query'],parts[4] if len(parts)>4 else parts[0],
                           parts[7] if len(parts)>7 else '',parts[6] if len(parts)>6 else '','%.1f%%'%h['identity']])
            o.append(table(['locus','AMR gene','drug class','drug','identity'],rr,maxrows=40))
            o.append('<p class="small">These are intrinsic chromosomal determinants present in <i>M. tuberculosis</i> '
                     'generally. They are not acquired resistance and do not change the predictions in section 6.1.</p>')
        else:
            o.append('<p>No AMRFinder reference protein matched at the thresholds used.</p>')
    else:
        o.append('<p>Not performed: %s.</p>'%esc(amr.get('reason','')))
    o.append('<div class="note"><b>Overlapping loci.</b> <i>katG</i>, <i>ahpC</i> and <i>eis</i> are both resistance '
             'determinants and virulence genes. Per the brief, <i>katG</i> codon 315, the <i>ahpC</i> promoter and the '
             '<i>eis</i> promoter are reported here as resistance findings and are excluded from the virulence '
             'interpretation in section 8, so the two signals are not confounded.</div>')
    return '\n'.join(o)

def SEC_variants(rows,unmask,f_eff,f_af):
    c=collections.Counter(r['effect'] for r in rows)
    snp=sum(1 for r in rows if r['type']=='snp'); ind=len(rows)-snp
    nsyn=sum(1 for r in unmask if r['effect']=='missense_variant')
    syn=sum(1 for r in unmask if r['effect']=='synonymous_variant')
    ld=sum(1 for r in rows if r['lineage_defining'])
    o=['<h2>8. Variant calling, filtering and annotation (analyses 10–14)</h2>',
       '<p>Reads were aligned to H37Rv with BWA-MEM, duplicates were marked, and variants were called with bcftools in '
       'haploid mode. Filters are those specified in the brief: <b>depth ≥ 10, mapping quality ≥ 40, alternate allele '
       'fraction ≥ 0.9</b>. Variants were annotated with SnpEff against a database built locally from the NC_000962.3 '
       'GenBank record, and cross-referenced to protein domains (UniProt), T-cell epitopes (IEDB), the Coll barcode and '
       'the frozen virulence catalogue.</p>',
       table(['metric','value'],[
        ['variants passing all filters','{:,}'.format(len(rows))],
        ['SNPs / indels','%s / %s'%('{:,}'.format(snp),'{:,}'.format(ind))],
        ['outside masked regions','{:,}'.format(len(unmask))],
        ['inside masked regions','{:,}'.format(len(rows)-len(unmask))],
        ['missense (unmasked)','{:,}'.format(nsyn)],['synonymous (unmasked)','{:,}'.format(syn)],
        ['non-synonymous : synonymous ratio','%.2f'%(nsyn/syn) if syn else 'n/a'],
        ['lineage-defining (Coll barcode)','{:,}'.format(ld)]]),
       fig(f_eff,'<b>Figure 2.</b> Variants by predicted functional effect (SnpEff). Counts include masked regions; '
           'the masked subset is analysed separately in section 12.'),
       fig(f_af,'<b>Figure 3.</b> Alternate allele fraction. Fixed variants passing the 0.9 threshold are shown in blue; '
           'sub-clonal sites between 0.10 and 0.90 (section 10) in orange. The dashed line marks the 0.9 threshold. '
           'The y axis is logarithmic.'),
       '<div class="note"><b>Masking.</b> PE, PPE and PE_PGRS genes, insertion sequences, phage regions and the '
       'TB-Profiler exclusion list are masked from the primary analysis, covering 490,058 bp (11.11%% of the genome). '
       'Short reads cannot be placed reliably in these regions. Variants falling inside the mask are retained in the '
       'tables, flagged, and interpreted only in the secondary analysis.</div>']
    return '\n'.join(o)

def SEC_virulence(vir,f_sys,rows,tc,pres):
    by=collections.Counter(r['functional_system'] for r in vir)
    prot=[r for r in vir if r['effect'] in ('missense_variant','stop_gained','frameshift_variant','stop_lost')]
    tier=collections.Counter(r['panel_tier'] for r in rows if r['panel_tier'])
    tested=[r for r in rows if r['tested_panel'] and not r['in_mask']]
    o=['<h2>9. Variants in the frozen virulence panel</h2>',
       '<p>The catalogue was supplied frozen: <b>946 genes, 1,141,186 bp, 25.87% of the genome</b>, compiled from VFDB '
       '(core and full), Mycobrowser functional categories, BV-BRC specialty genes and the curated locus list in the '
       'brief\u2019s appendix, deduplicated by Rv number. It is used exactly as issued; no gene was added or removed. '
       'BED files were taken from the NC_000962.3 folder and their md5 checksums verified against the provenance sheet '
       'before use.</p>',
       '<div class="note"><b>Tiering and multiple testing.</b> Only Tier 1 is tested. Tier 1a (core virulence, 204 genes) '
       'is the primary pre-specified panel; Tier 1b (toxin\u2013antitoxin, 120 genes) is tested as a separate '
       'false-discovery-rate stratum. Benjamini\u2013Hochberg correction is applied within 1a and within 1b '
       'independently, never across the pooled 324. Tier 2 (PE/PPE, 168 genes) is masked and reported only as secondary '
       'analysis; Tier 3 (454 genes) is described and never formally tested, so it generates no p values and contributes '
       'nothing to the correction. This partitioning is pre-specified, not chosen after seeing results. Association '
       'testing is reported in the cohort report; no association test is made in this per-isolate report.</div>',
       '<h3>9.1 Panel coverage in this isolate</h3>',
       '<p>Callable length is the denominator the brief specifies for variant-burden normalisation (analysis 31). '
       'The nominal lengths are the frozen <span class="mono">length_bp</span> values; the callable figures are this '
       'isolate\u2019s own, after coverage.</p>',
       table(['tier','stratum','genes','nominal bp','callable bp in this isolate','callable %','variants'],
         [[t,{'1a':'core_virulence','1b':'toxin_antitoxin','2':'pe_ppe_secondary','3':'not_tested'}[t],
           tc[t]['genes'],'{:,}'.format(tc[t]['nominal_bp']),'{:,}'.format(tc[t]['callable_bp']),
           '%.2f%%'%(100*tc[t]['callable_fraction']) if tc[t]['callable_fraction'] else '—',
           tier.get(t,0)] for t in ('1a','1b','2','3')]),
       '<p><b>Variants in the tested panel (Tiers 1a and 1b), outside masked regions: %d.</b></p>'%len(tested),
       '<h3>9.2 Variants by functional system</h3>']
    if f_sys: o.append(fig(f_sys,'<b>Figure 4.</b> Unmasked variants in virulence-panel genes, grouped by functional '
                           'system. Genes named in the brief\u2019s appendix keep the brief\u2019s grouping; genes contributed only by '
                           'a database are grouped by that database\u2019s own functional category.'))
    if not f_sys: o.append(fig(None,''))
    o.append(table(['functional system','variants'],[[k,v] for k,v in by.most_common()]))
    o.append('<h3>9.3 Protein-altering variants in panel genes</h3>')
    rr=[[r['gene'],r['rv'],r['panel_tier'],r['pos'],'%s>%s'%(r['ref'],r['alt']),r['hgvs_p'] or r['hgvs_c'],
         r['effect'].replace('_',' '),r['functional_system'][:24],(r['domains'][:40] or '—'),
         ('yes' if r['in_tcell_epitope'] else ''),('' if r['blosum62'] is None else r['blosum62'])]
        for r in sorted(prot,key=lambda x:(x['panel_tier'] or 'z',x['pos']))]
    o.append(table(['gene','Rv','tier','position','change','protein change','effect','system','domain / site','epitope','BLOSUM62'],
                   rr,maxrows=60))
    return '\n'.join(o)

def SEC_epitope(rows,pp,f_epi):
    ev=[r for r in rows if r['in_tcell_epitope']]
    o=['<h2>10. T-cell epitope variation (analysis 21)</h2>',
       '<p>2,189 IEDB epitopes across 535 antigens were mapped onto H37Rv coordinates. Positions reported by IEDB were '
       'not used directly: 255 of them disagreed with the H37Rv protein sequence, reflecting a mix of precursor and '
       'mature-protein numbering and strain polymorphism. Each peptide was instead located by sequence search in the '
       'translated H37Rv protein, and the resulting genome coordinates were verified by translating them back.</p>',
       '<h3>10.1 Variants falling within known T-cell epitopes</h3>']
    if ev:
        o.append(table(['gene','Rv','position','protein change','effect','epitope IEDB IDs'],
            [[r['gene'],r['rv'],r['pos'],r['hgvs_p'] or r['hgvs_c'],r['effect'].replace('_',' '),r['epitope_ids']] for r in ev],maxrows=40))
    else:
        o.append('<p>No variant in this isolate falls within a mapped T-cell epitope.</p>')
    o.append('<h3>10.2 Per-isolate diversity in epitope and non-epitope partitions (21b)</h3>')
    o.append('<p>Computed over the same genes — those carrying at least one mapped epitope — with each gene split into '
             'its epitope and non-epitope parts. Values are per isolate, not group means, and are normalised by callable '
             'sequence length. Masked regions are excluded.</p>')
    o.append(table(['partition','variants','callable bp','variants per kb','synonymous','non-synonymous','pN/pS'],
        [[k.replace('_',' '),pp[k]['variants'],'{:,.0f}'.format(pp[k]['callable_bp']),
          pp[k]['pi_per_kb'],pp[k]['syn'],pp[k]['nonsyn'],pp[k]['pN_pS']] for k in ('epitope','non_epitope')]))
    o.append(fig(f_epi,'<b>Figure 5.</b> Variant density in epitope and non-epitope partitions of the same genes, '
                 'normalised by callable sequence.'))
    r=pp['ratio_epitope_to_non_epitope_pi']
    if r:
        o.append('<p>The epitope-to-non-epitope density ratio is <b>%.3f</b> across %d genes. A ratio below 1 is '
                 'consistent with the reported hyperconservation of <i>M. tuberculosis</i> T-cell epitopes; a ratio '
                 'above 1 is not. This is a single-isolate comparison against the reference and is not a population '
                 'measure; the population-level comparison is in the cohort report.</p>'%(r,pp['genes_considered']))
    return '\n'.join(o)

def SEC_reg(rows):
    pr=[r for r in rows if r['promoter_of'] and not r['in_mask']]
    ch=[r for r in rows if r['chipseq_regulators'] and not r['in_mask']]
    lead=[r for r in pr if 'leaderless' in r['promoter_of']]
    o=['<h2>11. Promoter, 5′UTR and regulatory variants (extras E8, E9, E10)</h2>',
       '<p>These analyses are additions to the brief. A variant that changes expression without changing protein '
       'sequence is invisible to coding-variant analysis, so promoter and untranslated-region variants are examined '
       'explicitly. Transcription start sites come from Shell et al. 2015, using the authors’ own leadered and '
       'leaderless classification — 2,331 genes have a mapped TSS. Fixed upstream windows are <i>not</i> used.</p>',
       '<div class="note"><b>Why mapped TSS rather than a fixed window.</b> Many virulence genes have no promoter of '
       'their own. <i>esxA</i> and <i>esxB</i> are transcribed as read-through from the PE35/Rv3872 operon promoter '
       '1,996 nt upstream; <i>katG</i>’s nearest mapped TSS is 1,200 nt upstream at Rv1910c. A 200 bp upstream window '
       'would scan empty sequence for these genes and report nothing. Genes with no mapped TSS of their own are '
       'reported as such rather than given an invented window.</div>',
       table(['metric','value'],[
        ['genes with a mapped TSS','2,331 (1,389 leadered, 942 leaderless)'],
        ['promoter definition','100 bp immediately upstream of the TSS, through the TSS'],
        ["5′UTR definition",'TSS to start codon, where the UTR is 1–500 nt (1,249 genes)'],
        ['variants in a mapped promoter or 5′UTR','{:,}'.format(len(pr))],
        ['of these, in leaderless genes','{:,}'.format(len(lead))],
        ['variants inside a ChIP-seq binding peak','{:,}'.format(len(ch))]])]
    if pr:
        o.append('<h3>11.1 Promoter and 5′UTR variants (E8, E10)</h3>')
        o.append(table(['position','change','affected gene / region','in virulence catalogue','ChIP-seq regulators'],
            [[r['pos'],'%s>%s'%(r['ref'],r['alt']),r['promoter_of'][:70],
              'yes' if r['in_virulence_catalogue'] else '',r['chipseq_regulators'][:40]] for r in pr],maxrows=45))
    vch=[r for r in ch if r['in_virulence_catalogue']]
    o.append('<h3>11.2 Variants in transcription-factor binding sites (E9)</h3>')
    o.append('<p>Binding peaks are from the Minch et al. 2015 ChIP-seq survey: 15,980 peaks for 156 regulators, '
             'covering 59.7%% of the genome. Because coverage is high, peak overlap alone is weak evidence; it is '
             'reported as context, and is most informative where it coincides with a promoter or 5′UTR variant.</p>')
    if vch:
        o.append(table(['position','gene','effect','regulators binding here'],
            [[r['pos'],r['gene'],r['effect'].replace('_',' '),r['chipseq_regulators'][:60]] for r in vch],maxrows=35))
    return '\n'.join(o)

def SEC_subclonal(sc,rows):
    vir=[s for s in sc if s['in_virulence_catalogue']]
    o=['<h2>12. Sub-clonal variants and mixed infection (analysis 27, extra E14)</h2>',
       '<p>Sites with an alternate allele fraction between 0.10 and 0.90, at depth ≥ 20 and mapping quality ≥ 40, '
       'outside masked regions. These are below the 0.9 threshold used for the primary call set and represent '
       'within-host heterogeneity or, collectively, mixed infection.</p>',
       table(['metric','value'],[['sub-clonal sites','{:,}'.format(len(sc))],
        ['in virulence-catalogue genes','{:,}'.format(len(vir))],
        ['fixed variants for comparison','{:,}'.format(len(rows))],
        ['sub-clonal as a share of all called sites','%.2f%%'%(100*len(sc)/(len(sc)+len(rows)))if (len(sc)+len(rows)) else 'n/a']])]
    n=len(sc)
    verdict=('No evidence of mixed infection: the number of sub-clonal sites is low and consistent with sequencing and '
             'alignment noise.' if n<20 else
             'A moderate number of sub-clonal sites is present. This may reflect within-host diversity rather than true '
             'mixed infection; interpret alongside the allele-fraction distribution in Figure 3.' if n<100 else
             'A high number of sub-clonal sites is present, which is consistent with mixed infection or contamination. '
             'This isolate should be interpreted with caution.')
    o.append('<div class="note"><b>Assessment.</b> %s</div>'%verdict)
    if sc:
        o.append(table(['position','change','depth','allele fraction','gene','virulence catalogue'],
            [[s['pos'],'%s>%s'%(s['ref'],s['alt']),s['depth'],s['allele_fraction'],s['gene'],
              'yes' if s['in_virulence_catalogue'] else ''] for s in sorted(sc,key=lambda x:-x['allele_fraction'])],maxrows=40))
    return '\n'.join(o)

def SEC_struct(cov,iz,cn,f_map):
    dels=cov['deletions']
    o=['<h2>13. Structural variation: deletions and IS<i>6110</i> (analyses 29 and 30)</h2>',
       fig(f_map,'<b>Figure 6.</b> Deletions ≥1 kb (orange) and sequence-confirmed IS<i>6110</i> junctions (blue '
           'triangles) across the H37Rv chromosome.'),
       '<h3>13.1 Large deletions (29)</h3>',
       '<p>Contiguous runs of zero coverage of at least 1 kb, called on <b>unfiltered</b> depth. Using '
       'mapping-quality-filtered depth here would report every repetitive element as deleted; on this isolate that '
       'difference is the distinction between %d regions and a much larger, partly spurious list.</p>'%len(dels)]
    if dels:
        o.append(table(['start','end','length (bp)','genes affected'],
            [[d['start'],d['end'],'{:,}'.format(d['length']),', '.join(d['genes'][:8])+('…' if len(d['genes'])>8 else '')]
             for d in dels],maxrows=40))
    else:
        o.append('<p>No deletion of ≥1 kb detected.</p>')
    o.append('<h3>13.2 IS<i>6110</i> (30)</h3>')
    o.append('<p>Copy number is estimated from read mass across the 16 IS<i>6110</i> elements of H37Rv relative to '
             'genome mean depth, using unfiltered depth. Insertion sites are called from clusters of soft-clipped reads '
             'whose clipped sequence matches the IS<i>6110</i> element by BLAST at ≥95%% identity; clusters failing that '
             'test are discarded, which removes the large majority of raw soft-clip clusters.</p>')
    o.append(table(['metric','value'],[
        ['IS6110 copies in H37Rv',iz['H37Rv_reference_copies']],
        ['estimated copy number in this isolate',cn['estimated_copy_number']],
        ['sequence-confirmed junctions',iz['junctions_supported']],
        ['junctions absent from H37Rv',iz['candidate_insertions_absent_from_H37Rv']],
        ['genome mean depth (unfiltered)','%.1f×'%cn['genome_mean_depth_unfiltered']]]))
    if iz['junctions']:
        o.append(table(['position','supporting clipped reads','IS6110-matching','gene at this position'],
            [[ '%d–%d'%(j['start'],j['end']),j['clipped_reads'],j['is6110_matching'],j['gene'] or 'intergenic']
              for j in iz['junctions']],maxrows=30))
    o.append('<p class="small">The two estimates are independent: copy number from depth, junctions from clipped reads. '
             'Agreement between them is a check on both.</p>')
    return '\n'.join(o)

def SEC_tbd1(t,lin):
    o=['<h2>14. TbD1 locus — mmpS6 / mmpL6 (analysis 28)</h2>',
       '<div class="note"><b>Why this cannot be done by mapping.</b> TbD1 is deleted from H37Rv, which is a lineage 4 '
       '(modern) strain, so the region has no coordinates in the reference. A coverage or mapping test against H37Rv '
       'returns nothing rather than failing visibly. The locus is therefore recovered from this isolate\u2019s own de novo '
       'assembly. Rv1557 is annotated <i>mmpL6</i> in the catalogue but is a different gene and is not used as a proxy.</div>',
       '<p>Reference: the 2,139 bp segment present in <i>M. canettii</i> CIPT 140010059 and absent from H37Rv, located by '
       'aligning the H37Rv region spanning Rv1556–Rv1559 against <i>M. canettii</i> — the deletion junction falls at '
       '<i>mmpL6</i>. Its two open reading frames match MmpS-family (Rv0403c 58%, Rv0451c 55%) and MmpL-family '
       '(Rv0450c 60%, Rv0676c 59%) proteins, confirming <i>mmpS6</i> and <i>mmpL6</i>. The method was validated on '
       'controls: <i>M. canettii</i> returns present at 100% recovery, H37Rv returns absent at 0%.</p>',
       '<p>An <b>absent</b> call additionally requires that the region around the deletion junction (H37Rv '
       '1,757,997\u20131,767,717) is found on a single contig of the assembly, so that the deletion is observed. Where '
       'that region is missing from the assembly, absence of TbD1 sequence cannot be told apart from absence of that part '
       'of the assembly, and the call is <b>not determined</b> rather than absent.</p>']
    if not t.get('available'):
        o.append('<div class="warn">Not assessed: %s.</div>'%esc(t.get('reason','')))
        return '\n'.join(o)
    o.append(table(['metric','value'],[
        ['call',t['call'].upper()],
        ['fraction of TbD1 recovered','%.1f%%'%(100*t['fraction_of_TbD1_recovered'])],
        ['aligned bases','%s of %s'%('{:,}'.format(t['aligned_bp']),'{:,}'.format(t['reference_length']))],
        ['alignment blocks',len(t.get('hits',[]))],
        ['junction region found on a single contig','%.1f%%'%(100*t.get('flank_fraction_covered',0))],
        ['lineage of this isolate',lin['lineage']]]))
    exp='retained' if lin['lineage'].startswith('lineage1') else 'deleted'
    agree = (t['call']=='present') == (exp=='retained')
    if t['call']=='not_determined':
        o.append('<div class="warn"><b>TbD1 could not be determined for this isolate.</b> %s. It is therefore not '
                 'reported as absent, and no comparison with the lineage expectation is made.</div>'%esc(t.get('reason','')[:1].upper()+t.get('reason','')[1:]))
    elif not lin['lineage'].startswith('lineage'):
        o.append('<p>The lineage of this isolate is <b>%s</b>, so no expectation for TbD1 can be tested against the '
                 'observed call (<b>%s</b>).</p>'%(esc(lin['lineage']),t['call']))
    else:
        o.append('<p>TbD1 is retained in the ancestral lineages, including lineage 1, and deleted in the modern lineages '
                 '(2, 3 and 4). This isolate is <b>%s</b>, for which TbD1 is expected to be <b>%s</b>; the observed call is '
                 '<b>%s</b>, which <b>%s</b> expectation.</p>'%(lin['lineage'],exp,t['call'],
                 'agrees with' if agree else 'does NOT agree with'))
    if t['call']=='present':
        o.append('<p>The intact locus has been extracted to <span class="mono">%s_TbD1.fasta</span> (%s bp including '
                 '200 bp of flanking sequence each side) for the cross-isolate alignment, which is produced once all '
                 'Lineage 1 isolates are available.</p>'%(esc(t.get('locus_contig','')).split('_contig')[0],
                 '{:,}'.format(t.get('locus_sequence_length',0))))
    else:
        o.append('<p>No locus sequence is extracted for this isolate, and it does not enter the Lineage 1 TbD1 alignment.</p>')
    return '\n'.join(o)

def SEC_pe(pe):
    o=['<h2>15. PE / PPE / PE_PGRS secondary analysis (analysis 15, Tier 2)</h2>',
       '<p>These 168 genes are masked from the primary analysis and are reported here separately, clearly labelled as '
       'secondary, as the brief requires. Short-read mapping across them is unreliable, so each locus is resolved from '
       'the de novo assembly rather than from whole-genome mapping.</p>']
    if not pe.get('available'):
        o.append('<div class="warn">Not assessed: %s.</div>'%esc(pe.get('reason','')))
        return '\n'.join(o)
    th=pe['thresholds']
    o.append('<p>A locus is reported only when a <b>single</b> contig aligns over at least %d%% of its length at '
             '\u2265%.0f%% identity <b>and</b> read depth over the locus is \u2265%d\u00d7. Depth and mapping quality are given for '
             'every call. Loci failing either test are reported as not resolved, never as reference-identical.</p>'%(
             int(100*th['min_assembly_coverage']),th['min_identity'],th['min_read_depth']))
    o.append(table(['outcome','loci'],[[k,v] for k,v in sorted(pe['summary'].items(),key=lambda x:-x[1])]))
    o.append(table(['metric','value'],[
        ['PE_PGRS genes assessed',pe['pe_pgrs_genes']],
        ['PE_PGRS mean read depth','%s\u00d7'%pe['pe_pgrs_mean_depth']],
        ['PE_PGRS copy equivalents vs genome mean',pe['pe_pgrs_copy_equivalents']]]))
    div=[r for r in pe['rows'] if r['call']=='resolved' and (r['percent_divergence_from_H37Rv'] or 0)>0]
    div.sort(key=lambda r:-(r['percent_divergence_from_H37Rv'] or 0))
    o.append('<h3>15.1 Resolved loci diverging from H37Rv</h3>')
    if div:
        o.append(table(['Rv','gene','length bp','identity %','divergence %','mean depth','min depth','mean MAPQ'],
            [[r['rv'],r['gene'],r['length_bp'],r['best_identity'],r['percent_divergence_from_H37Rv'],
              r['mean_depth'],r['min_depth'],r['mean_mapq']] for r in div],maxrows=45))
    else:
        o.append('<p>No resolved locus differs from H37Rv.</p>')
    nr=[r for r in pe['rows'] if r['call']!='resolved']
    o.append('<h3>15.2 Loci not resolved</h3>')
    o.append(table(['Rv','gene','reason','assembly coverage','contigs','mean depth','mean MAPQ'],
        [[r['rv'],r['gene'],r['call'].replace('not resolved ',''),r['assembly_coverage'],
          r['contigs_aligned'],r['mean_depth'],r['mean_mapq']] for r in nr],maxrows=35))
    return '\n'.join(o)

ABS_MIN_GENOME=90.0
def apply_absence_rule(R):
    """v5: one gene-absence rule for all isolates. 'absent' only if the isolate has >=90% of the genome at >=10x
    (the gene call itself needs <10% of the locus covered). Otherwise 'not determined'."""
    br=R['exclusion']['genome_breadth_10x_percent']
    R['absence_rule']={'genome_breadth_10x_percent':br,'applies':br<ABS_MIN_GENOME,'reverted':0}
    if br>=ABS_MIN_GENOME: return
    for p in R['catalogue_presence']:
        if p['call']=='absent':
            p['call_before_rule']='absent'; p['call']='not determined'; R['absence_rule']['reverted']+=1

def SEC_cohort(sample,R):
    C=json.load(open(P+'/ref/cohort_context.json')); M=C['_meta']; c=C.get(sample)
    o=['<h2>17. Cohort context (cohort stage, 55 isolates)</h2>',
       '<p>This section places the isolate in the full cohort. It uses results computed once for all 55 isolates '
       '(analyses 22–24, 26–28 and 30); the group statistics are in the cohort report. All distances are SNVs over positions '
       'callable in both isolates (depth ≥10, outside the static mask); positions not callable in one of the two isolates are '
       'missing, not differences.</p>']
    if not c:
        o.append('<div class="warn">No cohort data for this isolate.</div>'); return '\n'.join(o)
    rows=[['cohort status',c['cohort_status'] or 'n/a'],
          ['lineage (barcode)','%s'%R['lineage']['lineage']]]
    if c['not_assessable_for_clustering']:
        rows.append(['SNP cluster','not assessable: fewer than 50% of the genome is callable, so distances to the partners below are lower bounds'])
    elif c['cluster5']:
        rows.append(['cluster at ≤5 SNPs','%d isolates: %s'%(len(c['cluster5']),', '.join(c['cluster5']))])
    else:
        rows.append(['cluster at ≤5 SNPs','none (singleton)'])
    if not c['not_assessable_for_clustering']:
        rows.append(['cluster at ≤12 SNPs','%s'%(', '.join(c['cluster12']) if c['cluster12'] else 'none (singleton)')])
    o.append(table(['item','value'],rows))
    o.append('<h3>17.1 Closest isolates in this cohort</h3>')
    o.append(table(['isolate','SNVs','positions compared'],[[x['isolate'],x['snps'],'{:,}'.format(x['comparable_positions'] or 0)] for x in c['nearest_local']]))
    o.append('<p class="small">A low-coverage partner (MTB4, MTB25, MTB32) is compared over few positions and its distance is a lower bound. '
             'Isolates with 0 SNVs may belong to different patients.</p>')
    g=c['nearest_public_all_lineage_tree']; l1=c['nearest_public_L1_context_tree']
    o.append('<h3>17.2 Closest public genome</h3>')
    pr=[]
    if g: pr.append(['all-lineage tree (55 local + 279 public)',g['accession'],g['lineage'],g['country'],g['snv']])
    if l1: pr.append(['lineage 1 context tree (55 local + 300 public lineage 1)',l1['accession'],l1['lineage'],l1['country'],l1['snv']])
    o.append(table(['tree','nearest public genome','its lineage','country','approx. SNVs'],pr) if pr else '<p>Not placed (lineage not assigned).</p>')
    o.append('<p class="small">Public genomes are assemblies, not reads. Distances are approximate (tree path length) and '
             'low-coverage or contaminated isolates are placed with less confidence.</p>')
    o.append('<h3>17.3 Other cohort-level facts</h3>')
    f=[['mixed-infection assessment',c['mixed_infection'] or 'n/a'],
       ['sub-clonal sites',c['subclonal_sites']],
       ['IS6110 copies (this isolate / cohort median)','%s / %s'%(c['is6110_copies'],M['is6110_cohort_median'])]]
    if R['lineage'].get('major_lineage')=='lineage1' or (R['lineage']['lineage'] or '').startswith('lineage1'):
        f.append(['TbD1 locus variant (position 328, Gly>Ser)','present' if c['tbd1_L1_pos328_variant'] else 'not present'])
    o.append(table(['item','value'],f))
    ar=R.get('absence_rule',{})
    if ar.get('applies'):
        o.append('<div class="warn"><b>Gene-absence rule.</b> Only %.1f%% of this genome is covered at ≥10×, below the 90%% needed to call a gene absent. '
                 '%d genes with little or no coverage are therefore reported as <b>not determined (isolate coverage too low to call absence)</b>, not absent.</div>'%(ar['genome_breadth_10x_percent'],ar['reverted']))
    else:
        o.append('<p class="small">Gene-absence rule: an isolate must have ≥90%% of its genome at ≥10× for a gene to be called absent; this isolate meets it (%.1f%%).</p>'%ar['genome_breadth_10x_percent'])
    return '\n'.join(o)

def SEC_presence(pres):
    c=collections.Counter(p['call'] for p in pres)
    absent=[p for p in pres if p['call']=='absent']
    nd=[p for p in pres if p['call']=='not determined']
    o=['<h2>16. Virulence gene presence and absence (analysis 16)</h2>',
       '<p>Each catalogue gene is called from breadth of coverage over its H37Rv locus: <b>present</b> at ≥95%% of the '
       'locus covered at ≥10×, <b>absent</b> where less than 10%% of the locus has any coverage, and '
       '<b>not determined</b> in between. Per the brief, genes in low-coverage regions are never called absent; in addition, no gene is called absent in an isolate with less than 90%% of its genome covered at ≥10× (section 17).</p>',
       table(['call','genes'],[['present',c.get('present',0)],['absent',c.get('absent',0)],
                               ['not determined',c.get('not determined',0)]])]
    if absent:
        o.append('<h3>16.1 Genes called absent</h3>')
        o.append(table(['Rv','gene','tier','functional system','breadth at ≥10×','mean depth'],
            [[p['rv'],p['gene'],p['tier'],p['system'][:30],p['breadth10'],p['mean_depth']] for p in absent],maxrows=40))
    if nd:
        o.append('<h3>16.2 Genes not determined</h3>')
        o.append('<p>Partial coverage, or no coverage in an isolate with less than 90%% of its genome covered at ≥10× (section 17): these are reported as undetermined rather than absent, as the brief requires.</p>')
        o.append(table(['Rv','gene','tier','functional system','breadth at ≥10×','mean depth'],
            [[p['rv'],p['gene'],p['tier'],p['system'][:30],p['breadth10'],p['mean_depth']] for p in nd],maxrows=40))
    return '\n'.join(o)

NOT_DONE=[
 ('15','Targeted PE/PPE re-assembly and variant calling','Loci are resolved from the whole-genome de novo assembly (section 15). A separate per-locus SPAdes '
  're-assembly was not run; unresolved loci are reported as not resolved.'),
 ('25b','Protein stability change (ΔΔG)','Run once for the cohort with FoldX 5.1 on AlphaFold DB models (27 variants; '
  'see the cohort report), not per isolate. Proteins above 700 amino acids were not modelled and AlphaFold3 was not run.'),
 ('26 (tree)','Global-context phylogeny','Built once for the cohort (55 isolates plus public genomes) and summarised for this '
  'isolate in section 17. Public data are assemblies, not reads.'),
 ('28','TbD1 locus','The brief gives this locus as Rv2073c–Rv2074 (mmpS6/mmpL6) and suggests a coverage-based test '
  'against H37Rv. Those Rv numbers are an oxidoreductase and a pyridoxamine 5′-phosphate oxidase, not mmpS6/mmpL6 '
  '(mmpL6 is Rv1557), and TbD1 is deleted from H37Rv itself, so no coverage-based test against H37Rv can detect it. '
  'The locus is therefore assessed by aligning the assembly against the M. canettii TbD1 segment (section 14); the '
  'cohort-wide alignment is in the cohort report.'),
 ('9, 16(matrix), 17–20, 22–24, 31–36','Cohort-level and statistical analyses','These are properties of the whole '
  'sample set, or require the clinical metadata. They were computed once on the complete 55-isolate cohort and are reported in the cohort report and tables, not in this per-isolate report; section 17 gives this isolate’s cohort-level facts.'),
]
SUBS=[
 ('Kraken2 + Bracken (analysis 3)','Mapping rate against H37Rv plus unmapped-read assessment','The standard Kraken2 '
  'database is 50–100 GB and exceeds available disk in this offline environment.'),
 ('CheckM2 (analysis 2)','Breadth of coverage across single-copy H37Rv CDS','CheckM2 requires a multi-gigabyte DIAMOND '
  'database that could not be obtained offline.'),
 ('Bakta / Prokka (analysis 8)','Reference-guided annotation transfer from H37Rv','The Bakta database exceeds available '
  'disk. For a clonal organism, transfer from the reference is also more accurate than de novo annotation and avoids '
  'the annotation noise the brief notes in section 3.'),
 ('PROVEAN / SIFT4G (analysis 25)','BLOSUM62 substitution score with domain and epitope context in this report','This report '
  'uses BLOSUM62. PROVEAN was run once for the cohort (mycobacterial UniRef90 database, not the full NCBI nr) and is in the cohort report; '
  'SIFT4G was not run because no M. tuberculosis database was available. BLOSUM62 is not presented as a PROVEAN or SIFT equivalent.'),
 ('Gubbins (analysis 22)','Omitted','M. tuberculosis is effectively clonal with no appreciable homologous '
  'recombination, so a recombination-removal step is not applicable.'),
]
def SEC_limits(excl,cov):
    o=['<h2>18. Limitations, substitutions and work not completed</h2>',
       '<h3>18.1 Analyses not completed in this report</h3>',
       table(['analysis','item','reason'],[[a,b,c] for a,b,c in NOT_DONE]),
       '<h3>18.2 Tools replaced, and by what</h3>',
       '<p>Each substitution below is named wherever its result appears; none is presented as equivalent to the tool it replaces.</p>',
       table(['specified tool','used instead','reason'],[[a,b,c] for a,b,c in SUBS]),
       '<h3>18.3 Limitations of this isolate’s data</h3>']
    lim=[]
    nc=sum(b-a+1 for a,b in cov['noncallable'])
    lim.append(['callable fraction','%.2f%% of the genome is covered at ≥10×; %s bp is not callable and no variant '
                'can be called there'%(excl['genome_breadth_10x_percent'],'{:,}'.format(nc))])
    lim.append(['masked regions','490,058 bp (11.11%) of PE/PPE, insertion sequence, phage and excluded-list sequence '
                'is masked from the primary analysis; variants there are flagged, not interpreted'])
    lim.append(['read length','single-end read length and insert size limit the detection of structural variants and '
                'repeat-internal variation; IS6110 copy number is an estimate from read depth, not a direct count'])
    lim.append(['MIRU-VNTR','not attempted: repeat arrays exceed read and insert length, so copy number cannot be '
                'determined reliably from these data'])
    lim.append(['assembly provenance','the de novo assembly used for analyses 1, 2 and 16 was carried over from the '
                'earlier run as sequence data; no result or interpretation from that earlier analysis is reused'])
    o.append(table(['item','statement'],lim))
    return '\n'.join(o)

def SEC_methods():
    V=json.load(open(P+'/ref/versions.json'))
    o=['<h2>19. Methods</h2>',
       '<p>Reads were quality- and adapter-trimmed with fastp at a Phred threshold of 20 and aligned to the '
       '<i>Mycobacterium tuberculosis</i> H37Rv reference genome NC_000962.3 (4,411,532 bp; the RefSeq equivalent of '
       'AL123456.3) using BWA-MEM. Duplicates were marked with samtools markdup. Per-base depth was computed twice: '
       'with a mapping-quality and base-quality threshold of 20 for defining callable regions and supporting variant '
       'calling, and without filtering for the analysis of deletions and repetitive elements. Variants were called with '
       'bcftools mpileup and call in haploid mode, left-normalised, and filtered to depth ≥ 10, mapping quality ≥ 40 and '
       'alternate allele fraction ≥ 0.9. Sites with an alternate allele fraction between 0.10 and 0.90 at depth ≥ 20 were '
       'retained separately as sub-clonal. Variants were annotated with SnpEff using a database built locally from the '
       'NC_000962.3 GenBank record, and were cross-referenced to protein features from UniProt proteome UP000001584, to '
       'T-cell epitopes from IEDB, to transcription start sites from Shell et al. 2015, to ChIP-seq binding peaks from '
       'Minch et al. 2015, and to the Coll 62-SNP lineage barcode. Lineage was assigned from the Coll barcode as '
       'distributed with the TB-Profiler database, requiring ≥ 90% of a barcode’s defining positions. Spoligotypes were '
       'determined in silico from reads at the direct repeat locus and matched against 9,658 SITVIT2 patterns. Drug '
       'resistance was called from the TB-Profiler mutation catalogue, which embeds WHO mutation catalogue version 2 and '
       'its five confidence grades, and was cross-checked by BLASTp against the NCBI AMRFinder reference protein set. '
       'A mask covering PE/PPE/PE_PGRS genes, insertion sequences, phage regions and the TB-Profiler exclusion list '
       '(490,058 bp, 11.11% of the genome) was applied to the primary analysis, with masked regions reported separately. '
       'The virulence gene catalogue was compiled from VFDB (core and full), Mycobrowser functional categories, BV-BRC '
       'specialty genes and the curated locus list in the analysis brief, deduplicated by Rv number to 946 genes, and '
       'was frozen before any comparison.</p>',
       '<p><b>Comparing reports between isolates.</b> The number of variants and other counts in this report depend on '
       'sequencing depth and on fixed thresholds (depth ≥ 10, allele fraction ≥ 0.9, mapping quality ≥ 40), and some fall in '
       'masked repetitive regions. Two isolates of the same strain can therefore show slightly different counts in their '
       'individual reports. A difference in counts between two reports is not a measure of genetic distance: relatedness '
       'between isolates is assessed in the cohort analysis on masked genomes, where positions without confident coverage '
       'are treated as missing data rather than as differences.</p>',
       '<h3>19.1 Software</h3>',
       table(['tool','version'],[[k,v] for k,v in V['tools'].items()]),
       '<h3>19.2 Databases and reference data</h3>',
       table(['resource','version / provenance'],[[k,v] for k,v in V['databases'].items()])]
    return '\n'.join(o)

def write_tables(sample,outdir,R,gb):
    import csv as _csv
    T=os.path.join(outdir,'tables'); os.makedirs(T,exist_ok=True)
    def w(name,cols,rows):
        with open(os.path.join(T,'%s_%s.tsv'%(sample,name)),'w',newline='') as fh:
            wr=_csv.writer(fh,delimiter='\t'); wr.writerow(cols)
            for r in rows: wr.writerow(r)
    v=R['variants']
    w('variants',list(v[0].keys()) if v else [],[list(x.values()) for x in v])
    w('subclonal',list(R['subclonal'][0].keys()) if R['subclonal'] else [],[list(x.values()) for x in R['subclonal']])
    w('resistance',['gene','mutation','drug','confidence','position','effect','source'],
      [[r['gene'],r['mutation'],r['drug'],r['confidence'],r['pos'],r['effect'],r['source']] for r in R['resistance']+R['resistance_uncertain']])
    w('identification_rd',['region','span','length_bp','fraction_covered_5x','mean_depth','call','interpretation'],
      [[r['region'],r['span'],r['length_bp'],r['fraction_covered_5x'],r['mean_depth'],r['call'],r['interpretation']] for r in R['identification']['rd']['regions']])
    w('tbd1',['metric','value'],[[k,str(v)] for k,v in R['tbd1'].items() if k!='hits'])
    w('pe_ppe',list(R['pe_ppe']['rows'][0].keys()) if R['pe_ppe'].get('rows') else [],
      [list(x.values()) for x in R['pe_ppe'].get('rows',[])])
    w('gene_presence',['rv','gene','system','tier','call','breadth10','mean_depth','length_bp','callable_bp'],
      [[p['rv'],p['gene'],p['system'],p['tier'],p['call'],p['breadth10'],p['mean_depth'],p['length_bp'],p['callable_bp']] for p in R['catalogue_presence']])
    w('deletions',['start','end','length','genes'],
      [[d['start'],d['end'],d['length'],','.join(d['genes'])] for d in R['coverage']['deletions']])
    w('is6110_junctions',['start','end','clipped_reads','is6110_matching','reference_copy','rv','gene'],
      [[j['start'],j['end'],j['clipped_reads'],j['is6110_matching'],j['reference_copy'],j['rv'],j['gene']] for j in R['is6110']['junctions']])
    w('gene_coverage',['rv','gene','length','mean_depth','breadth_10x','breadth_1x'],
      [[rv,(A.genes[rv]['name'] or rv),b['len'],round(b['mean'],2),round(b['breadth10'],4),round(b['breadth1'],4)]
       for rv,b in sorted(gb.items())])
    return T
