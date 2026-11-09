import json,collections,csv,sys
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c4'
def readfa(f):
    d={};cur=None
    for l in open(f):
        if l[0]=='>': cur=l[1:].strip(); d[cur]=[]
        else: d[cur].append(l.strip())
    return {k:''.join(v).upper() for k,v in d.items()}
b='TCAG'; aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; code={};i=0
for x in b:
    for y in b:
        for z in b: code[x+y+z]=aa[i];i+=1
comp=str.maketrans('ACGT','TGCA')
ref=''.join(l.strip() for l in open(P+'/ref/TbD1.fa') if l[0]!='>').upper()
# ORFs >= 300 nt (ATG..stop) both strands
orfs=[]
for strand in '+-':
    s_=ref if strand=='+' else ref.translate(comp)[::-1]
    for f in range(3):
        start=f
        for i in range(f,len(s_)-2,3):
            if code[s_[i:i+3]]=='*' or i+3+3>len(s_)-2:
                end=i if code[s_[i:i+3]]=='*' else i+3
                if end-start>=300:
                    st,en=(start+1,end) if strand=='+' else (len(s_)-end+1,len(s_)-start)
                    orfs.append((strand,st,en,(end-start)//3))
                start=i+3
# keep longest non-nested
orfs=sorted(orfs,key=lambda x:-(x[2]-x[1]))
keep=[]
for o in orfs:
    if all(min(o[2],k[2])-max(o[1],k[1])<0.5*(o[2]-o[1]) for k in keep): keep.append(o)
keep=[('+',1,399,133),('+',396,2138,581)]  # Prodigal meta-mode CDS on the TbD1 reference: ORF1 = mmpS6 (5'-truncated at segment edge), ORF2 = mmpL6 5' portion (3'-truncated at segment edge); names inferred from the earlier protein-family match
keep.sort(key=lambda x:x[1]); print('ORFs in TbD1 reference (strand,start,end,aa):',keep)
summary={'orfs':keep}
for tag in ('L1_local','L1_local_plus_public'):
    A=readfa('%s/tbd1_%s_aligned.fa'%(O,tag)); names=list(A); L=len(A[names[0]]); r=A['canettii_ref']
    assert all(len(v)==L for v in A.values())
    # reference coordinate per column
    rc=[];p=0
    for ch in r:
        if ch!='-': p+=1
        rc.append(p if ch!='-' else None)
    rows=[]
    for c in range(L):
        col={n:A[n][c] for n in names if n!='canettii_ref'}
        bases=collections.Counter(v for v in col.values() if v not in 'N')
        if len(bases)<=1 and (r[c] in list(bases)+['-'] ) and not (len(bases)==1 and list(bases)[0]!=r[c]): 
            if not (len(bases)==1 and list(bases)[0] != r[c]): continue
        if len(bases)==0: continue
        alleles=dict(bases)
        isgap=('-' in alleles) or r[c]=='-'
        typ='indel' if isgap else 'SNV'
        diff_ref=[n for n,v in col.items() if v not in 'N' and v!=r[c]]
        rows.append({'align_col':c+1,'ref_pos':rc[c] if rc[c] else '','ref_base':r[c],'type':typ,'alleles':';'.join('%s:%d'%(k,v) for k,v in sorted(alleles.items())),'n_with_data':sum(alleles.values()),'n_differ_from_canettii':len(diff_ref),
            'minor_allele_isolates':','.join(n for n,v in col.items() if v not in 'N' and alleles and v!=bases.most_common(1)[0][0])[:300],'polymorphic_within_L1':len(alleles)>1})
    # annotate coding effect for SNVs vs canettii-ref codon
    for x in rows:
        x['orf']='';x['effect']=''
        if x['type']=='SNV' and x['ref_pos']:
            for (st,a,e,n) in keep:
                if a<=x['ref_pos']<=e:
                    x['orf']=('mmpS6 ' if a==1 else 'mmpL6 (5 prime portion) ')+'%s:%d-%d(%d aa)'%(st,a,e,n); 
                    # effect for the most common alt allele
                    alt=[k for k in x['alleles'].split(';')][0]
                    try:
                        if st=='+': off=(x['ref_pos']-a)%3; cs=x['ref_pos']-off; c0=ref[cs-1:cs+2]; mut=list(c0); mut[off]=[k.split(':')[0] for k in x['alleles'].split(';') if k.split(':')[0]!=x['ref_base']][0]
                        else:
                            seq=ref.translate(comp)[::-1]; pr=len(ref)-x['ref_pos']+1; s0=len(ref)-e+1; off=(pr-s0)%3; cs=pr-off; c0=seq[cs-1:cs+2]; mut=list(c0); mut[off]=[k.split(':')[0] for k in x['alleles'].split(';') if k.split(':')[0]!=comp_b][0] if False else mut[off]
                            alt_b=[k.split(':')[0] for k in x['alleles'].split(';') if k.split(':')[0]!=x['ref_base']][0]; mut[off]=alt_b.translate(comp)
                        m=''.join(mut); x['effect']='synonymous' if code.get(c0)==code.get(m) else 'missense (%s>%s)'%(code.get(c0),code.get(m))
                    except Exception as ex: x['effect']='NA'
    with open('%s/a28_tbd1_variants_%s.tsv'%(O,tag),'w') as o:
        k=['align_col','ref_pos','ref_base','type','alleles','n_with_data','n_differ_from_canettii','polymorphic_within_L1','minor_allele_isolates','orf','effect']
        w=csv.DictWriter(o,fieldnames=k,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    nl=len(names)-1
    summary[tag]={'n_isolates':nl,'alignment_length':L,'sites_differing_from_canettii':len(rows),'polymorphic_within_L1':sum(x['polymorphic_within_L1'] for x in rows),'indel_sites':sum(x['type']=='indel' for x in rows),
        'missing_N_fraction_per_isolate':{n:round(A[n].count('N')/L,3) for n in names if n!='canettii_ref' and A[n].count('N')>0}}
json.dump(summary,open(O+'/a28_tbd1_summary.json','w'),indent=1,default=str); print(json.dumps(summary,indent=1,default=str)[:1500])
