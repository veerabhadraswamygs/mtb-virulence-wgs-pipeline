"""Analysis 8: uniform GFF3 re-annotation of all 55 de novo assemblies.
Gene calls = the same Prodigal (-p single) calls already used for every isolate (coordinates are in the .faa headers);
names transferred from H37Rv by DIAMOND blastp (>=60% identity and >=80% of both query and subject length)."""
import json,subprocess,os,re,sys,collections
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
OUT=P+'/cohort/a8'
G=json.load(open(P+'/ref/genes.json')); PR=json.load(open(P+'/ref/h37rv_prot.json'))
ref=OUT+'/H37Rv_proteome.faa'
if not os.path.exists(ref):
    with open(ref,'w') as o:
        for k,s in PR.items(): o.write('>%s\n%s\n'%(k,s))
    subprocess.run('%s/bin/diamond makedb --in %s -d %s/H37Rv_prot -p 1 --quiet'%(P,ref,OUT),shell=True,check=True)
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
summ=[]
for s in samples:
    faa='%s/work/%s/%s_proteome_p2.faa'%(P,s,s)
    tsv='%s/%s.dmnd.tsv'%(OUT,s)
    if os.path.getsize(faa)==0: open(tsv,'w').close()
    elif not os.path.exists(tsv):
        subprocess.run("%s/bin/diamond blastp -q %s -d %s/H37Rv_prot -o %s -p 1 --quiet --max-target-seqs 5 -e 1e-10 "
          "--outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"%(P,faa,OUT,tsv),shell=True,check=True)
    best={}
    for l in open(tsv):
        f=l.split('\t'); q=f[0]; pid=float(f[2]); L=int(f[3]); ql=int(f[4]); sl=int(f[5]); bs=float(f[11])
        if q not in best or bs>best[q][0]: best[q]=(bs,f[1],pid,(int(f[7])-int(f[6])+1)/ql,(int(f[9])-int(f[8])+1)/sl)
    genes=[]; seq=None
    for l in open(faa):
        if l.startswith('>'):
            h=l[1:].split(' # ')
            genes.append((h[0].strip(),int(h[1]),int(h[2]),'+' if h[3].strip()=='1' else '-',h[4]))
    n_asg=0; found=set(); ncl=collections.Counter()
    with open('%s/%s.reannot.gff3'%(OUT,s),'w') as o:
        o.write('##gff-version 3\n##annotation: Prodigal -p single + DIAMOND transfer from H37Rv NC_000962.3 (analysis 8)\n')
        for gid,st,en,sd,attr in genes:
            contig=gid.rsplit('_',1)[0]
            b=best.get(gid)
            if b and b[2]>=60 and b[3]>=0.8 and b[4]>=0.8:
                rv=b[1]; g=G.get(rv,{}); n_asg+=1; found.add(rv)
                a='ID=%s;Name=%s;locus_tag=%s;product=%s;ref_identity=%.1f;ref_qcov=%.2f;ref_scov=%.2f;category=%s'%(gid,g.get('name') or rv,rv,g.get('product','').replace(';',','),b[2],b[3],b[4],g.get('category',''))
            else:
                a='ID=%s;product=hypothetical protein (no H37Rv match >=60%% id, >=80%% coverage)'%gid
            m=re.search(r'partial=(\d+)',attr)
            if m and m.group(1)!='00': a+=';partial=%s'%m.group(1)
            o.write('%s\tProdigal+DIAMOND\tCDS\t%d\t%d\t.\t%s\t0\t%s\n'%(contig,st,en,sd,a))
    af=json.load(open('%s/work/%s/asm_filter.json'%(P,s)))
    tl=af.get('total_bp',0)
    summ.append({'sample':s,'assembly_total_bp':tl,'assembly_contigs':af.get('kept',0),'assessable_for_gene_content':bool(4000000<=tl<=4700000 and 3500<=len(genes)<=4600 and af.get('kept',0)<=200),'cds_called':len(genes),'cds_assigned_to_H37Rv':n_asg,'cds_unassigned':len(genes)-n_asg,
                 'distinct_H37Rv_genes_found':len(found),'H37Rv_genes_not_found':len(PR)-len(found)})
    json.dump({'sample':s,'found':sorted(found)},open('%s/%s.found.json'%(OUT,s),'w'))
json.dump(summ,open(OUT+'/a8_summary.json','w'),indent=1)
print(len(summ),'samples done')
