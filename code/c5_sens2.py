"""C5 sensitivity: one isolate per <=5-SNP cluster (connected components of the pairwise-complete SNP distance matrix; MTB4/25/32/35 are unclusterable singletons)."""
import re,sys
src=open('/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/code/c5_stats.py').read().split('# ================================================================== run diabetes (31-33)')[0]
exec(src)
m=pd.read_csv(P+'/cohort/c1/a24_snp_distance_matrix.tsv',sep='\t',index_col=0)
ids=[s for s in S]; par={s:s for s in ids}
def f(x):
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
for i,a in enumerate(ids):
    for b in ids[i+1:]:
        if a in m.index and b in m.columns and a not in LOW and b not in LOW and m.loc[a,b]<=5: par[f(a)]=f(b)
cl=collections.defaultdict(list)
for s in ids: cl[f(s)].append(s)
rep=[sorted(v,key=lambda x:int(x[3:]))[0] for v in cl.values()]; rep=sorted([r for r in rep if r not in LOW],key=lambda x:int(x[3:]))
rows=[]
for k,v in cl.items():
    v=sorted(v,key=lambda x:int(x[3:])); rows.append(dict(cluster_rep=v[0],n=len(v),members=','.join(v),lineage=LIN[v[0]],diabetes_status=','.join(sorted(set(C[s]['diabetes'] for s in v))),discordant_for_diabetes=len(set(C[s]['diabetes'] for s in v))>1))
cdf=pd.DataFrame(rows).sort_values('cluster_rep',key=lambda x:x.str[3:].astype(int)); pass
print('clusters',len(cl),'representatives',len(rep),'multi-member',int((cdf.n>1).sum()),'discordant',int(cdf.discordant_for_diabetes.sum()))
lab='strict: one isolate per <=5-SNP cluster, excluding 4 low-coverage/contaminated isolates (n=%d)'%len(rep)
b=burden_tests('diabetes','Yes','Diabetes mellitus',sel=rep,sens=lab)
pd.DataFrame(b).to_csv(O+'/a31_a33_burden_diabetes_strict_sensitivity.tsv',sep='\t',index=False)
A=assoc_tables('diabetes','Yes','Diabetes mellitus',sel=rep,sens=lab); A.to_csv(O+'/a32_a33_variant_gene_tests_diabetes_strict_sensitivity.tsv',sep='\t',index=False)
print('diabetes done',len(A)); sys.stdout.flush()
CAND2=[('sex','Male','Sex (male vs female)'),('smoking','Yes','Smoking (ADDITIVE HABITS includes Smoker)'),('ses','BPL','Socioeconomic status (BPL vs APL)'),('smear_3plus','Yes','Sputum smear 3+ vs lower'),('died','Yes','Clinical outcome: dead'),
 ('cough','Yes','Persistent cough'),('fever','Yes','Fever'),('weight_loss','Yes','Weight loss'),('chest_pain','Yes','Chest pain'),('breathlessness','Yes','Breathlessness')]
bb=[];vv=[]
for var,pos,label in CAND2:
    n1=sum(C[s][var]==pos for s in rep)
    if min(n1,len(rep)-n1)<5:
        print('skip',label,n1,len(rep)-n1); continue
    bb+=burden_tests(var,pos,label,sel=rep,sens=lab); vv.append(assoc_tables(var,pos,label,sel=rep,sens=lab)); print(label); sys.stdout.flush()
pd.DataFrame(bb).to_csv(O+'/a35_burden_other_clinical_strict_sensitivity.tsv',sep='\t',index=False)
pd.concat(vv).to_csv(O+'/a35_variant_gene_tests_other_clinical_strict_sensitivity.tsv',sep='\t',index=False)
print('done')

# burden only, excluding the 4 low-quality isolates (no cluster collapse), other clinical variables
sel=[s for s in S if s not in LOW]; bb=[]
for var,pos,label in CAND2+[]:
    n1=sum(C[s][var]==pos for s in sel)
    if min(n1,len(sel)-n1)>=5: bb+=burden_tests(var,pos,label,sel=sel,sens='excluding 4 low-coverage/contaminated isolates (n=51)')
pd.DataFrame(bb).to_csv(O+'/a35_burden_other_clinical_excl_lowcov.tsv',sep='\t',index=False); print('done2')
