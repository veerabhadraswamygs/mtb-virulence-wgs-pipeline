"""C6: FoldX on AlphaFold Server toxin-antitoxin complexes (vapBC36 Rv1982c+Rv1982A, vapBC6 Rv0656c+Rv0657c): interface energy (AnalyseComplex),
ddG of the cohort variant in the complex (BuildModel x3), and interface-distance of the mutated residue. Model 0 (top-ranked) of each AlphaFold Server job."""
import os,re,json,glob,shutil,subprocess,math,csv,sys
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; FX=glob.glob(P+'/src/foldx/foldx_*')[0]; MOL=P+'/src/foldx/molecules'
CX='/tmp/claude-0/cx'; W='/tmp/claude-0/cx/work'; os.makedirs(W,exist_ok=True)
JOBS={'fold_2026_11_05_08_02':('vapBC36','WT',None),'fold_2026_11_05_08_04':('vapBC36','AF_mutant_vapC36_L126H',None),
      'fold_2026_11_05_08_08':('vapBC6','WT',None),'fold_2026_11_05_08_15':('vapBC6','AF_mutant_vapC6_W65L',None)}
MUT={'vapBC36':('L',126,'H','Rv1982c vapC36 p.Leu126His'),'vapBC6':('W',65,'L','Rv0656c vapC6 p.Trp65Leu')}
AA={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def cif2pdb(cif,pdb):
    atoms=[]
    for l in open(cif):
        if l.startswith('ATOM'):
            t=l.split(); atoms.append(t)
    with open(pdb,'w') as o:
        last=None
        for i,t in enumerate(atoms,1):
            _,id_,el,name,alt,res,ch,ent,seq,ins,x,y,z,occ,b=t[:15]; ch=t[16] if len(t)>16 else ch
            if last and last!=ch: o.write('TER\n')
            last=ch
            nm=(' '+name) if len(name)<4 else name
            o.write('ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%(i,nm,res,ch,int(seq),float(x),float(y),float(z),float(occ),float(b),el))
        o.write('TER\nEND\n')
def run(d,args): return subprocess.run([FX]+args,cwd=d,capture_output=True,text=True,timeout=14400)
def analyse(d,pdb):
    run(d,['--command=AnalyseComplex','--pdb='+pdb,'--analyseComplexChains=A,B','--output-dir=.'])
    f=d+'/Summary_%s_AC.fxout'%pdb[:-4]; 
    for l in open(f):
        t=l.rstrip('\n').split('\t')
        if len(t)>5 and t[0].startswith(pdb[:-4]) or (len(t)>5 and 'pdb' in t[0].lower() and t[0].endswith('.pdb')): 
            try: return float(t[5])
            except: pass
    return None
def mindist(pdb,chres):
    A=[];B=[]
    for l in open(pdb):
        if l.startswith('ATOM'):
            ch=l[21]; xyz=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
            if ch=='A' and int(l[22:26])==chres and l[12:16].strip() not in('N','C','O','CA'): A.append(xyz)
            if ch=='B': B.append(xyz)
    return min(math.dist(a,b) for a in A for b in B) if A else None
def plddt(pdb,chres):
    for l in open(pdb):
        if l.startswith('ATOM') and l[21]=='A' and int(l[22:26])==chres and l[12:16].strip()=='CA': return float(l[60:66])
out=[]
for z,(name,kind,_) in JOBS.items():
    d=W+'/'+z; shutil.rmtree(d,ignore_errors=True); os.makedirs(d); shutil.copytree(MOL,d+'/molecules')
    cif=glob.glob('%s/%s/*model_0.cif'%(CX,z))[0]; sc=json.load(open(glob.glob('%s/%s/*summary_confidences_0.json'%(CX,z))[0]))
    cif2pdb(cif,d+'/m.pdb'); run(d,['--command=RepairPDB','--pdb=m.pdb','--output-dir=.'])
    rp='m_Repair.pdb' if os.path.exists(d+'/m_Repair.pdb') else 'm.pdb'
    ie=analyse(d,rp); wt,pos,mt,label=MUT[name]
    row={'complex':name,'model':kind,'variant':label,'iptm':sc.get('iptm'),'ptm':sc.get('ptm'),'ranking_score':sc.get('ranking_score'),
         'interaction_energy_kcal_mol':ie,'mutated_residue_min_sidechain_distance_to_partner_A':round(mindist(d+'/m.pdb',pos),2) if mindist(d+'/m.pdb',pos) else None,
         'residue_pLDDT':plddt(d+'/m.pdb',pos),'complex_ddG_kcal_mol':'','complex_ddG_sd':'','interaction_energy_after_mutation':''}
    if kind=='WT':
        open(d+'/individual_list.txt','w').write('%sA%d%s;\n'%(wt,pos,mt))
        run(d,['--command=BuildModel','--pdb='+rp,'--mutant-file=individual_list.txt','--numberOfRuns=3','--output-dir=.'])
        vals=[]
        for l in open(d+'/Dif_%s.fxout'%rp[:-4]):
            t=l.split('\t')
            if len(t)>1 and t[0].endswith('.pdb'):
                try: vals.append(float(t[1]))
                except: pass
        if vals:
            m=sum(vals)/len(vals); row['complex_ddG_kcal_mol']=round(m,2); row['complex_ddG_sd']=round((sum((v-m)**2 for v in vals)/len(vals))**.5,2)
        # interaction energy of mutant complex built by FoldX (first repetition)
        mp=glob.glob(d+'/%s_%s_1.pdb'%(rp[:-4],'*'))+glob.glob(d+'/%s_1_0.pdb'%rp[:-4])
        if mp: row['interaction_energy_after_mutation']=analyse(d,os.path.basename(mp[0]))
    out.append(row); print(row,flush=True)
os.makedirs(P+'/cohort/c4',exist_ok=True)
w=csv.DictWriter(open(P+'/cohort/c4/a25d_foldx_alphafold_complexes.tsv','w'),fieldnames=list(out[0]),delimiter='\t'); w.writeheader(); w.writerows(out)
