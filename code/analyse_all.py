"""Run every per-sample analysis once and cache to work/<sample>/results.json."""
import json,os,sys,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import analyse as A, mod_lineage as ML, mod_resist as MR, mod_tbd1 as MT, mod_pe as MP, mod_id as MI
P=A.P
def run(sample,force=False):
    W=P+'/work/'+sample; out=W+'/results.json'
    if os.path.exists(out) and not force: return json.load(open(out))
    t0=time.time(); R={}
    cov=json.load(open(W+'/cov.json'))
    R['coverage']={k:cov[k] for k in ('deletions','noncallable','mean_depth_mqfiltered','mean_depth_all')}
    R['variants']=A.variant_table(sample,W)
    R['subclonal']=A.subclonal(sample,W)
    R['purity']=A.purity(sample,W,W+'/%s_fastp.json'%sample)
    R['exclusion']=A.exclusion(cov,R['purity'])
    R['decontamination']=A.decontamination(cov,R['purity'])
    R['assembly']=A.assembly_stats(sample)
    R['completeness']=A.completeness(cov)
    R['epitope_partition']=A.partition_pi(sample,R['variants'],cov)
    R['is6110']=A.is6110(sample,W)
    R['is6110_copy_number']=A.is6110_copy_number(W)
    R['catalogue_presence']=A.catalogue_presence(sample,cov)
    R['tier_callable']=A.tier_callable(R['catalogue_presence'])
    R['amrfinder']=A.amrfinder(sample)
    R['lineage']=ML.lineage(sample,W)
    R['spoligotype']=ML.spoligotype(sample,W)
    h,v,g=MR.call(W+'/%s.ann.vcf'%sample)
    R['resistance']=h; R['resistance_uncertain']=v; R['resistance_genes']=len(g)
    R['fastp']=json.load(open(W+'/%s_fastp.json'%sample))['summary']
    R['fastp_filtering']=json.load(open(W+'/%s_fastp.json'%sample))['filtering_result']
    R['identification']=MI.run(sample)
    R['tbd1']=MT.run(sample)
    R['pe_ppe']=MP.run(sample)
    R['_seconds']=round(time.time()-t0,1)
    json.dump(R,open(out,'w'))
    return R
if __name__=='__main__':
    R=run(sys.argv[1],force='--force' in sys.argv)
    print('computed in %ss'%R['_seconds'])
    for k in ['variants','subclonal','resistance','resistance_uncertain']: print('  %-22s %d'%(k,len(R[k])))
    print('  lineage',R['lineage']['lineage'],'| SIT',R['spoligotype']['sit'],
          '| IS6110 cn',R['is6110_copy_number']['estimated_copy_number'],
          '| excluded',R['exclusion']['excluded'])
