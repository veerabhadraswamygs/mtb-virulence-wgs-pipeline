import subprocess,json,os,re
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
def first(cmd):
    r=subprocess.run(cmd,shell=True,capture_output=True,text=True)
    t=(r.stdout or r.stderr).strip().splitlines()
    return t[0] if t else ''
V={
 'fastp':first('%s/bin/fastp --version'%P),
 'bwa':re.sub(r'.*Version: ','bwa ',first('%s/bin/bwa 2>&1 | grep Version'%P)),
 'samtools':first('%s/bin/samtools --version'%P),
 'bcftools':first('%s/bin/bcftools --version'%P),
 'minimap2':'minimap2 '+first('%s/bin/minimap2 --version'%P),
 'freebayes':first('%s/bin/freebayes --version'%P),
 'BLAST+':first('%s/bin/blastn -version'%P),
 'DIAMOND':first('%s/bin/diamond --version'%P),
 'IQ-TREE':first('%s/bin/iqtree2 --version'%P),
 'Prodigal':first('%s/bin/prodigal -v 2>&1 | sed -n 2p'%P),
 'PAML codeml':'PAML 4.10.10 (codeml)',
 'SnpEff':first('java -jar %s/src/snpEff/snpEff.jar -version 2>/dev/null'%P) or 'SnpEff (snpEff_latest_core)',
 'Python':first('python3 -V'),
}
DB={
 'Reference genome':'H37Rv NC_000962.3 (RefSeq equivalent of AL123456.3), 4,411,532 bp',
 'Annotation':'NCBI RefSeq GFF + Mycobrowser H37Rv v5 (functional categories)',
 'SnpEff database':'H37Rv, built locally from the NC_000962.3 GenBank record',
 'tbdb':'TB-Profiler database (tbdb-master): Coll barcode 1,111 positions; 49,330 catalogued mutations; 9,658 SITVIT2 spoligotype patterns',
 'WHO catalogue':'WHO mutation catalogue v2 (via tbdb), five confidence grades',
 'VFDB':'core (setA) and full (setB), nucleotide and protein; database updated 23 Oct 2026, downloaded 27 Oct 2026',
 'BV-BRC':'release 3.59.14, genome 83332.12, specialty genes, downloaded 27 Oct 2026',
 'IEDB':'MTBC, linear peptides, positive T-cell assays, host Homo sapiens; queried 27 Oct 2026',
 'UniProt':'proteome UP000001584, feature columns (domains, regions, motifs, sites, TM, signal), downloaded 27 Oct 2026',
 'AMRFinder':'AMRProt reference protein set (NCBI)',
 'TSS reference':'Shell et al. 2015 (PLoS Genet) S3 Table, authors’ Leadered / Leaderless classification',
 'ChIP-seq':'Minch et al. 2015 (Nat Commun) Supplementary Table 1, 15,980 peaks, 156 regulators',
}
json.dump({'tools':V,'databases':DB},open(P+'/ref/versions.json','w'),indent=1)
for k,v in V.items(): print('%-14s %s'%(k,v))
