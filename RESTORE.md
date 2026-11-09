# Restore guide

This repository holds the code and the version records only. Sequence reads, VCFs, assemblies and result tables are deposited separately (see the manuscript's Availability of data and materials).

## What to install
Tool versions are listed in `env/versions.json` (fastp 1.4.0, BWA 0.7.18, samtools 1.24, bcftools 1.24, freebayes 1.3.10, minimap2 2.31, BLAST+ 2.17.0, DIAMOND 2.2.8, IQ-TREE 2.4.0, Prodigal 2.6.3, PAML 4.10.10, SnpEff 5.4c, Python 3.13) and the databases with download dates in the same file. Python packages: `env/python_packages.txt`.

## Layout
Place the repository as `<root>/code` and `<root>/ref`, then rebuild the reference files named in `env/versions.json` (H37Rv NC_000962.3 FASTA and GenBank, tbdb, VFDB, UniProt proteome UP000001584, AlphaFold DB models, IEDB export) under `<root>/ref`. Build the BWA index with `bwa index ref/H37Rv.fa`.

## Frozen panel check
`ref/panel/` holds the four frozen tier BED files (Tier 1a 211, 1b 120, 2 168, 3 447 genes). Their MD5 sums must be:

```
tier1a  a73ead619009c734f7d2a25cc64c04ed
tier1b  dff8d6ae0d20f7b4c4dddf5061608f03
tier2   75917c35aa6f89b09070547b1e7bf09f
tier3   c30cdbf19adfe295f17f627e6214b4e6
```

`code/check.py` re-verifies these for every sample. The pipeline fingerprint of release 1.0.0 is `74d12c396d65013d` (`env/pipeline_version.json`, `code/check_frozen.sh`).

## Run
Per sample: `code/run_sample.sh <sample>`; cohort analyses are the `c1_` to `c6_` scripts, in numerical order. Scripts contain the working paths of the environment where they were run; set `P` (project root) at the top of a script if your layout differs. In `build_reg.py`, `build_reg2.py` and `build_tss.py` the folder of published supplementary tables was renamed to `supplementary_tables` for this release.

## Integrity
`MANIFEST.md5` lists checksums of every file in this release: `md5sum -c MANIFEST.md5`.
