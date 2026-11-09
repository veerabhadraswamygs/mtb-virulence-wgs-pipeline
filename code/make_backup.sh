#!/bin/bash
# Rebuild the restore kit. Same filenames every time, so the copies in
# Downloads\PIPELINE_BACKUP are overwritten rather than accumulating.
set -e
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
B=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/backup
mkdir -p $B && cd $P
rm -f $B/pipeline_core.tar.gz.part-* $B/per_sample_results.tar.gz* $B/snpeff_db.tar.gz
tar -czf $B/.core.tgz --exclude='ref/af_pdb' --exclude='ref/H37Rv.fa.bwt' \
    --exclude='ref/H37Rv.fa.sa' --exclude='ref/H37Rv.fa.pac' \
    --exclude='ref/H37Rv.fa.amb' --exclude='ref/H37Rv.fa.ann' code ref
( cd $B && split -b 20m .core.tgz pipeline_core.tar.gz.part- && rm -f .core.tgz )
tar -czf $B/snpeff_db.tar.gz -C $P/src/snpEff data
# every sample that has a verified report
rm -rf $B/stage && mkdir -p $B/stage/results && cd $B/stage
for S in $(ls $P/out 2>/dev/null); do
  [ -f "$P/out/$S/${S}_complete_analysis_report_2.pdf" ] || continue
  mkdir -p results/$S
  cp /mnt/user-data/outputs/${S}_report2/${S}.filt.vcf.gz       results/$S/ 2>/dev/null || true
  cp /mnt/user-data/outputs/${S}_report2/${S}_denovo_p2.fasta   results/$S/ 2>/dev/null || true
  cp $P/work/$S/results.json results/$S/ 2>/dev/null || true
  cp $P/work/$S/cov.json     results/$S/ 2>/dev/null || true
  cp -r $P/out/$S/tables     results/$S/ 2>/dev/null || true
  cp $P/logs/$S.check        results/$S/ 2>/dev/null || true
done
tar -czf $B/per_sample_results.tar.gz results
( cd $B && split -b 20m per_sample_results.tar.gz per_sample_results.tar.gz.part- && rm -f per_sample_results.tar.gz )
cd $B && rm -rf stage
rm -f $B/precohort_results.tar.gz*
tar -czf $B/precohort_results.tar.gz -C $P --exclude="cal_*.npy" --exclude="*.mldist" --exclude="*.ckp.gz" --exclude="*.bionj" --exclude="*uniqueseq.phy" --exclude="*.splits.nex" --exclude="cohort/c3/aln" --exclude="cohort/c3/cml" --exclude="cohort/c3/hyphy" --exclude="cohort/c4/foldx" --exclude="cohort/c4/tree/*.mldist" --exclude="cohort/c1/treeX2.mldist" cohort pub/calls pub/acc_list.txt pub/ncbi_dataset/data/assembly_data_report.jsonl 2>/dev/null
( cd $B && split -b 20m precohort_results.tar.gz precohort_results.tar.gz.part- && rm -f precohort_results.tar.gz )
md5sum precohort_results.tar.gz.part-* pipeline_core.tar.gz.part-* snpeff_db.tar.gz per_sample_results.tar.gz.part-* RESTORE.md PIPELINE_CHANGELOG.md COHORT_NOTES.md > MANIFEST.md5
echo "samples covered: $(cat $B/per_sample_results.tar.gz.part-* | tar -tzf - | grep -c 'results/MTB[0-9]*/$')"
ls -la $B | awk '{print $5,$9}'
