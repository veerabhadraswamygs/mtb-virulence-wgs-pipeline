#!/bin/bash
# Package a finished sample into /mnt/user-data/outputs ready for delivery.
set -u
S=$1
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
O=/mnt/user-data/outputs
[ -f "$P/out/$S/${S}_complete_analysis_report_2.pdf" ] || exit 1
mkdir -p "$O/${S}_report2/tables"
cp $P/out/$S/${S}_complete_analysis_report_2.pdf $O/${S}_report2/
cp $P/out/$S/tables/*.tsv $O/${S}_report2/tables/ 2>/dev/null
for f in ${S}_denovo_p2.fasta ${S}.filt.vcf.gz ${S}.raw.vcf.gz ${S}_TbD1.fasta; do
  [ -f "$P/work/$S/$f" ] && cp "$P/work/$S/$f" $O/${S}_report2/
done
cd $O/${S}_report2 && python3 -I - "$S" <<'PY'
import csv,glob,os,sys
from openpyxl import Workbook
from openpyxl.styles import Font
S=sys.argv[1]
wb=Workbook(); wb.remove(wb.active)
nm={'variants':'Variants','subclonal':'Sub-clonal','resistance':'Resistance','gene_presence':'Panel genes',
    'deletions':'Deletions','is6110_junctions':'IS6110','gene_coverage':'Gene coverage','tbd1':'TbD1',
    'pe_ppe':'PE-PPE','identification_rd':'RD profile'}
for f in sorted(glob.glob('tables/*.tsv')):
    k=os.path.basename(f).replace(S+'_','').replace('.tsv','')
    ws=wb.create_sheet(nm.get(k,k)[:31])
    for r in csv.reader(open(f),delimiter='\t'): ws.append(r)
    if ws.max_row: 
        for c in ws[1]: c.font=Font(bold=True)
    ws.freeze_panes='A2'
wb.save('%s_all_tables_2.xlsx'%S)
PY
echo "packaged $S"
