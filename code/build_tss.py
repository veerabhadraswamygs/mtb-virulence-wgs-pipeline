"""Build the TSS / regulatory reference from Cortes 2013, Shell 2015 and Minch 2015."""
import json,sys,re
from openpyxl import load_workbook
U='/mnt/user-data/uploads/Downloads/limitation fix files'
R='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2/ref'
def sheets(path):
    wb=load_workbook(path,read_only=True,data_only=True); return wb
for tag,path in [('cortes_S1',U+'/supplementary_tables/1-s2.0-S2211124713006153-mmc2.xlsx'),
                 ('cortes_S6',U+'/supplementary_tables/1-s2.0-S2211124713006153-mmc6.xlsx'),
                 ('shell_S3', U+'/supplementary_tables/pgen.1005641.s010.xlsx'),
                 ('minch',    U+'/ncomms6829-s2.xlsx')]:
    wb=sheets(path)
    print('==',tag,'sheets:',wb.sheetnames[:6])
    ws=wb[wb.sheetnames[0]]
    for i,row in enumerate(ws.iter_rows(values_only=True)):
        print('   r%d:'%i,[str(x)[:22] if x is not None else '' for x in row[:9]])
        if i>=3: break
    wb.close()
