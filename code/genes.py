"""Build the H37Rv gene index: Rv locus -> coords, name, product, category, protein_id."""
import re, gzip, json, sys
REF='%s/ref'%sys.argv[1]
genes={}
# --- RefSeq GFF: coordinates + protein_id ---
for line in open(REF+'/H37Rv_refseq.gff'):
    if line.startswith('#'): continue
    f=line.rstrip('\n').split('\t')
    if len(f)<9 or f[2] not in ('CDS','gene'): continue
    a=dict(kv.split('=',1) for kv in f[8].split(';') if '=' in kv)
    lt=a.get('locus_tag')
    if not lt: continue
    g=genes.setdefault(lt,{'rv':lt,'start':None,'end':None,'strand':f[6],'name':None,
                           'product':None,'protein_id':None,'category':None,'type':None})
    s,e=int(f[3]),int(f[4])
    if g['start'] is None or s<g['start']: g['start']=s
    if g['end'] is None or e>g['end']:     g['end']=e
    g['strand']=f[6]
    if f[2]=='CDS':
        g['type']='CDS'
        if a.get('protein_id'): g['protein_id']=a['protein_id']
        if a.get('product'):    g['product']=a['product'].replace('%2C',',')
    if a.get('gene'): g['name']=a['gene']
print('RefSeq loci:',len(genes))
# --- Mycobrowser: functional category + gene name ---
mb=0
for line in open(REF+'/mycobrowser.gff'):
    if line.startswith('#'): continue
    f=line.rstrip('\n').split('\t')
    if len(f)<9: continue
    a={}
    for kv in f[8].split(';'):
        if '=' in kv:
            k,v=kv.split('=',1); a[k.strip()]=v.strip()
    lt=a.get('Locus')
    if not lt: continue
    g=genes.setdefault(lt,{'rv':lt,'start':int(f[3]),'end':int(f[4]),'strand':f[6],'name':None,
                           'product':None,'protein_id':None,'category':None,'type':f[2]})
    if a.get('Functional_Category'): g['category']=a['Functional_Category']; mb+=1
    if a.get('Name') and not g['name']: g['name']=a['Name']
    if a.get('Product') and not g['product']: g['product']=a['Product']
    if g['type'] is None: g['type']=f[2]
    if g['start'] is None: g['start']=int(f[3]); g['end']=int(f[4]); g['strand']=f[6]
print('after Mycobrowser:',len(genes),'| with category:',mb)
cds=[g for g in genes.values() if g['type']=='CDS']
print('CDS:',len(cds),'| with protein_id:',sum(1 for g in cds if g['protein_id']))
import collections
print('categories:',len(collections.Counter(g['category'] for g in genes.values() if g['category'])))
json.dump(genes,open(REF+'/genes.json','w'))
# protein_id -> rv map
pid={g['protein_id']:g['rv'] for g in genes.values() if g.get('protein_id')}
json.dump(pid,open(REF+'/protid2rv.json','w'))
print('protein_id map:',len(pid))
