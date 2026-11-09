import zipfile,gzip,sys,os,re,shutil
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
z=zipfile.ZipFile(sys.argv[1]); n=0; accs=[]
for name in z.namelist():
    m=re.match(r'ncbi_dataset/data/(GC[AF]_\d+\.\d)/.*\.(fna|fa)$',name)
    if not m: continue
    acc=m.group(1); out='%s/pub/fna/%s.fna.gz'%(P,acc)
    if os.path.exists(out): continue
    with z.open(name) as i, gzip.open(out+'.tmp','wb',compresslevel=3) as o: shutil.copyfileobj(i,o)
    os.rename(out+'.tmp',out)
    n+=1; accs.append(acc)
open(P+'/pub/new_accs.txt','w' if False else 'a').write('\n'.join(accs)+'\n')
print('ingested',n)
