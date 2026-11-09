"""C1 step 1: per-isolate callable-position arrays (depth>=10, outside static mask), packed."""
import numpy as np,subprocess,os,re,json
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'
O=P+'/cohort/c1'; L=4411532
mask=np.zeros(L+1,bool)
for l in open(P+'/ref/mask_static.bed'):
    f=l.split('\t'); mask[int(f[1])+1:int(f[2])+1]=True
print('masked bp',mask.sum())
np.save(O+'/mask.npy',np.packbits(mask))
samples=sorted([d for d in os.listdir(P+'/work') if re.match(r'MTB\d+$',d)],key=lambda x:int(x[3:]))
st={}
for s in samples:
    out=subprocess.run("zcat %s/work/%s/depth_all.txt.gz | cut -f3"%(P,s),shell=True,capture_output=True,text=True).stdout.split()
    d=np.array(out,dtype=np.int32); assert len(d)==L,(s,len(d))
    cal=np.zeros(L+1,bool); cal[1:]=d>=10
    st[s]={'mean_depth':round(float(d.mean()),1),'callable_unmasked':int((cal&~mask).sum()),'callable_frac_unmasked':round(float((cal&~mask).sum()/(L-mask.sum()+1)),4)}
    np.save('%s/cal_%s.npy'%(O,s),np.packbits(cal))
json.dump(st,open(O+'/callable_stats.json','w'),indent=1)
for s,x in sorted(st.items(),key=lambda t:t[1]['callable_frac_unmasked'])[:10]: print(s,x)
