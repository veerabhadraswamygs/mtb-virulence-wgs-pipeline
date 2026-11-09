import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, json
P='/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2'; O=P+'/cohort/c4'
A={};cur=None
for l in open(O+'/tbd1_L1_local_plus_public_aligned.fa'):
    if l[0]=='>': cur=l[1:].strip(); A[cur]=[]
    else: A[cur].append(l.strip())
A={k:''.join(v).upper() for k,v in A.items()}
ref=A.pop('canettii_ref'); L=len(ref)
loc=[k for k in A if k.startswith('MTB')]; pub=[k for k in A if not k.startswith('MTB')]
names=loc+pub
M=np.zeros((len(names),L),int)  # 0 same, 1 differ, 2 missing(N), 3 gap
for i,n in enumerate(names):
    for j,(a,b) in enumerate(zip(A[n],ref)):
        M[i,j]=2 if a=='N' else (3 if a=='-' else (0 if a==b else 1))
fig,ax=plt.subplots(2,1,figsize=(11,12),gridspec_kw={'height_ratios':[1,14]},sharex=True)
ax[0].set_ylim(0,1); ax[0].axis('off')
for st,en,lab,c in ((1,399,"mmpS6 (5'-truncated)",'#4C78A8'),(396,2138,"mmpL6 5' portion",'#72B7B2')):
    ax[0].add_patch(plt.Rectangle((st,0.2),en-st,0.5,color=c)); ax[0].text((st+en)/2,0.45,lab,ha='center',va='center',color='white',fontsize=9)
cols=matplotlib.colors.ListedColormap(['#EEEEEE','#D62728','#BBBBBB','#000000'])
ax[1].imshow(M,aspect='auto',cmap=cols,vmin=0,vmax=3,extent=[0.5,L+0.5,len(names)-0.5,-0.5],interpolation='nearest')
ax[1].set_yticks(range(len(names))); ax[1].set_yticklabels([n if n.startswith('MTB') else '' for n in names],fontsize=5)
ax[1].axhline(len(loc)-0.5,color='k',lw=0.8); ax[1].set_xlabel('Position in M. canettii TbD1 reference segment (bp)')
ax[1].text(L+10,len(loc)/2,'local L1 (n=%d)'%len(loc),rotation=90,va='center',fontsize=8); ax[1].text(L+10,len(loc)+len(pub)/2,'public L1 (n=%d)'%len(pub),rotation=90,va='center',fontsize=8)
ax[1].set_xlim(0,L+60)
import matplotlib.patches as mp
ax[0].legend(handles=[mp.Patch(color='#EEEEEE',label='identical to canettii reference'),mp.Patch(color='#D62728',label='differs'),mp.Patch(color='#BBBBBB',label='not recovered (N)')],loc='lower center',fontsize=7,ncol=3,framealpha=.9,bbox_to_anchor=(0.5,1.0),bbox_transform=ax[0].transAxes)
fig.suptitle('TbD1 locus across lineage 1: %d local and %d public assemblies'%(len(loc),len(pub)),fontsize=11)
fig.savefig(O+'/fig_tbd1_locus_alignment_L1.png',dpi=220,bbox_inches='tight'); fig.savefig(O+'/fig_tbd1_locus_alignment_L1.svg',bbox_inches='tight')
