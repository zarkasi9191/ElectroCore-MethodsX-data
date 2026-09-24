import numpy as np, matplotlib
import sys, csv, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ecore_io import ROOT, load_fit_export, load_measurement
from ecore_kk import kk_test
OUT = ROOT / "figures"; OUT.mkdir(exist_ok=True)

matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif'],'mathtext.fontset':'stix','font.size':8})
# (a) dataset B, from Hasil_Fitting_2_Puncak.xlsx (Auto_Fit_Ranking)
names=['M13 (true)','M14','M15','M9','M8']
from ecore_io import table
_h,_r=table('B','Auto_Fit_Ranking'); _rk={r[1].split('.')[0]:dict(zip(_h,r)) for r in _r}
raw=[_rk[m]['AIC'] for m in ('13','14','15','9','8')]
pen=[_rk[m]['AIC_penalized'] for m in ('13','14','15','9','8')]
# (b),(c) dataset B correlation matrices (M13, M14) from the software's uncertainty routine, 2 d.p.
lb=['$R_s$','$R_1$','$C_1$','$R_2$','$C_2$']
B=[[float(x) for x in row[1:i+1]] for i,row in enumerate(list(csv.reader(open(ROOT/'results'/'correlation_M13_datasetB.csv')))[1:])]
lc=['$R_s$','$R_1$','$Q_1$','$n_1$','$R_2$','$Q_2$','$n_2$']
C=[[float(x) for x in row[1:i+1]] for i,row in enumerate(list(csv.reader(open(ROOT/'results'/'correlation_M14_datasetB.csv')))[1:])]
fig=plt.figure(figsize=(6.6,4.25))
gs=fig.add_gridspec(2,2,height_ratios=[0.62,1.0],hspace=0.30,wspace=0.28,left=0.125,right=0.9,top=0.95,bottom=0.075)
a=fig.add_subplot(gs[0,:])
y=np.arange(len(names))[::-1]
for n_,r,p,yy in zip(names,raw,pen,y):
    if p!=r:
        a.annotate('',xy=(p,yy),xytext=(r,yy),arrowprops=dict(arrowstyle='->',color='#c0392b',lw=1.0))
        a.text((r+p)/2,yy+0.22,'+%.0f'%(p-r),ha='center',fontsize=7,color='#c0392b')
a.plot(raw,y,'o',mfc='white',mec='0.35',ms=6,label='raw AIC')
a.plot(pen,y,'o',color='#1e3a5f',ms=6,label='penalised AIC')
a.axvline(-1021.41,color='k',ls=':',lw=0.9)
a.text(-1019,4.25,'best penalised AIC',fontsize=6.8,va='center')
a.set_yticks(y); a.set_yticklabels(names); a.set_xlabel('AIC',labelpad=2); a.set_ylim(-0.6,4.7)
a.legend(loc='upper right',fontsize=7,frameon=True,borderpad=0.4,handletextpad=0.3)
a.set_title('(a) the penalty displaces the raw-AIC minimum (M9) past the true circuit',fontsize=8,pad=3)
def mat(ax,labels,L,title,solid=(),dotted=(),cbar=False):
    n=len(labels); M=np.full((n,n),np.nan)
    for i,row in enumerate(L):
        for j,v in enumerate(row): M[i,j]=v
    im=ax.imshow(M,cmap='RdBu_r',vmin=-1,vmax=1)
    for i,row in enumerate(L):
        for j,v in enumerate(row):
            ax.text(j,i,('%.2f'%v).replace('-','\u2212'),ha='center',va='center',fontsize=5.0 if n>7 else (6.0 if n>5 else 6.8),color='white' if abs(v)>0.7 else 'black')
    for (i,j) in solid: ax.add_patch(Rectangle((j-.5,i-.5),1,1,fill=False,lw=1.6,ec='k'))
    for (i,j) in dotted: ax.add_patch(Rectangle((j-.5,i-.5),1,1,fill=False,lw=1.3,ec='k',ls=':'))
    ax.set_xticks(range(n)); ax.set_xticklabels(labels); ax.set_yticks(range(n)); ax.set_yticklabels(labels)
    ax.tick_params(length=0,pad=1.5)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(title,fontsize=8,pad=3)
    return im
ab=fig.add_subplot(gs[1,0]); mat(ab,lb,B,'(b) selected M13, penalty 0')
ac=fig.add_subplot(gs[1,1]); im=mat(ac,lc,C,'(c) M14, penalty 7.9 (Tier 1)',solid=[(4,1)],dotted=[(3,2),(6,5)])
cb=fig.colorbar(im,ax=ac,fraction=0.046,pad=0.03); cb.set_label(r'$\rho_{ij}$',fontsize=8,labelpad=2); cb.ax.tick_params(labelsize=6.5)
fig.text(0.5,0.012,'Solid box: cross-group |ρ| > 0.90, penalized (soft).  Dotted: within-arc CPE pair, exempt.',ha='center',fontsize=6.8)
fig.savefig(OUT/'Fig7_identifiability.png',dpi=400); print('ok')
