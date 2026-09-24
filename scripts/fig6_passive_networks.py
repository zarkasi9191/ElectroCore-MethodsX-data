import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ecore_io import ROOT, load_fit_export, load_measurement
from ecore_kk import kk_test
OUT = ROOT / "figures"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif'],'font.size':8,'axes.titlesize':8.5,
  'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7,'mathtext.fontset':'stix'})
B1='#1f6fe0'; R1='#d62728'
fig,axs=plt.subplots(2,3,figsize=(7.2,4.3),gridspec_kw=dict(width_ratios=[1.5,1,1],hspace=0.5,wspace=0.30,left=0.07,right=0.99,top=0.94,bottom=0.1))
rows=[('A','Dataset A','M8',(45e3,72e3),'abc'),('B','Dataset B','M13',None,'def')]
for r,(fn,lab,model,mask,L) in enumerate(rows):
    f,Z,Zf=load_fit_export(fn)
    d1=100*(Z.real-Zf.real)/np.abs(Z); d2=100*(Z.imag-Zf.imag)/np.abs(Z)
    kr=kk_test(f,Z,20); k1,k2,chi=100*kr['delta_re'],100*kr['delta_im'],kr['chi2_kk']
    a=axs[r,0]
    a.plot(Z.real,-Z.imag,'o',mfc='none',mec=B1,ms=3,mew=0.8,label='data (%d pts)'%len(f))
    a.plot(Zf.real,-Zf.imag,'-',color=R1,lw=1.2,label='fit, %s'%model)
    x0,x1=Z.real.min(),Z.real.max(); pad=0.04*(x1-x0)
    a.set_xlim(x0-pad,x1+pad); a.set_ylim(0,(-Z.imag).max()*1.18); a.set_aspect('equal',adjustable='box')
    a.set_xlabel("$Z'$ (Ω)"); a.set_ylabel("$-Z''$ (Ω)"); a.set_title('(%s) %s, selected fit'%(L[0],lab))
    a.legend(loc='lower center',frameon=True,fontsize=6.5); a.grid(alpha=0.25,lw=0.5)
    lo=min(d1.min(),d2.min(),k1.min(),k2.min(),-2.3)-0.3; hi=max(d1.max(),d2.max(),k1.max(),k2.max(),2.3)+0.3
    for j,((r1,r2),t) in enumerate([((d1,d2),'(%s) fit residuals'%L[1]),((k1,k2),'(%s) KK residuals, $M$ = 20'%L[2])]):
        ax=axs[r,j+1]
        ax.semilogx(f,r1,'o',color=B1,ms=2.6,label=r'$\Delta_{re}$')
        ax.semilogx(f,r2,'s',mfc='none',mec=R1,ms=2.6,mew=0.8,label=r'$\Delta_{im}$')
        ax.axhline(0,color='k',lw=0.6)
        for y in (-2,-1,1,2): ax.axhline(y,color='orange',ls='--',lw=0.5)
        if mask:
            ax.axvspan(*mask,color='0.85',zorder=0,lw=0)
            ax.text(np.sqrt(mask[0]*mask[1]),hi-0.1,'masked',ha='center',va='top',fontsize=6.3,rotation=90)
        ax.set_ylim(lo,hi); ax.set_xlabel('Frequency (Hz)'); ax.set_title(t,pad=3); ax.grid(alpha=0.25,lw=0.5)
        if j==0: ax.set_ylabel('residual (%)'); ax.legend(loc='lower left',ncol=2,frameon=True,fontsize=6.5)
        else: ax.text(0.04,0.06,r'$\chi^2_{KK}$ = %.1f × 10$^{-5}$'%(chi*1e5),transform=ax.transAxes,fontsize=7,bbox=dict(fc='white',ec='0.7',lw=0.5,pad=2))
    print(lab,'chi2KK %.2e'%chi)
fig.savefig(OUT/'Fig6_passive_networks.png',dpi=400); print('ok')
