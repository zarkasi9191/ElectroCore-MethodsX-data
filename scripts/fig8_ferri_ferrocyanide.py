import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ecore_io import ROOT, load_fit_export, load_measurement
from ecore_kk import kk_test
OUT = ROOT / "figures"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif'],'mathtext.fontset':'stix','font.size':8})
cases=[(10,10,1e4,'#1f77b4'),(15,30,1e4,'#2ca02c'),(20,30,1e4,'#d62728')]   # c, certified fmin, fmax, colour
fig,(a,b)=plt.subplots(1,2,figsize=(7.0,3.1),gridspec_kw=dict(width_ratios=[1.0,1.3],wspace=0.26,left=0.075,right=0.99,top=0.91,bottom=0.15))
HF=[];HC=[]
for c,lo,hi,col in cases:
    f,Z=load_measurement(c)
    m=(f>=lo*0.999)&(f<=hi*1.001)
    a.plot(Z.real,-Z.imag,'o',ms=3,color=col,label='%d mM'%c)
    rf=kk_test(f,Z,20); rc=kk_test(f[m],Z[m],20)
    e=int(np.floor(np.log10(rc['chi2_kk']))); mant=rc['chi2_kk']/10**e
    hf,=b.semilogx(f,100*np.hypot(rf['delta_re'],rf['delta_im']),'-',color=col,lw=0.9,alpha=0.5,label='%d mM, full window'%c)
    hc,=b.semilogx(f[m],100*np.hypot(rc['delta_re'],rc['delta_im']),'o-',ms=2.6,lw=0.9,color=col,
                   label='%d mM, validated: $\\chi^2_{KK}$ = %.1f×10$^{%d}$'%(c,mant,e))
    HF.append(hf); HC.append(hc)
a.set_xlim(-60,2850); a.set_ylim(-60,1250); a.set_aspect('equal',adjustable='box'); a.set_xlabel("$Z'$ (Ω)"); a.set_ylabel("$-Z''$ (Ω)")
a.set_title('(a) measured spectra',fontsize=8)
a.legend(loc='lower right',fontsize=7,frameon=True)
b.axvspan(30,1e4,color='0.92',lw=0,zorder=0)
b.text(550,0.0125,'validated window (10 mM: from 10 Hz)',ha='center',va='bottom',fontsize=6.6,color='0.35')
b.set_yscale('log'); b.set_ylim(0.01,3000); b.set_xlabel('Frequency (Hz)'); b.set_ylabel('|Δ| (%)')
b.set_title('(b) Kramers–Kronig residual magnitude',fontsize=8)
b.legend(HF+HC,[h.get_label() for h in HF+HC],ncol=2,loc='upper left',fontsize=6.6,frameon=True,columnspacing=1.0,handlelength=1.8)
fig.savefig(OUT/'Fig8_ferri_ferrocyanide.png',dpi=400); print('ok')
