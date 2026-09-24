import sys, json, warnings, numpy as np
warnings.filterwarnings('ignore')
import tkstub
import plotting_export_sens as pes
sys.modules['plotting_export']=pes
import main_app_sens as ma, weighting
from kk import load
class V:
    def __init__(s,v): s.v=v
    def get(s): return s.v
def make_app():
    a=object.__new__(ma.ElectrocoreAnalyzer)
    a.excluded_indices=set(); a.user_initial_guess={}
    a.weighting_scheme=V('modulus'); a.irls_iterative=V(True); a.irls_max_iter=V(20)
    a.log_message=lambda *x,**k: None
    return a
out={}
for fn,mods in [('Hasil_Fitting_2_Puncak.xlsx',['13.','14.']),('Hasil_Fitting_1_Puncak.xlsx',['9.','13.','14.','15.','2.'])]:
    f,Z,_=load(fn); app=make_app()
    for mk in mods:
        k=[k for k in ma.EIS_MODELS if k.startswith(mk)][0]
        weighting.LAST_CORR.clear()
        r=app.fit_circuit(f,Z,ma.EIS_MODELS[k]['circuit'],ma.EIS_MODELS[k]['params'],k,n_restarts=12,tol=1e-8,max_nfev=1000)
        cs=r['circuit_string']; nm,C=weighting.LAST_CORR.get(cs,(None,None))
        print('==',fn[14],k[:40],'pen %.2f'%r['identifiability']['penalty'])
        for d in r['identifiability']['details']: print('    ',d[:160])
        if C is not None and fn.startswith('Hasil_Fitting_2'):
            out[mk]={'names':nm,'C':C.tolist()}
            print('   ',nm); print(np.round(C,3))
json.dump(out,open('corrB.json','w'),indent=1)
