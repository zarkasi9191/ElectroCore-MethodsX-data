import sys, copy, json, warnings, numpy as np
warnings.filterwarnings('ignore')
import tkstub
import plotting_export_sens as pes
sys.modules['plotting_export']=pes
import main_app_sens as ma
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
def autofit(f,Z,app):
    active={k:v for k,v in ma.EIS_MODELS.items() if int(k.split('.')[0])<=15}
    def stage(keys,nr,tol,nf):
        return [r for r in (app.fit_circuit(f,Z,active[k]['circuit'],active[k]['params'],k,n_restarts=nr,tol=tol,max_nfev=nf) for k in keys) if r['success']]
    coarse=stage(list(active),3,1e-6,300)
    top=[r['model_key'] for r in pes.rank_models_correctly(coarse)[:15]]
    fine=stage([k for k in active if k in top],12,1e-8,1000)
    fk={r['model_key'] for r in fine}
    return [r for r in coarse if r['model_key'] not in fk]+fine
configs=[('baseline',0.98,150,5,10),('ρsev 0.95',0.95,150,5,10),('ρsev 0.99',0.99,150,5,10),
         ('rel. error 100 %',0.98,100,5,10),('rel. error 200 %',0.98,200,5,10),
         ('window ×3',0.98,150,3,10),('window ×10',0.98,150,10,10),
         ('margin 5',0.98,150,5,5),('margin 15',0.98,150,5,15),
         ('all lenient',0.99,200,3,5),('all strict',0.95,100,10,15)]
out={}

for fn,true in [('Hasil_Fitting_1_Puncak.xlsx','7'),('Hasil_Fitting_2_Puncak.xlsx','13')]:
    f,Z,_=load(fn); app=make_app()
    res=autofit(f,Z,app)
    rows=[]
    for name,cs,rs,wf,mg in configs:
        ma._SENS.update(CS=cs,RS=rs,WF=wf); pes._SENS_MARGIN[0]=mg
        R=copy.copy(res)
        for r in R:
            idn=app.check_identifiability(r['params'],r['frequencies_used'],uncertainty=r.get('uncertainty'))
            r['identifiability']=idn
            am=r['advanced_metrics']; am['AIC_penalized']=am['AIC']+idn['penalty']
        rk=pes.rank_models_correctly(R)
        sel=rk[0]; tr=[(i+1,r) for i,r in enumerate(rk) if r['model_key'].split('.')[0]==true][0]
        n0=sum(1 for r in rk if r.get('_tier')==0)
        rows.append(dict(cfg=name,selected='M'+sel['model_key'].split('.')[0],sel_tier=sel.get('_tier'),
             sel_pen=round(sel['identifiability']['penalty'],2),true_rank=tr[0],true_tier=tr[1].get('_tier'),n_tier0=n0,
             top3=['M'+r['model_key'].split('.')[0]+'/T%s'%r.get('_tier') for r in rk[:3]]))
        print(fn[14],rows[-1],flush=True)
    out[fn]=rows
json.dump(out,open('sens_results.json','w'),indent=1,ensure_ascii=False)
