from pathlib import Path
import sys,json,numpy as np,time
r=Path(__file__).parent;sys.path.insert(0,str(r/'stage'))
from frozen_sigma_tangent import CachedFrozenSigmaResponse
models={}
for d in ['MnO_center','MnO_m_0.001','MnO_p_0.001']:
 m,_=CachedFrozenSigmaResponse.from_export(r/'native'/'MnO'/d/'response',diagnostic_gap_cutoff=1e-10,cache_mb=0,precontract=False)
 models[d]=m
refmu=models['MnO_m_0.001'].meta.mu
factor=14710.507641078/(3*148.11716074988342)
rows=[]
base={}
for shift in [0,-.001,.001,-.01,.01]:
 mu=refmu+shift/13.605693122994;values={}
 for d,m in models.items():
  dh=np.broadcast_to(-(mu-m.meta.mu)*np.eye(m.eps.shape[1]),(len(m.eps),m.eps.shape[1],m.eps.shape[1])).copy()
  scalar=m.lattice_scalar(dh,t=1,u=m.uc)+mu*21
  rho=m.density(dh,t=1,u=m.uc)
  values[d]=dict(K_lattice_plus_muN0_Ry=scalar,N=float(np.trace(rho,axis1=1,axis2=2).real.sum()))
 if shift==0:base=values
 delta=((values['MnO_p_0.001']['K_lattice_plus_muN0_Ry']-base['MnO_p_0.001']['K_lattice_plus_muN0_Ry'])-(values['MnO_m_0.001']['K_lattice_plus_muN0_Ry']-base['MnO_m_0.001']['K_lattice_plus_muN0_Ry']))/.002*factor
 row=dict(mu_shift_eV=shift,reference_mu_Ry=refmu,stress_sensitivity_GPa=delta,values=values)
 rows.append(row);print(json.dumps(row),flush=True)
(r/'charge_mu_probe.json').write_text(json.dumps(dict(scope='fixed saved potential/charge and frozen Sigma: sensitivity only, not a relaxed fixed-N correction',rows=rows),indent=2)+'\n')
