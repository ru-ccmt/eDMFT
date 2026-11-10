"""Use already converged endpoints ONLY as an independent response benchmark."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
import sys,json
from pathlib import Path
import numpy as np
from setup import ROOT,SOURCES
from test_kernel import load
from potential_probe import PotentialProbe
from wien_response_density import Density
from screen_response import assemble,charge

PATHS={
 'SrVO3':('/Users/haule/projects/materials/SrVO3','V_strain_m001','V_strain_p001'),
 'FeSe':('/Users/haule/projects/materials/FeSe_again/FeSe_stress','FeSe_RELA_strain_m001','FeSe_RELA_strain_p001')}
def run(case):
 f,d=load(case);p=PotentialProbe(case,d);base,mn,pn=PATHS[case];rows=[]
 for name in [pn,mn]:
  dr=Path(base)/name;cell=np.array([float(dr.joinpath(case+'.struct').read_text().splitlines()[3][j:j+10]) for j in [0,10,20]]);lam=float(np.mean(cell/p.cell))
  rho=Density(dr/(case+'.clmsum'),dr/(case+'.struct'));rows.append((rho.transport(lam,True),lam))
 (plus,lp),(minus,lm)=rows;ref=plus.add(minus,-1).scale(1/(lp-lm));ref.write(ROOT/case/'reference_q.density')
 op,core=p.apply_kernel(ref,'reference_induced',core=True)
 aff=np.load(ROOT/case/'affine_operator_0.0002500.npy')
 z=np.load('/private/tmp/projector_stationarity_probe_20261107/band_matrix_evidence.npz') if case=='SrVO3' else np.load('/private/tmp/projector_cross_material_20261107/FeSe/matrices.npz')
 R=z['response']+z['norm_response'];fac=14710.507641078/(3*p.vol)
 dr,info=f.response(op,True);val=assemble(case,f,d,dr,'reference_induced_density');response=p.rho.zeros().add(val).add(core)
 source=Density(ROOT/case/'source_0.0002500.density',p.structure);pred=source.copy().add(response);err=pred.copy().add(ref,-1)
 result={'kernel_correction_GPa':float(np.einsum('kij,kji->',R,op).real)*fac,'affine_plus_reference_kernel_GPa':float(np.einsum('kij,kji->',R,op+aff).real)*fac,'reference_charge':charge(ref,p),'source_norm':float(np.linalg.norm(source.vector())),'reference_norm':float(np.linalg.norm(ref.vector())),'fixed_point_residual_norm':float(np.linalg.norm(err.vector())),'resolvent':info}
 (ROOT/case/'reference_check.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True);pred.write(ROOT/case/'reference_predicted_q.density');np.save(ROOT/case/'reference_induced.npy',op)
if __name__=='__main__':run(sys.argv[1])
