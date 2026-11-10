"""Contract the portable driver's prediction with the archived adjoint."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json
import numpy as np
from pathlib import Path
from setup import ROOT
from test_kernel import load
from potential_probe import PotentialProbe
from wien_response_density import Density

def run(case,folder):
    workspace=Path(folder);f,d=load(case);p=PotentialProbe(case,d)
    q=Density(workspace/'response.density',p.structure)
    config=json.loads((workspace/'config.json').read_text())
    step=config.get('probe_step',1e-5);h=config.get('strain',.00025)
    kv,_=p.apply_kernel(q,workspace.name+'_contraction',step=step)
    aff=np.load(ROOT/case/f'affine_operator_{h:.7f}.npy')
    if case=='SrVO3':z=np.load('/private/tmp/projector_stationarity_probe_20261107/band_matrix_evidence.npz');analytic=5.67688415;numeric=2.933836407;refcor=-2.9867588495
    else:z=np.load('/private/tmp/projector_cross_material_20261107/FeSe/matrices.npz');analytic=-.1164916715;numeric=-1.0491512951;refcor=-1.0022154658
    R=z['response']+z['norm_response'];factor=14710.507641078/(3*p.vol)
    total=aff+kv;value=float(np.einsum('kij,kji->',R,total).real)*factor
    reference=Density(ROOT/case/'reference_q.density',p.structure)
    data={'case':case,'volume_Bohr3':p.vol,'correction_GPa':value,'correction_from_relaxed_reference_GPa':refcor,'difference_from_reference_correction_GPa':value-refcor,'baseline_GPa':analytic,'corrected_GPa':analytic+value,'numeric_GPa':numeric,'corrected_minus_numeric_GPa':analytic+value-numeric,'reference_density_error_norm':float(np.linalg.norm(q.copy().add(reference,-1).vector())),'convention':'hydrostatic mean tensile stress, already divided by 3','fixed_G':True}
    np.savez_compressed(workspace/'projector_correction_matrices.npz',potential=total,response=R,volume=p.vol)
    (workspace/'stress_correction.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data,indent=2),flush=True)
if __name__=='__main__':run(sys.argv[1],sys.argv[2])
