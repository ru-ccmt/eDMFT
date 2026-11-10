"""Diagnostic first-order fixed-point defect, never applied to published stress.

K(n*)-K(n) ~= lambda^T [F(n)-n], (I-J)^T lambda = dK/dn.
Use the CENTRAL adjoint at nearby endpoints in matching native coordinates.
This neglects endpoint variation of lambda and terms quadratic in the defect;
it is a convergence-sensitivity estimate, not a corrected finite difference.
"""
from pathlib import Path
import sys,json,numpy as np
root=Path(__file__).parent;sys.path.insert(0,str(root/'stage'))
from lapw_charge_response import NativeChargeNormalization,native_mixer_inversion
from lapw_lattice import native_structure,NativeLattice
from wien_response_density import Density
from wstruct import Struct
adj=np.load(root/'MnO_stationarity_adjoint.npz')['adjoint']
rows=[];reference=None
for name in ['MnO_center','MnO_m_0.001','MnO_p_0.001']:
 w=root/'native/MnO'/name;c='MnO';s=Struct();s.ReadStruct(str(w/(c+'.struct')));s=native_structure(s)
 inp=Density(w/(c+'.clmsum'),w/(c+'.struct'));raw=Density(w/(c+'.clmval'),w/(c+'.struct'))
 raw.add(Density(w/(c+'.clmcor'),w/(c+'.struct')))
 if reference is None:reference=inp;volume=NativeLattice.from_structure(s).volume
 assert tuple(reference.rad)==tuple(inp.rad) and tuple(reference.pw)==tuple(inp.pw)
 assert all(np.array_equal(a,b) for a,b in zip(reference.grids,inp.grids))
 normalization=NativeChargeNormalization(raw,s,float(np.dot(s.mult,s.Znuc)),inversion=native_mixer_inversion(s))
 residual=normalization.central-inp.vector()
 rows.append(dict(directory=name,residual_norm=float(np.linalg.norm(residual)),linearized_energy_defect_Ry=float(adj@residual)))
shift=(rows[2]['linearized_energy_defect_Ry']-rows[1]['linearized_energy_defect_Ry'])/.002*14710.507641078/(3*volume)
result=dict(scope=__doc__,rows=rows,estimated_change_of_numeric_mean_stress_GPa=shift)
(root/'stationarity_defect.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
