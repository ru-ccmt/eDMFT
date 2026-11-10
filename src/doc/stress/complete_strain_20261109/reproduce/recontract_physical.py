"""Exact recontraction after removing only the explicit-source projection.

J, the complete scalar covector c, and its adjoint are unchanged. Reuse them;
only s_nu changes. This is not a fit to finite-difference pressure.
"""
from pathlib import Path
import sys,json,numpy as np
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_total_stress import stress_from_contractions
from stress_response_adjoint import RY_BOHR3_TO_GPA
from lapw_potential_core_state import load_potential_core_state
case='FeSe';src=root/case/'optimized_total';dst=root/case/'physical_total';dst.mkdir(exist_ok=True)
r=json.loads((src/'report.json').read_text());sol=np.load(src/'adjoint.npz');delta=np.load(root/case/'coordinate_source.npz')['deltas'];lam=sol['adjoint'];correction=delta@lam;D=np.array(r['derivatives_Ry'])+correction
pc=load_potential_core_state(root/case/'refined/potential_core.npz');s=pc.structure;stress=stress_from_contractions(r['strains'],D,s.a*s.b*s.c)
r.update(stress_Ry_Bohr3=stress.tolist(),stress_GPa=(stress*RY_BOHR3_TO_GPA).tolist(),mean_tensile_stress_GPa=float(np.trace(stress)*RY_BOHR3_TO_GPA/3),compression_pressure_GPa=float(-np.trace(stress)*RY_BOHR3_TO_GPA/3),derivatives_Ry=D.tolist(),induced_Ry=(np.array(r['induced_Ry'])+correction).tolist(),coordinate_source_correction_Ry=correction.tolist(),includes_post_transport_normalization=False,execution_note='Exact source recontraction of the optimized total run: J, c and lambda unchanged. Timings belong to the saved full driver run before this final explicit-source-only correction.')
for i,d in enumerate(delta):
 b=np.load(src/f'source{i}.npz');np.savez_compressed(dst/f'source{i}.npz',strain=b['strain'],source=b['source']+d,geometry=b['geometry'],potential=b['potential'])
np.savez_compressed(dst/'adjoint.npz',**{k:sol[k] for k in sol.files});(dst/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(r['stress_GPa'],r['mean_tensile_stress_GPa'])
