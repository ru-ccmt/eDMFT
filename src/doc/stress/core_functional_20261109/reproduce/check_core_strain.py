"""Native full-mesh explicit hydrostatic core scalar, fixed RMT and inc shifts."""
from pathlib import Path
import json,sys
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_functional_cancellation import CoreEigenvalueFunctional
from lapw_potential_core_state import load_potential_core_state
from lapw_strain_potential import unequal_centered_derivative

def energies(path):
    result=[]
    with path.open('rb') as f:
        while True:
            header=np.fromfile(f,'i4',7)
            if not len(header):break
            result.append(float(np.fromfile(f,'f8',4)[0]))
            f.seek(int(4*header[3]*8),1)
    return np.array(result)
for case in ['SrVO3','FeSe']:
    pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    core=CoreEigenvalueFunctional.from_state(pc)
    base=Path('/private/tmp/edmft_screened_native_20261108/full_bz')/case
    source=np.load('/private/tmp/edmft_strain_charge_20261109/full_bz/'+case+'/analytic.npz')['potential_source']
    analytic=core.action(source)
    weights=np.array([pc.structure.mult[a]*st.occ for a,atom in pc.core.atoms.items() for st in atom.states])
    center=2*weights@energies(base/'central/core_native.bin')
    rows=[]
    for h in [.00025,.000125]:
        dirs=[base/f'source_{sign}_{h}' for sign in ['m','p']]
        offsets=[json.loads((p/'done.json').read_text())['actual_scale']-1 for p in dirs]
        vals=[2*weights@energies(p/'core_native.bin') for p in dirs]
        fd=float(unequal_centered_derivative(vals[0],center,vals[1],*offsets))
        rows.append(dict(step=h,offsets=offsets,fd_Ry=fd,error_Ry=fd-analytic))
    ratio=((rows[0]['offsets'][1]-rows[0]['offsets'][0])/(rows[1]['offsets'][1]-rows[1]['offsets'][0]))**2
    rich=(ratio*rows[1]['fd_Ry']-rows[0]['fd_Ry'])/(ratio-1)
    report=dict(analytic_Ry=analytic,steps=rows,richardson_Ry=rich,richardson_error_Ry=rich-analytic,
        scope='explicit hydrostatic source, not screened total stress',stress_applied=False)
    (root/case/'core_strain.json').write_text(json.dumps(report,indent=2)+'\n')
    print(case,json.dumps(report),flush=True)
