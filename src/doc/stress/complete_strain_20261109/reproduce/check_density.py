from pathlib import Path
import sys,json,re,time
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_potential_core_state import load_potential_core_state
from lapw_density_functional import PreparedDensityFunctional
from lapw_strain_potential import unequal_centered_derivative
from wstruct import Struct

def value(path):
    text=path.read_text();v=re.findall(r':DBG_EN_DEN_TOTAL:\s*(\S+)',text)
    return float(v[-1]) if v else float(re.findall(r':DEN.*?=\s*(\S+)',text)[-1])

for case in sys.argv[1:] or ['SrVO3','FeSe']:
    t=time.perf_counter()
    pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    c=PreparedDensityFunctional(pc)
    native=Path('/private/tmp/edmft_core_functional_20261109')/case
    center=value(native/'native_energy/0'/(case+'.scf0'))
    refs=json.loads((native/'density_energy_audit.json').read_text())
    report=dict(value_Ry=c.value_ry,native_value_Ry=center,central_error_Ry=c.value_ry-center,preparation_seconds=time.perf_counter()-t,probes={})
    for label,d in [('probe',np.load('/private/tmp/edmft_hartree_periodic_20261108/'+case+'_direction.npz')['density']),('screened_full_bz',np.load('/private/tmp/edmft_strain_charge_20261109/full_bz/'+case+'/analytic_density.npy'))]:
        analytic=c.action(d);numeric=refs[label]['sectors']['Hxc']['richardson_Ry']
        report['probes'][label]=dict(analytic_Ry=analytic,numeric_Ry=numeric,error_Ry=analytic-numeric)
    base=Path('/private/tmp/edmft_screened_native_20261108/full_bz')/case
    folders=[base/f'source_{sign}_0.000125' for sign in ['m','p']]
    cells=[]
    for folder in folders:
        st=Struct();st.ReadStruct(str(folder/(case+'.struct')));cells.append(np.array([st.a,st.b,st.c]))
    offsets=[json.loads((p/'done.json').read_text())['actual_scale']-1 for p in folders]
    strain=np.diag((cells[1]-cells[0])/(np.array([pc.structure.a,pc.structure.b,pc.structure.c])*(offsets[1]-offsets[0])))
    explicit,parts=c.explicit(strain)
    report['explicit']=dict(analytic_Ry=explicit,parts=parts,steps=[])
    base=Path('/private/tmp/edmft_screened_native_20261108/full_bz')/case
    for h in [.00025,.000125]:
        dirs=[base/f'source_{sign}_{h}' for sign in ['m','p']]
        off=[json.loads((d/'done.json').read_text())['actual_scale']-1 for d in dirs]
        vals=[value(d/(case+'.scf0')) for d in dirs]
        fd=float(unequal_centered_derivative(vals[0],center,vals[1],*off))
        report['explicit']['steps'].append(dict(step=h,offsets=off,fd_Ry=fd,error_Ry=fd-explicit))
    rows=report['explicit']['steps'];ratio=((rows[0]['offsets'][1]-rows[0]['offsets'][0])/(rows[1]['offsets'][1]-rows[1]['offsets'][0]))**2
    fd=(ratio*rows[1]['fd_Ry']-rows[0]['fd_Ry'])/(ratio-1)
    report['explicit'].update(richardson_Ry=fd,richardson_error_Ry=fd-explicit)
    dest=root/case;dest.mkdir(exist_ok=True)
    np.savez_compressed(dest/'density_covectors.npz',gradient=c.density_gradient,defect=c.variational_defect)
    (dest/'density.json').write_text(json.dumps(report,indent=2)+'\n')
    print(case,json.dumps(report),flush=True)
