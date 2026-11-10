"""Read-only native core and LAPW0-energy checks on archived endpoints."""
from pathlib import Path
import json, sys, time, re
import numpy as np

root = Path(__file__).resolve().parent
sys.path[:0] = [str(root / "stage"), "/Users/haule/projects/eDMFT_new/src/python"]
from lapw_potential_core_state import load_potential_core_state
from lapw_functional_cancellation import (
    CoreEigenvalueFunctional, SpatialPotentialPairing, stationary_potential_readout,
)
from wien_response_density import Density

old = Path('/private/tmp/edmft_xc_core_response_20261108')
readout = Path('/private/tmp/edmft_functional_readout_20261109')

def core_energies(path):
    energies = []
    with path.open('rb') as f:
        while True:
            h = np.fromfile(f, 'i4', 7)
            if not len(h): break
            assert len(h) == 7
            energies.append(float(np.fromfile(f, 'f8', 4)[0]))
            f.seek(int(4 * h[3] * 8), 1)
    return np.array(energies)

def density_energy(folder, case):
    text=(folder/(case+'.scf0')).read_text()
    debug=re.findall(r'^:DBG_EN_DEN_TOTAL:\s*(\S+)', text, re.M)
    return float(debug[-1]) if debug else float(re.findall(r'^:DEN.*?=\s*(\S+)',text,re.M)[-1])

for case in sys.argv[1:] or ['SrVO3','FeSe']:
    start=time.perf_counter()
    pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    c=CoreEigenvalueFunctional.from_state(pc)
    pair=SpatialPotentialPairing(pc.hartree,pc.structure.mult)
    weights=np.array([pc.structure.mult[a]*st.occ for a,atom in pc.core.atoms.items() for st in atom.states])
    shift=sum(pc.structure.mult[a]*sum(st.occ for st in atom.states)*pc.metadata['core'][a]['shift_Ry'] for a,atom in pc.core.atoms.items())
    energy=lambda p: float(2*weights@core_energies(p)+shift)
    scalar_error=c.value_ry-energy(old/case/'hxc_core/0/core_native.bin')
    assert abs(scalar_error)<2e-10
    direction=np.load('/private/tmp/edmft_hartree_periodic_20261108/'+case+'_direction.npz')['density']
    potential=pc.potential_operator.matvec(direction)
    core_density=pc.input_density.zeros()
    for a,atom in pc.core.atoms.items(): core_density.rad[a,0,0]=atom.density.copy()
    core_spatial=pair.gradient(pair.density_vector(core_density))
    input_spatial=pair.gradient(pair.density_vector(pc.input_density))
    report=dict(case=case,mode=pc.metadata['mode'],core_energy_Ry=c.value_ry,
        central_core_scalar_error_Ry=scalar_error,core_number=c.core_number,steps={},
        core_derivative_Ry=c.action(potential),core_exterior_derivative_Ry=float(c.exterior_gradient@potential),
        core_spatial_derivative_Ry=float(core_spatial@potential),
        core_minus_spatial_Ry=float((c.potential_gradient-core_spatial)@potential),
        stationary_DEN_derivative_Ry=float(-input_spatial@potential))
    for h in [.01,.005,.001,.0005]:
        cp=old/case/'hxc_core'/str(h)/'core_native.bin'
        cm=old/case/'hxc_core'/str(-h)/'core_native.bin'
        xp=old/case/'xc'/str(h);xm=old/case/'xc'/str(-h)
        entry={}
        if cp.exists() and cm.exists():
            numeric=(energy(cp)-energy(cm))/(2*h)
            entry.update(core_fd_Ry=numeric,core_error_Ry=numeric-c.action(potential))
        if (xp/(case+'.scf0')).exists():
            numeric=(density_energy(xp,case)-density_energy(xm,case))/(2*h)
            entry.update(DEN_fd_Ry=numeric,DEN_identity_defect_Ry=numeric+float(input_spatial@potential))
        report['steps'][str(h)]=entry
    target=root/case;target.mkdir(exist_ok=True)
    cv=np.load(readout/'three_k'/case/'covectors.npz')
    components=stationary_potential_readout(pc,cv['lattice_potential'])
    raw=pc.input_density.zeros().add(Density('/private/tmp/edmft_full_k_response_20261108/'+case+'/charge/0/'+case+'.clmval', '/private/tmp/edmft_full_k_response_20261108/'+case+'/charge/0/'+case+'.struct')).add(core_density)
    raw_spatial=pair.gradient(pair.density_vector(raw))
    report['potential_cancellation']={}
    induced=np.load('/private/tmp/edmft_strain_charge_20261109/three_k/'+case+'/analytic_induced_potential.npy')
    for label,v in [('probe',potential),('screened_induced',induced)]:
        report['potential_cancellation'][label]=dict(
            lattice_Ry=float(components['lattice']@v),core_Ry=c.action(v),
            negative_input_pairing_Ry=float(-input_spatial@v),
            stationary_sum_Ry=float(components['potential']@v),
            output_minus_input_pairing_Ry=float((raw_spatial-input_spatial)@v),
            core_minus_spatial_Ry=float((c.potential_gradient-core_spatial)@v),
            core_exterior_Ry=float(c.exterior_gradient@v))
    report.update(seconds=time.perf_counter()-start,stress_applied=False,full_klein_functional=False)
    np.savez_compressed(target/'covectors.npz',**{k:v for k,v in components.items() if isinstance(v,np.ndarray)},core_spatial=core_spatial,core_exterior=c.exterior_gradient,raw_output_spatial=raw_spatial)
    (target/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
