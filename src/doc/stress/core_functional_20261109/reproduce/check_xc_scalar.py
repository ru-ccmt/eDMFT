"""Compare the analytic XC scalar/gradient with native energy differences.

PBE minus LDA at identical density cancels all electrostatic contributions.
This independent native scalar check tests the actual angular/FFT energy
assembly, rather than assuming integration by parts or Vxc=dExc/dn.
"""
from pathlib import Path
from types import SimpleNamespace
import os,sys,json,shutil,subprocess,re,time
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python']
from lapw_functional_cancellation import PreparedXCEnergy
from lapw_potential_core_state import load_potential_core_state
from lapw_xc_response import MuffinTinXCResponse,InterstitialXCResponse,PeriodicXCResponse
old=Path('/private/tmp/edmft_xc_core_response_20261108')
def xc_state(pc,density,functional):
    mt=[]
    for a,source in enumerate(pc.xc.mt):
        center=np.column_stack([density.rad[k] for k in pc.hartree.keys if k[0]==a])
        extra={} if functional=='LDA' else dict(sine=source.sine.ravel(),cosine=source.cosine.ravel(),angular_derivatives=source.angular[1:])
        mt.append(MuffinTinXCResponse(source.r,source.ells,center,source.Y,source.weights,
            functional=functional,harmonics=source.harmonics,**extra))
    ir=InterstitialXCResponse(pc.hartree.reanalysis,np.array(list(density.pw.values())),
        [pc.structure.a,pc.structure.b,pc.structure.c],functional=functional,rho_min=pc.metadata['rho_min'])
    return SimpleNamespace(hartree=pc.hartree,structure=pc.structure,input_density=density,
        metadata=pc.metadata|dict(functional=functional),xc=PeriodicXCResponse(pc.hartree,mt,ir))
for case in ['SrVO3','FeSe']:
    pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
    e=PreparedXCEnergy(pc)
    # KXC differs from TOT only by requesting the separate E_xc scalar.
    src=old/case/'xc/0';dst=root/case/'xc_energy_central';dst.mkdir(exist_ok=True)
    for name in ['lapw0.def','xc_funcs.h',case+'.struct',case+'.in0',case+'.inm',case+'.clmsum']:
        shutil.copy2(src/name,dst/name)
    lines=(dst/(case+'.in0')).read_text().splitlines();lines[0]='KXC '+lines[0][4:]
    (dst/(case+'.in0')).write_text('\n'.join(lines)+'\n')
    env=os.environ.copy();env.update(WIEN_STRESS_DEBUG='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    with (dst/'run.log').open('w') as log:
        subprocess.run([str(old/'native/SRC_lapw0/lapw0'),'lapw0.def'],cwd=dst,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    native=float(re.findall(r':EXC.*?=\s*(\S+)',(dst/(case+'.scf0')).read_text())[-1])
    assert abs(native-e.value_ry)<5.1e-9,(case,native,e.value_ry)
    report=dict(central_XC_energy_Ry=e.value_ry,native_XC_energy_Ry=native,
        central_error_Ry=e.value_ry-native,directions={})
    lda=PreparedXCEnergy(xc_state(pc,pc.input_density,'LDA')) if case=='FeSe' else e
    ref=json.loads((root/case/'density_energy_audit.json').read_text())
    for label,d in [('probe',np.load('/private/tmp/edmft_hartree_periodic_20261108/'+case+'_direction.npz')['density']),('screened_full_bz',np.load('/private/tmp/edmft_strain_charge_20261109/full_bz/'+case+'/analytic_density.npy'))]:
        defect=float(e.variational_defect@d)
        diff=float((e.density_term_gradient-lda.density_term_gradient)@d)
        native_diff=ref[label]['sectors']['Hxc']['richardson_Ry']-ref[label]['sectors']['H+LDA']['richardson_Ry']
        entry=dict(analytic_XC_defect_Ry=defect,
            remaining_native_DEN_defect_Ry=ref[label]['sectors']['Hxc']['defect_Ry']-defect,
            analytic_XC_minus_LDA_density_term_Ry=diff,
            native_XC_minus_LDA_density_term_Ry=native_diff,
            XC_minus_LDA_error_Ry=diff-native_diff)
        assert abs(diff-native_diff)<2e-7,(case,label,entry)
        report['directions'][label]=entry
    np.savez_compressed(root/case/'xc_energy_covectors.npz',energy=e.energy_gradient,defect=e.variational_defect,density_term=e.density_term_gradient)
    (root/case/'xc_scalar.json').write_text(json.dumps(report,indent=2)+'\n')
    print(case,json.dumps(report),flush=True)
