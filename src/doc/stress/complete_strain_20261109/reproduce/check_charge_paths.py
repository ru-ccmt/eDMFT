from pathlib import Path
import sys,json,numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_potential_core_state import load_potential_core_state
from lapw_functional_cancellation import CoreEigenvalueFunctional
from lapw_density_functional import PreparedDensityFunctional
from wien_response_density import Density
from wstruct import Struct
for case,base,dirs in [('SrVO3',Path('/Users/haule/projects/materials/SrVO3'),['V_strain_m001','V_plus_force','V_strain_p001']),('FeSe',Path('/Users/haule/projects/materials/FeSe_again/FeSe_stress'),['FeSe_RELA_strain_m001','FeSe_RELA','FeSe_RELA_strain_p001'])]:
 pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
 C=PreparedDensityFunctional(pc);core=CoreEigenvalueFunctional.from_state(pc)
 prev=Path('/private/tmp/edmft_strain_charge_20261109/full_bz')/case
 q=np.load(prev/'analytic_density.npy');b=np.load(prev/'analytic.npz');strain=np.array(json.loads((prev/'REPORT.json').read_text())['strain'])
 cv=np.load(Path('/private/tmp/edmft_functional_readout_20261109/full_bz')/case/'covectors.npz')
 g=cv['lattice_charge']+pc.potential_operator.rmatvec(core.potential_gradient)+C.density_gradient
 explicit=float(cv['lattice_joint']@b['joint']+core.action(b['potential_source'])+C.explicit(strain)[0])
 cells=[];ns=[]
 for directory in dirs:
  path=base/directory;s=Struct();s.ReadStruct(str(path/(case+'.struct')))
  scale=float(np.mean(np.array([s.a,s.b,s.c])/np.array([pc.structure.a,pc.structure.b,pc.structure.c])))
  n=Density(path/(case+'.clmsum'),path/(case+'.struct')).transport(scale,True)
  ns.append(pc.input_density.zeros().add(n).vector());cells.append(scale)
 qfd=(ns[2]-ns[0])/(cells[2]-cells[0]);factor=14710.5077/(3*pc.hartree.volume)
 report=dict(case=case,explicit_Ry=explicit,response_norm=float(np.linalg.norm(q)),native_center_vs_live_input=float(np.linalg.norm(pc.input_density.vector()-ns[1])),fd_response_norm=float(np.linalg.norm(qfd)),response_relative_difference=float(np.linalg.norm(qfd-q)/np.linalg.norm(qfd)),mean_GPa_response=(explicit+float(g@q))*factor,mean_GPa_live_path=(explicit+float(g@qfd))*factor,readout_change_GPa=float(g@(qfd-q))*factor)
 print(json.dumps(report),flush=True)
 np.savez_compressed(root/case/'charge_paths.npz',q=q,q_live=qfd,gradient=g)
 (root/case/'charge_paths.json').write_text(json.dumps(report,indent=2)+'\n')
