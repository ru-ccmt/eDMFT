"""Physical source FD and complete readout without extra coordinate normalization."""
from pathlib import Path
import sys,json,numpy as np
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_potential_core_state import load_potential_core_state
from lapw_charge_response import NativeChargeNormalization,native_mixer_inversion
from lapw_strain_charge import charge_integration_geometry
from wien_response_operators import hydrostatic_transport_source
from wien_response_density import Density
from stress_response_adjoint import RY_BOHR3_TO_GPA
from wstruct import Struct
case='FeSe';pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/FeSe_potential_core.npz');s=pc.structure;inv=native_mixer_inversion(s);target=float(s.mult@s.Znuc)
raw=Density(root/case/'current/FeSe.clmval',root/case/'current/FeSe.struct');raw.add(Density('/private/tmp/edmft_xc_core_response_20261108/FeSe/hxc_core/0/FeSe.clmcor',root/case/'current/FeSe.struct'));n=NativeChargeNormalization(raw,s,target,inversion=inv);output=raw.from_vector(n.central);post=NativeChargeNormalization(output,s,target,inversion=inv);e=np.diag([0.,0.,1.]);tau=hydrostatic_transport_source(output).vector()/3;ig=charge_integration_geometry(output,s,e,inversion=inv);charge=-float(ig@n.central+(post.mt+post.ir)@tau);delta=post.p0/post.I*charge
analytic=np.load(root/case/'tensor/source1.npz')['source']+delta;rows=[]
for h in [.00025,.000125]:
 values=[];steps=[]
 for sign in [-1,1]:
  d=root/case/'zz_source'/(str(h)+('_m' if sign<0 else '_p'));ss=Struct();ss.ReadStruct(str(d/'FeSe.struct'));step=ss.c/s.c-1;steps.append(step);scale=(1+step)**(1/3)
  raw=Density(d/'FeSe.clmval',d/'FeSe.struct');raw.add(Density(d/'FeSe.clmcor',d/'FeSe.struct'));physical=NativeChargeNormalization(raw,ss,target,inversion=inv).central
  values.append(pc.input_density.zeros().add(raw.from_vector(physical)).transport(scale,True).vector())
 fd=(values[1]-values[0])/(steps[1]-steps[0]);rows.append(dict(step=h,relative_error=float(np.linalg.norm(fd-analytic)/np.linalg.norm(analytic)),absolute_error=float(np.linalg.norm(fd-analytic))))
# The physical relaxed density tangent is the inverse transport WITHOUT a
# second central-cell normalization. Contract it with the same scalar gradient.
values=[];steps=[]
for label in ['minus','plus']:
 d=root/case/'halfstep'/label;scale=json.loads((d/'RESULT.json').read_text())['actual_scale'];steps.append(scale-1);rho=Density(d/'FeSe.clmsum',d/'FeSe.struct');values.append(pc.input_density.zeros().add(rho).transport(scale,True).vector())
q=(values[1]-values[0])/(steps[1]-steps[0]);np.save(root/case/'halfstep/physical_density_tangent.npy',q)
r=json.loads((root/case/'refined_total/report.json').read_text());c=np.load(root/case/'refined_total/adjoint.npz')['charge_gradient'];D=float(r['explicit_Ry'][0]+c@q);mean=D*RY_BOHR3_TO_GPA/(3*s.a*s.b*s.c)
report=dict(physical_axial_source=rows,native_density_readout_Ry=D,native_density_readout_mean_GPa=mean,relaxed_numeric_mean_GPa=json.loads((root/case/'halfstep/REPORT.json').read_text())['mean_tensile_stress_GPa']);(root/case/'physical_sources.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
