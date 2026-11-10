"""Audit artificial central-coordinate charge projection after inverse transport."""
from pathlib import Path
import sys,json,numpy as np
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_potential_core_state import load_potential_core_state
from lapw_charge_response import NativeChargeNormalization,native_mixer_inversion
from lapw_strain_charge import charge_integration_geometry
from wien_response_operators import hydrostatic_transport_source
from wien_response_density import Density
from stress_response_adjoint import RY_BOHR3_TO_GPA
report={}
for case in ['FeSe','SrVO3']:
 folder=root/case/('refined' if case=='FeSe' else 'current');run=root/case/('refined_total' if case=='FeSe' else 'tensor')
 archive=folder/'potential_core.npz'
 if not archive.exists():archive=Path('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108')/(case+'_potential_core.npz')
 pc=load_potential_core_state(archive);s=pc.structure;inv=native_mixer_inversion(s)
 raw=Density(folder/(case+'.clmval'),folder/(case+'.struct'));corepath=folder/(case+'.clmcor');corepath=corepath if corepath.exists() else Path('/private/tmp/edmft_xc_core_response_20261108')/case/'hxc_core/0'/(case+'.clmcor');raw.add(Density(corepath,folder/(case+'.struct')))
 norm=NativeChargeNormalization(raw,s,float(np.dot(s.mult,s.Znuc)),inversion=inv)
 output=raw.from_vector(norm.central);post=NativeChargeNormalization(output,s,norm.target,inversion=inv)
 reportfile=run/('report.json' if case=='FeSe' else 'REPORT.json');r=json.loads(reportfile.read_text());sol=np.load(run/('adjoint.npz' if case=='FeSe' else 'solution.npz'));print(case,sol.files,flush=True)
 lam=sol['adjoint'];corrections=[];charges=[];deltas=[]
 for e in [np.eye(3),np.diag([0.,0.,1.])]:
  tau=np.trace(e)/3*hydrostatic_transport_source(output).vector()
  ig=charge_integration_geometry(output,s,e,inversion=inv)
  charge=-float(ig@norm.central+(post.mt+post.ir)@tau)
  delta=post.p0/post.I*charge
  corrections.append(float(lam@delta));charges.append(charge);deltas.append(delta)
 old=r.get('mean_tensile_stress_GPa',r.get('mean_GPa'))
 if old is None:old=np.trace(r['stress_GPa'])/3
 mean=old+corrections[0]*RY_BOHR3_TO_GPA/(3*s.a*s.b*s.c)
 report[case]=dict(omitted_source_charges=charges,derivative_corrections_Ry=corrections,old_mean_GPa=old,corrected_mean_GPa=mean)
 np.savez_compressed(root/case/'coordinate_source.npz',deltas=np.array(deltas))
 print(report[case],flush=True)
(root/'coordinate_charge.json').write_text(json.dumps(report,indent=2)+'\n')
