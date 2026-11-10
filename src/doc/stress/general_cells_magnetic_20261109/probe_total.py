from pathlib import Path
from contextlib import ExitStack
import sys,os,time,json,traceback,numpy as np
root=Path(__file__).parent;sys.path.insert(0,str(root/'stage'))
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_magnetic_response import combine_valence
from lapw_klein_functional import PreparedKleinFunctional
from lapw_energy_response import SearchedEnergyResponse
from lapw_strain_charge import hydrostatic_charge_source
from lapw_strain_field import diagonal_strain_potential
from wien_response_density import Density
results={};only=os.environ.get('TEST_MATERIAL','MnO')
for mat,d,c in [('MnO','MnO_center','MnO'),('MnTe','dmft','0_800_MnTe')]:
 if mat!=only:continue
 results[mat]=[]
 for cutoff in [1e-10,1e-4]:
  t=time.perf_counter();w=root/'native'/mat/d
  try:
   pc=load_potential_core_state(w/'potential_core.npz');core=Density(w/(c+'.clmcor'),w/(c+'.struct'))
   with ExitStack() as stack:
    branches=[stack.enter_context(load_native_valence_state(w,c,complex_vectors=False,energy_policy='fixed',matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so',diagnostic_gap_cutoff=cutoff,spin_channel=channel,cache_mb=128)) for channel in (['','dn'] if mat=='MnTe' else [''])]
    v=combine_valence(*branches) if len(branches)==2 else branches[0]
    state=NativeScreenedState(v,pc,core,target_charge=np.dot(pc.structure.mult,pc.structure.Znuc),fixed_number=(mat=='MnTe'),normalization_enabled=((w/(c+'.inm')).read_text().split()[2].upper()=='YES'))
    functional=PreparedKleinFunctional(state);energy=SearchedEnergyResponse(state,functional,w,c,complex_vectors=False,library='/Users/haule/dbin/liblapw_response.so')
    print(mat,cutoff,'prepared functional',time.perf_counter()-t,flush=True)
    e=np.eye(3);pot,_=diagonal_strain_potential(pc,e)
    source,parts,report=hydrostatic_charge_source(state,w,c,strain=e,potential_source=pot)
    source=energy.complete_source(source,pot)
    explicit,pieces=functional.explicit(parts['band_geometry'][:,None],pot[:,None],e[None]);explicit+=energy.explicit_addition(pot[:,None])
    vals,lam,solver=energy.screened.contractions(energy.charge_gradient,source[:,None],explicit=explicit,rtol=2e-12,atol=1e-11,restart=30,max_cycles=8)
    row=dict(diagnostic_only=True,cutoff_Ry=cutoff,mean_stress_GPa=float(vals[0]*14710.507641078/(3*pc.hartree.volume)),derivative_Ry=float(vals[0]),explicit_Ry=float(explicit[0]),induced_Ry=float(lam@source),pieces=pieces,metadata=state.metadata,solver=solver)
  except Exception as exc:row=dict(cutoff_Ry=cutoff,error=type(exc).__name__+': '+str(exc));traceback.print_exc()
  row['seconds']=time.perf_counter()-t;results[mat].append(row);print(json.dumps(row),flush=True);(root/('total_probe_'+mat+'.json')).write_text(json.dumps(results,indent=2)+'\n')
  if 'error' in row:break
