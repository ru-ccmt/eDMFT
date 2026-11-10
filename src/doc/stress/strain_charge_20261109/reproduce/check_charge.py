"""Analytic native material-charge source versus archived native endpoints."""
from pathlib import Path
import json,sys,time
import numpy as np
root=Path(__file__).resolve().parent
if root.name=='reproduce':root=root.parent
sys.path[:0]=[str(root/'stage'),'/Users/haule/projects/eDMFT_new/src/python',
             '/private/tmp/edmft_screened_native_20261108']
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_strain_charge import hydrostatic_charge_source
from lapw_strain_field import diagonal_strain_potential
from wstruct import Struct
from lapw_strain_potential import unequal_centered_derivative,NativePotentialFileLayout,read_potential_file
from lapw_hydrostatic_potential import hydrostatic_explicit_potential
from wien_response_density import Density
from native_reference import NativeReference


def run(case,mesh='three_k',solve=False):
 native=NativeReference(case,mesh);source=native.source
 work=root/mesh/case;work.mkdir(parents=True,exist_ok=True)
 with load_native_valence_state(source,case,complex_vectors=(case=='SrVO3'),energy_policy='fixed',
      matrix_quadrature='wien2k24',library='/Users/haule/dbin/liblapw_response.so') as v:
  print(case,mesh,'prepared',v.preparation_seconds,flush=True)
  pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
  core=Density(native.xcroot/case/'hxc_core/0'/(case+'.clmcor'),source/(case+'.struct'))
  state=NativeScreenedState(v,pc,core,target_charge=native.target,fixed_number=True)
  folders=[native.work/f'source_{label}_0.000125' for label in ['m','p']]
  cells=[]
  for folder in folders:
   s=Struct();s.ReadStruct(str(folder/(case+'.struct')));cells.append(np.array([s.a,s.b,s.c]))
  offsets=[json.loads((p/'done.json').read_text())['actual_scale']-1 for p in folders]
  strain=np.diag((cells[1]-cells[0])/(native.cell*(offsets[1]-offsets[0])))
  start=time.monotonic();potential,potential_parts=diagonal_strain_potential(pc,strain)
  analytic,parts,report=hydrostatic_charge_source(state,source,case,strain=strain,potential_source=potential)
  report.update(seconds=time.monotonic()-start,metadata=state.metadata,mesh=mesh,steps={})
  h0,_=hydrostatic_explicit_potential(pc);h1,_=diagonal_strain_potential(pc,np.eye(3))
  report['hydrostatic_field_relative']=float(np.linalg.norm(h0-h1)/np.linalg.norm(h0))
  layout=NativePotentialFileLayout(pc.hartree,read_potential_file(source/(case+'.vsp'),pc.structure.jrj),
       read_potential_file(source/(case+'.vns'),pc.structure.jrj),complex_interstitial=(case=='SrVO3'))
  def read(folder,suffix):
   return pc.input_density.zeros().add(Density(folder/(case+suffix),folder/(case+'.struct'))).vector()
  central=native.work/'central'
  centers={'.clmval':read(central,'.clmval'),'.clmcor':read(central,'.clmcor'),
           'normalized.npy':np.load(central/'normalized.npy')}
  nrad=sum(map(len,pc.input_density.rad.values()))
  np.savez_compressed(work/'analytic.npz',source=analytic,**parts)
  for step in [.00025,.000125]:
   folders=[native.work/f'source_{label}_{step}' for label in ['m','p']]
   offsets=[json.loads((p/'done.json').read_text())['actual_scale']-1 for p in folders]
   rows={};derivatives={}
   potentials=[layout.read(p/(case+'.vsp'),p/(case+'.vns')) for p in folders]
   numerical_potential=unequal_centered_derivative(potentials[0],layout.central,potentials[1],*offsets)
   selected_potential=layout.selector.matvec(potential)
   for suffix in centers:
    ends=[np.load(p/suffix) if suffix.endswith('.npy') else read(p,suffix) for p in folders]
    derivatives[suffix]=unequal_centered_derivative(ends[0],centers[suffix],ends[1],*offsets)
   pairs={'valence':(parts['valence_electronic']+parts['valence_radial']+parts['valence_geometry'],derivatives['.clmval']),
          'core':(parts['core'],derivatives['.clmcor']),
          'raw':(parts['raw_output'],derivatives['.clmval']+derivatives['.clmcor']),
          'source':(analytic,derivatives['normalized.npy'])}
   for key,(a,b) in pairs.items():
    rows[key]={}
    for sector,sl in [('all',slice(None)),('MT',slice(0,nrad)),('PW',slice(nrad,None))]:
     rows[key][sector]=dict(absolute=float(np.linalg.norm((a-b)[sl])),relative=float(np.linalg.norm((a-b)[sl])/max(np.linalg.norm(b[sl]),1e-30)))
   report['steps'][str(step)]=dict(actual_offsets=offsets,comparisons=rows,
       potential_relative=float(np.linalg.norm(selected_potential-numerical_potential)/np.linalg.norm(numerical_potential)))
   np.savez_compressed(work/f'native_{step}.npz',**derivatives)
   print(case,mesh,step,json.dumps(rows),flush=True)
  a,b=[np.load(work/f'native_{h}.npz') for h in [.00025,.000125]]
  steps=[report['steps'][str(h)]['actual_offsets'] for h in [.00025,.000125]]
  ratio=((steps[0][1]-steps[0][0])/(steps[1][1]-steps[1][0]))**2
  extrap=(ratio*b['normalized.npy']-a['normalized.npy'])/(ratio-1)
  report['richardson_source_relative']=float(np.linalg.norm(analytic-extrap)/np.linalg.norm(extrap))
  report['richardson_PW_relative']=float(np.linalg.norm((analytic-extrap)[nrad:])/np.linalg.norm(extrap[nrad:]))
  report['normalization_geometry_norm']=float(np.linalg.norm(parts['normalization_geometry']))
  report['post_normalization_change_norm']=float(np.linalg.norm(parts['post_normalization_change']))
  if solve:
   qold=np.load(native.work/'analytic_density.npy')
   result=state.solve_source(analytic,initial=qold,rtol=2e-8,atol=1e-9,restart=30,max_cycles=8)
   report['screened']=dict(result.report,relative_to_numerical_source_solution=float(np.linalg.norm(result.density-qold)/np.linalg.norm(qold)))
   np.save(work/'analytic_density.npy',result.density)
   np.save(work/'analytic_induced_potential.npy',result.induced_potential)
   print(case,mesh,'SCREENED',json.dumps(report['screened']),flush=True)
  (work/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
  np.save(work/'analytic_source.npy',analytic)
  print(case,mesh,'DONE',report['richardson_source_relative'],flush=True)

if __name__=='__main__':run(sys.argv[1],next((a for a in sys.argv[2:] if not a.startswith('--')),'three_k'),solve='--solve' in sys.argv)
