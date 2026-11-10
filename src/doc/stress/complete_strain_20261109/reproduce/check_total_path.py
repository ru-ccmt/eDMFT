"""Independent native full-functional directional test along the solved response."""
from pathlib import Path
import sys, json, shutil, re, time
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=['/private/tmp/edmft_screened_native_20261108','/private/tmp/edmft_functional_readout_20261109']
from native_reference import NativeReference
from check_material import native_scalar
sys.path.insert(0,str(root/'stage'))
from lapw_potential_core_state import load_potential_core_state
from lapw_functional_cancellation import CoreEigenvalueFunctional
from lapw_density_functional import PreparedDensityFunctional
from lapw_strain_potential import unequal_centered_derivative

def energy(path):
 vals=[]
 with path.open('rb') as f:
  while True:
   h=np.fromfile(f,'i4',7)
   if not len(h):break
   assert len(h)==7
   vals.append(float(np.fromfile(f,'f8',4)[0]));f.seek(int(4*h[3]*8),1)
 return np.array(vals)
def den(folder,case):
 text=(folder/(case+'.scf0')).read_text();m=re.findall(r':DBG_EN_DEN_TOTAL:\s*(\S+)',text)
 if not m:raise ValueError('high precision DEN absent')
 return float(m[-1])

case=sys.argv[1];hlist=[float(x) for x in sys.argv[2:]] or [.00025,.000125]
r=NativeReference(case);old=r.work;r.work=root/case/'total_path';r.work.mkdir(parents=True,exist_ok=True)
(r.work/'central').mkdir(exist_ok=True)
shutil.copy2(old/'central/lapw1.fixgvec',r.work/'central/lapw1.fixgvec')
r.env['WIEN_STRESS_DEBUG']='1'
pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
C=PreparedDensityFunctional(pc);core=CoreEigenvalueFunctional.from_state(pc)
weights=np.array([pc.structure.mult[a]*st.occ for a,atom in pc.core.atoms.items() for st in atom.states])
shift=sum(pc.structure.mult[a]*sum(st.occ for st in atom.states)*pc.metadata['core'][a]['shift_Ry'] for a,atom in pc.core.atoms.items())
prev=Path('/private/tmp/edmft_strain_charge_20261109/full_bz')/case
q=np.load(prev/'analytic_density.npy');parts=np.load(prev/'analytic.npz');strain=np.array(json.loads((prev/'REPORT.json').read_text())['strain'])
cv=np.load(Path('/private/tmp/edmft_functional_readout_20261109/full_bz')/case/'covectors.npz')
aL=float(cv['lattice_joint']@parts['joint']+cv['lattice_charge']@q)
aC=float(core.action(parts['potential_source']+pc.potential_operator.matvec(q)))
aD=float(C.explicit(strain)[0]+C.action(q))
analytic=dict(lattice=aL,core=aC,DEN=aD);analytic['total']=sum(analytic.values())
center=dict(lattice=native_scalar(r.source,r.base),core=core.value_ry,DEN=C.value_ry);center['total']=sum(center.values())
report=dict(case=case,mesh='full_bz',analytic_Ry=analytic,scope='directional native test along analytic q, not separately converged pressure',steps=[],timings=[])
for h in hlist:
 vals=[];offsets=[]
 for sign in [-1,1]:
  label=('m_' if sign<0 else 'p_')+str(h)
  actual=float(np.mean(np.array([float(f'{a*(1+sign*h):.6f}') for a in r.cell])/r.cell))
  rho=r.rho.from_vector(r.rho.vector()+(actual-1)*q)
  _,scale=r.evaluate(label,density=rho,scale=1+sign*h)
  folder=r.work/label;meta=json.loads((folder/'done.json').read_text())
  v=dict(lattice=native_scalar(folder,r.base,mu=meta['mu_Ry'],number=r.target_band),core=float(2*weights@energy(folder/'core_native.bin')+shift),DEN=den(folder,case))
  v['total']=sum(v.values());vals.append(v);offsets.append(scale-1)
 row=dict(h=h,offsets=offsets,values_Ry=vals)
 row['fd_Ry']={k:float(unequal_centered_derivative(vals[0][k],center[k],vals[1][k],*offsets)) for k in center}
 row['error_Ry']={k:row['fd_Ry'][k]-analytic[k] for k in center}
 report['steps'].append(row);report['timings']=r.timings
 if len(report['steps'])>1:
  a,b=report['steps'][-2:];ratio=((a['offsets'][1]-a['offsets'][0])/(b['offsets'][1]-b['offsets'][0]))**2
  rich={k:(ratio*b['fd_Ry'][k]-a['fd_Ry'][k])/(ratio-1) for k in center}
  report['richardson_Ry']=rich;report['richardson_error_Ry']={k:rich[k]-analytic[k] for k in center}
 (r.work/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
 print(case,row,flush=True)
print(json.dumps(report,indent=2),flush=True)
