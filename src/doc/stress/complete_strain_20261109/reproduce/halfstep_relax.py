"""Independent frozen-Sigma charge relaxation, seeded by the analytic tangent."""
from pathlib import Path
import sys,os,json,re,time,shutil,subprocess
import numpy as np
root=Path(__file__).resolve().parent;sys.path.insert(0,'/Users/haule/dbin')
from lapw_potential_core_state import load_potential_core_state
from wien_frozen_sigma_response import vector_signature
from wien_response_density import Density
case,label=sys.argv[1:3];h={'minus':-.0005,'plus':.0005,'center':0.0}[label]
out=root/case/'halfstep'/label;out.mkdir(parents=True,exist_ok=True)
source=root/case/'native_live/center';pc=load_potential_core_state('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz')
s=pc.structure;old=Path('/private/tmp/edmft_xc_core_response_20261108');shape0=np.array([s.a,s.b,s.c]);lines=(source/(case+'.struct')).read_text().splitlines();lines[3]=''.join(f'{x*(1+h):10.6f}' for x in shape0)+lines[3][30:];newcell=np.array([float(lines[3][i:i+10]) for i in [0,10,20]]);scale=float(np.mean(newcell/shape0))
if not (out/'prepared.json').exists():
 names=[case+'.'+z for z in ['in0','in1','in1c','in2','in2c','inm','inc','indmfl','indmfi','indmf2','klist','kgen','rotlm']]+['EF.dat','Eorb_imp.dat','projectorw.dat','lapw0.def','lapw1.def','dmft2.def','lcore.def']+[p.name for p in source.glob('sig.inp*')]
 for name in names:
  p=source/name
  if p.exists():shutil.copy2(p,out/name)
 (out/(case+'.struct')).write_text('\n'.join(lines)+'\n')
 q=np.load(root/case/'searched_readout/solution.npz')['density'];rho=pc.input_density.from_vector(pc.input_density.vector()+(scale-1)*q).transport(scale)
 template=Density(source/(case+'.clmsum'),source/(case+'.struct'))
 if set(template.rad)!=set(rho.rad) or set(template.pw)!=set(rho.pw):raise ValueError('native/archive density layouts differ')
 template.rad={k:rho.rad[k].copy() for k in template.rad};template.pw={k:rho.pw[k] for k in template.pw};template.write(out/(case+'.clmsum'))
 in1=out/(case+('.in1c' if case=='SrVO3' else '.in1'));text=in1.read_text().splitlines();v=text[1].split();v[0]=f'{float(v[0])/scale:.8f}';text[1]=' '.join(v);in1.write_text('\n'.join(text)+'\n')
 in2=out/(case+('.in2c' if case=='SrVO3' else '.in2'));text=in2.read_text().splitlines();text[0]='TOT'
 for i,line in enumerate(text):
  if 'GMAX' in line:text[i]=f'{float(line.split()[0])/scale:.8f} GMAX'
 in2.write_text('\n'.join(text)+'\n')
 shutil.copy2('/private/tmp/edmft_screened_native_20261108/full_bz/'+case+'/central/lapw1.fixgvec',out/'lapw1.fixgvec');(out/'lapw1.fixgvec.read').touch()
 (out/'prepared.json').write_text(json.dumps(dict(nominal=h,actual=scale-1,cell=newcell.tolist()),indent=2)+'\n')
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))};env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',WIEN_STRESS_DEBUG='1',WIENROOT='/Users/haule/wien2k_24',WIEN_DMFT_ROOT='/Users/haule/dbin',PYTHONDONTWRITEBYTECODE='1')
weights=np.array([s.mult[a]*st.occ for a,atom in pc.core.atoms.items() for st in atom.states]);shift=sum(s.mult[a]*sum(st.occ for st in atom.states)*pc.metadata['core'][a]['shift_Ry'] for a,atom in pc.core.atoms.items())
history=json.loads((out/'history.json').read_text()) if (out/'history.json').exists() else [];stable=0
for it in range(len(history)+1,41):
 start=time.perf_counter();(out/'core_native.bin').unlink(missing_ok=True)
 for name in ['SOlapm.dat','SOlapm.datup','SOlapm.datdn']:
  if (out/name).exists():raise ValueError('unexpected cached projector normalization')
 commands=[('lapw0',[str(old/'native/SRC_lapw0/lapw0'),'lapw0.def']),('lapw1',[str(Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1')),'lapw1.def']),('dmft2',['/Users/haule/dbin/dmft2_response','dmft2.def']),('lcore',[str(old/'native/SRC_lcore/lcore'),'lcore.def']),('merge',['/Users/haule/miniconda3-clean/envs/work311/bin/python','/Users/haule/dbin/merge_scf.py','--replace']),('mixer',['/Users/haule/wien2k_24/x','mixer','-f',case])]
 for program,cmd in commands:
  with (out/(program+'.log')).open('w') as f:subprocess.run(cmd,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
  if program in ['lapw1','dmft2'] and program.upper()+' END' not in (out/(program+'.log')).read_text():raise RuntimeError(str(out/(program+'.log')))
 den=float(re.findall(r'^:DBG_EN_DEN_TOTAL:\s*(\S+)',(out/(case+'.scf0')).read_text(),re.M)[-1]);vsum=float(re.findall(r'^:VSUM.*?=\s*(\S+)',(out/(case+'.scf2')).read_text(),re.M)[-1]);es=[]
 with (out/'core_native.bin').open('rb') as f:
  while True:
   head=np.fromfile(f,'i4',7)
   if not len(head):break
   es.append(float(np.fromfile(f,'f8',4)[0]));f.seek(int(4*head[3]*8),1)
 K=den+float(2*weights@np.array(es)+shift)+vsum
 dis=float(re.findall(r'^:DIS.*?([0-9.]+)\s*$',(out/(case+'.scfm')).read_text(),re.M)[-1]);delta=None if not history else K-history[-1]['K_Ry'];history.append(dict(iteration=it,K_Ry=K,DEN_Ry=den,VSUM_Ry=vsum,DIS=dis,delta_K_Ry=delta,seconds=time.perf_counter()-start))
 (out/'history.json').write_text(json.dumps(history,indent=2)+'\n');print(case,label,history[-1],flush=True)
 stable=stable+1 if dis<=1e-7 and delta is not None and abs(delta)<=2e-9 else 0
 if stable>=2:break
window=tuple(map(int,(out/(case+'.indmfl')).read_text().splitlines()[0].split()[:2]));same=vector_signature(out/(case+'.vector'),s.nat,window)==vector_signature(root/case/'current'/(case+'.vector'),s.nat,window)
report=dict(case=case,converged=stable>=2,same_vector_signature=same,nominal=h,actual_scale=scale,K_Ry=history[-1]['K_Ry'],iterations=len(history))
(out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
