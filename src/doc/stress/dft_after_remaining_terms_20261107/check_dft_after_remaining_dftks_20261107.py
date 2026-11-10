from pathlib import Path
import subprocess,shutil,os,concurrent.futures,re,json,hashlib,time
root=Path('/private/tmp/dft_after_remaining_dftks_20261107');root.mkdir(exist_ok=True)
wien=Path('/Users/haule/wien2k_24');old=Path('/Users/haule/projects/eDMFT_new/src/doc/stress/remaining_terms_audit_20261107/before')
sets=[('SrVO3','SrVO3',Path('/Users/haule/projects/materials/SrVO3'),['DFT_m001','DFT_center','DFT_p001']),('FeSe','FeSe',Path('/Users/haule/projects/materials/FeSe_again/FeSe_stress'),['DFT_RELA_m001','DFT_RELA','DFT_RELA_p001']),('MnTe','0_800_MnTe',Path('/Users/haule/projects/materials/MnTe'),['MnTe_m_0.001','0_800_MnTe','MnTe_p_0.001'])]
env=dict(os.environ,WIENROOT=str(wien),WIEN_DMFT_ROOT='/Users/haule/dbin',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
for k in list(env):
 if k.startswith(('WIEN_STRESS_','EDMFT_STRESS_')):del env[k]
env['WIEN_STRESS_DEBUG']='1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def tensor(p,tag):
 s=p.read_text();return [[float(x) for x in re.findall('^:'+tag+f'{i:03d}:'+r'\s*(.*)',s,re.M)[-1].split()[:3]] for i in [1,2,3]]
def run(spec):
 label,case,parent,names=spec;source=parent/names[1];d=root/label;d.mkdir(exist_ok=True);before=d/'before';before.mkdir(exist_ok=True);after=d/'after';after.mkdir(exist_ok=True)
 inputs={}
 for suf in ['struct','clmsum','clmup','clmdn','in0','in1','in1c','in2','in2c','inc','inm','klist','kgen','rotlm','vsp','vns','clmcor','r2v']:
  p=source/(case+'.'+suf)
  if p.exists():shutil.copy2(p,d/p.name);inputs[p.name]=sha(p)
 for suf in ['scfm','scf0','scfc','scf2']:
  p=source/(case+'.'+suf)
  if p.exists():shutil.copy2(p,before/('original_'+p.name))
 complex_case=(d/(case+'.in1c')).exists() and (d/(case+'.in1c')).stat().st_size>0
 for step in ['lapw0','lapw1','lapw2','lcore']:
  cmd=[str(wien/'x'),step,'-f',case]
  if complex_case and step in ['lapw1','lapw2']:cmd+=['-c']
  if step=='lapw1':cmd=['/Users/haule/dbin/x_dmft.py','lapw1']
  start=time.monotonic();print('START',label,step,flush=True)
  with (d/(step+'.log')).open('w') as f:subprocess.run(cmd,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
  print('END',label,step,round(time.monotonic()-start,1),flush=True)
  expected=case+'.'+dict(lapw0='scf0',lapw1='scf1',lapw2='scf2',lcore='scfc')[step]
  if not (d/expected).exists():raise RuntimeError((label,step,'no output'))
 for suf in ['scf0','scfc','scf2','scf1','vsp','vns','clmcor','energy']:
  p=d/(case+'.'+suf)
  if p.exists():shutil.copy2(p,after/p.name)
 # Re-evaluate the old density/core algorithms on exactly the same charge.
 for step in ['lapw0','lcore']:
  with (d/(step+'_old.log')).open('w') as f:subprocess.run([str(old/step),step+'.def'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 for suf in ['scf0','scfc','vsp','vns','clmcor']:shutil.copy2(d/(case+'.'+suf),before/(case+'.'+suf))
 import numpy as np
 v=np.array(tensor(after/(case+'.scf2'),'STR_VAL'));res={}
 for name,path in [('before',before),('after',after)]:
  dens=np.array(tensor(path/(case+'.scf0'),'STR_DENS'));core=np.array(tensor(path/(case+'.scfc'),'STR_CORE'));tot=(dens-core-v)*14710.5164057
  res[name]={'diagonal_GPa':np.diag(tot).tolist(),'mean_GPa':float(np.trace(tot)/3),'density_RyBohr3':dens.tolist(),'core_RyBohr3':core.tolist()}
 energies=[];charges=[];eorbs=[];basis=[]
 for folder in names:
  q=parent/folder;s=(q/(case+'.scfm')).read_text();energies.append(float(re.findall(r'^:ENE.*=\s*([-+0-9.]+)',s,re.M)[-1]));charges.append(float(re.findall(r'^:DIS.*\s([0-9.]+)\s*$',s,re.M)[-1]));ep=q/'Eorb_imp.dat';eorbs.append(sha(ep) if ep.exists() else None)
  basis.append({suf:sha(q/(case+'.'+suf)) for suf in ['struct','in0','in1c' if complex_case else 'in1','in2c' if complex_case else 'in2','klist']})
 vol=float(re.findall(r'^:VOL\s*:.*=\s*([-+0-9.]+)',(after/(case+'.scf0')).read_text(),re.M)[-1]);num=(energies[2]-energies[0])/(.006*vol)*14710.5164057
 res.update(case=case,source=str(source),input_hashes=inputs,reference_dirs=[str(parent/n) for n in names],energies_Ry=energies,charge_distances=charges,Eorb_hashes=eorbs,volume=vol,numerical_GPa=num,after_minus_numeric_GPa=res['after']['mean_GPa']-num,ordinary_outputs_identical={s:(before/(case+'.'+s)).read_bytes()==(after/(case+'.'+s)).read_bytes() for s in ['vsp','vns','clmcor']},reference_input_hashes=basis)
 assert all(sha(source/f)==h for f,h in inputs.items())
 (d/'results.json').write_text(json.dumps(res,indent=2)+'\n');print('RESULT',label,res['before']['mean_GPa'],res['after']['mean_GPa'],num,res['after_minus_numeric_GPa'],flush=True)
 return label,res
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:results=dict(pool.map(run,sets))
(root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
