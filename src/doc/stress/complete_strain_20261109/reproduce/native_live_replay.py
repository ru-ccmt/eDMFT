"""Same-input native energy replay in isolation; never overwrite material runs."""
from pathlib import Path
import sys,os,json,re,time,shutil,subprocess
import numpy as np
root=Path(__file__).resolve().parent;case=sys.argv[1]
base,dirs={'FeSe':(Path('/Users/haule/projects/materials/FeSe_again/FeSe_stress'),['FeSe_RELA_strain_m001','FeSe_RELA','FeSe_RELA_strain_p001']),'SrVO3':(Path('/Users/haule/projects/materials/SrVO3'),['V_strain_m001','V_plus_force','V_strain_p001'])}[case]
old=Path('/private/tmp/edmft_xc_core_response_20261108')
with np.load('/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/'+case+'_potential_core.npz') as z:
 meta=json.loads(str(z['metadata']));mult=z['mult'].copy();cell=z['cell'][:3].copy()
weights=np.array([mult[a]*st['occupation'] for a,c in enumerate(meta['core']) for st in c['states']]);shift=sum(mult[a]*c['shift_Ry']*sum(st['occupation'] for st in c['states']) for a,c in enumerate(meta['core']))
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))};env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',WIEN_STRESS_DEBUG='1')
out=root/case/'native_live';out.mkdir(exist_ok=True);rows=[]
for label,name in zip(['minus','center','plus'],dirs):
 src=base/name;d=out/label;d.mkdir(exist_ok=True)
 if (d/'RESULT.json').exists():rows.append(json.loads((d/'RESULT.json').read_text()));continue
 names=[case+'.'+s for s in ['struct','in0','in1','in1c','in2','in2c','inm','inc','indmfl','indmfi','indmf2','klist','kgen','rotlm','clmsum']]+['EF.dat','Eorb_imp.dat','projectorw.dat']+[p.name for p in src.glob('sig.inp*')]
 for name in names:
  p=src/name
  if p.exists():shutil.copy2(p,d/name)
 for program,folder in [('lapw0',old/case/'xc/0'),('lcore',old/case/'hxc_core/0'),('lapw1',root/case/'current'),('dmft2',root/case/'current')]:shutil.copy2(folder/(program+'.def'),d/(program+'.def'))
 in2=d/(case+('.in2c' if case=='SrVO3' else '.in2'));lines=in2.read_text().splitlines();lines[0]='TOT';in2.write_text('\n'.join(lines)+'\n')
 (d/'core_native.bin').unlink(missing_ok=True)
 commands=[('lapw0',old/'native/SRC_lapw0/lapw0'),('lapw1',Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1')),('dmft2',Path('/Users/haule/dbin/dmft2_response')),('lcore',old/'native/SRC_lcore/lcore')];timing={}
 for program,exe in commands:
  t=time.perf_counter()
  with (d/(program+'.log')).open('w') as f:subprocess.run([str(exe),program+'.def'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
  timing[program]=time.perf_counter()-t
  if program in ['lapw1','dmft2'] and program.upper()+' END' not in (d/(program+'.log')).read_text():raise ValueError(str(d))
 den=float(re.findall(r'^:DBG_EN_DEN_TOTAL:\s*(\S+)',(d/(case+'.scf0')).read_text(),re.M)[-1])
 vsum=float(re.findall(r'^:VSUM.*?=\s*(\S+)',(d/(case+'.scf2')).read_text(),re.M)[-1]);es=[]
 with (d/'core_native.bin').open('rb') as f:
  while True:
   h=np.fromfile(f,'i4',7)
   if not len(h):break
   es.append(float(np.fromfile(f,'f8',4)[0]));f.seek(int(4*h[3]*8),1)
 core=float(2*weights@np.array(es)+shift);row=dict(label=label,directory=str(d),DEN_Ry=den,core_Ry=core,VSUM_Ry=vsum,K_Ry=den+core+vsum,timings=timing)
 (d/'RESULT.json').write_text(json.dumps(row,indent=2)+'\n');rows.append(row);print(case,label,row,flush=True)
D=(rows[2]['K_Ry']-rows[0]['K_Ry'])/.002
report=dict(case=case,rows=rows,derivative_Ry=D,mean_tensile_GPa=D*14710.5077/(3*np.prod(cell)),scope='fresh energy replay at supplied densities, not reconvergence')
(out/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
