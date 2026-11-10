from pathlib import Path
import shutil,os,subprocess,time,json
root=Path(__file__).parent;r=root/'MnO';src=r/'center';d=r/'dftks_control';d.mkdir(exist_ok=True)
for p in src.iterdir():
 if p.is_file() and (p.suffix in ['.struct','.in0','.in1','.in1c','.in2','.inc','.inm','.klist','.kgen','.rotlm','.vsp','.vns','.clmsum','.clmcor','.scf0','.scfc','.indmfl','.indmfi','.indmf2','.def'] or p.name in ['sig.inp','sig.inp1','projectorw.dat','SOlapm.dat','EF.dat','Eorb_imp.dat']):shutil.copy2(p,d/p.name)
env=os.environ|{'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','EDMFT_STRESS_STERNHEIMER_DUMP':'1','EDMFT_STRESS_RESPONSE_DUMP':'1'}
# Serial dftKS produces one vector file; dmft2 distributes its k points.
for label,cmd in [('lapw1',['/Users/haule/dbin/lapw1','lapw1.def']),('dmft2',['/opt/homebrew/bin/mpirun','-n','4','/Users/haule/dbin/dmft2','dmft2.def'])]:
 t=time.monotonic()
 with (d/(label+'.log')).open('w') as f:p=subprocess.run(cmd,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 print(label,p.returncode,time.monotonic()-t,flush=True);assert p.returncode==0
import analyze_cases as a
analytic=float(__import__('numpy').trace(a.tensor(d/'MnO.scf0',':STR_DENS')-a.tensor(d/'MnO.scfc',':STR_CORE')-a.tensor(d/'MnO.scf2',':STR_VAL'))*14710.507641078/3)
print('dftKS analytic',analytic,flush=True);(d/'result.json').write_text(json.dumps({'analytic_GPa':analytic},indent=2))
