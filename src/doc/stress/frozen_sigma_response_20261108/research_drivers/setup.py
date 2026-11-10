from pathlib import Path
import shutil,re,json,os,subprocess,time
ROOT=Path(__file__).parent
SOURCES={'SrVO3':Path('/private/tmp/projector_stationarity_probe_20261107/SrVO3'),'FeSe':Path('/private/tmp/projector_cross_material_20261107/FeSe/center')}
ENV=os.environ|{'WIENROOT':'/Users/haule/wien2k_24','WIEN_DMFT_ROOT':'/Users/haule/dbin','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
for key in list(ENV):
 if key.startswith(('EDMFT_','WIEN_STRESS_','DFTKS_','WIEN_LAPW1_')): del ENV[key]
EXE=ROOT/'build/src/dmft2/dmft2'
def prepare(case,name):
 src=SOURCES[case];d=ROOT/case/name;d.mkdir(parents=True,exist_ok=True)
 for p in src.iterdir():
  if not p.is_file():continue
  if p.name.startswith(case+'.') and p.suffix in ['.struct','.in0','.in1','.in1c','.in2','.in2c','.inc','.inm','.klist','.kgen','.rotlm','.vsp','.vns','.clmsum','.clmcor','.scf0','.scfc','.indmfl','.indmfi','.indmf2','.energy'] or p.name in ['EF.dat','SOlapm.dat','projectorw.dat','Eorb_imp.dat','sig.inp','dmft2.def','_processes_'] or re.fullmatch(r'sig\.inp\d+',p.name):shutil.copy2(p,d/p.name)
 vp=d/(case+'.vector')
 if not vp.is_symlink() and not vp.exists() and (src/(case+'.vector')).exists():vp.symlink_to(src/(case+'.vector'))
 # The source SrVO3 uses dftKS complex vectors under the ordinary name.
 for suf in ['in2','in2c']:
  p=d/(case+'.'+suf)
  if p.exists():a=p.read_text().splitlines();a[0]='TOT '+a[0][4:];p.write_text('\n'.join(a)+'\n')
 p=d/(case+'.indmf2');a=p.read_text().splitlines();idx=next(i for i,l in enumerate(a) if l.strip() and not l.startswith('#'));mu=(src/'EF.dat').read_text().split()[0];a[idx]=f'c {mu} 0 1.0 .T. 0.02';p.write_text('\n'.join(a)+'\n')
 p=d/'dmft2.def'
 if case=='SrVO3':
  for z in src.glob(case+'.vector_*'):
   target=d/z.name
   if not target.exists():target.symlink_to(z)
  for z in src.glob(case+'.energy_*'):shutil.copy2(z,d/z.name)
 else:p.write_text(p.read_text().replace('_x',''))
 return d

def run(case,name,extra={}):
 d=prepare(case,name);t=time.monotonic()
 with (d/'run.log').open('w') as f:r=subprocess.run(['/opt/homebrew/bin/mpirun','-n','4',str(EXE),'dmft2.def'],cwd=d,env=ENV|extra,stdout=f,stderr=subprocess.STDOUT)
 result={'rc':r.returncode,'seconds':time.monotonic()-t};(d/'timing.json').write_text(json.dumps(result));print(case,name,result,flush=True)
 if r.returncode:raise RuntimeError(d)
 return d
if __name__=='__main__':
 import concurrent.futures
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  print(list(pool.map(lambda c:run(c,'central',{'EDMFT_RESPONSE_EXPORT':str(ROOT/c/'central/response')}),SOURCES)))
