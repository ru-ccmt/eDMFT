from pathlib import Path
import shutil,subprocess,os,json,hashlib,time
r=Path(__file__).resolve().parent;src=Path('/private/tmp/edmft_linear_response_20261108/portable_SrVO3_fixedg/central');d=r/'SrVO3_integration';d.mkdir(exist_ok=True)
for p in src.iterdir():
 if p.is_file() and (p.suffix in {'.struct','.in0','.in1','.in1c','.in2','.in2c','.inc','.inm','.klist','.kgen','.rotlm','.indmfl','.indmfi','.indmf2','.clmsum','.vsp','.vns','.vector','.energy','.nsh'} or p.name in ['dmft2.def','lapw1.def','EF.dat','projectorw.dat','Eorb_imp.dat','sig.inp','sig.inp1']):shutil.copy2(p,d/p.name)
p=d/'SrVO3.in2c';a=p.read_text().splitlines();a[0]='STRT';p.write_text('\n'.join(a)+'\n')
(d/'mpi_prefix.dat').write_text('/opt/homebrew/bin/mpirun -n 4\n')
(d/'SrVO3.stress_response.json').write_text(json.dumps({'enabled':True,'solver':{'mpi':['/opt/homebrew/bin/mpirun','-n','4'],'wien_root':'/Users/haule/wien2k_24','dmft2':'/Users/haule/dbin/dmft2_response','max_iterations':15}},indent=2)+'\n')
names=['SrVO3.in2c','SrVO3.indmf2','dmft2.def','lapw1.def','sig.inp1','SrVO3.clmsum','SrVO3.vector','SrVO3.energy','SrVO3.vsp']
def hashes():return {n:hashlib.sha256((d/n).read_bytes()).hexdigest() for n in names}
before=hashes();env=os.environ.copy()
for k in list(env):
 if k.startswith(('EDMFT_','WIEN_LAPW1_','WIEN_STRESS_','DFTKS_')):del env[k]
env.update(PYTHONPATH=str(r/'python')+':/Users/haule/projects/eDMFT_new/src/python:/Users/haule/dbin',WIENROOT='/Users/haule/wien2k_24',WIEN_DMFT_ROOT='/Users/haule/dbin',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
t=time.monotonic()
with (r/'SrVO3_integration.log').open('w') as log:result=subprocess.run(['/Users/haule/miniconda3-clean/envs/work311/bin/python',str(r/'python/x_dmft.py'),'--execute-only','dmft2'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
after=hashes();report={'returncode':result.returncode,'seconds':time.monotonic()-t,'prepared_and_physical_inputs_unchanged':before==after,'before':before,'after':after}
(r/'SrVO3_integration_result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if result.returncode:raise RuntimeError((r/'SrVO3_integration.log').read_text()[-3000:])
assert before==after
status=json.loads((d/'SrVO3.stress_response.status.json').read_text());assert status['status']=='complete';print(status['workspace'],flush=True)
