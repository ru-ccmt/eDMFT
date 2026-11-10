from pathlib import Path
import os,subprocess,shutil,json,hashlib
r=Path(__file__).resolve().parent;d=r/'SrVO3_installed_noresponse';src=r/'SrVO3_integration'
shutil.copytree(src,d,ignore=shutil.ignore_patterns('SrVO3.response_runs','SrVO3.stress_response.status.json','.SrVO3.dmft2-response.lock'))
p=d/'SrVO3.in2c';a=p.read_text().splitlines();a[0]='STRFT';p.write_text('\n'.join(a)+'\n');old=p.read_bytes()
env=os.environ.copy()
for k in list(env):
 if k.startswith(('EDMFT_','WIEN_STRESS_','WIEN_LAPW1_','DFTKS_')):del env[k]
env.update(WIENROOT='/Users/haule/wien2k_24',WIEN_DMFT_ROOT='/Users/haule/dbin',PYTHONPATH='/Users/haule/dbin',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
with (r/'installed_smoke.log').open('w') as log:rr=subprocess.run(['/Users/haule/miniconda3-clean/envs/work311/bin/python','/Users/haule/dbin/x_dmft.py','--execute-only','--no-response','dmft2'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
assert rr.returncode==0,(r/'installed_smoke.log').read_text()[-2000:]
assert p.read_bytes()==old
status=json.loads((d/'SrVO3.stress_response.status.json').read_text());assert status['status']=='disabled';assert not (d/'SrVO3.response_runs').exists()
report={'installed_execution_passed':True,'STRFT_preserved':True,'response_suppression_passed':True}
(r/'installed_smoke_result.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
