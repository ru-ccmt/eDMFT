"""Native mixer integration check in a disposable copied SrVO3 calculation."""
from pathlib import Path
import sys,os,json,shutil,hashlib,subprocess,re
root=Path(__file__).resolve().parent;case='SrVO3';live=Path('/Users/haule/projects/materials/SrVO3/V_plus_force');d=root/case/'publication';d.mkdir(exist_ok=True)
for suffix in ['struct','inm','inc','clmsum','clmval','clmcor','vsp','vns','scf0','scf1','scf2','scfc']:
 shutil.copy2(live/(case+'.'+suffix),d/(case+'.'+suffix))
for p in live.glob(case+'.in*'):
 if p.is_file():shutil.copy2(p,d/p.name)
shutil.copy2(live/'Eorb_imp.dat',d/'Eorb_imp.dat')
report=json.loads((root/case/'physical_total/report.json').read_text());report['converged']=True;p=d/'total_stress.json';p.write_text(json.dumps(report,indent=2)+'\n')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
status=dict(schema=1,backend='analytic',status='complete',stress_ready=True,stress_applied=False,converged=True,reference_id='native-publication-test',report=str(p),report_sha256=sha(p),parent_fingerprints={case+'.clmsum':sha(d/(case+'.clmsum'))})
(d/(case+'.stress_response.status.json')).write_text(json.dumps(status))
env=dict(os.environ,PYTHONPATH=str(root/'stage')+':/Users/haule/dbin',PYTHONDONTWRITEBYTECODE='1',WIENROOT='/Users/haule/wien2k_24',WIEN_DMFT_ROOT='/Users/haule/dbin',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
for label,command in [('merge',[sys.executable,str(root/'stage/merge_scf.py'),'--replace','-f',case]),('mixer',['/Users/haule/wien2k_24/x','mixer','-f',case])]:
 with (d/(label+'.log')).open('w') as f:subprocess.run(command,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 if label=='mixer' and 'MIXER END' not in (d/'mixer.log').read_text():raise ValueError(str(d/'mixer.log'))
rows=[]
text=(d/(case+'.scfm')).read_text()
for i in range(1,4):
 m=re.search(r'^:STRESS_GPa%03d:\s*(\S+)\s*(\S+)\s*(\S+)'%i,text,re.M)
 if m is None:raise ValueError('mixer total absent')
 rows.append(list(map(float,m.groups())))
import numpy as np
error=float(np.max(abs(np.array(rows)-report['stress_GPa'])))
if error>1e-5:raise ValueError('native mixer disagrees with response tensor')
result=dict(scope='SCF publication plumbing, copied native inputs',expected_GPa=report['stress_GPa'],native_mixer_GPa=rows,max_error_GPa=error,stress_applied=json.loads((d/(case+'.stress_response.status.json')).read_text())['stress_applied'])
(d/'REPORT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
