from pathlib import Path
import os,sys,subprocess,json,time
root=Path(__file__).parent
jobs=json.loads((root/'native_jobs.json').read_text())
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','WIEN_','DFTKS_'))}
env.update(WIENROOT='/Users/haule/wien2k_24',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
records=[]
for job in jobs:
 w=Path(job['directory']);c=job['case']
 def run(label,argv,extra=None,marker=None):
  t=time.perf_counter()
  with (w/(label+'.log')).open('w') as f:r=subprocess.run(argv,cwd=w,env=env|(extra or {}),stdout=f,stderr=subprocess.STDOUT)
  txt=(w/(label+'.log')).read_text(errors='replace')
  ok=r.returncode==0 and (marker is None or marker+' END' in ' '.join(txt.split())) and 'Fortran runtime error' not in txt
  row=dict(material=job['material'],name=job['name'],stage=label,seconds=time.perf_counter()-t,completed=ok,returncode=r.returncode)
  records.append(row); print(json.dumps(row),flush=True)
  (root/'native_runs.json').write_text(json.dumps(records,indent=2)+'\n')
  if not ok:print(txt[-2500:],flush=True);raise RuntimeError(label)
 for label in ['lapw0','lcore']:run(label+'_definition',['/Users/haule/wien2k_24/x',label,'-f',c,'-d'])
 run('lapw0',['/Users/haule/dbin/lapw0_response','lapw0.def'],marker='LAPW0')
 run('lapw1',['/Users/haule/dbin/lapw1_response','lapw1.def'],{'WIEN_LAPW1_MATRIX_DUMP':str(w/'native')},'LAPW1')
 for suf in ['']+(['dn'] if job['material']=='MnTe' else []):
  label='dmft2'+suf
  run(label,['/opt/homebrew/bin/mpirun','-np','4','/Users/haule/dbin/dmft2_response',label+'.def'],{'EDMFT_RESPONSE_EXPORT':str(w/('response' if not suf else 'response_dn'))},'DMFT2')
 run('lcore',['/Users/haule/dbin/lcore_response','lcore.def'],marker='CORE')
