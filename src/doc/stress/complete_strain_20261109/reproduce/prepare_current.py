"""Fresh native export using the user's current radial-energy setup, in scratch."""
from pathlib import Path
import os,sys,subprocess,shutil,time,json
root=Path(__file__).resolve().parent;case=sys.argv[1]
old=Path('/private/tmp/edmft_full_k_response_20261108/full_bz')/case/'charge/0'
live={'FeSe':Path('/Users/haule/projects/materials/FeSe_again/FeSe_stress/FeSe_RELA'),'SrVO3':Path('/Users/haule/projects/materials/SrVO3/V_plus_force')}[case]
d=root/case/'current';d.mkdir(parents=True,exist_ok=True)
for suffix in ['struct','in1','in1c','in2','in2c','klist','kgen','rotlm','indmfl','indmfi','indmf2','clmsum','vsp','vns']:
 p=old/(case+'.'+suffix)
 if p.exists():shutil.copy2(p,d/p.name)
for name in ['lapw1.def','dmft2.def','EF.dat','Eorb_imp.dat','projectorw.dat']+[p.name for p in old.glob('sig.inp*')]:
 shutil.copy2(old/name,d/name)
for suffix in ['in1c' if case=='SrVO3' else 'in1']:
 shutil.copy2(live/(case+'.'+suffix),d/(case+'.'+suffix))
# No historical overlap file is reused: the exported normalization matches
# precisely this central band basis and the complete current k mesh.
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))}
env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
exe=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1')
timings={}
for label,exe,extra in [('lapw1',exe,{'WIEN_LAPW1_MATRIX_DUMP':str(d/'native')}),('dmft2','/Users/haule/dbin/dmft2_response',{'EDMFT_RESPONSE_EXPORT':str(d/'response')})]:
 start=time.perf_counter()
 with (d/(label+'.log')).open('w') as f:subprocess.run([str(exe),label+'.def'],cwd=d,env=env|extra,stdout=f,stderr=subprocess.STDOUT,check=True)
 if label.upper()+' END' not in (d/(label+'.log')).read_text():raise RuntimeError(str(d/(label+'.log')))
 timings[label]=time.perf_counter()-start
print(case,timings,flush=True)
(d/'timings.json').write_text(json.dumps(timings,indent=2)+'\n')
