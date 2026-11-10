"""Fresh portable native snapshot after independent central reconvergence."""
from pathlib import Path
import sys,os,json,shutil,subprocess,time
root=Path(__file__).resolve().parent;sys.path[:0]=[str(root/'stage'),'/Users/haule/dbin']
from lapw_native_potential_export import prepare_native_potential_core
case=sys.argv[1];src=root/case/'halfstep/center';out=root/case/'refined';out.mkdir(exist_ok=True)
if json.loads((src/'RESULT.json').read_text())['converged'] is not True:raise ValueError('central run not converged')
names=[case+'.'+z for z in ['struct','in0','in1','in1c','in2','in2c','inm','inc','indmfl','indmfi','indmf2','klist','kgen','rotlm','clmsum']]+['EF.dat','Eorb_imp.dat','projectorw.dat','lapw0.def','lapw1.def','dmft2.def','lcore.def']+[p.name for p in src.glob('sig.inp*')]
for name in names:
 if (src/name).exists():shutil.copy2(src/name,out/name)
# The export worker requires the already converged central chemical
# potential. Its response is still solved at fixed N by the Python model.
p=out/(case+'.indmf2');lines=p.read_text().splitlines();j=next(i for i,l in enumerate(lines) if l.strip() and not l.lstrip().startswith('#'));words=lines[j].split();words[1]=(src/'EF.dat').read_text().split()[0];words[2]='0';lines[j]=' '.join(words);p.write_text('\n'.join(lines)+'\n')
env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))};env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',WIEN_STRESS_DEBUG='1')
old=Path('/private/tmp/edmft_xc_core_response_20261108/native');jobs=[('lapw0',old/'SRC_lapw0/lapw0',{}),('lapw1',Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1'),{'WIEN_LAPW1_MATRIX_DUMP':str(out/'native')}),('dmft2',Path('/Users/haule/dbin/dmft2_response'),{'EDMFT_RESPONSE_EXPORT':str(out/'response')}),('lcore',old/'SRC_lcore/lcore',{})]
(out/'core_native.bin').unlink(missing_ok=True);times={}
for name,exe,extra in jobs:
 t=time.perf_counter()
 with (out/(name+'.log')).open('w') as f:subprocess.run([str(exe),name+'.def'],cwd=out,env=env|extra,stdout=f,stderr=subprocess.STDOUT,check=True)
 if name in ['lapw1','dmft2'] and name.upper()+' END' not in (out/(name+'.log')).read_text():raise ValueError(str(out/(name+'.log')))
 times[name]=time.perf_counter()-t
prepare_native_potential_core(out,case,out/'potential_core.npz')
(out/'timings.json').write_text(json.dumps(times,indent=2)+'\n');print(case,times,flush=True)
