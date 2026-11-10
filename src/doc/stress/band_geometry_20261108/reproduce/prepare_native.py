"""Fixed-potential strained H/O and fixed-coefficient raw-projector endpoints.

Use native LAPW1 for matrices. Then replace only this scratch endpoint vector
with CENTRAL coefficients before DMFT2, to measure the raw operator derivative
without KS eigenvector/window motion. These replay densities/energies are not
physical outputs and must not be used for a stress/free-energy comparison.
"""
from pathlib import Path
import sys,os,shutil,json,time,subprocess,hashlib
import numpy as np
root=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
from lapw_native_snapshot import read_basis_matrix
from frozen_sigma_response import FrozenSigmaResponse
from wstruct import Struct

def run(case):
 source=Path('/private/tmp/edmft_full_k_response_20261108')/case/'charge/0'
 previous=Path('/private/tmp/edmft_screened_native_20261108/three_k')/case
 work=root/case;work.mkdir(exist_ok=True)
 central,records=FrozenSigmaResponse.from_export(source/'response')
 env={k:v for k,v in os.environ.items() if not k.startswith(('EDMFT_','DFTKS_','WIEN_'))}
 env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 s0=Struct();s0.ReadStruct(str(source/(case+'.struct')));cell0=np.array([s0.a,s0.b,s0.c])
 timings=[]
 def command(folder,exe,label,extra):
  start=time.monotonic()
  with (folder/(label+'.log')).open('w') as f:
   subprocess.run([str(exe),label+'.def'],cwd=folder,env=env|extra,stdout=f,stderr=subprocess.STDOUT,check=True)
  log=(folder/(label+'.log')).read_text()
  if any(x in log for x in ['Fortran runtime error','ERROR STOP','Error termination']):raise RuntimeError(log[-1000:])
  if label.upper()+' END' not in log:raise RuntimeError('missing native END: '+str(folder))
  for p in folder.glob(label+'.error*'):
   if p.read_text().strip():raise RuntimeError(p.read_text())
  timings.append(dict(folder=folder.name,program=label,seconds=time.monotonic()-start))
 for step in [.00025,.000125]:
  for label in ['m','p']:
   folder=work/f'{label}_{step}';folder.mkdir(exist_ok=True)
   names=[case+'.'+x for x in ['struct','vsp','vns','in1','in1c','in2','in2c','indmfl','indmfi','indmf2','klist','kgen','rotlm']]
   names+=['EF.dat','Eorb_imp.dat','projectorw.dat','lapw1.def','dmft2.def']+[p.name for p in source.glob('sig.inp*')]
   for name in names:
    if (source/name).exists():shutil.copy2(source/name,folder/name)
   shutil.copy2(previous/f'source_{label}_{step}'/(case+'.struct'),folder/(case+'.struct'))
   shutil.copy2(previous/'central/lapw1.fixgvec',folder/'lapw1.fixgvec');(folder/'lapw1.fixgvec.read').touch()
   exe=Path('/private/tmp/edmft_augmentation_response_20261108/native/SRC_lapw1')/('lapw1c' if case=='SrVO3' else 'lapw1')
   command(folder,exe,'lapw1',{'WIEN_LAPW1_MATRIX_DUMP':str(folder/'native')})
   for k,record in enumerate(records):
    o=read_basis_matrix(folder/f'native.k{k+1:08d}.bin',b'EDMFTOV1')
    np.testing.assert_array_equal(o.labels,record.labels)
   # Preserve the actual strained eigenvectors; the diagnostic replay uses
   # central coefficients and headers solely to export raw local overlaps.
   shutil.copy2(folder/(case+'.vector'),folder/'strained_eigen.vector')
   shutil.copy2(source/(case+'.vector'),folder/(case+'.vector'))
   command(folder,'/Users/haule/dbin/dmft2_response','dmft2',{'EDMFT_RESPONSE_EXPORT':str(folder/'response')})
   model,rr=FrozenSigmaResponse.from_export(folder/'response')
   for a,b in zip(records,rr):np.testing.assert_array_equal(a.coefficients,b.coefficients)
   s=Struct();s.ReadStruct(str(folder/(case+'.struct')))
   cell=np.array([s.a,s.b,s.c]);eta=float(np.mean(cell/cell0)-1)
   record=dict(case=case,eta=eta,cell=cell.tolist(),coefficient_gauge='central native vector replay',
     purpose='raw projector derivative only; replay energy and density are not physical',
     source_fingerprints={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in names if (source/name).exists()})
   (folder/'endpoint.json').write_text(json.dumps(record,indent=2)+'\n')
   print(case,folder.name,'done',flush=True)
 (work/'timings.json').write_text(json.dumps(timings,indent=2)+'\n')

if __name__=='__main__':run(sys.argv[1])
