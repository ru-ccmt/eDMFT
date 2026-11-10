from pathlib import Path
import json,os,shutil,struct,subprocess,sys
import numpy as np
from scipy.io import FortranFile,FortranEOFError
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'src/python'))
from check_oneparticle_stress_fd import read_basis_matrix

def vectors(path,nat,cplx):
 with FortranFile(path) as f:
  for _ in range(2*nat):f.read_record('u1')
  while True:
   try:h=f.read_record('u1').tobytes()
   except FortranEOFError:return
   n,nb=struct.unpack('<ddd10siid',h)[4:6]; labels=f.read_record('i4').reshape(3,n,order='F');E=[];C=[]
   for _ in range(nb):
    E.append(struct.unpack('<id',f.read_record('u1').tobytes())[1]);C.append(f.read_record('c16' if cplx else 'f8'))
   yield np.array(E),np.array(C).T,labels

def run(case,kind):
 nat,cplx=(4,True) if case=='SrVO3' else (2,False)
 src=Path('/private/tmp/radial_stress_response_20261108')/(case+'_'+kind)/'center'
 d=root/'tests'/(case+'_wien_'+kind);d.mkdir(exist_ok=True)
 for p in src.iterdir():
  if p.is_file() and (p.suffix in {'.struct','.vsp','.vns','.klist','.in1','.in1c'} or p.name=='lapw1.def'):shutil.copy2(p,d/p.name)
 p=d/(case+'.klist');p.write_text('\n'.join(p.read_text().splitlines()[:3])+'\nEND\n')
 env=os.environ.copy()
 for k in list(env):
  if k.startswith(('DFTKS_','WIEN_LAPW1_','WIEN_STRESS_','EDMFT_')):del env[k]
 env.update(OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1',WIEN_LAPW1_MATRIX_DUMP=str(d/'native'))
 exe=root/'wien/SRC_lapw1'/('lapw1c' if cplx else 'lapw1')
 with (d/'native.log').open('w') as out:r=subprocess.run([str(exe),'lapw1.def'],cwd=d,env=env,stdout=out,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError((d/'native.log').read_text())
 vals=[]
 for ik,(E,C,labels) in enumerate(vectors(d/(case+'.vector'),nat,cplx),1):
  O=read_basis_matrix(d/f'native.k{ik:08d}.bin',b'EDMFTOV1');H=read_basis_matrix(d/f'native.ham.k{ik:08d}.bin',b'EDMFTHV1')
  assert np.array_equal(labels,O.labels)
  vals.append({'k':ik,'overlap':float(np.max(abs(C.conj().T@O.matrix@C-np.eye(len(E))))),'H':float(np.max(abs(C.conj().T@H.matrix@C-np.diag(E)))),'H_conjugate':float(np.max(abs(C.conj().T@H.matrix.conj()@C-np.diag(E))))})
 print(case,kind,vals,flush=True)
 (d/'identities.json').write_text(json.dumps(vals,indent=2)+'\n')
 return vals
if __name__=='__main__':run(*sys.argv[1:])
