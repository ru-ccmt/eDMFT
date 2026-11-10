from material_augmentation import *
import os,shutil,time
case=sys.argv[1];mat=Material(case);s=mat.s
for ia in range(s.nat):
 for l in range(5):mat.radial(ia,l)
print('radial prepared',case,flush=True)
root=R/case/'vns_probe';root.mkdir(exist_ok=True)
env=os.environ.copy()
for k in list(env):
 if k.startswith(('WIEN_','EDMFT_','DFTKS_')):del env[k]
env.update(OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
exe=Path('/private/tmp/edmft_native_operator_20261108/native_probe/SRC_lapw1')/('lapw1c' if mat.complex else 'lapw1')
step=.01;ends=[]
for sign in [-1,1]:
 d=root/('minus' if sign<0 else 'plus');ends.append(d);d.mkdir(exist_ok=True)
 for p in mat.source.iterdir():
  if p.is_file() and (p.suffix in {'.struct','.vsp','.vns','.in1','.in1c','.klist'} or p.name=='lapw1.def'):shutil.copy2(p,d/p.name)
 p=d/(case+'.vns');text=p.read_text();patches=[]
 for block in re.finditer(r'ATOMNUMBER\s*=\s*(\d+)',text):
  ia=int(block[1])-1;end=text.find('ATOMNUMBER',block.end());end=len(text) if end<0 else end
  section=text[block.end():end]
  for channel in re.finditer(r'VLM\(R\) FOR L\s*(-?\d+)\s+M=\s*(\d+)',section):
   offset=block.end()+channel.end()
   for match in list(FLOAT.finditer(text[offset:end]))[:s.jrj[ia]]:
    patches.append((offset+match.start(),offset+match.end(),float(match.group().replace('D','E'))*(1+sign*step)))
 for a,b,value in reversed(patches):text=text[:a]+f'{value:{b-a}.12E}'+text[b:]
 p.write_text(text)
 if not (d/'done').exists():
  with (d/'run.log').open('w') as log:proc=subprocess.run([str(exe),'lapw1.def'],cwd=d,env=env|{'WIEN_LAPW1_MATRIX_DUMP':str(d/'native')},stdout=log,stderr=subprocess.STDOUT)
  proc.check_returncode();(d/'done').touch()
rows=[]
for ik,(ev,c,labels) in enumerate(mat.cv):
 targets=[np.zeros((len(ev),len(ev)),complex),np.zeros((len(ev),len(ev)),complex)]
 latom=0
 for ia in range(s.nat):
  for mu in range(s.mult[ia]):
   z,dz,*indices=mat.coefficients(ia,latom,ik)
   for transpose in [False,True]:targets[int(transpose)]+=z.conj().T@mat.kernel(ia,mat.vns[ia],indices,transpose)@z
   latom+=1
 hs=[read_basis_matrix(d/f'native.ham.k{ik+1:08d}.bin',b'EDMFTHV1').matrix for d in ends];dh=c.conj().T@((hs[1]-hs[0])/(2*step))@c
 row=dict(k=ik+1,scale=float(np.linalg.norm(dh)))
 for tr,h in enumerate(targets):row[str(tr)]=dict(relative=float(np.linalg.norm(h-dh)/np.linalg.norm(dh)),maximum=float(np.max(abs(h-dh))),hermitian=float(np.max(abs(h-h.conj().T))))
 rows.append(row);print(case,row,flush=True)
(root/'report.json').write_text(json.dumps(rows,indent=2)+'\n')
