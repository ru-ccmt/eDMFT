from pathlib import Path
import os,sys,re,struct,subprocess,shutil,json
import numpy as np
from scipy.io import FortranFile
from scipy.special import spherical_jn,eval_legendre
R=Path(__file__).resolve().parent
sys.path.insert(0,'/Users/haule/projects/eDMFT_new/src/python')
sys.path.insert(0,'/private/tmp/edmft_radial_all_l_20261108')
from check_oneparticle_stress_fd import read_basis_matrix
from test_wien_matrices import vectors
from wstruct import Struct
from lapw_response_assembly import NativeBandProjection
case=sys.argv[1];cplx=case=='SrVO3';nat=4 if cplx else 2
root=R/'native_probe'/case;root.mkdir(exist_ok=True)
src=Path('/private/tmp/edmft_radial_all_l_20261108/tests')/(case+'_all_l_restored')/'center'
env=os.environ.copy()
for key in list(env):
 if key.startswith(('WIEN_','EDMFT_','DFTKS_')):del env[key]
env.update(OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
exe=R/'native_probe/SRC_lapw1'/('lapw1c' if cplx else 'lapw1')

def run(name,step):
 d=root/name;d.mkdir(exist_ok=True)
 for p in src.iterdir():
  if p.is_file() and (p.suffix in {'.struct','.vsp','.vns','.klist','.in1','.in1c'} or p.name=='lapw1.def'):shutil.copy2(p,d/p.name)
 p=d/(case+'.klist');p.write_text('\n'.join(p.read_text().splitlines()[:3])+'\nEND\n')
 en=env|{'WIEN_LAPW1_MATRIX_DUMP':str(d/'native')}
 if step:
  text=(d/(case+'.vsp')).read_text();blocks=re.split(r'(ATOMNUMBER\s*=\s*\d+)',text)
  field=re.compile(r'[-+]?\d+\.\d+[DE][-+]\d+')
  matches=list(field.finditer(blocks[2]));assert len(matches)==len(v)
  for match,delta in reversed(list(zip(matches,step*dv))):
   a,b=match.span();value=float(match.group().replace('D','E'))+delta
   blocks[2]=blocks[2][:a]+f'{value:{b-a}.12E}'+blocks[2][b:]
  p=d/'basis.vsp';p.write_text(''.join(blocks));en['WIEN_LAPW1_BASIS_VSP']=str(p)
 if not (d/'done').exists():
  with (d/'run.log').open('w') as log:r=subprocess.run([str(exe),'lapw1.def'],cwd=d,env=en,stdout=log,stderr=subprocess.STDOUT)
  r.check_returncode();(d/'done').touch()
 return d
s=Struct();s.ReadStruct(str(src/(case+'.struct')));n=s.jrj[0];r0=s.r0[0];rmt=s.rmt[0];dx=np.log(rmt/r0)/(n-1);r=r0*np.exp(dx*np.arange(n))
block=re.split(r'ATOMNUMBER\s*=\s*\d+', (src/(case+'.vsp')).read_text())[1]
v=np.array([float(t.replace('D','E')) for t in re.findall(r'[-+]?\d+\.\d+[DE][-+]\d+',block)])
assert len(v)==n
dv=.07*r*np.sin(1.3*r);center=run('center',0)
with FortranFile(center/(case+'.vector')) as f:header=f.read_record('f8')
e=header[5];lapw=e<150;e=e if lapw else e-200
cv=list(vectors(center/(case+'.vector'),nat,cplx));tangents=[]
for ik,(ev,c,labels) in enumerate(cv,1):
 p=center/f'native.k{ik:08d}.bin';o=read_basis_matrix(p,b'EDMFTOV1');npw=o.nbasis-o.nlo
 with p.open('rb') as f:f.seek(44);kf=np.fromfile(f,'f8',3)
 q=2*np.pi*(labels[:,:npw].T+kf)/np.array([s.a,s.b,s.c]);qn=np.linalg.norm(q,axis=1);qu=q/qn[:,None]
 b=np.array([spherical_jn(5,qn*rmt),qn*spherical_jn(5,qn*rmt,derivative=True)])
 inp=root/f'channel{ik}.input'
 with inp.open('w') as f:
  f.write(f'{n} {npw} {int("RELA" in (center/(case+".struct")).read_text())} 5 {int(lapw)}\n')
  f.write(f'{r0:.17e} {dx:.17e} {rmt:.17e} {e:.17e} {s.Znuc[0]:.17e}\n')
  for a in [v,dv,b.ravel(order='F')]:f.write(' '.join(f'{x:.17e}' for x in a)+'\n')
 out=root/f'channel{ik}.bin';subprocess.run([str(R/'probe_radial'),str(inp),str(out)],check=True)
 rad=np.fromfile(out,'f8').reshape((2,npw,npw)).transpose(0,2,1)
 label=labels[:,:npw].T;phase=sum(np.exp(2j*np.pi*((label@pos)[:,None]-(label@pos)[None,:])) for pos in s.pos[:s.mult[0]])
 angular=4*np.pi/(s.a*s.b*s.c)*11*eval_legendre(5,np.clip(qu@qu.T,-1,1))*phase
 projection=NativeBandProjection(angular,c[:npw])
 transpose_error=projection.operator.check_transpose(samples=2)
 assert transpose_error<1e-12
 full=[]
 for a in rad:
  raw=a*angular;triangle=np.triu(raw);fold=triangle+np.triu(raw,1).conj().T
  aa=np.zeros((o.nbasis,o.nbasis),complex);aa[:npw,:npw]=fold;full.append(aa)
  np.testing.assert_allclose(projection.project(a),c.conj().T@aa@c,rtol=1e-11,atol=1e-15)
 tangents.append(full)
 print(case,ik,"native band projection transpose error",transpose_error,flush=True)
rows=[]
for step in [2e-3,1e-3,5e-4]:
 ends=[run(('minus' if sg<0 else 'plus')+f'_{step:.0e}',sg*step) for sg in [-1,1]]
 vals=[]
 for ik,((ev,c,labels),(dh,do)) in enumerate(zip(cv,tangents),1):
  row=dict(k=ik)
  for name,magic,target in [('H',b'EDMFTHV1',dh),('O',b'EDMFTOV1',do)]:
   suffix=f'native.{"ham." if name=="H" else ""}k{ik:08d}.bin'
   mats=[read_basis_matrix(d/suffix,magic) for d in ends]
   assert all(np.array_equal(m.labels,labels) for m in mats)
   fd=(mats[1].matrix-mats[0].matrix)/(2*step)
   row[name+'_relative']=float(np.linalg.norm(fd-target)/np.linalg.norm(target))
   row[name+'_max']=float(np.max(abs(fd-target)))
   band=c.conj().T@(fd-target)@c;row[name+'_band_max']=float(np.max(abs(band)))
  vals.append(row)
 rows.append(dict(step=step,kpoints=vals));print(case,rows[-1],flush=True)
(root/'report.json').write_text(json.dumps(rows,indent=2)+'\n')
