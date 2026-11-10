"""Validate actual cached G2 matrices and contractions with independent integrals.
Requires NumPy/SciPy and a Fortran compiler. Writes only to EDMFT_STRESS_TEST_DIR.
"""
from pathlib import Path
import os,sys,subprocess,shutil,re,json
import numpy as np
from scipy.special import sph_harm_y,gamma,gammainc
HERE=Path(__file__).resolve().parent
SRC=HERE.parent
OUT=Path(os.environ['EDMFT_STRESS_TEST_DIR']);OUT.mkdir(exist_ok=True,parents=True)
# Only expose the cache to this test program; all production arithmetic is unchanged.
s=(SRC/'stress_g2_matrices.f90').read_text().replace('  private\n','  private\n  public :: atoms\n',1)
(OUT/'g2_under_test.f90').write_text(s)
s=(SRC/'stress_kernels.f90').read_text()
m=re.search(r'^\s*subroutine C_factor\(.*?^\s*end subroutine C_factor',s,re.M|re.S|re.I)
(OUT/'c_factor.f90').write_text(m.group(0))
compiler=os.environ.get('FC',shutil.which('gfortran') or 'gfortran')
cmd=[compiler,'-O0','-g','-fcheck=all','-ffree-line-length-none',str(HERE/'g2_fixture.f90'),str(OUT/'g2_under_test.f90'),str(OUT/'c_factor.f90'),str(SRC/'ylm.f90'),str(SRC/'atpar/rint13.f90'),str(HERE/'g2_fixture_main.f90'),*(['-framework','Accelerate'] if sys.platform=='darwin' else ['-llapack','-lblas']),'-o',str(OUT/'g2_fixture')]
subprocess.run(cmd,cwd=OUT,check=True)
subprocess.run([str(OUT/'g2_fixture')],cwd=OUT,check=True)
with (OUT/'g2_fixture.bin').open('rb') as h:
 def read(dtype,shape):return np.fromfile(h,dtype=dtype,count=int(np.prod(shape))).reshape(shape,order='F')
 n,nb,nrad=read('i4',(3,));r0,dx=read('f8',(2,));rot=read('f8',(3,3))
 ly=read('i4',(n,))-1;l=read('i4',(n,));p=read('i4',(n,));e=read('f8',(n,));partner=read('i4',(n,))-1
 mats={}
 for key,shape in [('s0',(n,n)),('sr',(n,n,6)),('sa',(n,n,6)),('v0',(n,n)),('vr',(n,n,6)),('va',(n,n,6))]:mats[key]=read('c16',shape)
 a=read('c16',(9,nb));b=read('c16',(9,nb));c=read('c16',(9,nb,2))
 da=read('c16',(9,nb,6));db=read('c16',(9,nb,6));dc=read('c16',(9,nb,2,6));eb=read('f8',(nb,))
 kernels=[read('c16',(3,3,nb)) for _ in range(6)]
 assert not h.read(1)
alpha=.7+.2*p;beta=alpha[:,None]+alpha[None,:];power=l[:,None]+l[None,:]+2
R=r0*np.exp((nrad-1)*dx)
def integral(power):return gamma(power+1)*gammainc(power+1,beta*R)/beta**(power+1)
radS=integral(power);radD=l[:,None]*radS-alpha[:,None]*integral(power+1)
radV=.31*integral(power+2);radW=.31*(l[:,None]*integral(power+2)-alpha[:,None]*integral(power+3))
x,w=np.polynomial.legendre.leggauss(32);phi=np.arange(64)*2*np.pi/64
x,phi=np.meshgrid(x,phi,indexing='ij');w=np.repeat(w,64)*2*np.pi/64
v=np.stack([np.sqrt(1-x*x)*np.cos(phi),np.sqrt(1-x*x)*np.sin(phi),x],axis=-1).reshape(-1,3)
m=ly-l*(l+1)
def harmon(vec):
 th=np.arccos(np.clip(vec[:,2]/np.linalg.norm(vec,axis=1),-1,1));ph=np.arctan2(vec[:,1],vec[:,0])
 return sph_harm_y(l[None,:],m[None,:],th[:,None],ph[:,None])
y=harmon(v);vp=sph_harm_y(2,0,np.arccos(v[:,2]),np.arctan2(v[:,1],v[:,0]))
axes=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
refs={'s0':radS*(ly[:,None]==ly[None,:]),'v0':radV*(y.conj().T@((w*vp)[:,None]*y))}
for key in ['sr','sa','vr','va']:refs[key]=np.zeros((n,n,6),complex)
for q,(aa,bb) in enumerate(axes):
 gen=np.zeros((3,3));gen[aa,bb]+=.5;gen[bb,aa]+=.5
 dy=(harmon(v@(np.eye(3)+1e-6*gen))-harmon(v@(np.eye(3)-1e-6*gen)))/(2e-6)
 refs['sr'][:,:,q]=radD*(y.conj().T@((w*v[:,aa]*v[:,bb])[:,None]*y))
 refs['sa'][:,:,q]=radS*(dy.conj().T@(w[:,None]*y))
 refs['vr'][:,:,q]=radW*(y.conj().T@((w*vp*v[:,aa]*v[:,bb])[:,None]*y))
 refs['va'][:,:,q]=radV*(dy.conj().T@((w*vp)[:,None]*y))
errors={key:float(np.max(abs(mats[key]-ref))) for key,ref in refs.items()}
for key,value in errors.items():assert value<2e-7,(key,value)
coef=np.empty((n,nb),complex);deriv=np.empty((n,nb,6),complex)
for i in range(n):
 if p[i]==1:coef[i]=a[ly[i]];deriv[i]=da[ly[i]]
 elif p[i]==2:coef[i]=b[ly[i]];deriv[i]=db[ly[i]]
 else:coef[i]=c[ly[i],:,p[i]-3];deriv[i]=dc[ly[i],:,p[i]-3,:]
deriv[:,:,:3]+=.5*coef[:,:,None]
D=np.zeros((n,nb,3,3),complex)
for q,(aa,bb) in enumerate(axes):D[:,:,aa,bb]=deriv[:,:,q];D[:,:,bb,aa]=deriv[:,:,q]
D=np.einsum('ai,nkij,bj->nkab',rot,D,rot)
HC=e[:,None]*coef
for i,j in enumerate(partner):
 if j>=0:HC[j]+=coef[i]
HC-=coef*eb[None,:]
expected=[]
for q,(aa,bb) in enumerate(axes):
 expect=[np.sum(D[:,:,aa,bb].conj()*(refs['s0']@HC),axis=0),np.sum(coef.conj()*(refs['sr'][:,:,q]@HC),axis=0),np.sum(coef.conj()*(refs['sa'][:,:,q]@HC),axis=0),np.sum(D[:,:,aa,bb].conj()*(refs['v0']@coef),axis=0),np.sum(coef.conj()*(refs['vr'][:,:,q]@coef),axis=0),np.sum(coef.conj()*(refs['va'][:,:,q]@coef),axis=0)]
 for row,(actual,target) in enumerate(zip(kernels,expect)):
  error=float(np.max(abs(actual[aa,bb]-target)))
  assert error<2e-6,(q,row,error)
  expected.append(error)
result={'primitive_count':int(n),'two_LO_energies_per_l':True,'complex_pseudo_bands':int(nb),'nonidentity_local_axes':True,'matrix_max_errors':errors,'six_G2_contractions_max_error':max(expected)}
(OUT/'complete_g2_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2));print('PASS: complete all-LO radial matrices and all six complex G2 contractions')
