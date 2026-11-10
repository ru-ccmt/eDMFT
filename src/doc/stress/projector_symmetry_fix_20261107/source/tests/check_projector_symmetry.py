"""Independent checks of the production normalization adjoint and strain frame.

The forward symmetry and its adjoint are extracted from dmft2. Two equivalent
blocks use different orbital representations, so equivalent-atom averaging
and orbital averaging do not commute. Finite differences of the complete
inverse-square-root functional catch wrong operation order and conjugation.
Only the symmetry-operation provider is mocked; production averaging and
Sylvester routines are used unchanged. No material calculation is modified.
"""
from pathlib import Path
import ctypes as ct
import hashlib
import json
import os
import re
import subprocess
import tempfile
import sys
import shutil
import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[3]
# A scratch copy can explicitly select the real source tree.
ROOT = Path(os.environ.get('EDMFT_SOURCE_ROOT', str(ROOT)))
OUT = Path(os.environ.get('EDMFT_SYMMETRY_TEST_DIR') or
           os.environ.get('EDMFT_STRESS_TEST_DIR') or
           tempfile.mkdtemp(prefix='edmft_symmetry_'))
OUT.mkdir(parents=True, exist_ok=True)
def routine(file, name):
    text = (ROOT/'src/dmft2'/file).read_text()
    pattern = rf'^\s*SUBROUTINE\s+{name}\s*\(.*?^\s*END SUBROUTINE\s+{name}\s*$'
    match = re.search(pattern, text, re.M|re.S|re.I)
    if match is None: raise RuntimeError(name)
    return match.group(0)+'\n'

stubs = '''module sym2
integer,parameter :: iord=2
end module
module dmfts
integer,parameter :: maxdim=3,ncix=2,natom=2,lmaxp=0
end module
module com_mpi
integer :: myrank=1,master=0
end module
module test_rotations
complex*16 :: matrices(3,3,2,2)
end module
subroutine set_test_rotations(r)
use test_rotations
complex*16,intent(in) :: r(3,3,2,2)
matrices=r
end subroutine
subroutine CmpSymmetryTransformation(ig,RT,cfX,cix_orb,cixdim,iSx, &
 iorbital,nindo,norbitals,maxdim2,verbose)
use test_rotations
integer :: ig,norbitals,maxdim2,cix_orb(norbitals),cixdim(2), &
 iSx(maxdim2,norbitals),iorbital(2,1),nindo(norbitals)
complex*16 :: RT(3,3,2),cfX(maxdim2,maxdim2,norbitals,norbitals)
logical :: verbose
RT=matrices(:,:,:,ig)
end subroutine
subroutine FindUniqueCix(ng,sizes,groups)
integer,intent(out) :: ng,sizes(2),groups(2,2)
ng=1; sizes=0; groups=0
sizes(1)=2; groups(1,1)=1; groups(1,2)=2
end subroutine
'''
source=stubs
for name in ['SymmetrizeOverlap','SymmetrizeNormalizationAdjoint','AverageEquivalentCixMatrices']:
    source+=routine('symmetrize.f90',name)
source+=routine('l2mainSubs.f90','BuildNormalizationAdjoint')
source+=routine('lomain.f90','ProjectorGlobalAxes')
source+=routine('inverss.f90','inv_3x3')
(OUT/'symmetry_probe.f90').write_text(source)
compiler=os.environ.get('FC',shutil.which('gfortran') or 'gfortran')
library=OUT/('symmetry_probe.dylib' if sys.platform=='darwin' else 'symmetry_probe.so')
link_flags=['-framework','Accelerate'] if sys.platform=='darwin' else ['-llapack','-lblas']
subprocess.run([compiler,'-O0','-g','-fcheck=all','-ffree-line-length-none',
 '-shared','-fPIC',str(OUT/'symmetry_probe.f90'),*link_flags,
 '-o',str(library)],cwd=OUT,check=True)
lib=ct.CDLL(str(library))
F=lambda x:np.asfortranarray(x)
I=lambda x:np.asfortranarray(x,dtype=np.int32)
Z=lambda s:np.zeros(s,complex,order='F')
def call(name,*args):
    getattr(lib,name.lower()+'_')(*[x.ctypes.data_as(ct.c_void_p) for x in args])
def herm(a):return (a+a.conj().T)/2
def invsqrt(a):
    e,v=eigh(a)
    return (v/np.sqrt(e))@v.conj().T
rng=np.random.default_rng(71137)
def rand(*shape):return rng.normal(size=shape)+1j*rng.normal(size=shape)
rot=Z((3,3,2,2))
for b in range(2):
    v,_=np.linalg.qr(rand(3,3))
    rot[:,:,b,0]=np.eye(3)
    rot[:,:,b,1]=v@np.diag([1.,-1.,-1.])@v.conj().T
call('set_test_rotations',rot)
args=[Z((3,3,2,2)),I([1,2]),I([3,3]),I(np.tile(np.arange(1,4)[:,None],(1,2))),
 I([[1],[2]]),I([3,3]),I([2]),I([3])]
def forward(a):
    r=Z((3,3,2))
    for b in range(2):
        r[:,:,b]=sum(rot[:,:,b,g].conj().T@a[:,:,b]@rot[:,:,b,g] for g in range(2))/2
    r[:]=r.mean(axis=2)[:,:,None]
    return r
raw=Z((3,3,2));delta=Z((3,3,2));j=Z((3,3,2))
for b in range(2):
    t=rand(3,3);raw[:,:,b]=t.conj().T@t+np.eye(3)
    delta[:,:,b]=herm(rand(3,3));j[:,:,b]=herm(rand(3,3))
actual=raw.copy(order='F');call('symmetrizeoverlap',actual,*args)
assert np.max(abs(actual-forward(raw)))<1e-13
x=Z((3,3,2));moment=Z((3,3,2))
for b in range(2):
    x[:,:,b]=invsqrt(actual[:,:,b]);moment[:,:,b]=(j[:,:,b]@x[:,:,b]).T
qbar=Z((3,3,2))
call('buildnormalizationadjoint',qbar,moment,x,I(np.ones((3,3,2))),I([3,3]),I([3]),I([2]))
base=qbar.copy(order='F')
saved=os.environ.pop('EDMFT_STRESS_NORM_SYMMETRY',None)
call('symmetrizenormalizationadjoint',qbar,*args)
def value(t):
    s=forward(raw+t*delta)
    return sum(2*np.trace(j[:,:,b]@invsqrt(s[:,:,b])).real for b in range(2))
h=1e-5;fd=(value(h)-value(-h))/(2*h)
predicted=2*np.sum(qbar*delta).real
without=2*np.sum(base*delta).real
assert abs(fd-predicted)<2e-9,(fd,predicted)
assert abs(fd-without)>1e-4
for mode in ['OFF','DIAG','APPLY']:
    os.environ['EDMFT_STRESS_NORM_SYMMETRY']=mode
    result=base.copy(order='F');call('symmetrizenormalizationadjoint',result,*args)
    assert np.max(abs(result-(qbar if mode=='APPLY' else base)))<1e-14
if saved is None:os.environ.pop('EDMFT_STRESS_NORM_SYMMETRY',None)
else:os.environ['EDMFT_STRESS_NORM_SYMMETRY']=saved

# Fractional and Cartesian rotations are different objects. Fixed physical
# local axes require L=C B R_fractional B^-1, even for a nonorthogonal cell.
br=np.array([[1.3,.2,.1],[0.,.9,.17],[0.,0.,.71]])
c,_=np.linalg.qr(rng.normal(size=(3,3)))
r=np.diag([1.,-1.,-1.]);rf=np.linalg.solve(br,r@br)
axes=np.zeros((3,3),order='F')
call('projectorglobalaxes',axes,F(c),F(br),F(rf))
assert np.max(abs(axes-c@r))<1e-14
kf=np.array([.37,-.29,.61]);k=br@kf
errors=[];legacy=[]
for a,b in [(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]:
    e=np.zeros((3,3));e[a,b]+=.5;e[b,a]+=.5
    def local_vector(t):
        bt=(np.eye(3)-t*e)@br
        rt=bt@rf@np.linalg.inv(bt)
        ct=axes@np.linalg.inv(rt)
        return ct@bt@rf@kf
    fdq=(local_vector(1e-6)-local_vector(-1e-6))/2e-6
    exact=-axes@e@k
    errors.append(float(np.linalg.norm(fdq-exact)))
    legacy.append(float(np.linalg.norm((-c@e@r@k)-exact)))
assert max(errors)<5e-10
assert min(legacy[-2:])>.1
report={'normalization_functional_fd':fd,'production_adjoint':float(predicted),
 'without_symmetry_adjoint':float(without),'absolute_error':float(abs(fd-predicted)),
 'global_frame_fd_errors':errors,'legacy_frame_errors':legacy,
 'source_hashes':{f:hashlib.sha256((ROOT/'src/dmft2'/f).read_bytes()).hexdigest()
                  for f in ['symmetrize.f90','l2mainSubs.f90','lomain.f90','inverss.f90']}}
(OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
print('PASS: exact symmetry adjoint, default/diagnostic modes, fixed physical axes')
