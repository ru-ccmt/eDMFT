from pathlib import Path
import re, subprocess, hashlib, json

ROOT=Path(__file__).resolve().parents[3]
import os,tempfile
OUT=Path(os.environ.get('EDMFT_STRESS_TEST_DIR',tempfile.mkdtemp(prefix='edmft_stress_test_')))
OUT.mkdir(parents=True,exist_ok=True)


def routine(file,name):
    text=(ROOT/file).read_text()
    pat=rf'^\s*SUBROUTINE\s+{name}\s*\(.*?^\s*END SUBROUTINE\s+{name}\s*$'
    match=re.search(pat,text,re.M|re.S|re.I)
    if not match: raise RuntimeError(name)
    return match.group(0)+'\n'

stubs='''module fact
real*8 :: fct(100)
end module
module com_mpi
logical :: Qprint=.false.
integer :: myrank=1,master=0
end module
module dmft_stress
real*8 :: stress_norm_response_sign=1d0
end module
module param
integer,parameter :: lmax2n=8
end module
module defs
complex*16,parameter :: zeroc=(0d0,0d0)
end module
module w_atpar
type cache_t
 logical :: ready=.false.
 integer :: nterm=0
 integer,allocatable :: ih(:),ibra(:),iket(:)
 complex*16,allocatable :: radial(:,:,:),angular(:,:,:)
end type
type(cache_t) :: nsp_geometry(1)
complex*16 :: VPts(2,2,1,1),WPts(2,2,1,1)
! Match production one-based bounds, including the first-record test.
integer :: lpv(1,1)=0,lv(1,1)=0,llv(1,1)=0,mmv(1,1)=0,Nhmx(1)
end module
'''
names=['BuildWindowResponseMatrices','ComputeNormalizationDerivative',
       'BuildNormalizationAdjoint','BuildNormalizationWindowResponseMatrices']
sources=stubs+'\n'.join(routine('src/dmft2/l2mainSubs.f90',n) for n in names)
for n in ['AccumulateDirectProjectorPair','AccumulateDirectProjectorSpectral']:
    sources+=routine('src/dmft2/cmp_dmft_weights.f90',n)
for n in ['c_alpha_m','C_factor','build_nsp_geometry']:
    sources+=routine('src/dmft2/stress_kernels.f90',n)
sources+='''
subroutine init_audit_fact()
use fact
implicit none
integer :: i
fct=0d0
fct(1)=1d0
do i=3,99,2
 fct(i)=fct(i-2)*dble((i-1)/2)
enddo
end subroutine
subroutine angular_probe(lb,mb,lk,mk,lp,mp,rad,ang)
use w_atpar
implicit none
integer,intent(in) :: lb,mb,lk,mk,lp,mp
complex*16,intent(out) :: rad(3,3),ang(3,3)
integer :: it
if (allocated(nsp_geometry(1)%ih)) deallocate(nsp_geometry(1)%ih, &
 nsp_geometry(1)%ibra,nsp_geometry(1)%iket,nsp_geometry(1)%radial, &
 nsp_geometry(1)%angular)
nsp_geometry(1)%ready=.false.
lpv(1,1)=lb;lv(1,1)=lk;llv(1,1)=lp;mmv(1,1)=mp;Nhmx(1)=1
VPts=1d0;WPts=1d0
call build_nsp_geometry(1)
rad=0d0;ang=0d0
do it=1,nsp_geometry(1)%nterm
 if(nsp_geometry(1)%ibra(it)/=lb*(lb+1)+mb+1) cycle
 if(nsp_geometry(1)%iket(it)/=lk*(lk+1)+mk+1) cycle
 rad=nsp_geometry(1)%radial(:,:,it)
 ang=nsp_geometry(1)%angular(:,:,it)
enddo
end subroutine
'''
g2=(ROOT/'src/dmft2/stress_g2_matrices.f90').read_text()
newmodule='''module g2_angular_test_module
use iso_fortran_env,only: dp=>real64
implicit none
integer,parameter :: axis_a(6)=[1,2,3,2,1,1],axis_b(6)=[1,2,3,3,3,2]
complex(dp),parameter :: one=(1._dp,0._dp),zero=(0._dp,0._dp)
contains
'''
for name in ['gauss_legendre_stress','g2_angular_tables']:
    match=re.search(rf'^  subroutine {name}\(.*?^  end subroutine$',g2,re.M|re.S|re.I)
    if not match:raise RuntimeError(name)
    newmodule+=match.group(0)+'\n'
sources=newmodule+'end module\n'+sources
sources+='''
subroutine quadrature_probe(lb,mb,lk,mk,lp,mp,rad,ang)
use g2_angular_test_module
implicit none
integer,intent(in) :: lb,mb,lk,mk,lp,mp
complex(dp),intent(out) :: rad(3,3),ang(3,3)
complex(dp),allocatable :: ag(:,:),ar(:,:,:),aa(:,:,:)
integer :: lmax,nh,i,j,q,a,b
lmax=max(lb,lk);nh=(lmax+1)**2
allocate(ag(nh,nh),ar(nh,nh,6),aa(nh,nh,6))
call g2_angular_tables(lmax,lp,mp,one,zero,ag,ar,aa)
i=lb*(lb+1)+mb+1;j=lk*(lk+1)+mk+1
rad=zero;ang=zero
do q=1,6
 a=axis_a(q);b=axis_b(q)
 rad(a,b)=ar(i,j,q);rad(b,a)=ar(i,j,q)
 ang(a,b)=aa(i,j,q);ang(b,a)=aa(i,j,q)
enddo
deallocate(ag,ar,aa)
end subroutine
'''
(OUT/'extracted_routines.f90').write_text(sources)
files=[ROOT/'src/dmft2'/x for x in ['gaunt.f90','t3j.f90','notri.f90','ylm.f90']]
import shutil,sys
compiler=os.environ.get('FC',shutil.which('gfortran') or 'gfortran')
cmd=[compiler,'-O0','-g','-fcheck=all','-ffree-line-length-none',
     '-shared','-fPIC',str(OUT/'extracted_routines.f90'),*map(str,files),
     *(['-framework','Accelerate'] if sys.platform=='darwin' else ['-llapack','-lblas']),'-o',str(OUT/'audit_routines.dylib')]
subprocess.run(cmd,cwd=OUT,check=True)
manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
          [ROOT/'src/dmft2/l2mainSubs.f90',ROOT/'src/dmft2/stress_kernels.f90',ROOT/'src/dmft2/stress_g2_matrices.f90',ROOT/'src/dmft2/cmp_dmft_weights.f90',*files]}
(OUT/'source_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Compiled current production routines with bounds checking in',OUT)
