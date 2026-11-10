from pathlib import Path
import subprocess,shutil,os
r=Path(__file__).resolve().parent;d=r/'elect_fixture';d.mkdir(exist_ok=True);src=r/'SRC_lapw0'
(d/'stubs.f90').write_text('''module TestCases
 logical :: SimpInt=.false.,ConOpt=.true.
end module
module struct
 integer :: nat=1,ndif=1,iatnr(1)=-1,mult(1)=1,jri(1)=881
 real(8) :: zz(1)=0,vol=1,rmt(1)=1,r0(1)=1d-5,dx(1),rotij(3,3,1),rotloc(3,3,1)
 character(4) :: irel='NREL',lattic='P'
end module
module densit
 integer :: lmmaxx=264
 real(8) :: clm(881,264,1),clmsp(881,264,1,2)
end module
module work
 real(8) :: v(881,264,1),exc(881,264,1,2),vxc(881,264,1,2)
end module
module gga_str_kb
 real(8) :: gga_fin_str(9)=0
end module
module parallel
 integer :: myid=0,nstart(0:0)=1,nstop(0:0)=1
end module
module norm_kub
 real(8) :: c_kub(0:10,0:10)=0
end module
''')
(d/'main.f90').write_text('''program check
 use struct
 use densit
 use work
 implicit none
 integer :: lm(2,264,1),lmmax(1),mode,side,i,l,m
 real(8) :: r(881),xc(9),pi
 complex(8) :: inter(9)
 pi=acos(-1d0);dx=log(1d5)/880
 do i=1,881
 r(i)=r0(1)*exp((i-1)*dx(1))
 enddo
 rotloc=0;rotij=0
 do i=1,3
 rotloc(i,i,1)=1;rotij(i,i,1)=1
 enddo
 do side=1,3
 do mode=0,5
  lm=0;clm=0;v=0;exc=0;vxc=0;xc=0;inter=0
  clm(:,1,1)=4*pi*r*r
  v(:,1,1)=sqrt(4*pi)
  lmmax=1
  if(mode>0)then
   select case(mode)
   case(1);l=2;m=0
   case(2);l=2;m=1
   case(3);l=2;m=2
   case(4);l=-2;m=1
   case(5);l=-2;m=2
   end select
   lm(:,2,1)=[l,m];lm(:,3,1)=[l,-m];lmmax=2
   if(side==1.or.side==3)then
    v=0;v(:,2,1)=r*r
   endif
   if(side==2.or.side==3)then
    clm=0;clm(:,2,1)=r**4
   endif
  endif
  write(21,'(a,2i4)') 'CASE',side,mode
  call finl_elect_str(lm,lmmax,inter,0d0,xc)
 enddo
 enddo
end program
''')
shutil.copy2(src/'param.inc',d/'param.inc')
# Link the actual contraction and real-harmonic conversion, with identity geometry.
for mode,path in [('before',r/'before/lapw0/finl_elect_str.F'),('fixed',src/'finl_elect_str.F')]:
 code=path.read_text();(d/f'finl_{mode}.F').write_text(code)
 cmd=['/opt/homebrew/bin/gfortran','-O0','-g','-fcheck=all','-ffree-form','-ffree-line-length-none','-I.','stubs.f90',f'finl_{mode}.F','main.f90']+[str(src/(f+'.f')) for f in ['multfc','charge','spline','c_alpha_m','gaunt','notri','t3j','t3j0','ph']]+['-o',f'check_{mode}']
 subprocess.run(cmd,cwd=d,check=True,stdout=subprocess.DEVNULL)
 with (d/f'log_{mode}').open('w') as f:subprocess.run([str(d/f'check_{mode}')],cwd=d,env=os.environ|{'WIEN_STRESS_DEBUG':'1'},check=True,stdout=f,stderr=subprocess.STDOUT)
 shutil.copy2(d/'fort.21',d/f'output_{mode}.txt')
print('ELECT FIXTURE DONE')
