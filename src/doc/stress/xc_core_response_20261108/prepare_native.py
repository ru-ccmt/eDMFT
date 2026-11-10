"""Scratch-only native XC/core reference exports; no installed source edits."""
from pathlib import Path
import shutil
root=Path(__file__).resolve().parent
for name in ['SRC_lapw0','SRC_lcore']:
 shutil.copytree(Path('/Users/haule/wien2k_24')/name,root/'native'/name,dirs_exist_ok=True)

p=root/'native/SRC_lcore/resld.f';s=p.read_text()
mark=' 398  NSTOP=0'
assert s.count(mark)==1
s=s.replace(mark,""" 398  continue
      open(996,file='core_native.bin',access='stream',form='unformatted',position='append')
      write(996) nqn,nql,nk,np,imax,imat,nuc
      write(996) de,z,dvc,dfl
      write(996) dr(1:np),dv(1:np),dp(1:np),dq(1:np)
      close(996)
      NSTOP=0""")
p.write_text(s)

p=root/'native/SRC_lapw0/xcpot1.F';s=p.read_text()
mark='!     LOOP OVER ALL RADIAL MESH POINTS'
assert s.count(mark)==1
block="""      call response_xc_geometry(jatom,nrad,ncom,nspd,jri(jatom),lmmax(jatom),npt1,&
          r,yy,weight,qsin,qcos,dydt,dydt2,dydp,dydp2,dydtp,lm(:,:,jatom),GGA_SWITCH)
!
"""
s=s.replace(mark,block+mark)
# A stream of native invariants/potentials to check the new local formula.
mark='        CALL VXCLM2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU,&'
s=s.replace(mark,"""        if(GGA_SWITCH)then
          write(985) ir,k,ft,2*gmagu,2*g2u,gdggt
        else
          write(985) ir,k,ft,0.d0,0.d0,0.d0
        endif
"""+mark)
mark='!     if(ir.ge.780.and. k.lt.3) then'
s=s.replace(mark,'        write(985) vxcu\n'+mark)
s+='''
subroutine response_xc_geometry(a,nr,nc,na,n,lmmax,ng,r,yy,w,ss,cc,dt,dtt,dp,dpp,dtp,lm,gga)
implicit none
integer a,nr,nc,na,n,lmmax,ng,lm(2,nc+3)
real*8 r(nr),yy(nc,na),w(na),ss(na),cc(na),dt(nc,na),dtt(nc,na),dp(nc,na),dpp(nc,na),dtp(nc,na)
logical gga
character(80) fname
write(fname,'("xc_geometry_",i4.4,".bin")')a
open(984,file=trim(fname),access='stream',form='unformatted',status='replace')
write(984) n,lmmax,ng
write(984) lm(:,1:lmmax),r(1:n),yy(1:lmmax,1:ng),w(1:ng)
if(gga)then
write(984) ss(1:ng),cc(1:ng),dt(1:lmmax,1:ng),dtt(1:lmmax,1:ng),dp(1:lmmax,1:ng),dpp(1:lmmax,1:ng),dtp(1:lmmax,1:ng)
endif
close(984)
close(985)
write(fname,'("xc_points_",i4.4,".bin")')a
open(985,file=trim(fname),access='stream',form='unformatted',status='replace')
end subroutine
'''
p.write_text(s)
# Export raw XC harmonic potentials before adding Hartree.
p=root/'native/SRC_lapw0/lapw0.F';s=p.read_text();mark='      CALL CPUTIM(TIME4)'
assert s.count(mark)==1
s=s.replace(mark,"""      open(986,file='xc_native.bin',access='stream',form='unformatted',status='replace')
      write(986) nat,nkk,iff1,iff2,iff3,vol
      write(986) kzz(:,1:nkk)
      do jatom=1,nat
        write(986) jri(jatom),lmmax(jatom),lm(:,1:lmmax(jatom),jatom)
        write(986) vxc(1:jri(jatom),1:lmmax(jatom),jatom,1)
      enddo
      write(986) fft(:,:,:,1)
      close(986)
"""+mark)
p.write_text(s)
