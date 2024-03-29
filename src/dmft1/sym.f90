!..... program sym is used for calculation of transformation parameters
!..... of so wavefunctions. It transforms the symmetry operations
!..... to frame of reference z||M. Variable det=1 ... M not reversed
!..... det=-1 ... M reversed ( must by later coupled to time  inversion ).
!.... Angle phase(I) is used for spin part of transformation matrix 
!..... acting of spinor function.
SUBROUTINE SYM(THETA,FI)
  USE param
  USE com_mpi
  USE structure, ONLY: iz, tau
  USE sym2
  IMPLICIT REAL*8 (A-H,O-Z)
  !COMMON /GENER/  BR1(3,3),BR2(3,3)
  DIMENSION rot(3,3),trans(3,3),rotinv(3,3)
  REAL*8     :: Det, iz_proper(3,3), phi1, the1, psi1
  COMPLEX*16 :: R2(2,2)
  pi=acos(-1.d0)
  CALL symoper

!!! BRISI
  idet(:)=1
  phase(:)=0.d0
!!! BRISI

  return
  
  ct=cos(theta)
  cf=cos(fi)
  st=sin(theta)
  sf=sin(fi)
  rot(1,1)=cf*ct
  rot(1,2)=-sf
  rot(1,3)=cf*st
  rot(2,1)=sf*ct
  rot(2,2)=cf
  rot(2,3)=sf*st
  rot(3,1)=-st
  rot(3,2)=0
  rot(3,3)=ct
  if (myrank.EQ.master .OR. fastFilesystem) then
     write(6,*)'SPIN COORD SYSTEM:'
     write(6,5)((rot(i,j),j=1,3),i=1,3)
     write(6,*)
  endif
  !call INVERSSYMDEF(rot,rotinv)
  call inv_3x3(rot,rotinv)
  iordnew=0
  do  i=1,IORD   ! 100
     Det = detx(iz_cartesian(:,:,i))
     iz_proper = transpose( iz_cartesian(:,:,i)*Det )
     CALL Angles(phi1,the1,psi1, iz_proper )
     CALL Spin_Rotation(R2,phi1,the1,psi1)
     
     if (myrank.EQ.master .OR. fastFilesystem) then
        write(6,'(A,I2,A)')'-----*** Symmetry operation # ',i,' ***-----'
        write(6,*)'in global cartezian coordinates:'
        !!! Symmetry operation in global cartesian coordinates
        write(6,5)((iz_cartesian(ii,jj,i),jj=1,3),ii=1,3)
     endif
     !!! Symmetry operation in the coordianate system in which magnetization points in z-direction
     !!! Magnetization direction is specified in case.inso
     call transform(trans,rotinv,iz_cartesian(1,1,i),rot)
     if (myrank.EQ.master .OR. fastFilesystem) then
        write(6,*)'in spin coordinates:'
        write(6,5)((trans(ii,jj),jj=1,3),ii=1,3)
        write(6,*)'____________________________________'
        write(6,*)'spin rotation matrix:'
        write(6,'(2F10.6,3x,2F10.6)') (R2(1,jj),jj=1,2)
        write(6,'(2F10.6,3x,2F10.6)') (R2(2,jj),jj=1,2)
     endif
     !!! We can sort the lattice symmetry operations into three types:
     !!!    A) Do not invert magnetization, such as identity, inversion, rotation around the magnetization axis
     !!!    B) Operations which invert magnetization, such as reflection in xz or yz plane (mirror plane with magnetization in the mirror plane)
     !!!    C) Operatios which rotate magnetization in arbitrary direction are removed
     !!! Operations (A) -- diagonal in spin space, have det>0, (idet=1) and nonzero phase
     !!! Operations (B) -- equivalent to time reversal transformation, have det<0, idet=-1
     !!!
     !!! Both (A) and (B) operations in the spin coordinates must have this form:
     !!!
     !!!  (+-cos(t), sin(t),  0)
     !!!  (-+sin(t), cos(t),  0)
     !!!  (       0,     0, +-1)
     !!!  
     det=trans(1,1)*trans(2,2)-trans(1,2)*trans(2,1)
!.... calculation of the phase shift of so - functions under
!...  sym operations: for operations det=1 is the spin matrix 
!... diagonal: | exp(-i*phase/2),             0|
!.......       |               0,exp(i*phase/2)|
!.... for operations det=-1 is the spin matrix of the pure rotation
!.... off-diagonal (i.e. reverts spinor up->dn, dn->up) and must
!.... be combined with time inversion ( which reverts the spinor again )
!.... in this case the phase(i) defined in SYM includes also
!.... part of the time inversion ( the other part of time inversion
!.... - complex conjugation of wave function is included in XSPLT )

     
     DD=det*trans(3,3) !!! The full 3x3 determinant. This is because 2x2 determinant is +-unity, hence a13=0 and a23=0
     if (abs(trans(1,1)).gt.1d0) then
        a=trans(1,1)/abs(trans(1,1))
     else
        a=trans(1,1)
     end if
     
     b=trans(1,2)
     
     if (abs(b).gt.1d-8) then
        phase(i)=-DD*b/(abs(-DD*b))*acoss(DD*a)
     else
        phase(i)=acoss(DD*a)
     end if
     if (det.lt.0.) then
        phase(i)=phase(i)+pi
     end if
     if (det.gt.0.) then
        idet(i)=1
     else
        idet(i)=-1
     end if
     
     if (abs(1.-abs(det)).gt.1d-4) then
        if (myrank.EQ.master .OR. fastFilesystem) write(6,*)'symm. operation ',i,' removed !!!! so-det=',det
     else
        if (myrank.EQ.master .OR. fastFilesystem) then
           write(6,'(A,F8.5,A,F8.5,A,I2)') 'phase_up=', -phase(i)/pi/2.d0, '*pi  phase_dn=',phase(i)/pi/2.d0,'*pi  flip=', idet(i)
           write(6,*)
        endif
        iordnew=iordnew+1
        idet(iordnew)=idet(i)
        phase(iordnew)=phase(i)
        iz_cartesian(1:3,1:3,iordnew)=iz_cartesian(1:3,1:3,i)
        iz(1:3,1:3,iordnew)=iz(1:3,1:3,i)
        tau(1:3,iordnew)=tau(1:3,i)
     end if
     
  ENDDO  ! 100 !100  continue
  if ((myrank.EQ.master .OR. fastFilesystem) .and. (iord.ne.iordnew)) write(6,*) 'Original and SO-compatible symops changed:', iord, iordnew

  iord=iordnew
  
5 format(3f12.6)
END SUBROUTINE SYM
      
subroutine transform(T,Pinv,A,P)
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION T(3,3),P(3,3),A(3,3),Pinv(3,3)
  
  do i=1,3
     do j=1,3
        sum=0
        do k=1,3
           do l=1,3
              sum=sum+Pinv(i,k)*A(k,l)*P(l,j)
           end do
        end do
        T(i,j)=sum
     end do
     !       write(6,*)'Transf. matrix:',(T(i,k),k=1,3)
  end do
  return
end subroutine transform

real*8 function ACOSS(x)
  implicit real*8 (A-H,O-Z)
  
  if ((abs(x)-1.).gt.1d-4) then
     write(6,*)'x=',x
     stop 'ACOSS ERROR'
  end if
  if (x.ge.1) then
     acoss=0.0
  else if (x.le.-1.) then
     acoss=acos(-1.0)
  else
     acoss=acos(x)
  end if
end function ACOSS














