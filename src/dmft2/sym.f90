! Here is the right way of doing this:
!    Under symmetry operation Gamma, the Bloch wave function transforms as Gamma psi_{k,s} = psi_{Gamma k, Gamma s}, which 
!    not just gives another momentum point in the star, but also mixes spin-up with spin-down. For arbitrary group operation
!    you should use spin_rotate module, and calculate R_spin(Gamma) ==>  R2 = spin_rotation(angles(R_{3x3}))  
!    Then you should use R2 to rotate spin, namely Gamma Psi_{k,s} = psi_{Gamma k, R2 s}
!
!    The implementation below is a limited case of this, and takes into account only two case,
!       (A) Rz(f)=[[cos(f),sin(f),0],[-sin(f),cos(f),0],[0,0,+-1]]. This can also be multiplied by inversion, which does not do anything to the spin.
!           In this case R2=[[e^{i*f/2},0],[0,e^{-i*f/2}]] and hence up remains up and down remains down.
!       (B) Rx(pi)*Rz(f)=[[cos(f),sin(f),0],[sin(f),-cos(f),0],[0,0,+-1]].
!           In this case R2=Rx(pi)*Rz(f) = [[0,i*e^{i*(pi-f)/2}],[0,e^{-i(pi-f)/2}]], hence we have similar phase, but also up->down and down->up.
!           But we also know that u_{up} = u^*_{down} and u_{down}=-u^*_{up}, hence with idet and conjugation, we can work with single spin species.
!
!..... program sym is used for calculation of transformation parameters
!..... of so wavefunctions. It transforms the symmetry operations
!..... to frame of reference z||M. Variable det=1 ... M not reversed
!..... det=-1 ... M reversed ( must by later coupled to time  inversion ).
!.... Angle phase(I) is used for spin part of transformation matrix 
!..... acting of spinor function.
SUBROUTINE SYM(THETA,FI)
  USE param
  USE com_mpi
  USE struct
  USE sym2
  !IMPLICIT REAL*8 (A-H,O-Z)
  IMPLICIT NONE
  COMMON /GENER/  BR1(3,3),BR2(3,3)
  DIMENSION rot(3,3),trans(3,3),rotinv(3,3)
  pi=acos(-1.d0)

  !call INVERSSYMDEF(BR1,BR1inv)
  call inv_3x3(BR1,BR1inv)
  
  !CALL symoper
  ! local coordinate system is : rot=(e_theta, e_phi, e_r)
  ct=cos(theta)
  cf=cos(fi)
  st=sin(theta)
  sf=sin(fi)
  ! e_theta
  rot(1,1)=cf*ct
  rot(2,1)=sf*ct
  rot(3,1)=-st
  ! e_phi
  rot(1,2)=-sf
  rot(2,2)=cf
  rot(3,2)=0
  ! e_r
  rot(1,3)=st*cf
  rot(2,3)=st*sf
  rot(3,3)=ct
  
  if (myrank.EQ.master .OR. fastFilesystem) then
     write(6,*)'SPIN COORD SYSTEM:'
     write(6,5)((rot(i,j),j=1,3),i=1,3)
     write(6,*)
  endif
  !call INVERSSYMDEF(rot,rotinv)
  rotinv = transpose(rot)
  iordnew=0
  do  i=1,IORD   ! 100
     ! transforms symmetry matrices to cartesian system
     !   for .ortho. systems   opimat == iz, otherwise opimat == BR1 . iz . BR1^-1 
     opimat = matmul( matmul( BR1, iz(:,:,i)), BR1inv )
     
     if (myrank.EQ.master .OR. fastFilesystem) then
        write(6,*)'in global cartezian coordinates:'
        !!! Symmetry operation in global cartesian coordinates
        write(6,5)((opimat(ii,jj),jj=1,3),ii=1,3)
        write(6,*)'____________________________'
     endif
     !!! Symmetry operation in the coordianate system in which magnetization points in z-direction
     !!! Magnetization direction is specified in case.inso
     !call transform(trans,rotinv,opimat(1,1,i),rot)
     
     trans = matmul( matmul(rotinv, opimat, rot)
     
     if (myrank.EQ.master .OR. fastFilesystem) then
        write(6,*)'in spin coordinates:'
        write(6,5)((trans(ii,jj),jj=1,3),ii=1,3)
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
     
     if (myrank.EQ.master .OR. fastFilesystem) write(6,*)i,phase(i),det
     if (myrank.EQ.master .OR. fastFilesystem) write(6,*)
     if (det.gt.0.) then
        idet(i)=1
     else
        idet(i)=-1
     end if
     if (abs(1.-abs(det)).gt.1d-4) then
        if (myrank.EQ.master .OR. fastFilesystem) write(6,*)'symm. operation ',i,' removed !!!! so-det=',det
!!!	   STOP
     else
        iordnew=iordnew+1
        idet(iordnew)=idet(i)
        phase(iordnew)=phase(i)
        opimat(1:3,1:3,iordnew)=opimat(1:3,1:3,i)
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

!SUBROUTINE INVERSSYMDEF(A,AINV)
!  IMPLICIT NONE
!  !IMPLICIT REAL*8 (A-H,O-Z)
!  REAL*8, intent(in)  :: A(3,3)
!  REAL*8, intent(out) :: AINV(3,3)
!  ! locals
!  REAL*8 :: det
!  det= a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3)-a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
!  AINV(1,1) =(   A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / det
!  AINV(2,1) =( - A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / det
!  AINV(3,1) =(   A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / det
!  AINV(1,2) =( - A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / det
!  AINV(2,2) =(   A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / det
!  AINV(3,2) =( - A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / det
!  AINV(1,3) =(   A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / det
!  AINV(2,3) =( - A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / det
!  AINV(3,3) =(   A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / det
!  RETURN
!END SUBROUTINE INVERSSYMDEF
