SUBROUTINE symoper
  !     transforms symmetry matrices to cartesian system
  !     for .ortho. systems   iz_cartesian == iz
  !     for not.ortho. systems iz_cartesian == BR1 . iz . BR1^-1 
  USE param
  USE structure, ONLY : BR1, iz, iord
  USE sym2
  IMPLICIT NONE
  REAL*8   :: BR1in(3,3), tmp(3,3)
  INTEGER  :: i, j, k, ior
  !.........inverssymdef is a subroutine in sph-UP.frc and
  !.........calculates the inverse of an 3*3 matrix.......
  !CALL INVERSSYMDEF(BR1,BR1in)
  CALL inv_3x3(BR1,BR1in)
  !.......define symmetry matrices for cartesian system
  !
  ior=iord
  iz_cartesian(:,:,:)=0.d0
  do i=1,IORD
     tmp(:,:) = matmul(BR1(:,:),iz(:,:,i))
     iz_cartesian(:,:,i) = matmul(tmp(:,:),BR1in(:,:))
     do j=1,3
        do k=1,3
           if (abs(iz_cartesian(j,k,i)).lt.1e-6) iz_cartesian(j,k,i)=0.d0
           if (abs(iz_cartesian(j,k,i)-1.d0).lt.1e-6) iz_cartesian(j,k,i)=1.d0
           if (abs(iz_cartesian(j,k,i)+1.d0).lt.1e-6) iz_cartesian(j,k,i)=-1.d0
        end do
     end do
  end do
  
  return
!11 FORMAT(3(3I2,F10.5/))
!12 FORMAT(3(3f8.5/))
end SUBROUTINE symoper
