SUBROUTINE correct_rotij(rotij, lattic,BR2,ORTHO, mult,nat,ndif)
  IMPLICIT NONE
  REAL*8, intent(inout):: rotij(3,3,ndif)
  REAL*8, intent(in)   :: BR2(3,3)   ! BR2(:,i) is primitive reciprocal vector b_i
  LOGICAL, intent(in)  :: ORTHO
  CHARACTER*4, intent(in) :: LATTIC
  INTEGER, intent(in)  :: mult(nat), nat, ndif
  ! locals
  INTEGER :: index, jatom, j
  REAL*8  :: RBAS(3,3), GBAS(3,3), tmp(3,3)
  REAL*8  :: det, PI
  
  PI=ACOS(-1.0D0)
  
  IF(.NOT. ORTHO.and.lattic(1:3).ne.'CXZ')  then !     for  mon.CXZ type rotation skipped (caution!!!)
     ! BR2 is primitive reciprocal vector
     GBAS(:,:) = BR2(:,:)/(2.0d0*PI)

     !! RBAS becomes primitive direct vectors
     tmp(:,:) = transpose(GBAS)
     call INVERSSYMDEF(tmp,RBAS)
     RBAS(:,:) = transpose(RBAS)
     
     index = 0
     DO jatom=1,nat
        DO j=1,mult(jatom)
           index = index + 1
           !  in dmft1 we used : rotij_cartesian = BR1 * rotij * BR1inv^{-1}
           !  but here we use  : rotij_cartesian = BR2 * rotij * ((BR2^T)^{-1})^T
           tmp = matmul(GBAS, rotij(:,:,index))
           rotij(:,:,index) = matmul(tmp, RBAS)
        enddo
     enddo
  endif
  
END SUBROUTINE correct_rotij
