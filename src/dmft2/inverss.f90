SUBROUTINE inv_3x3(A,Ainv)
  IMPLICIT NONE
  REAL*8, intent(in)  :: A(3,3)
  REAL*8, intent(out) :: Ainv(3,3)
  ! locals
  REAL*8 :: det
  det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) &
      - A(3,1)*A(2,2)*A(1,3) - A(1,1)*A(3,2)*A(2,3) - A(2,1)*A(1,2)*A(3,3)
  Ainv(1,1) =(   A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / det
  Ainv(1,2) =( - A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / det
  Ainv(1,3) =(   A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / det
  Ainv(2,1) =( - A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / det
  Ainv(2,2) =(   A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / det
  Ainv(2,3) =( - A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / det
  Ainv(3,1) =(   A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / det
  Ainv(3,2) =( - A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / det
  Ainv(3,3) =(   A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / det
  RETURN
END SUBROUTINE inv_3x3
