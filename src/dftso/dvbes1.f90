SUBROUTINE DVBES1(FJ,DJ,SM,RI,NT)
  IMPLICIT NONE
  !        Arguments
  INTEGER, intent(in) :: NT
  REAL*8, intent(in)  :: RI, SM, FJ(NT)
  REAL*8, intent(out) :: DJ(NT)
  !..................................................................
  !   calculate the derivatives of the Bessel functions
  !      DJ = DFJ/DX where X = SM*RI
  !                                        D.D.KOELLING
  !..................................................................
  !   Local Parameters
  REAL*8, PARAMETER :: ZUP = 1.0D-5
  !        Local Scalars
  INTEGER            L, LM
  DOUBLE PRECISION   Q2, Q3, X
  !
  X = SM*RI
  IF (X .GT. ZUP) THEN
     Q2 = -1.d0/X
     Q3 = Q2
     DJ(1) = -FJ(2)
     LM = 1
     DO L = 2, NT
        Q3 = Q3 + Q2
        DJ(L) = FJ(LM) + Q3*FJ(L)
        LM = LM + 1
     ENDDO
  ELSE
     DJ(1) = 0.d0
     DJ(2) = 1.0D+0/3.0D+0
     DO L = 3, NT
        DJ(L) = 0.d0
     ENDDO
  ENDIF
  !
  RETURN
END SUBROUTINE DVBES1
