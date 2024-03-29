SUBROUTINE DVBES1(FJ,DJ,SM,RI,NT)
  !
  IMPLICIT NONE
  !
  !        Arguments
  !
  INTEGER            NT
  DOUBLE PRECISION   RI, SM
  DOUBLE PRECISION   DJ(NT), FJ(NT)
  !
  !     ..................................................................
  !
  !        calculate the derivatives of the Bessel functions
  !           DJ = DFJ/DX where X = SM*RI
  !                                             D.D.KOELLING
  !
  !     ..................................................................
  !
  !        Local Parameters
  !
  DOUBLE PRECISION   ONE, THIRD, ZUP, ZERO
  PARAMETER          (ONE   = 1.0D+0)
  PARAMETER          (THIRD = 1.0D+0/3.0D+0)
  PARAMETER          (ZERO  = 0.0D+0)
  PARAMETER          (ZUP   = 1.0D-5)
  !
  !        Local Scalars
  !
  INTEGER            L, LM
  DOUBLE PRECISION   Q2, Q3, X
  !
  X = SM*RI
  IF (X .GT. ZUP) THEN
     Q2 = -ONE/X
     Q3 = Q2
     DJ(1) = -FJ(2)
     LM = 1
     DO L = 2, NT
        Q3 = Q3 + Q2
        DJ(L) = FJ(LM) + Q3*FJ(L)
        LM = LM + 1
     ENDDO
  ELSE
     DJ(1) = ZERO
     DJ(2) = THIRD
     DO L = 3, NT
        DJ(L) = ZERO
     ENDDO
  ENDIF
  !
  RETURN
END SUBROUTINE DVBES1
