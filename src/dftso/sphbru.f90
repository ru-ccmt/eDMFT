REAL*8 FUNCTION SPHBRU(N,X)
  !IMPLICIT REAL*8(A-H,O-Z)
  IMPLICIT NONE
  INTEGER, intent(in) :: n
  REAL*8,  intent(in) :: x
  ! locals
  REAL*8, PARAMETER :: TOL = 1.D-14
  REAL*8  :: BB, DBB, FAC, FACD
  INTEGER :: i, j
  I=0
  BB=1.D0
  DO
     I=I+1
     FAC=1.D0
     DO J=1,I
        FAC=FAC/DFLOAT(J)
     ENDDO
     FACD=1.D0
     DO J=1,I
        FACD=FACD/(2.D0*DFLOAT(N)+3.D0+2.D0*DFLOAT(J-1))
     ENDDO
     DBB=(X*X/2.D0)**I*FAC*FACD
     BB=BB+DBB*(-1.D0)**I
     IF(DBB.LT.TOL) EXIT
  ENDDO
  SPHBRU=BB
  RETURN
END FUNCTION SPHBRU
