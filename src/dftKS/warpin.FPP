#ifdef _REAL
REAL*8     FUNCTION WARPIN(X,Y,Z)
#else  
COMPLEX*16 FUNCTION WARPIN(X,Y,Z)
#endif    
  use out, only      : kmax, WARP
  IMPLICIT NONE
  INTEGER            X,  Y,  Z
  !
  !INCLUDE 'param.inc'
  IF ((ABS(X) .LE. kmax(1)) .AND. (ABS(Y) .LE. kmax(2)) .AND. (ABS(Z) .LE. kmax(3))) THEN
     WARPIN = WARP(X,Y,Z)
  ELSE
     WARPIN = 0.d0
  ENDIF
  RETURN
END FUNCTION WARPIN
