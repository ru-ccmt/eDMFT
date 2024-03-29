SUBROUTINE CALI(AT,SUM,RX,DX,JWI)
  USE param, ONLY: nrad
  IMPLICIT NONE
  REAL*8, intent(in) :: AT(NRAD), RX(NRAD), DX
  REAL*8, intent(out):: SUM
  ! locals
  INTEGER :: jwi, i, jf, js
  REAL*8  :: sum1
  !---------------------------------------------------
  !DIMENSION AT(NRAD),RX(NRAD)
  !*****************************************************
  SUM=0.0D0
  SUM1=0.0D0
  JS=JWI-2*(JWI/2)+1
  JF=JWI-1
  DO I=JS,JF,2
     SUM1=SUM1+(2.D0*RX(I)*AT(I)+RX(I+1)*AT(I+1))
  ENDDO
  IF(JS.EQ.1) THEN
     SUM=(SUM1+SUM1-RX(JWI)*AT(JWI))*DX/3.D0
  ELSE IF(JS.EQ.2) THEN
     SUM=(SUM1+SUM1-RX(JWI)*AT(JWI)+RX(1)*AT(1))*DX/3.D0
  END IF
  RETURN
END SUBROUTINE CALI
