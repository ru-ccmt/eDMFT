SUBROUTINE CHARGE(RX,S,FPIV,FUNCT,DELTA,JEND,SUM)                 
  USE param, ONLY: NRAD
  IMPLICIT NONE
  REAL*8, intent(out):: SUM
  REAL*8, intent(in) :: RX(NRAD), FUNCT(NRAD)
  REAL*8, intent(in) :: S, FPIV, DELTA
  INTEGER, intent(in):: JEND
  ! locals
  INTEGER :: NPOINT, JS, JF, J
  
  SUM=0.0D0                                         
  NPOINT=JEND                                                       
  JS=2-MOD(NPOINT,2)                                            
  JF=JEND-2                                                     
  DO J=JS,JF,2                                                    
     SUM=SUM+FUNCT(J)*RX(J)+4.*FUNCT(J+1)*RX(J+1)+FUNCT(J+2)*RX(J+2)   
  ENDDO
  SUM=SUM*DELTA/3.0D0                                               
  SUM=SUM+(1-MOD(NPOINT,2))*( FUNCT(1)+FUNCT(JS) )/ 2.0D0*( RX(JS)-RX(1))
  SUM=SUM*FPIV*S                                                    
  RETURN                                                            
END SUBROUTINE CHARGE
