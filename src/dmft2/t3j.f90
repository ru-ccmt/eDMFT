REAL*8 FUNCTION T3J (J1,J2,J3,M1,M2,M3)                           
  !                                                                       
  !     ANGULAR MOMENTA ENTER WITH 2X ACTUAL VALUES                       
  !     ARRAY FCT MUST BE DECLARED AND FILLED GLOBALLY WITH               
  !     J-FACTORIAL IN THE 2J+1 POSITION                                  
  !     LARGEST  FACTORIAL NEEDED IS JA+JB+JC+1                           
  !     FUNCTION NOTRI(J,K,L) RETURNS 1 IF K,L,J FORM A TRIANGLE,-1 IF NOT
  !     FUNCTION PH IS NEEDED TO COMPUTE PHASES                           
  !     NOTRI AND M-TEST IN MAIN PROGRAM                                  
  !                                                                       
  use fact, only : FCT
  IMPLICIT REAL*8 (A-H,O-Z)
  !COMMON/FACT/FCT(100)                                              
  T3J=0.                                                            
  K0=J1+J2-J3                                                       
  K1=J1-M1                                                          
  K2=J2+M2                                                          
  L1=J2-J3-M1                                                       
  L2=J1-J3+M2                                                       
  KMAX=MIN0(K1,K2,K0)+1                                             
  KMIN=MAX0(0,L1,L2)+1                                              
  C=SQRT(FCT(K0+1)*FCT(J1-J2+J3+1)*FCT(J2-J1+J3+1)/FCT(J1+J2+J3+3)*FCT(J1+M1+1)*FCT(K1+1)*FCT(J2-M2+1)*FCT(K2+1)*FCT(J3+M3+1)*FCT(J3-M3+1))                                             
  DO I=KMIN,KMAX,2 ! 10
     K=I-1                                                             
     T3J=T3J+C*PH(K)/FCT(I)/FCT(K0-K+1)/FCT(K1-K+1)/FCT(K2-K+1)/FCT(I-L1)/FCT(I-L2)
  ENDDO
     
  T3J=PH(J1-J2-M3)*T3J                                              
  RETURN                                                            
END FUNCTION T3J

REAL*8 FUNCTION T3J0(J1,J2,J3)                                    
  use fact, only: FCT
  IMPLICIT REAL*8 (A-H,O-Z)
  !COMMON/FACT/FCT(100)                                              
  J=J1+J2+J3                                                        
  F1=SQRT(FCT(J-2*J1+1)*FCT(J-2*J2+1)*FCT(J-2*J3+1)/FCT(J+3))       
  F2=FCT(J/2+1)/FCT(J/2-J1+1)/FCT(J/2-J2+1)/FCT(J/2-J3+1)           
  T3J0=PH(J/2)*F1*F2                                                
  RETURN                                                            
END FUNCTION T3J0

REAL*8 FUNCTION PH(N)                                             
  IMPLICIT REAL*8 (A-H,O-Z)
  PH=1.0D0                                                          
  K=MOD(IABS(N),4)                                                  
  IF(K.EQ.0) GOTO 5                                                 
  PH=-1.0D0                                                         
5 RETURN                                                            
END FUNCTION PH
