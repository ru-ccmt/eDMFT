      SUBROUTINE DVBES1(FJ,DJ,SM,RI,NT)                                 
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
!-----X CALCULATE THE DERIVATIVES OF THE BESSEL FUNCTIONS.   X----X----X
!-----X   DJ=DFJ/DX WHERE X=SM*RI                                 X----X
!-----X                    D.D.KOELLING                      X----X----X
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DJ(*),FJ(*)                                             
      DATA ZUP/1.0D-5/,ZERO/0.0D0/,THIRD/0.3333333333333D0/,ONE/1.0D0/  
      X=SM                                                              
      IF(X.GT.ZUP) GOTO 20                                              
      DJ(1)=ZERO                                                        
      DJ(2)=THIRD                                                       
      DO 10 L=3,NT                                                      
   10 DJ(L)=ZERO                                                        
      RETURN                                                            
   20 Q2=-ONE/X                                                         
      Q3=Q2                                                             
      DJ(1)=-FJ(2)                                                      
      LM=1                                                              
      DO 30 L=2,NT                                                      
      Q3=Q3+Q2                                                          
      DJ(L)=FJ(LM)+Q3*FJ(L)                                             
      LM=LM+1                                                           
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
