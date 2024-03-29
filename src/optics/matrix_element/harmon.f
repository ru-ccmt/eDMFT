      SUBROUTINE HARMON(N,X,Y,Z,LMAX2,F,DF,RI)                           
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),Y(N),Z(N),F(LMAX2+1,N),DF(LMAX2+1,N)                           
      DIMENSION A(3)                                                    
      COMMON /GENER/  BR1(3,3)                                          
      LMX=LMAX2+1                                                        
      DO 1 I=1,N                                                        
      A(1)=X(I)*BR1(1,1)+Y(I)*BR1(1,2)+Z(I)*BR1(1,3)                    
      A(2)=X(I)*BR1(2,1)+Y(I)*BR1(2,2)+Z(I)*BR1(2,3)                    
      A(3)=X(I)*BR1(3,1)+Y(I)*BR1(3,2)+Z(I)*BR1(3,3)                    
      XM=SQRT(A(1)**2+A(2)**2+A(3)**2)                                  
      XA=RI*XM                                                          
      CALL SPHBES(LMAX2,XA,F(1,I))                                       
      CALL DVBES1(F(1,I),DF(1,I),XA,RI,LMX)                             
      DO 2 J=1,LMX                                                      
    2 DF(J,I)=XM*DF(J,I)                                                
    1 CONTINUE                                                          
      RETURN                                                            
      END                                                               
