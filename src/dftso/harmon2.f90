SUBROUTINE HARMON2(N,X,Y,Z,LMAX2,F,DF,RI)
  use structure, only: BR1
  IMPLICIT NONE
  INTEGER, intent(in) :: N, LMAX2
  REAL*8,  intent(in) :: X(N), Y(N), Z(N), RI
  REAL*8, intent(out) :: F(LMAX2+1,N), DF(LMAX2+1,N)
  ! locals
  REAL*8  :: A(3), xa, xm
  INTEGER :: i, lmx
  !
  LMX=LMAX2+1
  DO I=1,N
     A(1)=X(I)*BR1(1,1)+Y(I)*BR1(1,2)+Z(I)*BR1(1,3)
     A(2)=X(I)*BR1(2,1)+Y(I)*BR1(2,2)+Z(I)*BR1(2,3)
     A(3)=X(I)*BR1(3,1)+Y(I)*BR1(3,2)+Z(I)*BR1(3,3)
     XM=SQRT(A(1)**2+A(2)**2+A(3)**2)
     XA=RI*XM
     CALL SPHBES2(LMAX2,XA,F(:,I))
     CALL DVBES12(F(:,I),DF(:,I),XA,RI,LMX)
     DF(1:LMX,I)=XM*DF(1:LMX,I)
  ENDDO
  RETURN                                                            
END SUBROUTINE HARMON2

SUBROUTINE DVBES12(FJ,DJ,SM,RI,NT)                                 
  !-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
  !-----X CALCULATE THE DERIVATIVES OF THE BESSEL FUNCTIONS.   X----X----X
  !-----X   DJ=DFJ/DX WHERE X=SM*RI                                 X----X
  !-----X                    D.D.KOELLING                      X----X----X
  !-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
  IMPLICIT NONE
  REAL*8, intent(in) :: FJ(*)                                          
  REAL*8, intent(out):: DJ(*)
  REAL*8, intent(in) :: SM
  REAL*8, intent(in) :: RI
  INTEGER, intent(in):: NT
  ! locals
  INTEGER :: L, LM
  REAL*8  :: Q2, Q3, X
  REAL*8, PARAMETER :: ZUP = 1.0D-5
  X=SM
  IF (X.LE.ZUP) THEN
     DJ(1)=0.d0
     DJ(2)=1.0/3.d0
     DO L=3,NT                                                      
        DJ(L)=0.d0
     ENDDO
     RETURN
  END IF
  Q2=-1.d0/X                                                         
  Q3=Q2                                                             
  DJ(1)=-FJ(2)                                                      
  LM=1                                                              
  DO L=2,NT                                                      
     Q3=Q3+Q2                                                          
     DJ(L)=FJ(LM)+Q3*FJ(L)                                             
     LM=LM+1                                                           
  ENDDO
  RETURN                                                            
END SUBROUTINE DVBES12

SUBROUTINE SPHBES2(N,X,FJ)                                         
  IMPLICIT NONE
  REAL*8, intent(out) :: FJ(*)
  REAL*8, intent(in)  :: X
  INTEGER, intent(in) :: N
  ! locals
  REAL*8, PARAMETER :: XLIM = 0.1D0
  REAL*8, PARAMETER :: TNHF = 10.5D0
  REAL*8, PARAMETER :: FFT = 15.0D0
  REAL*8, PARAMETER :: T25 = 1.0D25
  REAL*8, PARAMETER :: TN25 = 1.0D-25
  REAL*8, PARAMETER :: TN50 = 1.0D-50
  REAL*8  :: CUFAC, FFN, FFO, FFP, FM, HFXSQ, sdr, ser, ta, twm, xi, xl, xlp
  INTEGER :: j, jj, m, mm, ns
  !***********************************************************************
  !***  VERSION III-UPPER LIMIT OBTAINED FROM THE EXPRESSIONS OF          
  !***             CORBATO AND URETSKY USING A LIMIT OF E SUB M OF        
  !***             2**-30 WHICH IS APPROXIMATELY 1.E-9                    
  !***            SUMMATION PROPERTY USED TO NORMALIZE RESULTS.           
  !***  ADDITIONAL FACTOR ADDED TO STARTING VALUE                         
  !***  N IS THE MAXIMUM L TO BE CALCULATED                               
  !***  X IS THE ARGUMENT                                                 
  !***  FJ IS THE ARRAY THAT THE SPHERICAL BESSEL FUNCTIONS ARE TO BE     
  !***  PLACED IN.                                                        
  !*****  MODIFIED TO NOT REQUIRE THE WORKING SPACE.                      
  !*****        29 MAY,1968                                               
  !***********************************************************************
  IF(N.LT.0) THEN
     WRITE(6,*) 'ERROR : N SHOULD NOT BE NEGATIVE'
     STOP 'sphbes2'
  ENDIF
  IF (X.LT.0.d0) THEN
     WRITE(6,*) 'ERROR : in sphbes2, X should not be negative'
     STOP 'sphbes2'
  ENDIF
  
  IF (X.LE.XLIM) THEN
     HFXSQ=0.5d0*X*X                                                      
     XL=1.d0                                                            
     TWM=1.d0                                                           
     M=0
     DO
        M=M+1
        TA=XL
        TWM=TWM+2.d0
        XL=XL/TWM
        TA=TA-XL*HFXSQ
        XLP=XL/(TWM+2.d0)
        FJ(M)=TA+0.5d0*XLP*HFXSQ*HFXSQ
        XL=XL*X
        IF (M.GT.N) EXIT
     ENDDO
     RETURN
  ENDIF
  !
  CUFAC=4.2D0                                                       
  IF (X.LT.(N-2)) CUFAC=TNHF/(N+0.5d0-X)                               
  NS = N+5+X*CUFAC
  !*******************  ADD ADDITIONAL FACTOR  ***************************
  NS=NS + (FFT/(1.d0+SQRT(X)))                                       
  FFO=0.d0                                                          
  FFN=TN25                                                          
  M=NS-1                                                            
  XI=1.d0/X                                                          
  FM=(M+M)+1.d0                                                      
  SDR=FM*TN50                                                       
314 CONTINUE
  FFP=FM*XI*FFN-FFO                                                 
  IF (ABS(FFP).LT.T25) GO TO 315                                    
  SDR=SDR*TN50                                                      
  FFP=FFP*TN25                                                      
  FFN=FFN*TN25                                                      
315 CONTINUE
  SDR=SDR + (FM-2.d0)*FFP*FFP                                        
  FFO=FFN                                                           
  FFN=FFP                                                           
  IF (M.LE.N) GO TO 316                                             
  M=M-1                                                             
  FM=FM-2.d0                                                         
  GO TO 314                                                         
316 CONTINUE
  FJ(M)=FFN                                                         
  FJ(M+1)=FFO                                                       
  GO TO 33                                                          
32 CONTINUE
  FJ(M)=FM*XI*FJ(M+1)-FJ(M+2)                                       
  IF(ABS(FJ(M)).GE.T25) GO TO 56                                    
  SDR=SDR + (FM-2.d0)*FJ(M)*FJ(M)                                    
  IF (M.LE.1) GO TO 34                                              
33 CONTINUE
  M = M-1                                                           
  FM=FM-2.d0                                                         
  GO TO 32                                                          
34 CONTINUE
  SER=1.d0/SQRT(SDR)                                                 
  MM = N+1
  DO M=1,MM                                                      
     FJ(M)=FJ(M)*SER                                                   
  ENDDO
  RETURN
56 CONTINUE
  JJ= M+1                                                           
  NS=N+1                                                            
  DO J = JJ,NS                                                   
     FJ(J)=FJ(J)*TN25                                                  
  ENDDO
  SDR=SDR*TN50                                                      
  GO TO 32                                                          
END SUBROUTINE SPHBES2

