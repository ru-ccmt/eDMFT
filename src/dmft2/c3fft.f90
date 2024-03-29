SUBROUTINE C3FFT(N1,N2,N3,C,LDC1,LDC2,ISIG,CWORK,csize,DWORK,dsize,IERR)
!
!     ..................................................................
! 1.     PROGRAM UNIT 'C3FFT'
!           3-dimensional Fast Fourier Transform of a COMPLEX*16 array
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           This routine performs either a forward or a backward FFT of
!           three-dimensional double complex (COMPLEX*16) data.
!           The primefactors of the number of datapoints along each
!           dimension of the data array should be small to yield better
!           performance.
!
! 3.     USAGE
!        PARAMETER-DESCRIPTION
!           C      - COMPLEX*16 array                     (input/output)
!                    Input: The (three-dimensional) array containing the
!                           data to be transformed.
!                    Output: the transformed (three-dimensional) data
!           LDC1   - INTEGER variable                            (input)
!                    The exact dimension of the first index of C as
!                    defined in the calling routine.
!           LDC2   - INTEGER variable                            (input)
!                    The exact dimension of the second index of C as
!                    defined in the calling routine.
!           N1     - INTEGER variable                            (input)
!                    The number of datapoints in the first dimension of
!                    array C.
!                    Constraint: N1 .LE. LDC1 (not checked here)
!           N2     - INTEGER variable                            (input)
!                    The number of datapoints in the second dimension of
!                    array C.
!                    Constraint: N2 .LE. LDC2 (not checked here)
!           N3     - INTEGER variable                            (input)
!                    The number of datapoints in the third dimension of
!                    array C.
!           IERR   - INTEGER variable                           (output)
!                    Error indicator.
!                    IERR .EQ. 0 ... no error
!                    IERR .NE. 0 ... currently not supported
!           ISIG   - INTEGER variable                            (input)
!                    Specifies the direction of the transform:
!                       ISIG .LT. 0 ... forward FFT
!                       ISIG .GT. 0 ... backward FFT
!                       ISIG = 0 ...... nothing will be done
!           CWORK  - COMPLEX*16 array                     (input/output)
!                    A work array which must be dimensioned at least
!                    max(N1,N2,N3) in the calling routine.
!           DWORK  - DOUBLE PRECISION array               (input/output)
!                    A work array which must be dimensioned at least
!                    4*(max(N1,N2,N3))+15 in the calling routine.
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           CFFTB  - backward Fast Fourier Transform (from FFTPACK)
!           CFFTF  - forward Fast Fourier Transform (from FFTPACK)
!           CFFTI  - initialization for CFFTB and CFFTF (from FFTPACK)
!
!        INDIRECTLY CALLED SUBROUTINES
!           PASSB  - (from FFTPACK)
!           PASSB2 - (from FFTPACK)
!           PASSB3 - (from FFTPACK)
!           PASSB4 - (from FFTPACK)
!           PASSB5 - (from FFTPACK)
!           PASSF  - (from FFTPACK)
!           PASSF2 - (from FFTPACK)
!           PASSF3 - (from FFTPACK)
!           PASSF4 - (from FFTPACK)
!           PASSF5 - (from FFTPACK)
!           PIMACH - (from FFTPACK)
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           The array containing the three-dimensional data to be
!           transformed is declared as COMPLEX*16.
!
! 4.     REMARKS
!           For repeated usage of this routine no separate set-up
!           routine is provided, because the time needed for this
!           set-up is neglectable compared to the actual transformation.
!
! 5.     METHOD
!           The three-dimensional FFT is done using one-dimensional
!           complex FFT routines (from FFTPACK). This is done by:
!              1. transforming data values along the 1. dimension of C
!              2. transforming data values along the 2. dimension of C
!              3. transforming data values along the 3. dimension of C
!           Because the data elements along the 2. and 3. dimension of
!           C are not stored in contiguous storage locations and
!           FFTPACK-routines don't provide means for specifying a
!           stride, the involved elements are copied into a temporary
!           one-dimensional field before transformation. Afterwards
!           the data values are stored back to their original
!           location in C.
!
!
  IMPLICIT NONE
  INTEGER, intent(in)  :: LDC1,  LDC2,  N1,  N2,  N3,  ISIG, csize, dsize
  INTEGER, intent(out) :: IERR
  DOUBLE PRECISION, intent(inout) :: DWORK(dsize)
  COMPLEX*16, intent(inout)       :: CWORK(csize)
  COMPLEX*16, intent(inout) :: C(LDC1,LDC2,N3)
  !
  INTEGER            I,  J,  K
  !
  IERR = 0
  IF (ISIG .LT. 0) THEN
     !
     !        forward transform in 1. dimension
     !
     CALL CFFTI(N1,DWORK)
     DO K = 1, N3      ! 20
        DO J = 1, N2   ! 10
           CALL CFFTF(N1,C(1,J,K),DWORK)
        ENDDO  !10         CONTINUE
     ENDDO     !20         CONTINUE
     !
     !        forward transform in 2. dimension
     !
     CALL CFFTI(N2,DWORK)
     DO K = 1, N3         ! 60
        DO I = 1, N1      ! 50
           DO J = 1, N2   ! 30
              CWORK(J) = C(I,J,K)
           ENDDO          ! 30          CONTINUE
           CALL CFFTF(N2,CWORK,DWORK)
           DO J = 1, N2   ! 40
              C(I,J,K) = CWORK(J)
           ENDDO          ! 40          CONTINUE
        ENDDO             ! 50          CONTINUE
     ENDDO                ! 60          CONTINUE
     !
     !        forward transform in 3. dimension
     !
     CALL CFFTI(N3,DWORK)
     DO J = 1, N2          ! 100
        DO I = 1, N1       ! 90
           DO K = 1, N3    ! 70
              CWORK(K) = C(I,J,K)
           ENDDO           ! 70          CONTINUE
           CALL CFFTF(N3,CWORK,DWORK)
           DO K = 1, N3    ! 80
              C(I,J,K) = CWORK(K)
           ENDDO           ! 80          CONTINUE
        ENDDO              ! 90          CONTINUE
     ENDDO                 ! 100         CONTINUE
  ELSEIF (ISIG .GT. 0) THEN
     !
     !        backward transform in 1. dimension
     !
     CALL CFFTI(N1,DWORK)
     DO K = 1, N3      ! 120
        DO J = 1, N2   ! 110
           CALL CFFTB(N1,C(1,J,K),DWORK)
        ENDDO          ! 110       CONTINUE
     ENDDO             ! 120       CONTINUE
     !
     !        backward transform in 2. dimension
     !
     CALL CFFTI(N2,DWORK)
     DO K = 1, N3         ! 160
        DO I = 1, N1      ! 150
           DO J = 1, N2   ! 130
              CWORK(J) = C(I,J,K)
           ENDDO          ! 130          CONTINUE
           CALL CFFTB(N2,CWORK,DWORK)
           DO J = 1, N2   ! 140
              C(I,J,K) = CWORK(J)
           ENDDO          ! 140          CONTINUE
        ENDDO             ! 150          CONTINUE
     ENDDO                ! 160          CONTINUE
     !
     !        backward transform in 3. dimension
     !
     CALL CFFTI(N3,DWORK)
     DO J = 1, N2         ! 200
        DO I = 1, N1      ! 190
           DO K = 1, N3   ! 170
              CWORK(K) = C(I,J,K)
           ENDDO          ! 170          CONTINUE
           CALL CFFTB(N3,CWORK,DWORK)
           DO K = 1, N3   ! 180
              C(I,J,K) = CWORK(K)
           ENDDO          ! 180          CONTINUE
        ENDDO             ! 190          CONTINUE
     ENDDO                ! 200          CONTINUE
  ENDIF
  !
  !        End of 'C3FFT'
  !
END SUBROUTINE C3FFT


!     SUBROUTINE CFFTB(N,C,WSAVE)                                               
!                                                                               
!     SUBROUTINE CFFTB COMPUTES THE BACKWARD COMPLEX DISCRETE FOURIER           
!     TRANSFORM (THE FOURIER SYNTHESIS). EQUIVALENTLY , CFFTB COMPUTES          
!     A COMPLEX PERIODIC SEQUENCE FROM ITS FOURIER COEFFICIENTS.                
!     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.                     
!                                                                               
!     A CALL OF CFFTF FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE             
!     SEQUENCE BY N.                                                            
!                                                                               
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE CFFTB MUST BE                 
!     INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).                         
!                                                                               
!     INPUT PARAMETERS                                                          
!                                                                               
!                                                                               
!     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS                
!            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.              
!                                                                               
!     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE            
!                                                                               
!     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15        
!             IN THE PROGRAM THAT CALLS CFFTB. THE WSAVE ARRAY MUST BE          
!             INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE) AND A            
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT             
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE               
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT           
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.                 
!             THE SAME WSAVE ARRAY CAN BE USED BY CFFTF AND CFFTB.              
!                                                                               
!     OUTPUT PARAMETERS                                                         
!                                                                               
!     C      FOR J=1,...,N                                                      
!                                                                               
!                C(J)=THE SUM FROM K=1,...,N OF                                 
!                                                                               
!                      C(K)*EXP(I*(J-1)*(K-1)*2*PI/N)                           
!                                                                               
!                            WHERE I=SQRT(-1)                                   
!                                                                               
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE            
!             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB              
!
SUBROUTINE CFFTB (N,C,WSAVE)                                              
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION       C(*)       ,WSAVE(*)                                      
  !                                                                               
  IF (N .EQ. 1) RETURN                                                      
  IW1 = N+N+1                                                               
  IW2 = IW1+N+N                                                             
  CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))                             
  RETURN                                                                    
END SUBROUTINE CFFTB

SUBROUTINE CFFTB1 (N,C,CH,WA,IFAC)                                        
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)               
  NF = IFAC(2)                                                              
  NA = 0                                                                    
  L1 = 1                                                                    
  IW = 1                                                                    
  DO K1=1,NF     ! 116                                                       
     IP = IFAC(K1+2)                                                        
     L2 = IP*L1                                                             
     IDO = N/L2                                                             
     IDOT = IDO+IDO                                                         
     IDL1 = IDOT*L1                                                         
     IF (IP .NE. 4) GO TO 103                                               
     IX2 = IW+IDOT                                                          
     IX3 = IX2+IDOT                                                         
     IF (NA .NE. 0) GO TO 101
     CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))                      
     GO TO 102                                                              
101  CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))                      
102  NA = 1-NA                                                              
     GO TO 115                                                              
103  IF (IP .NE. 2) GO TO 106                                               
     IF (NA .NE. 0) GO TO 104                                               
     CALL PASSB2 (IDOT,L1,C,CH,WA(IW))                                      
     GO TO 105                                                              
104  CALL PASSB2 (IDOT,L1,CH,C,WA(IW))                                      
105  NA = 1-NA                                                              
     GO TO 115                                                              
106  IF (IP .NE. 3) GO TO 109                                               
     IX2 = IW+IDOT                                                          
     IF (NA .NE. 0) GO TO 107                                               
     CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))                              
     GO TO 108                                                              
107  CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))                              
108  NA = 1-NA                                                              
     GO TO 115                                                              
109  IF (IP .NE. 5) GO TO 112                                               
     IX2 = IW+IDOT                                                          
     IX3 = IX2+IDOT                                                         
     IX4 = IX3+IDOT                                                         
     IF (NA .NE. 0) GO TO 110                                               
     CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))              
     GO TO 111                                                              
110  CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))              
111  NA = 1-NA                                                              
     GO TO 115                                                              
112  IF (NA .NE. 0) GO TO 113                                               
     CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))                    
     GO TO 114                                                              
113  CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))                   
114  IF (NAC .NE. 0) NA = 1-NA                                              
115  L1 = L2                                                                
     IW = IW+(IP-1)*IDOT                                                    
  ENDDO   ! 116  CONTINUE                                                                  
  IF (NA .EQ. 0) RETURN                                                     
  N2 = N+N                                                                  
  DO I=1,N2       ! 117                                                       
     C(I) = CH(I)                                                           
  ENDDO           ! 117     CONTINUE                                                                  
  RETURN                                                                    
END SUBROUTINE CFFTB1

!     SUBROUTINE CFFTF(N,C,WSAVE)                                               
!                                                                               
!     SUBROUTINE CFFTF COMPUTES THE FORWARD COMPLEX DISCRETE FOURIER            
!     TRANSFORM (THE FOURIER ANALYSIS). EQUIVALENTLY , CFFTF COMPUTES           
!     THE FOURIER COEFFICIENTS OF A COMPLEX PERIODIC SEQUENCE.                  
!     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.                     
!                                                                               
!     THE TRANSFORM IS NOT NORMALIZED. TO OBTAIN A NORMALIZED TRANSFORM         
!     THE OUTPUT MUST BE DIVIDED BY N. OTHERWISE A CALL OF CFFTF                
!     FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE SEQUENCE BY N.              
!                                                                               
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE CFFTF MUST BE                 
!     INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).                         
!                                                                               
!     INPUT PARAMETERS                                                          
!                                                                               
!                                                                               
!     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS                
!            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES. N            
!                                                                               
!     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE            
!                                                                               
!     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15        
!             IN THE PROGRAM THAT CALLS CFFTF. THE WSAVE ARRAY MUST BE          
!             INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE) AND A            
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT             
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE               
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT           
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.                 
!             THE SAME WSAVE ARRAY CAN BE USED BY CFFTF AND CFFTB.              
!                                                                               
!     OUTPUT PARAMETERS                                                         
!                                                                               
!     C      FOR J=1,...,N                                                      
!                                                                               
!                C(J)=THE SUM FROM K=1,...,N OF                                 
!                                                                               
!                      C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N)                          
!                                                                               
!                            WHERE I=SQRT(-1)                                   
!                                                                               
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE            
!             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB              
!                                                                               
SUBROUTINE CFFTF (N,C,WSAVE)                                              
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION       C(*)       ,WSAVE(*)                                      
  !                                                                               
  IF (N .EQ. 1) RETURN                                                      
  IW1 = N+N+1                                                               
  IW2 = IW1+N+N                                                             
  CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))                             
  RETURN                                                                    
END SUBROUTINE CFFTF

SUBROUTINE CFFTF1 (N,C,CH,WA,IFAC)                                        
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION    CH(*),C(*),WA(*),IFAC(*)               
  NF = IFAC(2)                                                              
  NA = 0                                                                    
  L1 = 1                                                                    
  IW = 1                                                                    
  DO K1=1,NF           ! 116
     IP = IFAC(K1+2)                                                        
     L2 = IP*L1                                                             
     IDO = N/L2                                                             
     IDOT = IDO+IDO                                                         
     IDL1 = IDOT*L1                                                         
     IF (IP .NE. 4) GO TO 103                                               
     IX2 = IW+IDOT                                                          
     IX3 = IX2+IDOT                                                         
     IF (NA .NE. 0) GO TO 101                                               
     CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))                      
     GO TO 102                                                              
101  CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))                      
102  NA = 1-NA                                                              
     GO TO 115                                                              
103  IF (IP .NE. 2) GO TO 106                                               
     IF (NA .NE. 0) GO TO 104                                               
     CALL PASSF2 (IDOT,L1,C,CH,WA(IW))                                      
     GO TO 105                                                              
104  CALL PASSF2 (IDOT,L1,CH,C,WA(IW))                                      
105  NA = 1-NA                                                              
     GO TO 115                                                              
106  IF (IP .NE. 3) GO TO 109                                               
     IX2 = IW+IDOT                                                          
     IF (NA .NE. 0) GO TO 107                                               
     CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))                              
     GO TO 108                                                              
107  CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))                              
108  NA = 1-NA                                                              
     GO TO 115                                                              
109  IF (IP .NE. 5) GO TO 112                                               
     IX2 = IW+IDOT                                                          
     IX3 = IX2+IDOT                                                         
     IX4 = IX3+IDOT                                                         
     IF (NA .NE. 0) GO TO 110                                               
     CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))              
     GO TO 111                                                              
110  CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))              
111  NA = 1-NA                                                              
     GO TO 115                                                              
112  IF (NA .NE. 0) GO TO 113                                               
     CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))                    
     GO TO 114                                                              
113  CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))                   
114  IF (NAC .NE. 0) NA = 1-NA                                              
115  L1 = L2                                                                
     IW = IW+(IP-1)*IDOT                                                    
  ENDDO !116  CONTINUE                                                                  
  IF (NA .EQ. 0) RETURN                                                     
  N2 = N+N                                                                  
  DO I=1,N2               ! 117                                              
     C(I) = CH(I)                                                           
  ENDDO                   ! 117     CONTINUE                                                                  
  RETURN                                                                    
END SUBROUTINE CFFTF1

!     SUBROUTINE CFFTI(N,WSAVE)                                                 
!                                                                               
!     SUBROUTINE CFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN             
!     BOTH CFFTF AND CFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH          
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND              
!     STORED IN WSAVE.                                                          
!                                                                               
!     INPUT PARAMETER                                                           
!                                                                               
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED                      
!                                                                               
!     OUTPUT PARAMETER                                                          
!                                                                               
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4*N+15            
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH CFFTF AND CFFTB          
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS            
!             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF           
!             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF CFFTF OR CFFTB.        
!                                                                               
SUBROUTINE CFFTI (N,WSAVE)                                                
  IMPLICIT REAL*8  (A-H,O-Z)
  DIMENSION        WSAVE(*)                                                  
  !                                                                               
  IF (N .EQ. 1) RETURN                                                      
  IW1 = N+N+1                                                               
  IW2 = IW1+N+N                                                             
  CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))                                     
  RETURN                                                                    
END SUBROUTINE CFFTI

SUBROUTINE CFFTI1 (N,WA,IFAC)
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
  DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
  NL = N
  NF = 0
  J = 0
101 J = J+1
  IF (J-4) 102,102,103
102 NTRY = NTRYH(J)
  GO TO 104
103 NTRY = NTRY+2
104 NQ = NL/NTRY
  NR = NL-NTRY*NQ
  IF (NR) 101,105,101
105 NF = NF+1
  IFAC(NF+2) = NTRY
  NL = NQ
  IF (NTRY .NE. 2) GO TO 107
  IF (NF .EQ. 1) GO TO 107
  DO I=2,NF     ! 106
     IB = NF-I+2
     IFAC(IB+2) = IFAC(IB+1)
  ENDDO         ! 106  CONTINUE
  IFAC(3) = 2
107 IF (NL .NE. 1) GO TO 104
  IFAC(1) = N
  IFAC(2) = NF
  TPI = 2.0D+0*PIMACH(DUM)
  ARGH = TPI/DFLOAT(N)
  I = 2
  L1 = 1
  DO K1=1,NF            ! 110
     IP = IFAC(K1+2)
     LD = 0
     L2 = L1*IP
     IDO = N/L2
     IDOT = IDO+IDO+2
     IPM = IP-1
     DO J=1,IPM         ! 109
        I1 = I
        WA(I-1) = 1.0D+0
        WA(I) = 0.0D+0
        LD = LD+L1
        FI = 0.0D+0
        ARGLD = DFLOAT(LD)*ARGH
        DO II=4,IDOT,2  ! 108
           I = I+2
           FI = FI+1.0D+0
           ARG = FI*ARGLD
           WA(I-1) = COS(ARG)
           WA(I) = SIN(ARG)
        ENDDO           ! 108        CONTINUE
        IF (IP .LE. 5) CYCLE !GO TO 109
        WA(I1-1) = WA(I-1)
        WA(I1) = WA(I)
     ENDDO              ! 109        CONTINUE
        
     L1 = L2
  ENDDO                 ! 110 CONTINUE
  RETURN
END SUBROUTINE CFFTI1

SUBROUTINE PASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)                  
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION       CH(IDO,L1,IP),CC(IDO,IP,L1),C1(IDO,L1,IP),WA(*),C2(IDL1,IP),CH2(IDL1,IP)
  IDOT = IDO/2                                                              
  NT = IP*IDL1                                                              
  IPP2 = IP+2                                                               
  IPPH = (IP+1)/2                                                           
  IDP = IP*IDO                                                              
  !                                                                               
  IF (IDO .LT. L1) GO TO 106                                                
  DO J=2,IPPH                  ! 103
     JC = IPP2-J                                                            
     DO K=1,L1                 ! 102
        DO I=1,IDO             ! 101
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)                                 
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)                                
        ENDDO                  ! 101 CONTINUE
     ENDDO                     ! 102 CONTINUE
  ENDDO                        ! 103 CONTINUE
  DO K=1,L1                    ! 105
     DO I=1,IDO                ! 104
        CH(I,K,1) = CC(I,1,K)                                               
     ENDDO                     ! 104 CONTINUE                                                               
  ENDDO                        ! 105 CONTINUE                                                                  
  GO TO 112                                                                 
106 CONTINUE
  DO J=2,IPPH                  ! 109
     JC = IPP2-J                                                            
     DO I=1,IDO                ! 108
        DO K=1,L1              ! 107
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)                                 
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)                                
        ENDDO                  ! 107        CONTINUE                                                            
     ENDDO                     ! 108        CONTINUE                                                               
  ENDDO                        ! 109        CONTINUE
  
  
  DO I=1,IDO                  ! 111
     DO K=1,L1                ! 110
        CH(I,K,1) = CC(I,1,K)                                               
     ENDDO                    ! 110 CONTINUE                                                               
  ENDDO                       ! 111 CONTINUE                                                                  
112 IDL = 2-IDO                                                               
  INC = 0                                                                   
  DO L=2,IPPH                 ! 116
     LC = IPP2-L                                                            
     IDL = IDL+IDO                                                          
     DO IK=1,IDL1             ! 113
        C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)                            
        C2(IK,LC) = WA(IDL)*CH2(IK,IP)                                      
     ENDDO                    ! 113    CONTINUE                                                               
     IDLJ = IDL                                                             
     INC = INC+IDO                                                          
     DO J=3,IPPH              ! 115
        JC = IPP2-J                                                         
        IDLJ = IDLJ+INC                                                     
        IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP                                  
        WAR = WA(IDLJ-1)                                                    
        WAI = WA(IDLJ)                                                      
        DO IK=1,IDL1          ! 114
           C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)                                
           C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)                             
        ENDDO                 ! 114 CONTINUE                                                            
     ENDDO                    ! 115 CONTINUE                                                               
  ENDDO                       ! 116 CONTINUE                                                                  
  DO J=2,IPPH                 ! 118
     DO IK=1,IDL1             ! 117
        CH2(IK,1) = CH2(IK,1)+CH2(IK,J)                                     
     ENDDO                    ! 117 CONTINUE                                                               
  ENDDO                       ! 118 CONTINUE                                                                  
  DO J=2,IPPH                 ! 120
     JC = IPP2-J                                                            
     DO IK=2,IDL1,2           ! 119
        CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)                                  
        CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)                                 
        CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)                                    
        CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)                                   
     ENDDO                    ! 119 CONTINUE
  ENDDO                       ! 120 CONTINUE
  NAC = 1                                                                   
  IF (IDO .EQ. 2) RETURN                                                    
  NAC = 0                                                                   
  DO IK=1,IDL1                ! 121
     C2(IK,1) = CH2(IK,1)                                                   
  ENDDO                       ! 121  CONTINUE                                                                  
  DO J=2,IP                   ! 123
     DO K=1,L1                ! 122
        C1(1,K,J) = CH(1,K,J)                                               
        C1(2,K,J) = CH(2,K,J)                                               
     ENDDO                    ! 122 CONTINUE                                                               
  ENDDO                       ! 123 CONTINUE                                                                  
  IF (IDOT .GT. L1) GO TO 127                                               
  IDIJ = 0                                                                  
  DO J=2,IP                 ! 126
     IDIJ = IDIJ+2                                                          
     DO I=4,IDO,2           ! 125
        IDIJ = IDIJ+2                                                       
        DO K=1,L1           ! 124                                                     
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)          
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)            
        ENDDO               ! 124 CONTINUE                                                            
     ENDDO                  ! 125 CONTINUE                                                               
  ENDDO                     ! 126 CONTINUE                                                                  
  RETURN                                                                    
127 IDJ = 2-IDO                                                               
  DO J=2,IP                ! 130
     IDJ = IDJ+IDO                                                          
     DO K=1,L1             ! 129
        IDIJ = IDJ                                                          
        DO I=4,IDO,2       ! 128
           IDIJ = IDIJ+2                                                    
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)          
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)            
        ENDDO              ! 128 CONTINUE                                                            
     ENDDO                 ! 129 CONTINUE                                                               
  ENDDO                    ! 130 CONTINUE                                                                  
  RETURN                                                                    
END SUBROUTINE PASSB



   
SUBROUTINE PASSB2 (IDO,L1,CC,CH,WA1)                                      
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION  CC(IDO,2,L1),CH(IDO,L1,2),WA1(*)                                                    
  IF (IDO .GT. 2) GO TO 102
  DO K=1,L1         ! 101
     CH(1,K,1) = CC(1,1,K)+CC(1,2,K)                                        
     CH(1,K,2) = CC(1,1,K)-CC(1,2,K)                                        
     CH(2,K,1) = CC(2,1,K)+CC(2,2,K)                                        
     CH(2,K,2) = CC(2,1,K)-CC(2,2,K)                                        
  ENDDO             ! 101  CONTINUE                                                                  
  RETURN                                                                    
102 CONTINUE
  DO K=1,L1         ! 104
     DO I=2,IDO,2   ! 103
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)                               
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)                                       
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)                                     
        TI2 = CC(I,1,K)-CC(I,2,K)                                           
        CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2                                 
        CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2                               
     ENDDO          ! 103    CONTINUE                                                               
  ENDDO             ! 104 CONTINUE                                                                  
  RETURN                                                                    
END SUBROUTINE PASSB2


SUBROUTINE PASSB3 (IDO,L1,CC,CH,WA1,WA2)
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION   CC(IDO,3,L1),CH(IDO,L1,3),WA1(*),WA2(*)
  DATA TAUR,TAUI /-0.5D+0,0.866025403784439D+0/
  IF (IDO .NE. 2) GO TO 102
  DO K=1,L1  ! 101
     TR2 = CC(1,2,K)+CC(1,3,K)
     CR2 = CC(1,1,K)+TAUR*TR2
     CH(1,K,1) = CC(1,1,K)+TR2
     TI2 = CC(2,2,K)+CC(2,3,K)
     CI2 = CC(2,1,K)+TAUR*TI2
     CH(2,K,1) = CC(2,1,K)+TI2
     CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
     CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
     CH(1,K,2) = CR2-CI3
     CH(1,K,3) = CR2+CI3
     CH(2,K,2) = CI2+CR3
     CH(2,K,3) = CI2-CR3
  ENDDO ! 101  CONTINUE
  RETURN
102 CONTINUE
  DO K=1,L1                ! 104
     DO I=2,IDO,2          ! 103
        TR2 = CC(I-1,2,K)+CC(I-1,3,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,2,K)+CC(I,3,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
        CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
     ENDDO                ! 103    CONTINUE
  ENDDO ! 104     CONTINUE
  RETURN
END SUBROUTINE PASSB3

SUBROUTINE PASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)                              
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION CC(IDO,4,L1),CH(IDO,L1,4),WA1(*),WA2(*),WA3(*)
  IF (IDO .NE. 2) GO TO 102                                                 
  DO K=1,L1     ! 101
     TI1 = CC(2,1,K)-CC(2,3,K)                                              
     TI2 = CC(2,1,K)+CC(2,3,K)                                              
     TR4 = CC(2,4,K)-CC(2,2,K)                                              
     TI3 = CC(2,2,K)+CC(2,4,K)                                              
     TR1 = CC(1,1,K)-CC(1,3,K)                                              
     TR2 = CC(1,1,K)+CC(1,3,K)                                              
     TI4 = CC(1,2,K)-CC(1,4,K)                                              
     TR3 = CC(1,2,K)+CC(1,4,K)                                              
     CH(1,K,1) = TR2+TR3                                                    
     CH(1,K,3) = TR2-TR3                                                    
     CH(2,K,1) = TI2+TI3                                                    
     CH(2,K,3) = TI2-TI3                                                    
     CH(1,K,2) = TR1+TR4                                                    
     CH(1,K,4) = TR1-TR4                                                    
     CH(2,K,2) = TI1+TI4                                                    
     CH(2,K,4) = TI1-TI4                                                    
  ENDDO ! 101 CONTINUE                                                                  
  RETURN                                                                    
102 CONTINUE
  DO K=1,L1            ! 104
     DO I=2,IDO,2      ! 103
        TI1 = CC(I,1,K)-CC(I,3,K)                                           
        TI2 = CC(I,1,K)+CC(I,3,K)                                           
        TI3 = CC(I,2,K)+CC(I,4,K)                                           
        TR4 = CC(I,4,K)-CC(I,2,K)                                           
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)                                       
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)                                       
        TI4 = CC(I-1,2,K)-CC(I-1,4,K)                                       
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)                                       
        CH(I-1,K,1) = TR2+TR3                                               
        CR3 = TR2-TR3                                                       
        CH(I,K,1) = TI2+TI3                                                 
        CI3 = TI2-TI3                                                       
        CR2 = TR1+TR4                                                       
        CR4 = TR1-TR4                                                       
        CI2 = TI1+TI4                                                       
        CI4 = TI1-TI4                                                       
        CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2                               
        CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2                                 
        CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3                               
        CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3                                 
        CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4                               
        CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4                                 
     ENDDO          ! 103 CONTINUE                                                               
  ENDDO             ! 104 CONTINUE                                                                  
  RETURN                                                                    
END SUBROUTINE PASSB4

SUBROUTINE PASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)                          
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION  CC(IDO,5,L1),CH(IDO,L1,5),WA1(*),WA2(*),WA3(*),WA4(*)
  DATA TR11,TI11,TR12,TI12 /0.309016994374947D+0,0.951056516295154D+0,-0.809016994374947D+0,0.587785252292473D+0/
  IF (IDO .NE. 2) GO TO 102                                                 
  DO K=1,L1       ! 101
     TI5 = CC(2,2,K)-CC(2,5,K)                                              
     TI2 = CC(2,2,K)+CC(2,5,K)                                              
     TI4 = CC(2,3,K)-CC(2,4,K)                                              
     TI3 = CC(2,3,K)+CC(2,4,K)                                              
     TR5 = CC(1,2,K)-CC(1,5,K)                                              
     TR2 = CC(1,2,K)+CC(1,5,K)                                              
     TR4 = CC(1,3,K)-CC(1,4,K)                                              
     TR3 = CC(1,3,K)+CC(1,4,K)                                              
     CH(1,K,1) = CC(1,1,K)+TR2+TR3                                          
     CH(2,K,1) = CC(2,1,K)+TI2+TI3                                          
     CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3                                      
     CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3                                      
     CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3                                      
     CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3                                      
     CR5 = TI11*TR5+TI12*TR4                                                
     CI5 = TI11*TI5+TI12*TI4                                                
     CR4 = TI12*TR5-TI11*TR4                                                
     CI4 = TI12*TI5-TI11*TI4                                                
     CH(1,K,2) = CR2-CI5                                                    
     CH(1,K,5) = CR2+CI5                                                    
     CH(2,K,2) = CI2+CR5                                                    
     CH(2,K,3) = CI3+CR4                                                    
     CH(1,K,3) = CR3-CI4                                                    
     CH(1,K,4) = CR3+CI4                                                    
     CH(2,K,4) = CI3-CR4                                                    
     CH(2,K,5) = CI2-CR5                                                    
  ENDDO ! 101 CONTINUE                                                                  
  RETURN                                                                    
102 CONTINUE
  DO K=1,L1               ! 104
     DO I=2,IDO,2         ! 103
        TI5 = CC(I,2,K)-CC(I,5,K)                                           
        TI2 = CC(I,2,K)+CC(I,5,K)                                           
        TI4 = CC(I,3,K)-CC(I,4,K)                                           
        TI3 = CC(I,3,K)+CC(I,4,K)                                           
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)                                       
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)                                       
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)                                       
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)                                       
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3                                   
        CH(I,K,1) = CC(I,1,K)+TI2+TI3                                       
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3                                 
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3                                   
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3                                 
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3                                   
        CR5 = TI11*TR5+TI12*TR4                                             
        CI5 = TI11*TI5+TI12*TI4                                             
        CR4 = TI12*TR5-TI11*TR4                                             
        CI4 = TI12*TI5-TI11*TI4                                             
        DR3 = CR3-CI4                                                       
        DR4 = CR3+CI4                                                       
        DI3 = CI3+CR4                                                       
        DI4 = CI3-CR4                                                       
        DR5 = CR2+CI5                                                       
        DR2 = CR2-CI5                                                       
        DI5 = CI2-CR5                                                       
        DI2 = CI2+CR5                                                       
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2                               
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2                                 
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3                               
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3                                 
        CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4                               
        CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4                                 
        CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5                               
        CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5                                 
     ENDDO               ! 103 CONTINUE                                                               
  ENDDO                  ! 104 CONTINUE                                                                  
  RETURN                                                                    
END SUBROUTINE PASSB5

SUBROUTINE PASSF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)                  
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION   CH(IDO,L1,IP),CC(IDO,IP,L1),C1(IDO,L1,IP),WA(*),C2(IDL1,IP),CH2(IDL1,IP)                                              
  IDOT = IDO/2                                                              
  NT = IP*IDL1                                                              
  IPP2 = IP+2                                                               
  IPPH = (IP+1)/2                                                           
  IDP = IP*IDO                                                              
  !                                                                               
  IF (IDO .LT. L1) GO TO 106                                                
  DO J=2,IPPH       ! 103
     JC = IPP2-J                                                            
     DO K=1,L1      ! 102
        DO I=1,IDO  ! 101                                                    
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)                                 
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)                                
        ENDDO       ! 101    CONTINUE                                                            
     ENDDO          ! 102    CONTINUE                                                               
  ENDDO             ! 103    CONTINUE
  
  DO K=1,L1         ! 105
     DO I=1,IDO     ! 104
        CH(I,K,1) = CC(I,1,K)                                               
     ENDDO          ! 104    CONTINUE                                                               
  ENDDO             ! 105    CONTINUE
  
  GO TO 112                                                                 
106 CONTINUE
  DO J=2,IPPH        ! 109
     JC = IPP2-J                                                            
     DO I=1,IDO      ! 108
        DO K=1,L1    ! 107
           CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)                                 
           CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)                                
        ENDDO       ! 107        CONTINUE                                                            
     ENDDO          ! 108        CONTINUE                                                               
  ENDDO             ! 109        CONTINUE
  
  
  DO I=1,IDO        ! 111
     DO K=1,L1      ! 110
        CH(I,K,1) = CC(I,1,K)                                               
     ENDDO          ! 110     CONTINUE                                                               
  ENDDO             ! 111     CONTINUE                                                                  
112 CONTINUE
  IDL = 2-IDO                                                               
  INC = 0                                                                   
  DO L=2,IPPH       ! 116
     LC = IPP2-L                                                            
     IDL = IDL+IDO                                                          
     DO IK=1,IDL1   ! 113
        C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)                            
        C2(IK,LC) = -WA(IDL)*CH2(IK,IP)                                     
     ENDDO          ! 113    CONTINUE                                                               
     IDLJ = IDL                                                             
     INC = INC+IDO                                                          
     DO J=3,IPPH   ! 115
        JC = IPP2-J                                                         
        IDLJ = IDLJ+INC                                                     
        IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP                                  
        WAR = WA(IDLJ-1)                                                    
        WAI = WA(IDLJ)                                                      
        DO IK=1,IDL1   ! 114
           C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)                                
           C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)                             
        ENDDO          ! 114  CONTINUE                                                            
     ENDDO             ! 115    CONTINUE                                                               
  ENDDO                ! 116 CONTINUE                                                                  
  DO J=2,IPPH          ! 118
     DO IK=1,IDL1      ! 117
        CH2(IK,1) = CH2(IK,1)+CH2(IK,J)                                     
     ENDDO             ! 117    CONTINUE                                                               
  ENDDO                ! 118 CONTINUE                                                                  
  DO J=2,IPPH          ! 120
     JC = IPP2-J                                                            
     DO IK=2,IDL1,2    ! 119
        CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)                                  
        CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)                                 
        CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)                                    
        CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)                                   
     ENDDO             ! 119     CONTINUE                                                               
  ENDDO                ! 120 CONTINUE                                                                  
  NAC = 1                                                                   
  IF (IDO .EQ. 2) RETURN                                                    
  NAC = 0                                                                   
  DO IK=1,IDL1         ! 121
     C2(IK,1) = CH2(IK,1)                                                   
  ENDDO                ! 121 CONTINUE                                                                  
  DO J=2,IP            ! 123
     DO K=1,L1         ! 122
        C1(1,K,J) = CH(1,K,J)                                               
        C1(2,K,J) = CH(2,K,J)                                               
     ENDDO             ! 122 CONTINUE
  ENDDO                ! 123 CONTINUE
  IF (IDOT .GT. L1) GO TO 127                                               
  IDIJ = 0                                                                  
  DO J=2,IP           ! 126
     IDIJ = IDIJ+2
     DO I=4,IDO,2     ! 125
        IDIJ = IDIJ+2
        DO K=1,L1     ! 124
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)          
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)            
        ENDDO         ! 124 CONTINUE
     ENDDO            ! 125 CONTINUE
  ENDDO               ! 126 CONTINUE
  RETURN                                                                    
127 IDJ = 2-IDO                                                               
  DO J=2,IP           ! 130
     IDJ = IDJ+IDO                                                          
     DO K=1,L1        ! 129
        IDIJ = IDJ
        DO I=4,IDO,2  ! 128
           IDIJ = IDIJ+2
           C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
           C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
        ENDDO         ! 128 CONTINUE                                                            
     ENDDO            ! 129 CONTINUE
  ENDDO               ! 130 CONTINUE
  RETURN                                                                    
END SUBROUTINE PASSF

SUBROUTINE PASSF2 (IDO,L1,CC,CH,WA1)                                      
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION CC(IDO,2,L1),CH(IDO,L1,2),WA1(*)                                                    
  IF (IDO .GT. 2) GO TO 102                                                 
  DO K=1,L1 ! 101
     CH(1,K,1) = CC(1,1,K)+CC(1,2,K)                                        
     CH(1,K,2) = CC(1,1,K)-CC(1,2,K)                                        
     CH(2,K,1) = CC(2,1,K)+CC(2,2,K)                                        
     CH(2,K,2) = CC(2,1,K)-CC(2,2,K)                                        
  ENDDO     ! 101 CONTINUE                                                                  
  RETURN                                                                    
102 CONTINUE
  DO K=1,L1        ! 104
     DO I=2,IDO,2  ! 103                                                     
        CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)                               
        TR2 = CC(I-1,1,K)-CC(I-1,2,K)                                       
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)                                     
        TI2 = CC(I,1,K)-CC(I,2,K)                                           
        CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2                                 
        CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2                               
     ENDDO         !  103    CONTINUE                                                               
  ENDDO            ! 104 CONTINUE                                                                  
  RETURN                                                                    
END SUBROUTINE PASSF2

SUBROUTINE PASSF3 (IDO,L1,CC,CH,WA1,WA2)
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION  CC(IDO,3,L1),CH(IDO,L1,3),WA1(*),WA2(*)
  DATA TAUR,TAUI /-0.5D+0,-0.866025403784439D+0/
  IF (IDO .NE. 2) GO TO 102

  ! !$OMP PARALLEL DO SHARED(CC,CH,TAUR,TAUI) PRIVATE(TR2,TI2,CR2,CI2,CR3,CI3) SCHEDULE(STATIC)
  DO K=1,L1
     TR2 = CC(1,2,K)+CC(1,3,K)
     CR2 = CC(1,1,K)+TAUR*TR2
     CH(1,K,1) = CC(1,1,K)+TR2
     TI2 = CC(2,2,K)+CC(2,3,K)
     CI2 = CC(2,1,K)+TAUR*TI2
     CH(2,K,1) = CC(2,1,K)+TI2
     CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
     CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
     CH(1,K,2) = CR2-CI3
     CH(1,K,3) = CR2+CI3
     CH(2,K,2) = CI2+CR3
     CH(2,K,3) = CI2-CR3
  ENDDO
  ! !$OMP END PARALLEL DO
  
  RETURN
102 CONTINUE
  ! !$OMP PARALLEL DO SHARED(CC,CH,TAUR,TAUI,WA1,WA2) PRIVATE(I,TR2,TI2,CR2,CI2,CR3,CI3,DR2,DI2,DR3,DI3) SCHEDULE(STATIC)
  DO K=1,L1         
     DO I=2,IDO,2   
        TR2 = CC(I-1,2,K)+CC(I-1,3,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,2,K)+CC(I,3,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
        CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
        CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
        CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
        CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
     ENDDO          
  ENDDO
  ! !$OMP END PARALLEL DO
  RETURN
END SUBROUTINE PASSF3


SUBROUTINE PASSF4 (IDO,L1,CC,CH,WA1,WA2,WA3)                              
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION  CC(IDO,4,L1),CH(IDO,L1,4),WA1(*),WA2(*),WA3(*)
  IF (IDO .NE. 2) GO TO 102                                                 
  DO K=1,L1   ! 101                                                          
     TI1 = CC(2,1,K)-CC(2,3,K)                                              
     TI2 = CC(2,1,K)+CC(2,3,K)                                              
     TR4 = CC(2,2,K)-CC(2,4,K)                                              
     TI3 = CC(2,2,K)+CC(2,4,K)                                              
     TR1 = CC(1,1,K)-CC(1,3,K)                                              
     TR2 = CC(1,1,K)+CC(1,3,K)                                              
     TI4 = CC(1,4,K)-CC(1,2,K)                                              
     TR3 = CC(1,2,K)+CC(1,4,K)                                              
     CH(1,K,1) = TR2+TR3                                                    
     CH(1,K,3) = TR2-TR3                                                    
     CH(2,K,1) = TI2+TI3                                                    
     CH(2,K,3) = TI2-TI3                                                    
     CH(1,K,2) = TR1+TR4                                                    
     CH(1,K,4) = TR1-TR4                                                    
     CH(2,K,2) = TI1+TI4                                                    
     CH(2,K,4) = TI1-TI4                                                    
  ENDDO       ! 101 CONTINUE                                                                  
  RETURN                                                                    
102 CONTINUE
  DO K=1,L1           ! 104
     DO I=2,IDO,2     ! 103
        TI1 = CC(I,1,K)-CC(I,3,K)                                           
        TI2 = CC(I,1,K)+CC(I,3,K)                                           
        TI3 = CC(I,2,K)+CC(I,4,K)                                           
        TR4 = CC(I,2,K)-CC(I,4,K)                                           
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)                                       
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)                                       
        TI4 = CC(I-1,4,K)-CC(I-1,2,K)                                       
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)                                       
        CH(I-1,K,1) = TR2+TR3                                               
        CR3 = TR2-TR3                                                       
        CH(I,K,1) = TI2+TI3                                                 
        CI3 = TI2-TI3                                                       
        CR2 = TR1+TR4                                                       
        CR4 = TR1-TR4                                                       
        CI2 = TI1+TI4                                                       
        CI4 = TI1-TI4                                                       
        CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2                               
        CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2                                 
        CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3                               
        CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3                                 
        CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4                               
        CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4                                 
     ENDDO          ! 103     CONTINUE                                                               
  ENDDO             ! 104     CONTINUE                                                                  
  RETURN                                                                    
END SUBROUTINE PASSF4

SUBROUTINE PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
  IMPLICIT REAL*8     (A-H,O-Z)
  DIMENSION       CC(IDO,5,L1), CH(IDO,L1,5), WA1(*), WA2(*), WA3(*), WA4(*)
  DATA TR11,TI11,TR12,TI12 /0.309016994374947D+0,-0.951056516295154D+0,-0.809016994374947D+0,-0.587785252292473D+0/
  IF (IDO .NE. 2) GO TO 102
  DO K=1,L1   ! 101
     TI5 = CC(2,2,K)-CC(2,5,K)
     TI2 = CC(2,2,K)+CC(2,5,K)
     TI4 = CC(2,3,K)-CC(2,4,K)
     TI3 = CC(2,3,K)+CC(2,4,K)
     TR5 = CC(1,2,K)-CC(1,5,K)
     TR2 = CC(1,2,K)+CC(1,5,K)
     TR4 = CC(1,3,K)-CC(1,4,K)
     TR3 = CC(1,3,K)+CC(1,4,K)
     CH(1,K,1) = CC(1,1,K)+TR2+TR3
     CH(2,K,1) = CC(2,1,K)+TI2+TI3
     CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
     CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
     CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
     CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
     CR5 = TI11*TR5+TI12*TR4
     CI5 = TI11*TI5+TI12*TI4
     CR4 = TI12*TR5-TI11*TR4
     CI4 = TI12*TI5-TI11*TI4
     CH(1,K,2) = CR2-CI5
     CH(1,K,5) = CR2+CI5
     CH(2,K,2) = CI2+CR5
     CH(2,K,3) = CI3+CR4
     CH(1,K,3) = CR3-CI4
     CH(1,K,4) = CR3+CI4
     CH(2,K,4) = CI3-CR4
     CH(2,K,5) = CI2-CR5
  ENDDO ! 101  CONTINUE
  RETURN
102 CONTINUE
  DO K=1,L1        ! 104
     DO I=2,IDO,2  ! 103
        TI5 = CC(I,2,K)-CC(I,5,K)
        TI2 = CC(I,2,K)+CC(I,5,K)
        TI4 = CC(I,3,K)-CC(I,4,K)
        TI3 = CC(I,3,K)+CC(I,4,K)
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
        CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
        CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
        CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
        CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
        CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
        CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
        CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
     ENDDO         ! 103     CONTINUE
  ENDDO            ! 104     CONTINUE
  RETURN
END SUBROUTINE PASSF5

FUNCTION PIMACH (DUM)
  IMPLICIT REAL*8 (A-H,O-Z)
  PIMACH = 4.0D+0*ATAN(1.0D+0)
  RETURN
END FUNCTION PIMACH
