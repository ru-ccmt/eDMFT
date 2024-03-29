SUBROUTINE SPHBES(N,X,FJ)
  IMPLICIT NONE
  INTEGER            N
  DOUBLE PRECISION   X
  DOUBLE PRECISION   FJ(*)
  !..................................................................
  !   Version III-upper limit obtained from the expressions of
  !              Corbato and Uretsky using a limit of E sub M of
  !              2**-30 which is approximately 1.E-9
  !             summation property used to normalize results.
  !   additional factor added to starting value
  !   N is the maximum L to be calculated, min N=2
  !   X is the argument
  !   FJ is the array where the spherical bessel functions are to be
  !   placed in.
  !..................................................................
  !        Parameters
  DOUBLE PRECISION   XLIM,  HF,  ZERO,  ONE,  TWO,  TNHF,  FFT
  DOUBLE PRECISION   T25,  TN25,  TN50
  PARAMETER          (XLIM = 0.1D+0,  HF = 0.5D+0,  ZERO = 0.0D+0)
  PARAMETER          (ONE = 1.0D+0,  TWO = 2.0D+0,  TNHF = 10.5D+0)
  PARAMETER          (FFT = 15.0D+0,  T25 = 1.0D+25)
  PARAMETER          (TN25 = 1.0D-25,  TN50 = 1.0D-50)
  !        Local Scalars
  INTEGER            J, JJ, M, MM, NS
  DOUBLE PRECISION   CUFAC, FFN, FFO, FFP, FM, HFXSQ, SDR, SER, TA
  DOUBLE PRECISION   TWM, XI, XL, XLP
  !        External Subroutines
  EXTERNAL           OUTERR
  !        Intrinsic Functions
  INTRINSIC          ABS, SQRT
  IF (N .LT. 2) THEN
     WRITE (6,6000)
  ELSEIF (X .LT. ZERO) THEN
     WRITE (6,6010)
  ELSEIF (X .GT. XLIM) THEN
     CUFAC = 4.2D+0
     IF (X .LT. (N-2)) CUFAC = TNHF/(N+HF-X)
     NS = N + 5 + int(X*CUFAC)
     !        add additional factor
     NS = NS + int(FFT/(ONE+SQRT(X)))
     !
     IF ((NS + 2) .GT. 100) THEN
        WRITE (6,6020) X, NS
        IF (X .GT. FFT) GOTO 900
        NS = 98
     ENDIF
     FFO = ZERO
     FFN = TN25
     M = NS - 1
     XI = ONE/X
     FM = (M + M) + ONE
     SDR = FM*TN50
     !        LOOP ... IF (M .LE. N) EXIT ... ENDLOOP
     
10   CONTINUE
     FFP = FM*XI*FFN - FFO
     IF (ABS(FFP) .GE. T25) THEN
        SDR = SDR*TN50
        FFP = FFP*TN25
        FFN = FFN*TN25
     ENDIF
     SDR = SDR + (FM-TWO)*FFP*FFP
     FFO = FFN
     FFN = FFP
     IF (M .LE. N) GOTO 15
     M = M - 1
     FM = FM - TWO
     !        END LOOP
     GOTO 10
15   CONTINUE
     
     FJ(M) = FFN
     FJ(M+1) = FFO
     
     !        REPEAT ... UNTIL (M .LE. 1)
20   CONTINUE
     M = M - 1
     FM = FM - TWO
30   CONTINUE
     FJ(M) = FM*XI*FJ(M+1) - FJ(M+2)
     IF (ABS(FJ(M)) .GE. T25) THEN
        JJ = M + 1
        NS = N + 1
        DO J = JJ, NS
           FJ(J) = FJ(J)*TN25
        ENDDO
        SDR = SDR*TN50
        GOTO 30
     ENDIF
     SDR = SDR + (FM-TWO)*FJ(M)*FJ(M)
     IF (M .GT. 1) GOTO 20
     !        END REPEAT
     SER = ONE/SQRT(SDR)
     MM = N + 1
     DO M = 1, MM
        FJ(M) = FJ(M)*SER
     ENDDO
  ELSE
     HFXSQ = HF*X**2
     XL = ONE
     TWM = ONE
     M = 0
     !        REPEAT ... UNTIL (M .GT. N)
60   CONTINUE
     M = M + 1
     TA = XL
     TWM = TWM + TWO
     XL = XL/TWM
     TA = TA - XL*HFXSQ
     XLP = XL/(TWM+TWO)
     FJ(M) = TA + HF*XLP*HFXSQ*HFXSQ
     XL = XL*X
     IF (M .LE. N) GOTO 60
     !        END REPEAT
  ENDIF
  !
  RETURN
  !
  !        Error messages
  !
900 CALL OUTERR('SPHBES','Error')
  STOP 'SPHBES - Error'
  !
6000 FORMAT (1H1,' ERROR, N SHOULD BE ge 2 ')
6010 FORMAT (1H1,' ERROR, X SHOULD NOT BE NEGATIVE ')
6020 FORMAT (1X,'WARNING: FOR X =',G13.6,' BESSH WANTS TO START AT N=',I5)
END SUBROUTINE SPHBES
