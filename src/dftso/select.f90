SUBROUTINE SELECT(L,EI,DE,EMAIN,REL,VR,RNOT,DX,JRI,ZZ)
  use MPI, only: Qprint
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  !        Arguments
  INTEGER            JRI, L
  DOUBLE PRECISION   DE, DX, EI, RNOT, ZZ
  DOUBLE PRECISION   VR(*)
  LOGICAL            REL
  CHARACTER*4        EMAIN
  !..................................................................
  !   automatic energy search:
  !      finds energy parameter by looking where value and slope of
  !      wavefunction changes sign  (E-top:    u  changes sign)
  !                                 (E-bottom: u' changes sign)
  !..................................................................
  !   Local Scalars
  INTEGER            NODE
  DOUBLE PRECISION   DDNOLD, DNOLD, DU, DUDN, DUPOLD, DUTST, DUUP
  DOUBLE PRECISION   E1, E2, EIDN, EIUP, FL, U, UDN, UPOLD, UTST
  DOUBLE PRECISION   UUP
  CHARACTER*67       ERRMSG
  !        External Subroutines
  EXTERNAL           OUTERR
  !        Intrinsic Functions
  INTRINSIC          MAX
  !
  E1 = -100.0D+0
  E2 = -100.0D+0
  if (de.lt.1.d-8) goto 600
  FL = L
  EI = EI*0.5D+0
  DE = DE*0.5D+0
  CALL DIRACOUT_old(REL,VR,RNOT,DX,JRI,EI,L,1,U,DU,NODE,ZZ)
  EIUP = EI
  EIDN = EI
  UPOLD = U
  DUPOLD = DU
  DNOLD = U
  DDNOLD = DU
10 CONTINUE
      EIUP = EIUP + DE
      CALL DIRACOUT_old(REL,VR,RNOT,DX,JRI,EIUP,L,1,UUP,DUUP,NODE,ZZ)
      UTST = UPOLD*UUP
      DUTST = DUPOLD*DUUP
      !      WRITE(6,99)L,EIUP,UTST,DUTST
      !   99 FORMAT(' L,EIUP,UTST,DUTST',I2,F8.4,2E12.3)
      UPOLD = UUP
      DUPOLD = DUUP
      IF (UTST .LT. 0.0D+0) THEN
         E2 = EIUP
         IF (E1 .GT. -30.0D+0) GOTO 30
      ENDIF
      IF (DUTST .LT. 0.0D+0) THEN
         E1 = EIUP
         IF (E2 .GT. -30.0D+0) GOTO 30
      ENDIF
      IF ((E1 .LT. -30.0D+0) .AND. (EIUP .LT. 2.5D+0)) THEN
   20     CONTINUE
          EIDN = EIDN - DE
          CALL DIRACOUT_old(REL,VR,RNOT,DX,JRI,EIDN,L,1,UDN,DUDN,NODE,ZZ)
          UTST = DNOLD*UDN
          DUTST = DDNOLD*DUDN
          !         WRITE(6,98)L,EIDN,UTST,DUTST
          !   98    FORMAT(' L,EIDN,UTST,DUTST',I2,F8.4,2E12.3)
          DNOLD = UDN
          DDNOLD = DUDN
          IF (UTST .LT. 0.0D+0) THEN
             E2 = EIDN
             IF (E1 .GT. -30.0D+0) GOTO 30
          ENDIF
          IF (DUTST .LT. 0.) THEN
             E1 = EIDN
             IF (E2 .GT. -30.0D+0) GOTO 30
          ENDIF
          IF (E2 .LT. -30.0D+0) THEN
             GOTO 10
          ELSEIF (EIDN .GT. (EI - 3.0D+0)) THEN
             GOTO 20
          ENDIF
       ELSEIF (EIUP .LT. 0.5D+0) THEN
          GOTO 10
       ENDIF
30     CONTINUE
       EI = EI*2.0D+0
       E1 = E1*2.0D+0
       E2 = E2*2.0D+0
       IF ((E1 .LT. -60.0D+0) .AND. (E2 .LT. -60.0D+0)) THEN
          IF (EMAIN .EQ. 'STOP') THEN
             !        Error: no energy limits found
             GOTO 900
          ELSE
            EI = 2.0D+0
         ENDIF
      ELSEIF (E2 .LT. -60.0D+0) THEN
         IF (EMAIN .EQ. 'STOP') THEN
            !        Error: no energy limits found
            GOTO 900
         ELSE
            EI = MAX(E1,EI)
         ENDIF
      ELSEIF (E1 .LT. -60.0D+0) THEN
         !        Error: no energy limits found
         GOTO 900
      ELSE
         EI = (E1+E2)*0.5D+0
      ENDIF
      !
      RETURN
      !
600   continue
      if (Qprint) WRITE (21,6000) L, EI, E1, E2
      return
      !        Error messages
900   WRITE (ERRMSG,9000) L
      CALL OUTERR('SELECT',ERRMSG)
      WRITE (ERRMSG,9010) E1, E2
      CALL OUTERR('SELECT',ERRMSG)
      STOP 'SELECT - Error'
      !
6000  FORMAT (10X,'E(',I2,')=',F10.4,3X,'E(BOTTOM)=',F9.3,3X,'E(TOP)=',F9.3)
9000  FORMAT ('no energy limits found for L=',I2)
9010  FORMAT ('E-bottom ',F10.5,3X,'E-top ',F10.5)
      !        End of 'SELECT'
      END
