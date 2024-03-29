SUBROUTINE SELECT(L,EI,DE,EMAIN,REL,VR,RNOT,DX,JRI,ZZ,jatom,efermi)
  use radial, only: A,B
  use param
  use mpi, only: Qprint, stop_MPI
  IMPLICIT NONE
  !        Arguments
  !include 'param.inc'
  REAL*8, intent(out)    :: Ei
  REAL*8, intent(inout)  :: dE
  INTEGER, intent(in)    :: L, jri, jatom
  REAL*8, intent(in)     :: Rnot, dx, ZZ, efermi
  CHARACTER*4, intent(in):: EMAIN
  LOGICAL, intent(in)    :: REL
  REAL*8, intent(in)     :: VR(*)
  ! locals
  REAL*8   :: Rnet(NRAD)
  LOGICAL  :: UP, DEBUG, ShowIter
  REAL*8   :: DExit, EI_High, EI_low, EIout, Reduce, Upper
  INTEGER  :: ifail, iiij, numcalls, NumIter
  !..................................................................
  !   automatic energy search:
  !      finds energy parameter by looking where value and slope of
  !      wavefunction changes sign  (E-top:    u  changes sign)
  !                                 (E-bottom: u' changes sign)
  !..................................................................
  !        Local Scalars
  INTEGER            NODE, NODE1, NODE2
  DOUBLE PRECISION   DDNOLD, DNOLD, DU, DUDN, DUPOLD, DUUP
  DOUBLE PRECISION   E1, E2, EIDN, EIUP, FL, U, UDN, UPOLD
  DOUBLE PRECISION   UUP
  CHARACTER*67       ERRMSG
  !        External Subroutines
  EXTERNAL           OUTERR, OUTWINC
  !
  REAL*8 :: SearchUP  ! External function
  REAL*8 :: SearchDN  ! External function
  !
  DEBUG=.false.
  !DEBUG=.true.
  ShowIter=.false.
  numcalls=0
  !     Setup RMT mesh. Moved out of loop for efficiency
  do iiij=1,JRI
     RNET(iiij)=RNOT*(exp(DX*(iiij-1.d0)))
  enddo

  ! Be careful with these hard coded hacks. Debug....
  Upper = 0.5  
  if(efermi.gt.0.8) Upper = efermi*0.5D0 + 0.25D0 
  if(efermi.gt.0.8.and.l.gt.1) Upper = efermi*0.5D0 + 0.5D0 

  
  DE = max(DE,0.02D0)
  NODE1 = -1
  NODE2 = -1
  Reduce = 0.5D0
  FL = L
  EI = EI*0.5D+0
  DE = DE*0.5D+0
  CALL OUTWINC(A,B,NODE,U,DU,REL,VR,Rnet,DX,JRI,EI,FL,ZZ,nrad)        ! value at Ei
  
  numcalls=numcalls+1
  EIUP = EI
  EIDN = EI
  E1 = -100.0D+0
  E2 = -100.0D+0
  UPOLD = U
  DUPOLD = DU
  DNOLD = U
  DDNOLD = DU
  DO while (EIup.LT.upper)
10   CONTINUE
     EIUP = EIUP + DE
     CALL OUTWINC(A,B,NODE,UUP,DUUP,REL,VR,Rnet,DX,JRI,EIUP,FL,ZZ,nrad)  ! value at EI+dE
     numcalls=numcalls+1
     IF ((UPOLD*UUP) .LT. 0.0D+0) THEN   ! change of sign for U(R) occured between EIUP and EIUP+DE
        E2 = EIUP
        NODE2=NODE
        IF (E1 .GT. -50.0D+0) EXIT       ! if zero of U(R) is bracketed and E1 is not too small, we can go and optimize
     ENDIF
     UPOLD = UUP
     IF ((DUPOLD*DUUP) .LT. 0.0D+0) THEN ! change of sign for dU/dR occured between EIUP and EIUP+DE
        E1 = EIUP
        NODE1 = NODE
        IF (E2 .GT. -50.0D+0) EXIT      ! if zero of dU/dr(R) is bracketed and E2 is not too small, we can go and optimize
     ENDIF
     DUPOLD = DUUP
     
     IF ((E1 .LT. -50.0D+0) .AND. (EIUP .LT. 2.5D+0)) THEN
20      CONTINUE
         
        EIDN = EIDN - DE
        CALL OUTWINC(A,B,NODE,UDN,DUDN,REL,VR,Rnet,DX,JRI,EIDN,FL,ZZ,nrad)
        numcalls=numcalls+1
        IF ((DNOLD*UDN) .LT. 0.0D+0) THEN  ! change of sign for U(R) occured between EIDN and EIDN+DE
           E2 = EIDN
           NODE2 = NODE
           IF (E1 .GT. -50.0D+0) EXIT
        ENDIF
        DNOLD = UDN
        IF ((DDNOLD*DUDN) .LT. 0.) THEN    ! change of sign for dU/dR occured between EIDN and EIDN+DE
           E1 = EIDN
           NODE1 = NODE
           IF (E2 .GT. -50.0D+0) EXIT
        ENDIF
        DDNOLD = DUDN
        
        IF (E2 .LT. -50.0D+0) GOTO 10
        IF (EIDN .GT. (EI - 3.0D+0)) GOTO 20
     ELSEIF (EIUP .LT. upper) THEN
        GOTO 10
     ENDIF
  ENDDO

  ! Here we bracket the zero of U and dU/dr
  if (showiter) write(6,66)0,numcalls,2*E1,2*E2,2*DE
  DExit   = 1D-8
  NumIter = 20 
  DO IIIJ=1,NumIter
     !
     IF( E2 .gt. -100.D0 ) then
        ! Searching for U(R)=0 from (E2-DE,E2+5*DE)
        EI_Low = E2-DE
        EI_High= E2+DE*5.D0
        Up     =.true.
        EIout = SearchUP(A,B,REL,VR,RNET,DX,JRI,FL,ZZ ,EI_Low, EI_High, Up, DE, numcalls, ifail, debug, nrad)
        IF(ifail.eq.0) E2 = EIOUT     ! Search succeded. We have U(R)=0 at E2
        if(debug)write(6,*)'Returned UP ',E2,ifail,EIOUT
     ENDIF
     IF( E1 .gt. -100.D0 ) then
        ! Searching for dU/dr(R)=0 from (E1-DE,E1+5*DE)
        EI_Low = E1-DE
        EI_High= E1+DE*5.D0
        Up     =.true.
        EIout = SearchDN(A,B,REL,VR,RNET,DX,JRI,FL,ZZ,EI_Low, EI_High, Up, DE, numcalls, ifail, debug, nrad)
        IF(ifail.eq.0) E1=EIOUT     ! Search succeded. We have dU/dr(R)=0 at E1
        if(debug)write(6,*)'Returned DN ',E1,EIOUT,ifail
     ENDIF
     !     End of Search. Now reduce the step size and iterate
     if(showiter) write(6,66)iiij,numcalls,2*E1,2*E2,2*DE
66   format('Pass ',i4,' Iterations ',i4,' E1 & E2 ',2ES13.5,' Step ',ES12.4,' (Ryd)')
     !
     DE=DE*Reduce
     if(DE.lt.DExit)exit                     ! Termination condition
  ENDDO
  
  EI = EI*2.0D+0
  E1 = E1*2.0D+0
  E2 = E2*2.0D+0
  IF ((E1 .LT. -160.0D+0) .AND. (E2 .LT. -160.0D+0)) THEN
     IF (EMAIN .EQ. 'STOP') THEN
        !        Error: no energy limits found
        GOTO 900
     ELSE
        EI = 2.0D+0                    ! both E1 and E2 very negative and Emain='C' => Ei=2.0
     ENDIF
  ELSEIF (E2 .LT. -160.0D+0) THEN
     IF (EMAIN .EQ. 'STOP') THEN
        !        Error: no energy limits found
        GOTO 900
     ELSE
        EI = MAX(E1,EI)               ! E2 very negative, but E1 not, and Emain='C' => Ei=max(E1,Ei)
     ENDIF
  ELSEIF (E1 .LT. -160.0D+0) THEN     ! E1 very negative but E2 not (case of dU/dr(E1)=0) => problem!
     !        Error: no energy limits found
     GOTO 900
  ELSE
     EI = (E1+E2)*0.5D+0              ! both E1 and E2 are normal => set the average
  ENDIF
  if (Qprint) then
     WRITE (6,6000)  l,jatom,L, EI, E1, E2, Node1, Node2, numcalls
     WRITE (21,6000) l,jatom,L, EI, E1, E2, Node1, Node2, numcalls
  endif
  
  RETURN
  !        Error messages
900 WRITE (ERRMSG,9000) jatom, L
  CALL OUTERR('SELECT',ERRMSG)
  WRITE (ERRMSG,9010) E1, E2
  CALL OUTERR('SELECT',ERRMSG)
  call stop_MPI
  STOP 'SELECT - Error'
  !
6000 FORMAT (':E',i1,'_',i4.4,':',1X,'E(',I2,')=',F10.4,3X,'E(BOTTOM)=',F9.3,3X,'E(TOP)=',F9.3,2i3,i6)
9000 FORMAT ('no energy limits found for atom ',i3,'  L=',I2)
9010 FORMAT ('E-bottom ',F10.5,3X,'E-top ',F10.5)
END SUBROUTINE SELECT


REAL*8 Function SearchDN(A,B,REL,VR,RNET,DX,JRI,FL,ZZ ,EI_Low, EI_High, Up, DE, numcalls, ifail, debug, nrad)
  ! Finds the sign chnage of the radial derivative of U
  IMPLICIT NONE
  DOUBLE PRECISION, intent(in) :: VR(*), RNET(NRAD)
  LOGICAL,          intent(in) :: REL, Up, debug
  REAL*8,           intent(in) :: DE, DX, EI_low, EI_high, FL, ZZ
  INTEGER,          intent(in) :: nrad, jri
  INTEGER,          intent(inout):: numcalls
  REAL*8,           intent(out):: A(nrad), B(nrad)
  INTEGER,          intent(out):: ifail
  REAL*8 :: DDE, EA, EB, DDNOLD, DUDN, EI, UDN, EIout
  INTEGER:: node
  IFAIL=0
  IF( UP )then
     DDE=DE
     EA=EI_LOW
     EB=EI_High
  ELSE
     DDE=-DE
     EA=EI_High
     EB=EI_Low
  ENDIF
  CALL OUTWINC(A,B,NODE,UDN,DDNOLD,REL,VR,Rnet,DX,JRI,EA,FL,ZZ,nrad)
  if(debug)write(6,*)'DN',EA,DDNOLD
  EI = EA
  DO WHILE(EI.LT.EB)
     EI = EI + DDE
     CALL OUTWINC(A,B,NODE,UDN,DUDN,REL,VR,Rnet,DX,JRI,EI,FL,ZZ,nrad)
     numcalls=numcalls+1
     if(debug)write(6,*)'DN',EI,DUDN
     !
     IF ((DDNOLD*DUDN) .LT. 0.0D+0) then  !       Has the sign changed ?
        EIOUT=min(EI,EI-DDE)
        if(debug)write(6,*)'Sign Dn ',EIOUT,DDNOLD
        SearchDN=EIout
        return
     endif
     DDNOLD = DUDN
  enddo
  SearchDN=0.d0
  IFAIL=1
  return
end Function SearchDN

REAL*8 Function SearchUP(A,B,REL,VR,RNET,DX,JRI,FL,ZZ,EI_Low, EI_High, Up, DE, numcalls, ifail, debug, nrad )
  ! Finds the sign chnage of the function U
  IMPLICIT NONE
  DOUBLE PRECISION, intent(in) :: VR(*), RNET(NRAD)
  LOGICAL,          intent(in) :: REL, Up, debug
  REAL*8,           intent(in) :: DE, DX, EI_low, EI_high, FL, ZZ
  INTEGER,          intent(in) :: nrad, jri
  INTEGER,          intent(inout):: numcalls
  REAL*8,           intent(out):: A(nrad), B(nrad)
  INTEGER,          intent(out):: ifail
  ! locals
  REAL*8 :: DDE, EA, EB, DDNOLD, DUDN, EI, UPOLD, UUP, EIOUT
  INTEGER:: node
  IFAIL=0
  IF( UP )then  ! Energy will be increasing from Ei_low to Ei_high in increments DE
     DDE=DE
     EA=EI_LOW
     EB=EI_High
  ELSE
     DDE=-DE   ! Energy will be decreasing from Ei_high to Ei_low in decrements DE
     EA=EI_High
     EB=EI_Low
  ENDIF
  CALL OUTWINC(A,B,NODE,UPOLD,DDNOLD,REL,VR,Rnet,DX,JRI,EA,FL,ZZ,nrad)
  if(debug)write(6,*)'UP',EA,UPOLD
  EI = EA
  DO WHILE(EI.LT.EB)
     EI = EI+DDE 
     CALL OUTWINC(A,B,NODE,UUP,DUDN,REL,VR,Rnet,DX,JRI,EI,FL,ZZ,nrad)
     numcalls=numcalls+1
     if(debug)write(6,*)'UP',EI,UUP
     !
     IF ( (UPOLD*UUP) .LT. 0.0D+0) then  ! Has the sign changed ?
        EIOUT=min(EI,EI-DDE)
        if(debug)write(6,*)'Sign Up',EIOUT,UUP
        SearchUP = EIOUT
        RETURN
     endif
     UPOLD = UUP
  ENDDO
  SearchUp=0.d0
  IFAIL=1
  return
end Function SearchUP
