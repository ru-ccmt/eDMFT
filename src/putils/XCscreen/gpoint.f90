! @Copyright 2007 Kristjan Haule

SUBROUTINE tfpoint(thetas,phis,tweights,ln)
  !************************************************!
  !     GENERATES POINTS ON A UNIT SPHERE          !
  !************************************************!
  IMPLICIT NONE
  !
  REAL*8, intent(out)  :: thetas(ln+1), phis(2*ln+1), tweights(ln+1)! spt(2,(ln+1)*(2*ln+1)),weight((ln+1)*(2*ln+1))
  INTEGER, intent(in)  :: ln
  !
  REAL*8 :: dsum
  REAL*8 :: fake_shift, phi, pi
  INTEGER :: i, j, lg
  !
  PI=4.D0*ATAN(1.D0)
  lg = ln+1
  !
  CALL GRULE(lg,thetas,tweights) ! Gauss points in the interval [1...0]
  !
  do i=1,(lg+1)/2                ! adding symmetrical points to [0,...-1]
     j = lg+1-i
     thetas(j)=-thetas(i)
     tweights(j)=tweights(i)
  enddo
  do i=1,lg
     thetas(i) = dacos(thetas(i))
  enddo
  !
  fake_shift=EXP(-4.D0)  ! small shift such that we do not start at phi=0
  DO i=0,2*lg-2       ! 2*l-1 points in phi direction
     PHI=PI*(2*dble(i)/dble(2*lg-1) + fake_shift)
     phis(i+1)=PHI
  ENDDO
  !
  dsum = sum(tweights)
  tweights(:) = tweights(:)*2./dsum
  !write(6,*)'Gauss-Legendre grid of ',lg,'x',2*lg-1
  return
end SUBROUTINE tfpoint

SUBROUTINE apoint(spt,weight,ln)
  !************************************************!
  !     GENERATES POINTS ON A UNIT SPHERE          !
  !************************************************!
  IMPLICIT NONE
  !
  REAL*8, intent(out)  :: spt(2,(ln+1)*(2*ln+1)),weight((ln+1)*(2*ln+1))
  INTEGER, intent(in)  :: ln
  !
  REAL*8 :: XX((ln+1)*(2*ln+1))
  REAL*8 :: dsum
  REAL*8 :: fake_shift, phi, pi
  REAL*8 :: weight0(ln+1)
  INTEGER :: i, j, index, lg
  !
  PI=4.D0*ATAN(1.D0)
  lg = ln+1
  !
  CALL GRULE(lg,XX,WEIGHT) ! Gauss points in the interval [0...1]
  !
  DO i=1,(lg+1)/2
     j = lg+1-i
     XX(j)=-XX(i)
     WEIGHT(j)=WEIGHT(i)
  ENDDO
  weight0(:) = WEIGHT(:lg)
  !
  fake_shift=EXP(-4.D0)  ! small shift such that we do not start at phi=0
  INDEX=0
  DO j=1,lg            ! over all theta points
     DO i=0,2*lg-2       ! 2*l-1 points in phi direction
        PHI=PI*(2*dble(i)/dble(2*lg-1) + fake_shift)
        !INDEX=lg*i+j
        INDEX=INDEX+1
        WEIGHT(INDEX)=weight0(j)
        SPT(1,INDEX)=dacos(XX(j))
        SPT(2,INDEX)=PHI
     ENDDO
  ENDDO
  !
  dsum = sum(weight)
  weight(:) = weight(:)*4.D0*PI/dsum
  write(6,*)'Gauss-Legendre grid of ',lg,'x',2*lg-1
  return
end SUBROUTINE apoint

SUBROUTINE gpoint(spt,weight,ln)
  !************************************************!
  !     GENERATES POINTS ON A UNIT SPHERE          !
  !************************************************!
  IMPLICIT NONE
  !
  REAL*8, intent(out)  :: spt(3,(ln+1)*(2*ln+1)),weight((ln+1)*(2*ln+1))
  INTEGER, intent(in)  :: ln
  !
  REAL*8 :: XX((ln+1)*(2*ln+1))
  REAL*8 :: dsum
  REAL*8 :: fake_shift, phi, pi, rxy, z
  INTEGER :: i, j, index, lg
  !
  PI=4.D0*ATAN(1.D0)
  lg = ln+1
  !
  CALL GRULE(lg,XX,WEIGHT) ! Gauss points in the interval [0...1]
  !
  DO i=1,(lg+1)/2
     j = lg+1-i
     XX(j)=-XX(i)
     WEIGHT(j)=WEIGHT(i)
  ENDDO
  !
  fake_shift=EXP(-4.D0)  ! small shift such that we do not start at phi=0
  DO j=1,lg            ! over all theta points
     Z=XX(j)             ! z   = cos(theta)
     RXY=sqrt(1.D0-Z*Z)  ! rxy = sin(theta)
     DO i=0,2*lg-2       ! 2*l-1 points in phi direction
        PHI=PI*(2*dble(i)/dble(2*lg-1) + fake_shift)
        INDEX=lg*i+j
        WEIGHT(INDEX)=WEIGHT(j)
        SPT(1,INDEX)=COS(PHI)*RXY
        SPT(2,INDEX)=SIN(PHI)*RXY
        SPT(3,INDEX)=Z
     ENDDO
  ENDDO
  !
  dsum = sum(weight)
  weight(:) = weight(:)*4.D0*PI/dsum
  write(6,*)'Gauss-Legendre grid of ',lg,'x',2*lg-1
  return
end SUBROUTINE gpoint

SUBROUTINE GRULE(N,X,W)
  IMPLICIT NONE
  !
  !     DETERMINES THE (N+1)/2 NONNEGATIVE POINTS X(I) AND
  !     THE CORRESPONDING WEIGHTS W(I) OF THE N-POINT
  !     GAUSS-LEGENDRE INTEGRATION RULE, NORMALIZED TO THE
  !     INTERVAL \-1,1\. THE X(I) APPEAR IN DESCENDING ORDER.
  !
  !     THIS ROUTINE IS FROM 'METHODS OF NUMERICAL INTEGRATION',
  !     P.J. DAVIS AND P. RABINOWITZ, PAGE 369.
  !
  INTEGER, intent(in) :: N
  REAL*8, intent(out) :: X(*), W(*)
  !
  REAL*8  :: D1, d2pn, d3pn, d4pn, den, dp, dpn, e1, fx, h, p, pi, pk, pkm1, pkp1, t, t1, u, v, x0
  INTEGER :: i, it, k, m
  !
  PI=4.D0*ATAN(1.D0)
  M=(N+1)/2
  E1=N*(N+1)
  DO I=1,M         ! 1
     T=(4*I-1)*PI/(4*N+2)
     X0=(1.D0-(1.D0-1.D0/N)/(8.D0*N*N))*COS(T)
     !--->    ITERATE ON THE VALUE  (M.W. JAN. 1982)
     DO IT=1,3     ! 2
        PKM1=1.D0
        PK=X0
        DO K=2,N   ! 3
           T1=X0*PK
           PKP1=T1-PKM1-(T1-PKM1)/K+T1
           PKM1=PK
           PK=PKP1
        ENDDO      ! 3
        DEN=1.D0-X0*X0
        D1=N*(PKM1-X0*PK)
        DPN=D1/DEN
        D2PN=(2.D0*X0*DPN-E1*PK)/DEN
        D3PN=(4.D0*X0*D2PN+(2.D0-E1)*DPN)/DEN
        D4PN=(6.D0*X0*D3PN+(6.D0-E1)*D2PN)/DEN
        U=PK/DPN
        V=D2PN/DPN
        H=-U*(1.D0+.5D0*U*(V+U*(V*V-U*D3PN/(3.D0*DPN))))
        P=PK+H*(DPN+.5D0*H*(D2PN+H/3.D0*(D3PN+.25D0*H*D4PN)))
        DP=DPN+H*(D2PN+.5D0*H*(D3PN+H*D4PN/3.D0))
        H=H-P/DP
        X0=X0+H
     ENDDO
     X(I)=X0
     FX=D1-H*E1*(PK+.5D0*H*(DPN+H/3.D0*(D2PN+.25D0*H*(D3PN+.2D0*H*D4PN))))
     W(I)=2.D0*(1.D0-X(I)*X(I))/(FX*FX)
  ENDDO
  IF(M+M.GT.N) X(M)=0.D0
  RETURN
END SUBROUTINE GRULE
