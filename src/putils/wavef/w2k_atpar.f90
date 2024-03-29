! @Copyright 2007 Kristjan Haule
SUBROUTINE OUTWIN(A, B, REL, V, Rnet, Nr, dh, ER, FL, Z, VAL,SLO,Nodes) 
  !         Integration of skalarrel. Schroedinger Equation
  ! ----------------------------------------------------------------
  !  Input:
  !    ER    Energy in Rydbergs
  !    FL    Angular momentum 
  !    Z     Nuclear charge
  !    V     rad.sym. Potential in Hartree
  !    RNOT  first radial grid point
  !    DH    log. increment
  !    Nr    Number of radial mesh points
  !  Output:
  !    VAL,SLO:  Wave function and gradient at the spherical boundary
  !    Nodes:    number of nodes
  !
  !    Rydberg units
  !
  ! ----------------------------------------------------------------
  IMPLICIT NONE
  REAL*8, intent(out) :: A(Nr), B(Nr)    ! solution with large and small component
  LOGICAL, intent(in) :: rel             ! real of complex
  REAL*8, intent(in)  :: V(Nr)           ! radial symmetric potential
  REAL*8, intent(in)  :: Rnet(Nr), dh, ER, Z ! radial grid-points, log. increment, Energy in Rydbergs, nuclear charge
  INTEGER, intent(in) :: Nr, fl          ! Number of radial mesh points, Angular momentum l
  INTEGER, intent(out):: nodes           ! number of nodes
  REAL*8, intent(out) :: val, slo        ! value at the spherical boundary:
  !    for large component u(r)=g(r,1) with psi=(u/r)*ylm
  !           val is u(rmax), i.e. (rmax * radial w.f.)
  !           slo is radial derivative of g at rmax (rmax * radial w.f.)'
  !
  ! Note: if r = b(exp(a*z)-1) then Jacobian dr=a(r+b) dz
  !
  ! locals
  INTEGER :: K, ij
  REAL*8  :: D(2,3), C, FLLP1, R83SQ, R1, R2, R3, H83, G0, AA, X, Y, U, ZZ, DRDI
  REAL*8  :: DG1, DG2, DG3, DF1, DF2, DF3, F0, S, SF, B1, B2, DET, PHI, R
  !
  !     Hartree in Ryd
  !E=EH*2.d0
  !
  A(:)=0.0
  B(:)=0.0
  Nodes = 0
  ZZ = Z + Z
  C = 2.d0*137.0359895d0
  !
  if(.not.rel) C=1e10  ! set to infinity
  ! angular momentum part
  FLLP1 = FL*(FL + 1.d0)
  ! some coefficients for integration
  R83SQ = 64.D0/9.D0
  R1 = 1.D0/9.D0
  R2 = -5.D0*R1
  R3 = 19.D0*R1
  H83 = 8.D0/3.D0
  !
  G0 = 1.d0
  IF (Z .LT. 0.9D0) THEN
     S = FL+1.d0
     SF = FL
     F0 = FL/C
  ELSE
     AA = ZZ/C
     S = DSQRT(FLLP1 + 1.D0 - AA*AA)
     SF = S
     F0 = G0*(S - 1.D0)/AA
  ENDIF
  DO  K = 1,3 ! 
     R = RNET(K)
     DRDI = DH*R
     A(K) = (R**S)*G0
     B(K) = (R**SF)*F0
     D(1,K) = DRDI*A(K)*S/R
     D(2,K) = DRDI*B(K)*SF/R
  ENDDO
  !
  DG1 = D(1,1)
  DG2 = D(1,2)
  DG3 = D(1,3)
  DF1 = D(2,1)
  DF2 = D(2,2)
  DF3 = D(2,3)
  DO K = 4, Nr
     R = RNET(K)
     DRDI = DH*R
     !       Factor two because V is in Hartree-Rydberg !
     PHI = (ER - 2.d0*V(K)/R)*DRDI/C
     U = DRDI*C + PHI
     X = -DRDI/R
     Y = -FLLP1*X*X/U + PHI
     DET = R83SQ - X*X + U*Y
     B1 = A(K-1)*H83 + R1*DG1 + R2*DG2 + R3*DG3
     B2 = B(K-1)*H83 + R1*DF1 + R2*DF2 + R3*DF3
     A(K) = (B1*(H83-X) + B2*U)/DET
     B(K) = (B2*(H83+X) - B1*Y)/DET
     IF (A(K)*A(K-1) .LT. 0D0) Nodes = Nodes + 1
     DG1 = DG2
     DG2 = DG3
     DG3 = U*B(K) - X*A(K)
     DF1 = DF2
     DF2 = DF3
     DF3 = X*B(K) - Y*A(K)
  ENDDO
  !
  DO ij=1,Nr 
     B(ij)=B(ij)*c/2.d0
  ENDDO 
  !
  VAL = A(Nr)/RNET(Nr)
  SLO = DG3/(DH*RNET(Nr))
  SLO = (SLO-VAL)/RNET(Nr) 
  RETURN
END SUBROUTINE OUTWIN

REAL*8 FUNCTION RINT13(rel,A1,B1,A2,B2,Nr,DX,R0)
  !PERFORMS RADIAL INTEGRALS
  !coded by D.D.KOELLING
  IMPLICIT NONE
  logical, intent(in):: rel
  REAL*8, intent(in) :: A1(Nr),B1(Nr),A2(Nr),B2(Nr)
  REAL*8, intent(in) :: DX, R0
  INTEGER, intent(in):: Nr
  !
  ! temporary
  REAL*8 :: D, S, CIN, R, R1, Z2, Z4, P1, P2
  INTEGER:: J, J1
  !
  CIN=1.d0/137.0359895d0**2
  IF(.NOT.REL) CIN=1E-22
  !
  D=EXP(DX)
  J=3-MOD(Nr,2)
  R=r0*(D**(J-1))
  J1=J-1
  R1=R/D
  Z4=0
  Z2=0
  DO
     Z4=Z4+R*(A1(J)*A2(J)+CIN*B1(J)*B2(J))                                 
     R=R*D
     J=J+1
     IF(J.GE.Nr) EXIT
     Z2=Z2+R*(A1(J)*A2(J)+CIN*B1(J)*B2(J))
     R=R*D
     J=J+1
  ENDDO
  P1=r0*(A1(1)*A2(1)+CIN*B1(1)*B2(1))
  P2=R1*(A1(J1)*A2(J1)+CIN*B1(J1)*B2(J1))
  S=2*Z2+4*Z4+R*(A1(J)*A2(J)+CIN*B1(J)*B2(J))+P2
  S=(DX*S+P1)/3.0D0
  IF(J1.GT.1) S=S+0.5D0*DX*(P1+P2)
  RINT13=S
  RETURN
END FUNCTION RINT13

REAL*8 FUNCTION RINT13g(C1,C2,A1,B1,A2,B2,Nr,DX_,R0_)
  ! C1=1, C2=[1.d0/137.0359895d0**2 or 0]
  IMPLICIT NONE
  !        but here C1 and C2 are arbitrary
  REAL*8, intent(in) :: C1, C2
  REAL*8, intent(in) :: A1(Nr), B1(Nr), A2(Nr), B2(Nr)
  REAL*8, intent(in) :: DX_, R0_
  INTEGER, intent(in):: Nr
  ! locals
  INTEGER::            J, J1
  REAL*8 ::   D, S, P1, P2, R, R1, Z2, Z4
  D = EXP(DX_)
  J = 3-MOD(Nr,2)
  J1=J-1
  R = r0_*(D**(J-1))
  R1 = R/D
  Z4 = 0
  Z2 = 0
  DO
     Z4=Z4+R*(C1*A1(J)*A2(J) + C2*B1(J)*B2(J))     ! Z4 ~ R*A*X (r)
     R=R*D
     J=J+1
     IF (J .GE. Nr) EXIT
     Z2 = Z2 + R*(C1*A1(J)*A2(J) + C2*B1(J)*B2(J))  ! Z2 ~ R*A*X (r+dr)   
     R=R*D
     J=J+1
  enddo
  P1 = R0_*(C1*A1(1)*A2(1) + C2*B1(1)*B2(1))
  P2 = R1*(C1*A1(J1)*A2(J1) + C2*B1(J1)*B2(J1))
  S = 2*Z2 + 4*Z4 + R*(C1*A1(J)*A2(J) + C2*B1(J)*B2(J)) + P2
  S = (DX_*S+P1)/3.0D+0
  IF (J1.GT.1) S=S+0.5D0*DX_*(P1+P2)
  RINT13g=S
  RETURN
END FUNCTION RINT13g

! Nr -> jrj(jatom)
SUBROUTINE ReadPotential(VR, filename, maxNr, Nr, nat)
!!! Reading potential
  REAL*8, intent(out)       :: VR(maxNr,nat)
  CHARACTER*200, intent(in) :: filename
  INTEGER, intent(in)       :: nat, maxNr, Nr(nat)
  !! temporary
  INTEGER :: jtape, iscf, idummy, j, jatom
  !
  jtape=18
  open(jtape,file=filename,status='old')
  READ(jtape, '(50X,I2,//)') ISCF

  do jatom=1,nat
     READ(jtape,1980)
     READ(jtape,2000) idummy
     READ(jtape,2031)
     READ(jtape,2022)(VR(J,jatom),J=1,Nr(jatom))
     READ(jtape,2031)
     READ(jtape,2030)
  enddo
  
  VR(:,:)=VR(:,:)/2.0D0  ! To Hartree's?
  close(jtape)
2022 FORMAT(3X,4E19.12)
2000 FORMAT(16X,I2//)
1980 FORMAT(3X)
2031 FORMAT(/)
2030 FORMAT(///)
END SUBROUTINE ReadPotential


SUBROUTINE ReadLinearizationE(E, filename, nat, lmax2)
!!! Reading linearization energies
  CHARACTER*200, intent(in) :: filename
  INTEGER, intent(in) :: LMAX2
  REAL*8, intent(out) :: E(0:LMAX2,nat)
  ! locals
  INTEGER :: itape, i
  REAL*8  :: elo
  !jatom = jatom+1  ! to convert from python to fortran
  itape=9
  open(itape,file=filename,status='old',form='unformatted')
  DO jatom=1,nat
     READ(itape) E(:,jatom)     ! linearization energies
     READ(itape) elo            ! local
  ENDDO
  close(itape)
END SUBROUTINE ReadLinearizationE


  ! NORM OF V(0,0)=V00(R)*R/SQRT(4.D0*PI)
  !DO i=1,jri
  !   Rx(i)=r0*(exp(dx*(i-1.d0)))
  !ENDDO
  ! dx <- dx(jatom)
  ! El <- E(l)

  ! DX <- DX(jatom)
  ! jrj <- jrj(jatom)
  ! r0 <- r0(jatom)
  ! Nrad, JRJ


SUBROUTINE ATPAR(A,B,Ae,Be,Aee,Bee,Pei,rel,l,VR,Rx,El,r0,dx,Zn,Nr)
  IMPLICIT NONE
  REAL*8, intent(out)  :: A(Nr), B(Nr)   ! large and small component of Sch. equation: u(r)
  REAL*8, intent(out)  :: Ae(Nr), Be(Nr) ! large and small component of derivative of Sch. equation: udot(r)
  REAL*8, intent(out)  :: Aee(Nr), Bee(Nr) ! large and small component of the second derivative of Sch. equation: udot(r)
  REAL*8, intent(out)  :: Pei            ! <udot|udot>
  !INTEGER, intent(out) :: Nodel          ! number of nodes of u(r)
  !REAL*8, intent(out)  :: dUv, dUve, dUvee  ! derivatives at the Rmt: dUv,dUve,dUvee of ul, dot{ul}, dot_dot{ul}
  LOGICAL, intent(in) :: rel             ! relativistic or not
  INTEGER, intent(in) :: l               ! orbital momentum
  REAL*8, intent(in)  :: VR(Nr),Rx(Nr)   ! radial potenial, radial mesh
  REAL*8, intent(in)  :: El              ! linearization energy for this l
  REAL*8, intent(in)  :: r0, dx, Zn      ! first point of the radial grid, dr for the grid, nuclear charge
  INTEGER, intent(in) :: Nr              ! number of radial points
  ! External
  INTERFACE 
     FUNCTION RINT13(rel,A1,B1,A2,B2,Nr,DX,R0)
       REAL*8 :: RINT13
       logical, intent(in):: rel
       REAL*8, intent(in) :: A1(Nr),B1(Nr),A2(Nr),B2(Nr)
       REAL*8, intent(in) :: DX, R0
       INTEGER, intent(in):: Nr
     END FUNCTION RINT13
  END INTERFACE
  ! temporary
  REAL*8  :: A1(Nr), B1(Nr), A2(Nr), B2(Nr)
  REAL*8  :: Uv1, dUv1, Uv2, dUv2, Uvee, Uv, dElR, ovlp, cross, Uve, dUv, dUve, dUvee
  INTEGER :: Nodel, i
  !
  !  CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE
  !  dElR IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN Rydbers
  dElR=4.0D-3
  
  CALL OUTWIN(A1, B1, rel, VR, Rx, Nr, dx, El-dElR, l, Zn, Uv1, dUv1, Nodel)
  ! Normalization
  ovlp = RINT13(rel, A1, B1, A1, B1, Nr, dx, r0)
  A1(:) = A1(:)/SQRT(ovlp)
  B1(:) = B1(:)/SQRT(ovlp)
  Uv1=Uv1/SQRT(ovlp)         ! Normalized value at Rmt and energy El-dElR
  dUv1=dUv1/SQRT(ovlp)
  
  CALL OUTWIN(A2, B2, rel, VR, Rx, Nr, dx, El+dElR, l, Zn, Uv2, dUv2, Nodel)
  ! Normalization
  ovlp = RINT13(rel, A2, B2, A2, B2, Nr, dx, r0)
  A2(:) = A2(:)/SQRT(ovlp)
  B2(:) = B2(:)/SQRT(ovlp)
  Uv2=Uv2/SQRT(ovlp)         ! Normalized value at Rmt and energy El+dElR
  dUv2=dUv2/SQRT(ovlp)

  ! Finite difference
  UvE = (Uv2-Uv1)/(2.*dElR)       ! energy derivative of the value at Rmt, i.e., udot(Rmt)
  dUvE= (dUv2-dUv1)/(2.*dElR)     ! d/dr udot(Rmt)
  Ae(:) = (A2(:)-A1(:))/(2.*dElR)
  Be(:) = (B2(:)-B1(:))/(2.*dElR)
  
  ! Schroedinger equation solver at El
  CALL OUTWIN(A, B, rel, VR, Rx, Nr, dx, El, l, Zn, Uv, dUv, Nodel)
  ovlp = RINT13(rel, A, B, A, B, Nr, dx, r0)
  A(:) = A(:)/SQRT(ovlp)
  B(:) = B(:)/SQRT(ovlp) 
  Uv=Uv/SQRT(ovlp)         ! Normalized value at Rmt and energy El+dElR
  dUv=dUv/SQRT(ovlp)
 
  ! Orthogonalization between u(r) and udot(r)
  CROSS = RINT13(rel, A, B, Ae, Be, Nr, dx, r0)
  Ae(:) = Ae(:)-CROSS*A(:)
  Be(:) = Be(:)-CROSS*B(:)
  UvE = UvE-CROSS*Uv
  dUvE=dUvE-CROSS*dUv
  
  ! PEI = <udor(r)|udot(r)>
  PEI = RINT13(rel, Ae, Be, Ae, Be, Nr, dx, r0)

  ! Finite difference
  Aee(:) = (A2(:)-2*A(:)+A1(:))/(2.*dElR**2)
  Bee(:) = (B2(:)-2*B(:)+B1(:))/(2.*dElR**2)
  Uvee = (Uv2-2*Uv+Uv1)/(2.*dElR**2)
  dUvee= (dUv2-2*dUv+dUv1)/(2.*dElR**2)
  ! Orthogonalization between udotdot(r) and u(r)
  CROSS = RINT13(rel, Aee, Bee, A, B, Nr, dx, r0)
  Aee(:) = Aee(:)-CROSS*A(:)
  Bee(:) = Bee(:)-CROSS*B(:)
  Uvee = Uvee - CROSS*Uv
  dUvee=dUvee - CROSS*dUv
  ! Orthogonalization between udotdot(r) and udot(r)
  CROSS = RINT13(rel, Aee, Bee, Ae, Be, Nr, dx, r0)
  CROSS = CROSS/PEI
  Aee(:) = Aee(:)-CROSS*Ae(:)
  Bee(:) = Bee(:)-CROSS*Be(:)
  Uvee = Uvee - CROSS*UvE
  dUvee=dUvee - CROSS*dUvE

  !print *, RINT13(rel, A, B, Ae, Be, Nr, dx, r0), RINT13(rel, A, B, Aee, Bee, Nr, dx, r0), RINT13(rel, Ae, Be, Aee, Bee, Nr, dx, r0)
  
  RETURN
END SUBROUTINE ATPAR
