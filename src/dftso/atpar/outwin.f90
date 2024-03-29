SUBROUTINE OUTWIN(A,B,Nodes,VAL,SLO,REL,V,Rnot,dh,jri,EH,FL,Z,NRAD)
  !         Integration of skalarrel. Schroedinger Equation
  ! ----------------------------------------------------------------
  !  Input:
  !    EH    Energy in Hartree
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
  !    A and B are the solution of the Schroedinger-like equation -a''+ (l(l+1)/r^2 -2*U(r)/r-enu)a=0
  !    Here a is normalize like Int[ a(r)**2, dr]=1
  !    VAL is u(r)/r at R_{MT} and SLO is derivative d(a/r) = ( da/dr - a/r)/r at R_MT
  ! ----------------------------------------------------------------
  IMPLICIT NONE
  !
  REAL*8,  intent(out) :: A(NRAD), B(NRAD)    ! solution with large and small component
  INTEGER, intent(out) :: nodes               ! number of nodes
  REAL*8,  intent(out) :: val, slo            ! value at the spherical boundary:
  LOGICAL, intent(in)  :: REL
  REAL*8,  intent(in)  :: V(NRAD)
  REAL*8,  intent(in)  :: Rnot, dh, EH, Z, FL ! radial grid-points, log. increment, Energy in Rydbergs, nuclear charge
  INTEGER, intent(in)  :: jri, Nrad           ! Number of radial mesh points, Angular momentum l
  !    for large component u(r)=g(r,1) with psi=(u/r)*ylm                                                                                                                         
  !           val is u(rmax), i.e. (rmax * radial w.f.)                                                                                                                           
  !           slo is radial derivative of g at rmax (rmax * radial w.f.)'                                                                                                         
  !                                                                                                                                                                               
  ! Note: if r = b(exp(a*z)-1) then Jacobian dr=a(r+b) dz                                                                                                                         
  !                                                                                                                                                                               
  ! locals
  INTEGER :: iiij, K
  REAL*8  :: Rnet(NRAD)
  REAL*8  :: D(2,3), E, C, FLLP1, R83SQ, R1, R2, R3, H83, G0, AA, X, Y, U, ZZ, DRDI
  REAL*8  :: DG1, DG2, DG3, DF1, DF2, DF3, F0, S, SF, B1, B2, DET, PHI, R
  !print *, 'Inside OUTWIN r0=', Rnot, 'Eh=', EH, 'FL=', FL, 'Z=', Z, 'NRAD=', NRAD
  !
  !print *, '**********outwin******'
  !print *, 'rel=', REL
  !print *, 'V(0)=', V(1)
  !print *, 'Rnot=', Rnot
  !print *, 'DH=', DH
  !print *, 'JRI=', JRI
  !print *, 'EH=', EH
  !print *, 'FL=',FL
  !print *, 'Z=', Z
  !print *, 'nrad=', nrad
  !print *, 'VAL=', VAL
  !print *, 'SLO=', SLO
  !print *, 'Nodes=', Nodes
  !
  E=EH*2.d0  ! Converts from Hartree to Rydbergs
  A(:)=0.d0
  B(:)=0.d0
  !
  DO iiij=1,JRI
     Rnet(iiij)=RNOT*(exp(DH*(iiij-1.d0)))
  ENDDO
  !C
  Nodes = 0
  ZZ = Z + Z
  C = 2.d0*137.0359895d0
  !
  if(.not.rel) C=1e10   !  set to infinity
  !
  FLLP1 = FL*(FL + 1.d0)
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
  DO K = 1,3
     R = RNET(K)
     DRDI = DH*R
     A(K) = (R**S)*G0
     B(K) = (R**SF)*F0
     D(1,K) = DRDI*A(K)*S/R
     D(2,K) = DRDI*B(K)*SF/R
  ENDDO
  !
  !
  DG1 = D(1,1)
  DG2 = D(1,2)
  DG3 = D(1,3)
  DF1 = D(2,1)
  DF2 = D(2,2)
  DF3 = D(2,3)
  DO  K = 4, JRI
     R = RNET(K)
     DRDI = DH*R
     !
     !       Factor two in front of V because of Hartree-Rdyberg units!
     !
     PHI = (E - 2.d0*V(K)/R)*DRDI/C
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
  !
  DO iiij=1,JRI
     B(iiij)=B(iiij)*c/2.d0
  ENDDO
  !
  VAL = A(JRI)/RNET(JRI)
  SLO = DG3/(DH*RNET(JRI))
  SLO = (SLO-VAL)/RNET(JRI) 
  RETURN
END SUBROUTINE OUTWIN
