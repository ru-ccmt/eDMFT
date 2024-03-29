SUBROUTINE OUTWIN(REL,V,RNOT,DH,JRI,EH,FL,VAL,SLO,Nodes,Z) 
  !         Integration der skalarrel. Schroedingergleichung
  ! ----------------------------------------------------------------
  !  Input:
  !    EH    Energie in Hartree 
  !    FL    Drehimpuls 
  !    Z     Kernladung
  !    V     rad.sym. Potential in Hartree
  !    RNOT  erster radialer Netzpunkt
  !    DH    log. Schrittweite
  !    JRI   Anzahl radialer Netzpunkte 
  !  Output:
  !    VAL,SLO:  Wellenfunktion und Steigung am Kugelrand
  !    Nodes:    Anzahl Knoten 
  !
  !    Rydberg Einheiten 
  !
  ! ----------------------------------------------------------------
  IMPLICIT REAL*8 (A-H,O-Z)
  !
  include 'param.inc'
  !
  logical rel
  DIMENSION D(2,3),V(NRAD),RNET(NRAD)
  COMMON /UHELP/   A(NRAD),B(NRAD),AP(NRAD),BP(NRAD),AE(NRAD),BE(NRAD)
  !C
  !     Hartree in Ryd
  E=EH*2.d0
  !
  DO iiij=1,JRI
     RNET(iiij)=RNOT*(exp(DH*(iiij-1.d0)))
  ENDDO
!C
  Nodes = 0
  ZZ = Z + Z
  C = 274.074D0
!
  if(.not.rel) C=1.D+10
!
  FLLP1 = FL*(FL + 1.d0)
  R83SQ = 64.D0/9.D0
  R1 = 1.D0/9.D0
  R2 = -5.D0*R1
  R3 = 19.D0*R1
  H83 = 8.D0/3.D0
!
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
  DO  2  K = 1,3
     R = RNET(K)
     DRDI = DH*R
     A(K) = (R**S)*G0
     B(K) = (R**SF)*F0
     D(1,K) = DRDI*A(K)*S/R
     D(2,K) = DRDI*B(K)*SF/R
    2 CONTINUE
!
!
      DG1 = D(1,1)
      DG2 = D(1,2)
      DG3 = D(1,3)
      DF1 = D(2,1)
      DF2 = D(2,2)
      DF3 = D(2,3)
      DO  4  K = 4, JRI
        R = RNET(K)
        DRDI = DH*R
!
!       Faktor zwei vor V wegen Hartree-Rydberg !
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
    4 CONTINUE
!
!
      DO 5001 iiij=1,JRI
         B(iiij)=B(iiij)*c/2.d0
 5001 CONTINUE
!
      VAL = A(JRI)/RNET(JRI)
      SLO = DG3/(DH*RNET(JRI))
      SLO = (SLO-VAL)/RNET(JRI) 
      RETURN
      END
