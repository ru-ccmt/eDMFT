SUBROUTINE forfhs(jatom)
  use lolog, only : loor, lapw,ilo
  use atspdt, only  : LQIND, LM, LQNS, VNS1, VNS2, VNS3, GFAC
  use loint, only : VNS1LO, VNS2LO, VNS3LO
  use param
  IMPLICIT NONE
  INTEGER, intent(in):: JATOM
  !INCLUDE 'param.inc'
  !..................................................................
  !   FORFHS writes the radial non muffin-tin integrals uvu, uevu, ...
  !   and the Gaunt-coefficients for the forces on file 71.
  !..................................................................
  !   Local Scalars
  INTEGER            L0, Ll, LP, LM1
  INTEGER            LQX, M, MP
  COMPLEX*16         TUU, TUD, TDU, TDD
  COMPLEX*16         TUU21, TUU12, TUD21, TDU12, TUU22
  integer jlo,klo
  !
  DO LQX = 1, LQIND(JATOM)
     L0 = LQNS(1,LQX,JATOM) - 1
     LL = LQNS(2,LQX,JATOM)
     LP = LQNS(3,LQX,JATOM) - 1
     LM1 = LQNS(4,LQX,JATOM)
     M = LQNS(5,LQX,JATOM)
     MP = LQNS(6,LQX,JATOM)
     tuu=gfac(lqx,jatom)*vns1(l0+1,lm1,lp+1,jatom)
     tdd=gfac(lqx,jatom)*vns2(l0+1,lm1,lp+1,jatom)
     tud=gfac(lqx,jatom)*vns3(l0+1,lm1,lp+1,jatom)
     tdu=gfac(lqx,jatom)*vns3(lp+1,lm1,l0+1,jatom)
     write(71,'(6i3,4es19.12,/,6(3x),4es19.12)') l0,ll,lp,m,lm(2,lm1,jatom),mp,tuu,tdd,tud,tdu
     IF ((L0 .LE. LOMAX).OR.(LP .LE. LOMAX)) THEN
        tuu21=(0.0d0,0.0d0)
        tuu12=(0.0d0,0.0d0)
        tuu22=(0.0d0,0.0d0)
        tud21=(0.0d0,0.0d0)
        tdu12=(0.0d0,0.0d0)
        IF (L0 .LE. LOMAX) THEN
           do jlo=1,ilo(l0,jatom)
              if (.not.lapw(l0,jatom).and.jlo.eq.1) cycle
              IF (LOOR(L0,jatom)) THEN
                 tuu21=gfac(lqx,jatom)*vns1lo(l0+1,lm1,lp+1,jlo,jatom)
                 tud21=gfac(lqx,jatom)*vns2lo(l0+1,lm1,lp+1,jlo,jatom)
                 IF (LP .LE. LOMAX) THEN
                    do klo=1,ilo(lp,jatom)
                       if (.not.lapw(lp,jatom).and.klo.eq.1) cycle
                       IF (LOOR(LP,jatom)) THEN
                          tuu22= gfac(lqx,jatom)*vns3lo(l0+1,lm1,lp+1,jlo,klo,jatom)
                       ENDIF
                    enddo
                 ENDIF
              ENDIF
           enddo
        ENDIF
        IF (LP .LE. LOMAX) THEN
           do klo=1,ilo(lp,jatom)
              if (.not.lapw(lp,jatom).and.klo.eq.1) cycle
              tuu12=gfac(lqx,jatom)*vns1lo(lp+1,lm1,l0+1,klo,jatom)
              tdu12=gfac(lqx,jatom)*vns2lo(lp+1,lm1,l0+1,klo,jatom)       
           enddo
        ENDIF
        write(71,'(6(3x),6es19.12,/,6(3x),4es19.12)') tuu21,tuu12,tuu22,tud21,tdu12
     ENDIF
  ENDDO
  RETURN
END SUBROUTINE forfhs
