SUBROUTINE forfhs(jatom)
  use lolog, only : loor,lapw,ilo
  use atspdt, only : LQIND,LM,LQNS,LQMM,VNS1,VNS2,VNS3, &
       WNS1,WNS2,WNS3,WNS4,GFAC
  use loint, only : VNS1LO,VNS2LO,VNS3LO,WNS1LO,WNS2LO,WNS3LO, &
       WNS4LO,WNS5LO
  use structure, only : IATNR
  use param, only : LOMAX,nloat
  IMPLICIT NONE
  INTEGER, intent(in) :: JATOM
  INTEGER :: L0,LL,LLL,LP,LM1,LQX,M,MP,MM,jlo,klo
  COMPLEX*16 :: TUU,TUD,TDU,TDD,WUU,WUD,WDU,WDD
  COMPLEX*16 :: HPHASE,VPUU,VPDD,VPUD,VPDU,WPUU,WPDD,WPUD,WPDU
  COMPLEX*16, allocatable :: TUU21(:),TUU12(:),TUD21(:),TDU12(:),TUU22(:,:)
  COMPLEX*16, allocatable :: WUU21(:),WUU12(:),WUD21(:),WDU12(:),WUU22(:,:)
  COMPLEX*16, allocatable :: VPUU21(:),VPUU12(:),VPUD21(:),VPDU12(:),VPUU22(:,:)
  COMPLEX*16, allocatable :: WPUU21(:),WPUU12(:),WPUD21(:),WPDU12(:),WPUU22(:,:)

  ! case.nsh is a sequential-unformatted interchange file.  The leading
  ! V block is the historical force input.  Appended W and phase-weighted
  ! bare blocks are used by WIEN2k-24 LAPW2 stress.  Prefix reads remain
  ! backward compatible because the Fortran runtime skips the record tail.
  allocate(TUU21(nloat),TUU12(nloat),TUD21(nloat),TDU12(nloat),TUU22(nloat,nloat))
  allocate(WUU21(nloat),WUU12(nloat),WUD21(nloat),WDU12(nloat),WUU22(nloat,nloat))
  allocate(VPUU21(nloat),VPUU12(nloat),VPUD21(nloat),VPDU12(nloat),VPUU22(nloat,nloat))
  allocate(WPUU21(nloat),WPUU12(nloat),WPUD21(nloat),WPDU12(nloat),WPUU22(nloat,nloat))

  DO LQX=1,LQIND(JATOM)
     L0=LQNS(1,LQX,JATOM)-1
     LL=LQNS(2,LQX,JATOM)
     LP=LQNS(3,LQX,JATOM)-1
     if (LL.eq.0) cycle
     LM1=LQNS(4,LQX,JATOM)
     M=LQNS(5,LQX,JATOM)
     MP=LQNS(6,LQX,JATOM)
     LLL=LM(1,LM1,JATOM)
     MM=LQMM(LQX,JATOM)

     ! Phase converting a complex spherical harmonic to the real/cubic
     ! potential harmonic stored in VLM.  This lets LAPW2 recombine the bare
     ! radial matrix with a different angular factor without dividing by a
     ! possibly vanishing original Gaunt coefficient.
     HPHASE=(1.0d0,0.0d0)
     IF (IATNR(JATOM).GT.0) THEN
        IF (LL.EQ.3 .OR. LL.EQ.5 .OR. LL.EQ.7 .OR. LL.EQ.9) THEN
           IF (MM.LT.0) THEN
              HPHASE=(0.0d0,1.0d0)
           ELSE
              HPHASE=(0.0d0,-1.0d0)
           ENDIF
        ENDIF
     ELSEIF (MM.NE.0) THEN
        IF (LLL.LT.0) THEN
           HPHASE=(0.0d0,1.0d0)
           IF (MOD(MM,2).EQ.0 .AND. MM.GT.0) HPHASE=-HPHASE
        ELSE
           IF (MOD(MM,2).EQ.1 .AND. MM.GT.0) HPHASE=-HPHASE
        ENDIF
        HPHASE=HPHASE/SQRT(2.0d0)
     ENDIF

     TUU=GFAC(LQX,JATOM)*VNS1(L0+1,LM1,LP+1,JATOM)
     TDD=GFAC(LQX,JATOM)*VNS2(L0+1,LM1,LP+1,JATOM)
     TUD=GFAC(LQX,JATOM)*VNS3(L0+1,LM1,LP+1,JATOM)
     TDU=GFAC(LQX,JATOM)*VNS3(LP+1,LM1,L0+1,JATOM)
     WUU=GFAC(LQX,JATOM)*WNS1(L0+1,LM1,LP+1,JATOM)
     WDD=GFAC(LQX,JATOM)*WNS2(L0+1,LM1,LP+1,JATOM)
     WUD=GFAC(LQX,JATOM)*WNS3(L0+1,LM1,LP+1,JATOM)
     WDU=GFAC(LQX,JATOM)*WNS4(L0+1,LM1,LP+1,JATOM)
     VPUU=HPHASE*VNS1(L0+1,LM1,LP+1,JATOM)
     VPDD=HPHASE*VNS2(L0+1,LM1,LP+1,JATOM)
     VPUD=HPHASE*VNS3(L0+1,LM1,LP+1,JATOM)
     VPDU=HPHASE*VNS3(LP+1,LM1,L0+1,JATOM)
     WPUU=HPHASE*WNS1(L0+1,LM1,LP+1,JATOM)
     WPDD=HPHASE*WNS2(L0+1,LM1,LP+1,JATOM)
     WPUD=HPHASE*WNS3(L0+1,LM1,LP+1,JATOM)
     WPDU=HPHASE*WNS4(L0+1,LM1,LP+1,JATOM)
     write(71) L0,LL,LP,M,MM,MP,TUU,TDD,TUD,TDU,WUU,WDD,WUD,WDU, &
          VPUU,VPDD,VPUD,VPDU,WPUU,WPDD,WPUD,WPDU

     IF ((L0.LE.LOMAX).OR.(LP.LE.LOMAX)) THEN
        TUU21=0.d0; TUU12=0.d0; TUD21=0.d0; TDU12=0.d0; TUU22=0.d0
        WUU21=0.d0; WUU12=0.d0; WUD21=0.d0; WDU12=0.d0; WUU22=0.d0
        VPUU21=0.d0; VPUU12=0.d0; VPUD21=0.d0; VPDU12=0.d0; VPUU22=0.d0
        WPUU21=0.d0; WPUU12=0.d0; WPUD21=0.d0; WPDU12=0.d0; WPUU22=0.d0

        IF (L0.LE.LOMAX) THEN
           do jlo=1,ilo(l0,jatom)
              if (.not.lapw(l0,jatom).and.jlo.eq.1) cycle
              if (.not.loor(l0,jatom)) cycle
              TUU21(jlo)=GFAC(LQX,JATOM)*VNS1LO(L0+1,LM1,LP+1,jlo,JATOM)
              TUD21(jlo)=GFAC(LQX,JATOM)*VNS2LO(L0+1,LM1,LP+1,jlo,JATOM)
              WUU21(jlo)=GFAC(LQX,JATOM)*WNS1LO(L0+1,LM1,LP+1,jlo,JATOM)
              WUD21(jlo)=GFAC(LQX,JATOM)*WNS2LO(L0+1,LM1,LP+1,jlo,JATOM)
              VPUU21(jlo)=HPHASE*VNS1LO(L0+1,LM1,LP+1,jlo,JATOM)
              VPUD21(jlo)=HPHASE*VNS2LO(L0+1,LM1,LP+1,jlo,JATOM)
              WPUU21(jlo)=HPHASE*WNS1LO(L0+1,LM1,LP+1,jlo,JATOM)
              WPUD21(jlo)=HPHASE*WNS2LO(L0+1,LM1,LP+1,jlo,JATOM)
              IF (LP.LE.LOMAX) THEN
                 do klo=1,ilo(lp,jatom)
                    if (.not.lapw(lp,jatom).and.klo.eq.1) cycle
                    if (.not.loor(lp,jatom)) cycle
                    TUU22(jlo,klo)=GFAC(LQX,JATOM)*VNS3LO(L0+1,LM1,LP+1,jlo,klo,JATOM)
                    WUU22(jlo,klo)=GFAC(LQX,JATOM)*WNS3LO(L0+1,LM1,LP+1,jlo,klo,JATOM)
                    VPUU22(jlo,klo)=HPHASE*VNS3LO(L0+1,LM1,LP+1,jlo,klo,JATOM)
                    WPUU22(jlo,klo)=HPHASE*WNS3LO(L0+1,LM1,LP+1,jlo,klo,JATOM)
                 enddo
              ENDIF
           enddo
        ENDIF

        IF (LP.LE.LOMAX) THEN
           do klo=1,ilo(lp,jatom)
              if (.not.lapw(lp,jatom).and.klo.eq.1) cycle
              TUU12(klo)=GFAC(LQX,JATOM)*VNS1LO(LP+1,LM1,L0+1,klo,JATOM)
              TDU12(klo)=GFAC(LQX,JATOM)*VNS2LO(LP+1,LM1,L0+1,klo,JATOM)
              WUU12(klo)=GFAC(LQX,JATOM)*WNS4LO(L0+1,LM1,LP+1,klo,JATOM)
              WDU12(klo)=GFAC(LQX,JATOM)*WNS5LO(L0+1,LM1,LP+1,klo,JATOM)
              VPUU12(klo)=HPHASE*VNS1LO(LP+1,LM1,L0+1,klo,JATOM)
              VPDU12(klo)=HPHASE*VNS2LO(LP+1,LM1,L0+1,klo,JATOM)
              WPUU12(klo)=HPHASE*WNS4LO(L0+1,LM1,LP+1,klo,JATOM)
              WPDU12(klo)=HPHASE*WNS5LO(L0+1,LM1,LP+1,klo,JATOM)
           enddo
        ENDIF

        write(71) (TUU21(jlo),jlo=1,nloat),(WUU21(jlo),jlo=1,nloat), &
             (VPUU21(jlo),jlo=1,nloat),(WPUU21(jlo),jlo=1,nloat)
        write(71) (TUU12(jlo),jlo=1,nloat),(WUU12(jlo),jlo=1,nloat), &
             (VPUU12(jlo),jlo=1,nloat),(WPUU12(jlo),jlo=1,nloat)
        write(71) (TUD21(jlo),jlo=1,nloat),(WUD21(jlo),jlo=1,nloat), &
             (VPUD21(jlo),jlo=1,nloat),(WPUD21(jlo),jlo=1,nloat)
        write(71) (TDU12(jlo),jlo=1,nloat),(WDU12(jlo),jlo=1,nloat), &
             (VPDU12(jlo),jlo=1,nloat),(WPDU12(jlo),jlo=1,nloat)
        write(71) ((TUU22(jlo,klo),jlo=1,nloat),klo=1,nloat), &
             ((WUU22(jlo,klo),jlo=1,nloat),klo=1,nloat), &
             ((VPUU22(jlo,klo),jlo=1,nloat),klo=1,nloat), &
             ((WPUU22(jlo,klo),jlo=1,nloat),klo=1,nloat)
     ENDIF
  ENDDO

  deallocate(TUU21,TUU12,TUD21,TDU12,TUU22,WUU21,WUU12,WUD21,WDU12,WUU22)
  deallocate(VPUU21,VPUU12,VPUD21,VPDU12,VPUU22,WPUU21,WPUU12,WPUD21,WPDU12,WPUU22)
END SUBROUTINE forfhs
