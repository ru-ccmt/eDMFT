#ifdef _REAL
#define Ttype REAL*8
#else
#define Ttype COMPLEX*16
#endif
SUBROUTINE PRTKPT(kpoint, KZZ, Eigval, Z, HSrows, NE, Qprint)
  use lapw_timer, only: READ_CPU_TIME, READ_WALL_TIME, time_pkp, START_TIMER, STOP_TIMER
  use lolog, only   : nlo
  use comc, only    : IPGR, KNAME
  use comi, only    : NBELW
  use coml, only    : SPRSWF, PRNTWF, WFTAPE
  use comr, only    : ELow, WEIGHT
  use kpts, only    : K3
  use lstapw, only  : RKMT
  use nmr
  use structure, only  : LATTIC, ortho
  use param
  use mpi, only: myrank, master, vector_para
  IMPLICIT NONE
  !        Arguments
  INTEGER, intent(in) :: kpoint, HSrows, NE
  INTEGER, intent(in) :: KZZ(3,HSrows)
  REAL*8,  intent(in) :: Eigval(NE)
  Ttype,   intent(in) :: Z(HSrows,NE)
  LOGICAL, intent(in) :: Qprint
  !        Local Scalars
  INTEGER :: j, jj
  REAL*8  :: ssx,ssy,ssz
  !
  !        External Subroutines
  !
  EXTERNAL           TAPEWF, WFPNT
  !
  !        setup timing
  !
  CALL START_TIMER(time_pkp)
  !
  !        write out time needed to calculate matrix-elements
  !
  !
  !        write out results of variational calculation
  !
  !  transformation into primitiv monoclinic basis
  IF(.not.ORTHO.and.lattic(1:3).eq.'CXZ') then
     ssx=K3(1,kpoint)+K3(3,kpoint)
     ssy=K3(2,kpoint)
     ssz=K3(3,kpoint)-K3(1,kpoint)   ! fixed CXZ bug
  else
     ssx=K3(1,kpoint)
     ssy=K3(2,kpoint)
     ssz=K3(3,kpoint)
  endif
  if (Qprint) then
     WRITE(6,6000) SSX, SSY, SSZ, KNAME(KPOINT), HSrows, WEIGHT(KPOINT), IPGR(KPOINT), (EIGVAL(J),J=1,NE)
     WRITE(6,6010) NBELW, ELow
     WRITE(6,6030)

     IF (KPOINT .EQ. 1) THEN
        WRITE(21,6001) K3(1,KPOINT), K3(2,KPOINT), K3(3,KPOINT),KNAME(KPOINT),HSrows,nlo,RKMT,WEIGHT(KPOINT),IPGR(KPOINT)
        do j=1,ne,5
           write(21,6002) j, (EIGVAL(JJ),JJ=j,min(j+4,NE))
           if(mod(j,50).eq.46)  write(21,*)
        enddo
     ENDIF
     IF (NBELW .GT. 0) THEN
        WRITE(21,6011) NBELW, ELow
     ENDIF
     IF (KPOINT .EQ. 1) THEN
        WRITE(21,6030)
     ENDIF
  endif
  
  IF (.NOT. SPRSWF) THEN
     IF (PRNTWF) THEN
        if (Qprint) Call WFPNT(KZZ, Z, HSrows,NE,NUME)
        if (myrank.eq.master) CALL WFPNT1(K3(:,kpoint),KZZ,Eigval,Z,HSrows,NE,NUME)
     ENDIF
     IF (WFTAPE) THEN
        if (vector_para .or. myrank.eq.master) CALL TAPEWF(KNAME(kpoint),K3(:,kpoint),KZZ,Eigval,Z,HSrows,NE,WEIGHT(kpoint),IPGR(kpoint))
     ENDIF
  ENDIF
  !        perfom timings
  CALL STOP_TIMER(time_pkp)
  !
  RETURN
  !
6000 FORMAT(/5X,'K=',3F10.5,3X,A10/5X,' MATRIX SIZE',I6,'  WEIGHT=',F5.2, '  PGR: ',A3,/5X,'EIGENVALUES ARE:',8(/2X,5F13.7))
6001 FORMAT(/7X,'K=',3F10.5,3X,A10/':RKM  :',1X,'MATRIX SIZE',I6,'LOs:',I4,'  RKM=',F5.2,'  WEIGHT=',F5.2, '  PGR: ',A3,/7X,'EIGENVALUES ARE:')
6002 format(':EIG',i5.5,':',3X,5F13.7)
6010 FORMAT(I13,' EIGENVALUES BELOW THE ENERGY ',F10.5)
6011 FORMAT(':WAR   : ',I13,' EIGENVALUES BELOW THE ENERGY ',F10.5)
6030 FORMAT(7X,14('****')/)
END SUBROUTINE PRTKPT
