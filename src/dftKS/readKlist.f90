SUBROUTINE findKlistLength(ITAPE,numkpt,ERRMSG,info)
  IMPLICIT NONE
  INTEGER,      intent(in)  :: ITAPE
  INTEGER,      intent(out) :: info, numkpt
  CHARACTER*67, intent(out) :: ERRMSG
  ! locals
  CHARACTER*10 :: KNAME
  INTEGER :: kindex, ios, n_k_p
  ERRMSG=''
  if((ITAPE.lt.1).or.(ITAPE.gt.100)) then ! is ITAPE valid unit?
     INFO = 8
     WRITE (ERRMSG,'("Invalid k-point file on unit ",i3)') ITAPE
     RETURN
  endif
  !
  n_k_p = 0
  kindex = 0
  DO
     kindex = kindex + 1
     READ (ITAPE,'(A10)',IOSTAT=ios) KNAME
     IF(ios.ne.0) THEN  ! Something wrong....
        INFO = 8
        WRITE (ERRMSG,'("Invalid k-point file on unit ",i3)') ITAPE
        RETURN
     ENDIF
     IF (KNAME .EQ. 'END       ') then  ! reached the end
        INFO = 0
        NUMKPT = kindex - 1
        REWIND(ITAPE)
        return
     ENDIF
     ! Modification to for the k.p method: Do not calculate LAPW1 for those kpoints with kname=' kp '; these states will instead be calculated in LAPWKP 
     IF(KNAME(1:4).EQ.' kp ') THEN
        kindex=kindex-1  ! these points will instead be calculated in LAPWKP
        n_k_p = n_k_p+1
     ENDIF
  END DO
END SUBROUTINE findKlistLength


SUBROUTINE readKlist(ITAPE,numkpt,KNAME,K3,WEIGHT,IPGR,E1,E2,ERRMSG,info)
  IMPLICIT NONE
  INTEGER,      intent(in)  :: ITAPE, numkpt
  CHARACTER*10, intent(out) :: KNAME(numkpt)
  REAL*8,       intent(out) :: WEIGHT(numkpt), K3(3,numkpt) !SX(numkpt), SY(numkpt), SZ(numkpt)
  REAL*8,       intent(out) :: E1, E2
  CHARACTER*3,  intent(out) :: IPGR(numkpt)
  CHARACTER*67, intent(out) :: ERRMSG
  INTEGER,      intent(out) :: info
  ! locals
  INTEGER :: kindex, ios, ISX, ISY, ISZ, IDV
  REAL*8  :: E1_, E2_
  LOGICAL :: newform
  !
  WEIGHT(:)=0.0D0
  K3(:,:)=0.d0
  
  if (numkpt.le.0) then
     INFO = 8
     WRITE (ERRMSG,'("Empty k-point file on unit ",i3)') ITAPE
     return
  endif
  
  kindex = 1
  READ (ITAPE,5101,IOSTAT=ios) KNAME(kindex), ISX, ISY, ISZ, IDV, WEIGHT(kindex), E1, E2, IPGR(kindex)
  
  IF(ios==0) THEN   ! This is the new form of case.klist
     newform=.true.
  ELSE              ! Must be the old form of case.klist, which had smaller #digits
     REWIND(ITAPE)
     READ (ITAPE,5100,IOSTAT=ios) KNAME(kindex), ISX, ISY, ISZ, IDV,WEIGHT(kindex), E1, E2, IPGR(kindex)
     if (ios.ne.0) then ! It is not even the old format. Something wrong....
        INFO = 8
        WRITE (ERRMSG,'("Invalid k-point file (2) on unit ",i3)') ITAPE
        return
     endif
     newform=.false.
  ENDIF
  IF (KNAME(kindex) .EQ. 'END       ') then  ! empty k-point file
     INFO = 8
     WRITE (ERRMSG,'("Empty k-point file on unit ",i3)') ITAPE
     return
  ENDIF

  K3(1,kindex) = DBLE(ISX)/DBLE(IDV)
  K3(2,kindex) = DBLE(ISY)/DBLE(IDV)
  K3(3,kindex) = DBLE(ISZ)/DBLE(IDV)
  IF(KNAME(kindex)(1:4).EQ.' kp ') kindex=kindex-1  ! Modification to for the k.p method
  
  DO 
     kindex = kindex + 1
     IF(kindex.GT.numkpt) THEN
        INFO=0
        return
     END IF
     !        read in K-point to be calculated and energy range
     IF(newform) THEN
        READ (ITAPE,5101) KNAME(kindex), ISX, ISY, ISZ, IDV,WEIGHT(kindex), E1_, E2_, IPGR(kindex)
     ELSE
        READ (ITAPE,5100) KNAME(kindex), ISX, ISY, ISZ, IDV,WEIGHT(kindex), E1_, E2_, IPGR(kindex)
     ENDIF
     IF (KNAME(kindex) .EQ. 'END       ') THEN
        INFO = 8
        WRITE (ERRMSG,'("Reached the end of the k-point file on second reading ",i3)') ITAPE
        return
     ENDIF
     K3(1,kindex) = DBLE(ISX)/DBLE(IDV)
     K3(2,kindex) = DBLE(ISY)/DBLE(IDV)
     K3(3,kindex) = DBLE(ISZ)/DBLE(IDV)
     IF(KNAME(kindex)(1:4).EQ.' kp ') kindex=kindex-1  ! Modification to for the k.p method
  END DO
  return
5101 FORMAT(A10,4I10,3F5.2,A3)
5100 FORMAT(A10,4I5,3F5.2,A3)
END SUBROUTINE readKlist
