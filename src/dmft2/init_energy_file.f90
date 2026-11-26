SUBROUTINE INIT_ENERGY_FILE(sumw0, nkpt, numkpt, nmat, nume)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use com_mpi, only: check_input_scan,mesh_failure,vectors
  USE com_mpi,only: nvector, vector_para, fvectors, Reduce_MPI2, myrank, master, VECFN, Qprint
  USE xa2, only: init_xa2
  USE com, only : nspin,nat
  USE param, only: nemin0, nemax0
  USE dmfts, only: DM_Emin, DM_Emax, projector, Qcomplex
  IMPLICIT NONE
  !
  REAL*8, intent(out) :: sumw0
  INTEGER, intent(out):: nkpt, numkpt, nmat, nume
  !
  CHARACTER*200 :: energy_filename, vector_filename
  LOGICAL       :: energy_file_exists
  CHARACTER*10  :: KNAME
  character(512) :: scan_error
  integer :: start_k
  INTEGER       :: k, ivector, itape, ios
  INTEGER       :: N, NEn, NUM, I, ii
  REAL*8        :: EMIST, SS, TT, ZZ, wgh, E1
  INTEGER       :: wkx, wky, wkz
  REAL*8        :: wAsr
  COMPLEX*16    :: wAsc
  !--------------------------------------------------------------------- 
  !find nkpt, nmat and nume in energy file
  if (abs(projector).ge.4) then
     DM_Emin=1000.
     DM_Emax=-1000.
  endif
  nspin=1                                                         
  scan_error=''
  k=0
  nmat=0
  nume=0
  sumw0=0
  DO ivector=1,nvector
     start_k=k
     if (vector_para) then
        energy_filename = fvectors(ivector, 3)
        vector_filename = fvectors(ivector, 1)
     else
        energy_filename = VECFN(3)
        vector_filename = VECFN(1)
     endif
     !INQUIRE(FILE=energy_filename, EXIST=energy_file_exists)
     energy_file_exists=.False.
     if (energy_file_exists) then
        itape=30
        open(itape,FILE=energy_filename,STATUS='old',FORM='formatted')
        DO I=1,NAT                                                  
           READ(itape,'(f9.5)') EMIST
           READ(itape,'(f9.5)') EMIST
        ENDDO
        ios=0
        DO WHILE (ios == 0)
           READ(itape,'(3e19.12,a10,2i6,F5.1)',IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
           IF (ios /= 0) CYCLE
           k=k+1
           nmat=MAX(n,nmat)
           nume=MAX(nen,nume)
           DO ii=1,NEn
              READ(itape,*) NUM,E1
              if (abs(projector).ge.4) then
                 if (ii.eq.nemin0 .and. E1.lt.DM_Emin) DM_Emin=E1
                 if (ii.eq.nemax0 .and. E1.gt.DM_Emax) DM_Emax=E1
              endif
           ENDDO
           sumw0=sumw0+wgh
        ENDDO
        close(itape)
     else
        itape=9
        if (vector_para) then
           open(itape,FILE=vector_filename,STATUS='old',FORM='unformatted',iostat=ios)
           if (ios/=0) then
              scan_error='Cannot open mapped vector '//trim(vector_filename)
              exit
           endif
        else
           rewind(itape)
        endif
        DO i=1,nat
           READ(itape,err=900,end=900) emist
           READ(itape,err=900,end=900) emist
        enddo
        ios=0
        DO WHILE (ios == 0)
           READ(itape,IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
           !print *, 'GOT', SS, TT, ZZ, KNAME, ios
           if (ios<0) exit
           if (ios/=0) goto 900
           if (N<=0.or.NEn<=0.or..not.ieee_is_finite(wgh)) goto 900
           k=k+1
           nmat=MAX(n,nmat)
           nume=MAX(nen,nume)
           READ(itape,err=900,end=900) (wkx,wky,wkz,i=1,N)
           DO ii=1,NEn 
              READ(itape,err=900,end=900) NUM,E1
              if (abs(projector).ge.4) then
                 if (ii.eq.nemin0 .and. E1.lt.DM_Emin) DM_Emin=E1
                 if (ii.eq.nemax0 .and. E1.gt.DM_Emax) DM_Emax=E1
              endif
              if (Qcomplex) then
                 READ(itape,err=900,end=900) (wAsc,i=1,N)
              else
                 READ(itape,err=900,end=900) (wAsr,i=1,N)
              endif
           ENDDO
           sumw0=sumw0+wgh
        ENDDO
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     endif
     if (vector_para) then
        if (k-start_k/=vectors(ivector,2)) then
           scan_error='Map k-point count disagrees with vector records in '//trim(vector_filename)
           exit
        endif
     endif
  END DO
  goto 910
900 scan_error='Truncated or invalid LAPW vector record in '//trim(vector_filename)
910 call check_input_scan(scan_error)
  nkpt=k

  CALL Reduce_MPI2(numkpt, nkpt, sumw0, DM_Emin, DM_Emax)

  if (Qprint) WRITE(6,*) 'nume=', nume, 'numkpt=', numkpt

  if (numkpt<=0.or..not.ieee_is_finite(sumw0).or.sumw0<=0.d0) then
     if (vector_para) then
        call mesh_failure('Empty or zero-weight k mesh: check _processes_ (dftKS) or _processes_<rank> (WIEN2k), and matching vectors.')
     else
        call mesh_failure('Empty or zero-weight k mesh in the shared vector '//trim(VECFN(1)))
     endif
  endif
END SUBROUTINE INIT_ENERGY_FILE
