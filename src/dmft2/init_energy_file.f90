SUBROUTINE INIT_ENERGY_FILE(sumw0, nkpt, numkpt, nmat, nume)
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
  k=0
  nmat=0
  nume=0
  sumw0=0
  DO ivector=1,nvector
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
           open(itape,FILE=vector_filename,STATUS='old',FORM='unformatted')
        else
           rewind(itape)
        endif
        DO i=1,nat
           READ(itape) emist
           READ(itape) emist
        enddo
        ios=0
        DO WHILE (ios == 0)
           READ(itape,IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
           !print *, 'GOT', SS, TT, ZZ, KNAME, ios
           IF (ios /= 0) CYCLE
           k=k+1
           nmat=MAX(n,nmat)
           nume=MAX(nen,nume)
           READ(itape) (wkx,wky,wkz,i=1,N)
           DO ii=1,NEn 
              READ(itape) NUM,E1
              if (abs(projector).ge.4) then
                 if (ii.eq.nemin0 .and. E1.lt.DM_Emin) DM_Emin=E1
                 if (ii.eq.nemax0 .and. E1.gt.DM_Emax) DM_Emax=E1
              endif
              if (Qcomplex) then
                 READ(itape) (wAsc,i=1,N)
              else
                 READ(itape) (wAsr,i=1,N)
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
  END DO
  nkpt=k

  CALL Reduce_MPI2(numkpt, nkpt, sumw0, DM_Emin, DM_Emax)

  if (Qprint) WRITE(6,*) 'nume=', nume, 'numkpt=', numkpt

END SUBROUTINE INIT_ENERGY_FILE
