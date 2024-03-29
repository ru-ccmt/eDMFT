module error
  use com_mpi, ONLY: stop_MPI, myrank, master
  IMPLICIT NONE
  INTEGER, parameter :: STDERR=99
contains
  
  SUBROUTINE OUTERR(SRNAME,ERRMSG)
    IMPLICIT NONE
    CHARACTER*(*), intent(in)  ::   SRNAME, ERRMSG
    WRITE(STDERR,'(A" - "A)') SRNAME, ERRMSG
    call stop_MPI
    STOP 'ERROR in dmft1'
    RETURN
  END SUBROUTINE OUTERR
  
  SUBROUTINE ERRFLG(FNAME,MSG)
    IMPLICIT NONE
    CHARACTER*(*), intent(in) :: FNAME, MSG
    if (myrank.eq.master) then
       OPEN (99,FILE=FNAME,ERR=900)
       WRITE (99,9000) MSG
    else
       !! Nothing for now. We do not want to have many error files at this point
       !call FilenameMPI2(FNAME)
       !OPEN (99,FILE=FNAME,ERR=900)
       !WRITE (99,9000) MSG
    endif
    RETURN
900 write(*,*)'Cannot open error-file'
    call stop_MPI
    STOP 'ERRFLG - couldn''t open errorflag-file.'
    !
9000 FORMAT (A)
  END SUBROUTINE ERRFLG

  SUBROUTINE ERRCLR
    IMPLICIT NONE
    if (myrank.eq.master) then
       close(99, status='delete')
    else
       !! Nothing for now, as we do not create this files, we should not delete them
       !close(99, status='delete')
    endif
    RETURN
  END SUBROUTINE ERRCLR

end module error
