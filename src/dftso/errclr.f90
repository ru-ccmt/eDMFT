SUBROUTINE ERRCLR(FNAME)
  use mpi, ONLY: myrank, master
  IMPLICIT NONE
  CHARACTER*(*)      FNAME
  if (myrank.eq.master) then
     close(99, status='delete')
  else
     !! Nothing for now, as we do not create this files, we should not delete them
     !close(99, status='delete')
  endif
  RETURN
END SUBROUTINE ERRCLR
