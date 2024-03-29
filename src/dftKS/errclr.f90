SUBROUTINE ERRCLR(FNAME)
  use mpi, ONLY: myrank, master
  IMPLICIT NONE
  CHARACTER*(*)      FNAME
  if (myrank.eq.master) then
     close(99, status='delete')
  endif
  RETURN
END SUBROUTINE ERRCLR
