SUBROUTINE GTFNAM(DEFFN,ERRFN)
  USE mpi, only:  FilenameMPI, nargs, argv, stop_MPI, myrank, master, FilenameMPI2
  IMPLICIT NONE
  CHARACTER*(*)      DEFFN, ERRFN
  CHARACTER*5        ERREXT
  PARAMETER          (ERREXT = 'error')
  !
  !        Local Scalars
  !
  INTEGER            I
  !
  !        extract the command-line argument
  !
  if(nargs.lt.1) then
     call stop_MPI
     STOP 'GTFNAM - Exactly one commandline argument has to be given.'
  endif
  DEFFN = argv(1)
  DO I = LEN(DEFFN), 1, -1
     IF (DEFFN(I:I) .EQ. '.') THEN
        IF (LEN(ERRFN) .LT. (I+LEN(ERREXT))) THEN
           call stop_MPI
           STOP 'GTFNAM - string ERRFN too short to hold filename.'
        ENDIF
        ERRFN(1:I) = DEFFN(1:I)
        ERRFN(I+1:LEN(ERRFN)) = ERREXT
        if (myrank.ne.master) call FilenameMPI2(ERRFN)
        RETURN
     ENDIF
  ENDDO
  !
  !        the name of the definition file contains no '.', it is assumed
  !        that this name contains no extension - append the extension
  !        '.error' to get a name for the error file.
  DO I = LEN(DEFFN), 1, -1
     IF (DEFFN(I:I) .NE. ' ') THEN
        IF (LEN(ERRFN) .LT. (I+1+LEN(ERREXT))) THEN
           call stop_MPI
           STOP 'GTFNAM - string ERRFN too short to hold filename.'
        ENDIF
        ERRFN(1:I) = DEFFN(1:I)
        ERRFN(I+1:LEN(ERRFN)) = '.' // ERREXT
        if (myrank.ne.master) call FilenameMPI2(ERRFN)
        RETURN
     ENDIF
  ENDDO
  !        filename contains only spaces
END SUBROUTINE GTFNAM
