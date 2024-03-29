SUBROUTINE GTFNAM(DEFFN,ERRFN)
  USE mpi, only:  FilenameMPI, nargs, argv, stop_MPI
  IMPLICIT NONE
  CHARACTER*(*)      DEFFN, ERRFN
  !        Local Parameters
  !
  CHARACTER*5        ERREXT
  PARAMETER          (ERREXT = 'error')
  !
  !        Local Scalars
  !
  INTEGER            I
  !
  !        extract the command-line argument
  !
  !iarg=iargc()
  !modifications for pgi-mpi compiler
  if(nargs.lt.1) then
     call stop_MPI
     STOP 'GTFNAM - Exactly one commandline argument has to be given.'
  endif
  DEFFN = argv(1)
  !iarg=1
  !CALL GETARG(iarg,DEFFN)
  !        generate a name for the error-message file
  DO I = LEN(DEFFN), 1, -1
     IF (DEFFN(I:I) .EQ. '.') THEN
        IF (LEN(ERRFN) .LT. (I+LEN(ERREXT))) THEN
           call stop_MPI
           STOP 'GTFNAM - string ERRFN too short to hold filename.'
        ENDIF
        ERRFN(1:I) = DEFFN(1:I)
        ERRFN(I+1:LEN(ERRFN)) = ERREXT
        RETURN
     ENDIF
  ENDDO
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
        RETURN
     ENDIF
  ENDDO
  !        filename contains only spaces
  RETURN
END SUBROUTINE GTFNAM
