      SUBROUTINE GTFNAM(DEFFN,ERRFN)
      CHARACTER*(*)      DEFFN, ERRFN
!
!     ..................................................................
!
! 1.     PROGRAM UNIT 'GTFNAM'
!           Provide names for definition- and error-files for LAPW.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Read the commandline-argument (exactly one has to be given)
!           specifying the name of the definition file ('lapw1.def')
!           and generate the name of the error file by replacing the
!           extension of the definition filename with '.error'. If no
!           extension can be found '.error' is appended.
!
! 3.     USAGE
!           CHARACTER*80 DEFFN, ERRFN
!           CALL GTFNAM(DEFFN,ERRFN)
!
!        ARGUMENT-DESCRIPTION
!           DEFFN  - CHARACTER*(*) string                       (output)
!                    on exit contains the filename of the lapw2-
!                    definition file 'lapw2.def'.
!           ERRFN  - CHARACTER*(*) string                       (output)
!                    on exit contains the filename of the file where
!                    error messages are stored (derived from the file-
!                    name of the definition-file).
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           - Subroutine GETARG for the extraction of command-line
!             arguments is used (the index for referencing command-
!             line argument starts with 0 on some machines and with
!             1 on other machines).
!           - A compiler directive to enable the use of the extension
!             subroutine GETARG is used (HP-Version)
!
! 4.     REMARKS
!           It is assumed that filename-extensions are separated by
!           character '.'.
!
! 5.     METHOD
!           - get commandline-argument (taking into account that
!             comandline-arguments are referenced starting with index
!             0 on some machines and starting with index 1 on others).
!           - 'lapw1.def' := commandline-argument
!           - look for the last occurence of character '.' in the
!             commandline-argument
!           - if found replace all characters after that '.' with
!             'error' giving the error filename
!           - otherwise append '.error' giving the error filename
!           
! 6.     DATE
!           25. August 1993                                 Version 1.01
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
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
      iarg=iargc()
      if(iarg.ne.1) goto 900
      CALL GETARG(1,DEFFN)
!      write(*,*) 'deffn',deffn,'ende'   
!      CALL GETARG(0,DEFFN)
!      IF (DEFFN .EQ. '      ') THEN
!         CALL GETARG(3,DEFFN)
!         IF (DEFFN .NE. '      ') GOTO 900
!         CALL GETARG(2,DEFFN)
!         IF (DEFFN .EQ. '      ') GOTO 900
!      ELSE
!         deffn=''
!         CALL GETARG(2,DEFFN)
!         IF (DEFFN .NE. '      ') GOTO 900
!         CALL GETARG(1,DEFFN)
!         IF (DEFFN .EQ. '      ') GOTO 900
!      ENDIF
!
!        generate a name for the error-message file
!
      DO 10 I = LEN(DEFFN), 1, -1
         IF (DEFFN(I:I) .EQ. '.') THEN
            IF (LEN(ERRFN) .LT. (I+LEN(ERREXT))) GOTO 910
            ERRFN(1:I) = DEFFN(1:I)
            ERRFN(I+1:LEN(ERRFN)) = ERREXT
            GOTO 30
         ENDIF
   10 CONTINUE
!
!        the name of the definition file contains no '.', it is assumed
!        that this name contains no extension - append the extension
!        '.error' to get a name for the error file.
!
      DO 20 I = LEN(DEFFN), 1, -1
         IF (DEFFN(I:I) .NE. ' ') THEN
            IF (LEN(ERRFN) .LT. (I+1+LEN(ERREXT))) GOTO 910
            ERRFN(1:I) = DEFFN(1:I)
            ERRFN(I+1:LEN(ERRFN)) = '.' // ERREXT
            GOTO 30
         ENDIF
   20 CONTINUE
!
!        filename contains only spaces
!
      GOTO 900
   30 CONTINUE
!
      RETURN
!
!        Errors
!
  900 STOP 'GTFNAM - Exactly one commandline argument has to be given.'
  910 STOP 'GTFNAM - string ERRFN too short to hold filename.'
!
!        End of 'GTFNAM'
!
      END
