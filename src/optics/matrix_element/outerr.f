      SUBROUTINE OUTERR(SRNAME,ERRMSG)
      CHARACTER*(*)      SRNAME
      CHARACTER*(*)      ERRMSG
!
!     ..................................................................
!
! 1.     PROGRAM UNIT 'OUTERR'
!           Prints out an (error) message on the standard error device.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           OUTERR is a gemeral purpose message output routine for
!           FORTRAN 77 subroutines. A message given as an argument
!           is printed on standard error and execution of the calling
!           routine resumes.
!
! 3.     USAGE
!           SUBROUTINE CALLER
!           ...
!           EXTERNAL OUTERR
!           ...
!           CHARACTER*67 ERRMSG
!           WRITE (ERRMSG,FMT='(''error number:'',I3)') 10
!           CALL OUTERR('CALLER',ERRMSG)
!           ...
!           END
!
!        ARGUMENT-DESCRIPTION
!           - SRNAME     (input) CHARACTER*(*)
!                        The name of the routine which called 'OUTERR'.
!                        The length of 'SRNAME' should not exceed
!                        6 characters. Any characters beyond that
!                        limit will not be printed.
!           - ERRMSG     (input) CHARACTER*(*)
!                        The (error) message to be printed.
!                        The length of 'ERRMSG' should not exceed
!                        67 characters. Any characters beyond that
!                        limit will not be printed.
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
!           If OUTERR is called by
!              CALL OUTERR('FOO','couldn''t open file.')
!           the resulting (output) error message will look like:
!
!           'FOO' - couldn't open file.
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           The actual value of the standard error device (STDERR)
!           has to be setup for different brand's of computers.
!           For instance: STDERR = 7 on HP9000 (800/700) under HPUX
!                                = 0 on IBM RS6000 under AIX
!                                = 0 on Siemens/Fujitsu S100 under UXP/M
!           
! 4.     REMARKS
!           In the program LAPW1 (part of the WIEN93 software package),
!           error messages are written to a special error-file, accessed
!           as logical filenumber 99.
!
! 5.     METHOD
!           Write error message on output device 'standard error'.
!           The used format depends on the length of the input strings.
!
! 6.     DATE
!           24. August 1993                                 Version 1.21
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!
!     ..................................................................
!
      INTEGER            STDERR
      PARAMETER          (STDERR = 99)
!
      IF ((LEN(ERRMSG) .GT. 67) .AND. (LEN(SRNAME) .GT. 6)) THEN
         WRITE (STDERR,9010) SRNAME, ERRMSG
      ELSEIF (LEN(ERRMSG) .GT. 67) THEN
         WRITE (STDERR,9020) SRNAME, ERRMSG
      ELSEIF (LEN(SRNAME) .GT. 6) THEN
         WRITE (STDERR,9030) SRNAME, ERRMSG
      ELSE
         WRITE (STDERR,9040) SRNAME, ERRMSG
      ENDIF
!
      RETURN
!
 9010 FORMAT (' ''', A6, ''' - ',A67)
 9020 FORMAT (' ''', A,  ''' - ',A67)
 9030 FORMAT (' ''', A6, ''' - ',A)
 9040 FORMAT (' ''', A,  ''' - ',A)
!
!     End of 'OUTERR'
!
      END
