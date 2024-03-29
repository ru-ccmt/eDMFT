      SUBROUTINE ERRCLR(FNAME)
      CHARACTER*(*)      FNAME
!
!     ..................................................................
!
! 1.     PROGRAM UNIT 'ERRCLR'
!           Clears the contents of a (error-flag) file.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Clears the contents of a file indicating an error condition,
!           remove the error flagging. ERRCLR works in conjunction with
!           the routine ERRFLG.
!
! 3.     USAGE
!           CALL ERRCLR('lapw1.error')
!
!        ARGUMENT-DESCRIPTION
!           FNAME  - CHARACTER*(*) string                        (input)
!                    The name of the file acting as error-flag.
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           ERRFLG - create an errorflag-file (notify error conditions)
!
!        INPUT/OUTPUT (READ/WRITE)
!           The contents of the file given as argument FNAME is deleted.
!           File FNAME is created if not existing and otherwise
!           overwritten.
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           none
!
! 4.     REMARKS
!           The best way to use this routine is to call ERRFLG at the
!           start of a program writing some message to the errorflag-
!           file and ERRCLR before a successful exit of the program. By
!           checking the contents of the errorflag-file it is possible
!           to determine if the program was successfully completed.
!
!           This method has the advantage of working even if some
!           runtime-error occurs which is not taken care of in the
!           program.
!
! 5.     METHOD
!           - close the errorflag-file (it was left opened by 'ERRFLG')
!           - open the errorflag-file
!           - clear the contents of the file by writing an end-of-file
!             marker at the beginning
!           - close the errorflag-file
!
! 6.     DATE
!           24. August 1993                                 Version 1.01
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
      CLOSE (99)
      OPEN (99,FILE=FNAME,STATUS='REPLACE')
!      ENDFILE (99)
!      CLOSE (99)
!
      RETURN
!
!        End of 'ERRCLR'
!
      END
