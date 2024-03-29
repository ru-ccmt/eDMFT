      SUBROUTINE CPUTIM(DSEC)
      IMPLICIT NONE
!
!        Scalar Arguments
!
      DOUBLE PRECISION   DSEC
!
!        Locals
!
      INTEGER            ISEC
      INTEGER            DUMMY(4)
!
!        External Functions
!
      INTEGER            TIMES
      EXTERNAL           TIMES
!
!        Intrinsic Functions
!
      INTRINSIC          DBLE
!
      ISEC = TIMES(DUMMY)
      DSEC = DBLE(DUMMY(1))/100.0D0
!
      RETURN
!
!        End of 'CPUTIM'
!
      END

      SUBROUTINE WALLTIM(DSEC)
      IMPLICIT NONE
      DOUBLE PRECISION DSEC
      DSEC=0.0D0
      RETURN
      END
