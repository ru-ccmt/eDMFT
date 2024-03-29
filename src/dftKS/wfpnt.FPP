#ifdef _REAL
#define Ttype REAL*8
#else
#define Ttype COMPLEX*16
#endif
SUBROUTINE WFPNT (KZZ, Z, HSrows,NE,NUME)
  !use matrices, only: KZZ, Z
  IMPLICIT NONE
  INTEGER, intent(in) :: NE, HSrows, NUME
  INTEGER, intent(in) :: KZZ(3,HSrows)
  Ttype,   intent(in) :: Z(HSrows,NE)
  !
  !..................................................................
  !  print output editor for non-relativistic wavefunctions
  !  (cubic version)                     P.I.SORANTIN
  !..................................................................
  INTEGER            NCOLS
  PARAMETER          (NCOLS = 9)
#ifndef _REAL  
  COMPLEX*16         IMAG
  PARAMETER          (IMAG = (0.0D+0,1.0D+0))
#endif  
  INTEGER            I, ILOWER, IUPPER, J, K
  !
  IF (NE .GT. NUME) THEN
     WRITE(6,9000) NUME
     GOTO 999
  ENDIF
  !        print out all computed Eigenvectors stored in 'Z'
  !        formatted in 'NCOLS' columns
  WRITE(6,6030)
  DO ILOWER = 1, NE, NCOLS
     IUPPER = MIN(NE,ILOWER+NCOLS-1)
     WRITE(6,6000)
     WRITE(6,6040) (I,I=ILOWER,IUPPER)
     DO J=1,HSrows
        WRITE(6,6010) KZZ(1,J), KZZ(2,J), KZZ(3,J)
#ifdef _REAL           
        WRITE(6,6020) (Z(J,K),K=ILOWER,IUPPER)
#else
        !        print the real-part of the eigenvectors
        WRITE(6,6020) (DBLE(Z(J,K)),K=ILOWER,IUPPER)
        !        print the imaginary-part of the eigenvectors
        WRITE(6,6021) (DBLE(Z(J,K)*IMAG),K=ILOWER,IUPPER)
#endif
     ENDDO
  ENDDO
  !
999 RETURN
  !
6000 FORMAT(/,/,/,3X,'RECIPROCAL LATTICE VECTORS',/)
6010 FORMAT(1X,I3,2I4)
6020 FORMAT(18X,9(F9.6,2X),5X,'REALPART')
#ifndef _REAL  
6021 FORMAT(18X,9(F9.6,2X),5X,'IMAGPART')
#endif  
6030 FORMAT(1X)
6040 FORMAT(18X,9(I3,'.ENERGY '))
9000 FORMAT(' NUMBER OF ENERGIES GT ',I4,' NUME IN BASISO')
END SUBROUTINE WFPNT
