#ifdef _REAL
#define Ttype REAL*8
#else
#define Ttype COMPLEX*16
#endif
SUBROUTINE WFPNT1 (K3,KZZ,Eigval,Z,HSrows,NE,NUME)
  !use matrices, only: eigval, kzz, Z
  IMPLICIT NONE
  REAL*8,  intent(in) :: K3(3)
  INTEGER, intent(in) :: HSrows, NE, NUME
  INTEGER, intent(in) :: KZZ(3,HSrows)
  REAL*8,  intent(in) :: Eigval(NE)
  Ttype,   intent(in) :: Z(HSrows,NE)
!..................................................................
!   print output editor for non-relativistic wavefunctions
!   (cubic version)                     P.I.SORANTIN
!..................................................................
  INTEGER            NCOLS
  PARAMETER          (NCOLS = 1)
#ifndef _REAL  
  COMPLEX*16         IMAG
  PARAMETER          (IMAG = (0.0D+0,1.0D+0))
#endif  
  INTEGER            ILOWER, IUPPER, J
  IF (NE .GT. NUME) THEN
     WRITE(6,9000) NUME
     GOTO 999
  ENDIF
  !        print out all computed Eigenvectors stored in 'Z'
  !        formatted in 'NCOLS' columns
  WRITE(55,6030)
  write(55,6000) K3(1),K3(2),K3(3),NE,HSrows
  DO J=1,HSrows
     WRITE(55,6010) KZZ(1,J), KZZ(2,J), KZZ(3,J)
  enddo
  DO ILOWER = 1, NE, NCOLS
     IUPPER = MIN(NE,ILOWER+NCOLS-1)
     WRITE(55,6005) ilower, eigval(ilower)
     WRITE(55,6020) (Z(j,ilower),j=1,HSrows)
  ENDDO
  !
999 RETURN
  !
6000 FORMAT(/,3f10.6,2i5)
6010 FORMAT(1X,3i5)
6005 format(i6,f20.8)
6020 FORMAT(8F9.6)
!#ifndef _REAL  
!6021 FORMAT(18X,9(F9.6,2X),5X,'IMAGPART')
!#endif  
6030 FORMAT(1X)
!6040 FORMAT(18X,9(I3,'.ENERGY '))
9000 FORMAT(' NUMBER OF ENERGIES GT ',I4,' NUME IN BASISO')
END SUBROUTINE WFPNT1
