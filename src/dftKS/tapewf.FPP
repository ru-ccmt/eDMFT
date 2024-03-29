#ifdef _REAL
#define Ttype REAL*8
#else
#define Ttype COMPLEX*16
#endif
SUBROUTINE TAPEWF(KNAME,K3,KZZ,Eigval,Z,HSrows,NE,WEIGHT,IPGR)
  !use matrices, only: kzz, Z, EIGVAL
  use structure, only: LATTIC, ortho
  IMPLICIT NONE
  !
  !        Scalar Arguments   
  !
  INTEGER, intent(in) :: NE, HSrows
  REAL*8,  intent(in) :: K3(3), WEIGHT
  INTEGER, intent(in) :: KZZ(3,HSrows)
  REAL*8,  intent(in) :: Eigval(NE)
  Ttype,   intent(in) :: Z(HSrows,NE)
  CHARACTER*3,  intent(in) :: IPGR
  CHARACTER*10, intent(in) :: KNAME
  !..................................................................
  !   write the wavefunction (unformatted) to a file (vector-file)
  !..................................................................
  !   Local Scalars   
  INTEGER :: I, J
  REAL*8  :: SX, SY, SZ
  !       
  !  transformation into primitiv monoclinic basis
  IF(.not.ORTHO.and.lattic(1:3).eq.'CXZ') then
     SX=K3(1)+K3(3)
     SY=K3(2)
     SZ=K3(3)-K3(1)    ! fixed CXZ bug
  else
     SX=K3(1)
     SY=K3(2)
     SZ=K3(3)
  endif
  !                                                                       
  WRITE(10) SX, SY, SZ, KNAME, HSrows, NE, WEIGHT
  WRITE(10) (KZZ(1,I),KZZ(2,I),KZZ(3,I),I=1,HSROWS)
#ifdef Extended
  WRITE(11,'(3D27.20,a10,2i6,f5.1,a3)') SX, SY, SZ, KNAME, HSROWS, NE, WEIGHT, IPGR   
#else
  WRITE(11,'(3es19.12,a10,2i6,f5.1,a3)') SX, SY, SZ, KNAME, HSROWS, NE, WEIGHT, IPGR      
#endif
  DO I = 1, NE
     WRITE(11,*) I, Eigval(I)
     WRITE(10) I, Eigval(I)
#ifdef _REAL     
     WRITE(10) (Z(J,I),J=1,HSROWS)
#else     
     WRITE(10) (CONJG(Z(J,I)),J=1,HSROWS)
#endif
  END DO
  RETURN
END SUBROUTINE TAPEWF
