! @Copyright 2007 Kristjan Haule
! 

FUNCTION CountSelfenergy(fh_sig, ncorr)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh_sig, ncorr
  INTEGER :: CountSelfEnergy
  ! locals
  REAL*8 :: sigma_vec(2*ncorr), ome
  INTEGER:: nom, ios, i
  nom=0
  ios=0
  DO
     ! read single complex number for each correlation index
     READ(fh_sig,*,IOSTAT=ios) ome,(sigma_vec(i),i=1,2*ncorr) 
     IF (ios.NE.0) EXIT ! If the read statement was not successful, exit
     nom=nom+1
  ENDDO
  IF (ios.GT.0) print *, 'ERROR Readings self energy file', fh_sig, ios, nom
  CountSelfenergy = nom
  REWIND(fh_sig)
  RETURN
END FUNCTION CountSelfenergy

SUBROUTINE ReadSelfenergy(fh_sig, sigma, omega, s_oo, gammac, ncorr, nom, maxsize, Sigind, maxdim)
  IMPLICIT NONE
  INTEGER, intent(in)     :: fh_sig, ncorr, nom, maxsize, maxdim
  INTEGER, intent(in)     :: Sigind(maxdim,maxdim)
  REAL*8, intent(out)     :: omega(nom)
  REAL*8, intent(in)      :: gammac
  COMPLEX*16, intent(out) :: sigma(maxsize,nom), s_oo(maxsize)
  ! locals
  REAL*8 :: sigma_vec(2*ncorr), ome
  INTEGER:: ios, iom, i, j, label
  LOGICAL :: broaden_component(ncorr)
  ! Preserve the historical floor only for exclusively diagonal components.
  ! A shared diagonal/off-diagonal label cannot be broadened independently.
  broaden_component = .FALSE.
  DO i=1,maxdim
     label = ABS(Sigind(i,i))
     IF (label.GT.0.AND.label.LE.ncorr) broaden_component(label) = .TRUE.
  ENDDO
  DO j=1,maxdim
     DO i=1,maxdim
        IF (i.EQ.j) CYCLE
        label = ABS(Sigind(i,j))
        IF (label.GT.0.AND.label.LE.ncorr) broaden_component(label) = .FALSE.
     ENDDO
  ENDDO
  sigma=0
  ios=0
  DO iom=1,nom
     ! read single complex number for each correlation index
     READ(fh_sig,*,IOSTAT=ios) ome,(sigma_vec(i),i=1,2*ncorr) 
     IF(ios.NE.0) then
        print *, 'ios=', ios, 'in reading self-energy' !EXIT ! If the read statement was not successful, exit
        exit
     endif
     omega(iom) = ome
     DO i=1,ncorr
        IF (broaden_component(i).AND.ABS(sigma_vec(2*i)).LT.gammac) sigma_vec(2*i) = -gammac
        sigma(i,iom) = dcmplx(sigma_vec(2*i-1),sigma_vec(2*i))
     ENDDO
  ENDDO
  
  READ(fh_sig,*,IOSTAT=ios) ome,(sigma_vec(i),i=1,2*ncorr) 
  IF(ios.NE.0) print *, 'ERROR Readings s_oo', fh_sig, ios
  DO i=1,ncorr
     s_oo(i) = dcmplx(sigma_vec(2*i-1),sigma_vec(2*i))
  ENDDO
  IF (ios.GT.0) print *, 'ERROR Readings self energy file', fh_sig, ios
END SUBROUTINE ReadSelfenergy
