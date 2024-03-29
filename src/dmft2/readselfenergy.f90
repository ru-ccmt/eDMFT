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

SUBROUTINE ReadSelfenergy(fh_sig, sigma, omega, s_oo, gammac, ncorr, nom, maxsize)
  IMPLICIT NONE
  INTEGER, intent(in)     :: fh_sig, ncorr, nom, maxsize
  REAL*8, intent(out)     :: omega(nom)
  REAL*8, intent(in)      :: gammac
  COMPLEX*16, intent(out) :: sigma(maxsize,nom), s_oo(maxsize)
  ! locals
  real*8, PARAMETER       :: Ry2eV = 13.60569253d0
  REAL*8 :: sigma_vec(2*ncorr), ome
  INTEGER:: ios, iom, i
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
        IF (ABS(sigma_vec(2*i)).LT.gammac) sigma_vec(2*i) = -gammac ! Set minimum broadening for correlated orbitals
        sigma(i,iom) = dcmplx(sigma_vec(2*i-1),sigma_vec(2*i))
     ENDDO
  ENDDO
  
  READ(fh_sig,*,IOSTAT=ios) ome,(sigma_vec(i),i=1,2*ncorr) 
  IF(ios.NE.0) print *, 'ERROR Readings s_oo', fh_sig, ios
  DO i=1,ncorr
     s_oo(i) = dcmplx(sigma_vec(2*i-1),sigma_vec(2*i))
  ENDDO
  IF (ios.GT.0) print *, 'ERROR Readings self energy file', fh_sig, ios

  ! Transforming to Rydbergs
  omega(:) = omega(:)/Ry2eV
  sigma(:,:) = sigma(:,:)/Ry2eV
  s_oo(:) = s_oo(:)/Ry2eV

  
END SUBROUTINE ReadSelfenergy

SUBROUTINE FitParabola(a, b, c, chi2, xw, fw, n)
  IMPLICIT NONE
  REAL*8, intent(out) :: a, b, c, chi2
  REAL*8, intent(in)  :: xw(n), fw(n)
  INTEGER, intent(in) :: n
  ! locals
  REAL*8 :: S, Sx, Sxx, Sxxx, Sxxxx, Sy, Sxy, Sxxy, x, y, den
  INTEGER :: i
  DO i=1,n
     x = xw(i)
     y = fw(i)
     S     = S    + 1
     Sx    = Sx   + x
     Sxx   = Sxx  + x*x
     Sxxx  = Sxxx + x*x*x
     Sxxxx = Sxxxx+ x*x*x*x
     Sy    = Sy   + y
     Sxy   = Sxy  + x*y
     Sxxy  = Sxxy + x*x*y
  ENDDO
  den = Sxx**3 + S*Sxxx**2 + Sx**2*Sxxxx - Sxx*(2*Sx*Sxxx + S*Sxxxx)
  a = (Sxx**2*Sxxy - Sx*Sxxx*Sxxy + Sx*Sxxxx*Sxy + Sxxx**2*Sy -  Sxx*(Sxxx*Sxy + Sxxxx*Sy))/den
  b = (-Sx*Sxx*Sxxy + S*Sxxx*Sxxy + Sxx**2*Sxy - S*Sxxxx*Sxy - Sxx*Sxxx*Sy +  Sx*Sxxxx*Sy)/den
  c = (Sx**2*Sxxy - S*Sxx*Sxxy + S*Sxxx*Sxy + Sxx**2*Sy -  Sx*(Sxx*Sxy + Sxxx*Sy))/den
  
  chi2=0
  DO i=1,n
     chi2 = chi2 + (a+b*xw(i)+c*xw(i)*xw(i) - fw(i))**2
  ENDDO
  chi2 = chi2/n
END SUBROUTINE FitParabola

SUBROUTINE FitLine(a, b, chi2, xw, fw, n)
  IMPLICIT NONE
  REAL*8, intent(out) :: a, b, chi2
  REAL*8, intent(in)  :: xw(n), fw(n)
  INTEGER, intent(in) :: n
  ! locals
  REAL*8 :: S, Sx, Sxx, Sy, Sxy, x, y, den
  INTEGER :: i
  DO i=1,n
     x = xw(i)
     y = fw(i)
     S   = S   + 1
     Sx  = Sx  + x
     Sxx = Sxx + x*x
     Sy  = Sy  + y
     Sxy = Sxy + x*y
  ENDDO
  den = S*Sxx-Sx*Sx
  a = (Sxx*Sy-Sx*Sxy)/den
  b = (S*Sxy-Sx*Sy)/den
  chi2=0
  DO i=1,n
     chi2 = chi2 + (a+b*xw(i) - fw(i))**2
  ENDDO
  chi2 = chi2/n
END SUBROUTINE FitLine

SUBROUTINE ApproximateSelfenergy(LowE, sigma, omega, ncorr, nom, maxsize, matsubara, s_oo, WL)
  IMPLICIT NONE
  REAL*8, intent(out)    :: LowE(5,maxsize)
  INTEGER, intent(in)    :: ncorr, nom, maxsize
  REAL*8, intent(in)     :: omega(nom), WL
  COMPLEX*16, intent(in) :: sigma(maxsize,nom)
  LOGICAL, intent(in)    :: matsubara
  COMPLEX*16, intent(in) :: s_oo(maxsize)
  ! locals
  REAL*8 :: fr(4), fi(4)
  INTEGER :: n_fitr, n_fiti, i, iom, imin, iom0
  REAL*8  :: ar, br, cr, chi2r, ai, bi, chi2i, Amg, Gmg, Cmg, mins, fr0, fr1, fi0, fi1, fix
  REAL*8,PARAMETER :: Ry2eV= 13.60569253d0 
  n_fitr = 4
  n_fiti = 3
  WRITE(*,*) 'The pole approximation for the self-energy:'

  LowE=0.0
  if (matsubara) then
     DO i=1,ncorr
        DO iom=1,n_fitr
           fr(iom) = dble(1/(sigma(i,iom)-s_oo(i)))
        ENDDO
        DO iom=1,n_fiti
           fi(iom) = aimag(1/(sigma(i,iom)-s_oo(i)))
        ENDDO
        CALL FitParabola(ar, br, cr, chi2r, omega(:n_fitr), fr, n_fitr )
        CALL FitLine(ai, bi, chi2i, omega(:n_fiti), fi, n_fiti)
     
        Amg = 1/bi
        Cmg = -Amg*ar
        Gmg = Amg*ai
        ! Trying to approximate low energy self-energy by Amg/(omega-Cmg+i*Gmg)+s_oo
        if (Amg.GT.0) then
           LowE(1,i) = 1.0
           LowE(2,i) = Amg
           LowE(3,i) = Cmg
           LowE(4,i) = Gmg
           LowE(5,i) = s_oo(i)
        else
           LowE(5,i) = sigma(i,1) ! Sigma(omega=0)
        endif
        WRITE(*,'(A,F12.7,1x,A,F12.7,1x,A,F12.7,1x,A,F12.7)') 'A=', Amg*Ry2eV**2, 'G=', Gmg*Ry2eV, 'C=', Cmg*Ry2eV, 's0=', LowE(5,i)
     ENDDO
  else
     DO i=1,ncorr
        Amg=0.
        Gmg=0.
        mins=1e10
        imin=0
        DO iom=1,nom-1
           fr0 = dble(1/(sigma(i,iom)-s_oo(i)))
           fr1 = dble(1/(sigma(i,iom+1)-s_oo(i)))
           fi0 = aimag(1/(sigma(i,iom)-s_oo(i)))
           fi1 = aimag(1/(sigma(i,iom+1)-s_oo(i)))
           fix = fi0 - (fi1-fi0)/(fr1-fr0)*fr0

           !print *, iom, (sigma(i,iom)-s_oo(i))*Ry2eV, fr0/Ry2eV, fr1/Ry2eV, fi0/Ry2eV, fi1/Ry2eV, fix/Ry2eV
           !print *, 'abs(om)', abs(omega(iom)).LT.WL  , fr0*fr1 .LT. 0, abs(fix).LT.mins, 'mins=', mins/Ry2eV

           if ((abs(omega(iom)).LT.WL) .and. (fr0*fr1 .LT. 0) .and. (abs(fix).LT.mins)) then ! om inside interval [-WL,WL], real-part changes sign, and imaginary part the smallest
              imin = iom
              Cmg = omega(iom) - (omega(iom+1)-omega(iom))*fr0/(fr1-fr0)
              mins = abs(fix)
              Amg = 1./((fr1-fr0)/(omega(iom+1)-omega(iom)))
              Gmg = Amg*mins
              !print *, omega(iom), Cmg, Amg
           endif
        ENDDO
        

        if (Amg.GT.0 .and. abs(Gmg).LT.10.) then
           LowE(1,i) = 1.0
           LowE(2,i) = Amg
           LowE(3,i) = Cmg
           LowE(4,i) = Gmg
           LowE(5,i) = s_oo(i)
        else
           iom0=1
           DO iom=1,nom
              if ( abs(omega(iom)) .LT. abs(omega(iom0)) ) iom0 = iom
           ENDDO
           LowE(5,i) = sigma(i,iom0)
        endif
        WRITE(*,'(A,F12.7,1x,A,F12.7,1x,A,F12.7,1x,A,F12.7)') 'A=', Amg*Ry2eV**2, 'G=', Gmg*Ry2eV, 'C=', Cmg*Ry2eV, 's0=', LowE(5,i)*Ry2eV
     ENDDO
  endif
END SUBROUTINE ApproximateSelfenergy
