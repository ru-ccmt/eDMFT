!-*- f90 -*- - file contains F90 code in free format

SUBROUTINE occupr(charge, mu, Ek, omega, L, nbands, nomega)
  IMPLICIT NONE
  REAL*8, intent(out)    :: charge
  REAL*8, intent(in)     :: mu
  COMPLEX*16, intent(in) :: Ek(nbands,nomega)
  REAL*8, intent(in)     :: omega(nomega)
  REAL*8, intent(in)     :: L
  INTEGER, intent(in)    :: nbands, nomega
  !f2py integer intent(hide), depend(omega)  :: nomega=shape(omega,0)
  !f2py integer intent(hide), depend(Ek)     :: nbands=shape(Ek,0)
  ! local
  REAL*8  :: wdens(nbands), w1, w2, pi
  INTEGER :: ip, iom
  pi = acos(-1.0)
  !print *, 'nbands=', nbands
  !print *, 'nom=', nomega
  !print *, 'pi=', pi
  wdens=0
  ! For the lowest frequency, integral from -infinity to (om_1+om2_2)/2
  w1 = L
  w2 = 0.5*(omega(2)+omega(1))
  do ip=1,nbands
     wdens(ip) = aimag(log(w2+mu-Ek(ip,1))-log(w1+mu-Ek(ip,1)))
  enddo
  ! The rest of the frequencies
  do iom=2,nomega-1
     w1 = 0.5*(omega(iom)+omega(iom-1))
     w2 = 0.5*(omega(iom+1)+omega(iom))
     if (omega(iom+1) .gt. 0.0) w2 = 0.0
     do ip=1,nbands
        wdens(ip) = wdens(ip) + aimag(log(w2+mu-Ek(ip,iom))-log(w1+mu-Ek(ip,iom)))
     enddo
     if (omega(iom+1) .gt. 0.0) exit
  enddo
  charge = sum(wdens)*(-1/pi)
  return 
END SUBROUTINE occupr
