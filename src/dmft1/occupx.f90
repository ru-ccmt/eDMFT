FUNCTION ferm(x)
  IMPLICIT NONE
  REAL*8 :: ferm, x
  if (abs(x).lt.100.) then
     ferm = 1/(exp(x)+1.)
     RETURN
  endif
  if (x.lt.0) then
     ferm = 1
     RETURN
  endif
  ferm = 0
  RETURN
END FUNCTION ferm

SUBROUTINE occupi(charge, dcharge1, dcharge2, mu, Ek, Ek0, wk, omega, nbands, nomega, nkp, maxbands)
  IMPLICIT NONE
  REAL*8, intent(out)    :: charge(nomega)
  REAL*8, intent(out)    :: dcharge1, dcharge2
  REAL*8, intent(in)     :: mu
  COMPLEX*16, intent(in) :: Ek(maxbands,nkp,nomega)
  REAL*8, intent(in)     :: Ek0(maxbands,nkp)
  REAL*8, intent(in)     :: wk(nkp)
  REAL*8, intent(in)     :: omega(nomega)
  INTEGER, intent(in)    :: nbands(nkp), nomega, nkp, maxbands
  !f2py integer intent(hide), depend(omega)  :: nomega=shape(omega,0)
  !f2py integer intent(hide), depend(Ek)     :: nkp=shape(Ek,1)
  !f2py integer intent(hide), depend(Ek)     :: maxbands=shape(Ek,0)
  ! local
  REAL*8  :: ferm!, cc
  REAL*8  :: pi, T, Om !, dcharge1, dcharge2
  INTEGER :: iom, ikp, i, Ntot
  complex*16  :: ci, csum, wsum
  DATA ci/(0.0D0,1.0D0)/
  pi = acos(-1.0)

  Ntot = 0.5*(omega(nomega)/omega(1)-1) + 0.4999
  !T = omega(nomega)/(2*Ntot+1)/pi

  T = omega(1)/pi 
  !print *, 'T=', T, 'n=', Ntot
  !WRITE(*,*)
  !cc=0
  do iom=1,nomega
     csum=0
     do ikp=1,nkp
        wsum=0
        do i=1,nbands(ikp)
           wsum = wsum + 1/(omega(iom)*ci+mu-Ek(i,ikp,iom)) - 1/(omega(iom)*ci+mu-Ek0(i,ikp))
        enddo
        csum = csum + wsum*wk(ikp)
     enddo
     charge(iom) = 2*T*dble(csum)
     !WRITE(*,'(A,1x,I5,1x,A,1x,I3,1x,A,1x,I3,1x,F25.15)') 'iom=', iom, 'ikp=', nkp, 'i=', 1, dble(csum)*4*T
     !cc = cc+2*T*dble(csum)
  enddo
  !WRITE(*,'(A,1x,F25.15)') 'Final results', cc*2

  dcharge1=0
  dcharge2=0
  Om = omega(nomega)+pi*T ! The bosonic frequency needed to account for the sum of the remaining fermionic frequencies
  do ikp=1,nkp
     do i=1,nbands(ikp)
        ! correction due to subtracting a constant energy 1/(iom+mu-Ek0)
        dcharge1 = dcharge1 + ferm((Ek0(i,ikp)-mu)/T)*wk(ikp)
        ! correction due to finite number of points we sum up
        dcharge2 = dcharge2 - (atan((dble(Ek(i,ikp,nomega))-mu)/Om)-atan((Ek0(i,ikp)-mu)/Om))*wk(ikp)/pi 
     enddo
  enddo
  !print *, 'dcharge1=', dcharge1
  !print *, 'dcharge2=', dcharge2
  !dcharge = dcharge1 + dcharge2
  return 
END SUBROUTINE occupi

