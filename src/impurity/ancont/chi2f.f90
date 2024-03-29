! @Copyright 2007 Kristjan Haule
! 
subroutine chi2(chi, nrm, chi4, gweigh, vary, gwfix, fixed, sqmc, ifunr, ifuni, ders, expand, expand_sig,&
                & nvary, nfix, nal, nom, nds)
  IMPLICIT NONE
  REAL*8, intent(out):: chi, nrm, chi4
  REAL*8, intent(in) :: gweigh(nvary)
  !REAL*8, intent(in) :: vary(nvary)
  INTEGER, intent(in):: vary(nvary)
  REAL*8, intent(in) :: gwfix(nfix)
  !REAL*8, intent(in) :: fixed(nfix)
  INTEGER, intent(in) :: fixed(nfix)
  REAL*8, intent(in) :: sqmc(nom,4)
  REAL*8, intent(in) :: ifunr(nal,nom)
  REAL*8, intent(in) :: ifuni(nal,nom)
  REAL*8, intent(in) :: ders(nal,nds)
  REAL*8, intent(in) :: expand(nds)
  REAL*8, intent(in) :: expand_sig(nds)
  INTEGER, intent(in) :: nvary, nfix, nal, nom, nds
  !f2py integer intent(hide), depend(gweigh)  :: nvary=shape(gweigh,0)
  !f2py integer intent(hide), depend(gwfix)   :: nfix=shape(gwfix,0)
  !f2py integer intent(hide), depend(sqmc)    :: nom=shape(sqmc,0)
  !f2py integer intent(hide), depend(ifunr)   :: nal=shape(ifunr,0)
  !f2py integer intent(hide), depend(expand)  :: nds=shape(expand,0)
  REAL*8 :: sig0, gi, gr
  REAL*8 :: tders(nds)
  INTEGER :: im, i

  sig0 = sqrt(0.5*(sqmc(1,3)**2+sqmc(1,4)))
  chi = 0
  do im=1,nom
     gi = 0.0
     gr = 0.0
     do i=1,nvary
        gi = gi + ifuni(vary(i),im)*gweigh(i)
        gr = gr + ifunr(vary(i),im)*gweigh(i)
     enddo
     do i=1,nfix
        gi = gi + ifuni(fixed(i),im)*gwfix(i)
        gr = gr + ifunr(fixed(i),im)*gwfix(i)
     enddo
     chi = chi + ((sqmc(im,1)-gr)*sig0/sqmc(im,3))**2+((sqmc(im,2)-gi)*sig0/sqmc(im,4))**2
  enddo

  nrm = sum(gweigh)+sum(gwfix)

  tders = 0
  do i=1,nvary
     tders = tders + ders(vary(i),:)*gweigh(i)
  enddo
  do i=1,nfix
     tders = tders + ders(fixed(i),:)*gwfix(i)
  enddo
  chi4 = 0
  do i=1,nds
     chi4 = chi4 + ((tders(i)-expand(i))/expand_sig(i))**2
  enddo

end subroutine chi2
  

subroutine matsum(gr, gi, En, iom, x0, dh, wb, nom, nx)
  IMPLICIT NONE
  REAL*8, intent(out) :: gr(nom), gi(nom)
  REAL*8, intent(in)  :: En
  REAL*8, intent(in)  :: iom(nom)
  REAL*8, intent(in)  :: x0(nx)
  REAL*8, intent(in)  :: dh(nx)
  REAL*8, intent(in)  :: wb(nx)
  INTEGER, intent(in) :: nom, nx
  !f2py integer intent(hide), depend(iom)  :: nom=shape(iom,0)
  !f2py integer intent(hide), depend(x0)   :: nx=shape(x0,0)
  !REAL*8 :: simps
  REAL*8 :: omn, w1, sumr, sumi
  INTEGER :: n, k
  
  gr=0
  gi=0
  do n=1,nom
     omn = iom(n)
     sumr=0
     sumi=0
     do k=1,nx
        w1 = wb(k)/(omn**2+(x0(k)+En)**2)
        sumr = sumr + w1*(x0(k)+En)*dh(k)
        sumi = sumi + w1*omn*dh(k)
     enddo
     gr(n) = -sumr
     gi(n) = -sumi
  enddo
end subroutine matsum



subroutine KramsKron(Frc, om, F0, wb, x0, dhx, En, nom, nx)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Frc(nom)
  REAL*8, intent(in)  :: om(nom)
  REAL*8, intent(in)  :: F0(nom)
  REAL*8, intent(in)  :: wb(nx)
  REAL*8, intent(in)  :: x0(nx)
  REAL*8, intent(in)  :: dhx(nx)
  REAL*8, intent(in)  :: En
  INTEGER, intent(in) :: nom, nx
  !f2py integer intent(hide), depend(om)  :: nom=shape(om,0)
  !f2py integer intent(hide), depend(x0)  :: nx=shape(x0,0)
  REAL*8, PARAMETER :: pi = 3.14159
  INTEGER :: j, k
  REAL*8  :: omj, wj, sums, Fre, Fri
  
  do j=1,nom
     omj = om(j)
     wj = F0(j)
     sums=0
     do k=1,nx
        sums = sums + (wb(k)-wj)*dhx(k)/(omj-(x0(k)+En))
     enddo
     Fre = sums - wj*log(abs((x0(nx)+En-omj)/(omj-x0(1)-En)))
     Fri = -pi*F0(j)
     Frc(j) = dcmplx(Fre, Fri)
  enddo
end subroutine KramsKron

