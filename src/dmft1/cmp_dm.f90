! @Copyright 2007 Kristjan Haule
! 

! mweight_ikp = (mweight(ikp)/tweight/nsymop)
  
SUBROUTINE cmp_DensityMatrix(g_inf, g_ferm, EF, E, s_oo, STrans, DMFTU, omega, nl, ll, iSx, cix, mweight_ikp, iso, iorbital, nbands, nume, nemin, maxdim, maxdim2, ncix, nomega, csize, maxsize, natom, lmax2, norbitals)
  use param, ONLY: gamma
  IMPLICIT NONE
  complex*16, intent(inout) :: g_inf(maxdim,maxdim,ncix,nomega), g_ferm(maxdim,maxdim,ncix)
  real*8,     intent(in)    :: EF, E(nume), omega(nomega), mweight_ikp
  complex*16, intent(in)    :: s_oo(maxsize,ncix), STrans(maxsize,ncix,nbands,nbands), DMFTU(nbands,maxdim2,norbitals) 
  integer,    intent(in)    :: nl(natom), ll(natom,4), cix(natom,4), iso, iorbital(natom,lmax2+1), iSx(maxdim2, norbitals)
  integer,    intent(in)    :: maxdim, maxdim2, ncix, nomega, nume, maxsize, nbands, natom, lmax2, nemin, csize, norbitals
  ! External function
  Interface
     FUNCTION ferm(x) result( fm )
       REAL*8 :: fm, x
     end Function ferm
  end interface
  !
  real*8, allocatable     :: e_inf(:)
  complex*16, allocatable :: A_inf(:,:), cA_inf(:,:,:), tmp(:,:), tmp2(:,:), gmk(:,:,:)
  complex*16, parameter   :: imag = (0.d0, 1.d0)
  integer     :: iband, icase, jcase, l1case, icix, l1, nind1, iorb1, iom, l2case, l2, nind2, iorb2, i, ind1, ind2
  real*8      :: beta, pi
  complex*16  :: xomega
  REAL*8,PARAMETER       :: Ry2eV= 13.60569253d0

  pi=acos(-1.0D0)
  beta = pi/omega(1)
  if (beta.LT.0) then
     ! You are probably on real axis and should not select matsubara
     return
  endif

  ALLOCATE( e_inf(nbands), A_inf(nbands,nbands) )
  A_inf = 0.d0
  DO i=1,nbands
     A_inf(i,i) = E(i+nemin-1)
  ENDDO

  CALL AddSigma_optimized2(A_inf, s_oo, STrans, csize, 1, nbands, ncix, maxsize)

  CALL heigsys(A_inf, e_inf, nbands) ! hermitian eigenvalues and eigenvectors
  
  allocate( cA_inf(maxdim2,nbands,norbitals) )
  cA_inf(:,:,:) = 0.0d0
  DO icase=1,natom      
     do l1case=1,nl(icase) 
        icix = cix(icase,l1case)
        if ( icix.EQ.0 ) CYCLE
        l1 = ll(icase,l1case)
        nind1 = (2*l1+1)*iso
        iorb1 = iorbital(icase,l1case)
        ! ham = A_inf * ek * A_inf^C
        ! 1/(iom+mu-ham) = A_inf * 1/(iom+mu-ek) * A_inf^C
        ! U^+ 1/(iom+mu-ham) U = U^+ A_inf * 1/(iom+mu-ek) A (U^+ A_inf)^C = cA_inf * 1/(iom+mu-ek) * cA_inf^C
        !
        call zgemm('C','N', nind1, nbands, nbands, (1.d0,0.d0), DMFTU(:,:,iorb1), nbands, A_inf(:,:), nbands, (0.d0,0.d0), cA_inf(:,:,iorb1), maxdim2)
     enddo
  ENDDO

  allocate( tmp(maxdim2,nbands), tmp2(maxdim2,maxdim2) )
  allocate( gmk(maxdim,maxdim,ncix) )

  do iom=1,nomega
     xomega = omega(iom)*IMAG
     gmk(:,:,:)=0.d0
     DO icase=1,natom      
        do l1case=1,nl(icase) 
           icix = cix(icase,l1case)
           if ( icix.EQ.0 ) CYCLE
           l1 = ll(icase,l1case)
           nind1 = (2*l1+1)*iso
           iorb1 = iorbital(icase,l1case)
           DO jcase=1,natom
              do l2case=1,nl(jcase)
                 if ( cix(jcase,l2case).NE.icix ) CYCLE
                 l2 = ll(jcase,l2case)
                 nind2 = (2*l2+1)*iso
                 iorb2 = iorbital(jcase,l2case)
                 ! gf_inf <- tmp2(ind1,ind2) += cA_inf(ind1,iband,icix) * 1/(iom+mu-e_inf(iband)) * conjg(cA_inf(ind2,iband,icix))
                 do iband=1,nbands
                    tmp(:nind1,iband) = cA_inf(:nind1,iband,iorb1)*1/(xomega+EF-e_inf(iband)+IMAG*gamma)
                 enddo
                 call zgemm('N','C', nind1, nind2, nbands, (1.d0,0.d0), tmp, maxdim2, cA_inf(:,:,iorb2), maxdim2, (0.d0,0.d0), tmp2, maxdim2)
                 do ind1=1,nind1
                    do ind2=1,nind2
                       gmk(iSx(ind1,iorb1), iSx(ind2,iorb2), icix) = tmp2(ind1,ind2)
                    enddo
                 enddo
              enddo
           enddo
        ENDDO
     enddo
     g_inf(:,:,:,iom) = g_inf(:,:,:,iom) + gmk(:,:,:)*mweight_ikp
  ENDDO

  gmk(:,:,:)=0.d0
  DO icase=1,natom      
     do l1case=1,nl(icase) 
        icix = cix(icase,l1case)
        if ( icix.EQ.0 ) CYCLE
        l1 = ll(icase,l1case)
        nind1 = (2*l1+1)*iso
        iorb1 = iorbital(icase,l1case)
        DO jcase=1,natom
           do l2case=1,nl(jcase)
              if ( cix(jcase,l2case).NE.icix ) CYCLE
              l2 = ll(jcase,l2case)
              nind2 = (2*l2+1)*iso
              iorb2 = iorbital(jcase,l2case)
              ! gf_inf <- tmp2(ind1,ind2) += cA_inf(ind1,iband,icix) * ferm(e_inf(iband)-mu) * conjg(cA_inf(ind2,iband,icix))
              do iband=1,nbands
                 tmp(:nind1,iband) = cA_inf(:nind1,iband,iorb1)*ferm((e_inf(iband)-EF)*beta)
              enddo
              call zgemm('N','C', nind1, nind2, nbands, (1.d0,0.d0), tmp, maxdim2, cA_inf(:,:,iorb2), maxdim2, (0.d0,0.d0), tmp2, maxdim2)
              do ind1=1,nind1
                 do ind2=1,nind2
                    gmk(iSx(ind1,iorb1), iSx(ind2,iorb2), icix) = tmp2(ind1,ind2)
                 enddo
              enddo
           enddo
        enddo
     ENDDO
  enddo
  g_ferm(:,:,:) = g_ferm(:,:,:) + gmk(:,:,:)*mweight_ikp

  DEALLOCATE( A_inf, e_inf )
  deallocate( tmp, tmp2 )
  deallocate( gmk )
  deallocate( cA_inf )
END SUBROUTINE cmp_DensityMatrix

SUBROUTINE Print_DensityMatrix(gmloc, g_inf, g_ferm, omega, maxdim, ncix, Sigind, nomega, cixdim, iso)
  USE splines, ONLY: zspline3
  IMPLICIT NONE
  COMPLEX*16, intent(in) :: g_inf(maxdim,maxdim,ncix,nomega), g_ferm(maxdim,maxdim,ncix), gmloc(maxdim, maxdim, ncix, nomega)
  REAL*8,     intent(in) :: omega(nomega)
  INTEGER,    intent(in) :: maxdim, ncix, nomega, cixdim(ncix), iso, Sigind(maxdim,maxdim,ncix)
  ! Functions
  Interface
     Function NumericSecondDeriv(func,omega,ii,nomega) result(df2)
       INTEGER, intent(in)    :: nomega, ii
       COMPLEX*16, intent(in) :: func(nomega)
       REAL*8, intent(in)    :: omega(nomega)
       COMPLEX*16 :: df2
     End Function NumericSecondDeriv
  end interface
  ! locals
  COMPLEX*16, allocatable :: DM(:,:,:), ncorr(:,:)
  COMPLEX*16 :: dg(nomega)
  REAL*8, allocatable :: omega_all(:)
  COMPLEX*16, allocatable :: dg_all(:)
  INTEGER, allocatable :: cind(:), cini(:)
  INTEGER :: iom, n0_om, nom_all, cixdm, icix, ip, iq, fh_DM, cixdms, it
  REAL*8  :: pi, beta, nf, wgh
  COMPLEX*16 :: df2, df3
  real*8, PARAMETER       :: Ry2eV = 13.60569193
  !
  allocate( DM(maxdim,maxdim,ncix) )
  pi = 2*acos(0.d0)
  beta = pi/omega(1)
  if (beta.LT.0) then
     ! You are probably on real axis and should not select matsubara
     return
  endif
  do iom=1,nomega
     if (abs((omega(iom)*beta/pi+1)/2. - iom).GT.1e-1) EXIT
  enddo
  n0_om = iom-1
  nom_all = int((omega(nomega)*beta/pi+1)/2.+0.1)
  allocate(omega_all(nom_all))
  do iom=1,nom_all
     omega_all(iom) = (2*iom-1)*pi/beta
  enddo
  allocate( dg_all(nom_all) )
  
  DM(:,:,:)=0.d0
  DO icix=1,ncix
     cixdm = cixdim(icix)
     DO ip=1,cixdm
        do iq=1,cixdm
           ! Interpolating gloc on all Matsubara points
           dg(:) = gmloc(ip,iq,icix,:)-g_inf(ip,iq,icix,:)
           !print *, 'n0_om=', n0_om, 'nomega=', nomega
           dg_all(:n0_om) = dg(:n0_om)
           if (n0_om.lt.nomega) then
              df2 = NumericSecondDeriv(dg,omega,n0_om,nomega)
              df3 = NumericSecondDeriv(dg,omega,nomega-1,nomega)
              dg_all(n0_om:) = zspline3( omega(n0_om:), dg(n0_om:), omega_all(n0_om:), df2, df3) ! The tail needs to be interpolated
           endif
           DM(ip,iq,icix) = sum(dg_all)/beta
        enddo
     ENDDO
     DM(:,:,icix) = transpose(conjg(DM(:,:,icix))) + DM(:,:,icix)
     DM(:,:,icix) = DM(:,:,icix) + g_ferm(:,:,icix)
  ENDDO
  deallocate( dg_all )
  deallocate( omega_all )

  wgh = 2.d0/iso
  fh_DM=21

  do icix=1,ncix
     cixdm = cixdim(icix)
     allocate( cind(cixdm) )
     ! If Sigind(i,i)=0, we eliminate the i-th column and rown, because such orbital should be treated as non-correlated                                                                                    
     cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated                                                                                    
     cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block                                                                                         
     DO ip=1,cixdm
        it = Sigind(ip,ip,icix)
        if (it.gt.0) then
           cixdms = cixdms + 1
           cind(ip) = cixdms
        endif
     ENDDO
     allocate( cini(cixdms))
     do ip=1,cixdm
        if (cind(ip).gt.0) cini(cind(ip))=ip
     enddo

     allocate( ncorr(cixdms,cixdms) )
     nf=0     ! total nd  
     ncorr=0  ! Density matrix
     DO ip=1,cixdms
        DO iq=1,cixdms
           ncorr(ip,iq) = DM(cini(ip),cini(iq),icix)
        ENDDO
        nf = nf + dble(ncorr(ip,ip))
     ENDDO

     WRITE(fh_DM,'(A,f12.6,I3,2x,I3,2x,A)') ':NCOR', (nf*wgh), icix, cixdms, '# nf, icix, cixdm'
     DO ip=1,cixdms
        do iq=1,cixdms
           WRITE(fh_DM,'(f12.7,1x,f12.7,3x)',advance='no') dble(ncorr(ip,iq))*wgh, aimag(ncorr(ip,iq))*wgh
        end do
        WRITE(fh_DM,*)
     END DO
     WRITE(fh_DM,*)

     WRITE(6,'(A,f12.6,I3,2x,I3,2x,A)') ':NCOR', (nf*wgh), icix, cixdms, '# nf, icix, cixdm'
     DO ip=1,cixdms
        do iq=1,cixdms
           WRITE(6,'(f12.7,1x,f12.7,3x)',advance='no') dble(ncorr(ip,iq))*wgh, aimag(ncorr(ip,iq))*wgh
        end do
        WRITE(6,*)
     END DO
     WRITE(6,*)

     deallocate( ncorr )
     deallocate( cind, cini )

     DO ip=1,cixdm
        do iq=1,cixdm
           WRITE(fh_DM,'(f12.7,1x,f12.7,3x)',advance='no') dble(DM(ip,iq,icix))*wgh, aimag(DM(ip,iq,icix))*wgh
        end do
        WRITE(fh_DM,*)
     END DO
     WRITE(fh_DM,*)
     
  ENDDO
  deallocate( DM )
END SUBROUTINE Print_DensityMatrix


FUNCTION ferm(x) result( fm )
  IMPLICIT NONE
  REAL*8 :: fm, x
  if (abs(x).lt.100.) then
     fm = 1/(exp(x)+1.)
     RETURN
  endif
  if (x.lt.0) then
     fm = 1
     RETURN
  endif
  fm = 0
  RETURN
END FUNCTION ferm

Function NumericSecondDeriv(funcy,omega,ii,nomega) result(df2)
  ! Numerically determines the derivative of function func(omega) at point ii
  IMPLICIT NONE
  COMPLEX*16 :: df2
  INTEGER, intent(in)   :: nomega, ii
  COMPLEX*16, intent(in):: funcy(nomega)
  REAL*8, intent(in)    :: omega(nomega)
  ! locals
  COMPLEX*16 :: f0, f1, f2
  REAL*8     :: dx0, dx1
  if (ii.le.1 .or. ii.ge.nomega) then
     print *, 'Can not determine derivative outside the range of a function'
     print *, 'nomega=', nomega, 'and ii=', ii
  endif
  f0 = funcy(ii-1)
  f1 = funcy(ii)
  f2 = funcy(ii+1)
  dx0 = omega(ii)-omega(ii-1)
  dx1 = omega(ii+1)-omega(ii)
  df2 = ((f2-f1)/dx1 - (f1-f0)/dx0)/(0.5*(dx0+dx1))
End Function NumericSecondDeriv
