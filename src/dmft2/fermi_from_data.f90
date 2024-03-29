!------------------------------------------!
! Data to determine the chemical potential !
!------------------------------------------!
MODULE muzerox
  REAL*8                :: w_sum, w_gamma, w_norm1, w_beta
  INTEGER               :: MF, w_nomega, w_npomega
  INTEGER               :: n0_om, nkp, max_nbands
  LOGICAL               :: wprint, wmatsubara, wprint1, Qcheckbands
  REAL*8                :: elecn, DM_EF
  REAL*8,    ALLOCATABLE:: Ek(:,:), wgh(:)
  COMPLEX*16,ALLOCATABLE:: zEk(:,:,:,:)
  INTEGER,   ALLOCATABLE:: nemm(:,:)
  REAL*8,   ALLOCATABLE :: abom(:,:), w_omega(:)
  INTEGER,  ALLOCATABLE :: nomq(:), iomq(:,:), jomq(:)
  REAL*8,   ALLOCATABLE :: womq(:,:)
  REAL*8,     PARAMETER :: Ry2eV= 13.60569253d0 
  !
CONTAINS
  !
  SUBROUTINE read_muzero(filenames)
    IMPLICIT NONE
    CHARACTER*180, allocatable, intent(in) :: filenames(:)
    ! 
    INTEGER :: fh(2)
    INTEGER :: ikp, iikp, iom, j, qmax, max_nbands, ios, fi, i
    !
    REAL*8  :: elecn_, DM_EF_, w_sum, w_gamma_, w_norm1, w_beta_
    INTEGER :: nkp_, max_nbands_, qmax_, w_nomega_, n0_om_, w_npomega_
    LOGICAL :: wprint_, wmatsubara_, wprint1_, Qcheckbands_
    !
    MF = size(filenames)
    fh(1)=15
    fh(2)=16
    do i=1,MF
       open(fh(i),FILE=filenames(i),STATUS='old',FORM='unformatted',IOSTAT=ios)
       if (ios.NE.0) then
          WRITE(6,*) 'ERROR: reading ', filenames(i)
          STOP 'ERROR: Fermi'
       endif
    enddo
    
    i=1
    READ(fh(i)) elecn, DM_EF, nkp, max_nbands
    READ(fh(i)) w_nomega, n0_om, w_beta, qmax, wmatsubara
    READ(fh(i)) w_npomega, w_gamma
    READ(fh(i)) wprint,wprint1,Qcheckbands
    
    if (MF>1) then
       i=2
       READ(fh(i)) elecn_, DM_EF_, nkp_, max_nbands_
       READ(fh(i)) w_nomega_, n0_om_, w_beta_, qmax_, wmatsubara_
       READ(fh(i)) w_npomega_, w_gamma_
       READ(fh(i)) wprint_,wprint1_,Qcheckbands_
       if (elecn.ne.elecn_) print *, 'ERROR: Total charge different in up/dn fermi'
       if (nkp.ne.nkp_) then
          print *, 'ERROR: up and dn have different number of k-points'
          STOP 'fermi: # of k-points different'
       endif
       max_nbands = max(max_nbands,max_nbands_)
       if (w_npomega.ne.w_npomega_ .or. qmax.ne.qmax_) then
          print *, 'ERROR up and dn have different number of frequency points'
          STOP 'fermi: # of frequency points different'
       endif
    endif
    
    ALLOCATE( nomq(w_nomega), jomq(w_nomega), iomq(w_nomega,qmax), womq(w_nomega,qmax), w_omega(w_nomega) )
    ALLOCATE( abom(2,w_nomega) )
    
    do i=MF,1,-1
       do iom=1,w_nomega
          READ(fh(i)) w_omega(iom), abom(1,iom), abom(2,iom), nomq(iom), jomq(iom)
          do j=1,nomq(iom)
             READ(fh(i)) womq(iom,j), iomq(iom,j)
          enddo
       enddo
    enddo
    
    allocate( wgh(nkp), nemm(3,nkp) )
    allocate( zEk(max_nbands,w_nomega,nkp,MF) )
    zEk(:,:,:,:)=1000.d0
    
    do i=MF,1,-1
       do ikp=1,nkp
          READ(fh(i)) iikp, wgh(ikp), nemm(1,ikp), nemm(2,ikp), nemm(3,ikp)
          do iom=1,w_nomega
             do j = 1,nemm(3,ikp)-nemm(2,ikp)+1
                READ(fh(i)) zEk(j,iom,ikp,i)
             enddo
          enddo
       enddo
    enddo
  END SUBROUTINE read_muzero
  
  SUBROUTINE deallocate_muzero
    deallocate( zEk )
    deallocate( wgh, nemm )
    DEALLOCATE( abom )
    DEALLOCATE( nomq, jomq, iomq, womq, w_omega )
  END SUBROUTINE deallocate_muzero


  SUBROUTINE print_muzero(fh,ii)
    IMPLICIT NONE
    INTEGER, intent(in) :: fh, ii
    !
    INTEGER :: ikp, iom, j
    REAL*8  :: wkp
    WRITE(fh,*) elecn, nkp, size(zEk,1)
    WRITE(fh,*) w_nomega, n0_om, w_beta, size(iomq,2), wmatsubara
    WRITE(fh,*) w_npomega, w_gamma
    WRITE(fh,*) wprint,wprint1,Qcheckbands
    do iom=1,w_nomega
       WRITE(fh,*) w_omega(iom), abom(1,iom), abom(2,iom), nomq(iom), jomq(iom)
       do j=1,nomq(iom)
          WRITE(fh,*) '  ', womq(iom,j), iomq(iom,j)
       enddo
    enddo
    do ikp=1,nkp
       WRITE(fh,*) ikp, wgh(ikp), nemm(1,ikp), nemm(2,ikp), nemm(3,ikp)
       do iom=1,w_nomega
          do j = 1,nemm(3,ikp)-nemm(2,ikp)+1
             WRITE(fh,'(2f25.15)',advance='no') zEk(j,iom,ikp,ii)
          enddo
          WRITE(fh,*)
       enddo
    enddo
  END SUBROUTINE print_muzero
  
  REAL*8 FUNCTION Density(EF)
    IMPLICIT NONE
    REAL*8, intent(in)    :: EF
    ! functions
    REAL*8 :: ferm
    ! locals
    REAL*8    :: cc, zw0, wkp, ek0, Om, drest, ci, fer, Temperature
    INTEGER   :: ikp, num, num0, iom, nemin, DM_nemin, DM_nemax, j0, j1, j2, i, nbands, ii
    COMPLEX*16:: cE, ccn, cek, e0, e1, e2, ca, cb, omn, csum1, csum
    REAL*8     :: pi
    REAL*8     :: xx0, xx1, xx2, x2x0, x2x1
    COMPLEX*16 :: imag
    pi = acos(-1.)
    imag=dcmplx(0.d0,1.d0)

    IF (wmatsubara) THEN
       if (wprint1) WRITE(*,'(A2,1x,A2,1x,A4,1x,A4,1x,A4,1x,A14,1x,A40,1x,A20,1x,A20)') 'ik','i','iom','j1','iomq','omn', 'cek', '1./(omn-cek)', 'dsum'
       zw0=0
       !$OMP PARALLEL DO PRIVATE(nemin,DM_nemin,DM_nemax,wkp,cc,num0,ek0,Om,drest,ci,iom,j0,j1,j2,e0,e1,e2,i,omn,cek,csum,csum1,xx0,xx1,xx2,x2x0,x2x1,ii) SCHEDULE(STATIC) REDUCTION(+:zw0)
       do ikp=1,nkp
          do ii=1,MF
             nemin    = nemm(1,ikp)
             DM_nemin = nemm(2,ikp)
             DM_nemax = nemm(3,ikp)
             wkp = wgh(ikp)
             if (DM_nemin.GT.nemin) zw0 = zw0 + wkp*(DM_nemin-nemin) ! These bands are fully filled

             cc=0
             DO num=DM_nemin,DM_nemax
                num0 = num-DM_nemin+1
                ek0 = dble(zEk(num0,w_nomega,ikp,ii))-EF
                csum1=0
                csum=0
                do iom=1,n0_om-1
                   cek = zEk(num0,iom,ikp,ii)-EF
                   omn = dcmplx(0, w_omega(iom) )
                   csum1 = csum1 + 1./(omn-cek) - 1./(omn-ek0)
                   !if (wprint1) WRITE(*,'(I2,1x,I2,1x,I4,1x,I4,1x,I4,1x,f14.10,1x,2f20.15,1x,f20.15,1x,f20.15)') ikp,num,iom,iom, iom, aimag(omn), cek, dble(1./(omn-cek)-1./(omn-ek0)), dble(csum1)
                enddo

                do iom=n0_om,w_nomega-1
                   j0 = jomq(iom-1)
                   j1 = jomq(iom)
                   j2 = jomq(iom+1)
                   e0 = zEk(num0,iom-1,ikp,ii)-EF
                   e1 = zEk(num0,iom,  ikp,ii)-EF
                   e2 = zEk(num0,iom+1,ikp,ii)-EF
                   do i=1,nomq(iom)
                      omn = dcmplx(0, (2*iomq(iom,i)-1)*pi/w_beta)
                      if (.False.) then
                         if (iomq(iom,i).LT.j1) then
                            cek = e1 + (e0-e1)*(iomq(iom,i)-j1)/(j0-j1) ! linear interpolation
                         else if (iomq(iom,i).GT.j1) then
                            cek = e1 + (e2-e1)*(iomq(iom,i)-j1)/(j2-j1) ! linear interolation
                         else
                            cek = e1
                         endif
                      else
                         xx0 = iomq(iom,i)-j0
                         xx1 = iomq(iom,i)-j1
                         xx2 = iomq(iom,i)-j2
                         x2x0 = j2-j0
                         x2x1 = j2-j1
                         cek = xx2/(j1-j0)*(e0*xx1/x2x0-e1*xx0/x2x1) + e2*xx0*xx1/(x2x0*x2x1)
                      endif
                      csum = csum + womq(iom,i)*(1./(omn-cek) - 1./(omn-ek0))
                      !if (wprint1) WRITE(*,'(I2,1x,I2,1x,I4,1x,I4,1x,I4,1x,f14.10,1x,2f20.15,1x,f20.15,1x,f20.15)') ikp,num,iom,j1,iomq(iom,i),aimag(omn), cek, dble(womq(iom,i)*(1./(omn-cek) - 1/(omn-ek0))), dble(csum)
                   enddo
                enddo
                iom = w_nomega
                j0 = jomq(iom-1)
                j1 = jomq(iom)
                e0 = zEk(num0,iom-1,ikp,ii)-EF
                e1 = zEk(num0,iom,  ikp,ii)-EF
                do i=1,nomq(iom)
                   omn = dcmplx(0, (2*iomq(iom,i)-1)*pi/w_beta)
                   cek = e1 + (e0-e1)*(iomq(iom,i)-j1)/(j0-j1)
                   csum = csum + womq(iom,i)*(1./(omn-cek) - 1./(omn-ek0))
                   !if (wprint1) WRITE(*,'(I2,1x,I2,1x,I4,1x,I4,1x,I4,1x,f14.10,1x,2f20.15,1x,f20.15,1x,f20.15)') ikp,num,iom,j1,iomq(iom,i),aimag(omn), cek, dble(1./(omn-cek)), dble(csum)
                enddo

                Om = (2*jomq(w_nomega)*pi/w_beta)
                drest = (atan(ek0/Om)-atan(dble(cek)/Om))/pi
                ci = 2*dble(csum+csum1)/w_beta + drest + ferm(ek0*w_beta)

                if (Qcheckbands) then
                   if (ci.GT.1.0) ci=1.0
                   if (ci.LT.0.0) ci=0.0
                endif
                cc = cc + ci
             ENDDO
             zw0 = zw0 + wkp * cc
          enddo
       enddo
       !$OMP END PARALLEL DO
    ELSE
       Temperature = 1./w_beta
       zw0=0
       do ii=1,MF
          do ikp=1,nkp
             nemin    = nemm(1,ikp)
             DM_nemin = nemm(2,ikp)
             DM_nemax = nemm(3,ikp)
             nbands = DM_nemax-DM_nemin+1
             wkp = (wgh(ikp)/w_sum*2.0)*w_norm1
             if (DM_nemin.GT.nemin) zw0 = zw0 + wkp*(DM_nemin-nemin) ! These bands are fully filled
             do iom=1,w_npomega
                fer = ferm(w_omega(iom)/Temperature)
                ca = abom(1,iom) + w_gamma*IMAG
                cb = abom(2,iom) + w_gamma*IMAG
                DO num=1,nbands
                   cE = zEk(num,iom,ikp,ii)
                   ccn = log(cE-cb-EF)-log(cE-ca-EF)
                   zw0 = zw0 - wkp * fer * aimag(ccn)/pi
                   !if (wprint1) WRITE(991,*) ikp, iom, num, cE, ca, cb, EF, aimag(ccn)/pi!, zw0 
                ENDDO
             enddo
          enddo
       enddo
    ENDIF
    if (wprint) WRITE(*,'(A,f13.8,A,f14.8,A,f10.4,A,f14.8)') 'EF[eV]=', EF*Ry2eV, ' Density=', zw0, ' N=', elecn, ' EF[Ry]=', EF
    Density = zw0/MF-elecn
    RETURN
  END FUNCTION Density
END MODULE muzerox


REAL*8 FUNCTION ferm(x)
  IMPLICIT NONE
  REAL*8, intent(in) :: x
  if (x.LT.-100.) then
     ferm = 1.0
     RETURN
  endif
  if (x.GT.100.) then
     ferm = 0.0
     RETURN
  endif
  ferm = 1/(exp(x)+1)
  RETURN
END FUNCTION ferm

program Fermi
  use muzerox, ONLY: read_muzero, deallocate_muzero, print_muzero, Density, DM_EF, elecn, wprint, Ry2eV
  IMPLICIT NONE
  CHARACTER*180 :: filename
  REAL*8 :: dd, dd2, def, tEF, Ea, Eb, fatol, xatol, xrtol
  INTEGER :: itmax, i, MF, fi
  CHARACTER(len=180), allocatable :: argv(:), filenames(:)

  allocate( argv(iargc()) )
  DO i = 1,iargc()
     CALL getarg(i, argv(i))
     !WRITE (*,*) argv(i)
  END DO
  if (size(argv)<1) then
     print *, 'No argumens in fermi. Need filename of band energies'
     STOP 'ERROR: fermi, no filename given'
  endif
  
  MF = size(argv)
  if (MF>2) MF=2
  allocate( filenames(MF) )
  filenames(1:MF) = argv(1:MF)
  call read_muzero(filenames)
  
  deallocate(argv)
  
  !call print_muzero(6,1)
  
  itmax = 500
  fatol = 1e-10
  xatol = 1e-9
  xrtol = 1e-9

  
  tEF = DM_EF  ! Start at previous EF
  dd = Density(tEF)
  def = abs(dd)/elecn*1.
  if (def.LT.0.00001) def=0.00001
  if (def.GT.0.05)  def=0.05
  
  !print *, 'dd=', dd, 'de=', def
  dd2 = Density(tEF-def*sign(1.0d0,dd))
  !print *, 'dd2=', dd2, 'de=', def
  def = 0.5 * abs(dd)*def/abs(dd2-dd)  ! 0.5*|n-n0|*1/(dn/de)
  !print *, 'new-def=', def
  if (def.LT.0.00001) def=0.00001
  ! new
  
  if (abs(dd).gt.1e-10) then
     if (dd.lt.0) then
        do while (dd.lt.0)
           tEF = tEF + def
           dd = Density(tEF)
        end do
        Ea = tEF-def
        Eb = tEF
     else
        do while (dd.gt.0)
           tEF = tEF - def
           dd = Density(tEF)
        enddo
        Ea = tEF
        Eb = tEF+def
     endif
     wprint = .FALSE.
     CALL brent ( Density, fatol, itmax, Ea, Eb, xatol, xrtol )
     tEF = Ea
  endif

  DM_EF = tEF
  print *, 'EF set to', DM_EF*Ry2eV

  fi=10
  open(fi,FILE='EF.dat',status='unknown')
  WRITE(fi,*) DM_EF*Ry2eV
  close(fi)
  
  call deallocate_muzero
  deallocate( filenames )
end program
