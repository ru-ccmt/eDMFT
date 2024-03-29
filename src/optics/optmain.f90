! @Copyright 2007 Kristjan Haule

SUBROUTINE Cmp_Optics(fenergy, fUdmft, fsymop, fh_m, fhb, ommax, Temperature_, delta, Nd, gammac, gamma, alphaV, Ndirection, Qsym, Qsimple, nom, nat, iso, norbitals, ncix, natom, nkpt, nmat, nume, Qcomplex, lmax2, maxdim2, maxdim, maxsize, dwindow, InterbandOnly)
  USE com_mpi, ONLY: nprocs, myrank, master, Reduce_MPI, Reduce_MPI_dos
  IMPLICIT NONE
  CHARACTER*100, intent(in) :: fenergy, fUdmft, fsymop!, fmommat
  REAL*8, intent(in)  :: ommax, delta, gammac, gamma, alphaV(3,3,Ndirection), dwindow, Temperature_
  INTEGER, intent(in) :: Nd, Ndirection, fh_m
  LOGICAL, intent(in) :: Qsym, Qsimple
  INTEGER, intent(in) :: fhb, nom, nat, iso, norbitals, ncix, natom, nkpt, nmat, nume
  LOGICAL, intent(in) :: Qcomplex
  INTEGER, intent(in) :: lmax2, maxdim2, maxdim, maxsize
  LOGICAL, intent(in) :: InterbandOnly
  ! Functions
  !COMPLEX*16 :: Opt_Inside_Simple, Opt_Inside, Opt_Noncorrelated
  interface
     REAL*8 Function ferm(x)
       REAL*8, intent(in):: x
     end function ferm
  end interface
  ! locals
  COMPLEX*16 :: ctmp(Ndirection), ctmp2(Ndirection), copt(Ndirection)
  CHARACTER*10:: skii
  LOGICAL :: Tcompute, pform
  INTEGER :: ios, pr_proc
  INTEGER :: itape, i, j, k, l, is, N, nemin, nemax, ikp, iikp, ip, iq, NE, nb_min, nb_max, dir, n_min, n_max
  REAL*8  :: EF, VOL, EMIST, renorm_wgh, wommax, wdelta, wgamma, wdwindow, Ek(nume), wgh(nkpt), womega, eps_p, eps_m, ferm_factors, Temperature
  integer :: fh_p, fh_s, fh_d, fh_o!, fh_m
  integer :: numk, nsymop, nbands, isym, ddi, iw, istart, iend, ix, iw_p, iw_m
  integer :: iorb, wkii, wkis, iom, iband, nord, N0, Nw, Nsymw, icix, tnorbitals
  integer :: nindo(norbitals), csize(ncix), tnindo(norbitals), cixdim(ncix), nl(natom), ll(natom,4)
  INTEGER :: cix(natom,4), iorbital(natom,4), nind(natom,4), iSx(maxdim2,norbitals), Sigind(maxdim,maxdim,ncix)
  INTEGER,    ALLOCATABLE :: amesh(:), bandind(:)
  REAL*8,     ALLOCATABLE :: zomega(:), conduc(:,:), opimat(:,:,:), omega(:)
  COMPLEX*16, ALLOCATABLE :: STrans(:,:,:,:), DMFTU(:,:,:), gc(:)
  COMPLEX*16, ALLOCATABLE :: sigma(:,:,:), wsigma(:,:), gij(:,:), zekw(:,:), Al(:,:,:), Ar(:,:,:)
  COMPLEX*16, ALLOCATABLE :: Cv(:,:,:), Dv(:,:,:), Vel(:,:,:), Ve(:,:,:), Cv0(:,:,:)
  COMPLEX*16, ALLOCATABLE :: L1(:,:,:), L2(:,:,:), L3(:,:,:), L4(:,:,:), Ltmp(:,:)
  COMPLEX*16, ALLOCATABLE :: lg_q(:), lg_p(:)
  LOGICAL,    ALLOCATABLE :: Qinside(:)
  ! constants
  REAL*8, PARAMETER      :: Ry2eV = 13.60569193
  REAL*8, PARAMETER      :: elresunit = 2.173985e-5  ! [ aBohr*hbar/ecoul^2 ] in units of Ohm*cm !!! Natural unit for electrical resistivity: [ aBohr*hbar/ecoul^2 ]
  REAL*8, PARAMETER      :: PI = 3.14159265358979
  COMPLEX*16, PARAMETER  :: IMAG = (0.0D0,1.0D0)
  REAL*8 :: COmega(2*Nd+1)
  
  ! Converting to Rydbergs
  wommax = ommax/Ry2eV
  wdelta = delta/Ry2eV
  wgamma = gamma/Ry2eV
  wdwindow = dwindow/Ry2eV
  Temperature = Temperature_/Ry2eV
  
  fh_p = 399 ! file for DMFT transformation UDMFT

  ! Logarithmic frequency mesh for integration
  N0 = int(log((wommax/wdelta)/Nd)/log(2.) + 1 + 0.5) ! Number of different logarithmic intervals
  if (N0 .LT. 1) N0=1 ! should not be zero
  Nw = (N0+1)*Nd/2   ! Number of all points in the logarithmic mesh  
  ALLOCATE ( amesh(2*Nd+1), zomega(2*Nw), conduc(Nw,Ndirection), gc(2*Nw) )
  gc=0
  print *, 'Number of logarithmic intervals=', N0, 'Number of all points in the logarithmic mesh=', Nw
  
  ! Many index arrays are imported from DMFT
  CALL Read_The_Rest_Basic_Arrays(fhb, nindo, cixdim, nl, ll, iorbital, cix, nind, csize, iSx, Sigind, EF, VOL, norbitals, ncix, natom, maxdim, maxdim2)
  
  ! Reading the self-energy for all orbitals
  ALLOCATE( sigma(nom,maxsize,ncix), omega(nom), wsigma(maxsize,ncix) )
  do icix=1,ncix
     itape = 80+icix
     WRITE(skii,fmt='(I2)') icix
     open(itape, file=TRIM('sig.inp')//ADJUSTL(TRIM(skii)), status='old', form='formatted')
     CALL ReadSelfenergy(itape, sigma(:,:,icix), omega, gammac, csize(icix), nom, maxsize)
  enddo
  ! Self-energy was in eV -> transform to Ry
  omega(:) = omega(:)/Ry2eV
  sigma(:,:,:) = sigma(:,:,:)/Ry2eV
  
  ! Kohn-Sham eigenvalues from the file
  open (59, file=fenergy, status='old', form='formatted')

  ! Reading symmetrization for matrix elements. It is used to transform from an irreducible k-point to general k-point.
  fh_s = 98
  open(fh_s, file=fsymop, status='old')
  READ(fh_s, fmt='(I6)') nord
  ALLOCATE( opimat(3,3,nord) )
  DO k=1,nord
     READ(fh_s, fmt='(3(3f8.5/))') ((opimat(j,i,k),i=1,3),j=1,3)
  END DO

  !! Reading linearization energies in Energy file
  itape=59
  do i=1,nat
     READ(itape, fmt='(f9.5)') EMIST !---- At the beginninge we have linearization energies --!
     READ(itape, fmt='(f9.5)') EMIST !---- We just skip them ---------------------------------!
  enddo
  
  wkis=0 ! DMFT transofmration might be in many files. Starting by the first.
  ! DMFT transformation
  WRITE(skii,fmt='(I2)') wkis

  !pform = .true.
  pform = .false.
  if (pform) then
     open(fh_p, file=TRIM(fUdmft)//ADJUSTL(TRIM(skii)), status='old',form='formatted')
  else
     open(fh_p, file=TRIM(fUdmft)//ADJUSTL(TRIM(skii)), status='old',form='unformatted')
  endif

  if (pform) then
     READ(fh_p,*) numk, nsymop, tnorbitals
     READ(fh_p,*) (tnindo(iorb), iorb=1,norbitals)
  else
     READ(fh_p) numk, nsymop, tnorbitals
     READ(fh_p) (tnindo(iorb), iorb=1,norbitals)
  endif

  ! optics matrix elements
  !open(fh_m, file=fmommat, status='old')
  READ(fh_m,*) 

  COmega=0

  
!!! ---------- Preparation of arrays for paralel executaion --------------
  pr_proc  = floor(nkpt/DBLE(nprocs)+0.999)  ! The maximum number of points calculated per processor                                          
  WRITE(6,'(A,I4,2x,A,I4)') 'pr_proc=', pr_proc, 'tot-k=', nkpt


  wgh=0
  conduc = 0
  wkii=0
  DO ikp=1,nkpt   ! over all irreducible k-points
     !--- We need to go over all k-points even though we will compute only some of them on this processor.
     !--- This is because we need to read vector file sequentially.
     iikp = ikp-myrank*pr_proc            ! The index of the point to compute. If negative, do not compute! 
     if (iikp.GT.pr_proc) EXIT            ! Processor finished. The rest of the points will be taken care of by the other processors.
     Tcompute=.FALSE.
     if (iikp.gt.0) Tcompute=.TRUE.       ! If Tcompute is true, the point needs to be computed.
     ! Reading band energies from energy file
     CALL ReadEnergiesK(Ek, wgh, NE, ikp, nkpt, nume)
     ! Start reading DMFT transformation
     CALL ReadDMFT_TransK_Outside(nsymop, nbands, nemin, numk, wkii, wkis, fUdmft, fh_p, pform, ikp)
     nemax = nemin+nbands-1
     
     ALLOCATE( bandind(NE) )
     
     DO ip=1,NE
        bandind(ip)=ip
     ENDDO
     DO ip=1,NE
        do iq=ip+1,NE
           if (abs(Ek(ip)-Ek(iq)).LT.1e-11) then
              bandind(iq)=bandind(ip)
           endif
        enddo
     ENDDO
     
     ALLOCATE( Qinside(NE) )
     Qinside=.FALSE.
     do ip=1,nbands
        Qinside(ip+nemin-1) = .TRUE.
     enddo


     ALLOCATE ( DMFTU(nbands,maxdim2,norbitals) )
     ALLOCATE ( Vel(NE,NE,3) )
     IF (Tcompute) THEN
        ALLOCATE( STrans(maxsize,ncix,nbands,nbands) )
        allocate( gij(nbands,nbands), zekw(NE,2*Nd+1), Al(nbands,nbands,2*Nd+1), Ar(nbands,nbands,2*Nd+1))
        ALLOCATE( Ve(NE,NE,3), Cv0(NE,NE,Ndirection) )
        ALLOCATE( Cv(NE,NE,Ndirection), Dv(NE,NE,Ndirection) )
        ALLOCATE( L1(NE,NE,3), L2(NE,NE,3), L3(NE,NE,3), L4(NE,NE,3), Ltmp(nbands,nbands) )
        ALLOCATE( lg_q(NE), lg_p(NE) )
     ELSE
        wgh(ikp)=0
     ENDIF
     
     CALL ReadVelocityK(Vel, NE, fh_m, nb_min, nb_max) !nemin, nbands)
          
     if (Qsym) then 
        Nsymw = nsymop  ! Over all reducible k-points (all group transformations)
     else
        Nsymw = 1       ! Only irreducible k-points needed
     endif
     
     DO isym=1, nsymop
        
        ! Continuing reading DMFT transformation
        CALL ReadDMFT_TransK_Inside(DMFTU, fh_p, pform, isym, nindo, nbands, norbitals, maxdim2)
        
        if (isym.GT.Nsymw .OR. .NOT.Tcompute) CYCLE ! If only the average optics is needed, we can go over the irreducible BZ only.

        !print *, ikp, isym
      
        ! For more efficient transformation of self-energy, we create special array of transformation STrans
        CALL CompressSigmaTransformation2(STrans, DMFTU, Sigind, iSx, cix, csize, iorbital, ll, nl, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize)
        
        ! Changing matrix elements to general k-point
        CALL Transform_To_Reducible(Ve, isym, Vel, opimat, NE, nord)
        
        !!!! For non-correlated
        L1 = Ve
        L2 = Ve
        L3 = Ve
        L4 = Ve
        DO ip=nb_min,nb_max
           DO iq=nb_min,nb_max
              ctmp=0
              DO dir=1,Ndirection
                 DO j=1,3
                    DO l=1,3
                       ctmp(dir) = ctmp(dir) + Ve(ip,iq,j)*Ve(iq,ip,l)*alphaV(l,j,dir)
                    ENDDO
                 ENDDO
              ENDDO
              Cv0(ip,iq,:) = ctmp
              Cv(ip,iq,:) = ctmp
              Dv(ip,iq,:) = ctmp
           ENDDO
        ENDDO
        ! For non-correlated
        do i=1,2*Nd+1
           zekw(:,i) = Ek(:)
        enddo

        iw=1
        ddi = 1
        do j=1,N0
           ! The algoritm for the frequency mesh is the following:
           ! We start with Nd negative and Nd positive points of equidistant mesh. Incluing zero, we have 2*Nd+1 points.
           ! The points are [-Nd,-Nd+1,...0, 1,.., Nd]*wdelta
           ! With this mesh, we can compute all integrals up to frequency Nd*wdelta, where wdelta is the distance between the points
           ! In the second step, we keep every second point: [-Nd,-Nd+2,...0,2,...Nd]*wdelta
           ! and we add additional points, so that we have again 2*Nd+1 points. We thus have the following points
           ! [-2*Nd, -2*Nd+2,...0,2,..2*Nd]*wdelta.
           ! The old points (one half of all the points) do not need to be recomputed. Only the new points are recomputed.
           ! The interval is doubles N0 times. N0 depends on the cutoff ommax.
           if (j.EQ.1) then
              ix=1
              do i=-Nd,Nd
                 amesh(i+Nd+1) = ddi*i
                 CALL FindEigsys(zekw(nemin:nemax,i+Nd+1), Al(:,:,i+Nd+1), Ar(:,:,i+Nd+1), ix, wsigma, gij, amesh(i+Nd+1)*wdelta, omega, STrans, Ek, nom, ncix, nbands, nemin, maxsize, nume, csize, sigma)
              enddo
           else
              do i=1,Nd/2 ! These points were computed in previous step. Just reshuffle them into the right place.
                 amesh( i+Nd+1) = amesh( 2*i+Nd+1)
                 zekw(nemin:nemax,i+Nd+1) = zekw(nemin:nemax,2*i+Nd+1)
                 Al(:,:,i+Nd+1) = Al(:,:,2*i+Nd+1)
                 Ar(:,:,i+Nd+1) = Ar(:,:,2*i+Nd+1)
              enddo
              do i=-1,-Nd/2,-1 ! These points were computed in previous step. Just reshuffle them into the right place.
                 amesh( i+Nd+1) = amesh( 2*i+Nd+1)
                 zekw(nemin:nemax,i+Nd+1) = zekw(nemin:nemax,2*i+Nd+1)
                 Al(:,:,i+Nd+1) = Al(:,:,2*i+Nd+1)
                 Ar(:,:,i+Nd+1) = Ar(:,:,2*i+Nd+1)
              enddo
              do i=-Nd,-Nd/2-1 ! These points are new. Need to compute eigensystem for them.
                 amesh(i+Nd+1) = ddi*i
                 CALL FindEigsys(zekw(nemin:nemax,i+Nd+1), Al(:,:,i+Nd+1), Ar(:,:,i+Nd+1), ix, wsigma, gij, amesh(i+Nd+1)*wdelta, omega, STrans, Ek, nom, ncix, nbands, nemin, maxsize, nume, csize, sigma)
              enddo
              do i=Nd/2+1,Nd,1 ! These points are new. Need to compute eigensystem for them.
                 amesh(i+Nd+1) = ddi*i
                 CALL FindEigsys(zekw(nemin:nemax,i+Nd+1), Al(:,:,i+Nd+1), Ar(:,:,i+Nd+1), ix, wsigma, gij, amesh(i+Nd+1)*wdelta, omega, STrans, Ek, nom, ncix, nbands, nemin, maxsize, nume, csize, sigma)
              enddo
           endif
        
           if (j.EQ.1) then
              istart=1
              iend=Nd
           else 
              istart=Nd/2+1
              iend=Nd
           endif
           if (Temperature.NE.0.0) then
              istart=1
              iend=Nd-4*Temperature/wdelta
              if (iend<1) iend=1
           endif
           WRITE(*,'(A,I4,1x,A,I3,1x,A,I3,1x,A,I3,1x,A,I4,1x,A,I4,1x,A,F10.3)') 'ikp=', ikp, 'isym=', isym, 'j=', j, 'istart=', istart, 'iend=', iend, 'Nd=', Nd, '4*T/d=', 4*Temperature/wdelta

           ix=1
           do k=istart,iend
              ! Frequency mesh for printing
              zomega(Nw+iw)   = amesh(k+Nd+1)*wdelta ! positive frequencies
              zomega(Nw-iw+1) = amesh(Nd-k+1)*wdelta ! negative frequencies
              ! For Density of state
              do ip=1,nbands
                 iband = ip+nemin-1
                 gc(Nw+iw)   = gc(Nw+iw)   + wgh(ikp)/(zomega(Nw+iw)  +EF-zekw(iband,Nd+k+1)+wgamma*IMAG)
                 gc(Nw-iw+1) = gc(Nw-iw+1) + wgh(ikp)/(zomega(Nw-iw+1)+EF-zekw(iband,Nd-k+1)+wgamma*IMAG)
              enddo                            
              
              womega = amesh(k+Nd+1)*wdelta
              
              n_min = nb_min
              n_max = nb_max
              do ip=nb_min,nb_max
                 if ( (Ek(ip)-EF).LT.(-womega-wdwindow) ) n_min = ip
              enddo
              do ip=nb_max,nb_min,-1
                 if ( (Ek(ip)-EF).GT.(womega+wdwindow) ) n_max = ip
              enddo

              ctmp=0  ! integration over epsilon=[0,omega]
              if (Temperature.NE.0.0) then !! FINITE TEMPERATURE
                 ! !$OMP PARALLEL DO SHARED(COmega,amesh,Al,Ar,Ve,alphaV) PRIVATE(iw_p,iw_m,eps_p,eps_m,ferm_factors,Cv,Dv,L1,L2,L3,L4,Ltmp,copt) SCHEDULE(STATIC) &
                 ! !$OMP& REDUCTION(+:ctmp)
                 do l=k+1-Nd,Nd+1
                    iw_p = Nd+l                 ! eps
                    iw_m = Nd-k+l               ! eps-omega
                 
                    eps_p = amesh(iw_p)*wdelta    ! eps
                    eps_m = amesh(iw_m)*wdelta    ! eps-omega
                    ferm_factors = ferm(eps_m/Temperature)-ferm(eps_p/Temperature)
                 
                    CALL cmp_C_and_D(Cv, Dv, Al(:,:,iw_p), Ar(:,:,iw_p), Al(:,:,iw_m), Ar(:,:,iw_m), Ve, alphaV, Ndirection, L1, L2, L3, L4, Ltmp, Qinside, nbands, NE, nemin, nemax, n_min, n_max)
                    
                    if (Qsimple) then
                       CALL Opt_Inside_Simple(copt, Cv, Dv, zekw, k, -10000, amesh, iw_p, iw_m, wdelta, wgamma, EF, nbands, Nd, bandind, Qinside, nume, Ndirection, NE, nemin, n_min, n_max, InterbandOnly)
                       ctmp = ctmp + copt*ferm_factors
                    else
                       CALL Opt_Inside(copt, Cv, Dv, zekw, k, -10000, amesh, iw_p, iw_m, wdelta, wgamma, EF, lg_q, lg_p, nbands, Nd, bandind, Qinside, nume, Ndirection, NE, nemin, n_min, n_max, InterbandOnly)
                       ctmp = ctmp + copt*ferm_factors
                    endif

                    !if (k.eq.1) WRITE(*,'(A,I3,1x,A,F12.6,1x,A,I3,1x,A,F12.6,1x,A,F12.6,1x,A,F12.6)') 'iw_p=', iw_p, 'eps_p=', eps_p/Temperature, 'iw_m=', iw_m, 'eps_m=', eps_m/Temperature, 'ferm=', ferm_factors, 'opt=', -dble(copt(1))
                    if (k.eq.1) then
                       COmega(iw_p) = COmega(iw_p) - dreal(copt(1))*wgh(ikp)*(4*pi)/(2*pi**2*elresunit*VOL)
                    endif
                 enddo
                 ! !$OMP END PARALLEL DO
              else
                 do l=1,k+1
                    iw_p = Nd+l                 ! eps
                    iw_m = Nd-k+l               ! eps-omega
                 
                    CALL cmp_C_and_D(Cv, Dv, Al(:,:,iw_p), Ar(:,:,iw_p), Al(:,:,iw_m), Ar(:,:,iw_m), Ve, alphaV, Ndirection, L1, L2, L3, L4, Ltmp, Qinside, nbands, NE, nemin, nemax, n_min, n_max)
                    
                    if (Qsimple) then
                       CALL Opt_Inside_Simple(copt, Cv, Dv, zekw, k, l, amesh, iw_p, iw_m, wdelta, wgamma, EF, nbands, Nd, bandind, Qinside, nume, Ndirection, NE, nemin, n_min, n_max, InterbandOnly)
                       ctmp = ctmp + copt
                    else
                       CALL Opt_Inside(copt, Cv, Dv, zekw, k, l, amesh, iw_p, iw_m, wdelta, wgamma, EF, lg_q, lg_p, nbands, Nd, bandind, Qinside, nume, Ndirection, NE, nemin, n_min, n_max, InterbandOnly)
                       ctmp = ctmp + copt
                    endif
                 enddo
              endif
              
              !!! For non-correlated
              !womega = amesh(k+Nd+1)*wdelta ! omega
              CALL Opt_Noncorrelated(copt,Cv0,Ek,womega,EF,wgamma,nb_min,nb_max,nemin,nemax,nume,Qinside,NE,Ndirection,InterbandOnly)
              ctmp2 = copt
              
              DO dir=1,Ndirection
                 ! We add proper dimension to optics
                 conduc(iw,dir) = conduc(iw,dir) - dreal(ctmp(dir) + ctmp2(dir))*wgh(ikp)*(4*pi)/(2*pi**2*elresunit*VOL)
              ENDDO
              iw = iw + 1
           enddo        
           ddi = ddi * 2
        enddo 
     ENDDO
     DEALLOCATE( bandind )
     DEALLOCATE( Qinside )
     DEALLOCATE( DMFTU, Vel )
     IF (Tcompute) THEN
        DEALLOCATE( lg_q, lg_p )
        DEALLOCATE( L1, L2, L3, L4, Ltmp )
        DEALLOCATE( Cv, Dv, Ve, Cv0 )
        DEALLOCATE( gij, zekw, Al, Ar)
        DEALLOCATE( STrans )
     ENDIF
  ENDDO

  ! k-points weights are not normalized. Need to normalize at the end!
  renorm_wgh=0
  DO ikp=1,nkpt
     renorm_wgh = renorm_wgh + wgh(ikp)
  ENDDO

  CALL Reduce_MPI(conduc, COmega, Nw, Ndirection, renorm_wgh, Nd)
  CALL Reduce_MPI_dos(gc,Nw)
  
  IF (myrank.EQ.master) THEN
     
     ! Normalization over irreducible k-points and over symmetrization (all k-points)
     gc = gc/(renorm_wgh*Nsymw)
     conduc = conduc/(renorm_wgh*Nsymw) * (2.0/iso)  ! changed in 2019. Needs to be multiplied by 2 due to spin, except if spin-orbit is included
     
     ! Output of total DOS and optics
     fh_d = 501
     open(fh_d, file='optdos.dat', status='unknown', form='formatted')
     WRITE(fh_d,*) '# omega[eV]   DOS[1/eV]'
     DO iom=1,2*Nw
        if (zomega(iom).NE.0) WRITE(fh_d,*) zomega(iom)*Ry2eV, -aimag(gc(iom))/PI/Ry2eV
     ENDDO
     close(fh_d)
     
     fh_o = 502
     open(fh_o, file='optics.dat', status='unknown', form='formatted')
     WRITE(fh_o,*) '# omega[eV]   sigma[1/(Ohm*cm)]'
     DO iom=1,Nw
        if (zomega(Nw+iom).GT.0.0) then
           WRITE(fh_o, '(f14.10,2x)', advance='no') zomega(Nw+iom)*Ry2eV
           DO dir=1,Ndirection
              WRITE(fh_o, '(f20.5,1x)', advance='no') conduc(iom,dir)/(zomega(Nw+iom))
           ENDDO
           WRITE(fh_o,*)
        endif
     ENDDO

     if (Temperature.NE.0.0) then
        COmega = COmega/(renorm_wgh*Nsymw)
        fh_d = 503
        open(fh_d, file='what_we_integrate.dat', status='unknown', form='formatted')
        WRITE(fh_d,*) '# omega[eV]   OMEGA[1/eV]'
        DO iom=2,2*Nd+1
           WRITE(fh_d,*) (iom-Nd-1)*wdelta*Ry2eV, COmega(iom)/wdelta
        ENDDO
        close(fh_d)
     endif
     
  ENDIF
  
998 CONTINUE

  close(fh_s)
  close(fh_p)
  close(59)
  close(60)
  DEALLOCATE( amesh, zomega, conduc, gc )
  DEALLOCATE( opimat )
  DEALLOCATE( sigma, omega, wsigma )
9010 FORMAT(/,2X,' KP:',I6,' NEMIN NEMAX : ',2I5, ' dE:',2f5.2,' K:',a10 /)
9040 FORMAT(3X,2I4,6E13.6,F13.8)
END SUBROUTINE Cmp_Optics



SUBROUTINE ReadVelocityK(Vel, NE, fh_m, nb_min, nb_max)!, nemin, nbands)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Vel(NE,NE,3)
  INTEGER, intent(out) :: nb_min, nb_max
  INTEGER, intent(in) :: fh_m, NE!, nemin, nbands
  ! locals
  INTEGER      :: iikp, i1, i2, ib1, ib2, ixyz
  REAL*8       :: wEmin, wEmax, dE
  CHARACTER*10 :: kname
  COMPLEX*16   :: O(3)
  
  Vel = 0
!!! Reading optics matrix elements
  READ(fh_m,9010) iikp,nb_min,nb_max,wEmin,wEmax,kname
  do i1=nb_min,nb_max
     do i2=i1,nb_max
        READ(fh_m, 9040) ib1,ib2, O(1), O(2), O(3), dE
        if (i1.ne.ib1 .or. i2.ne.ib2) then
           print *, 'ERROR', i1, ib1, i2, ib2
        endif
        !if (i1.GE.nemin .AND. i1.LT.nemin+nbands .AND. i2.GE.nemin .AND. i2.LT.nemin+nbands) then
        if (i1.LE.NE .and. i2.LE.NE) then
           DO ixyz=1,3
              Vel(i1,i2,ixyz) = O(ixyz)
              if (i2.NE.i1) then
                 Vel(i2,i1,ixyz) = conjg(O(ixyz))
              endif
           ENDDO
        endif
     enddo
  enddo
  if (nb_max.GT.NE) nb_max = NE
!9010 FORMAT(/,2X,' KP:',I6,' NEMIN NEMAX : ',2I5, ' dE:',2f5.2,' K:',a10 /)
9010 FORMAT(/,2X,4X,I6,15X,2I5,4X,2f5.2,3X,a10 /)
9040 FORMAT(3X,2I4,6E13.6,F13.8)
END SUBROUTINE ReadVelocityK

SUBROUTINE ReadEnergiesK(Ek, wgh, NE, ikp, nkpt, nume)
  IMPLICIT NONE
!!!! Reading band energies from energy file
  REAL*8, intent(out) :: wgh(nkpt), Ek(nume)
  INTEGER, intent(out):: NE
  INTEGER, intent(in) :: ikp, nkpt, nume
  ! locals
  INTEGER :: is, ios, itape, ii, NUM
  REAL*8  :: S, T, Z, WG, E1
  CHARACTER*10 :: KNAME
  INTEGER :: N
  itape=59
  READ(itape,'(3e19.12,a10,2i6,e19.12)',IOSTAT=ios) S,T,Z,KNAME,N,NE,WG
  wgh(ikp) = WG
  DO ii=1,NE
     READ(itape,*) NUM, E1
     Ek(NUM)=E1
  ENDDO
END SUBROUTINE ReadEnergiesK
  

SUBROUTINE ReadDMFT_TransK_Outside(nsymop, nbands, nemin, numk, wkii, wkis, fUdmft, fh_p, pform, ikp)  
!!! Start reading DMFT transformation
  IMPLICIT NONE
  INTEGER, intent(out)   :: nsymop, nbands, nemin
  INTEGER, intent(inout) :: numk, wkii, wkis
  CHARACTER*100, intent(in) :: fUdmft
  INTEGER, intent(in)    :: fh_p, ikp
  LOGICAL, intent(in)    :: pform
  ! locals
  CHARACTER*10 :: skii
  INTEGER :: tnorbitals, iikp, tmaxdim2
  if (wkii==numk) then ! There are many udmfile's because of parallel execution! Need to read one after the other.
     wkii=0
     wkis = wkis + 1
     WRITE(skii,fmt='(I2)') wkis
     close(fh_p)
     if (pform) then
        open(fh_p, file=TRIM(fUdmft)//ADJUSTL(TRIM(skii)), status='old', form='formatted')
        READ(fh_p,*) numk, nsymop, tnorbitals
     else
        open(fh_p, file=TRIM(fUdmft)//ADJUSTL(TRIM(skii)), status='old', form='unformatted')
        READ(fh_p) numk, nsymop, tnorbitals
     endif
  endif

  if (pform) then
     READ(fh_p,*) iikp, nbands, tmaxdim2, tnorbitals, nemin
  else
     READ(fh_p) iikp, nbands, tmaxdim2, tnorbitals, nemin
  endif
  
  if (iikp.NE.ikp) print *, 'ERROR: ikp and iikp=', ikp, iikp
  
  wkii = wkii + 1
END SUBROUTINE ReadDMFT_TransK_Outside

SUBROUTINE ReadDMFT_TransK_Inside(DMFTU, fh_p, pform, isym, nindo, nbands, norbitals, maxdim2)
  ! Continuing reading DMFT transformation
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: DMFTU(nbands,maxdim2,norbitals)
  INTEGER, intent(in)     :: fh_p, isym, nindo(norbitals), nbands, norbitals, maxdim2
  LOGICAL, intent(in)     :: pform
  ! local
  INTEGER :: iisym, iorb, ind, i
  if (pform) then
     READ(fh_p,*) iisym
  else
     READ(fh_p) iisym
  endif
  if (iisym.NE.isym) print *, 'ERROR: isym and iisym=', isym, iisym
  DO iorb=1,norbitals
     DO ind=1,nindo(iorb)
        if (pform) then
           !READ(fh_p,*) (DMFTU(i,ind,iorb),i=1,nbands)
           DO i=1,nbands
              READ(fh_p,'(f16.9,1x,f16.9,3x)',advance='no') DMFTU(i,ind,iorb)
           ENDDO
           READ(fh_p,*)
        else
           READ(fh_p) (DMFTU(i,ind,iorb),i=1,nbands)
        endif
     ENDDO
  ENDDO
  !DO iorb=1,norbitals
  !   DO ind=1,nindo(iorb)
  !      DO i=1,nbands
  !         WRITE(6,'(I2,1x,I2,1x,I2,1x,I3,1x,2f14.9)') iisym,iorb,ind,i,DMFTU(i,ind,iorb)
  !      ENDDO
  !   ENDDO
  !ENDDO
END SUBROUTINE ReadDMFT_TransK_Inside

SUBROUTINE Read_The_Rest_Basic_Arrays(fhb, nindo, cixdim, nl, ll, iorbital, cix, nind, csize, iSx, Sigind, EF, VOL, norbitals, ncix, natom, maxdim, maxdim2)
  IMPLICIT NONE
  INTEGER, intent(in)  :: fhb
  INTEGER, intent(out) :: nindo(norbitals), cixdim(ncix), nl(natom)
  INTEGER, intent(out) :: ll(natom,4), iorbital(natom,4), cix(natom,4), nind(natom,4)
  INTEGER, intent(out) :: csize(ncix), iSx(maxdim2,norbitals), Sigind(maxdim,maxdim,ncix)
  REAL*8, intent(out)  :: EF, VOL
  INTEGER, intent(in)  :: norbitals, ncix, natom, maxdim, maxdim2
  ! locals
  character*100 :: STR
  INTEGER :: iorb, icix, icase, lcase, ip, iq, ip1
  !!! ReadBasicArrays
  READ(fhb, *) STR ! 'nindo'
  READ(fhb, *) (nindo(iorb),iorb=1,norbitals)
  READ(fhb, *) STR ! 'cidim'
  READ(fhb, *) (cixdim(icix),icix=1,ncix)
  READ(fhb, *) STR ! 'nl'
  READ(fhb, *) (nl(icase), icase=1,natom)
  READ(fhb, *) STR ! 'll, iorbital, cix, nind'
  do icase=1,natom
     do lcase=1,nl(icase)
        READ(fhb, *) ll(icase,lcase)
        READ(fhb, *) iorbital(icase,lcase)
        READ(fhb, *) cix(icase,lcase)
        READ(fhb, *) nind(icase,lcase)          ! nind
     enddo
  enddo

  READ(fhb, *) STR ! csize
  READ(fhb, *) (csize(icix), icix=1,ncix)
  READ(fhb, *) STR ! 'iSx'
  do iorb=1,norbitals
     do ip1=1,nindo(iorb)
        READ(fhb, *) iSx(ip1,iorb)
     enddo
  enddo
  READ(fhb, *) STR ! 'Sigind'
  DO icix=1,ncix
     DO ip=1,cixdim(icix)
        DO iq=1,cixdim(icix)
           READ(fhb, *) Sigind(ip,iq,icix)
        ENDDO
     ENDDO
  ENDDO
  READ(fhb, *) STR ! EF
  READ(fhb, *) EF, VOL
END SUBROUTINE Read_The_Rest_Basic_Arrays

SUBROUTINE Transform_To_Reducible(Ve, isym, Vel, opimat, NE, nord)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Ve(NE,NE,3)
  INTEGER, intent(in)     :: isym
  COMPLEX*16, intent(in)  :: Vel(NE,NE,3)
  REAL*8, intent(in)      :: opimat(3,3,nord)
  INTEGER, intent(in)     :: NE, nord
  ! locals
  INTEGER    :: ip, iq, i, ii
  COMPLEX*16 :: osm
  ! Changing matrix elements to general k-point
  DO ip=1,NE
     DO iq=1,NE
        do i=1,3
           osm=0
           do ii=1,3
              osm = osm + Vel(ip,iq,ii)*opimat(ii,i,isym)
           end do           
           Ve(ip,iq,i) = osm
        end do
     ENDDO
  ENDDO
END SUBROUTINE Transform_To_Reducible


SUBROUTINE Opt_Inside_Simple(copt, Cv, Dv, zekw, k, l, amesh, iw_p, iw_m, wdelta, wgamma, EF, nbands, Nd, bandind, Qinside, nume, Ndirection, NE, nemin, nb_min, nb_max, InterbandOnly)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: copt(Ndirection)
  COMPLEX*16, intent(in) :: Cv(NE,NE,Ndirection), Dv(NE,NE,Ndirection), zekw(NE,2*Nd+1)
  INTEGER, intent(in) :: k, l, amesh(2*Nd+1), iw_p, iw_m, Ndirection
  REAL*8, intent(in)  :: wdelta, wgamma, EF
  INTEGER, intent(in) :: nbands, Nd, bandind(NE), nume, NE, nemin, nb_min, nb_max
  LOGICAL, intent(in) :: InterbandOnly, Qinside(NE)
  ! locals
  REAL*8  :: dx, wdx
  INTEGER :: p, q, dir
  COMPLEX*16  :: eps_p, eps_m, womega, IMAG  
  !------------------------------------------------------------------
  DATA IMAG/(0.0D0,1.0D0)/
  
  eps_p = amesh(iw_p)*wdelta    ! eps
  eps_m = amesh(iw_m)*wdelta    ! eps-omega
  womega = amesh(k+Nd+1)*wdelta ! omega
  
  dx = womega/k
  !dx=(amesh(Nd+2)-amesh(Nd+1))*wdelta
  wdx=dx
  if (l.EQ.1 .or. l.EQ.k+1) wdx=0.5*dx ! Trapezoid rule integration: the first and last interval has 1/2 of the integration weight
  copt = 0
  DO dir=1,Ndirection
     DO p=nb_min,nb_max
        DO q=nb_min,nb_max
           if (InterbandOnly .and. bandind(p).eq.bandind(q)) CYCLE
           if (.NOT.Qinside(p) .AND. .NOT.Qinside(q) ) CYCLE  !! both p and q are outside the window -> non-interacting problem computed outside
           copt(dir) = copt(dir) + wdx*Cv(p,q,dir)/(     (eps_p+EF-zekw(p,iw_p)+wgamma*IMAG)*(eps_m+EF-zekw(q,iw_m)+wgamma*IMAG))
           copt(dir) = copt(dir) - wdx*Dv(p,q,dir)/(conjg(eps_p+EF-zekw(p,iw_p)+wgamma*IMAG)*(eps_m+EF-zekw(q,iw_m)+wgamma*IMAG))
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE Opt_Inside_Simple

SUBROUTINE Opt_Inside(copt, Cv, Dv, zekw, k, l, amesh, iw_p, iw_m, wdelta, wgamma, EF, lg_q, lg_p, nbands, Nd, bandind, Qinside, nume, Ndirection, NE, nemin, nb_min, nb_max, InterbandOnly)
  ! Function computes: 
  !          Cv(ip,iq)*wes - Dv(ip,iq)*wos
  !          where wes and wos are appropraite weights obtained by integrating over frequency
  !          in a small interval of frequency. Frequency points are given by amesh(:)*wdelta
  IMPLICIT NONE
  COMPLEX*16, intent(out)   :: copt(Ndirection)
  COMPLEX*16, intent(in)    :: Cv(NE,NE,Ndirection), Dv(NE,NE,Ndirection), zekw(NE,2*Nd+1)
  INTEGER, intent(in)       :: k, l, amesh(2*Nd+1), iw_p, iw_m
  REAL*8, intent(in)        :: wdelta, wgamma, EF
  COMPLEX*16, intent(inout) :: lg_q(NE), lg_p(NE)
  INTEGER, intent(in)       :: nbands, Nd, bandind(NE), NE, nume, Ndirection, nb_min, nb_max, nemin
  LOGICAL, intent(in)       :: InterbandOnly, Qinside(NE)
  ! locals
  REAL*8  :: dx, wdx, a, b, womega, eps_p, eps_m, big
  INTEGER :: p, q, dir
  COMPLEX*16 :: denom_es, denom_os, wes, wos
  COMPLEX*16 :: IMAG
  !------------------------------------------------------------------
  DATA IMAG/(0.0D0,1.0D0)/
  
  big = 1e6 
  
  eps_p = amesh(iw_p)*wdelta    ! eps
  eps_m = amesh(iw_m)*wdelta    ! eps-omega
  womega = amesh(k+Nd+1)*wdelta ! omega
  

  if (l.EQ.1 .or. iw_p.EQ.1) then        ! The first interval to integrate is smaller
     a=amesh(iw_p)*wdelta
  else 
     a=0.5*(amesh(iw_p-1)+amesh(iw_p))*wdelta
  endif
  if (l.eq.k+1 .or. iw_p.EQ.2*Nd+1) then ! The last interval to integrate is smaller
     b=amesh(iw_p)*wdelta
  else                    ! The rest of the intervals are the same
     b=0.5*(amesh(iw_p)+amesh(iw_p+1))*wdelta
  endif

  DO p=nb_min,nb_max
     lg_q(p) = log(b+EF-womega-zekw(p,iw_m)+wgamma*IMAG) - log(a+EF-womega-zekw(p,iw_m)+wgamma*IMAG)
     lg_p(p) = log(b+EF-       zekw(p,iw_p)+wgamma*IMAG) - log(a+EF-       zekw(p,iw_p)+wgamma*IMAG)
  ENDDO

  copt=0
  DO dir=1,Ndirection
     DO p=nb_min,nb_max
        DO q=nb_min,nb_max
           if (InterbandOnly .and. bandind(p).eq.bandind(q)) CYCLE
           if (.NOT.Qinside(p) .AND. .NOT.Qinside(q) ) CYCLE  !! both p and q are outside the window -> non-interacting problem computed outside
           denom_es = womega -       zekw(p,iw_p)  + zekw(q,iw_m)
           denom_os = womega - conjg(zekw(p,iw_p)) + zekw(q,iw_m) -2*wgamma*IMAG
           wes = (lg_q(q)-      lg_p(p)) /denom_es
           wos = (lg_q(q)-conjg(lg_p(p)))/denom_os
           if (abs(wes).gt.big) then 
              print *, 'Warning: the weightes is very big at ', womega, eps_p, p, q, wes
              wes=0
           endif
           if (abs(wos).gt.big) then
              print *, 'Warning: the weightos is very big at ', womega, eps_p, p, q, wos
              wos=0
           endif
           copt(dir) = copt(dir) + Cv(p,q,dir)*wes - Dv(p,q,dir)*wos
        ENDDO
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE Opt_Inside


SUBROUTINE FindEigsys(zek, Al, Ar, ix, wsigma, gij, womega, omega, STrans, Ek, nom, ncix, nbands, nemin, maxsize, nume, csize, sigma)
  IMPLICIT NONE
  COMPLEX*16, intent(out)   :: zek(nbands), Al(nbands,nbands), Ar(nbands, nbands)
  INTEGER, intent(inout)    :: ix
  COMPLEX*16, intent(inout) :: wsigma(maxsize, ncix)
  COMPLEX*16, intent(inout) :: gij(nbands,nbands)
  REAL*8, intent(in)     :: womega
  REAL*8, intent(in)     :: omega(nom)
  COMPLEX*16, intent(in) :: STrans(maxsize,ncix,nbands,nbands)
  REAL*8, intent(in)     :: Ek(nume)
  INTEGER, intent(in)    :: nom, ncix, nbands, nemin, maxsize, nume, csize(ncix)
  COMPLEX*16, intent(in) :: sigma(nom, maxsize, ncix)
  ! locals
  INTEGER :: icix, ind, i
  
  gij=0
  DO i=1,nbands
     gij(i,i) = Ek(i+nemin-1)  !-----  g^-1 of the LDA part in band representation ----!
  ENDDO
  CALL findNext(womega, ix, omega, nom)              
  do icix=1,ncix
     do ind=1,csize(icix)
        CALL interp(wsigma(ind,icix), sigma(:,ind,icix), womega, ix, omega, nom)
     enddo
  enddo
  CALL AddSigma_optimized2(gij, wsigma, STrans, csize, 1, nbands, ncix, maxsize)           
  CALL eigsys(gij, zek, Al, Ar, nbands)
END SUBROUTINE FindEigsys


SUBROUTINE cmp_C_and_D(Cv, Dv, Al_p, Ar_p, Al_m, Ar_m, Ve, alphaV, Ndirection, L1, L2, L3, L4, Ltmp, Qinside, ndim, NE, nemin, nemax, nb_min, nb_max)
  !--------------------------------------
  IMPLICIT NONE
  !---------- Passed variables ----------     
  COMPLEX*16, intent(out)   :: Cv(NE,NE,Ndirection), Dv(NE,NE,Ndirection)
  COMPLEX*16, intent(in)    :: Al_p(ndim,ndim), Ar_p(ndim,ndim), Al_m(ndim,ndim), Ar_m(ndim,ndim)
  COMPLEX*16, intent(in)    :: Ve(NE,NE,3)
  REAL*8, intent(in)        :: alphaV(3,3,Ndirection)
  COMPLEX*16, intent(inout) :: L1(NE,NE,3), L2(NE,NE,3), L3(NE,NE,3), L4(NE,NE,3), Ltmp(ndim,ndim)
  INTEGER, intent(in)       :: ndim, Ndirection, NE, nb_min, nb_max, nemin, nemax
  LOGICAL, intent(in)       :: Qinside(NE)
  !---------- Local variables -----------
  COMPLEX*16 :: Al_p_c(ndim,ndim), Ar_p_c(ndim,ndim)
  COMPLEX*16 :: c1tmp, c2tmp
  INTEGER    :: p,q, i,j,l, dir

  ! Calculation of matrices C and D that are multiplied to the double denominator
  ! C is for equal sign and D for opposite sign
  Al_p_c = conjg(transpose(Al_p))
  Ar_p_c = conjg(transpose(Ar_p))
  
  DO i=1,3
     ! inside the small window
     Ltmp = matmul(Al_m, Ve(nemin:nemax,nemin:nemax,i))
     L1(nemin:nemax,   nemin:nemax,   i) = matmul(Ltmp, Ar_p)
     L3(nemin:nemax,   nemin:nemax,   i) = matmul(Ltmp, Al_p_c)
     
     ! cross terms, where one index is inside window, and one is outside
     if (nb_min.LT.nemin) then
        L1(nb_min:nemin-1,nemin:nemax,   i) = matmul(Ve(nb_min:nemin-1,nemin:nemax,i),Ar_p)
        L1(nemin:nemax,   nb_min:nemin-1,i) = matmul(Al_m, Ve(nemin:nemax,nb_min:nemin-1,i))
        L3(nb_min:nemin-1,nemin:nemax,   i) = matmul(Ve(nb_min:nemin-1,nemin:nemax,i), Al_p_c)
        L3(nemin:nemax,   nb_min:nemin-1,i) = L1(nemin:nemax,nb_min:nemin-1,i)
     endif
     if (nb_max.GT.nemax) then
        L1(nemax+1:nb_max,nemin:nemax,   i) = matmul(Ve(nemax+1:nb_max,nemin:nemax,i),Ar_p)
        L1(nemin:nemax,   nemax+1:nb_max,i) = matmul(Al_m, Ve(nemin:nemax,nemax+1:nb_max,i))
        L3(nemax+1:nb_max,nemin:nemax,   i) = matmul(Ve(nemax+1:nb_max,nemin:nemax,i), Al_p_c)
        L3(nemin:nemax,nemax+1:nb_max,   i) = L1(nemin:nemax,nemax+1:nb_max,i)
     endif
     
     ! inside the small window
     Ltmp = matmul(Ve(nemin:nemax,nemin:nemax,i), Ar_m)
     L2(nemin:nemax,   nemin:nemax,   i) = matmul(Al_p, Ltmp)
     L4(nemin:nemax,nemin:nemax,i) = matmul(Ar_p_c, Ltmp)
     
     ! cross terms, where one index is inside window, and one is outside
     if (nb_min.LT.nemin) then
        L2(nb_min:nemin-1,nemin:nemax,   i) = matmul(Ve(nb_min:nemin-1,nemin:nemax,i),Ar_m)
        L2(nemin:nemax,   nb_min:nemin-1,i) = matmul(Al_p, Ve(nemin:nemax,nb_min:nemin-1,i))
        L4(nb_min:nemin-1,nemin:nemax,   i) = L2(nb_min:nemin-1,nemin:nemax,   i)
        L4(nemin:nemax,   nb_min:nemin-1,i) = matmul(Al_p_c, Ve(nemin:nemax,nb_min:nemin-1,i))
     endif
     if (nb_max.GT.nemax) then
        L2(nemax+1:nb_max,nemin:nemax,   i) = matmul(Ve(nemax+1:nb_max,nemin:nemax,i),Ar_m)
        L2(nemin:nemax,   nemax+1:nb_max,i) = matmul(Al_p, Ve(nemin:nemax,nemax+1:nb_max,i))
        L4(nemax+1:nb_max,nemin:nemax,   i) = L2(nemax+1:nb_max,nemin:nemax,   i)
        L4(nemin:nemax,   nemax+1:nb_max,i) = matmul(Al_p_c, Ve(nemin:nemax,nemax+1:nb_max,i))
     endif
  ENDDO
  
  DO p=nb_min,nb_max
     DO q=nb_min,nb_max
        ! When both p and q are outside, we do not need to update Cv and Dv, because they 
        ! are equal to v^2.
        if (.NOT.Qinside(p) .AND. .NOT.Qinside(q) ) CYCLE
        DO dir=1,Ndirection
           c1tmp=0
           c2tmp=0
           DO j=1,3
              DO l=1,3
                 c1tmp = c1tmp + L2(p,q,j)*L1(q,p,l)*alphaV(l,j,dir)
                 c2tmp = c2tmp + L4(p,q,j)*L3(q,p,l)*alphaV(l,j,dir)
              ENDDO
           ENDDO
           Cv(p,q,dir) = c1tmp
           Dv(p,q,dir) = c2tmp
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE cmp_C_and_D

SUBROUTINE cmp_C_and_D_old(Cv, Dv, Al_p, Ar_p, Al_m, Ar_m, Ve, alphaV, Ndirection, L1, L2, L3, L4, Ltmp, ndim)
  !--------------------------------------
  IMPLICIT NONE
  !---------- Passed variables ----------     
  COMPLEX*16, intent(out)   :: Cv(ndim,ndim,Ndirection), Dv(ndim,ndim,Ndirection)
  COMPLEX*16, intent(in)    :: Al_p(ndim,ndim), Ar_p(ndim,ndim), Al_m(ndim,ndim), Ar_m(ndim,ndim)
  COMPLEX*16, intent(in)    :: Ve(ndim,ndim,3)
  REAL*8, intent(in)        :: alphaV(3,3,Ndirection)
  COMPLEX*16, intent(inout) :: L1(ndim,ndim,3), L2(ndim,ndim,3), L3(ndim,ndim,3), L4(ndim,ndim,3), Ltmp(ndim,ndim)
  INTEGER, intent(in)       :: ndim, Ndirection
  !---------- Local variables -----------
  COMPLEX*16 :: c1tmp, c2tmp
  INTEGER    :: p,q, i,j,l, dir

  ! Calculation of matrices C and D that are multiplied to the double denominator
  ! C is for equal sign and D for opposite sign
  DO i=1,3
     Ltmp = matmul(Al_m, Ve(:,:,i))
     L1(:,:,i) = matmul(Ltmp, Ar_p)
     L3(:,:,i) = matmul(Ltmp, conjg(transpose(Al_p)))
     Ltmp = matmul(Ve(:,:,i), Ar_m)
     L2(:,:,i) = matmul(Al_p, Ltmp)
     L4(:,:,i) = matmul(conjg(transpose(Ar_p)), Ltmp)
  ENDDO
  Cv = 0
  Dv = 0
  DO p=1,ndim
     DO q=1,ndim
        DO dir=1,Ndirection
           c1tmp=0
           c2tmp=0
           DO j=1,3
              DO l=1,3
                 c1tmp = c1tmp + L2(p,q,j)*L1(q,p,l)*alphaV(l,j,dir)
                 c2tmp = c2tmp + L4(p,q,j)*L3(q,p,l)*alphaV(l,j,dir)
              ENDDO
           ENDDO
           Cv(p,q,dir) = Cv(p,q,dir) + c1tmp
           Dv(p,q,dir) = Dv(p,q,dir) + c2tmp
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE cmp_C_and_D_old


SUBROUTINE Opt_Noncorrelated(copt, Cv0, Ek, womega, EF, wgamma,nb_min, nb_max, nemin, nemax, nume, Qinside, NE, Ndirection, InterbandOnly)
!!! Add non-correlated bands !!!  
  IMPLICIT NONE
!!! For non-correlated
  COMPLEX*16, intent(out) :: copt(Ndirection)
  COMPLEX*16, intent(in) :: Cv0(NE,NE,Ndirection)
  REAL*8, intent(in)     :: Ek(nume)
  REAL*8, intent(in)  :: womega, EF
  REAL*8, intent(in)  :: wgamma
  INTEGER, intent(in) :: nb_min, nb_max, nemin, nemax, nume, NE, Ndirection
  LOGICAL, intent(in) :: InterbandOnly, Qinside(NE)
  ! local variables
  REAL*8     :: aw, bw
  INTEGER    :: ip, iq, dir
  COMPLEX*16 :: denom_es, denom_os, wes, wos
  COMPLEX*16, ALLOCATABLE :: lg_p(:), lg_q(:)
  COMPLEX*16, PARAMETER  :: IMAG = (0.0D0,1.0D0)

  ALLOCATE( lg_p(nb_max), lg_q(nb_max) )
  
  aw=0.0
  bw=womega
  DO ip=nb_min,nb_max
     if (Qinside(ip)) CYCLE ! These bands were taken into account in DMFT part
     lg_q(ip) = log(bw+EF-womega-Ek(ip)+wgamma*IMAG) - log(aw+EF-womega-Ek(ip)+wgamma*IMAG)
     lg_p(ip) = log(bw+EF-       Ek(ip)+wgamma*IMAG) - log(aw+EF-       Ek(ip)+wgamma*IMAG)
  ENDDO
  copt=0
  DO ip=nb_min,nb_max

     !if (InterbandOnly .OR. (ip.GE.nemin .and.ip.LE.nemax)) CYCLE ! These bands were taken into account in DMFT part
     DO iq=nb_min,nb_max

        if ( Qinside(ip) .or. Qinside(iq) ) CYCLE ! These bands were taken into account in DMFT part

        denom_es = womega - Ek(ip) + Ek(iq)
        denom_os = womega - Ek(ip) + Ek(iq) -2*wgamma*IMAG
        wes = (lg_q(iq)-      lg_p(ip)) /denom_es
        wos = (lg_q(iq)-conjg(lg_p(ip)))/denom_os
        DO dir=1,Ndirection
           copt(dir) = copt(dir) + Cv0(ip,iq,dir)*(wes - wos)
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE( lg_p, lg_q )
END SUBROUTINE Opt_Noncorrelated



program read_optm
  USE com_mpi, ONLY: start_MPI, stop_MPI
  IMPLICIT NONE
  ! functions
  INTEGER :: CountSelfenergy
  ! variables
  CHARACTER*100 :: case, fmommat, fsymop, fUdmft, fenergy, fbasicArrays, fsigname
  CHARACTER*2   :: updn, dnup, so
  INTEGER :: fhb, fhi, fh_m
  CHARACTER*100 :: STR
  INTEGER :: nom, nk, iso, nat, norbitals, ncix, natom, nkpt, nmat, nume, lmax2, maxdim2, maxdim, maxsize
  LOGICAL :: Qsym, Qsimple, Qcomplex, InterbandOnly, ProjectToCorrelated
  REAL*8  :: gamma, gammac, ommax, delta, dwindow, Temperature
  INTEGER :: Nd, Nd0, i, j, l
  CHARACTER*5 :: adum
  INTEGER :: Ndirection
  LOGICAL :: file_exists
  REAL*8, ALLOCATABLE :: alphaV(:,:,:)

  
  CALL start_MPI()

  Qsimple = .False.
  updn=''
  fhi = 995
  open(fhi, file='dmftopt.in', status='old', form='formatted')
  READ(fhi, *) Temperature
  READ(fhi, *) case
  READ(fhi, *) gamma
  READ(fhi, *) gammac
  READ(fhi, *) ommax
  READ(fhi, *) delta
  READ(fhi, *) Nd0
  READ(fhi, *) Qsym
  READ(fhi, *) InterbandOnly
  READ(fhi, *) dwindow
  READ(fhi, *) Ndirection
  print *, 'Ndirection=', Ndirection
  ALLOCATE( alphaV(3,3,Ndirection) )
  do l=1,Ndirection
     READ(fhi, *) alphaV(1,1,l), alphaV(1,2,l), alphaV(1,3,l)
     READ(fhi, *) alphaV(2,1,l), alphaV(2,2,l), alphaV(2,3,l)
     READ(fhi, *) alphaV(3,1,l), alphaV(3,2,l), alphaV(3,3,l)
  enddo
  close(fhi)

  !if (.NOT.(updn.EQ.'up' .OR.updn.EQ.'dn')) updn=''
  !!
  !dnup='dn'
  !IF (updn.EQ.'dn') THEN
  !   dnup = 'up'
  !ELSEIF (updn.EQ.'up') THEN
  !   dnup = 'dn'
  !ELSE
  !   updn = ''
  !ENDIF
  !!
  
  !!! Counting the number of frequency points in self-energy
  fsigname = 'sig.inp'//TRIM(updn)//'1'
  print *, 'fsigname=', fsigname
  open(81, file=fsigname, status='old', form='formatted')
  nom = CountSelfenergy(81, 1)-1 !--- How many frequency point exists in the input file? ---!
  close(81)
  
  fbasicArrays = 'BasicArrays.dat'
  !!! Reading the first part of the header file to get dimensions
  fhb=996
  open(fhb,file=fbasicArrays,status='old', form='formatted')
  READ(fhb, *) STR ! 'nat, iso, norbitals, ncix, natom'
  READ(fhb, *) nat, iso, norbitals, ncix, natom
  READ(fhb, *) STR ! nkpt, nmat, nume
  READ(fhb, *) nkpt, nmat, nume
  READ(fhb, *) STR ! Qcomplex
  READ(fhb, *) Qcomplex
  READ(fhb, *) STR ! 'lmax2, maxdim2, maxdim, maxsize'
  READ(fhb, *) lmax2, maxdim2, maxdim, maxsize

  
  so = ''
  IF (iso.EQ.2) so='so'
  fenergy = TRIM(case)//'.energy'//TRIM(so)//TRIM(updn)
  fmommat = TRIM(case)//'.mommat'//TRIM(updn)
  fsymop = TRIM(case)//'.symop'
  fUdmft = 'Udmft'//TRIM(updn)//'.'

  fh_m = 99  ! file for velocity matrix elements
  
  INQUIRE(file=fmommat, exist=file_exists)
  if (.not. file_exists) then
     fmommat = TRIM(case)//'.mommat2'//TRIM(updn)
     INQUIRE(file=fmommat, exist=file_exists)
     if (.not. file_exists) WRITE(*,*) 'ERROR File not found: Expecting that ', fmommat, 'exists.'
  endif
  open(fh_m, file=fmommat, status='old')

  if (Temperature.NE.0.0) then
     ! We will use linear mesh only. We need to set very large Nd, so that
     ! cutoff >=4*T, and ommax/(2*delta)=Nd
     if (Nd0 .LT. int(2*Temperature/Delta+1.)) Nd0=int(2*Temperature/Delta+1.)
     if (ommax.LT.2*Nd0*delta) then
        ommax=2*Nd0*delta
     else
        Nd0 = ommax/(2*delta)
     endif
  endif

  
  WRITE(*,'(A,A)') 'case=', TRIM(case)
  WRITE(*,'(A,f10.7,1x,A,f10.7)') 'gamma=', gamma, 'gammac=', gammac
  WRITE(*,'(A,f10.3,1x,A,f14.10,1x,A,I4)') 'ommax=', ommax, 'delta=', delta, 'Nd=', Nd0
  WRITE(*,'(A,L)') 'Qsym=', Qsym
  WRITE(*, '(A)') 'alphaV='
  WRITE(*,*) (((alphaV(j,i,l), i=1,3), j=1,3),l=1,Ndirection)
  WRITE(*,'(A,A)') 'fmommat=', fmommat
  WRITE(*,'(A,A)') 'fenergy=', fenergy
  WRITE(*,'(A,A)') 'fUdmft=', fUdmft
  WRITE(*,'(A,A)') 'fsymop=', fsymop
  WRITE(*,'(A,L)') 'InterbandOnly=', InterbandOnly
  WRITE(*,'(A,f10.7)') 'dwindow=', dwindow
  WRITE(*,'(A,f10.7)') 'Temperature=', Temperature


  Nd = 2*Nd0
  CALL Cmp_Optics(fenergy, fUdmft, fsymop, fh_m, fhb, ommax, Temperature, delta, Nd, gammac, gamma, alphaV, Ndirection, Qsym, Qsimple, nom, nat, iso, norbitals, ncix, natom, nkpt, nmat, nume, Qcomplex, lmax2, maxdim2, maxdim, maxsize, dwindow, InterbandOnly)
  
  close(fh_m)
  DEALLOCATE( alphaV )
  CALL stop_MPI()
  
end program read_optm

REAL*8 Function ferm(x)
  IMPLICIT NONE
  REAL*8, intent(in):: x
  if (x.gt.100) then
     ferm=0.0
     return
  endif
  if (x.lt.-100) then
     ferm=1.0
     return
  endif
  ferm = 1./(1.+exp(x))
  return
end function ferm
