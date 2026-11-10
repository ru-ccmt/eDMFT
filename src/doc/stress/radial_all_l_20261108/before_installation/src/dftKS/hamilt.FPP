#ifdef _REAL
#define tpcast(X) dble(X)  !! typecasting in fortran
#define XCOPY DCOPY
#define Ttype REAL*8
#define tsize 8
#ifdef _INTEL_VML
! intel vmlib & _REAL
#define vcos_sin(res,X,n) vdcos(n,X,res)
#else
! non-intel & _REAL
#define vcos_sin(res,X,n) vcos(res,X,n)
#endif
#else
#define tpcast(X) X
#define XCOPY ZCOPY
#define Ttype COMPLEX*16
#define tsize 16
#ifdef _INTEL_VML
! intel vmlib & complex
#define vcos_sin(res,X,n) vzcis(n,X,res)
#else
! non-intel & complex
#define vcos_sin(res,X,n) vcosisin(res,X,n)
#endif
#endif

#ifdef _INTEL_VML
#define vsin_cos(res_sin,res_cos,X,n) vdsincos(n,X,res_sin,res_cos)
#else
#define vsin_cos(res_sin,res_cos,X,n) vsincos(res_sin,res_cos,X,n)
#endif

SUBROUTINE HAMILTON(LMX,NAT,NV,U,UP,DU,DUP,MULT,POS,V,Vol,pei,E)
  use matrices, only : HS, HSDIAG, HSROWS, KZZ, RK, Kn
  use loabc, only : ALO, BLO, CLO, ELO, DPLO, PI12LO, PE12LO,pilolo,nloat_new
  use lapw_timer, only : time_albl, time_loop260, time_phase,time_legendre, time_step_function, time_h,time_us, time_loop210, time_loop230,time_loop240, time_overlap, time_lo, START_TIMER, STOP_TIMER, READ_CPU_TIME, READ_wall_TIME,init_all_timer,time_distiouter,time_hamilt
  use lolog, only : nlo, lapw, ilo, nlov
  use structure, only: ndf, ROTLOC, rotij, Rmt
  use albl, ONLY: init_albl, end_albl, Al_c, Al_r, Bl_c, Bl_r
  use param, ONLY: lmax, lomax, nloat, blocksize
  use mpi, ONLY: Qprint
  IMPLICIT NONE
  !        Arguments
  INTEGER, intent(in) :: LMX, NAT, NV  ! max-l, n-sorts, N-plane waves
  REAL*8, intent(in)  :: U(0:LMAX-1,NAT), UP(0:LMAX-1,NAT), DU(0:LMAX-1,NAT), DUP(0:LMAX-1,NAT) ! radial wave functions
  INTEGER, intent(in) :: MULT(NAT) ! how many atoms of the same sort
  REAL*8, intent(in)  :: POS(3,NDF) ! atom positions
  REAL*8, intent(in)  :: V(NAT), Vol  ! Relative volume of MT-sphere, Volume of the unit cell
  REAL*8, intent(in)  :: pei(0:LMAX-1,NAT) ! <dotu|dotu>    ! STRANGE: why not get from module??
  REAL*8, intent(in)  :: E(0:LMAX-1,NAT)  ! linearization energies , STRANGE Elo is taken from module
  !
  !        External Subroutines
  EXTERNAL  :: WARPIN, CPUTIM, DVBES1, ROTATE, SPHBES, YLM
  ! Constants
  DOUBLE PRECISION   TOL
  PARAMETER          (TOL=1.0D-8)
  !        Locals
  Ttype, allocatable :: Hamilt(:,:), Sk_all(:,:)
  Ttype, allocatable :: PHASE(:), Sk_ints(:,:), help_exp(:)
  Ttype      :: phs, WARPIN
  REAL*8     :: V3, atmp0, tmp3, stmp, ARGX, ATMP, BTMP, PI4R2V, PI4R2V2, PI, VIFPR4, akinlo, C1, C11, C12, C2, atmpa, atmpb, atmpab, atmpba, wtmpa, wtmpb, den_j, x0
  REAL*8     :: DTIMH, DTIMHw, DTIMLG, DTIMMA, DTIMPH,DTIMDIST, time_albl_tot, time_albl_totw, DTIMLGw, DTIMMAw, DTIMPHw,DTIMDISTw,DTIMUSw, DTIMS, DTIMUS
  INTEGER    :: jlo_i,jlop,jlo_j, jnlo, m_j, mu_j, lomaxlomax, m, inlo, i_g, j_g, multmax
  INTEGER    :: i, j, mup, mu_i, mp, mu, jatom, l, lmp, m_i, ihelp, iouter, imin, imax, latom, lfirst
  COMPLEX*16 :: C11UP, C12UP, cYp, AcYp, BcYp, CcYp, C11UP_kinlo, C12UP_kinlo, CLYLN, C1UP, C2UP, akinloup
  REAL*8     :: rotloc_x_BR2_x_rotij(3,3,ndf), ROTV2(3), XL(lmax-2), atom_pos(3)
  COMPLEX*16, allocatable :: YL(:,:,:), C3UP(:), C13UP(:), PHSC(:)
  REAL*8,     allocatable :: C3(:), C13(:), PPLX(:), PPLY(:), SLEN(:,:), TMP1(:), TMP2(:)!, DFJ(:), FJ(:), 
  REAL*8,     allocatable :: KxK(:,:), x(:), Pl(:,:,:), help_x(:), help_cos(:), j1x(:), help_sin(:), tmp_y(:)
  INTEGER,    allocatable :: iind(:,:,:)

  CALL START_TIMER(time_hamilt)
  nloat=nloat_new
  allocate( Hamilt( HSrows, BLOCKSIZE) )
  allocate( Sk_all( HSrows, BLOCKSIZE) )
  ! monsters
  allocate( Pl(HSrows,0:LMX-1,blocksize) )
  allocate( SLEN(HSrows,blocksize) )
  allocate( KxK(HSrows,blocksize) )
  allocate( Sk_ints(HSrows,blocksize))
  ! large
  allocate( PPLX(HSrows), PPLY(HSrows), TMP1(HSrows), TMP2(HSrows), X(HSrows))
  allocate( help_x(HSrows), help_cos(HSrows), help_sin(HSrows), help_exp(HSrows) )
  allocate( j1x(HSrows), tmp_y(HSrows), PHASE(HSrows) )

  CALL init_albl(lmx-1,'HAM')

  multmax = maxval(mult)
  if (Qprint) then
     write(6,'(a20,f12.1,a,a,2i6)') 'allocate Sk_all',HSrows*blocksize*tsize/(1024.0*1024.0),' MB ','         dimensions',HSrows,blocksize
     write(6,'(a20,f12.1,a,a,2i6)') 'allocate hamilt',HSrows*blocksize*tsize/(1024.0*1024.0),' MB ','         dimensions',HSrows,blocksize
     write(6,'(a20,f12.1,a,a,2i6)') 'allocate Sk_ints',HSrows*blocksize*tsize/(1024.0*1024.0),' MB ','         dimensions',HSrows,blocksize
     write(6,'(a20,f12.1,a,a,2i6)') 'allocate slen',HSrows*blocksize*8.0/(1024.0*1024.0),' MB ','         dimensions',HSrows,blocksize
     write(6,'(a20,f12.1,a,a,2i6)') 'allocate KxK',HSrows*blocksize*8.0/(1024.0*1024.0),' MB ','         dimensions',HSrows,blocksize
     write(6,'(a20,f12.1,a,a,3i6)') 'allocate legendre',dble((HSrows)*LMAX*blocksize*8)/(1024.0*1024.0),' MB ','         dimensions',HSrows,lmax,blocksize
     write(6,'(a20,f12.1,a,a,2i6)') 'allocate al,bl (row)',HSrows*lmx*2.0*2*8.0/(1024.0*1024.0),' MB ','         dimensions',HSrows,lmx
     write(6,'(a20,f12.1,a,a,2i6)') 'allocate al,bl (col)',blocksize*lmx*2.0*2*8.0/(1024.0*1024.0),' MB ','         dimensions',blocksize,lmx
     write(6,'(a20,f12.1,a,a,3i6)') 'allocate YL',dble((LOMAX+1)*(LOMAX+1)*HSrows*multmax*16)/(1024.0*1024.0),' MB ','         dimensions',(LOMAX+1)*(LOMAX+1)-1,HSrows,multmax
  endif
  
  PI = 4.d0*atan(1.d0)
  DO I=1,LMAX-2
     XL(I) = DBLE(I)/DBLE(I+1)
  enddo
  DTIMPH = 0
  DTIMUS = 0
  DTIMLG = 0
  DTIMMA = 0
  DTIMS = 0
  DTIMH = 0
  DTIMHw = 0
  DTIMDIST = 0
  time_albl_tot = 0
  DTIMPHw = 0
  DTIMUSw = 0
  DTIMLGw = 0
  DTIMMAw = 0
  DTIMDISTw = 0
  time_albl_totw= 0
  CALL START_TIMER(time_loop260)

  DO IOUTER=0, NV/(BLOCKSIZE)

     imin = IOUTER*BLOCKSIZE + 1
     imax = MIN(NV, (IOUTER + 1)*BLOCKSIZE)
     !        precompute Legendre-Polynomials
     CALL START_TIMER(time_legendre)
     DO i_g=imin,imax
        ihelp = i_g-imin+1
        ! Notice that Kn(:,i_g) is k+K in cartesian global coordinate system of conventional unit cell, i.e., Kn = BR2*ik
        DO j=1,i_g  ! This loop build xi = (K_i+k)(K_j+k)/(|K_i+k||K_j+k|)
           KxK(j,ihelp) =  dot_product(Kn(:,i_g),Kn(:,j))  ! KxK = (K_i+k)(K_j+k)
           SLEN(j,ihelp) = sqrt((Kn(1,j)-Kn(1,i_g))**2+(Kn(2,j)-Kn(2,i_g))**2+(Kn(3,j)-Kn(3,i_g))**2)  ! slen=|K_j-K_i|
           X(j)=KxK(j,ihelp)   ! X = (K_i+k)(K_j+k)
           den_j=RK(i_g)*RK(j) ! den_j = |K_i+k||K_j+k|
           IF (den_j.GT.tol) X(j) = X(j)/den_j ! X = (K_i+k)(K_j+k)/(|K_i+k||K_j+k|)
        enddo
        PPLX(:i_g) = X(:i_g)*X(:i_g)
        PPLY(:i_g) = PPLX(:i_g) - 1.d0
        Pl(:i_g,0,ihelp) = 1.d0
        Pl(:i_g,1,ihelp) = X(:i_g)
        ! Recursion for Legendre Polynomials
        DO l=2,lmx-2
           Pl(:i_g,l,ihelp) = PPLX(:i_g) + XL(l-1)*PPLY(:i_g)
           PPLX(:i_g) = X(:i_g)*Pl(:i_g,l,ihelp)
           PPLY(:i_g) = PPLX(:i_g) - Pl(:i_g,l-1,ihelp)
        enddo
        Pl(:i_g,lmx-1,ihelp) = PPLX(:i_g) + XL(lmx-2)*PPLY(:i_g)
        ! correcting for x~1.0 where recursion does not work very well
        do j=1,i_g
           IF (ABS(x(j)) .GT. 0.9999999D+0) THEN
              do l=2,lmx-1
                 Pl(j,l,ihelp) = x(j)*Pl(j,l-1,ihelp)
              enddo
           ENDIF
        enddo
     enddo
     
     CALL STOP_TIMER(time_legendre)
     DTIMLG = DTIMLG + READ_CPU_TIME(time_legendre)
     DTIMLGw = DTIMLGw + READ_wall_TIME(time_legendre)

     Sk_all(:,:) = 0.d0
     Sk_ints(:,:)=0.d0
     Hamilt(:,:) = 0.d0
     DO jatom=1,nat
        lfirst = 1+sum(mult(1:jatom-1))
        VIFPR4 = 4.d0*pi*Rmt(jatom)**4/Vol

        CALL START_TIMER(time_albl)
        !        precompute al(kn) and bl(kn)
        call make_albl(jatom,Rmt(jatom),lmx-1,NV,'HAM',imin,imax)
        CALL STOP_TIMER(time_albl)
        time_albl_tot=time_albl_tot+READ_CPU_TIME(time_albl)
        time_albl_totw=time_albl_totw+READ_wall_TIME(time_albl)
        imin = IOUTER*BLOCKSIZE + 1
        imax = MIN(NV, (IOUTER + 1) * BLOCKSIZE)

        DO i_g =imin,imax
           ihelp = i_g-imin+1
           !        calculate phase-factors
           CALL START_TIMER(time_phase)
           PHASE(:) = 0.d0
           DO mu=1,mult(jatom)
              latom = lfirst + mu -1
              atom_pos(:) = 2.d0*Pi*POS(:,latom)
              DO j=1,i_g
                 help_x(j)= (KZZ(1,j)-KZZ(1,i_g))*atom_pos(1)+(KZZ(2,j)-KZZ(2,i_g))*atom_pos(2)+(KZZ(3,j)-KZZ(3,i_g))*atom_pos(3)
              enddo
              call vcos_sin(help_exp,help_x,i_g)    ! exp(i*help_x)
              PHASE(1:i_g) = PHASE(1:i_g) + help_exp(1:i_g)
           ENDDO
           
!           !$OMP PARALLEL DO PRIVATE(latom,atom_pos,x0) REDUCTION(+:PHASE)
!           DO mu=1,mult(jatom)
!              latom = lfirst + mu -1
!              atom_pos(:) = 2.d0*Pi*POS(:,latom)
!              do j=1,i_g
!                 x0 = (KZZ(1,j)-KZZ(1,i_g))*atom_pos(1)+(KZZ(2,j)-KZZ(2,i_g))*atom_pos(2)+(KZZ(3,j)-KZZ(3,i_g))*atom_pos(3)
!#ifdef _REAL
!                 PHASE(j) = PHASE(j) + cos(x0)
!#else                 
!                 PHASE(j) = PHASE(j) + dcmplx(cos(x0),sin(x0))
!#endif                 
!              enddo
!           ENDDO
!           !$OMP END PARALLEL DO
           
           CALL STOP_TIMER(time_phase)
           DTIMPH = DTIMPH + READ_CPU_TIME(time_phase)
           DTIMPHw = DTIMPHw + READ_wall_TIME(time_phase)
           CALL START_TIMER(time_us)
           CALL START_TIMER(time_step_function)
           
           ! Interstitial part of Hamiltonian and Overlap.
           help_x(1:i_g-1) = slen(1:i_g-1,ihelp)*Rmt(jatom)  ! |K_i-K_j|*Rmt
           tmp_y(1:i_g-1) = 1.d0/help_x(1:i_g-1)**3          ! 1/(|K_i-K_j|*Rmt)**3
           call vsin_cos(help_sin,help_cos,help_x,i_g-1)     ! sin(|K_i-K_j|*Rmt), cos(|K_i-K_j|*Rmt)
           v3=3.d0*V(jatom)   ! 3.*V_MT/V_cell
           ! j1(x)/x = (x*cos(x)-sin(x))/x**3 where x=|K_i-K_j|*Rmt 
           ! O += 3 * V_MT/V_cell * phase * j1(x)/x |x=(|K_i-K_j|*Rmt)
           j1x(1:i_g-1) = -(help_x(1:i_g-1)*help_COS(1:i_g-1)-help_SIN(1:i_g-1))*tmp_y(1:i_g-1)
           Sk_ints(1:i_g-1,ihelp) = Sk_ints(1:i_g-1,ihelp) - v3*PHASE(1:i_g-1)*j1x(1:i_g-1)
           Sk_ints(i_g,ihelp) = Sk_ints(i_g,ihelp) + 1.d0/nat - V(jatom)*mult(jatom)  ! if j==i => K_i-K_j==0 => j1(x)/x=1/3.
           
           CALL STOP_TIMER(time_step_function)
           DTIMS = DTIMS + READ_CPU_TIME(time_step_function)
           
           CALL STOP_TIMER(time_us)
           DTIMUS = DTIMUS + READ_CPU_TIME(time_us)
           DTIMUSw = DTIMUSw + READ_wall_TIME(time_us)
           
           CALL START_TIMER(time_overlap)
           !        calculate Overlap
           ! Muffin-thin part of the symmetric Hamiltonian
           tmp1(:) = 0.d0
           tmp2(:) = 0.d0
           DO l=0,lmx-1
              if (.not.lapw(l,jatom)) then
                 atmp = AL_c(ihelp,l)
                 atmp0 = U(l,jatom)*DU(l,jatom)*Rmt(jatom)**2
		 !  !$OMP PARALLEL DO PRIVATE(stmp)
                 DO j=1,i_g
                    !  surface part of kinetic energy 
                    stmp = (2*l+1)*Pl(j,l,ihelp)*atmp*AL_r(j,l)
                    tmp2(j) = tmp2(j) + stmp
                    tmp1(j) = tmp1(j) + (E(l,jatom) + atmp0)*stmp
                 end do
		 !  !$OMP END PARALLEL DO
              else
                 atmpa = u(l,jatom)*du(l,jatom)*Rmt(jatom)**2
                 atmpb = up(l,jatom)*dup(l,jatom)*Rmt(jatom)**2
                 atmpab = u(l,jatom)*dup(l,jatom)*Rmt(jatom)**2
                 atmpba = up(l,jatom)*du(l,jatom)*Rmt(jatom)**2
                 atmp = Al_c(ihelp,L)
                 btmp = Bl_c(ihelp,L)*pei(L,JATOM)
                 wtmpa = (E(l,jatom)+atmpa)*atmp + atmpba*Bl_c(ihelp,l)
                 wtmpb = (E(l,jatom)*pei(l,jatom) + atmpb)*Bl_c(ihelp,L) + (1.d0+atmpab)*atmp
		 !  !$OMP PARALLEL DO PRIVATE(tmp3)
                 DO j=1,i_g
                    tmp3 = (2*l+1)*Pl(j,l,ihelp)
                    tmp2(j) = tmp2(j) + ( atmp*Al_r(j,l) +  btmp*Bl_r(j,l))*tmp3
                    tmp1(j) = tmp1(j) + (wtmpa*Al_r(j,l) + wtmpb*Bl_r(j,l))*tmp3
                 end do
		 !  !$OMP END PARALLEL DO
              endif
           enddo
	   !	   !$OMP PARALLEL DO PRIVATE(phs)
           DO j=1,i_g
              phs = phase(j)*VIFPR4
              Sk_all(j,ihelp) = Sk_all(j,ihelp) + tmp2(j)*phs
              Hamilt(j,ihelp) = Hamilt(j,ihelp) + tmp1(J)*phs
           enddo
           !	!$OMP END PARALLEL DO
	   
           CALL STOP_TIMER(time_overlap)
           DTIMMA = DTIMMA + READ_CPU_TIME(time_overlap)
           DTIMMAw = DTIMMAw + READ_wall_TIME(time_overlap)
        END DO  ! over i_g
     enddo      ! over jatom

     ! Finally finish the interstitials
     CALL START_TIMER(time_h)
     DO i_g=imin,imax
        ihelp = i_g-imin+1
        Sk_all(1:i_g,ihelp) = Sk_all(1:i_g,ihelp) + Sk_ints(1:i_g,ihelp)
        Hamilt(1:i_g,ihelp) = Hamilt(1:i_g,ihelp) + KxK(1:i_g,ihelp)*Sk_ints(1:i_g,ihelp)
        DO j=1,i_g
           Hamilt(j,ihelp) = Hamilt(j,ihelp) + WARPIN(KZZ(1,j)-KZZ(1,i_g),KZZ(2,j)-KZZ(2,i_g),KZZ(3,j)-KZZ(3,i_g))
        enddo
     ENDDO
     CALL STOP_TIMER(time_h)
     DTIMH  = DTIMH + READ_CPU_TIME(time_h)
     DTIMHw = DTIMHw + READ_wall_TIME(time_h)
     
     CALL start_TIMER(time_distiouter)

     imin = IOUTER*BLOCKSIZE + 1
     imax = MIN(NV, (IOUTER+1)* BLOCKSIZE)
     DO i_g=imin,imax
        ihelp = i_g-imin+1
        CALL XCOPY(i_g-1, Hamilt(1,ihelp), 1, HS(i_g, 1), HSrows )  ! HS(i_g,:) <- Hamilt(:,ihelp)
        CALL XCOPY(i_g, Sk_all(1,ihelp), 1, HS(1,i_g), 1 )          ! HS(:,i_g)   <- Sk_all(:,ihelp)
        HSDIAG(i_g)=Hamilt(i_g,ihelp)                               ! HSDIAG(i_g) <- Hamilt(i_g,i_g)
     END DO
     CALL STOP_TIMER(time_distiouter)
     DTIMDIST = DTIMDIST + READ_CPU_TIME(time_distiouter)
     DTIMDISTw = DTIMDISTw + READ_wall_TIME(time_distiouter)
  END DO
  
  CALL STOP_TIMER(time_loop260)
  
  call end_albl('HAM')

  deallocate( Pl )
  deallocate( SLEN)
  deallocate( KxK )
  deallocate( Sk_ints )
  deallocate( PPLX, PPLY, TMP1,TMP2, X)
  deallocate( help_x, help_cos, help_sin, help_exp )
  deallocate( j1x, tmp_y, PHASE )

  if (Qprint) then
     WRITE (6,'(a,2f12.1)') 'Time for al,bl    (hamilt, cpu/wall) : ', time_albl_tot, time_albl_totw
     WRITE (6,'(a,2f12.1)') 'Time for legendre (hamilt, cpu/wall) : ', DTIMLG, DTIMLGw
     WRITE (6,'(a,2f12.1)') 'Time for phase    (hamilt, cpu/wall) : ', DTIMPH, DTIMPHw
     WRITE (6,'(a,2f12.1)') 'Time for us       (hamilt, cpu/wall) : ', DTIMUS, DTIMUSw
     WRITE (6,'(a,2f12.1)') 'Time for overlaps (hamilt, cpu/wall) : ', DTIMMA, DTIMMAw
     WRITE (6,'(a,2f12.1)') 'Time for distrib  (hamilt, cpu/wall) : ', DTIMDIST, DTIMDISTw
     WRITE (6,'(a,2f12.1)') 'Time for warpin   (hamilt, cpu/wall) : ', DTIMH, DTIMHw
     WRITE (6,'(a,2f12.1)') 'Time sum iouter   (hamilt, cpu/wall) : ', READ_CPU_TIME(time_loop260), READ_wall_TIME(time_loop260)

     WRITE (6,'(a,i8)') ' number of local orbitals, nlo (hamilt) ', NLO
  endif
  IF (NLO .NE. 0) THEN

     ! This is really confusing and please be careful:
     !   In lapw1 we are using   rotloc * BR2 * rotij *(i1+j1/N1,i2+j2/N2,i3+j3/N3)
     !   more precisely          rotloc * BR2 * rotij * BR2^T^{-1}*T*BR2 * (i1+j1/N1,i2+j2/N2,i3+j3/N3)
     !   but in dmft1 and dmft2 we use
     !                           crotloc * BR1 * rotij *(i1+j1/N1,i2+j2/N2,i3+j3/N3)
     !   which meads that rotloc rotates in cartesian coordinates of the conventional unit cell.
     !   while crotloc rotates in cartesian coordinates of the primitive unit cell.
     DO jatom=1,nat
        lfirst = 1+sum(mult(1:jatom-1))
        DO mu=1,mult(jatom)
           latom = lfirst + mu -1
           rotloc_x_BR2_x_rotij(:,:,latom) = matmul(ROTLOC(:,:,jatom), ROTIJ(:,:,latom))
        ENDDO
     ENDDO
     
     lomaxlomax = (lomax+1)*(lomax+1)
     allocate( C3UP(nloat), C13UP(nloat), C3(nloat), C13(nloat) )
     allocate( iind(lomaxlomax,nloat,multmax) )
     allocate( PHSC(HSrows), YL(0:lomaxlomax-1,HSrows,multmax) )
     if (Qprint) then
        write(6,'(a20,f12.1,a,a,3i6)') 'allocate YL  ',dble(lomaxlomax*HSrows*multmax*16)/(1024.0*1024.0),' MB ','         dimensions',lomaxlomax-1,HSrows,multmax
        write(6,'(a20,f12.1,a,a,i6)') 'allocate phsc',dble(HSrows*16)/(1024.0*1024.0),' MB ','         dimensions',HSrows
     endif
     CALL init_albl(lomax,'HAM')
     !------ calculate LO-Overlap-matrix elements ---------------------------
     CALL START_TIMER(time_lo)
     
     deallocate( Hamilt,Sk_all)
     allocate( Hamilt(1:HSrows, 1:NLO) )
     allocate( Sk_all(1:HSrows, 1:NLO) )
     Hamilt(:,:)=0.d0
     Sk_all(:,:)=0.d0
     
     DO jatom=1,nat
        lfirst = 1+sum(mult(1:jatom-1))
        ! Again precomputing alm,blm to make chi continuous accross MT-sphere
        call make_albl(jatom,Rmt(jatom),lomax,NV,'HAM',1,0)
        ! index for local orbital is precomputed for simpler lookup
        iind(:,:,:)=0
        i=nlov(jatom)
        DO l=0,lomax
           do jlo_i=1,ilo(l,jatom)
              do mu=1,mult(jatom)
                 DO m=-l,l
                    i=i+1
                    iind(l*(l+1)+m+1,jlo_i,mu)=i
                 END DO
              end do
           end do
        END DO
        ! Precalculate spherical-harmonics for Alm, Blm
        ! We compute Ylm's only up to NV (at the moment missing local orbitals tems)
        DO mu=1,mult(jatom)
           latom = lfirst + mu -1
           DO j_g=1,NV
              !-----  rotij transformation from the first atom to an equivalent atom 
              !-----  Here it has already been transformed to cartesian coordinates by BR2
              !-----  Rotloc rotates to local coordinate system from struc file
              ROTV2 = matmul(rotloc_x_BR2_x_rotij(:,:,latom),Kn(:,j_g))
              CALL YLM(ROTV2,LOMAX,Yl(0,j_g,mu))
           ENDDO
        ENDDO
        pi4r2v = 4.d0*pi*Rmt(jatom)**2/sqrt(Vol)
        pi4r2v2 = pi4r2v**2
        DO l=0,lomax
           do jlo_i=1,ilo(l,jatom)
              !        C1,C2,C3 are needed for the Overlap-matrix update
              !        C11,C12,C13 are needed for the Hamilton-matrix update
              !    The local orbital, which vanishes at the MT-boundary is called u_{real}. Then we have:
              !    u_{real} = alo*u + blo*dotu + clo*u_lo
              !
              !    (-nabla^2+V_{sym})|u_{real}> = alo*el*u + blo*(el*dotu + u) + clo*elo*u_lo = (alo*el+blo)*u + blo*el*dotu + clo*elo*u_lo
              !
              !    C1  = <u_{real}|u>    = alo + clo*<u|u_lo>
              !    C2  = <u_{real}|dotu> = blo*<dotu|dotu> + clo*<dotu|u_lo>
              !    C3  = <u_{real}|u_lo> = clo + blo*<dotu|u_lo> + alo*<u|u_lo>
              !    C11 = 0.5*(<u_{real}|H_{sym}|u>+<u|H_{sym}|u_{real}>) = alo*el + 0.5*( blo + clo*<u|u_lo>*(el+elo))
              !    C12 = 0.5*(<u_{real}|H_{sym}|u>+<u|H_{sym}|u_{real}>) = blo*<dotu|dotu>*el + 0.5*(alo+clo*<u|u_lo>+clo*<dotu|u_lo>*(el+elo))
              !    C13 = 0.5*(<u_{real}|H_{sym}|u_lo>+<u_lo|H|u_{real}>) = clo*elo + 0.5*(blo*<u|u_lo> + (el+elo)*(alo*<u|u_lo>+blo*<dotu|u_lo>))
              !
              C1  = ALO(l,jlo_i,jatom) + pi12lo(l,jlo_i,jatom)*CLO(l,jlo_i,jatom)                       
              C2  = BLO(l,jlo_i,jatom)*pei(l,jatom) + pe12lo(l,jlo_i,jatom)*CLO(l,jlo_i,jatom)
              C11 = ALO(l,jlo_i,jatom)*E(l,jatom) + 0.5d0*(BLO(l,jlo_i,jatom) + CLO(l,jlo_i,jatom)*pi12lo(l,jlo_i,jatom)*(Elo(l,jlo_i,jatom) + E(l,jatom)))
              C12 = BLO(l,jlo_i,jatom)*pei(l,jatom)*E(l,jatom) + 0.5d0*(ALO(l,jlo_i,jatom) + CLO(l,jlo_i,jatom)*pi12lo(l,jlo_i,jatom) + CLO(l,jlo_i,jatom)*pe12lo(l,jlo_i,jatom)*(Elo(l,jlo_i,jatom) + E(l,jatom)))
              do jlop=1,jlo_i
                 C3(jlop)  = CLO(l,jlo_i,jatom)*pilolo(l,jlop,jlo_i,jatom) + pe12lo(l,jlop,jatom)*BLO(l,jlo_i,jatom) + ALO(l,jlo_i,jatom)*pi12lo(l,jlop,jatom)
                 C13(jlop) = CLO(l,jlo_i,jatom)*0.5d0*(Elo(l,jlo_i,jatom)+Elo(l,jlop,jatom))*pilolo(l,jlop,jlo_i,jatom) + 0.5d0*BLO(l,jlo_i,jatom)*pi12lo(l,jlop,jatom) + 0.5d0*(Elo(l,jlop,jatom) + E(l,jatom))*(ALO(l,jlo_i,jatom)*pi12lo(l,jlop,jatom)+BLO(l,jlo_i,jatom)*pe12lo(l,jlop,jatom))
              enddo
              !
              ! Precalculate surface kinetic energy factor
              ! akinlo = 1/2 * R^2 * du_{real}/dr(r=R)
              !
              akinlo = 0.5d0*Rmt(jatom)**2*(ALO(l,jlo_i,jatom)*du(l,jatom) + BLO(l,jlo_i,jatom)*dup(l,jatom) + CLO(l,jlo_i,jatom)*dplo(l,jlo_i,jatom))
              
              DO mu_i=1,mult(jatom)
                 DO m_i=-l,l
                    inlo = iind(l*(l+1)+m_i+1,jlo_i,mu_i)
                    ihelp = inlo
                    i_g = NV + inlo
                    
                    DO mup=1,mult(jatom)
                       latom = lfirst + mup -1
                       !  Above we calculated Ylm's for j<NV, i.e., only for reciprocal vectors in the absence of local orbitals
                       !  Here we add the j=i_g term, so that all j<=i_g terms are available to use.
                       ROTV2 = matmul(rotloc_x_BR2_x_rotij(:,:,latom), Kn(:,i_g) )
                       CALL YLM(ROTV2,l,Yl(0,i_g,mup))
                       !        determine phase factor
                       do j_g=1,i_g
                          ARGX = 2.d0*PI*dot_product(POS(:,latom), KZZ(:,j_g)-KZZ(:,i_g) )
                          PHSC(j_g) = dcmplx(dcos(argx),dsin(argx))
                       enddo
                       DO mp=-l,l
                          lmp = l*(l+1) + mp
                          CLYLN = conjg(Yl(lmp,i_g,mup))*PI4R2V2
                          C1UP = C1*CLYLN
                          C2UP = C2*CLYLN
                          do jlop=1,jlo_i
                             C3UP(jlop) = C3(jlop)*CLYLN
                          enddo
                          C11UP = C11*CLYLN
                          C12UP = C12*CLYLN
                          do jlop=1,jlo_i
                             C13UP(jlop) = C13(jlop)*CLYLN
                          enddo
                          akinloup = akinlo*CLYLN
                          C11UP_kinlo = C11UP + akinloup * u(l,jatom)
                          C12UP_kinlo = C12UP + akinloup * up(l,jatom)
                          ! This is the cross tem <chi_\vK| \chi_{local}>
                          DO j_g=1,NV            
                             cYp = Yl(lmp,j_g,mup) * PHSC(j_g)
                             AcYp = AL_r(j_g,l) * cYp
                             BcYp = BL_r(j_g,L) * cYp
                             Sk_all(j_g,ihelp) = Sk_all(j_g,ihelp) + tpcast( AcYp *  C1UP        + BcYp * C2UP )
                             Hamilt(j_g,ihelp) = Hamilt(j_g,ihelp) + tpcast( AcYp *  C11UP_kinlo + BcYp * C12UP_kinlo )
                          enddo
                          
                          do jlo_j=1,jlo_i  ! This is the final term <chi_{local}|chi_{local}>
                             do mu_j=1,mult(jatom)    ! The following loops l,jlo_j,mu_j,m_j give index j for the local orbital
                                do m_j=-l,l
                                   jnlo = iind(l*(l+1)+m_j+1,jlo_j,mu_j)
                                   j_g = NV + jnlo
                                   if (j_g.gt.i_g) CYCLE      ! Since H is hermitian, we just compute j<=i
                                   cYp = Yl(lmp,j_g,mup) * PHSC(j_g)
                                   AcYp = ALO(L,jlo_j,jatom) * cYp
                                   BcYp = BLO(L,jlo_j,jatom) * cYp
                                   CcYp = CLO(L,jlo_j,jatom) * cYp
                                   Sk_all(j_g,ihelp) = Sk_all(j_g,ihelp) + tpcast( AcYp * C1UP  + BcYp * C2UP  + CcYp * C3UP(jlo_j) )
                                   Hamilt(j_g,ihelp) = Hamilt(j_g,ihelp) + tpcast( AcYp * C11UP + BcYp * C12UP + CcYp * C13UP(jlo_j))
                                enddo
                             enddo
                          enddo   ! jlo_j loop =>
                       enddo      ! mp   loop  => set of local orbitals
                    enddo         ! mup  loop  =>
                 enddo     ! m_i loop  =>
              enddo        ! mu_i      =>  set of local orbitals
           enddo           ! jlo_i     =>
        enddo   ! l
     enddo      ! jatom

     !***** Finally store the results for local orbitals
     DO i_g=1,NLO
        do j_g=1,NV+i_g-1
           HS(NV+i_g,j_g) = Hamilt(j_g,i_g)
        enddo
        do j_g=1,NV+i_g
           HS(j_g,NV+i_g) = Sk_all(j_g,i_g)
        enddo
        HSDIAG(NV+i_g)=Hamilt(NV+i_g,i_g)
     END DO
     
     CALL STOP_TIMER(time_lo)
     if (Qprint) WRITE (6,'(a,2f12.1)') 'Time for los      (hamilt, cpu/wall) : ', READ_CPU_TIME(time_lo),READ_wall_TIME(time_lo)
     deallocate( PHSC, YL )
     call end_albl('HAM')
     deallocate( C3UP, C13UP, C3, C13 )
     deallocate( iind )
  ENDIF

  deallocate( hamilt,Sk_all)
  CALL STOP_TIMER(time_hamilt)
  
  RETURN
END SUBROUTINE HAMILTON
