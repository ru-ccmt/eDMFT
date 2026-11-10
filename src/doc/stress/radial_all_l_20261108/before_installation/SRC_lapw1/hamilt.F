      SUBROUTINE HAMILT(LTOP,NAT,NV,U,UP,DU,DUP,MULT,R,POS,V,VOLUME,CHN,E,INS)
!
#ifdef Parallel
      use matrices, only : H, S,HSROWS, KZZ, XK, YK, ZK, RK
#else
      use matrices, only : HS, HSDIAG, HSROWS, KZZ, XK, YK, ZK, RK
#endif
      use loabc, only : ALO, BLO, CLO, ELO, PLO, DPLO, PELO, &
                        DPELO, PEILO, PI12LO, PE12LO,pilolo,nloat_new,LOSECDER
      use parallel, only : ICTXTALL, DESCHS,DESCZ, NPE, MYID,  &
                           NPROW, NPCOL, MYROWHS, MYCOLHS, &
                           BLOCKSIZE, LDHS, LCOLHS, LCOLZ, &
                           INIT_PARALLEL, BARRIER, CLOSE_PARALLEL
      use lapw_timer, only : time_albl, time_loop260, time_phase,   &
                        time_legendre, time_step_function, time_h,   &
        	       time_us, time_loop210, time_loop230,   &
      		       time_loop240, time_overlap, time_lo,   &
      		       START_TIMER, STOP_TIMER, READ_CPU_TIME, READ_wall_TIME, &
                       init_all_timer,time_distiouter
      use lolog, only : nlo, loor, lapw, ilo
      use rotmat, only: ROTIJ, ROTLOC
      use struk, only : NDF,multmax
      use coml, only   : mgga
      use albl
      use stdmath
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Parameters
!
      DOUBLE PRECISION   ONE, TOL, UTOP, ZERO, HALF
      DOUBLE PRECISION   TWO, THREE, FOUR
      PARAMETER          (ONE=1.0D+0, TOL=1.0D-8)
      PARAMETER          (TWO=2.0D+0, THREE=3.0D0, FOUR=4.0D0)
      PARAMETER          (HALF=0.5D0)
      PARAMETER          (UTOP=0.9999999D+0, ZERO=0.0D+0)
      DOUBLE COMPLEX         CIMAG
      PARAMETER          (CIMAG = (0.0D0,1.0D0))
      DOUBLE COMPLEX         CZERO
      PARAMETER          (CZERO=(0.0D0,0.0D0))
      DOUBLE COMPLEX     CONE
      PARAMETER          (CONE=(1.0d0,0.0d0))
!
!        Arguments
!
      INTEGER            INS, LTOP
      INTEGER, INTENT(IN) :: NAT, NV
      INTEGER            MULT(NAT)
      DOUBLE PRECISION   VOLUME
      DOUBLE PRECISION, INTENT(IN) :: R(NAT)
      DOUBLE PRECISION   CHN(0:LMAX-1,NAT), DU(0:LMAX-1,NAT)
      DOUBLE PRECISION   DUP(0:LMAX-1,NAT), E(0:LMAX-1,NAT)
      DOUBLE PRECISION, INTENT(IN) ::  POS(3,NDF), V(NAT)
      DOUBLE PRECISION   POSTwoPi(3,NDF)
      DOUBLE PRECISION   U(0:LMAX-1,NAT), UP(0:LMAX-1,NAT)

!        YL(lm,:,i) - temporary storage of spherical harmonics for
!                     local orbitals
!
      DOUBLE COMPLEX, allocatable :: YL(:,:,:)
!
!        Locals
!
      DOUBLE PRECISION   V3
      INTEGER            jlo,jlop
      DOUBLE PRECISION   atmp1, tmp3, stmp, htmp, btmp1
      DOUBLE PRECISION   akinlor,akinlorup,akinloiup,akinbr,akinbi
      INTEGER            I, INDATM, indatm0, J, JEQ, JEQO, JNEQ, K, L
      INTEGER            L0M0, MO1, MO2, N, NNLO, IATM, IHELP, IOUTER
      INTEGER            IMIN, IMAX, II
      DOUBLE PRECISION   ARGX, ARGY, ARGZ, ATMP, BESR, BTMP
      DOUBLE PRECISION   ARGXX, ARGYY, ARGZZ
      DOUBLE PRECISION   C1R, C1I, C2R, C2I, C3I
      DOUBLE PRECISION   C11R, C11I, C12R, C12I, C13I
      DOUBLE PRECISION   C1RUP, C1IUP, C2RUP, C2IUP
      DOUBLE PRECISION   C11RUP,C11IUP,C12RUP,C12IUP
      DOUBLE PRECISION   C6R, C6I, C7R, C7I, C8R, C8I
      DOUBLE PRECISION   CABR, CABI, CABCR, CABCI, CLR, CLI
      DOUBLE PRECISION   CLYLR, CLYLI, CLYLNR, CLYLNI
      DOUBLE PRECISION   DLL1, PI4R2V, DEN
      DOUBLE PRECISION   DTIME1, DTIME2, DTIMH, DTIMLG, DTIMMA, DTIMPH,DTIMDIST
      DOUBLE PRECISION   DTIME3, dtime240, dtime241, dtime242
      DOUBLE PRECISION   dtime260, dtime211, dtime231, dtime210,dtime230,time_albl_tot
      DOUBLE PRECISION   time_albl_totw, DTIMLGw, DTIMMAw, DTIMPHw,DTIMDISTw,DTIMUSw
      DOUBLE PRECISION   DTIMS, DTIMUS, PI, Q, RKN, TMP,  VIFPR4, rk_tmp1

      DOUBLE PRECISION, allocatable :: C3R(:),C13R(:),C3RUP(:), C3IUP(:),C13RUP(:),C13IUP(:)   
      DOUBLE PRECISION, allocatable ::   rk_tmp(:)
      DOUBLE PRECISION, allocatable ::   PPLX(:), PPLY(:)
      DOUBLE PRECISION, allocatable ::   SLEN(:,:), TMP1(:), TMP2(:)
      DOUBLE PRECISION, allocatable ::   X2(:,:), XL(:),x(:)
      DOUBLE PRECISION, allocatable ::   P(:,:,:)
      DOUBLE PRECISION, dimension(3) ::  FK, ROTV1, ROTV2, VEC
      DOUBLE PRECISION, allocatable ::   xk_tmp_c(:),yk_tmp_c(:),zk_tmp_c(:)
      DOUBLE PRECISION, allocatable ::   xk_tmp_r(:),yk_tmp_r(:),zk_tmp_r(:)
#if defined (INTEL_VML)
!_COMPLEX      DOUBLE COMPLEX, allocatable ::     help_exp(:)
!_REAL      DOUBLE PRECISION ::     help_exp
      DOUBLE PRECISION, allocatable ::   help1(:), help_cos(:)
      DOUBLE PRECISION, allocatable ::   help_sin(:), tmp_y(:)
      DOUBLE PRECISION, allocatable ::   tmp_ytmp(:)
#else
      DOUBLE PRECISION              :: help_sin, help1, help_cos,tmp_y, tmp_ytmp
      DOUBLE COMPLEX                :: help_exp
#endif
!_REAL      DOUBLE PRECISION ::   PHS
!_REAL      DOUBLE PRECISION, allocatable ::   PHASE(:),spanelus(:,:)
!_COMPLEX      DOUBLE COMPLEX ::     PHS
!_COMPLEX      DOUBLE COMPLEX, allocatable ::     PHASE(:),spanelus(:,:)
      DOUBLE COMPLEX, allocatable ::         PHSC(:)
!_REAL      DOUBLE PRECISION, allocatable ::  HPANEL(:,:)
!_REAL      DOUBLE PRECISION, allocatable ::  SPANEL(:,:)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: HPANEL(:, :)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: SPANEL(:, :)
            integer, allocatable:: kzz_tmp(:,:), kzz_ihelp(:,:)
            integer nv1,nv2
!
!        External Functions
!
!_REAL      DOUBLE PRECISION   USTPHX, WARPIN, vtau_warp_func
!_COMPLEX      DOUBLE COMPLEX         USTPHX, WARPIN, vtau_warp_func
      EXTERNAL           USTPHX, WARPIN
!
!        External Subroutines
!
      EXTERNAL           CPUTIM, ROTATE, YLM
!
!        Intrinsic Functions
!
      INTRINSIC          ABS, ATAN, DCONJG, COS, DBLE, DCMPLX, EXP, SIN
      INTRINSIC          SQRT

      integer i_g,j_g,hsdim_r,hsdim_c,i_col,ipr,ipc,ii_g
      integer global2local, jnnlo
      integer, allocatable :: j_end(:), j_end_1(:)

!     For precomputation of bessel functions and derivatives
      double precision, allocatable :: FJ1(:,:), DFJ1(:,:)
      double precision, allocatable :: FJ2(:,:), DFJ2(:,:)
      integer ltop2, imin_old, imax_old, jend_rk! , jend_rk_all
      double precision  rneq_old
      logical recalc1, recalc2

#ifdef Barrier
      CALL BARRIER
#endif

!
      allocate( XL(LMAX-2) )

#ifdef Parallel
      HSdim_r=LDHS
      HSdim_c=LCOLHS
#else
      HSdim_c=HSROWS
      HSdim_r=HSROWS
#endif

      allocate( rk_tmp(HSdim_r) )
      allocate( kzz_TMP(3,HSdim_r) )
      allocate(xk_tmp_r(HSdim_r),yk_tmp_r(HSdim_r),zk_tmp_r(HSdim_r))

!_REAL     write(6,'(a20,f12.1,a,a,2i6)') 'allocate spanel',HSdim_r*blocksize*8.0/(1024.0*1024.0),' MB ',&
!_REAL           '         dimensions',HSdim_r,blocksize
!_REAL     write(6,'(a20,f12.1,a,a,2i6)') 'allocate hpanel',HSdim_r*blocksize*8.0/(1024.0*1024.0),' MB ',&
!_REAL           '         dimensions',HSdim_r,blocksize
!_COMPLEX  write(6,'(a20,f12.1,a,a,2i6)') 'allocate spanel',HSdim_r*blocksize*16.0/(1024.0*1024.0),' MB ',&
!_COMPLEX        '         dimensions',HSdim_r,blocksize
!_COMPLEX  write(6,'(a20,f12.1,a,a,2i6)') 'allocate hpanel',HSdim_r*blocksize*16.0/(1024.0*1024.0),' MB ',&
!_COMPLEX        '         dimensions',HSdim_r,blocksize
!_REAL      write(6,'(a20,f12.1,a,a,2i6)') 'allocate spanelus',HSdim_r*blocksize*8.0/(1024.0*1024.0),' MB ',&
!_REAL           '         dimensions',HSdim_r,blocksize
!_COMPLEX   write(6,'(a20,f12.1,a,a,2i6)') 'allocate spanelus',HSdim_r*blocksize*16.0/(1024.0*1024.0),' MB ',&
!_COMPLEX           '         dimensions',HSdim_r,blocksize
      write(6,'(a20,f12.1,a,a,2i6)') 'allocate slen',HSdim_r*blocksize*8.0/(1024.0*1024.0),' MB ',&
           '         dimensions',HSdim_r,blocksize
      write(6,'(a20,f12.1,a,a,2i6)') 'allocate x2',HSdim_r*blocksize*8.0/(1024.0*1024.0),' MB ',&
           '         dimensions',HSdim_r,blocksize
      write(6,'(a20,f12.1,a,a,3i6)') 'allocate legendre',dble((HSdim_r)*LMAX*blocksize*8)/(1024.0*1024.0),' MB ',&
           '         dimensions',HSdim_r,lmax,blocksize
      write(6,'(a20,f12.1,a,a,2i6)') 'allocate al,bl (row)',HSdim_r*ltop*2.0*2*8.0/(1024.0*1024.0),' MB ',&
           '         dimensions',HSdim_r,ltop
      write(6,'(a20,f12.1,a,a,2i6)') 'allocate al,bl (col)',blocksize*ltop*2.0*2*8.0/(1024.0*1024.0),' MB ',&
           '         dimensions',blocksize,ltop
      write(6,'(a20,f12.1,a,a,3i6)') 'allocate YL',dble((LOMAX+1)*(LOMAX+1)*HSdim_r*multmax*16)/(1024.0*1024.0),' MB ',&
           '         dimensions',(LOMAX+1)*(LOMAX+1)-1,HSdim_r,multmax

      PI = FOUR*ATAN(ONE)
      DO I = 1, LMAX - 2
         XL(I) = DBLE(I)/DBLE(I+1)
      enddo
      DTIMPH = ZERO
      DTIMUS = ZERO
      DTIMLG = ZERO
      DTIMMA = ZERO
      DTIMS = ZERO
      DTIMH = ZERO
      DTIMDIST=ZERO
      time_albl_tot=zero
      DTIMPHw = ZERO
      DTIMUSw = ZERO
      DTIMLGw = ZERO
      DTIMMAw = ZERO
      DTIMDISTw=ZERO
      time_albl_totw=zero
!
!        compute Overlap-matrix
!
      CALL START_TIMER(time_loop260)

      j=0
      DO J_g = 1,NV
#ifdef Parallel
         ipr=mod((j_g-1)/blocksize,nprow)
         if (ipr.ne.myrowhs) cycle
         j=j+1
#else
         j=j_g
#endif
         kzz_tmp(1:3,j)=KZZ(1:3,J_g) 
         xk_tmp_r(j)=XK(J_g)
         yk_tmp_r(j)=YK(J_g)
         zk_tmp_r(j)=ZK(J_g)
         rk_tmp(j)=RK(J_g)
      enddo
      jend_rk=j
!print*,myid,jend_rk, hsdim_r

      POSTwoPi  = POS*TWO*PI

!$omp parallel private(atmp,atmp1,btmp,btmp1,dll1,j,tmp3,stmp,htmp, &
!$omp& tmp1,tmp2,ihelp, phase, &
!$omp& jeq, indatm,argxx,argx, help1, help_cos, help_sin, help_exp, &
!$omp& v3,tmp_y,tmp_ytmp,besr,ipr,i,i_g,j_g,fk, ipc, imin, imax, &
!$omp& j_end_1,j_end, kzz_ihelp,xk_tmp_c,yk_tmp_c,zk_tmp_c, &
!$omp& rk_tmp1,x2,slen,x,den,pplx,pply,p,k,spanel,spanelus,hpanel, &
!$omp& jneq,INDATM0,vifpr4,recalc1,recalc2,imin_old,imax_old,rneq_old,FJ1,FJ2,DFJ1,DFJ2) &
!$omp& reduction(max:DTIMPH,DTIMS,DTIMH,DTIMUS,DTIMMA,time_albl_tot,DTIMLG, &
!$omp& DTIMPHw,DTIMUSw,DTIMMAw,time_albl_totw,DTIMLGw,DTIMDIST,DTIMDISTw) 

!!!$omp& jneq,INDATM0,vifpr4,imin_old,imax_old,rneq_old,jend_rk,FJ1,FJ2,DFJ1,DFJ2) &

      allocate( HPANEL( HSdim_r, BLOCKSIZE) )
      allocate( SPANEL( HSdim_r, BLOCKSIZE) )
      allocate( xk_tmp_c(blocksize),yk_tmp_c(blocksize),zk_tmp_c(blocksize) )
      allocate( kzz_ihelp(3,blocksize) )
      allocate( j_end(blocksize),j_end_1(blocksize) )
      allocate( PHASE(HSdim_r) )
      allocate( spanelus(HSdim_r,blocksize))
      allocate( SLEN(HSdim_r,blocksize), TMP1(HSdim_r) )
      allocate( TMP2(HSdim_r) )
      allocate( X2(HSdim_r,blocksize), X(HSdim_r) )
      allocate( PPLX(HSdim_r), PPLY(HSdim_r) )
      allocate( P(HSdim_r,0:LTOP-1,blocksize) )

#if defined (INTEL_VML)
      allocate( help1(HSdim_r), help_cos(HSdim_r), help_sin(HSdim_r) )
      allocate( tmp_ytmp(HSdim_r) )
      allocate( tmp_y(HSdim_r) )
!_COMPLEX      allocate( help_exp(HSdim_r) )
#endif

!     Allocate enough values for bessel functions and derivatives
      ltop2=max(ltop-1,2)       ! Warning, make_albl is called with ltop-1 (not ltop)
      allocate ( FJ1(0:ltop2,HSdim_r) , FJ2 (0:ltop2,Blocksize) )
      allocate ( DFJ1(0:ltop2,HSdim_r), DFJ2(0:ltop2,Blocksize) )
!
!     setup counters
      imin_old=-1 ; imax_old=-1
      rneq_old=-1

      CALL init_albl(ltop-1,'HAM')
#if defined (_OPENMP)
      if (DTIMPH .lt. -HUGE(DTIMPH)) then
         DTIMPH = ZERO
         DTIMUS = ZERO
         DTIMLG = ZERO
         DTIMMA = ZERO
         DTIMS = ZERO
         DTIMH = ZERO
         DTIMDIST=ZERO
         time_albl_tot=zero
         DTIMPHw = ZERO
         DTIMUSw = ZERO
         DTIMLGw = ZERO
         DTIMMAw = ZERO
         DTIMDISTw=ZERO
         time_albl_totw=zero
      endif
#endif

!$omp do schedule(dynamic,1)
      iouter_loop: DO IOUTER = 0, NV/(BLOCKSIZE)

! To properly count the CPU and wallclock time with OpenMP we need to use
! the max reduction, however we need the values to be initialized to 0 instead of -inf.
! Proper solution is to bake our own with omp declare reduction, but this is currently
! impossible due to missing ifort support (2018).

#ifdef Parallel
         ipc=mod(iouter,npcol)
         if (mycolhs.ne.ipc) cycle
#endif
         IMIN = IOUTER*BLOCKSIZE + 1
         IMAX = MIN(NV, (IOUTER + 1)*BLOCKSIZE)
!
!        precompute Legendre-Polynomials
!
         CALL START_TIMER(time_legendre)

         IHELP = 0
!!DEC$ PARALLEL
         legrende: DO I_g = IMIN, IMAX
            IHELP = IHELP + 1

#ifdef Parallel
            j=0
            DO J_g = 1, I_g-1
               ipr=mod((j_g-1)/blocksize,nprow)
               if (ipr.ne.myrowhs) cycle
               j=j+1
            enddo
            j_end_1(ihelp)=j
            ipr=mod((i_g-1)/blocksize,nprow)
            if (ipr.eq.myrowhs) j=j+1
            j_end(ihelp)=j
#else
            j_end(ihelp)=i_g
            j_end_1(ihelp)=i_g-1
#endif

            kzz_ihelp(1:3,ihelp)=KZZ(1:3,I_g) 
            xk_tmp_c(ihelp)=xk(i_g)
            yk_tmp_c(ihelp)=yk(i_g)
            zk_tmp_c(ihelp)=zk(i_g)
            rk_tmp1=RK(I_g)

            DO J = 1,j_end(ihelp)
               X2(J,ihelp) =  XK_tmp_c(ihelp)*XK_tmp_r(J) +  &
                              YK_tmp_c(ihelp)*YK_tmp_r(J) +  &
                              ZK_tmp_c(ihelp)*ZK_tmp_r(J)
               SLEN(J,ihelp) = SQRT((XK_tmp_r(J)-XK_tmp_c(ihelp))*(XK_tmp_r(J)-XK_tmp_c(ihelp))+ &
                                    (YK_tmp_r(J)-YK_tmp_c(ihelp))*(YK_tmp_r(J)-YK_tmp_c(ihelp))+ &
                                    (ZK_tmp_r(J)-ZK_tmp_c(ihelp))*(ZK_tmp_r(J)-ZK_tmp_c(ihelp)))
               x(j)=X2(J,ihelp)
               den=rk_tmp1*rk_tmp(j)
               IF (DEN .GT. TOL) X(J) = X(J)/DEN               
               PPLX(J) = X(J)*X(J)
               PPLY(J) = PPLX(J) - ONE
               P(J,0,ihelp) = ONE
               P(J,1,ihelp) = X(J)
            enddo

            DO K = 2, LTOP - 2
               DO J = 1, j_end(ihelp)
                  IF (ABS(X(J)) .LT. UTOP) THEN
                     P(J,K,ihelp) = PPLX(J) + XL(K-1)*PPLY(J)
                     PPLX(J) = X(J)*P(J,K,ihelp)
                     PPLY(J) = PPLX(J) - P(J,K-1,ihelp)
                     P(J,LTOP-1,ihelp) = PPLX(J) + XL(LTOP-2)*PPLY(J)
                  ELSE
                     P(J,K,ihelp) = X(J)*P(J,K-1,ihelp)
                  ENDIF
               enddo
            enddo
            
            DO J = 1, j_end(ihelp)
               IF (ABS(X(J)) .LT. UTOP) THEN
                  P(J,LTOP-1,ihelp) = PPLX(J) + XL(LTOP-2)*PPLY(J)
               ELSE
                  P(J,LTOP-1,ihelp) = X(J)*P(J,LTOP-2,ihelp)
               ENDIF
            enddo

         enddo legrende

         CALL STOP_TIMER(time_legendre)
         DTIMLG = DTIMLG + READ_CPU_TIME(time_legendre)
         DTIMLGw = DTIMLGw + READ_wall_TIME(time_legendre)

!_REAL          SPANEL(:,:) = ZERO
!_COMPLEX       SPANEL(:,:) = CZERO            
!_REAL          spanelus(:,:)=ZERO
!_COMPLEX       spanelus(:,:)=CZERO
!_REAL          HPANEL(:,:) = ZERO
!_COMPLEX       HPANEL(:,:) = CZERO

         jneq_loop: DO JNEQ = 1, NAT

            if (jneq.eq.1) then
               INDATM0 = 0
            else
               indatm0= indatm0+mult(jneq-1)
            endif

            VIFPR4 = FOUR*PI*VOLUME*R(JNEQ)**4

!
!        precompute al(kn) and bl(kn)
!
            CALL START_TIMER(time_albl)

!           New version
!
!           Have imin/imax changed ?
            if( (imin .ne. imin_old) .or. (imax .ne. imax_old) )then
                recalc2=.true.
            else
                recalc2=.false.
            endif
!
!           Has r(jneq) changed ?
            if( abs(r(jneq) - rneq_old) .gt. 1D-10)then
                recalc1=.true.
                recalc2=.true.
            else
                recalc1=.false.
            endif
!
!           New version, where spherical bessels and derivatives 
!           are only recalculated if needed
!if(recalc1) print*,myid,recalc1,jend_rk,HSdim_r,rneq_old
!if(recalc2) print*,myid,recalc2,imin,imin_old,imax,imax_old
            imin_old=imin ; imax_old=imax ; rneq_old=r(jneq)
            call make_alblH(jneq,r(jneq),ltop-1,NV,'HAM',imin,imax, &
                FJ1,FJ2,DFJ1,DFJ2, HSdim_r, BlockSize, &
                recalc1,recalc2,ltop2,rk_tmp, jend_rk) !HSdim_r) !jend_rk)
!            call make_albl(jneq,r(jneq),ltop-1,NV,'HAM',imin,imax)
!if(recalc1)print*, myid,hsdim_r,al_r(1:HSdim_r,0)
            CALL STOP_TIMER(time_albl)
            time_albl_tot=time_albl_tot+READ_CPU_TIME(time_albl)
            time_albl_totw=time_albl_totw+READ_wall_TIME(time_albl)

            IMIN = IOUTER*BLOCKSIZE + 1
            IMAX = MIN(NV, (IOUTER + 1) * BLOCKSIZE)
            IHELP = 0
!            POSTwoPi  = POS*TWO*PI

!DEC$ PARALLEL
            imin1: DO I_g = IMIN, IMAX

               IHELP = IHELP + 1
!
!        calculate phase-factors
!
               CALL START_TIMER(time_phase)

!_REAL         PHASE(1:j_end(ihelp)) = ZERO
!_COMPLEX      PHASE(1:j_end(ihelp)) = CZERO

               DO JEQ = 1, MULT(JNEQ)
            
                  INDATM = INDATM0 + jeq
                  
                  ARGXX = POSTwoPi(1,INDATM)*KZZ_ihelp(1,ihelp) &
                        + POSTwoPi(2,INDATM)*KZZ_ihelp(2,ihelp) &
                        + POSTwoPi(3,INDATM)*KZZ_ihelp(3,ihelp)

#if !defined (INTEL_VML)
!$omp simd private(argx,help1,help_cos,help_exp)
#endif
                  DO J = 1,j_end(ihelp)
                     ARGX = POSTwoPi(1,INDATM)*KZZ_tmp(1,J) - ARGXX
                     ARGX = POSTwoPi(2,INDATM)*KZZ_tmp(2,J) + ARGX
                     ARGX = POSTwoPi(3,INDATM)*KZZ_tmp(3,J) + ARGX
#if defined (INTEL_VML)
                     help1(j)=(ARGX)
                  enddo

!_REAL            call vdcos(j_end(ihelp),help1,help_cos)
!_COMPLEX         call vzcis(j_end(ihelp),help1,help_exp)

! in case of older mk/vml you may miss vzcis (IPO Error: unresolved : vzcis_)
! Either upgrade the libraries (ifort+mkl) or:
! Comment the line below and uncomment the other ones. 
! Uncomment also all other lines above where   help_tmpcos occurs
!!_COMPLEX         call vdcos(j_end(ihelp),help1,help_tmpcos)
!!_COMPLEX         call vdsin(j_end(ihelp),help1,help_tmpsin)
!!_COMPLEX         do j = 1, j_end(ihelp)
!!_COMPLEX            help_exp(j) = dcmplx(help_tmpcos(j),help_tmpsin(j))
!!_COMPLEX         end do
#else
                     help1 = (ARGX)
!_REAL               help_cos = dcos(help1)
!_COMPLEX            help_exp = dcmplx(dcos(help1),dsin(help1))
#endif
#if defined (INTEL_VML)
                  do j = 1,j_end(ihelp)
!_REAL               PHASE(j) = PHASE(j) + help_cos(j)
!_COMPLEX            PHASE(j) = PHASE(j) + help_exp(j)
#else
!_REAL               PHASE(j) = PHASE(j) + help_cos
!_COMPLEX            PHASE(j) = PHASE(j) + help_exp
#endif
                  enddo
               enddo 

               CALL STOP_TIMER(time_phase)
               DTIMPH = DTIMPH + READ_CPU_TIME(time_phase)
               DTIMPHw = DTIMPHw + READ_wall_TIME(time_phase)

               CALL START_TIMER(time_us)
!
!        initialize Overlap-matrixelements with step-function values
!
ins=0    !fixed PB 20.5.08. because alternative (rm case.vns) does not work anymore
               IF (INS .EQ. 0) THEN

                  CALL START_TIMER(time_step_function)
                  v3=THREE*V(JNEQ)

#if !defined (INTEL_VML)
!$omp simd private(help1,tmp_y,help_sin,help_cos,besr)
#endif

                  do j=1,j_end_1(ihelp)
#if defined (INTEL_VML)
                     help1(j) = slen(j,ihelp)*r(jneq)
                     tmp_y(j) = help1(j)*help1(j)*help1(j)
                  end do
                  call vdsincos(j_end_1(ihelp),help1,help_sin,help_cos)
#else
                     help1 = slen(j,ihelp)*r(jneq)
                     tmp_y = help1*help1*help1
                     help_sin = dsin(help1) 
                     help_cos = dcos(help1)
                     tmp_y = 1.0D0 / tmp_y
#endif

#if defined (INTEL_VML)
                  call vdinv(j_end_1(ihelp), tmp_y, tmp_ytmp)
!!LDM                  call dcopy(j_end_1(ihelp), tmp_ytmp, 1, tmp_y, 1)
                  DO J = 1, j_end_1(ihelp)
!!LDM                     BESR = (help1(j)*help_COS(j)-help_SIN(j))*tmp_y(j)
                     BESR = (help1(j)*help_COS(j)-help_SIN(j))*tmp_ytmp(j)
                     spanelus(j,ihelp)=spanelus(j,ihelp)+v3*PHASE(j)*BESR
#else
                     BESR = (help1*help_COS-help_SIN)*tmp_y
                     spanelus(j,ihelp)=spanelus(j,ihelp)+v3*PHASE(j)*BESR
#endif
                  enddo
#ifdef Parallel
                  ipr=mod((i_g-1)/blocksize,nprow)
                  if (ipr.eq.myrowhs) then
                     i=global2local(i_g,nprow,blocksize)                     
                     spanelus(i,ihelp)=spanelus(i,ihelp)+ONE/DBLE(NAT) - V(JNEQ)*DBLE(MULT(JNEQ))
                  endif
#else
                  spanelus(i_g,ihelp)=spanelus(i_g,ihelp)+ONE/DBLE(NAT) - V(JNEQ)*DBLE(MULT(JNEQ))
#endif
                  CALL STOP_TIMER(time_step_function)
                  DTIMS = DTIMS + READ_CPU_TIME(time_step_function)
                  
                  CALL START_TIMER(time_h)

                  if (jneq.eq.nat) then

                     j=0
                     DO J = 1, j_end(ihelp)
                        SPANEL(J,IHELP) = SPANEL(J,IHELP)+spanelus(j,ihelp)
                        HPANEL(J,IHELP) = HPANEL(J,IHELP)+X2(j,ihelp)*spanelus(j,ihelp)
                        HPANEL(J,IHELP) = HPANEL(J,IHELP) + &
                            WARPIN(KZZ_tmp(1,J)-KZZ_ihelp(1,ihelp), &
                                       KZZ_tmp(2,J)-KZZ_ihelp(2,ihelp), &
                                       KZZ_tmp(3,J)-KZZ_ihelp(3,ihelp))
                     enddo
                     if(mgga) then
                        do j = 1, j_end(ihelp)
                           HPANEL(J,IHELP) = HPANEL(J,IHELP) + & ! Add (G+k).(G'+k)*VTAU
                                vtau_warp_func(KZZ_tmp(1,J)- KZZ_ihelp(1,ihelp), &  ! for MGGA calculation 
                                KZZ_tmp(2,J)-KZZ_ihelp(2,ihelp), &               ! KZZ_ihelp = G'
                                KZZ_tmp(3,J)-KZZ_ihelp(3,ihelp)) * &
                                (xk_tmp_r(j)*xk_tmp_c(ihelp) + &                  ! CHECK THIS
                                 yk_tmp_r(j)*yk_tmp_c(ihelp) + & ! seems ok...
                                 zk_tmp_r(j)*zk_tmp_c(ihelp))
                        enddo
                     endif    
                  endif

                  CALL STOP_TIMER(time_h)
                  DTIMH = DTIMH + READ_CPU_TIME(time_h)

! else if (ins .ne. 0)
! Should this branch not better be completely removed? It is hardcoded
! out by ins=0 above
! and array sizes are even incompatible?? (PHASE(HSROWS,NAT) vs.
! PHASE(HSROWS) above!
! jdoumont 5/10/20
       ELSE
                     
                  if (jneq.eq.nat) then
                     j=0
                     DO J_g = 1, I_g
#ifdef Parallel
                        ipr=mod((j_g-1)/blocksize,nprow)
                        if (ipr.ne.myrowhs) cycle
#endif
                        j=j+1                        
                        FK(1) = XK(J_g) - XK(I_g)
                        FK(2) = YK(J_g) - YK(I_g)
                        FK(3) = ZK(J_g) - ZK(I_g)
                        SPANEL(J,IHELP) = USTPHX(FK,PHASE,j,NAT,V,R,MULT)
                        HPANEL(J,IHELP) = X2(j,ihelp)*SPANEL(J,IHELP)
                     enddo
                  endif
                  
       ENDIF
!
               CALL STOP_TIMER(time_us)
               DTIMUS = DTIMUS + READ_CPU_TIME(time_us)
               DTIMUSw = DTIMUSw + READ_wall_TIME(time_us)

               CALL START_TIMER(time_overlap)
!
!        calculate Overlap-matrixelement
!
               TMP1(1:j_end(ihelp)) = ZERO
               TMP2(1:j_end(ihelp)) = ZERO

               DO L = 0, LTOP - 1
                  DLL1 = DBLE(L+L+1)*VIFPR4 
                  if (.not.lapw(l,jneq)) then
                     ATMP = AL_c(L,ihelp) *DLL1
                     atmp1 = ATMP*U(L,JNEQ)*DU(L,JNEQ)*R(JNEQ)**2
!                     DLL1 = DBLE(L+L+1)
!$omp simd private(tmp3,stmp,htmp)
                     DO J = 1, j_end(ihelp)
                        TMP3 = P(J,L,ihelp) !*DLL1
                        STMP = ATMP*AL_r(J,L)*TMP3
                        TMP2(J) = TMP2(J) + STMP
!                        HTMP = STMP*E(L,JNEQ)
                        TMP1(J) = TMP1(J) + STMP*E(L,JNEQ) !+ HTMP
!  surface part of kinetic energy 
                        tmp1(j) = tmp1(j)+TMP3*ATMP1*AL_r(J,L)
                     end do
                  else

!....for lapw all BL lines must be kept
                     ATMP = AL_c(L,ihelp)
                     BTMP = BL_c(L,ihelp)*CHN(L,JNEQ)
                     atmp1 = R(JNEQ)**2* (ATMP*U(L,JNEQ)*DU(L,JNEQ) + &
                          BL_c(L,ihelp)*UP(L,JNEQ)*DU(L,JNEQ))
                     btmp1 = R(JNEQ)**2*(BL_c(L,Ihelp)*UP(L,JNEQ)*DUP(L,JNEQ) + &
                          AL_c(L,ihelp)*U(L,JNEQ)*DUp(L,JNEQ))
!                     DLL1 = DBLE(L+L+1)
                     ATMP = ATMP*DLL1
                     BTMP = BTMP*DLL1
                     ATMP1= ATMP1*DLL1
                     BTMP1= BTMP1*DLL1
!$omp simd private(tmp3,stmp,htmp)
                     DO J = 1, j_end(ihelp)
                        TMP3 = P(J,L,ihelp) !*DLL1
                        STMP = ATMP*AL_r(J,L)
                        STMP = (STMP + BTMP*BL_r(J,L))*TMP3
!                        STMP = STMP*TMP3
                        TMP2(J) = TMP2(J) + STMP
                        HTMP = ATMP*BL_r(J,L)*TMP3
!                        HTMP = HTMP*TMP3
!                        HTMP = HTMP + STMP*E(L,JNEQ)
                        TMP1(J) = TMP1(J) + HTMP + STMP*E(L,JNEQ)
!  surface part of kinetic energy 
                        tmp1(j) = tmp1(j)+&
                             TMP3*(ATMP1*AL_r(J,L)+Btmp1*BL_r(J,L))
                     end do
                  endif
                  
               enddo
!$omp simd
               DO J = 1, j_end(ihelp)
!                  TMP1(J) = TMP1(J)*VIFPR4
!                  TMP2(J) = TMP2(J)*VIFPR4
                  SPANEL(J,IHELP) = SPANEL(J,IHELP) + PHASE(J)*TMP2(J)
                  HPANEL(J,IHELP) = HPANEL(J,IHELP) + PHASE(J)*TMP1(J)
               enddo

               CALL STOP_TIMER(time_overlap)

               DTIMMA = DTIMMA + READ_CPU_TIME(time_overlap)
               DTIMMAw = DTIMMAw + READ_wall_TIME(time_overlap)

             END DO imin1
         enddo jneq_loop

         CALL start_TIMER(time_distiouter)

         IHELP = 0
         IMIN = IOUTER*BLOCKSIZE + 1
         IMAX = MIN(NV, (IOUTER + 1)* BLOCKSIZE)
         
         DO I_g = IMIN, IMAX
            IHELP = IHELP + 1
#ifdef Parallel
            i=global2local(i_g,npcol,blocksize)
!_REAL      CALL DCOPY(HSdim_r, HPANEL(1,IHELP), 1, H(1,I), 1 )
!_REAL      CALL DCOPY(HSdim_r, SPANEL(1,IHELP), 1, S(1,I), 1 )
!_COMPLEX   CALL ZCOPY(HSdim_r, HPANEL(1,IHELP), 1, H(1,I), 1 )
!_COMPLEX   CALL ZCOPY(HSdim_r, SPANEL(1,IHELP), 1, S(1,I), 1 )            
#else
!_REAL      CALL DCOPY(I_g-1, HPANEL(1,IHELP), 1, HS(I_g, 1), HSrows )
!_REAL      CALL DCOPY(I_g, SPANEL(1,IHELP), 1, HS(1,I_g), 1 )
!_COMPLEX   CALL ZCOPY(I_g-1, HPANEL(1,IHELP), 1, HS(I_g, 1), HSrows )
!_COMPLEX   CALL ZCOPY(I_g, SPANEL(1,IHELP), 1, HS(1,I_g), 1 )
            HSDIAG(I_g)=HPANEL(I_g,IHELP)
#endif
         END DO

         CALL STOP_TIMER(time_distiouter)
         DTIMDIST = DTIMDIST + READ_CPU_TIME(time_distiouter)
         DTIMDISTw = DTIMDISTw + READ_wall_TIME(time_distiouter)

        END DO iouter_loop
!$omp end do nowait

      call end_albl('HAM')

      deallocate( PHASE )
      deallocate( spanelus )
      deallocate( j_end, j_end_1)
      deallocate( SLEN, TMP1,TMP2)
      deallocate( X, x2 )
      deallocate( kzz_ihelp )
      deallocate( PPLX, PPLY )
      deallocate( P )
      deallocate ( FJ1 , FJ2,  DFJ1, DFJ2 )
#if defined (_OPENMP)
      deallocate( HPANEL, SPANEL )
#endif

#if defined (INTEL_VML)
      deallocate( help1, help_cos, help_sin )
!_COMPLEX      deallocate( help_exp )
      deallocate( tmp_ytmp )
      deallocate( tmp_y )
#endif
!!$omp end parallel
!$omp single

      CALL STOP_TIMER(time_loop260)

!this is needed ATM since the first alloc was in a parallel section so we can't just reuse it
#if defined (_OPENMP)
      allocate( HPANEL( HSdim_r, BLOCKSIZE) )
      allocate( SPANEL( HSdim_r, BLOCKSIZE) )
#endif

#ifdef Barrier      
      CALL BARRIER
#endif
      WRITE (6,'(a,i8)') ' number of local orbitals, nlo (hamilt) ', NLO

       IF (NLO .NE. 0) THEN

      allocate( PHSC(HSrows) )
      allocate( YL(0:(LOMAX+1)*(LOMAX+1)-1,HSrows,multmax) )
      nloat=nloat_new
      allocate(C3R(nloat),C13R(nloat),C3RUP(nloat),C3IUP(nloat),C13RUP(nloat),C13IUP(nloat))
      write(6,'(a20,f12.1,a,a,3i6)') 'allocate YL  ',dble((LOMAX+1)*(LOMAX+1)*HSrows*multmax*16)/(1024.0*1024.0),' MB ',&
           '         dimensions',(LOMAX+1)*(LOMAX+1)-1,HSrows,multmax
      write(6,'(a20,f12.1,a,a,i6)') 'allocate phsc',dble(HSrows*16)/(1024.0*1024.0),' MB ',&
           '         dimensions',HSrows

      CALL init_albl(lomax,'HAM')
!
!------ calculate LO-Overlap-matrixelements ---------------------------
!
         CALL START_TIMER(time_lo)
         INDATM = 0
         IATM = 0
         NNLO = 0
!
!        initialize matrixelements
!
!_REAL       HPANEL(1:HSdim_r, 1:BLOCKSIZE)=ZERO
!_REAL       SPANEL(1:HSdim_r, 1:BLOCKSIZE)=ZERO
!_COMPLEX       HPANEL(1:HSdim_r, 1:BLOCKSIZE)=CZERO
!_COMPLEX       SPANEL(1:HSdim_r, 1:BLOCKSIZE)=CZERO


         I_g = NV + 1
         IMIN = I_g
         IMAX = MIN(NV + NLO,((IMIN-1)/(BLOCKSIZE*NPCOL)+1)*BLOCKSIZE*NPCOL)

         jneq_lo_loop: DO JNEQ = 1, NAT

!
!        precompute al(kn) and bl(kn)
!
            call make_albl(jneq,r(jneq),lomax,NV,'HAM',1,0)

            DO JEQ = 1, MULT(JNEQ)
               IATM = IATM + 1
!
!        precalculate spherical-harmonics for Alm, Blm
!
               DO J_g = 1, NV
#ifdef Parallel
                  ipr=mod((j_g-1)/blocksize,nprow)
                  if (ipr.ne.myrowhs) cycle
#endif
                  VEC(1) = XK(J_g)
                  VEC(2) = YK(J_g)
                  VEC(3) = ZK(J_g)
                  CALL ROTATE(VEC,ROTIJ(1,1,IATM),ROTV1)
                  CALL ROTATE(ROTV1,ROTLOC(1,1,JNEQ),ROTV2)
                  CALL YLM(ROTV2,LOMAX,YL(0,J_g,JEQ))

               enddo
            enddo

            PI4R2V = FOUR*PI*SQRT(VOLUME)*R(JNEQ)*R(JNEQ)
            CLR = ONE
            CLI = ZERO
            
            l_loop: DO L = 0, LOMAX
               jlo_loop:  do  jlo=1,ilo(l,jneq)

!        precalculate factors of a,b,c for LO's
!        C1,C2,C3 are needed for the Overlap-matrix update
!        C11,C12,C13 are needed for the Hamilton-matrix update
!
                  C1R  = ALO(L,jlo,JNEQ) + &
                         PI12LO(L,jlo,JNEQ)*CLO(L,jlo,JNEQ)
                  C1I  = ZERO
                  C2R  = BLO(L,jlo,JNEQ)*CHN(L,JNEQ) + &
                         PE12LO(L,jlo,JNEQ)*CLO(L,jlo,JNEQ)
                  C2I  = ZERO
                  if (LOSECDER(L,jlo,JNEQ)) then
                     C11R = ALO(L,jlo,JNEQ)*E(L,JNEQ) + &
                            HALF*BLO(L,jlo,JNEQ)  + &
                            E(L,JNEQ)*CLO(L,jlo,JNEQ)*PI12LO(L,jlo,JNEQ) !non orthogonal term
                     C12R = HALF*ALO(L,jlo,JNEQ) + &
                                 E(L,JNEQ)*BLO(L,jlo,JNEQ)*CHN(L,JNEQ) + &
                                 E(L,JNEQ)*CLO(L,jlo,JNEQ)*PE12LO(L,jlo,JNEQ) + &
                                 CLO(L,jlo,JNEQ)*CHN(L,JNEQ) + &
                                 HALF*CLO(L,jlo,JNEQ)*PI12LO(L,jlo,JNEQ) !non orthogonal term
                  else
                     C11R = ALO(L,jlo,JNEQ)*E(L,JNEQ) + &
                            HALF*(BLO(L,jlo,JNEQ) + &
                            CLO(L,jlo,JNEQ)*PI12LO(L,jlo,JNEQ)*(ELO(L,jlo,JNEQ) + &
                            E(L,JNEQ)))
                     C12R = BLO(L,jlo,JNEQ)*CHN(L,JNEQ)*E(L,JNEQ) +&
                            HALF*(ALO(L,jlo,JNEQ) + &
                            CLO(L,jlo,JNEQ)*PI12LO(L,jlo,JNEQ) + &
                            CLO(L,jlo,JNEQ)*PE12LO(L,jlo,JNEQ)*(ELO(L,jlo,JNEQ) +&
                            E(L,JNEQ)))
                  endif
                  C11I = ZERO
                  C12I = ZERO
                  do jlop=1,jlo
                     C3R(jlop)  = CLO(L,jlo,JNEQ)*pilolo(l,jlop,jlo,jneq) +&
                                  BLO(L,jlo,JNEQ)*PE12LO(L,jlop,JNEQ) + &
                                  ALO(L,jlo,JNEQ)*PI12LO(L,jlop,JNEQ)
                  enddo
                  C3I  = ZERO
                  do jlop=1,jlo
                     if (LOSECDER(L,jlo,JNEQ)) then
                        if (LOSECDER(L,jlop,JNEQ)) then
                           C13R(jlop) = E(L,JNEQ)*BLO(L,jlo,JNEQ)*PE12LO(L,jlop,JNEQ) + &
                                        BLO(L,jlo,JNEQ)*CHN(L,JNEQ) + &
                                        E(L,JNEQ)*CLO(L,jlo,JNEQ)*pilolo(l,jlop,jlo,jneq) + &
                                        TWO*CLO(L,jlo,JNEQ)*PE12LO(L,jlop,JNEQ) + &
                                        E(L,JNEQ)*ALO(L,jlo,JNEQ)*PI12LO(L,jlo,JNEQ) + & !non orthogonal term
                                        BLO(L,jlo,JNEQ)*PI12LO(L,jlo,JNEQ) !non orthogonal term
                        else
                           C13R(jlop) = HALF*(E(L,JNEQ)+ELO(L,jlop,JNEQ))*ALO(L,jlo,JNEQ)*PI12LO(L,jlop,JNEQ) + &
                                        HALF*(E(L,JNEQ)+ELO(L,jlop,JNEQ))*BLO(L,jlo,JNEQ)*PE12LO(L,jlop,JNEQ) + &
                                        HALF*(E(L,JNEQ)+ELO(L,jlop,JNEQ))*CLO(L,jlo,JNEQ)*pilolo(l,jlop,jlo,jneq) + &
                                        HALF*BLO(L,jlo,JNEQ)*PI12LO(L,jlop,JNEQ) + &
                                        CLO(L,jlo,JNEQ)*PE12LO(L,jlop,JNEQ)
                        endif
                     else
!                        if (LOSECDER(L,jlop,JNEQ)) then
!                           print *,"Something wrong in case.in1 file, second derivative LO must be last"
!                        endif
                        C13R(jlop) =   HALF*CLO(L,jlo,JNEQ)*(ELO(L,jlo,JNEQ)+ELO(L,jlop,JNEQ))*pilolo(l,jlop,jlo,jneq)+ &
                                       HALF*(BLO(L,jlo,JNEQ)*PI12LO(L,jlop,JNEQ) + &
                                            (ELO(L,jlop,JNEQ) + E(L,JNEQ))*(ALO(L,jlo,JNEQ)*PI12LO(L,jlop,JNEQ)+&
                                            BLO(L,jlo,JNEQ)*PE12LO(L,jlop,JNEQ)))
                     endif
                  enddo
                  C13I = ZERO

!
! Precalculate surface kinetic energy factor
!
                  akinloR=HALF*r(jneq)**2*(alo(l,jlo,jneq)*DU(L,JNEQ)+ &
                       blo(l,jlo,jneq)*DUP(L,JNEQ)+ &
                       CLO(l,jlo,jneq)*DPLO(l,jlo,jneq))

                  jeq0_loop: DO JEQO = 1, MULT(JNEQ)
                     mo1_loop:   DO MO1 = -L, L

                        NNLO = NNLO + 1

!    compute only local columns
!
                        IHELP=MOD(I_g,BLOCKSIZE)
                        IF (IHELP.EQ.0) IHELP=BLOCKSIZE

                        jeq_loop: DO JEQ = 1, MULT(JNEQ)
                           INDATM = INDATM + 1

!
!        calculate spherical-harmonics for Alm, Blm
!
                           VEC(1) = XK(NV+NNLO)
                           VEC(2) = YK(NV+NNLO)
                           VEC(3) = ZK(NV+NNLO)
                           CALL ROTATE(VEC,ROTIJ(1,1,INDATM),ROTV1)
                           CALL ROTATE(ROTV1,ROTLOC(1,1,JNEQ),ROTV2)
                           jnnlo=NV+NNLO
                           CALL YLM(ROTV2,L,YL(0,jnnlo,JEQ))

#ifdef Parallel
                           ipc=mod((i_g-1)/blocksize,npcol)
           local_column:   if (ipc.eq.mycolhs) then
#endif                           
!
!        determine phase factor
!
#if defined(_OPENMP) && !defined (Parallel)
! cycle is not allowed in the openmp simd loop
!$omp simd private(argx,argy,argz)
#endif
                              DO J_g = 1, NV + NNLO
#ifdef Parallel
                                    ipr=mod((j_g-1)/blocksize,nprow)
                                    if (ipr.ne.myrowhs) cycle
#endif                           

                                 ARGX = POSTwoPi(1,INDATM)*(KZZ(1,NV+NNLO)-KZZ(1,J_g))
                                 ARGY = POSTwoPi(2,INDATM)*(KZZ(2,NV+NNLO)-KZZ(2,J_g))
                                 ARGZ = POSTwoPi(3,INDATM)*(KZZ(3,NV+NNLO)-KZZ(3,J_g))
                                 ARGX=(ARGX+ARGY+ARGZ)
                                 PHSC(J_g) = DCMPLX(DCOS(ARGX),-DSIN(ARGX) )
                              enddo

                              mo2_loop: DO MO2 = -L, L

                                 L0M0 = L*(L+1) + MO2
!
!        CLYLN = (CIMAG**L)*DCONJG(YL(L0M0,NV+NNLO,JEQ))
!
                                 CLYLNR = CLR*DBLE(YL(L0M0,jNNLO,JEQ)) + &
                                      CLI*DIMAG(YL(L0M0,jNNLO,JEQ))
                                 CLYLNI = CLI*DBLE(YL(L0M0,jNNLO,JEQ)) - &
                                      CLR*DIMAG(YL(L0M0,jNNLO,JEQ))
!
!        update C1,C2,C3 for the Overlap-matrix 
!        C11,C12,C13 for the Hamilton-matrix 
!
                                 C1RUP = C1R*CLYLNR*PI4R2V
                                 C1IUP = C1R*CLYLNI*PI4R2V
                                 C2RUP = C2R*CLYLNR*PI4R2V
                                 C2IUP = C2R*CLYLNI*PI4R2V
                                 do jlop=1,jlo
                                    C3RUP(jlop) = C3R(jlop)*CLYLNR*PI4R2V
                                    C3IUP(jlop) = C3R(jlop)*CLYLNI*PI4R2V
                                 enddo
                                 C11RUP = C11R*CLYLNR*PI4R2V
                                 C11IUP = C11R*CLYLNI*PI4R2V
                                 C12RUP = C12R*CLYLNR*PI4R2V
                                 C12IUP = C12R*CLYLNI*PI4R2V
                                 do jlop=1,jlo
                                    C13RUP(jlop) = C13R(jlop)*CLYLNR*PI4R2V
                                    C13IUP(jlop) = C13R(jlop)*CLYLNI*PI4R2V
                                 enddo
                                 akinlorup = akinloR*CLYLNR*PI4R2V
                                 akinloiup = akinloR*CLYLNI*PI4R2V

                                 j=0
                                 DO J_g = 1, NV
#ifdef Parallel
                                    ipr=mod((j_g-1)/blocksize,nprow)
                                    if (ipr.ne.myrowhs) cycle
#endif
                                    j=j+1

!
!        CLYL = (CIMAG**L)*DCONJG(YL(L0M0,J,JEQ))
!
                                    CLYLR = CLR*DBLE(YL(L0M0,J_g,JEQ)) + &
                                         CLI*DIMAG(YL(L0M0,J_g,JEQ))
                                    CLYLI = CLI*DBLE(YL(L0M0,J_g,JEQ)) - &
                                         CLR*DIMAG(YL(L0M0,J_g,JEQ))
!
!        C6 ... Alm, C7 ... Blm
!
                                    C6R = CLYLR*AL_r(J,L)*PI4R2V
                                    C6I = CLYLI*AL_r(J,L)*PI4R2V
                                    C7R = CLYLR*BL_r(J,L)*PI4R2V
                                    C7I = CLYLI*BL_r(J,L)*PI4R2V
!
!        CAB = DCONJG(C6)*C1 + DCONJG(C7)*C2
!
                                    CABR =  C6R*C1RUP + C7R*C2RUP + &
                                         C6I*C1IUP + C7I*C2IUP
                                    CABI = -C6I*C1RUP - C7I*C2RUP + &
                                         C6R*C1IUP + C7R*C2IUP
!
!        update Overlap-matrix
!
                                    SPANEL(J,IHELP) = SPANEL(J,IHELP) + &
                                         (CABR*DBLE(PHSC(J_g)) - CABI*DIMAG(PHSC(J_g)))
!_COMPLEX                           SPANEL(J,IHELP) = SPANEL(J,IHELP) +  &
!_COMPLEX                                    CIMAG *  (CABR*DIMAG(PHSC(J_g)) + &
!_COMPLEX                                             CABI*DBLE(PHSC(J_g)))
!
!        CAB = DCONJG(C6)*C11 + DCONJG(C7)*C12
!
                                    CABR =  C6R*C11RUP + C7R*C12RUP + &
                                         C6I*C11IUP + C7I*C12IUP
                                    CABI = -C6I*C11RUP - C7I*C12RUP + &
                                         C6R*C11IUP + C7R*C12IUP
                                    akinbr = (c6r*u(l,jneq)+c7r*up(l,jneq)) * akinlorup &
                                         +(c6i*u(l,jneq)+c7i*up(l,jneq)) * akinloiup
                                    akinbi = (c6r*u(l,jneq)+c7r*up(l,jneq)) * akinloiup &
                                         -(c6i*u(l,jneq)+c7i*up(l,jneq)) * akinlorup
!
!        update Hamilton-matrix
!
                                    HPANEL(J,IHELP) = HPANEL(J,IHELP) + &
                                         (CABR*DBLE(PHSC(J_g)) - CABI*DIMAG(PHSC(J_g)))
!_COMPLEX                           HPANEL(J,IHELP) = HPANEL(J,IHELP) + &
!_COMPLEX                                       CIMAG *(CABR*DIMAG(PHSC(J_g)) + &
!_COMPLEX                                               CABI*DBLE(PHSC(J_g)))
                                    HPANEL(J,IHELP) = HPANEL(J,IHELP) + &
                                         (akinbr*DBLE(PHSC(J_g)) - &
                                         akinbi*DIMAG(PHSC(J_g)))
!_COMPLEX                           HPANEL(J,IHELP) = HPANEL(J,IHELP) +  &
!_COMPLEX                                     CIMAG *  (akinbr*DIMAG(PHSC(J_g)) + &
!_COMPLEX                                               akinbi*DBLE(PHSC(J_g)))
                                 enddo
                                 
                                 jlop_loop: do  jlop=1,jlo
                                       
!                                    write(0,'(7(a,i3))') 'jlop=',jlop,'  jlo=',jlo,'  mo1=',mo1,'  jeg=',jeqo,'  nnlo=',NNLO,'  L=',L,'  IHELP=',ihelp
                                    nv1=NV + NNLO - (MO1+L) - (2*L+1)*(JEQO-1) - (jlo-jlop)*(2*l+1)*mult(jneq)
                                    nv2=nv1+(2*l+1)*mult(jneq)-1
                                    nv2=min(nv2,I_g)

!                                    write(0,*) nv1,nv2,i_g

!                                    DO J_g = NV + NNLO - (MO1+L) - (2*L+1)*(JEQO-1) - (jlo-jlop)*(2*l+1)*mult(jneq),&
!                                             NV + NNLO - (jlo-jlop)*(MO1+L+1)-(jlo-jlop)*(JEQO-1)*(2*l+1)

                                    do j_g=nv1,nv2

#ifdef Parallel
                                       ipr=mod((j_g-1)/blocksize,nprow)
                                       if (ipr.ne.myrowhs) cycle
                                       j=global2local(j_g,nprow,blocksize)
#else
                                       j=j_g
#endif

!
!        CLYL = (CIMAG**L)*DCONJG(YL(L0M0,J,JEQ))
!
                                       CLYLR = CLR*DBLE(YL(L0M0,J_g,JEQ)) + &
                                            CLI*DIMAG(YL(L0M0,J_g,JEQ))
                                       CLYLI = CLI*DBLE(YL(L0M0,J_g,JEQ)) - &
                                            CLR*DIMAG(YL(L0M0,J_g,JEQ))
!
!        C6 ... CONJG(ALOlm), C7 ... CONJG(BLOlm), C8 ... CONJG(CLOlm)
!
                                       C6R = CLYLR*ALO(L,jlop,JNEQ)*PI4R2V
                                       C6I = CLYLI*ALO(L,jlop,JNEQ)*PI4R2V
                                       C7R = CLYLR*BLO(L,jlop,JNEQ)*PI4R2V
                                       C7I = CLYLI*BLO(L,jlop,JNEQ)*PI4R2V
                                       C8R = CLYLR*CLO(L,jlop,JNEQ)*PI4R2V
                                       C8I = CLYLI*CLO(L,jlop,JNEQ)*PI4R2V
!                                                         
!        CABC = DCONJG(C6)*C1 + DCONJG(C7)*C2 + DCONJG(C8)*C3
!
                                       CABCR =  C6R*C1RUP + C6I*C1IUP + &
                                            C7R*C2RUP + C7I*C2IUP +    &
                                            C8R*C3RUP(jlop) + C8I*C3IUP(jlop)    
                                       CABCI = -C6I*C1RUP + C6R*C1IUP - &
                                            C7I*C2RUP + C7R*C2IUP - &
                                            C8I*C3RUP(jlop) + C8R*C3IUP(jlop) 
!
!        update Overlap-matrix
!
                                       SPANEL(J,IHELP) = SPANEL(J,IHELP) + &
                                            (CABCR*DBLE(PHSC(J_g)) - &
                                            CABCI*DIMAG(PHSC(J_g)))

!_COMPLEX                              SPANEL(J,IHELP) = SPANEL(J,IHELP) +  &
!_COMPLEX                                     CIMAG *  (CABCR*DIMAG(PHSC(J_g)) + &
!_COMPLEX                                               CABCI*DBLE(PHSC(J_g)))
!
!        CABC = DCONJG(C6)*C11 + DCONJG(C7)*C12 + DCONJG(C8)*C13
!
                                       CABCR =  C6R*C11RUP + C6I*C11IUP + &
                                            C7R*C12RUP + C7I*C12IUP +    &
                                            C8R*C13RUP(jlop) + C8I*C13IUP(jlop)    
                                       CABCI = -C6I*C11RUP + C6R*C11IUP - &
                                            C7I*C12RUP + C7R*C12IUP - &
                                            C8I*C13RUP(jlop) + C8R*C13IUP(jlop) 
!
!        update Hamilton-matrix
!
                                       HPANEL(J,IHELP) = HPANEL(J,IHELP) + &
                                            (CABCR*DBLE(PHSC(J_g)) - &
                                            CABCI*DIMAG(PHSC(J_g)))
!_COMPLEX                                 HPANEL(J,IHELP) = HPANEL(J,IHELP) +  &
!_COMPLEX                                      CIMAG * (CABCR*DIMAG(PHSC(J_g)) + &
!_COMPLEX                                               CABCI*DBLE(PHSC(J_g)))
                                    enddo

                                 enddo jlop_loop
                              enddo mo2_loop

#ifdef Parallel
                     endif local_column
#endif
                              
                        enddo jeq_loop

                        INDATM = INDATM - MULT(JNEQ)

                        IF(I_g.EQ.IMAX) THEN

                           ii=0
                           DO II_g = IMIN, IMAX
#ifdef Parallel
                              ipc=mod((ii_g-1)/blocksize,npcol)
                              if (ipc.ne.mycolhs) cycle
                              ii=global2local(ii_g,npcol,blocksize)
                              IHELP = MOD(II_g,BLOCKSIZE)
                              IF(IHELP.EQ.0) IHELP=BLOCKSIZE
!_REAL                              CALL DCOPY(HSdim_r, HPANEL(1,IHELP), 1, H(1,II), 1 )
!_REAL                              CALL DCOPY(HSdim_r, SPANEL(1,IHELP), 1, S(1,II), 1 )
!_COMPLEX                           CALL ZCOPY(HSdim_r, HPANEL(1,IHELP), 1, H(1,II), 1 )
!_COMPLEX                           CALL ZCOPY(HSdim_r, SPANEL(1,IHELP), 1, S(1,II), 1 )
#else
                              ii=ii+1
                              ii=ii_g
                              IHELP = MOD(II,BLOCKSIZE)
                              IF(IHELP.EQ.0) IHELP=BLOCKSIZE
!_REAL                              CALL DCOPY(II-1, HPANEL(1,IHELP), 1, HS(II, 1), HSROWS )
!_REAL                              CALL DCOPY(II, SPANEL(1,IHELP), 1, HS(1,II), 1 )
!_COMPLEX                           CALL ZCOPY(II-1, HPANEL(1,IHELP), 1, HS(II, 1), HSROWS )
!_COMPLEX                           CALL ZCOPY(II, SPANEL(1,IHELP), 1, HS(1,II), 1 )
                              HSDIAG(II)=HPANEL(II,IHELP)
#endif
                           END DO

!_REAL                           HPANEL(:,:)=ZERO
!_REAL                           SPANEL(:,:)=ZERO
!_COMPLEX                        HPANEL(:,:)=CZERO
!_COMPLEX                        SPANEL(:,:)=CZERO
                           IMIN = IMAX+1
                           IMAX = MIN(NV + NLO,((IMIN-1)/(BLOCKSIZE*NPCOL)+1)*BLOCKSIZE*NPCOL)
                              
                        END IF

                        I_g = I_g + 1
                        
                     enddo mo1_loop
                  enddo jeq0_loop
                  
               enddo jlo_loop
            enddo l_loop
            TMP = CLR
            CLR = -CLI
            CLI = TMP
            INDATM = INDATM + MULT(JNEQ)

         enddo jneq_lo_loop

         nevertrue: if(imax.gt.imin) THEN
            ii=0
            DO II_g = IMIN, IMAX
#ifdef Parallel
               ipc=mod((ii_g-1)/blocksize,npcol)
               if (ipc.ne.mycolhs) cycle
               ii=global2local(ii_g,npcol,blocksize)
               IHELP = MOD(II_g,BLOCKSIZE)
               IF(IHELP.EQ.0) IHELP=BLOCKSIZE
!_REAL         CALL DCOPY(HSdim_r, HPANEL(1,IHELP), 1, H(1,II), 1 )
!_REAL         CALL DCOPY(HSdim_r, SPANEL(1,IHELP), 1, S(1,II), 1 )
!_COMPLEX      CALL ZCOPY(HSdim_r, HPANEL(1,IHELP), 1, H(1,II), 1 )
!_COMPLEX      CALL ZCOPY(HSdim_r, SPANEL(1,IHELP), 1, S(1,II), 1 )
#else
               ii=ii+1
               ii=ii_g
               IHELP = MOD(II,BLOCKSIZE)
               IF(IHELP.EQ.0) IHELP=BLOCKSIZE
!_REAL         CALL DCOPY(II-1, HPANEL(1,IHELP), 1, HS(II, 1), HSROWS )
!_REAL         CALL DCOPY(II, SPANEL(1,IHELP), 1, HS(1,II), 1 )
!_COMPLEX      CALL ZCOPY(II-1, HPANEL(1,IHELP), 1, HS(II, 1), HSROWS )
!_COMPLEX      CALL ZCOPY(II, SPANEL(1,IHELP), 1, HS(1,II), 1 )
               HSDIAG(II)=HPANEL(II,IHELP)
#endif
            END DO
         END IF nevertrue
	 CALL STOP_TIMER(time_lo)
!         WRITE (6,'(a,2f12.1)') 'Time for los      (hamilt, cpu/wall) : ', &
!                             READ_CPU_TIME(time_lo),READ_wall_TIME(time_lo)

         deallocate( PHSC )
         deallocate( YL )
         call end_albl('HAM')

      ENDIF
!print*,hsdiag
      deallocate( hpanel,spanel)
!$omp end single
!$omp end parallel
      WRITE (6,'(a,2f12.2)') 'Time for al,bl    (hamilt, cpu/wall) : ', time_albl_tot, time_albl_totw
      WRITE (6,'(a,2f12.2)') 'Time for legendre (hamilt, cpu/wall) : ', DTIMLG, DTIMLGw
      WRITE (6,'(a,2f12.2)') 'Time for phase    (hamilt, cpu/wall) : ', DTIMPH, DTIMPHw
      WRITE (6,'(a,2f12.2)') 'Time for us       (hamilt, cpu/wall) : ', DTIMUS, DTIMUSw
      WRITE (6,'(a,2f12.2)') 'Time for overlaps (hamilt, cpu/wall) : ', DTIMMA, DTIMMAw
      WRITE (6,'(a,2f12.2)') 'Time for distrib  (hamilt, cpu/wall) : ', DTIMDIST, DTIMDISTw
!      WRITE (6,'(a,2f12.2)') 'Time sum iouter   (hamilt, cpu/wall) : ', &
!                     READ_CPU_TIME(time_loop260), READ_wall_TIME(time_loop260)
         WRITE (6,'(a,2f12.2)') 'Time sum iouter   (hamilt, cpu/wall) : ', &
                     READ_CPU_TIME(time_loop260), READ_wall_TIME(time_loop260)
         WRITE (6,'(a,2f12.2)') 'Time for los      (hamilt, cpu/wall) : ', &
                             READ_CPU_TIME(time_lo),READ_wall_TIME(time_lo)


      deallocate( rk_tmp )
      deallocate( XL )
      deallocate( kzz_TMP )
      deallocate( xk_tmp_r,yk_tmp_r,zk_tmp_r)

#ifdef Barrier
      CALL BARRIER
#endif

      RETURN
!
!         End of 'HAMILT'
!
      END
