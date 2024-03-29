! @Copyright 2007 Kristjan Haule
! 

SUBROUTINE cmp_overlap(projector,Olapm, SOlapm, Qcomplex, nsymop, csort, iorbital, cix_orb, nindo, cixdim, iSx, noccur, cfX, crotloc_x_rotij, maxucase, maxdim2, norbitals, pr_proc, SIMPLE, Qsymmetrize)
  USE com_mpi,  ONLY: myrank, master, vectors, nvector, vector_para, fvectors, AllReduce_MPI, Qprint
  USE com,      ONLY: MINWAV, MAXWAV, iso, emin, emax
  USE abc,      ONLY: KX, KY, KZ, BK3, E, A, ALM, ALML, aK
  USE param,    ONLY: nmat , LMAX2, iblock, nrf, NRAD, LOMAX, nloat, nemin0, nemax0, max_nl
  USE w_atpar,  ONLY: w_FJ, w_DFJ, w_jatom, w_alo, w_nlo, w_nlov, w_nlon, w_ilo, w_lapw, w_ri_mat, w_P, w_DP, nnlo
  USE structure,ONLY: RMT, VOL, pos, tau, nat, rotij, tauij, iord, iz, BR1, rot_spin_quantization
  USE case,     ONLY: iatom, isort, nl, ll, crotloc, shft, maxdim, ncix, natom, Sigind, csize, maxsize, cix, ifirst
  USE sym2,     ONLY: tmat, idet, iz_cartesian
  USE kpts,     ONLY: mweight, tweight, numkpt
  USE p_project,ONLY: P_rfi, rix_mat, n_al_ucase, max_lcase, al_ucase, l_al_ucase, j_al_ucase, kNri, rix, dri, al_interstitial
  USE matpar,   only: atpar, alo, nlo, nlov, nlon, ilo, lapw, RI_MAT, P, DP
  IMPLICIT NONE
  INTEGER, intent(in)     :: projector
  COMPLEX*16, intent(out) :: Olapm(maxdim,maxdim,ncix), SOlapm(maxdim,maxdim,ncix)
  LOGICAL, intent(in)     :: Qcomplex
  INTEGER, intent(in)     :: nsymop, csort(nat), pr_proc
  INTEGER, intent(in)     :: iorbital(natom,lmax2+1), cix_orb(norbitals), nindo(norbitals)
  INTEGER, intent(in)     :: cixdim(ncix), iSx(maxdim2, norbitals), noccur(maxsize,ncix)
  COMPLEX*16, intent(in)  :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  INTEGER, intent(in)     :: maxucase, maxdim2, norbitals
  REAL*8,  intent(in)     :: crotloc_x_rotij(3,3,max_nl,natom)
  LOGICAL, intent(in)     :: SIMPLE, Qsymmetrize
  ! interfaces
  interface
     REAL*8 Function Romb(y,N,dx)
       IMPLICIT NONE
       REAL*8, intent(in) :: y(N)
       REAL*8, intent(in) :: dx
       INTEGER, intent(in):: N
     end Function Romb
     REAL*8 FUNCTION detx(a)
       IMPLICIT NONE
       REAL*8, intent(in) :: a(3,3)
     end FUNCTION detx
  end interface
  ! locals
  real*8,     allocatable :: a_real(:)
  complex*16, allocatable :: DMFTU(:,:,:), olp(:,:)
  complex*16, allocatable :: Olapmk(:,:,:)
  COMPLEX*16, allocatable :: URx(:,:), tmp(:,:)
  !COMPLEX*16, allocatable :: Uu(:,:), Ud(:,:)
  COMPLEX*16, allocatable :: a_interstitial(:,:,:), h_interstitial(:,:,:)
  REAL*8,     allocatable  :: phi_jl(:,:)
  REAL*8,     allocatable :: aKR(:), jlr(:), jlp(:)
  COMPLEX*16, allocatable :: olocf(:,:), work(:),tmp1(:,:), solap(:,:)
  REAL*8,     allocatable :: ws(:), rwork(:)
  INTEGER,    allocatable :: iwork(:)
  !
  INTEGER, allocatable :: cind(:), cini(:)
  INTEGER :: cixdms
  !
  INTEGER      :: ikp, iikp, iks, nkp, ivector, N, NE, NEMIN, NEMAX, nbands, NUM, isym, icase, lcase, lfirst
  INTEGER      :: latom, jatom, iucase, i, IND_YL, ibb, lda, ldb, ldc, l1
  INTEGER      :: nind, idt, icix, iorb, m1, lms1, ind1, ind2, is1, iorb1, iorb2, nind1, nind2
  INTEGER      :: cixdm, ip, iq, it, is, itape, L, M, ii, i3, N1, num1, irf1, irf2
  logical      :: Tcompute
  character*10 :: KNAME
  real*8       :: S, T, Z, exxx, FAC, PI, TWOPI, ARG123, ARG2, ARGT, EMIST, rx
  complex*16   :: h_yl(2*LMAX2+1,iblock)
  complex*16   :: h_alyl(2*LMAX2+1,iblock,2), h_blyl(2*LMAX2+1,iblock,2)
  complex*16   :: YL((LMAX2+1)*(LMAX2+1))
  complex*16   :: PHSHEL, CFAC, IMAG, cc, corc, sm1, sm2, ff
  real*8       :: Olapc(maxsize,ncix), ARGT2, tmp3(3,3)
!  COMPLEX*16   :: tuu, tud, tdu, tdd
  CHARACTER*200:: FNAME
  !
  real*8       :: BK(3), BKROT(3), BKROT2(3), BKRLOC(3)
  REAL*8       :: Kn(3), crotloc_x_BR1(3,3)
  INTEGER      :: lwork, lrwork, liwork, info, iind, ir, Nri, isize, isize_
  LOGICAL      :: nonzero_interst
  REAL*8       :: Det, phi1, the1, psi1
  COMPLEX*16   :: Rispin(2,2)
  INTEGER      :: Nstar
  REAL*8       :: k_star(3,iord), Trans3(3,3)
  INTEGER      :: gind(iord), nsymop2, igi, j
  !------------------------------------------------------------------

  DATA IMAG/(0.0D0,1.0D0)/
  
  PI=ACOS(-1.0D0)
  TWOPI=2.D0*PI

  if (Qprint) WRITE(6,'(A,I3,2x,A,I3)') 'pr_proc=', pr_proc, 'tot-k=', numkpt
  
  ALLOCATE( a_real(nmat) ) !---  for eigenvectors -------------!
  allocate( olp(maxdim2,maxdim2) )
  allocate( Olapmk(maxdim,maxdim,ncix) )
  Olapm=0

  iikp=0
  DO ivector=1,nvector
     
     DO is=1,iso    !------ over up/dn ---------------------!
        itape=8+is
        if (vector_para) then
           FNAME = fvectors(ivector,is)
           open(itape,FILE=FNAME,STATUS='old',FORM='unformatted')
        else
           !------------ vector files need to be read from the beginning -------------------!           
           rewind(itape) !--- both vector files: 9 and 10 rewind -------------------------!
        endif
        DO I=1,NAT
           READ(itape) EMIST !---- At the beginninge we have linearization energies --!
           READ(itape) EMIST !---- We just skip them ---------------------------------!
        ENDDO
     END DO
     
     if (vector_para) then
        nkp = vectors(ivector,2)
     else
        nkp= numkpt
     endif
     
     DO iks=1,nkp   !------ Over all irreducible k-points ------!
        if (vector_para) then
           ikp = vectors(ivector,3)+iks  ! successive index in k-point table from case.klist
           iikp = iikp+1                 ! successive index in k-point table on this processor
        else
           ikp = iks                     ! successive index in k-point table from case.klist
           !--- We need to go over all k-points even though we will compute only some of them on this processor.
           !--- This is because we need to read vector file sequentially.
           iikp = ikp-myrank*pr_proc            ! The index of the point to compute. If negative, do not compute! 
        endif
        
        Tcompute=.FALSE.
        if (iikp.gt.0) Tcompute=.TRUE.       ! If Tcompute is true, the point needs to be computed.
        if (iikp.GT.pr_proc) EXIT            ! Processor finished. The rest of the points will be taken care of by the other processors.
        
        if (Tcompute .and. abs(projector).eq.5) allocate( phi_jl(nmat, n_al_ucase) )
              
        !------- reading from vector for both spins -----------!
        DO is=1,iso    !------ over up/dn ---------------------!
           itape=8+is
           READ(itape,END=998) S,T,Z,KNAME,N,NE  !--- Here we can jupm out of loop 4 and stop at 998 -----------------------------------!
           IF(N.GT.MAXWAV) MAXWAV=N                                         
           IF(N.LT.MINWAV) MINWAV=N
           READ(itape) (KX(I),KY(I),KZ(I),I=1,N) !--- Reads all reciprocal vectors -----------------------------------------------------!
           NEMIN=1
           NEMAX=0
           NUM=0
           DO WHILE (NUM.NE.NE)
              READ(itape) NUM,exxx                !----- eigenvalues read -----!
              E(NUM)=exxx
              if (.not.Qcomplex) then
                 READ(itape) (A_real(I),I=1,N)    !----- eigenvector read -----!
                 a(:,num,is) = a_real(:)
              else
                 READ(itape) (A(I,NUM,is),I=1,N)  !--- eigenvector complex of size n (where n is number of reciprocal vectors) --!
              end if
              if (abs(projector).LT.4) then
                 IF(E(NUM).LT.EMIN) NEMIN=NEMIN+1
                 IF(E(NUM).LT.EMAX) NEMAX=NEMAX+1
              !else
              !   nemin=nemin0
              !   nemax=nemax0
              endif
           ENDDO
           
           if (abs(projector).GE.4) then
              nemin=nemin0
              nemax=min(nemax0,NE)
           endif
           
           IF (Tcompute .and. is.eq.1) THEN
              isize=N-nnlo
              DO I=1,N                                !----- Over reciprocal vectors --------------------------------------------------------!
                 BK3(1,I)=(S+KX(I))                   !----  Creates K+k where K is reciprocal vector and k is irreducible k-point ----------!
                 BK3(2,I)=(T+KY(I))
                 BK3(3,I)=(Z+KZ(I))
                 Kn(:) = matmul(BR1,BK3(:,I))
                 aK(I) = sqrt(Kn(1)**2+Kn(2)**2+Kn(3)**2)
              ENDDO
              DO iucase=1,maxucase
                 jatom = w_jatom(iucase)
                 !--- Computes Bessel functions for all needed l's and atoms. --------------------------------!
                 !---  For optimization purposes, done only for irreducible k-points -------------------------!
                 !---  output: FJ, DFJ / Spherical Bessel functions and derivative, uses common block BR1/general ---!
                 CALL HARMON(isize,aK(:isize),LMAX2,w_FJ(:,:isize,iucase),w_DFJ(:,:isize,iucase),RMT(jatom))
              ENDDO
              ! New
              if (abs(projector).EQ.5) then
                 Nri=2**kNri+1    ! Number of radial points in the interstitials
                 phi_jl(:,:)=0
                 allocate( aKR(Nri), jlr(Nri), jlp(Nri) )
                 DO iind=1,n_al_ucase
                    if (abs(dri(iind)).gt.1e-10) then
                       l     = l_al_ucase(iind)
                       jatom = j_al_ucase(iind)
                       DO i=1,isize         ! over all reciprocal vectors K
                          rx = RMT(jatom)
                          do ir=1,Nri
                             aKR(ir)=rx*aK(i)    ! |k+K|*r
                             rix(ir)=rx          !  r
                             rx = rx + dri(iind)
                          enddo
                          CALL sphbes2(l,Nri,aKR,jlr)  ! spherical bessel : j_l(|k+K|*r)
                          jlr(:) = jlr(:)*rix(:)       ! r*j_l . We need to do that, because the integral is Int[ (phi(r)/r)*j_l(r)*r^2]=Int[phi(r)*j_l(r)*r]
                          jlp(:) = jlr(:)*P_rfi(:,iind)
                          phi_jl(i,iind) = romb(jlp, Nri, dri(iind))  ! Integration over r on rix mesh
                       ENDDO
                    endif
                 ENDDO
                 deallocate( aKR, jlr, jlp )
              endif
              ! New
           ENDIF
        ENDDO  !
        if (.not.Tcompute) CYCLE  ! This k-points was read, but will not be computed on this processor
        !---  Finished reading eigenvectors and eigenvalues for both spins -----!
        !---   E(iband) contains eigenvalues   ---------------------------------!
        !---   A(:,iband,is) contains eigenvectors -----------------------------!
        nbands = nemax-nemin+1
        allocate( DMFTU(nbands,maxdim2,norbitals) )
        allocate( URx(nbands,maxdim2), tmp(nbands,maxdim2) )
        
        if (abs(projector).eq.5) then
           allocate( h_interstitial(2*LMAX2+1,iblock,2) )
           allocate( a_interstitial(2*LMAX2+1,nbands,iso), al_interstitial(2*LMAX2+1,nbands,iso,max_lcase) )
        endif

        if (.False.) then
           ! Finds all star members and corresponding group operations
           ! This is switched off, as there are some cases of selection of G vectors, which are not completely symmetric
           ! Taking only the star members can then lead to less symmetry as expected.
           BK(1)=S
           BK(2)=T
           BK(3)=Z
           call STERN(BK, Nstar, k_star, gind, iz, tau, iord, Qcomplex)
           nsymop2 = min(nsymop,Nstar)
        else
           nsymop2=nsymop
           do igi=1,iord
              gind(igi)=igi
           enddo
        endif
        !print *, 'tweight=', tweight, 'nsymop=', nsymop2, 'mweight=', mweight(ikp)

        !----------- sum over all k-points in the star of the irreducible k-point ---------------------------------!
        DO igi=1,nsymop2  !-- Over all symmetry operarations -> all k-points in the star of irreducible k point --!
           isym = gind(igi)
           DMFTU=0
           DO icase=1,natom  !--------------- over all atoms requested in the input ------------------------!
              latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
              jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )
              lfirst = ifirst(latom)
              iucase = csort(jatom)  ! The renumbert sorts, such that the required atoms from the input give continuous index
              !----  Setting all values of common-blocks which depend on atom and were computed by atpar ---!
              alo(:,:,:,:) = w_alo(:,:,:,:,iucase)
              nlo = w_nlo(iucase)                 
              nlov = w_nlov(iucase)               
              nlon = w_nlon(iucase)               
              ilo(:) = w_ilo(:,iucase)            
              lapw(:) = w_lapw(:,iucase)          
              ri_mat(:,:,:,:) = w_ri_mat(:,:,:,:,iucase)
              P(:,:,:) = w_p(:,:,:,iucase)        
              DP(:,:,:) = w_dp(:,:,:,iucase)      
              if ((nlo+nlon+nlov).NE.nnlo) then
                 WRITE(6,*) 'ERROR: nlo+nlon+nlov should be equal to nnlo but is not', nlo+nlon+nlov, nnlo
                 STOP
              endif
              isize = N-(nlo+nlon+nlov)
              
              if (abs(projector).eq.5) al_interstitial(:,:,:,:)=0.d0

              FAC=4.0D0*PI*RMT(jatom)**2/SQRT(VOL)
              do lcase=1,nl(icase) !----------- loop over L(jatom) requested in the ionput ---------------!
                 l=ll(icase,lcase) !------ current L --!

                 if (iso.eq.2) then
                    !!  local_axis_defined_by_locrot  <- local_axis_of_equivalent_atom <- group_operation_symmetry <- from_spin_quantization_to_global_cartesian
                    !!* Trans3 = crotloc(:,:,icase) * rotij_cartesian * iz_cartesian(:,:,isym) * rot_spin_quantization
                    tmp3 = matmul(iz_cartesian(:,:,isym), rot_spin_quantization)
                    Trans3 = matmul(crotloc_x_rotij(:,:,lcase,icase),tmp3)
                    Det = detx(Trans3)
                    Trans3 = transpose(Trans3*Det)
                    CALL Angles_zxz(phi1,the1,psi1, Trans3 )
                    CALL Spin_Rotation(Rispin,phi1,the1,psi1)
                 endif
                 crotloc_x_BR1(:,:) = matmul( crotloc(:,:,lcase,icase),BR1 )
                 
                 ALM = 0.0         !------  ALM(m,band,nrf,is) will hold product of eigenvectors and a/b expansion coefficients --!
                 nonzero_interst=.false.
                 if (abs(projector).eq.5) then
                    a_interstitial=0
                    iind=al_ucase(icase,lcase)
                    if (abs(dri(iind)).gt.1e-10) nonzero_interst=.true.
                 endif
                 !--------- blocks are for efficiency. Matrix is multiplied in block form. This must be important in the past, while modern BLAS should do that better. I think it is obsolete.
                 DO ii=1,isize,iblock !------ iblock is 128 for 32-bit system -------!
                    !-------- nlo-number of local orbitals -----!
                    i3=0
                    do i=ii,min(ii+iblock-1,isize)  ! 121
                       !---------  rotates ylm(k+K) to ylm(k'+K) where k' is in star of irreducible k. ------------!
                       i3=i3+1
                       ! BKROT = R_a.(k+K) transforms to the reducible k-point
                       BKROT = matmul(TMAT(:,:,isym), BK3(:,I))
                       ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom 
                       BKROT2 = matmul(rotij(:,:,latom), BKROT)
                       !---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
                       ! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
                       !BKROT3 = matmul(BR1,BKROT2)
                       !!---- BKRLOC = crotloc.R_n.R_a.(k+K),  rotates according to the user specified local coordinate system.
                       !BKRLOC = matmul(crotloc(:,:,icase), BKROT3)
                       BKRLOC = matmul(crotloc_x_BR1, BKROT2)
                       !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
                       CALL YLM (BKRLOC,LMAX2,YL)  ! 
                       ! ARG123 = (R_g.(k+K)) *  R(iatom) * 2pi
                       ARG123 = dot_product( BKROT2, POS(:,lfirst) )*TWOPI
                       ! ARGT = (k+K)*tau(isym) * 2pi
                       ARGT = dot_product(BK3(:3,i), TAU(:3,isym))*TWOPI
                       ! ARGT2 = (R_a.(k+K)).tau_n * 2pi
                       ARGT2 = dot_product(BKROT, tauij(:3,latom))*TWOPI
                       ! ARG2 = (R_a.(k+K)) *  shft * 2pi
                       ARG2  = dot_product(BKROT, shft(:3,latom) )*TWOPI
                       ! PHSEHL = e^{I*2pi*( (R_a.(k+K))*tau_n + (K+k)*tau(isym) + (R_n.R_a.(k+K)*R(first)))}
                       PHSHEL=EXP(IMAG*(ARG123+ARG2+ARGT+ARGT2))
                       DO  M=1,2*L+1
                          IND_YL=M+L*L
                          h_yl(M,i3)=conjg(yl(ind_yl))*phshel !----- h_yl is rotated yl when k is rotated to k' -----!
                       END DO
                    enddo
                 
                    DO is=1,iso  !--- over both spins
                       i3=0
                       do i=ii,min(ii+iblock-1,isize)
                          i3=i3+1
                          DO M=1,2*L+1
                             if (lapw(l)) then
                                h_ALYL(m,i3,is)=(w_DFJ(L,I,iucase)*P(l,is,2)-w_FJ(L,I,iucase)*DP(l,is,2))* h_yl(M,i3)!*ph_spin(is) ! derivatives of bessel functions and spheric harmonics
                                h_BLYL(m,i3,is)=(w_FJ(L,I,iucase)*DP(l,is,1)-w_DFJ(L,I,iucase)*P(l,is,1))* h_yl(M,i3)!*ph_spin(is) 
                             else
                                h_ALYL(m,i3,is)=w_FJ(L,I,iucase)/P(l,is,1)/RMT(jatom)**2*h_yl(M,i3)!*ph_spin(is)
                                h_BLYL(m,i3,is)=(0.d0,0.d0)
                             end if
                          END DO
                       enddo
                       !!!  New
                       if (abs(projector).eq.5 .and. nonzero_interst) then
                          do i=ii,min(ii+iblock-1,isize)
                             DO M=1,2*L+1
                                !    h_interstitial = <P_phi|j_l> Y_L^*(R(k+K))*exp(i*(k+K)*r_latom)*exp(i*spin_phase)/Rmt^2
                                h_interstitial(M,i-ii+1,is) = phi_jl(i,iind)*h_yl(M,i-ii+1)/RMT(jatom)**2!*ph_spin(is)
                             ENDDO
                          enddo
                       endif
                       !!!  New
                       ibb=min(iblock,isize-ii+1)
                       ldb=nmat                    
                       lda=2*LMAX2+1
                       ldc=lda
                       isize_=2*l+1
                       !---- h_alyl(2*lmax+1,iblock,is)  contains rotated Apw's, such that chi(r) = (Apw*u(r) + Bpw*udot(r))*Ylm(r)
                       !---- h_blyl(2*lmax+1,iblock,is)  contains rotated Bpw's
                       !---- A(:N,iband,is)              contains eigenvectors, also named C(k+G,inad,is)
                       !---- alm[2*l+1,iband][irf=1,is] += sum_{iK\in block} h_alyl[2*l+1,iK][is]*A[iK,iband][is]
                       !---- alm[2*l+1,iband][irf=2,is] += sum_{iK\in block} h_blyl[2*l+1,iK][is]*A[iK,iband][is]
                       !---- 
                       !---- The results is:
                       !---- alm[lm,iband][1,is] = sum_G Apw(lm,is,K+G) * C(k+G,iband,is)
                       !---- alm[lm,iband][2,is] = sum_G Bpw(lm,is,K+G) * C(k+G,iband,is)
                       !---- Where C(k+G,iband,is) are eigenvectors, and Apw and Bpw are expansion coefficients defined in Shick et.al., PRB 60, 10763 (1999).
                       call zgemm('N','N',isize_,nemax-nemin+1,ibb,(1.d0,0.d0), h_alyl(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), alm(1,nemin,1,is),ldc)
                       call zgemm('N','N',isize_,nemax-nemin+1,ibb,(1.d0,0.d0), h_blyl(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), alm(1,nemin,2,is),ldc)
                       if (abs(projector).eq.5 .and. nonzero_interst) then !! The new
                          call zgemm('N','N',isize_,nbands,ibb,(1.d0,0.d0),h_interstitial(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), a_interstitial(1,1,is),ldc)
                       endif !! The new
                    ENDDO !----------- over both spins ---------!
                 ENDDO    !----------- over iblock -------------!

                 !-------------- Adds localized orbitals to alm. -------------------!
                 if (nlo.ne.0) then
                    call lomain(nemin,nemax,lfirst,latom,n,jatom,isym,L,iso,crotloc_x_BR1)
                 end if

                 !do m=1,2*l+1
                 !   do num=nemin,nemax
                 !      do irf=1,nrf
                 !         WRITE(6,'(3I4,2x,4F10.6)') m,num,irf, alm(m,num,irf,1), alm(m,num,irf,2)
                 !      enddo
                 !   enddo
                 !enddo
                 
                 !--------- FACCFAC = (i)^l*4.*PI*RMT(JATOM)**2/SQRT(VOL) --------!
                 CFAC = FAC*IMAG**L      !------  (i)^l*fac -----!
                 
                 if (iso.eq.2) then
                    ALML(l,:(2*L+1),nemin:nemax,:nrf,1) = ALM(:(2*L+1),nemin:nemax,:nrf,1)*(CFAC*Rispin(1,1)) + ALM(:(2*L+1),nemin:nemax,:nrf,2)*(CFAC*Rispin(2,1))
                    ALML(l,:(2*L+1),nemin:nemax,:nrf,2) = ALM(:(2*L+1),nemin:nemax,:nrf,1)*(CFAC*Rispin(1,2)) + ALM(:(2*L+1),nemin:nemax,:nrf,2)*(CFAC*Rispin(2,2))
                    !ALML(l,:(2*L+1),nemin:nemax,:nrf,:iso) = ALM(:(2*L+1),nemin:nemax,:nrf,:iso)*CFAC
                    if (abs(projector).eq.5 .and. nonzero_interst) then
                       al_interstitial(:(2*L+1),:nbands,1,lcase) = a_interstitial(:(2*L+1),:nbands,1)*(CFAC*Rispin(1,1)) + a_interstitial(:(2*L+1),:nbands,2)*(CFAC*Rispin(2,1))
                       al_interstitial(:(2*L+1),:nbands,2,lcase) = a_interstitial(:(2*L+1),:nbands,1)*(CFAC*Rispin(1,2)) + a_interstitial(:(2*L+1),:nbands,2)*(CFAC*Rispin(2,2))
                       !al_interstitial(:(2*L+1),:nbands,:iso,lcase) = a_interstitial(:(2*L+1),:nbands,:iso)*CFAC
                    endif
                 else
                    ALML(l,:(2*L+1),nemin:nemax,:nrf,1) = ALM(:(2*L+1),nemin:nemax,:nrf,1)*CFAC   
                    if (abs(projector).eq.5 .and. nonzero_interst) al_interstitial(:(2*L+1),:nbands,1,lcase) = a_interstitial(:(2*L+1),:nbands,1)*CFAC 
                 endif
              enddo  !--------- end of atom L loop (alm)  ----------------!
              !----------------  ALML(l,m,iband,ifr,ispin) contains (A*C,B*C) for all l,m,iband,ifr=(1,2),ispin --------!
              !---------- Computing the DMFT transformation --------------!
              do lcase=1,nl(icase)
                 
                 l1=ll(icase,lcase)
                 nind=(2*l1+1)*iso
                 idt = idet(isym)
                 icix = cix(icase,lcase)
                 iorb = iorbital(icase,lcase)
                 
                 if (icix.gt.0) then
                    URx=0
                    N1=2*l1+1

                    if (abs(projector).eq.1) then
                       ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U
                       !$OMP PARALLEL DO SHARED(URx) PRIVATE(num1,is1,m1,lms1,ind1,cc,irf1) SCHEDULE(STATIC)
                       do num1=1,nbands      !nemin,nemax     ! over bands
                          do is1=1,iso            ! over spin-1
                             do m1=-l1,l1           ! over m-1
                                lms1=l1+1+m1
                                ind1=l1+1+idt*m1+N1*(is1-1)
                                cc=0
                                do irf1=1,nrf
                                   cc = cc + alml(l1,lms1,num1+nemin-1,irf1,is1)*ri_mat(irf1,1,l1,is1)*(idt)**m1
                                enddo
                                if (idt.gt.0) then
                                   URx(num1,ind1) = cc
                                else
                                   URx(num1,ind1) = conjg(cc)
                                endif
                             enddo
                          enddo
                       enddo
                       !$OMP END PARALLEL DO
                    elseif (abs(projector).eq.2 .or. abs(projector).eq.4) then
                       ! L -> l1
                       ! nrfmax -> nrf
                       ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U
                       !$OMP PARALLEL DO SHARED(URx) PRIVATE(num1,is1,m1,lms1,ind1,cc,sm1,sm2,irf1,irf2,ff) SCHEDULE(STATIC)
                       do num1=1,nbands      !nemin,nemax     ! over bands
                          do is1=1,iso            ! over spin-1
                             do m1=-l1,l1           ! over m-1
                                lms1=l1+1+m1
                                ind1=l1+1+idt*m1+N1*(is1-1)
                                cc=0
                                sm1=0
                                sm2=0
                                do irf1=1,nrf
                                   do irf2=1,nrf
                                      ff = alml(l1,lms1,num1+nemin-1,irf1,is1)*conjg(alml(l1,lms1,num1+nemin-1,irf2,is1))
                                      sm1 = sm1 + ff * ri_mat(irf1,irf2,l1,is1)
                                      sm2 = sm2 + ff * ri_mat(irf1,1,l1,is1) * ri_mat(1,irf2,l1,is1)
                                   enddo
                                   cc = cc + alml(l1,lms1,num1+nemin-1,irf1,is1)*ri_mat(irf1,1,l1,is1)*(idt)**m1
                                enddo
                                cc = cc * sqrt(abs(sm1/sm2))
                                if (idt.gt.0) then
                                   URx(num1,ind1) = cc
                                else
                                   URx(num1,ind1) = conjg(cc)
                                endif
                             enddo
                          enddo
                       enddo
                       !$OMP END PARALLEL DO
                    elseif (abs(projector).eq.3) then
                       ! For correlated orbitals, we will also create transformation matrix U, such that P(LL'ij)=U^\dagger*U
                       !$OMP PARALLEL DO SHARED(URx) PRIVATE(num1,is1,m1,lms1,ind1,cc,corc,irf1,irf2) SCHEDULE(STATIC)
                       do num1=1,nbands      !nemin,nemax     ! over bands
                          do is1=1,iso            ! over spin-1
                             do m1=-l1,l1           ! over m-1
                                lms1=l1+1+m1
                                ind1=l1+1+idt*m1+N1*(is1-1)
                                irf1 = 1
                                cc = alml(l1,lms1,num1+nemin-1,irf1,is1)*sqrt(ri_mat(irf1,irf1,l1,is1))*(idt)**m1
                                ! Adding correction due to \cdot{u} and local orbitals
                                ! The results is  A(i,L,irf=1)*sqrt(o(L)) * sqrt( 1 + \sum_{irf>1} |A(i,L,irf)|^2*o(L,irf))
                                if ( abs(cc)>1e-6 ) then
                                   corc = 0
                                   do irf1=1,nrf
                                      do irf2=1,nrf
                                         if (irf1.eq.1 .and. irf2.eq.1) cycle
                                         corc = corc + alml(l1,lms1,num1+nemin-1,irf1,is1)*conjg(alml(l1,lms1,num1+nemin-1,irf2,is1))*ri_mat(irf1,irf2,l1,is1)
                                      enddo
                                   enddo
                                   cc = cc*sqrt(1 + abs(corc)/abs(cc)**2)
                                endif
                                if (idt.gt.0) then
                                   URx(num1,ind1) = cc
                                else
                                   URx(num1,ind1) = conjg(cc)
                                endif
                             enddo
                          enddo
                       enddo
                       !$OMP END PARALLEL DO
                    elseif (abs(projector).eq.5) then  ! fixed projector
                       ii = al_ucase(icase,lcase)
                       !$OMP PARALLEL DO SHARED(URx,rix_mat) PRIVATE(num1,is1,m1,lms1,ind1,cc,irf1) SCHEDULE(STATIC)
                       do num1=1,nbands  
                          do is1=1,iso   
                             do m1=-l1,l1 
                                lms1=l1+1+m1
                                ind1=l1+1+idt*m1+N1*(is1-1)
                                cc=0
                                do irf1=1,nrf
                                   cc = cc + alml(l1,lms1,num1+nemin-1,irf1,is1)*rix_mat(irf1,ii)*(idt)**m1
                                enddo
                                !!! New
                                cc = cc + al_interstitial(lms1,num1,is1,lcase)*(idt)**m1
                                !!! New
                                if (idt.gt.0) then
                                   URx(num1,ind1) = cc
                                else
                                   URx(num1,ind1) = conjg(cc)
                                endif
                             enddo
                          enddo
                       enddo
                       !$OMP END PARALLEL DO
                    else
                       print *, 'Only projector=[1,2,3,4,5,-1,-2,-3,-4,-5] is allowed!'
                       stop
                    endif

                    nind1 = nind
                    iorb1=iorb
                    
                    !if (iso.eq.2) then
                    !   no1 = nind1/iso
                    !   ALLOCATE( Uu(nbands,no1), Ud(nbands,no1))
                    !   tuu = Rspin(1,1,iorb1)
                    !   tud = Rspin(1,2,iorb1)
                    !   tdu = Rspin(2,1,iorb1)
                    !   tdd = Rspin(2,2,iorb1)
                    !   Uu(:,:) = URx(:,1:no1)
                    !   Ud(:,:) = URx(:,no1+1:2*no1)
                    !   URx(:,1:no1)       = Uu*tuu + Ud*tdu
                    !   URx(:,no1+1:2*no1) = Uu*tud + Ud*tdd
                    !   DEALLOCATE( Uu, Ud )
                    !endif
                    
                    DO iorb2=1,norbitals
                       if ( cix_orb(iorb2).NE. icix ) CYCLE
                       nind2 = nindo(iorb2)
                       call zgemm('N','C', nbands, nind2, nind1, (1.d0,0.d0), URx,nbands, cfX(:,:,iorb2,iorb1),maxdim2, (0.d0,0.d0), tmp,nbands)
                       DMFTU(:,:nind2,iorb2) = DMFTU(:,:nind2,iorb2) + conjg(tmp(:,:nind2))
                    ENDDO
                    
                 endif
              enddo  ! over L-case
           ENDDO     ! over atom-cases
           !-------- computing correlated green's function --------!

           Olapmk=0
           DO iorb1=1,norbitals
              icix = cix_orb(iorb1)
              if ( icix.EQ.0 ) CYCLE
              nind1 = nindo(iorb1)
              DO iorb2=1,norbitals
                 if ( icix.NE.cix_orb(iorb2) ) CYCLE
                 nind2 = nindo(iorb2)
        
                 call zgemm('C','N', nind1, nind2, nbands, (1.d0,0.d0), DMFTU(:,:,iorb1),nbands, DMFTU(:,:,iorb2),nbands, (0.d0,0.d0), olp,maxdim2)
                 
                 do ind1=1,nind1
                    do ind2=1,nind2
                       Olapmk( iSx(ind1,iorb1), iSx(ind2,iorb2), icix ) = olp(ind1,ind2)
                    enddo
                 enddo
              enddo
           ENDDO
           Olapm(:,:,:) = Olapm(:,:,:) + Olapmk(:,:,:) * (mweight(ikp)/tweight/nsymop2) ! proper weight for the reducible k-point
        ENDDO !------  isym: over star of the irreducible k-point ----!
     
        deallocate( DMFTU)
        deallocate( URx, tmp)

        if (abs(projector).eq.5) then
           deallocate( phi_jl )
           deallocate( h_interstitial )
           deallocate( a_interstitial, al_interstitial )
        endif
        
        WRITE(*,'(I3,A,1x,I5,1x,A,I4,A,I2)') myrank, ') Finished k-point number', ikp, 'with #bands=', nbands, ' Nstar=', nsymop2
        call flush(6)
     ENDDO !---- over reducible k-point: ikp
  
998  CONTINUE !---- irreducible k-points end (jump) from reading
     !--- end k-points

     DO is=1,iso    !------ over up/dn ---------------------!
        itape=8+is
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     ENDDO
     
  ENDDO ! over different vector files


  CALL AllReduce_MPI(Olapm, maxdim, ncix)
  if (myrank.eq.master) then
     WRITE(*,*) 'Renormalizing Gloc to account for the interstitials'
  endif

  if (Qsymmetrize) Call SymmetrizeOverlap(Olapm, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2)
  
  if (SIMPLE) then
     Olapc=0
     DO icix=1,ncix
        cixdm = cixdim(icix)
        !---- Olap to vector form, and s_oo to matrix form
        DO ip=1,cixdm
           do iq=1,cixdm
              it = abs(Sigind(ip,iq,icix))
              if (it.gt.0) then
                 Olapc(it, icix) =  Olapc(it, icix) + real(Olapm(ip,iq,icix))
              endif
           enddo
        ENDDO
        DO it=1,csize(icix)
           Olapc(it, icix) = Olapc(it, icix)/noccur(it,icix)
        ENDDO
     ENDDO
     if (myrank.eq.master) then
        WRITE(*,'(A)',advance='no') '  Z due to finite E-window='
        do icix=1,ncix
           do it=1,csize(icix)
              WRITE(*,'(F15.12)',advance='no') Olapc(it,icix)
           enddo
        enddo
        WRITE(*,*)
     endif
  else if (.True.) then
     !!! We copute SOlapm = 1/sqrt(Olapm) for later use
     lwork  = 4*maxdim + maxdim*maxdim
     lrwork = 10*maxdim + 2*maxdim*maxdim + 1
     liwork = 3+10*maxdim
     allocate( work(lwork), rwork(lrwork), iwork(liwork) )

     SOlapm(:,:,:) = 0.d0
     DO icix=1,ncix
        cixdm = cixdim(icix)
        
        allocate( cind(cixdm) )
        cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated                                                                                    
        cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block                                                                                         
        DO ip=1,cixdm
           if (Sigind(ip,ip,icix) .ne. 0) then
              cixdms = cixdms + 1
              cind(ip) = cixdms
           endif
        ENDDO
        allocate( cini(cixdms))
        do ip=1,cixdm
           if (cind(ip).gt.0) cini(cind(ip))=ip
        enddo
        deallocate( cind )

        allocate( olocf(cixdms,cixdms), ws(cixdms), tmp1(cixdms,cixdms), solap(cixdms,cixdms) )
        do ip=1,cixdms
           do iq=1,cixdms
              olocf(ip,iq) = Olapm(cini(ip),cini(iq),icix)
           enddo
        enddo
        CALL ZHEEVD('V','U', cixdms, olocf, cixdms, ws, work, lwork, rwork, lrwork, iwork, liwork, info )
        if (info .ne. 0) then
           print *, 'Diagonalization of renormalization factor failed. Info-zheevd=', info
        endif
        do ip=1,cixdms
           tmp1(:,ip) = olocf(:,ip)*(1./sqrt(abs(ws(ip))))
        enddo
        call zgemm('N','C', cixdms, cixdms, cixdms, (1.d0,0.d0), tmp1, cixdms, olocf, cixdms, (0.d0,0.d0), solap, cixdms)


        do ip=1,cixdm
           SOlapm(ip,ip,icix)=1.0d0
        enddo
        do ip=1,cixdms
           do iq=1,cixdms
              SOlapm(cini(ip),cini(iq),icix) = solap(ip,iq)
           enddo
        enddo
        
        if (myrank.eq.master) then
           !WRITE(*,'(A)',advance='no') '  Z due to finite E-window='
           !do it=1,cixdm
           !   WRITE(*,'(F15.12)',advance='no') ws(it)
           !enddo
           !WRITE(*,*)
           WRITE(*,'(A)',advance='no') 'Z to renormalize='
           do ip=1,cixdms
              WRITE(*,'(F16.10)',advance='no') 1./dble(SOlapm(cini(ip),cini(ip),icix))
           enddo
           WRITE(*,*)
        endif
        
        deallocate( cini )
        deallocate( olocf, ws, tmp1, solap )
     ENDDO
     
     deallocate( work, rwork, iwork )

     ! Save for later for dmft2 step
     open(996, FILE='SOlapm.dat', status='unknown')
     DO icix=1,ncix
        cixdm = cixdim(icix)
        WRITE(996,*) icix
        do ind1=1,cixdm
           do ind2=1,cixdm
              WRITE(996,'(2F20.16)',advance='no')  SOlapm(ind1,ind2,icix)
           enddo
           WRITE(996,*)
        enddo
     ENDDO
     close(996)
     
  else
     !!! We copute SOlapm = 1/sqrt(Olapm) for later use
     lwork  = 4*maxdim + maxdim*maxdim
     lrwork = 10*maxdim + 2*maxdim*maxdim + 1
     liwork = 3+10*maxdim
     allocate( work(lwork), rwork(lrwork), iwork(liwork) )
     allocate( olocf(maxdim,maxdim), ws(maxdim), tmp1(maxdim,maxdim) )
     
     DO iorb1=1,norbitals
        icix = cix_orb(iorb1)
        if ( icix.EQ.0 ) CYCLE
        nind1 = nindo(iorb1)
        olocf=0
        do ind1=1,nind1
           do ind2=1,nind1
              olocf(ind1,ind2) = Olapm( iSx(ind1,iorb1), iSx(ind2,iorb1), icix )
           enddo
        enddo
        
        CALL ZHEEVD('V','U', nind1, olocf, maxdim, ws, work, lwork, rwork, lrwork, iwork, liwork, info )
        if (info .ne. 0) then
           print *, 'Diagonalization of renormalization factor failed. Info-zheevd=', info
        endif
        
        do ip=1,nind1
           tmp1(:,ip) = olocf(:,ip)*(1./sqrt(abs(ws(ip))))
        enddo
        call zgemm('N','C', nind1, nind1, nind1, (1.d0,0.d0), tmp1, maxdim, olocf,maxdim, (0.d0,0.d0), SOlapm(:,:,icix),maxdim)

        if (myrank.eq.master) then
           WRITE(*,*) 'Z due to finite E-window=', ws(:nind1)
        endif
        
     ENDDO
     deallocate( work, rwork, iwork )
     deallocate( olocf, ws, tmp1 )
  endif


  DEALLOCATE( a_real ) !---  for eigenvectors -------------!
  deallocate( olp )
  deallocate( Olapmk )

  !CALL cputim(t2c)
  !CALL walltim(t2w)
  !time1c = time1c + t2c-t1c
  !time1w = time1w + t2w-t1w

  !print *, 'total time overlap:', time1c, time1w
  !print *, 'overlap interst. proj:', time2c, time2w
  
END SUBROUTINE cmp_overlap


