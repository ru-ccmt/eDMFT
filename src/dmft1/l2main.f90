SUBROUTINE L2MAIN(Qcomplex,nsymop,mode,projector,Qrenormalize,fUdmft)
  !------------------------------------------------------------------------------
  !-- Main routine for setting up DMFT transformation to local basis set       --
  !-- The key object, DMFTransform can transform self-energy to Kohn-Sham base --
  !-- as well as project the Kohn-Sham green's function to localized base      --
  !------------------------------------------------------------------------------  
  !--  mode: 'g' -> computes the green's function
  !--        'e' -> eigenvalues are printed
  !--        'u' -> printing transformation
  !--        'x' -> computes xqtl2, just for debugging
  USE param,    ONLY: LMAX2, LOMAX, NRAD, nloat, nrf, nmat, nume, gamma, gammac, aom_default, bom_default, nom_default, IBLOCK, matsubara, Cohfacts, ComputeLogGloc, maxbands, nemin0, nemax0, max_nl, cmp_partial_dos, Hrm, lmaxr
  USE structure,ONLY: VOL, RMT, iord, iz, nat, mult, pos, tau, rotij, tauij, jri, BR1, rot_spin_quantization
  USE case,     ONLY: cf, csize, maxsize, nl, cix, ll, Sigind, iatom, legend, maxdim, natom, ncix, shft, isort, crotloc, ifirst
  USE sym2,     ONLY: tmat, idet, iz_cartesian
  USE com,      ONLY: MINWAV, MAXWAV, iso, iso2, emin, emax, EF
  USE abc,      ONLY: KX, KY, KZ, BK3, E, A, ALM, ALML, aK
  USE kpts,     ONLY: mweight, tweight
  USE w_atpar,  ONLY: w_alo, w_nlo, w_nlov, w_nlon, w_ilo, w_lapw, w_ri_mat, w_P, w_DP, w_jatom, w_FJ, w_DFJ, w_allocate, w_deallocate, w_rfk, w_allocate_rfk, w_deallocate_rfk, nnlo, el_store, elo_store
  USE p_project,ONLY: p_allocate, p_deallocate, p_cmp_rixmat, P_rfi, n_al_ucase, max_lcase, al_ucase, al_ucase_, l_al_ucase, j_al_ucase, nl_case, kNri, rix, dri, al_interstitial
  USE com_mpi,  ONLY: nprocs, myrank, master, cpuID, vectors, vector_para, VECFN, fvectors, FilenameMPI, FilenameMPI2, Gather_MPI, Reduce_MPI, FindMax_MPI, Gather_procs, nvector, Reduce_DM_MPI, Qprint
  USE kpts,     only: numkpt
  USE matpar,   only: atpar, alo, nlo, nlov, nlon, ilo, lapw, RI_MAT, P, DP, RF1, RF2
  use DMFTProjector, only : CmpDMFTrans
  use PrintG, only : PrintGloc
  IMPLICIT NONE
  !------- Input variables
  LOGICAL, intent(in)      :: Qcomplex, Qrenormalize
  INTEGER, intent(in)      :: nsymop, projector
  CHARACTER*1, intent(in)  :: mode
  CHARACTER*200, intent(in):: fUdmft
  !------ External Functions
  INTEGER :: CountSelfenergy
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
  !------- allocatable local varaibles
  real*8, PARAMETER       :: Ry2eV = 13.60569193d0
  complex*16, allocatable :: DMFTrans(:,:,:,:,:), STrans(:,:,:,:), GTrans(:,:,:), DMFTU(:,:,:)
  complex*16, allocatable :: xqtl2(:,:,:,:), gk(:), gij(:,:), gloc(:,:), gtot(:), gmk(:,:,:), gmloc(:,:,:,:), Ekp(:,:,:,:)
  complex*16, allocatable :: g_inf(:,:,:,:), g_ferm(:,:,:)!, DM(:,:,:)
  complex*16, allocatable :: sigma(:,:,:), Deltac(:,:,:), Glc(:,:,:), zek(:), s_oo(:,:)
  real*8,     allocatable :: a_real(:), omega(:), Eimpc(:,:), Olapc(:,:), dek(:)
  real*8,     allocatable :: Cohf0(:)
  complex*16, allocatable :: Eimpm(:,:,:), Eimpmk(:,:,:), Cohf(:,:), evl(:,:), evr(:,:)
  complex*16, allocatable :: Olapm(:,:,:), Olapmk(:,:,:), Olapm0(:,:,:), SOlapm(:,:,:)
  integer,    allocatable :: noccur(:,:), nbandsk(:), nemink(:), n_ik(:), tn_ik(:), korder(:)
  integer,    allocatable :: cixdim(:), iSx(:,:), nindo(:), csort(:), cix_orb(:)
  complex*16, allocatable :: tEk(:,:,:,:), cfX(:,:,:,:)!, Rspin(:,:,:)
  integer,    allocatable :: tnbands(:), tnemin(:)     
  INTEGER,    allocatable :: pr_procs(:)
  complex*16, allocatable :: A2(:,:)
  ! For p_interstitial
  COMPLEX*16, allocatable :: a_interstitial(:,:,:), h_interstitial(:,:,:)
  REAL*8, allocatable     :: phi_jl(:,:), crotloc_x_rotij(:,:,:,:)
  REAL*8, allocatable     :: aKR(:), jlr(:), jlp(:)
  !------- other local arrays
  complex*16  :: h_yl(2*LMAX2+1,iblock)
  complex*16  :: h_alyl(2*LMAX2+1,iblock,2)
  complex*16  :: h_blyl(2*LMAX2+1,iblock,2)
  complex*16  :: YL((LMAX2+1)*(LMAX2+1))
  integer     :: iorbital(natom,lmax2+1), csizes(ncix)
  real*8      :: BK(3), BKROT(3), BKROT2(3), BKRLOC(3), el_read(lmaxr)
  !------ other local variables
  logical     :: Tcompute, read_overlap, Qsymmetrize
  character*10:: KNAME
  complex*16  :: PHSHEL, CFAC, xomega
  real*8      :: PI, TWOPI, EMIST, S, T, Z, exxx, ARG123, ARG2, ARGT, ARGT2, FAC, logGloc(ncix), Kn(3), rx, mweight_ikp
  real*8      :: crotloc_x_BR1(3,3)
  integer     :: norbitals, cnemin, cnbands
  integer     :: l1, idt, num, NUM_MAX, iscf, maxdim2
  integer     :: nomega, icix
  integer     :: i, ii, jj, j, iom, icase, lcase, is, itape, jtape, jatom, i3, lda, ldb, ldc, irf, iind, iat
  integer     :: iorb, L, nind, it, ikp, iks, iikp, ivector, N, NE, NEMIN, NEMAX, nbands, isym, igi, M, ibb, LATOM
  integer     :: fh_sig, fh_dos, fh_gc, fh_dt, fh_eig, fh_Eimp
  integer     :: iucase, pr_proc, pr_procr, max_bands, natm, maxucase
  integer     :: fh_p, ind, lfirst, nkp, il, nr0, ir, Nri, isize, info, iband
  character*100:: filename
  logical     :: pform, nonzero_shft
  CHARACTER*200 :: FNAME
  !
  REAL*8      :: Det, phi1, the1, psi1, WGH
  COMPLEX*16  :: Rispin(2,2)
  complex*16  :: IMAG, gtc
  LOGICAL    :: RENORM_SIMPLE, cmp_DM
  INTEGER    :: Nstar
  REAL*8     :: k_star(3,iord), rotij_cartesian(3,3), BR1inv(3,3), tmp3(3,3), Trans3(3,3)
  INTEGER    :: gind(iord), nsymop2,  iom_zero
  REAL*8     :: time0, time1, time2, time3, read_time, read_time2, atpar_time, overlap_time, bess_time, zgemm_time, zgem2_time, alm_time, trans_time, comprs_time, imp_time, inv_time, gc_time, dm_time
  LOGICAL    :: vector_out
  INTEGER    :: vcn_out(2)

  RENORM_SIMPLE=.False.
  cmp_DM = matsubara
  read_overlap = nsymop.ne.iord .and. mode.NE.'g'  ! We do not intent to go over all symmetry operations (like plotting), hence need overlap from before
  Qsymmetrize  = nsymop.ne.iord .and. mode.EQ.'g'   ! We also do not go over all symmetry operations, but we will symmetrize

!!!! Here you should set ComputeLogGloc to true if you want to compute it for free energy!!! 
  ComputeLogGloc = .False.
  !------------------------------------------------------------------
  DATA IMAG/(0.0D0,1.0D0)/
  !------------------------------------------------------------------     
  fh_sig = 80
  fh_dos = 500
  fh_gc  = 180
  fh_dt  = 280
  fh_eig = 380
  fh_Eimp = 380
  !------ Some parameters which should be read from the input! Not yet implemented! -----!
  !----------  Some old stuff ---------------------
  PI=ACOS(-1.0D0)
  TWOPI=2.D0*PI
  MAXWAV=0
  MINWAV=100000

  tmat(:,:,:iord) =  iz(:,:,:iord)

  call cputim(time0)

  natm = sum(mult) 
  ALLOCATE( csort(nat) )
  CALL Create_Atom_Arrays(csort, maxucase, isort, iatom, nat, natm, natom)

  !----------- Find index to atom/L named iorbital(icase,lcase) -------------------------------!
  CALL Create_Orbital_Arrays(iorbital, norbitals, maxdim2, nl, ll, cix, natom, lmax2, iso)

  allocate( noccur(maxsize,ncix) )
  ALLOCATE( cixdim(ncix), iSx(maxdim2, norbitals), nindo(norbitals), cix_orb(norbitals), cfX(maxdim2,maxdim2,norbitals,norbitals) )
  cfX=0

  CALL Create_Other_Arrays(cixdim, iSx, noccur, nindo, cix_orb, cfX, CF, nl, ll, cix, iorbital, csize, csizes, Sigind, iso, natom, maxdim, lmax2, ncix, maxsize, norbitals, maxdim2)

  !ALLOCATE( Rspin(2,2,norbitals) )
  !CALL GetSpinRotation(Rspin,rotij,crotloc,BR1,norbitals,natom,natm,iso,lmax2,iatom,nl,cix,ll,iorbital)


  filename = 'BasicArrays.dat'
  if (mode.EQ.'u') CALL PrintSomeArrays(filename, nat, iso, norbitals, ncix, natom, numkpt, nmat, nume, Qcomplex, lmax2, maxdim2, maxdim, maxsize, nindo, cixdim, nl, ll, cix, iorbital, csize, iSx, Sigind, EF, VOL)

  if (abs(projector).ge.5) then
     filename='projectorw.dat'
     call p_allocate(filename, maxucase,csort)
     call w_allocate_rfk(n_al_ucase) !------ contains all variables which depend on atom / which are changed by atpar -----------!
  endif

  call w_allocate(maxucase,nat,iso2) !------ contains all variables which depend on atom / which are changed by atpar -----------!

  vector_out = Hrm(1:2).eq.'HW' .and. .not.vector_para 
  if ( vector_out ) then
     do is=1,iso
        vcn_out(is) = 1973+is
        FNAME = TRIM(VECFN(is))//'_dmft'
        open(vcn_out(is),FILE=FNAME,STATUS='unknown',FORM='unformatted')
     enddo
  end if
  
  do is=1,iso2
     itape=8+is
     if (vector_para) then
        if (nvector.ge.1) then
           FNAME = fvectors(1,is)
        else   ! No k-point to compute, but still need linearization energies
           FNAME = TRIM(VECFN(1))//'_1'
        endif
        open(itape,FILE=FNAME,STATUS='old',FORM='unformatted')
     else
        !------------ vector files need to be read from the beginning -------------------!           
        rewind(itape) 
     endif

     DO iat=1,nat
        READ(itape) el_read(1:lmaxr)   ! linearization energies
        el_store(0:lmax2,iat,is) = el_read(1:lmax2+1)
        READ(itape) elo_store(0:lomax,1:nloat,iat,is) ! local
        if ( vector_out ) then  
           write(vcn_out(is)) (el_read(j),j=1,lmaxr)   ! linearization energies
           write(vcn_out(is)) ((elo_store(j,ii,iat,is),j=0,lomax),ii=1,nloat) ! local
           !print *, 'el=', el_read(1:lmaxr)
           !print *, 'elo=', elo_store(0:lomax,1:nloat,iat,is)
        endif
     ENDDO
  end do

  call cputim(time1)
  read_time = time1-time0

  time0=time1
  !-------------------------------------------------------------------------------------!
  !--------------- Reading potential parameters for all atoms.   -----------------------!
  !-- It stores  necessary parameters and uses them below, when looping over atoms -----!
  !-- The step is necessary because green's function of each orbital depends on the ----!
  !--  self-energy of all l's and all atoms. If more than one atom is correlated, we ---!
  !--  can not avoid looping over all correlated atoms in the inside loop --------------!
  !-------------------------------------------------------------------------------------!
  do is=1,iso2 !------ preparing input files case.vsp and case.vspn containing potential parameters -----!
     jtape=17+is
     rewind(jtape)
     READ(jtape,2032) ISCF
  end do
  DO jatom=1,nat    !------  over all atoms (all sorts) even if we do not require qtl for them ------------------!
     iucase = csort(jatom) !--- this particular atom is requested in the input under index iucase --!
     if(iucase.gt.0) then                 !--- requested in the inut ----!
        if (Qprint) write(6,'(A,1x,I2,1x,A,1x,I2)')' Calculation for atom', jatom, 'which has iucase=', iucase
        w_jatom(iucase) = jatom
     else
        if (Qprint) write(6,'(A,1x,I2,1x,A)')' ATOM=', jatom, ' left out ' !-- not requested in the input ----!
     endif
     !----  mult(jatom) stands for multiplicity of each atom -> number of equivalent atoms of this sort --------!
     !--- Reading radial functions U(R), UE(R), for all atoms including those that are left out ----------------!
     do is=1,iso2
        jtape=17+is
        CALL ATPAR(jatom,jtape,el_store(:,:,is),elo_store(:,:,:,is),is,iso2)
        if (vector_para) then
           close(itape)
        else
           rewind(itape) 
        endif
     end do
     nnlo=nlo+nlon+nlov  ! sum of all local orbitals. It does not depend on jatom

     if (iucase.gt.0) then ! save values of potential parameters
        if (Qprint) write(6,'(A,I2,A,I2)')'Saving potential parameters for atom',jatom, ' into ', iucase
        w_nlo(iucase)  = nlo
        w_nlov(iucase) = nlov
        w_nlon(iucase) = nlon
        !w_loor(:,iucase) = loor(:)
        w_ilo(:,iucase)  = ilo(:) 
        w_lapw(:,iucase) = lapw(:)
        w_alo(:,:,:,1,iucase) = alo(:,:,:,1)                 ! alo(l,jlo,irf,is)
        w_p(:,1,:,iucase) = P(:,1,:)                         ! P(lmax,is,irf)
        w_dp(:,1,:,iucase) = DP(:,1,:)                       ! DP(lmax,is,irf)
        w_ri_mat(:,:,:,1,iucase) = ri_mat(:,:,:,1)           ! ri_mat(nrf,nrf,0:lmax,is)

        if (abs(projector).ge.5) then
           nr0 = jri(jatom)
           do irf=1,nrf
              do il=1,nl_case(iucase)
                 ii = al_ucase_(iucase,il)
                 l  = l_al_ucase(ii)
                 w_rfk(:nr0,1,irf,ii) = rf1(:nr0,l,1,irf)
                 w_rfk(:nr0,2,irf,ii) = rf2(:nr0,l,1,irf)
              enddo
           enddo
        endif

        if (iso2.eq.2) then
           w_alo(:,:,:,2,iucase) = alo(:,:,:,2)                 ! alo(l,jlo,irf,is)
           w_p(:,2,:,iucase) = P(:,2,:)                         ! P(lmax,is,irf)
           w_dp(:,2,:,iucase) = DP(:,2,:)                       ! DP(lmax,is,irf)
           w_ri_mat(:,:,:,2,iucase) = ri_mat(:,:,:,2)           ! ri_mat(nrf,nrf,0:lmax,is)
        else if (iso.eq.2 .and. iso2.eq.1)then 
           w_alo(:,:,:,2,iucase) = alo(:,:,:,1)                 ! alo(l,jlo,irf,is)
           w_p(:,2,:,iucase) = P(:,1,:)                         ! P(lmax,is,irf)
           w_dp(:,2,:,iucase) = DP(:,1,:)                       ! DP(lmax,is,irf)
           w_ri_mat(:,:,:,2,iucase) = ri_mat(:,:,:,1)           ! ri_mat(nrf,nrf,0:lmax,is)
        endif
     endif
  ENDDO

  if (abs(projector).ge.5) then
     CALL p_cmp_rixmat(w_rfk)
     call w_deallocate_rfk()
  endif

  call cputim(time1)
  atpar_time = time1-time0

  !---------------- Reading self-energy ------------------------------------------------------------------!
  !---------------- Input self-energy must be in eV units. Will convert to Rydbergs below ----------------!
  !---------------- Resulting gloc written to file is also in eV units -----------------------------------!
  if (ncix.gt.0) then
     nomega = CountSelfenergy(fh_sig+1, csize(1))-1 !--- How many frequency point exists in the input file? ---!
  else                                           !--- For non-correlated case we do not have self-energy ---!
     nomega = nom_default                        !--- using default ----------------------------------------!
  endif

  allocate( sigma(maxsize,ncix,nomega), omega(nomega), s_oo(maxsize,ncix) ) !----- allocating necessary arrays for self-energy ---!

  if (ncix.eq.0) then
     do iom=1,nomega
        omega(iom) = aom_default + (bom_default-aom_default)*iom/(nomega-1.)     !-- default frequency mesh --!
     enddo
  endif

  !----------- Actual reading of the self-energy -------------!
  if (mode.eq.'g' .or. mode.eq.'e') then
     if (Qprint) then
        WRITE(6,*) '------------------------------------------------------------------------------'
        WRITE(6,*) '----------------- SELF ENERGY AND ITS STRUCTURE ------------------------------'
     endif
     do icix=1,ncix
        CALL ReadSelfenergy(fh_sig+icix, sigma(:,icix,:), omega, s_oo(:,icix), gammac, csize(icix), nomega, maxsize)
     enddo
     !--------------  change from eV to Rydbergs --------------------------!
     omega(:) = omega(:)/Ry2eV
     sigma(:,:,:) = sigma(:,:,:)/Ry2eV
     s_oo(:,:) = s_oo(:,:)/Ry2eV
  endif
  !---------------- Reading self-energy --------------------------------!


!!! ---------- Preparation of arrays for paralel executaion --------------
  if (vector_para) then
     pr_proc = sum(vectors(:,2))
     pr_procr = pr_proc
     !numk = pr_proc
     if (Qprint) WRITE(6,'(A,I3,2x,A,I3)') 'pr_proc=', pr_proc, 'tot-k=', numkpt
  else
     pr_proc  = floor(numkpt/DBLE(nprocs)+0.999)  ! The maximum number of points calculated per processor     
     ! correcting pr_proc is the number of k-points is not dividable
     pr_procr = pr_proc
     if ((myrank+1)*pr_proc .gt. numkpt) then
        if (numkpt-myrank*pr_proc.gt.0) then
           pr_procr = numkpt-myrank*pr_proc
        else
           pr_procr = 0
        endif
     endif
     if (Qprint) WRITE(6,'(A,I3,2x,A,I3)') 'pr_proc=', pr_proc, 'pr_procr=', pr_procr, 'tot-k=', numkpt
  endif


  ! frequency closer to zero
  iom_zero=1
  do iom=2,nomega
     if (abs(omega(iom)).LT.abs(omega(iom_zero))) then
        iom_zero = iom
     end if
  end do
  if (vector_out .and. myrank.eq.master) then
     WRITE(6,*) 'iom_zero=', iom_zero, 'omega(iom_zero)=', omega(iom_zero)
  endif
  
  ALLOCATE( pr_procs(nprocs) )
  CALL Gather_procs(pr_procr, pr_procs, nprocs)

  if (myrank.eq.master .AND. SUM(pr_procs).NE.numkpt) then
     WRITE(6,*) 'ERROR: sum(pr_procs) should be numkpt, but is not', sum(pr_procs), numkpt
  endif

  !if (iso.eq.2) then
  allocate( crotloc_x_rotij(3,3,max_nl,natom) )
  !call INVERSSYMDEF(BR1,BR1inv)
  call inv_3x3(BR1,BR1inv)
  DO icase=1,natom  !--------------- over all atoms requested in the input ------------------------!
     latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
     tmp3 = matmul(BR1,rotij(:,:,latom))
     rotij_cartesian = matmul(tmp3,BR1inv)
     do il=1,nL(icase)
        crotloc_x_rotij(:,:,il,icase) = matmul(crotloc(:,:,il,icase),rotij_cartesian)
     end do
  END DO
  !endif

  call cputim(time0)
  if (Qrenormalize .and. abs(projector).le.5) then
     allocate(Olapm0(maxdim,maxdim,ncix), SOlapm(maxdim,maxdim,ncix) )
     if (read_overlap) then
        call read_overlap_from_file(info, SOlapm, cixdim, maxdim, ncix)
        if (info.ne.0) then
           WRITE(6,*) 'ERROR : file SOlapm.dat could not be found. Run regular dmft1 first to produce this file.'
           STOP
        endif
     else
        CALL cmp_overlap(projector,Olapm0, SOlapm, Qcomplex, nsymop, csort, iorbital, cix_orb, nindo, cixdim, iSx, noccur, cfX, crotloc_x_rotij, maxucase, maxdim2, norbitals, pr_proc, RENORM_SIMPLE, Qsymmetrize)
     endif
  endif
  call cputim(time1)
  overlap_time = time1-time0

  !------------ allocating some important arrays --------------!
  ALLOCATE( a_real(nmat) ) !---  for eigenvectors -------------!


  if (mode.EQ.'g') then
     allocate( gmloc(maxdim, maxdim, ncix, nomega) )
     gmloc=0
     allocate(Olapm(maxdim,maxdim,ncix), Olapmk(maxdim,maxdim,ncix))
     Olapm=0
     Olapmk=0
     allocate(Eimpm(maxdim,maxdim,ncix), Eimpmk(maxdim,maxdim,ncix))
     Eimpm=0
     Eimpmk=0
     if (cmp_partial_dos) then
        allocate( gloc(norbitals,nomega) ) !--- local green's function and momentum dependent green's function --!
        gloc=0
     endif
     allocate(gtot(nomega))
     gtot=0
     if (cmp_DM) then
        allocate( g_inf(maxdim,maxdim,ncix,nomega) )
        allocate( g_ferm(maxdim,maxdim,ncix) )
        g_inf(:,:,:,:)=0.d0
        g_ferm(:,:,:)=0.d0
     endif
  elseif (mode.EQ.'e') then
     CALL FindMax_MPI(max_bands, maxbands)
     !print *, 'Maximum bands over all cores=', max_bands
     allocate(Ekp(nomega,nsymop,max_bands,pr_procr), nbandsk(pr_procr), nemink(pr_procr), n_ik(pr_procr))
     Ekp=0
     nbandsk=0
     nemink=0
     n_ik=0
  elseif (mode.EQ.'u') then
     print *, 'cpuID=', cpuID
     filename = TRIM(fUdmft)//ADJUSTL(TRIM(cpuID))  ! fUdmft comes from def file
     fh_p = 501
     pform = .False.
     if (pform) then
        open(fh_p,file=TRIM(filename),status='unknown', form='formatted')
     else 
        open(fh_p,file=TRIM(filename),status='unknown', form='unformatted')
     endif


     if (pform) then
        WRITE(fh_p,*) pr_procr, nsymop, norbitals
        WRITE(fh_p,*) (nindo(iorb), iorb=1,norbitals)
     else
        WRITE(fh_p) pr_procr, nsymop, norbitals
        WRITE(fh_p) (nindo(iorb), iorb=1,norbitals)
     endif

  endif
  if (Cohfacts) ALLOCATE( Cohf0(nume) )

  bess_time=0
  zgemm_time=0
  alm_time=0
  trans_time=0
  comprs_time=0
  imp_time=0
  inv_time=0
  gc_time=0
  dm_time=0
!!! start-kpoints
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

        if (Tcompute .and. abs(projector).ge.5) allocate( phi_jl(nmat, n_al_ucase) )

        read_time2=0
        call cputim(time0)
        !------- reading from vector for both spins -----------!
        DO is=1,iso    !------ over up/dn ---------------------!
           call cputim(time2)
           itape=8+is
           READ(itape,END=998) S,T,Z,KNAME,N,NE,WGH  !--- Here we can jupm out of loop 4 and stop at 998 -----------------------------------!
           IF(N.GT.MAXWAV) MAXWAV=N                                         
           IF(N.LT.MINWAV) MINWAV=N
           READ(itape) (KX(I),KY(I),KZ(I),I=1,N) !--- Reads all reciprocal vectors -----------------------------------------------------!

           if (vector_out) then
              write(vcn_out(is)) S,T,Z,KNAME,N,NE,WGH
              write(vcn_out(is)) (KX(I),KY(I),KZ(I),I=1,N)
           endif

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
           NUM_MAX=NUM
           call cputim(time3)
           read_time2 = read_time2 + time3-time2

           if (abs(projector).ge.4) then
              nemin=nemin0
              nemax=min(nemax0,NE)
           endif

           IF (Tcompute .and. is.eq.1) THEN         !----- No need to compute this for both spins ----!
              isize=N-nnlo
              DO I=1,N                              !----- Over reciprocal vectors --------------------------------------------------------!
                 BK3(1,I)=(S+KX(I))                 !----  Creates K+k where K is reciprocal vector and k is irreducible k-point ----------!
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
                 ! You could optimize here and compute up to N-nnlo only. Because the rest is for local orbitals.
              ENDDO
              ! New
              if (abs(projector).ge.5) then
                 Nri=2**kNri+1    ! Number of radial points in the interstitials
                 phi_jl(:,:)=0.d0
                 allocate( aKR(Nri), jlr(Nri), jlp(Nri) )
                 DO iind=1,n_al_ucase
                    if (abs(dri(iind)).gt.1e-10) then
                       l     = l_al_ucase(iind)
                       jatom = j_al_ucase(iind)
                       DO i=1,isize         ! over all reciprocal vectors K. You could go up to N-nnlo only.
                          rx = RMT(jatom)
                          do ir=1,Nri
                             aKR(ir)=rx*aK(i)    ! |k+K|*r
                             rix(ir)=rx          !  r
                             rx = rx + dri(iind)
                          enddo
                          CALL sphbes2(l,Nri,aKR,jlr)  ! spherical bessel : j_l(|k+K|*r)
                          jlr(:) = jlr(:)*rix(:)       ! r*j_l . We need to do that, because the integral is Int[ (phi(r)/r)*j_l(r)*r^2]=Int[phi(r)*j_l(r)*r]
                          jlp(:) = jlr(:)*P_rfi(:,iind)  ! projector from file projectorw.dat
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
        if (cmp_partial_dos) allocate( DMFTrans(nbands,nbands,maxdim2,maxdim2,norbitals) )
        allocate( DMFTU(nbands,maxdim2,norbitals) )

        if (abs(projector).ge.5) then
           allocate( h_interstitial(2*LMAX2+1,iblock,2) )
           allocate( a_interstitial(2*LMAX2+1,nbands,iso), al_interstitial(2*LMAX2+1,nbands,iso,max_lcase) )
        endif

        if (mode.EQ.'e') then
           nbandsk(iikp) = nbands
           nemink(iikp) = NEMIN
           n_ik(iikp) = ikp
        elseif (mode.EQ.'x') then
           allocate( xqtl2(nbands,maxdim2,maxdim2,norbitals) )
           xqtl2=0
        elseif (mode.EQ.'u') then
           if (pform) then
              WRITE(fh_p,*) ikp, nbands,maxdim2,norbitals, nemin
           else
              WRITE(fh_p) ikp, nbands,maxdim2,norbitals, nemin
           endif
        endif
        call cputim(time1)
        bess_time = bess_time + time1-time0 - read_time2
        read_time = read_time + read_time2
        time0=time1

        if (Cohfacts) then
           ALLOCATE( Cohf(nbands,nbands) )
           READ(1001,*) jj, cnemin, cnbands
           DO i=1,cnbands
              READ(1001,*) (Cohf0(j),j=1,cnbands)
              if (i+cnemin-1.ge.nemin .and. i+cnemin-1.le.nemax) then
                 Cohf(i+cnemin-nemin,:) = Cohf0(nemin-cnemin+1:nemax-cnemin+1)
              endif
           ENDDO
           OPEN(1002,FILE='cohfactorsd.dat',status='unknown')
           call cputim(time1)
           read_time = read_time + time1-time0
        endif


        allocate( STrans(maxsize,ncix,nbands,nbands) )  ! WWW
        if (cmp_partial_dos) allocate( GTrans(nbands,nbands,norbitals) )

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
        !----------- sum over all k-points in the star of the irreducible k-point ---------------------------------!
        DO igi=1,nsymop2  !-- Over all symmetry operarations -> all k-points in the start of irreducible k point --!
           isym = gind(igi)
           DMFTU=0
           DO icase=1,natom  !--------------- over all atoms requested in the input ------------------------!
              call cputim(time0)
              zgem2_time=0
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
              if (abs(projector).ge.5) al_interstitial(:,:,:,:)=0.d0

              nonzero_shft = sum(abs(shft(:3,latom))) .GT. 1e-10
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
                 if (abs(projector).ge.5) a_interstitial=0
                 !--------- blocks are for efficiency. Matrix is multiplied in block form. This must be important in the past, while modern BLAS should do that better. I think it is obsolete.
                 DO ii=1,isize,iblock !------ iblock is 128 for 32-bit system -------!
                    !-------- nlo-number of local orbitals -----!
                    do i=ii,min(ii+iblock-1,isize)  ! 121
                       !---------  rotates ylm(k+K) to ylm(k'+K) where k' is in star of irreducible k. ------------!
                       ! BK=BK3(:,i) -----  reciprocal vector and irreducible vector: G=K+k ----!
                       ! BKROT = R_a.(k+K) transforms to the reducible k-point
                       BKROT = matmul(TMAT(:,:,isym), BK3(:,I))
                       ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom 
                       BKROT2 = matmul(rotij(:,:,latom), BKROT)
                       !---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
                       ! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
                       !BKROT3 = matmul(BR1, BKROT2)
                       !!---- BKRLOC = crotloc.R_n.R_a.(k+K),  rotates according to the user specified local coordinate system.
                       !BKRLOC = matmul(crotloc(:,:,icase), BKROT3)
                       BKRLOC = matmul(crotloc_x_BR1, BKROT2)
                       !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
                       CALL YLM (BKRLOC,LMAX2,YL)  ! 
                       ! (R_n.R_a.(k+K)) *  R(first) * 2pi
                       ARG123 = dot_product(BKROT2, POS(:,lfirst))*TWOPI
                       ! ARGT = (k+K)*tau(isym) * 2pi
                       ARGT = dot_product(BK3(:3,I), TAU(:3,isym))*TWOPI
                       ! ARGT2 = (R_a.(k+K)).tau_n * 2pi
                       ARGT2= dot_product(BKROT, tauij(:3,latom))*TWOPI
                       ! ARG2 = (R_a.(k+K)) *  shft * 2pi
                       ARG2 = 0.d0
                       if (nonzero_shft) ARG2 = dot_product(BKROT,shft(:3,latom))*TWOPI ! Before user rotation, but already on 
                       ! PHSEHL = e^{I*2pi*( (R_a.(k+K))*tau_n + (K+k)*tau(isym) + (R_n.R_a.(k+K)*R(first)))}
                       PHSHEL=EXP(IMAG*(ARG123+ARG2+ARGT+ARGT2))
                       i3=i-ii+1
                       DO  m=1,2*l+1
                          h_yl(m,i3)=conjg(yl(l*l+m))*phshel !----- h_yl is rotated yl when k is rotated to k' -----!
                       END DO
                    enddo

                    DO is=1,iso  !--- over both spins
                       do i=ii,min(ii+iblock-1,isize)
                          i3 = i-ii+1
                          DO M=1,2*L+1
                             if (lapw(l)) then
                                ! P(l,is,1) = ul(Rmt)
                                ! P(l,is,2) = dot{ul}(Rmt)
                                ! DP(l,is,1) = d/dr ul(Rmt)
                                ! DP(l,is,2) = d/dr dot{ul}(Rmt)
                                ! FJ(l,G,iucase) = jl(|k+G|Rmt)
                                ! DFJ(l,G,iucase)=d/dr jl(|k+G|Rmt)
                                h_ALYL(m,i3,is)=(w_DFJ(L,I,iucase)*P(l,is,2)-w_FJ(L,I,iucase)*DP(l,is,2))* h_yl(M,i3) ! derivatives of bessel functions and spheric harmonics
                                h_BLYL(m,i3,is)=(w_FJ(L,I,iucase)*DP(l,is,1)-w_DFJ(L,I,iucase)*P(l,is,1))* h_yl(M,i3)
                             else
                                h_ALYL(m,i3,is)=w_FJ(L,I,iucase)/P(l,is,1)/RMT(jatom)**2*h_yl(M,i3)
                                h_BLYL(m,i3,is)=(0.d0,0.d0)
                             end if
                          END DO
                       enddo
!!! New
                       if (abs(projector).ge.5) then
                          iind=al_ucase(icase,lcase)
                          do i=ii,min(ii+iblock-1,isize)
                             DO M=1,2*L+1
                                !    h_interstitial = <P_phi|j_l> Y_L^*(R(k+K))*exp(i*(k+K)*r_latom)*exp(i*spin_phase)/Rmt^2
                                h_interstitial(M,i-ii+1,is) = phi_jl(i,iind)*h_yl(M,i-ii+1)/RMT(jatom)**2
                             ENDDO
                          enddo
                       endif
!!! New
                       ibb=min(iblock,isize-ii+1)
                       lda=2*LMAX2+1
                       ldc=lda
                       ldb=nmat

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
                       call cputim(time2)
                       call zgemm('N','N',2*l+1,nemax-nemin+1,ibb,(1.d0,0.d0), h_alyl(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), alm(1,nemin,1,is),ldc)
                       call zgemm('N','N',2*l+1,nemax-nemin+1,ibb,(1.d0,0.d0), h_blyl(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), alm(1,nemin,2,is),ldc)
                       if (abs(projector).ge.5) then !! The new
                          call zgemm('N','N',2*l+1,nbands,ibb,(1.d0,0.d0),h_interstitial(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), a_interstitial(1,1,is),ldc)
                       endif !! The new
                       call cputim(time3)
                       zgem2_time = zgem2_time + time3-time2
                    ENDDO !----------- over both spins ---------!
                 ENDDO    !----------- over iblock -------------!

                 !-------------- Adds localized orbitals to alm. -------------------!
                 if (nlo.ne.0) then
                    call lomain (nemin,nemax,lfirst,latom,n,jatom,isym,L,iso,crotloc_x_BR1)
                 end if

                 CFAC = FAC*IMAG**L      !------  (i)^l*fac -----!
                 if (iso.eq.2) then
                    ALML(l,:(2*L+1),nemin:nemax,:nrf,1) = ALM(:(2*L+1),nemin:nemax,:nrf,1)*(CFAC*Rispin(1,1)) + ALM(:(2*L+1),nemin:nemax,:nrf,2)*(CFAC*Rispin(2,1))
                    ALML(l,:(2*L+1),nemin:nemax,:nrf,2) = ALM(:(2*L+1),nemin:nemax,:nrf,1)*(CFAC*Rispin(1,2)) + ALM(:(2*L+1),nemin:nemax,:nrf,2)*(CFAC*Rispin(2,2))
                    !ALML(l,:(2*L+1),nemin:nemax,:nrf,:iso) = ALM(:(2*L+1),nemin:nemax,:nrf,:iso)*CFAC
                    if (abs(projector).ge.5) then
                       ! This implements spin rotation depending on the group operation. The orbitals are rotated by rotation of momenta, while spin rotation needs to be added
                       ! to Gamma * chi_{\vK,s} = chi_{Gamma^{-1}*\vK,\Gamma^{-1}s}
                       al_interstitial(:(2*L+1),:nbands,1,lcase) = a_interstitial(:(2*L+1),:nbands,1)*(CFAC*Rispin(1,1)) + a_interstitial(:(2*L+1),:nbands,2)*(CFAC*Rispin(2,1))
                       al_interstitial(:(2*L+1),:nbands,2,lcase) = a_interstitial(:(2*L+1),:nbands,1)*(CFAC*Rispin(1,2)) + a_interstitial(:(2*L+1),:nbands,2)*(CFAC*Rispin(2,2))
                       !al_interstitial(:(2*L+1),:nbands,:iso,lcase) = a_interstitial(:(2*L+1),:nbands,:iso)*CFAC
                    endif
                 else
                    ALML(l,:(2*L+1),nemin:nemax,:nrf,1) = ALM(:(2*L+1),nemin:nemax,:nrf,1)*CFAC
                    if (abs(projector).ge.5) al_interstitial(:(2*L+1),:nbands,1,lcase) = a_interstitial(:(2*L+1),:nbands,1)*CFAC
                 endif
              enddo  !--------- end of atom L loop (alm)  ----------------!
              !----------------  ALML(l,m,iband,ifr,ispin) contains (A*C,B*C) for all l,m,iband,ifr=(1,2),ispin --------!

              call cputim(time1)
              alm_time = alm_time + time1-time0
              zgemm_time = zgemm_time + zgem2_time
              time0=time1
              !---------- Computing the DMFT transformation --------------!
              do lcase=1,nl(icase)
                 l1=ll(icase,lcase)
                 nind=(2*l1+1)*iso
                 idt = idet(isym)
                 icix = cix(icase,lcase)
                 iorb = iorbital(icase,lcase)
                 CALL CmpDMFTrans(DMFTrans, DMFTU, alml, ri_mat, projector, cfX, l1, iorb, cix_orb, nindo, iso, nemin, nbands, nind, maxdim2, icix, idt, icase, lcase, lmax2, nume, nrf, nrf, norbitals)
                 if (mode.EQ.'x' .and. cmp_partial_dos) then
                    CALL CmpXqtl(xqtl2(:,:,:,iorb), DMFTrans(:,:,:,:,iorb), nbands, nind, maxdim2)
                 endif
              enddo
              call cputim(time1)
              trans_time = trans_time + time1-time0
           ENDDO  ! over atom-cases
           !-------- computing correlated green's function --------!

           call cputim(time0)
           if (Qrenormalize) then
              if (abs(projector).le.5) then
                 CALL RenormalizeTrans(DMFTU, Olapm0, SOlapm, cix_orb, cixdim, nindo, iSx, projector, nbands, maxdim2, norbitals, maxdim, ncix, RENORM_SIMPLE)
              else
                 CALL RenormalizeTransK(DMFTU, cix_orb, cixdim, nindo, iSx, Sigind, projector, nbands, maxdim2, norbitals, maxdim, ncix)
              endif
           endif

           !call Debug_Print_Projector(DMFTU,nindo,norbitals,nbands,maxdim2)

           CALL CompressSigmaTransformation2(STrans, DMFTU, Sigind, iSx, cix, iorbital, ll, nl, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize) ! WWW
           call cputim(time1)
           comprs_time = comprs_time + time1-time0
           time0=time1

           IF (mode.EQ.'g') THEN
              if (cmp_partial_dos) CALL CompressGcTransformation4(GTrans, DMFTrans, nbands, nindo, norbitals, maxdim2)
              CALL GetImpurityLevels2(Eimpmk, Olapmk, DMFTU, STrans, E, EF, s_oo, nemin, nume, iSx, iorbital, ll, nl, csize, cix, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize)
              Olapm(:,:,:) = Olapm(:,:,:) + Olapmk(:,:,:) * (mweight(ikp)/tweight/nsymop2) ! proper weight for the reducible k-point
              Eimpm(:,:,:) = Eimpm(:,:,:) + Eimpmk(:,:,:) * (mweight(ikp)/tweight/nsymop2) ! proper weight for the reducible k-point
           ELSE IF (mode.EQ.'u') THEN
              if (pform) then
                 WRITE(fh_p,*) isym
              else
                 WRITE(fh_p) isym
              endif
              DO iorb=1,norbitals
                 DO ind=1,nindo(iorb)
                    if (pform) then
                       do i=1,nbands
                          write(fh_p,'(f16.9,1x,f16.9,3x)',advance='no') dble(DMFTU(i,ind,iorb)), aimag(DMFTU(i,ind,iorb))
                       enddo
                       write(fh_p,*)
                    else
                       write(fh_p) (DMFTU(i,ind,iorb),i=1,nbands)
                    endif
                 ENDDO
              ENDDO
           ENDIF
           call cputim(time1)
           imp_time = imp_time + time1-time0

           if (Cohfacts) WRITE(1002,'(A,5I5,2x,A)') '#', ikp, isym, nbands, nemin, nomega, ': ikp, isym, nbands nemin, nomega, kweight'
           if (vector_out) ALLOCATE( A2(nmat,nbands) )

           !$OMP PARALLEL DO SHARED(gtot,gmloc,gloc,Ekp,sigma,STrans,cohf) PRIVATE(gij,i,zek,evl,evr,xomega,gk,gmk,gtc) SCHEDULE(STATIC)
           do iom=1,nomega                !-------------------- over frequency ---------------------------------!
              ALLOCATE( gij(nbands,nbands) )
              if (mode.EQ.'e') allocate( zek(nbands) )
              if (mode.EQ.'g') allocate( gmk(maxdim,maxdim,ncix) )
              if (mode.EQ.'g' .and. cmp_partial_dos) allocate( gk(norbitals) )

              if (Cohfacts) allocate( evl(nbands,nbands), evr(nbands,nbands))
              gij=0                       !-------- band structure part of g^-1 --------------------------------!           
              ! Building of G^-1 or H_LDA
              IF (mode.EQ.'e') THEN         !--- We will only compute eigenvalues and eigenvectors
                 DO i=1,nbands
                    gij(i,i) = E(i+nemin-1)  !---   preparing hamiltonian: epsk+sigma
                 ENDDO
                 ! Adding self-energy
                 if (ncix.gt.0) then         !-------- adding self-energy if correlated orbitals exist ------------!
                    CALL AddSigma_optimized2(gij, sigma(:,:,iom), STrans, csize, 1, nbands, ncix, maxsize) ! WWW
                 endif
                 if (Hrm(1:1).eq.'H') then ! For Hermitian mode force Hamiltonian to be Hermitian
                    gij(:,:) = 0.5*(gij(:,:)+conjg(transpose(gij(:,:))))
                 endif
                 if (Cohfacts) then
                    if (Hrm(1:1).eq.'H') then
                       CALL zheigsys(gij, zek, evl, evr, nbands)
                    else
                       CALL eigsys(gij, zek, evl, evr, nbands)
                    endif
                    call zgemm('N','N', nbands,  nbands,  nbands, (1.d0,0.d0), evl, nbands, Cohf, nbands, (0.d0,0.d0), gij, nbands)
                    call zgemm('N','N', nbands,  nbands,  nbands, (1.d0,0.d0), gij, nbands, evr, nbands, (0.d0,0.d0), evl, nbands)
                    !$OMP CRITICAL 
                    WRITE(1002,'(F14.8,1x)',advance='no') omega(iom)*Ry2eV
                    do i=1,nbands
                       WRITE(1002,'(f14.8,1x,f14.8,2x)',advance='no') evl(i,i)
                    enddo
                    WRITE(1002,*)
                    !$OMP END CRITICAL 
                 else if ( vector_out .and. iom.eq.iom_zero) then
                    allocate( evr(nbands,nbands), dek(nbands) )
                    ! in HV mode we write out new Hermitian vector file, which can be used to determine irreducible representations of every band
                    CALL dheigsys(gij, dek, evr, nbands)
                    !$OMP CRITICAL 
                    write(6,6000) S,T,Z,KNAME,nbands, 1.0,' ', (dek(i),i=1,nbands)
                    write(6,6010) nemin, E(nemin)
                    write(6,6030)
                    !$OMP END CRITICAL
                    do is=1,iso
                       call zgemm('N','N', N,  nbands,  nbands, (1.d0,0.d0), A(1,nemin,is), nmat, evr, nbands, (0.d0,0.d0), A2, nmat)
                       !$OMP CRITICAL 
                       do num=1,nemin-1
                          write(vcn_out(is))  num,E(num)
                          if (.not.Qcomplex) then
                             write(vcn_out(is)) (dble(A(i,num,is)),i=1,N)
                          else
                             write(vcn_out(is)) (A(i,num,is),i=1,N)
                          endif
                       enddo
                       do num=nemin,nemin+nbands-1
                          iband = num-nemin+1
                          write(vcn_out(is))  num, dek(iband)
                          if (.not.Qcomplex) then
                             write(vcn_out(is)) (dble(A2(i,iband)),i=1,N)
                          else
                             write(vcn_out(is)) (A2(i,iband),i=1,N)
                          endif
                       enddo
                       do num=nemin+nbands,NE
                          write(vcn_out(is))  num,E(num)
                          if (.not.Qcomplex) then
                             write(vcn_out(is)) (dble(A(i,num,is)),i=1,N)
                          else
                             write(vcn_out(is)) (A(i,num,is),i=1,N)
                          endif
                       enddo
                       !$OMP END CRITICAL 
                    enddo
                    deallocate( evr, dek )
                 else
                    if (Hrm(1:1).eq.'H') then
                       CALL zheigvals(gij, zek, nbands)
                    else
                       CALL eigvals(gij, zek, nbands)
                    endif
                 endif
                 Ekp(iom,isym,:nbands,iikp) = zek(:nbands)
              ELSEIF (mode.EQ.'g') THEN
                 IF (matsubara) THEN
                    xomega = omega(iom)*IMAG
                 ELSE
                    xomega = omega(iom)
                 ENDIF
                 DO i=1,nbands
                    gij(i,i) = xomega+EF-E(i+nemin-1) + IMAG*gamma  !-----  g^-1 of the LDA part in band representation ----!
                 ENDDO
                 call cputim(time0)
                 ! Adding self-energy
                 if (ncix.gt.0) then         !-------- adding self-energy if correlated orbitals exist ------------!
                    CALL AddSigma_optimized2(gij, sigma(:,:,iom), STrans, csize, -1, nbands, ncix, maxsize) ! WWW
                 endif
                 call cputim(time1)
                 comprs_time = comprs_time + time1-time0
                 time0=time1
                 if (ncix.gt.0) then
                    CALL zinv(gij,nbands)    !-------- inversion of matrix to get g -------------------------------!
                 else
                    DO i=1,nbands
                       gij(i,i) = 1/gij(i,i) !-------- g is diagonal if no correlations ---------------------------!
                    ENDDO
                 endif
                 call cputim(time1)
                 inv_time = inv_time + time1-time0
                 time0=time1

                 CALL CmpGkc2(gmk, gij, DMFTU, iSx, iorbital, ll, nl, cix, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands)
                 gmloc(:,:,:,iom) = gmloc(:,:,:,iom) + gmk(:,:,:)*(mweight(ikp)/tweight/nsymop2) ! proper weight for the reducible k-point

                 ! G not resolved in m-quantum number. It is computed by using wave functions, u,\cdot{u},...
                 if (cmp_partial_dos) then
                    CALL CmpGknc(gk, gij, GTrans, nindo, norbitals, ncix, nbands)
                    gloc(:,iom) = gloc(:,iom) + gk(:)*(mweight(ikp)/tweight/nsymop2) ! proper weight for the reducible k-point
                 endif

                 ! Total g
                 gtc=0
                 DO i=1,nbands       ! over bands-1
                    gtc = gtc + gij(i,i)
                 ENDDO
                 ! The rest of the bands added to total g!
                 DO i=1,nemin-1 
                    gtc = gtc + 1/(xomega+EF-E(i) + (0.d0, 1.d0)*gamma*5)
                 ENDDO
                 ! nemax
                 DO i=nemax+1,NUM_MAX
                    gtc = gtc + 1/(xomega+EF-E(i) + (0.d0, 1.d0)*gamma*5)
                 ENDDO
                 gtot(iom) = gtot(iom) + gtc*(mweight(ikp)/tweight/nsymop2) 

                 call cputim(time1)
                 gc_time = gc_time + time1-time0

              ELSEIF (mode.EQ.'x') THEN
                 ! do nothing
              ELSEIF (mode.EQ.'u') THEN
                 ! do nothing
              ELSE
                 print *, 'mode not yet implemented 1!'
              ENDIF
              deallocate( gij )
              if (mode.EQ.'e') deallocate( zek )
              if (mode.EQ.'g') deallocate( gmk )
              if (mode.EQ.'g' .and. cmp_partial_dos) deallocate( gk )
              if (Cohfacts) DEALLOCATE( evl, evr )
           enddo  !--  iom loop
           !$OMP END PARALLEL DO

           if (mode.EQ.'g' .and. cmp_DM) then

              call cputim(time0)
              mweight_ikp = (mweight(ikp)/tweight/nsymop2)
              call cmp_DensityMatrix(g_inf, g_ferm, EF, E, s_oo, STrans, DMFTU, omega, nl, ll, iSx, cix, mweight_ikp, iso, iorbital, nbands, nume, nemin, maxdim, maxdim2, ncix, nomega, csize, maxsize, natom, lmax2, norbitals)
              call cputim(time1)
              dm_time = dm_time + time1-time0
           endif
           if (vector_out) DEALLOCATE( A2 )
        ENDDO !------  isym: over star of the irreducible k-point ----!

        if (Cohfacts) then
           DEALLOCATE(Cohf)
        endif
        deallocate( DMFTU, STrans )
        if (cmp_partial_dos) deallocate( DMFTrans, GTrans )

        if (abs(projector).ge.5) then
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

  END DO ! over different vector files

  !---- Deallocating some arrays before MPI_Gather
  CALL w_deallocate()
  if (abs(projector).ge.5) call p_deallocate()
  DEALLOCATE(a_real)
  if (mode.EQ.'g') then
     deallocate( Olapmk, Eimpmk )
  endif

  if (Cohfacts) THEN
     DEALLOCATE(Cohf0)
     CLOSE(1002)
     CLOSE(1001)
  endif

  if ( vector_out ) then
     do is=1,iso
        close(vcn_out(is))
     enddo
  endif
  
  !if (iso.eq.2) then
  deallocate( crotloc_x_rotij )
  !endif

  !DEALLOCATE( Rspin )

  !***************** MPI CALLS ***************************
  if (mode.eq.'g') then
     CALL Reduce_MPI(gloc, gtot, gmloc, Olapm, Eimpm, norbitals, nomega, maxdim, ncix)
     if (cmp_DM) then
        CALL Reduce_dm_MPI(g_inf, g_ferm, maxdim, ncix, nomega)
     endif

  elseif (mode.eq.'e') then
     if (myrank.eq.master) then
        ALLOCATE( tEk(nomega,nsymop,max_bands,numkpt) )
        ALLOCATE( tnbands(numkpt), tnemin(numkpt), tn_ik(numkpt) )
     endif

     CALL Gather_MPI(tEk, tnbands, tnemin, tn_ik, Ekp, nbandsk, nemink, n_ik, pr_procr, pr_procs, nprocs, numkpt, max_bands, nomega, nsymop) 
  endif
  !***************** MPI CALLS ***************************
  DEALLOCATE( pr_procs )

  !------------- Printing the calculated local green's function for all requested orbitals in requested shape --------!
  IF (mode.EQ.'g' .and. myrank.eq.master) THEN

     if (Qsymmetrize) then
        WRITE(6,*) 'Symmetrizing the Greens function since the group symmetrization is turned off.'
        Call SymmetrizeLocalQuantities(gmloc, Eimpm, Olapm, s_oo, sigma, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, noccur, maxdim2, nomega)
        if (cmp_DM) Call SymmetrizeLocalQuantities2(g_inf, g_ferm, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2, nomega)
     endif

     if (cmp_DM) then
        call Print_DensityMatrix(gmloc, g_inf, g_ferm, omega, maxdim, ncix, Sigind, nomega, cixdim, iso)
        deallocate( g_inf )
        deallocate( g_ferm )
     end if

     allocate( Deltac(maxsize, ncix, nomega), Glc(maxsize, ncix, nomega), Eimpc(maxsize, ncix), Olapc(maxsize,ncix) )

     CALL GetDelta2(Deltac, Glc, Eimpc, Olapc, logGloc, matsubara, omega, sigma, s_oo, gamma, gmloc, Eimpm, Olapm, Sigind, csize, cixdim, noccur, ncix, maxsize, maxdim, nomega, projector, ComputeLogGloc, Qsymmetrize)

     CALL PrintGloc(fh_dos, fh_gc, fh_dt, Glc, gloc, gtot, Deltac, omega, csize, csizes, nl, ll, legend, iatom, ncix, nomega, natom, norbitals, maxsize, Ry2eV)

     WRITE(6,'(A)',advance='no') 'Eimp='
     do icix=1,ncix
        do it=1,csize(icix)
           WRITE(6,'(F15.8)',advance='no') Eimpc(it,icix)*Ry2eV
        enddo
     enddo
     WRITE(*,*)
     WRITE(*,'(A)',advance='no') 'final Z due to E-window='
     do icix=1,ncix
        do it=1,csize(icix)
           WRITE(*,'(F15.8)',advance='no') Olapc(it,icix)
        enddo
        WRITE(*,*)
     enddo

     do icix=1,ncix
        write(fh_Eimp+icix,'(A)',advance='no') 'Ed=['
        do it=1,csizes(icix)
           write(fh_Eimp+icix,'(f14.8,A)',advance='no') Eimpc(it,icix)*Ry2eV, ','
        enddo
        write(fh_Eimp+icix,'(A)') ']'
        write(fh_Eimp+icix,'(A)',advance='no') 'Olap=['
        do it=1,csizes(icix)
           write(fh_Eimp+icix,'(f14.8,A)',advance='no') Olapc(it,icix), ','
        enddo
        write(fh_Eimp+icix,'(A)') ']'
        write(fh_Eimp+icix,'(A,f14.8)') 'logGloc=',sum(logGloc)*Ry2eV
     enddo
     deallocate( Deltac , Glc, Eimpc, Olapc)

  ELSEIF (mode.EQ.'e' .and. myrank.eq.master) THEN
     ! Print Eigenvalues

     ALLOCATE( korder(numkpt) )
     do ikp=1,numkpt
        korder(tn_ik(ikp))=ikp
     enddo

     DO j=1,numkpt
        ikp = korder(j)
        nbands = tnbands(ikp)
        nemin = tnemin(ikp)
        DO isym=1,nsymop
           WRITE(fh_eig,'(A,5I5,2x,F23.20,2x,A)') '#', j, isym, nbands, nemin, nomega, mweight(ikp)/tweight/nsymop, ': ikp, isym, nbands nemin, nomega, kweight'
           DO iom=1,nomega
              WRITE(fh_eig,'(F19.14,2x)',advance='no') omega(iom)*Ry2eV
              DO i=1,nbands
                 WRITE(fh_eig,'(2E24.16,1x)',advance='no') tEk(iom,isym,i,ikp)*Ry2eV
              ENDDO
              WRITE(fh_eig,*)
           ENDDO
        ENDDO
     ENDDO
     DEALLOCATE( korder )
     DEALLOCATE( tEk, tnbands, tnemin, tn_ik )
  ENDIF
  !------ Deallocation of memory ---------------

  DEALLOCATE( sigma, omega, s_oo)
  deallocate( noccur )
  deallocate( cixdim, iSx, nindo, cix_orb, cfX )
  deallocate( csort )

  if (mode.EQ.'g') then
     if (cmp_partial_dos) DEALLOCATE( gloc )
     DEALLOCATE( gmloc, gtot )
     deallocate( Olapm, Eimpm )
  elseif (mode.EQ.'e') then
     DEALLOCATE(Ekp, nbandsk, nemink, n_ik)
  endif

  if (Qrenormalize .and. abs(projector).le.5) deallocate( Olapm0, SOlapm )
  !------ Deallocation of memory ---------------

  if (Qprint) then
     write(6,'(//,3X,"=====>>> CPU TIME SUMMARY",/)')
     write(6,*)
     write(6,'(a,f10.2)') '   read files        :', read_time
     write(6,'(a,f10.2)') '   atpar             :', atpar_time
     write(6,'(a,f10.2)') '   cmp overlap       :', overlap_time
     write(6,'(a,f10.2)') '   bessel and Ylm    :', bess_time
     write(6,'(a,f10.2)') '   zgemm inside      :', zgemm_time
     write(6,'(a,f10.2)') '   cmp alms          :', alm_time
     write(6,'(a,f10.2)') '   cmp projector     :', trans_time
     write(6,'(a,f10.2)') '   compress projectr :', comprs_time
     write(6,'(a,f10.2)') '   impurity levels   :', imp_time
     write(6,'(a,f10.2)') '   matrix inversion  :', inv_time
     write(6,'(a,f10.2)') '   greens function   :', gc_time
     write(6,'(a,f10.2)') '   density matrix    :', dm_time
  endif



  RETURN                    
2032 FORMAT(50X,I2,//)
6000 FORMAT(/5X,'K=',3F10.5,3X,A10/5X,' MATRIX SIZE',I6,'  WEIGHT=',F5.2, '  PGR: ',A3,/5X,'EIGENVALUES ARE:',8(/2X,5F13.7))
6010 FORMAT(I13,' EIGENVALUES BELOW THE ENERGY ',F10.5)
6030 FORMAT(7X,14('****')/)
END SUBROUTINE L2MAIN


SUBROUTINE Debug_Print_Projector(DMFTU,nindo,norbitals,nbands,maxdim2)
  IMPLICIT NONE
  COMPLEX*16, intent(in) :: DMFTU(nbands,maxdim2,norbitals)
  INTEGER, intent(in)    :: nindo(norbitals), norbitals, nbands, maxdim2
  !locals
  INTEGER :: iorb, ind, iband
  open(988,FILE='U_1.dat',form='formatted',status='unknown',access='append')
  do iorb=1,norbitals
     do ind=1,nindo(iorb)
        WRITE(988,*) 'iorb=', iorb, 'ind=', ind
        do iband=1,nbands
           WRITE(988,'(2F20.14,1x)',advance='no') DMFTU(iband,ind,iorb) 
        enddo
        WRITE(988,*)
     enddo
  enddo
  close(988)
END SUBROUTINE Debug_Print_Projector
         

              

