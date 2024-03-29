SUBROUTINE L2MAIN(coord,NSPIN1,sumw,tclm,tclm_w,tfour,tfour_w)
  USE defs,  ONLY: PI, ZERO
  USE param, ONLY: IBLOCK, LMAX2, LOMAX, NCOM, NLOAT, NRAD, NMAT, NUME, NSYM, NKPT, fh_vec, numkpt, filename_V_nsh, filename_V_vns, filename_V_sph, nwave, iff1, iff2, iff3, kmax!, NGAU
  USE atspdt,ONLY: e_store!, EL, P, DP, PE, DPE, PEI
  USE com,   ONLY: NAT, EMIN, EMAX, ELECN, XWT!, rel, NSPIN, NK, EF, NB, weigh, MAXWAV, MINWAV,  NBAND
  USE charp, ONLY: init_charp, zero_charp, fini_charp
  USE chard, ONLY: init_chard, zero_chard, fini_chard
  USE charf, ONLY: init_charf, zero_charf, fini_charf
  USE lo,    ONLY: elo_store!, nlov, nlon, nlo, loor, lapw, ilo, alo, blo, clo, pi12lo, pe12lo, pr12lo!, a1lo, b1lo
  USE lohelp,ONLY: init_lohelp, fini_lohelp !, u21, ue21, u12, ue12, u22, zerotc_lohelp, !, sum12, sum21, sum22, sume21, sume12
  USE structure, ONLY: iatnr, mult, isplit, jri, pos, rotij, rmt, rotloc, r0, dx, vol, br1, natm, cform, ZZ !, v, tauij, br2
  USE xa,    ONLY: LM, R, fj, dfj, E, WEIGHT, init_xa, fini_xa !, TC100, TCA100, TCB100, SUMA, SUMB, SUMAB, SUMBA, AVEC, BK, BKRLOC, BKROT
  USE xa3,   ONLY: bk3, aK, bk3lo, init_xa3, fini_xa3, k3, k3lo, As, As_lo!, A, A_lo
  USE reclat,ONLY: inst, tauk, kzz
  USE muzero,ONLY: nomq, jomq, iomq, womq, abom, n0_om, init_muzero, fini_muzero
  USE com_mpi,ONLY: nprocs,myrank,master,FilenameMPI,Reduce_MPI,Gather_procs,FindMax_MPI,pr_proc,pr_procr,pr_procs,vector_para,vectors,nvector,fvectors,VECFN,FindMaxK_MPI,Qprint, stop_MPI
  USE p_project,ONLY: p_allocate, p_deallocate, p_cmp_rixmat, phi_jl, rix, P_rfi, l_al_ucase, kNri, n_al_ucase, dri, j_al_ucase
  USE w_atpar, ONLY: w_allocate0,  w_allocate, w_deallocate1, w_deallocate2, w_lm, w_lmmax, w_RHOLM,w_xwt1,w_xwt1l,w_xwt1h,w_xwteh,w_xwtel, precmp_w_atpar, Read_nsh, Deallocate_nsh
  USE Forces, ONLY: Force4, Force4_mine, Force_surface, forcea, Qforce, init_forces, fini_forces, Force_surface_extra
  USE readPotential, ONLY: CheckFilesExist, close_V_vns, close_V_vsp, init_V_vns, init_V_vsp, read_V_vns, read_V_vsp
  USE char, ONLY: modus
  USE dmfts  !, ONLY: recomputeEF, wl
  USE sym2,   ONLY: tau, iz, iord
  USE xa2, ONLY: lda_weight, lda_nbmax
  USE fftw3_omp, ONLY: fft_init, fft_fini, fft_run, fft, fft_init_step1, fft_init_step2, fft_fini_step2
  USE fact, only: fct
  USE Fermi, only: cmp_EF
  IMPLICIT NONE
  CHARACTER*5, intent(in)  :: coord
  REAL*8, intent(in)       :: sumw
  INTEGER, intent(in)      :: NSPIN1
  REAL*8, intent(out)      :: TCLM,TCLM_w,TFOUR,TFOUR_w
  ! Common blocks
  !REAL*8     :: FCT
  !COMMON /FACT/   FCT(100)
  ! Functions
  Interface
     integer FUNCTION CountSelfenergy(fh_sig, ncorr)
       integer, intent(in) :: fh_sig, ncorr
     end FUNCTION CountSelfenergy
     real*8 FUNCTION FreeE0(Energy, Temp)
       real*8, intent(in) :: Energy, Temp
     end FUNCTION FreeE0
  end interface
  interface
     REAL*8 Function Romb(y,N,dx)
       IMPLICIT NONE
       REAL*8, intent(in) :: y(N)
       REAL*8, intent(in) :: dx
       INTEGER, intent(in):: N
     end Function Romb
  end interface
  ! interfaces
  ! Locals
  real*8, PARAMETER       :: Ry2eV = 13.60569253d0
  !
  !CHARACTER*10  :: KNAME
  CHARACTER*100 :: CDUMMY, filename
  CHARACTER*200 :: FNAME
  LOGICAL    :: Rho_Renormalize, Tcompute
  LOGICAL    :: debug, more_kpoints
  INTEGER    :: I, J, N, ISCF, LFIRST, jatom, ilm1, ilm3, n0, jlm, lmmax, imax, nemin, nemax, DM_nemin, DM_nemax, DM_nemaxx, label, i1
  INTEGER    :: nnlo, isize, num, ilm2, ilm, li, mi, it, ip, iq, ndim
  REAL*8     :: emist, fac, ETOT2, FTOT2, FTOT0, SQRT2, SQFP, Y, TDE1, TDE2!, S,T,Z
  REAL*8     :: time_bl, time_bl_w, time_reduc, time_reduc_w, time_write, time_writeclm, time_writescf, time4_w, time3, time4, time_force, time_force_w
  REAL*8     :: time_write_w, time_writeclm_w, time_ilm, time_ilm_w, time_radprod, time_radprod_w, time_m, time_m_w, time_rad, time_rad_w, time3_w
  REAL*8     :: time_rd_w, time_rd_c, time_atpar_w, time_atpar_c, time_writescf_w, time_r_w, time_r_c, time_lo, time_lo_w, time_lnG,time_lnGw
  REAL*8     :: time_writeefg, time_writeefg_w, t1c, t1w, t2c, t2w, t3c, t3w, time1, time1_w, time2, time2_w, time_mt, time_mt_w
  REAL*8     :: time_dmf, time_dmfw, time_int, time_intw, time_ed, time_edw, time_dmf0, time_dmf0w, time_dmf1, time_dmf1w, time_dmfwgh, time_dmfwgh_w
  !
  INTEGER    :: LM_MAX, ikp, is, itape, icix, iom, max_nbands, norbitals, maxdim2, nbands, nbandsx, nbands_dft, nomega, npomega, fh_sig, fh_dos, nip, nipc
  INTEGER    :: ift1, ift2, ift3, iff1t, iff2t, iff3t, ia1, ia2, jx, qmax, iikp, iks, itape1, ivector, iind, ir, l, Nri, icase, nedim, info
  REAL*8     :: OVRDIV, TC, EFGFACT, EFG20, EFG22, EFG2M, QXX, QYY, QZZ, QXY, VXX, VYY, VZZ, vnorm1, twgh, wgamma, rx
  REAL*8     :: volin, beta, logG, logGloc, dens, eimp_nd, eimp_nd2, DeltaG, dEtot, EF, TrGSigVdc
  LOGICAL    :: RENORM_SIMPLE, Qforce_j
  complex*16 :: wkp
  COMPLEX*16, ALLOCATABLE :: cfX(:,:,:,:), DMFTU(:,:,:,:), STrans(:,:,:,:), Aweight(:,:), Olapm0(:,:,:), SOlapm(:,:,:), AEweight(:,:)
  COMPLEX*16, ALLOCATABLE :: sumfft(:,:,:), tloc(:,:,:), rho1(:), rhok(:), Rspin(:,:,:)!, vrho1(:), vrhok(:)
  COMPLEX*16, ALLOCATABLE :: ekin1(:), ekink(:), Edimp0(:,:,:)
  REAL*8,     ALLOCATABLE :: zw2(:), zwe(:), LowE(:,:,:), omega(:), Edimp(:)
  INTEGER,    ALLOCATABLE :: nindo(:), cix_orb(:), cixdim(:), iSx(:,:), noccur(:,:), iorbital(:,:)
  COMPLEX*16, ALLOCATABLE :: sigma(:,:,:), s_oo(:,:), gloc(:), Gdloc(:,:,:), lgTrans(:,:,:,:), wEpsw(:,:), Gdloc0(:,:,:,:)
  REAL*8,     ALLOCATABLE :: aKR(:), jlr(:), jlp(:), Vlm(:,:), Vr(:), rh1_tmp(:,:), rh2_tmp(:,:), rh3_tmp(:,:)
  REAL*8,     ALLOCATABLE :: Nds(:,:), dNds(:,:), Nds_new(:)
  INTEGER,    allocatable :: lg_deg(:) 
  !INTEGER    :: keigen(3,nmat)
  REAL*8     :: zw0(nume), Kn(3)
  COMPLEX*16 :: DM(maxdim,maxdim,ncix), dDM(maxdim,maxdim,ncix) ! Density matrix
  REAL*8     :: fsph(3,natm)
  REAL*8     :: fsph2(3,natm)
  REAL*8     :: fnsp(3,natm)
  REAL*8     :: fvdrho(3,nat), fdvrho(3,nat)
  REAL*8     :: forb(3,natm), fextra(3,nat), fsur(3,nat), fkdc(3), fsur2(9), fsur_norm
  LOGICAL    :: Qforce2

  Qforce2=.True.
  
  RENORM_SIMPLE = .False.
  Rho_Renormalize = .FALSE.
  debug = .FALSE.
  wgamma = 1e-10
  time_dmf1=0
  time_dmf1w=0
  fh_sig = 80
  fh_dos = 500

  ETOT2=0.0
  FTOT2=0.0
  FTOT0=0.0
  XWT=0.0                                                           
  SQRT2=SQRT(2.0D0)                                                 
  SQFP=SQRT(4.D0*PI)                                                

  CALL Get_A_Dimensions(norbitals, maxdim2, lmaxp, iso, natom, nl, ll)
  
  CALL init_lohelp
  CALL init_xa(nat)

  allocate( iorbital(natom,lmaxp+1) )
  allocate( nindo(norbitals), cix_orb(norbitals), cixdim(ncix) )
  allocate( iSx(maxdim2, norbitals) )
  allocate( cfX(maxdim2,maxdim2,norbitals,norbitals) )
  allocate( noccur(maxsize,ncix) )
  CALL Set_A_Arrays(iorbital, nindo, cix_orb, iSx, noccur, cixdim, cfX, CF, norbitals, iso, natom, nl, ll, cix, ncix, maxdim, maxdim2, lmaxp, Sigind, csize, maxsize)
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
  do icix=1,ncix
     CALL ReadSelfenergy(fh_sig+icix, sigma(:,icix,:), omega, s_oo(:,icix), gammac, csize(icix), nomega, maxsize)
  enddo
  
  !--------------  change from eV to Rydbergs --------------------------!
  gamma = gamma/Ry2eV
  !---------------- Reading self-energy --------------------------------!
  
  ALLOCATE( LowE(5,maxsize,ncix) )
  LowE = 0.0
  do icix=1,ncix
     IF (recomputeEF.gt.1) THEN
        CALL ApproximateSelfenergy(LowE(:,:,icix), sigma(:,icix,:), omega, csize(icix), nomega, maxsize, matsubara, s_oo(:,icix), WL)
     ENDIF
  enddo
  

  CALL CountFrequencyArrays2(n0_om, qmax, beta, omega, nomega, matsubara)
  
  if (matsubara) Temperature = 1/beta

  CALL init_muzero(nomega,qmax)

  CALL CreateFrequencyArrays2(nomq, jomq, iomq, womq, abom, npomega, n0_om, qmax, beta, omega, nomega, matsubara, Temperature)

  ALLOCATE( Rspin(2,2,norbitals) )
  Rspin=0.0

  if (iso.EQ.2) CALL GetSpinRotation(Rspin,rotij,crotloc,norbitals,natom,natm,iso,lmaxp,iatom,nl,cix,ll,iorbital)
  
  
  vnorm1=1.d0  
  IF((iso.EQ.2).AND.(nspin1.EQ.1)) vnorm1=0.5d0

  Y=1.0D0                                                           
  DO I=1,49
     J=2*I-1
     FCT(J)=Y  ! FCT(i) = ((i+1)/2)!
     Y=Y*I
  END DO
  
  LFIRST=1
  nnlo=0
  time_bl=zero; time_bl_w=zero; time_reduc=zero; time_reduc_w=zero
  time_write=zero; time_writeclm=zero; time_writescf=zero; time_write_w=zero
  time_writeclm_w=zero; time_writescf_w=zero; time_ilm=zero; time_ilm_w=zero
  time_radprod=zero; time_radprod_w=zero; time_m=zero; time_m_w=zero; time_rad=zero; time_rad_w=zero; time_rd_w=zero; time_rd_c=zero
  time_r_w=zero; time_r_c=zero; time_atpar_w=zero; time_atpar_c=zero; time_writeefg=zero; time_writeefg_w=zero
  time_dmf=zero; time_dmfw=zero; time_int=zero; time_intw=zero; time_ed=zero; time_edw=zero
  time_dmfwgh=zero; time_dmfwgh_w=zero; time_mt=zero; time_mt_w=0; time_lo=0; time_lo_w=0;
  time_force=zero; time_force_w=zero
  
  allocate (e_store(0:lmax2,nat),elo_store(0:lomax,1:nloat,nat))
  
  itape=9
  if (vector_para) then
     if (nvector.ge.1) then
        FNAME = fvectors(1,1)
     else   ! No k-point to compute, but still need linearization energies
        FNAME = TRIM(VECFN(1))//'_1'
     endif
     open(itape,FILE=FNAME,STATUS='old',FORM='unformatted')
  else
     !------------ vector files need to be read from the beginning -------------------!           
     rewind(itape) 
  endif
  DO i=1,nat
     READ(itape) e_store(0:lmax2,i)
     READ(itape) elo_store(0:lomax,1:nloat,i)
  enddo
  if (vector_para) then
     close(itape)
  else
     rewind(itape) 
  endif

  !print *, 'CheckFileExist=', CheckFilesExist(filename_V_nsh, filename_V_vns)
  !print *, 'modus=', modus, 'True=', (modus.eq.'FOR ')
  Qforce = CheckFilesExist(filename_V_nsh, filename_V_vns) .AND. (modus.eq.'FOR ')
  if (myrank.EQ.master) then
     write(6,*) '    FORCE-CALCULATION:', Qforce
  end if
  
  CALL w_allocate0(nat)

  call init_forces(nat)
  if (Qforce) then
     fsph(:,:) = 0.d0
     fnsp(:,:) = 0.d0
     fsph2(:,:) = 0.d0
     fvdrho(:,:) = 0.d0
  endif
  allocate( lg_deg(ntcix) )
  lg_deg(:) = 1    ! Here was a bug in 2017. When averaging over different cores, some might have no k-points and lg_deg might be zero, causing NaN's after averaging.

  ! Reads case.in2 to determin LM array
  ! and maximum LM -> LM_MAX=max(LM(1,:))
  LM_MAX=0
  do jatom=1,nat
     ! neg L MEANS NEGATIVE SPHERICAL HARMONIC COMB (SEE KURKI-SUONIO)
     ! (+-l,m) == (l,m,+-) ==> +- is (cos(phi),sin(phi))
     read(5,1003) ( (lm(j,jlm),j=1,2), jlm=1,ncom )
     ! Find maximum LMMAX
     DO JLM=2,NCOM
        if (Qforce) then
           ! 1 x, -1 x  ==>  forcea(0)=.true.   either x or y or z
           ! 1 1        ==>  forcea(1)=.true.   only x
           !-1 1        ==>  forcea(2)=.true.   only y
           ! 1 0        ==>  forcea(3)=.true.   only z
           if(abs(lm(1,jlm)).eq.1)                    forcea(0,jatom)=.true.
           if((lm(1,jlm).eq.1).and.(lm(2,jlm).eq.1))  forcea(1,jatom)=.true.
           if((lm(1,jlm).eq.-1).and.(lm(2,jlm).eq.1)) forcea(2,jatom)=.true.
           if((lm(1,jlm).eq.1).and.(lm(2,jlm).eq.0))  forcea(3,jatom)=.true.
        endif
        IF(LM(1,JLM).EQ.0) EXIT
     ENDDO
     LMMAX=JLM-1
     if (LMMAX.GT.LM_MAX) LM_MAX = LMMAX
     
     !if (Qprint) then
     !   write(6,*) ' atom',jatom,' ncomu',0,' lmmax',LMMAX
     !   WRITE(6, 1005) LMMAX,((LM(J,JLM),J=1,2), JLM=1,LMMAX)             
     !endif
     !if (myrank.EQ.master) then
     !   WRITE(21,13)   JATOM,IATNR(JATOM),(POS(I,lfirst),I=1,3),MULT(JATOM)
     !   WRITE(21,1005) LMMAX,((LM(J,JLM),J=1,2), JLM=1,lmmax)             
     !endif
     w_lm(1:2,1:NCOM,jatom) = lm(1:2,1:NCOM)
     w_lmmax(jatom) = lmmax
  enddo

!!! BRISI
  !forcea(:,:) = .true.
!!! BRISI

  ! allocates memory to store radial wave functions and 
  ! all results of subroutine "atpar", which is called outside
  ! the k - loop.
  CALL w_allocate(LM_MAX,nat)

  CALL precmp_w_atpar(cform,zz,nnlo,ISCF)

  if (Qforce) CALL Read_nsh()
  nip=1
  if (Qforce) nip=4

  nipc=1 ! nip for correlated atoms
  do icase=1,natom ! over all correlated atoms
     Qforce_j = Qforce .AND. forcea(0,isort(iatom(icase)))
     if (Qforce_j) nipc=4
!!! BRISI
     !nipc=4
!!! BRISI
  enddo
  if (myrank.EQ.master) WRITE(6,*) 'nipc=', nipc, 'For force calculation on correlated atom it should be 4'
  
  if (abs(projector).ge.5) then
     filename="projectorw.dat"
     call p_allocate(filename)
     call p_cmp_rixmat()
  endif
  
  !---------------------------------
  ! START LOOP OVER ALL k-POINTS
  !---------------------------------
  CALL init_xa3(nmat,nnlo,nume,iso)

  
  if (vector_para) then
     pr_proc = sum(vectors(:,2))
     pr_procr = pr_proc
     if (Qprint) WRITE(6,'(A,I3,2x,A,I3)') 'pr_proc=', pr_proc, 'tot-k=', numkpt
     !WRITE(6,'(A,I3,2x,A,I3)') 'pr_procr=', pr_procr, 'tot-k=', numkpt
  else
     pr_proc  = floor(nkpt/DBLE(nprocs)+0.999)  ! The maximum number of points calculated per processor     
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

  ALLOCATE( pr_procs(nprocs) )
  CALL Gather_procs(pr_procr, pr_procs, nprocs)
  if (myrank.EQ.master .and. SUM(pr_procs).NE.numkpt) then
     WRITE(6,*) 'ERROR: sum(pr_procs) should be numkpt, but is not', sum(pr_procs), numkpt
  endif

  max_nbands=0
  ! If DMFT-transformation need to be normalized, we compute the renormalization coefficients
  if (Qrenormalize .and. abs(projector).le.5) then
     allocate(Olapm0(maxdim,maxdim,ncix), SOlapm(maxdim,maxdim,ncix) )
     call read_overlap_from_file(info, Olapm0, SOlapm, cixdim, maxdim, ncix)
     time_dmf0=0; time_dmf0w=0
     if (info.ne.0) CALL cmp_overlap(Olapm0, SOlapm, max_nbands, nnlo, norbitals, natom, maxdim2, emin, emax, iorbital, cix_orb, iSx, nindo, noccur, cixdim, cfX, Rspin, sumw, lmaxp, time_dmf0, time_dmf0w, RENORM_SIMPLE)
  endif
  
  kmax(:)=0
  ! Chemical potential is recomputed on DMFT-self-energy
  if (recomputeEF.GT.0 .or. mode.EQ.'e') then
     CALL cmp_EF(DM_EF, recomputeEF, max_nbands, Olapm0, SOlapm, omega, sigma, nnlo, norbitals, natom, maxdim2, nomega, npomega, emin, emax, iorbital, cix_orb, cixdim, iSx, nindo, cfX, Rspin, sumw, wgamma, beta, vnorm1, LowE, lmaxp, RENORM_SIMPLE, kmax, time_dmf1, time_dmf1w)
  else if (recomputeEF.LT.0) then  ! For recomputeEF=-1 we use DFT EF and DFT weights (just for debugging purposes)
     allocate( lda_weight(nkpt*1,nume) )
     CALL LDA_fermi(EF,lda_weight,lda_nbmax,nkpt,nume,1,nat)
     DM_EF=EF
  endif
  
  IF (mode.EQ.'e') goto 6170 !!! mode='e' : Computes only eigenvalues and prints them to the disc
  
  ! The new chemical potential
  IF (Qprint) WRITE(6,1060)  ELECN,DM_EF
  if (myrank.EQ.master) WRITE(21,1060) ELECN,DM_EF    

  IF (mode.EQ.'m') goto 6170 !!! mode==m : Compute only the chemical potential!



  ! ---- for interstitial ----
  ! ..... FIND MAX RECIPROCAL LATTICE VECTORS
  if (sum(kmax).eq.0) then
     kmax(:)=0
     DO ivector=1,nvector
        if (vector_para) then
           open(fh_vec,FILE=fvectors(ivector,1),STATUS='old',FORM='unformatted')
        else
           rewind(fh_vec) !--- both vector files: 9 and 10 rewind -------------------------!
        endif
        DO I=1,NAT
           READ(fh_vec) EMIST !---- At the beginninge we have linearization energies --!
           READ(fh_vec) EMIST !---- We just skip them ---------------------------------!
        ENDDO
        DO iks=1,vectors(ivector,2)   ! kpoint loop begin
           CALL Read_Vec_Spin_DontStore(nemin, DM_nemin, DM_nemax, more_kpoints, n0, emin, emax, nnlo, kmax)
           IF (.NOT.more_kpoints) EXIT
        ENDDO
        if (vector_para) then
           close(fh_vec)
        else
           rewind(fh_vec)
        endif
     END DO
     CALL FindMaxK_MPI(kmax)
  endif

  IF (Qprint) WRITE(6,*) 'max k indices:',kmax(1),kmax(2),kmax(3)
  kmax(1)=kmax(1)+1
  kmax(2)=kmax(2)+1
  kmax(3)=kmax(3)+1
  ! -----  set iff values for FFT calculation
  iff1=2*kmax(1)
  iff2=2*kmax(2)
  iff3=2*kmax(3)
  iff1=(iff1+1)*2
  iff2=(iff2+1)*2
  iff3=(iff3+1)*2
  iff1t=iff1
  iff2t=iff2
  iff3t=iff3
  CALL IFFLIM(iff1t,iff1)
  CALL IFFLIM(iff2t,iff2)
  CALL IFFLIM(iff3t,iff3)
  IF (Qprint) WRITE(6,*) 'n,iff1,iff2,iff3',n-nnlo,iff1,iff2,iff3
  !
  if (Qforce) then
     ift1 = iff1
     ift2 = iff2
     ift3 = iff3
  else
     ift1 = 1
     ift2 = 1
     ift3 = 1
  endif
  allocate( tloc(ift1,ift2,ift3) )
  allocate ( sumfft(iff1,iff2,iff3) )
  sumfft=0.0
  tloc=0.0
  ! --- finished for interstitial
  
  do jatom=1,nat
     lfirst = 1+sum(mult(1:jatom-1))
     if (Qprint) then
        write(6,*) ' atom',jatom,' ncomu',0,' lmmax',LMMAX
        write(6, 1005) LMMAX,((w_lm(j,jlm,jatom),j=1,2), jlm=1,lmmax)
     endif
     if (myrank.EQ.master) then
        WRITE(21,13)   jatom, iatnr(jatom), (POS(i,lfirst),i=1,3),mult(jatom)
        WRITE(21,1005) LMMAX,((w_lm(j,jlm,jatom),j=1,2), jlm=1,lmmax)
     endif
  end do

  if (myrank.eq.master) then
     WRITE(8,787) ISCF                        
     WRITE(8,*) '   NORM OF CLM(R) =          '                        
     WRITE(8,*)                                                 
     WRITE(6,800)
  endif
  
  w_RHOLM(:,:,:)=0.0
  w_xwt1(:,:)=0.0
  w_xwt1l(:,:)=0.0
  w_xwt1h(:,:)=0.0
  w_xwteh(:,:)=0.0
  w_xwtel(:,:)=0.0
     
  REWIND 5
  READ(5,fmt='(//,1A)') CDUMMY
  
  CALL init_charp(nume)
  CALL init_chard(nume)
  CALL init_charf(nume)

  allocate( gloc(nomega) )
  
  allocate( Gdloc(ntcix,nomega,nipc) )
  allocate( Gdloc0(maxdim,maxdim,ncix,nomega), Edimp0(maxdim,maxdim,ncix) )
  Gdloc=0
  Edimp0=0
  Gdloc0=0
!!! ---------- Preparation of arrays for paralel executaion --------------                                                     
  !pr_proc  = floor(nkpt/DBLE(nprocs)+0.999)  ! The maximum number of points calculated per processor                                          
  IF (Qprint) WRITE(6,'(A,I3,2x,A,I3)') 'pr_proc=', pr_proc, 'tot-k=', nkpt
  IF (Qprint) print *, 'mpr_proc=', pr_proc, 'myrank=', myrank, 'nprocs=', nprocs, 'nkpt=', nkpt

  allocate( Nds(nipc,ntcix), dNds(nipc,ntcix) )
  Nds=0.d0
  DM=0

  gloc=0
  dens=0
  iikp=0
  DO ivector=1,nvector
     DO is=1,iso    !------ over up/dn ---------------------!
        itape1=8+is
        !itape2=29+is
        if (vector_para) then
           open(itape1,FILE=fvectors(ivector,is),STATUS='old',FORM='unformatted')
           !open(itape2,FILE=fvectors(ivector,is+2),STATUS='old',FORM='formatted')
        else
           rewind(itape1) !--- both vector files: 9 and 10 rewind -------------------------!
           !rewind(itape2) !--- both vector files: 9 and 10 rewind -------------------------!
        endif
        DO I=1,NAT
           READ(itape1) EMIST !---- At the beginninge we have linearization energies --!
           READ(itape1) EMIST !---- We just skip them ---------------------------------!
        ENDDO
     END DO

     DO iks=1,vectors(ivector,2) ! kpoint loop begin
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
        
        call cputim(t1c)
        call walltim(t1w)
        
        !!! Need to read both spins!
        CALL Read_Vec_Spin(E, twgh, As, As_lo, k3, k3lo, bk3, bk3lo, nemin, nemax, DM_nemin, DM_nemax, more_kpoints, n0, emin, emax, DM_Emin, DM_Emax, iso, nmat, nume, nnlo)

        IF (.NOT.Tcompute) CYCLE
        IF (.NOT.more_kpoints) EXIT
        
        call cputim(t2c)
        call walltim(t2w)
        time_rd_c=time_rd_c+t2c-t1c
        time_rd_w=time_rd_w+t2w-t1w

        isize=n0-nnlo
        DO I=1,isize                                !----- Over reciprocal vectors -----------------------!
           Kn(:) = matmul(BR1,BK3(:,I))             !----  Creates K+k in cartesian coordinates ----------!
           aK(I) = sqrt(Kn(1)**2+Kn(2)**2+Kn(3)**2) !----  calculates |K+k| for use in bessel functions --!
        ENDDO
        do jatom=1,nat
           CALL HARMON(isize,aK(:isize),lmax2,fj(:,:isize,jatom),dfj(:,:isize,jatom),rmt(jatom))
        enddo
!!! BRISI
        !call Hamiltonian(n0, nnlo,ikp)
        !!if (iks.eq.5)
        !STOP
!!! BRISI
        ! New
        if (abs(projector).ge.5) then
           allocate( phi_jl(nmat, n_al_ucase) )
           
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
     
        call cputim(t1c)
        call walltim(t1w)
     
        nbands = DM_nemax-DM_nemin+1
        allocate( DMFTU(nbands,maxdim2,norbitals,nipc) )

        ! DMFT-projection is build

        CALL Build_DMFT_Projector(DMFTU, cfX, Rspin, iorbital, norbitals, nipc, n0, nnlo, nbands, cix_orb, nindo, DM_nemin, DM_nemax, maxdim2, lmaxp)
        
        if (abs(projector).ge.5) deallocate( phi_jl )

        if (Qrenormalize) then
           if (abs(projector).le.5) then
              CALL RenormalizeTrans(DMFTU, Olapm0, SOlapm, cix_orb, cixdim, nindo, iSx, nbands, nipc, maxdim2, norbitals, maxdim, ncix, RENORM_SIMPLE)
           else
              CALL RenormalizeTransK(DMFTU, cix_orb, cixdim, nindo, iSx, Sigind, projector, nbands, nipc, maxdim2, norbitals, maxdim, ncix)
           endif
        endif
        allocate( STrans(maxsize,ncix,nbands,nbands) )
        allocate( lgTrans(nbands,nbands,ntcix,nipc) )
        
        CALL CompressSigmaTransformation2(STrans, lgTrans, lg_deg, DMFTU, Sigind, Sigind_orig, iSx, cix, iorbital, ll, nl, natom, iso, ncix, maxdim, maxdim2, lmaxp, norbitals, nbands, maxsize, ntcix, nipc)

        call cputim(t2c)
        call walltim(t2w)

        time_dmf=time_dmf+t2c-t1c
        time_dmfw=time_dmfw+t2w-t1w
        
        wkp = (twgh/sumw)*2.0
        zw0(:)=0
        DO num=nemin,DM_nemin-1
           zw0(num) = dble(wkp) ! These bands are fully filled
        ENDDO

        call cputim(t1c)
        call walltim(t1w)

        nedim=1
        if (Qforce) nedim=nbands
        allocate( Aweight(nbands,nbands) )
        allocate( AEweight(nedim,nedim) )

        CALL cmp_dmft_weights2(Aweight, gloc, logG, dens,  Gdloc, Gdloc0, DMFTU, dNds, Edimp0, AEweight, nedim, STrans, sigma, s_oo, wkp, omega, vnorm1, DM_EF, Temperature, wgamma, gamma, csize, nbands, DM_nemin, DM_nemax, maxsize, ncix, nomega, npomega, nume, matsubara, lgTrans, lg_deg, ntcix, nipc, ikp, maxdim2, iorbital, iSx, norbitals)

        Nds = Nds + dNds
        
        call cputim(t2c)
        call walltim(t2w)
        time_dmfwgh=time_dmfwgh+t2c-t1c
        time_dmfwgh_w=time_dmfwgh_w+t2w-t1w

        CALL GetLocalDensityMatrix(dDM, DMFTU(:,:,:,1), Aweight, vnorm1, iorbital, iSx, nbands, maxdim2, norbitals, lmaxp)
        DM = DM + dDM

        DEALLOCATE( STrans )
        deallocate( lgTrans )
        DEALLOCATE( DMFTU )
        
        ! We use the following formula to implement the total energy:
        ! dE = \sum_{k,i} E_{k,i}^{LDA} * n_{i,i}^{DMFT}
        ! where E_{k,i} are Kohn-Sham eigenvalues
        ! and n_{i,i}^{DMFT}= T*sum_{iom}g_{k,ii}^{DMFT} is the DMFT density matrix in KS base
        DO i=nemin,DM_nemin-1
           ETOT2 = ETOT2 + E(i)*dble(wkp)*vnorm1
        ENDDO
        DO i=1,nbands
           ETOT2 = ETOT2 + E(i+DM_nemin-1)*dble(Aweight(i,i))*vnorm1
        ENDDO
        
        TDE1=0
        DO i=nemin,DM_nemin-1
           dens = dens + dble(wkp)*vnorm1
           TDE1 = TDE1 + FreeE0(E(i)-DM_EF,Temperature)*dble(wkp)*vnorm1
        ENDDO
        TDE2=0
        DO i=1,nbands
           TDE2 = TDE2 + FreeE0(E(i+DM_nemin-1)-DM_EF,Temperature)*dble(wkp)*vnorm1
        ENDDO
        FTOT2 = FTOT2 + TDE1+TDE2+logG
        FTOT0 = FTOT0 + TDE1+TDE2
        !WRITE(6,'(A,f11.6,1x,A,f11.6,1x,A,f11.6,1x,A,f11.6,1x,A,f11.6,1x,A,f11.6)') 'Density=', dens, 'ETOT(band)=', ETOT2, 'trlgG+mu*N=', FTOT2+DM_EF*dens, 'logG=',logG,'TDE=', TDE1+TDE2, 'trlgG=', FTOT2
        
        allocate( zw2(nbands) )
        CALL Diagonalize_DMFT_WEIGHTS(zw2, Aweight, nbands, DM_nemin, DM_nemaxx, .True.)

        ! Weights for all bands, which are needed below
        weight=0
        weight(:DM_nemin-1) = zw0(:DM_nemin-1)  ! These are fully occupied
        weight(DM_nemin:DM_nemax) = zw2(:)      ! These have DMFT weights
        if (recomputeEF.LT.0) then ! This is just for debugging of DMFT
           ! We take DFT weights are replace DMFT weights with those computed on LDA level.
           weight(:nemax) = lda_weight(ikp,:nemax)
           DM_nemaxx = lda_nbmax
           !Aweight(:,:) = 0.0d0
           !do num=1,nbands
           !   Aweight(num,num)=1.d0
           !enddo
        endif

        nbands_dft = 1
        if (Qforce) nbands_dft = nemax-nemin+1
        allocate( wEpsw(nbands_dft, nbands_dft) )
        allocate( zwe(nedim) )
        if (Qforce) then
           CALL GetAEweight(wEpsw, Aweight, AEweight, zw2, nbands, nbands_dft, nemin, DM_nemin, DM_nemaxx)
           CALL Diagonalize_DMFT_WEIGHTS(zwe, AEweight, nbands, DM_nemin, DM_nemaxx, .False.)
           
           !do i=1,DM_nemaxx-DM_nemin+1
           !   WRITE(6,'(I4,1x,6F20.10)') i, zw2(i)/zw2(1), zwe(i)/zw2(1), zwe(i)/zw2(i), E(i+DM_nemin-1), zw2(i), dble(wEpsw(i+DM_nemin-nemin,i+DM_nemin-nemin))
           !enddo
           
        endif
        
        call cputim(t3c)
        call walltim(t3w)
        time_ed=time_ed+t3c-t2c
        time_edw=time_edw+t3w-t2w
        
        ! we need to have alm up to DM_nemax
        nemax     = max(DM_nemax,nemax)
        do num=nemin,DM_nemaxx
           xwt=xwt+weight(num)*vnorm1
        enddo

        WRITE(*,'(I3,A,1x,I3,1x,I4,1x,A,I4,1x,A,I4,1x,A,I4,1x,A,I4,1x,A,I4,1x)') myrank, ') Finished k-point', ikp, iikp, 'with nemin=', nemin, 'nemax=', nemax, 'DM_nemin=', DM_nemin, 'DM_nemax=', DM_nemax, 'DM_nemaxx=', DM_nemaxx
     
        deallocate( zw2 )
        
        nbandsx = DM_nemaxx-DM_nemin+1

        CALL cmp_MT_density(w_RHOLM, Aweight, wEpsw, weight, fsph, fnsp, fsph2, nemin, nemax, DM_nemin, DM_nemax, DM_nemaxx, nbands, nbandsx, nbands_dft, lm_max, n0, nnlo, coord, ikp, DM_EF, time_bl,time_bl_w, time_reduc, time_reduc_w, time_radprod, time_radprod_w, time_m, time_m_w, time_rad, time_rad_w, time_ilm, time_ilm_w, time_lo, time_lo_w, time_force, time_force_w)
        
        DEALLOCATE( wEpsw )
        CALL cputim(time1)
        CALL walltim(time1_w)

        time_mt = time_mt + time1-t3c
        time_mt_w = time_mt_w + time1_w-t3w
        
        CALL cmp_interstitial(sumfft, tloc, Aweight, AEweight, zwe, nedim, iff1, iff2, iff3, ift1, ift2, ift3, n0, nnlo, nbands, nbandsx, nemin, DM_nemin, DM_nemaxx, Qforce)
        
        DEALLOCATE( zwe )
        DEALLOCATE( AEweight )
        DEALLOCATE( Aweight )
        
        CALL cputim(time2)
        CALL walltim(time2_w)
        time_int =time_int +time2-time1
        time_intw =time_intw +time2_w-time1_w
     ENDDO     !  kpoint loop end

     DO is=1,iso    !------ over up/dn ---------------------!
        itape1=8+is
        if (vector_para) then
           close(itape1)
        else
           rewind(itape1)
        endif
     ENDDO
     
  ENDDO

  deallocate( dNds )
  
  if (abs(projector).ge.5) call p_deallocate()
  if (Qforce) CALL Deallocate_nsh()

  allocate( Edimp(ntcix) )
  call SymmetrizeLocalQuantities(Gdloc0, Edimp0, Gdloc(:,:,1), Edimp, cfX, cix_orb, cixdim, iSx, iorbital, nindo, lg_deg, norbitals, maxdim2, nomega)
  deallocate( Edimp0 )
  
  CALL Reduce_MPI(xwt, ETOT2, FTOT2, FTOT0, gloc, w_RHOLM, DM, Gdloc, fsph, fnsp, fsph2, Edimp, Nds, w_xwt1, w_xwteh, w_xwtel, w_xwt1h, w_xwt1l, sumfft, tloc, nomega, NRAD, LM_MAX, nat, natm, iff1, iff2, iff3, ift1, ift2, ift3, maxdim, ncix, ntcix, nipc, Qforce)
  
  if (myrank.eq.master) then

     CALL cputim(time1)
     CALL walltim(time1_w)
     
     call SymmetrizeDensityMatrix(DM, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2)

     write(6,*) 'Sum of eigenvalues=', ETOT2*Ry2eV, 'eV'
     write(6,*) 'TrLogG=', (FTOT2+DM_EF*elecn)*Ry2eV, 'eV'
     if (matsubara) then
        ! THIS SHOULD BE EXTENDED TO REAL AXIS
        CALL cmpLogGdloc(logGloc, eimp_nd, eimp_nd2, DeltaG, forb, TrGSigVdc, Gdloc, Edimp, s_oo, DM, Nds, Temperature, Sigind, Sigind_orig, cixdim, ncix, maxdim, maxsize, ntcix, sigma, nomega, csize, fh_sig, nipc, SOlapm)
     else
        logGloc=0.
        DeltaG=0.
        eimp_nd=0.
        eimp_nd2=0.
     endif
     
     CALL cputim(time2)
     CALL walltim(time2_w)
     time_lnG = time2-time1
     time_lnGw =time2_w-time1_w
     
     logGloc = logGloc*2./iso
     DeltaG = DeltaG*2./iso
     WRITE(6,*) 'TrLog(Gloc)[eV]=', logGloc*Ry2eV
     WRITE(6,*) 'Eimp*Nd[eV]=', eimp_nd*Ry2eV
     WRITE(6,*) '(Eimp+s_oo)*Nd[eV]=', eimp_nd2*Ry2eV
     WRITE(6,*) 'Tr((Delta-w*dDelta)G)=', DeltaG*Ry2eV
     WRITE(6,'(A,f15.12,1x,A,f15.10)') 'Ratio to renormalize=', (elecn/xwt), 'rho-rho_expected=', (xwt-elecn)
     if (Rho_Renormalize) then
        w_RHOLM = w_RHOLM * (elecn/xwt)
     endif
     
     do iom=1,nomega
        WRITE(fh_dos,'(f18.10,1x,2f18.10)') omega(iom)*Ry2eV, -dimag(gloc(iom))/(pi*Ry2eV)
     enddo

     open(1999,file='gdloc.dat',status='unknown')
     WRITE(1999,'(A)',advance='no') '#'
     do i=1,ntcix
        WRITE(1999,'(f10.5,1x)',advance='no') Edimp(i)
     enddo
     WRITE(1999,*)
     do i=1,nomega
        WRITE(1999,'(f18.12,2x)',advance='no') omega(i)*Ry2eV
        do it=1,ntcix
           WRITE(1999,'(f16.12,1x,f16.12,2x)',advance='no') dble(Gdloc(it,i,1))/Ry2eV, aimag(Gdloc(it,i,1))/Ry2eV
        enddo
        WRITE(1999,*)
     enddo
     close(1999)
     
  endif

  deallocate( lg_deg )
  deallocate( Nds )
  deallocate( Gdloc, Edimp )
  DEALLOCATE( gloc )
  deallocate( Gdloc0 )

  if (myrank.eq.master) then
     if (Qforce) then
        CALL init_V_vns(filename_V_vns, 9919)  ! Start reading non-spherical potential
        CALL init_V_vsp(filename_V_sph, 18, ISCF)    ! Start reading spherical potential
        allocate( Vlm(nrad,0:lm_max-1) )
        allocate( Vr(nrad) )
     endif
     if (Qforce2 .and. Qforce) then
        allocate( rh3_tmp(ncom,nat) )
        !allocate( drh3_tmp(ncom,nat) )
     endif
     DO jatom=1,nat
        lm(1:2,1:NCOM) = w_lm(1:2,1:NCOM,jatom)
        LMMAX = w_lmmax(jatom)
        R(:)=0.d0
        DO i=1,jri(jatom)
           R(i)=r0(jatom)*Exp((i-1.)*dx(jatom)) ! Radial mesh
        ENDDO
        IMAX=JRI(JATOM)
        
        CALL cputim(time1)
        CALL walltim(time1_w)
        
        DO ILM=1,LMMAX
           OVRDIV=1.0
           Li=LM(1,ILM)
           Mi=LM(2,ILM)
           IF(Li.NE.0) THEN
              IF(Mi.NE.0) w_RHOLM(:IMAX,ILM,jatom)=w_RHOLM(:IMAX,ILM,jatom)*SQRT2  ! all but m=0 term is normalized by sqrt(2)
              CYCLE
           ENDIF
           ! IF IT COMES HERE, Li==0
           CALL CHARGE(R,OVRDIV,SQFP,w_RHOLM(1,ILM,jatom),DX(JATOM),JRI(JATOM),TC)
           WRITE(6,207)  JATOM,TC
           WRITE(21,720) jatom,jatom,TC
           IF(ISPLIT(JATOM).EQ.1) THEN
              WRITE(6,250) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,250) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.2) THEN                               
              WRITE(6,251) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,251) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.3) THEN                               
              WRITE(6,252) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,252)jatom,JATOM, JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.-2) THEN                              
              WRITE(6,253) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,253) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.4) THEN                               
              WRITE(6,254) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,254) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.5) THEN                               
              WRITE(6,255) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,255) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.6) THEN                               
              WRITE(6,256) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,256) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.8) THEN                               
              WRITE(6,257) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
              WRITE(21,257) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,14)
           ELSE IF(ISPLIT(JATOM).EQ.15) THEN                               
              WRITE(6,258) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,21)
              WRITE(21,258) jatom,JATOM,JATOM,(w_XWT1(I,jatom),I=0,3),(w_XWT1(I,jatom),I=7,21)
           ENDIF
           do i=0,3
              if(w_xwt1l(i,jatom).lt.0.00005d0) then
                 w_XWT1l(i,jatom)=0.00001d0
                 w_XWTel(i,jatom)=w_XWT1l(i,jatom)*10.00001d0
              endif
              if(w_xwt1h(i,jatom).lt.0.00005d0) then
                 w_XWT1h(i,jatom)=0.00001d0
                 w_XWTeh(i,jatom)=w_XWT1h(i,jatom)*10.00001d0
              endif
           enddo
           write(6,260)
           WRITE(6,261) jatom,(w_XWT1l(I,jatom),w_XWTel(I,jatom)/w_XWT1l(I,jatom),I=0,3)
           write(6,262)
           WRITE(6,263) jatom,(w_XWT1h(I,jatom),w_XWTeh(I,jatom)/w_XWT1h(I,jatom),I=0,3)
           write(21,260)
           WRITE(21,261) jatom,(w_XWT1l(I,jatom),w_XWTel(I,jatom)/w_XWT1l(I,jatom),I=0,3)
           write(21,262)
           WRITE(21,263) jatom,(w_XWT1h(I,jatom),w_XWTeh(I,jatom)/w_XWT1h(I,jatom),I=0,3)
           !
           ! if li==0 and mi!=0 .... should this ever happen?
           IF(Mi.NE.0) w_RHOLM(:IMAX,ILM,jatom)=w_RHOLM(:IMAX,ILM,jatom)*SQRT2
        enddo
     
        call cputim(time3)
        call walltim(time3_w)
        
        ! Writting case.clmval
        WRITE(8,1990) JATOM
        WRITE(8,2001) LMMAX
        DO ILM=1,LMMAX
           Li=LM(1,ILM)
           Mi=LM(2,ILM)
           WRITE(8,2011) Li,Mi
           WRITE(8,2022) ( w_RHOLM(I,ILM,jatom), I=1,IMAX )
           WRITE(8,2031)
        end DO
        WRITE(8,2030)
        
        call cputim(time4)
        call walltim(time4_w)
        time_writeclm=time_writeclm+time4-time3
        time_writeclm_w=time_writeclm_w+time4_w-time3_w

        if (Qforce) then
           ! Old FOMAI2 
           ! calculates integral of Veff*grad(rholm)
           CALL read_V_vns(Vlm(:,0:lmmax-1),jatom,nrad,imax,lmmax)
           CALL read_V_vsp(Vr(1:imax),jatom,imax)                  
           Vlm(1:imax,0) = Vr(1:imax)/R(1:imax)*(sqrt(4.d0*pi)) ! No need for conversion back to Rydbergs
           fvdrho(:,jatom) = Force4_mine(w_rholm(:,1:lmmax,jatom),Vlm(:,0:lmmax-1),R,lmmax,lm(:,1:lmmax),imax,dx(jatom),nrad,lmax2,lmmax)

           allocate( rh1_tmp(nrad,ncom), rh2_tmp(nrad,ncom) )
           rh1_tmp(:,:)=0.d0
           rh2_tmp(:,:)=0.d0
           do ilm1=1,lmmax
              rh1_tmp(1:imax,ilm1) = Vlm(1:imax,ilm1-1)*R(1:imax)**2
              rh2_tmp(1:imax,ilm1) = w_rholm(1:imax,ilm1,jatom)/R(1:imax)**2
           enddo
           fdvrho(:,jatom) = Force4_mine(rh1_tmp,rh2_tmp,R,lmmax,lm(:,1:lmmax),imax,dx(jatom),nrad,lmax2,lmmax)
           deallocate( rh1_tmp, rh2_tmp )
        end if
        if (Qforce2 .and. Qforce) then
           ! It turns out that w_rholm is going to be changed below, hence we will save its value at the MT-boundary for later
           rh3_tmp(1:lmmax,jatom) = w_rholm(imax,1:lmmax,jatom)/R(imax)**2  ! Density at the MT-boundary
        endif
        

        !print*, 'LM before LABEL='
        !do ilm1=1,lmmax
        !   print*, lm(1,ilm1), lm(2,ilm1)
        !end do

        DO ilm1=1,lmmax
           IF (LM(1,ILM1).EQ.2.AND.LM(2,ILM1).EQ.0) THEN  ! This changes quadrupolar moment y_{2,0} component of rho_lm
              w_RHOLM(1:IMAX,ILM1,jatom)=w_RHOLM(1:IMAX,ILM1,jatom)/R(1:IMAX)**3
              LABEL=1
              DO ILM2=ILM1+1,LMMAX
                 IF (LM(1,ILM2).EQ.2.AND.LM(2,ILM2).EQ.2) THEN ! This changes y_{2,2,+}
                    w_RHOLM(1:IMAX,ILM2,jatom)=w_RHOLM(1:IMAX,ILM2,jatom)/R(1:IMAX)**3
                    LABEL=2
                    DO ILM3=ILM2+1,LMMAX
                       IF (LM(1,ILM3).EQ.-2.AND.LM(2,ILM3).EQ.2) THEN  ! This changes y_{2,2,-}
                          w_RHOLM(1:IMAX,ILM3,jatom)=w_RHOLM(1:IMAX,ILM3,jatom)/R(1:IMAX)**3
                          LABEL=3
                          GOTO 7372
                       ENDIF
                    ENDDO
                    GOTO 7372
                 ENDIF
              ENDDO
              GOTO 7372
           ENDIF
        ENDDO
        
        GOTO 7371
     
7372 CONTINUE 
     
        EFGFACT=0.02997925
        FAC=0.8D0*PI*324.14D0
        DO I1=0,jri(jatom)-5,1000   ! 7374
           ! integrates: TC = integrate(R*RHOLM)
           CALL CHARGE(R,1.D0,1.D0,w_RHOLM(1,ILM1,jatom),DX(JATOM),JRI(JATOM)-I1,TC) 
           TC=TC*FAC
           EFG20=TC
           IF (LABEL.GT.1) THEN
              CALL CHARGE(R,1.D0,1.D0,w_RHOLM(1,ILM2,jatom),DX(JATOM),JRI(JATOM)-I1,TC) 
              TC=TC*FAC
              EFG22=TC
              IF (LABEL.EQ.3) THEN
                 CALL CHARGE(R,1.D0,1.D0,w_RHOLM(1,ILM3,jatom),DX(JATOM),JRI(JATOM)-I1,TC)
                 TC=TC*FAC
                 EFG2M=TC
                 QXX=(-EFG20/SQRT(3.D0)+EFG22)*SQRT(15.D0/4.D0/PI)
                 QYY=(-EFG20/SQRT(3.D0)-EFG22)*SQRT(15.D0/4.D0/PI)
                 QZZ=(EFG20*2.D0/SQRT(3.D0))*SQRT(15.D0/4.D0/PI)
                 QXY=EFG2M*SQRT(15.D0/4.D0/PI)
                 QXX=-QXX*EFGFACT
                 QXY=-QXY*EFGFACT
                 QYY=-QYY*EFGFACT
                 QZZ=-QZZ*EFGFACT
                 IF (I1.EQ.0) WRITE(6,3019)
                 IF (I1.EQ.0) WRITE(21,3019)
                 WRITE(6,3121) jatom,QXX,QXY,QYY,QZZ,R(JRI(JATOM)-I1)
                 IF (I1.EQ.0) WRITE(21,3121) jatom,QXX,QXY,QYY,QZZ,R(JRI(JATOM)-I1)      
              ELSE
                 VXX=(-EFG20/SQRT(3.D0)+EFG22)*SQRT(15.D0/4.D0/PI)
                 VYY=(-EFG20/SQRT(3.D0)-EFG22)*SQRT(15.D0/4.D0/PI)
                 VZZ=(EFG20*2.D0/SQRT(3.D0))*SQRT(15.D0/4.D0/PI)
                 IF (I1.EQ.0) WRITE(6,3017)
                 VXX=-VXX*EFGFACT
                 VYY=-VYY*EFGFACT
                 VZZ=-VZZ*EFGFACT
                 IF (I1.EQ.0) WRITE(21,3017)
                 WRITE(6,3021) jatom,VXX,VYY,VZZ,R(JRI(JATOM)-I1)
                 IF (I1.EQ.0) WRITE(21,3021) jatom,VXX,VYY,VZZ,R(JRI(JATOM)-I1)
              ENDIF
           ELSE
              EFG20=EFG20*2.D0/SQRT(3.D0)*SQRT(15.D0/4.D0/PI)
              EFG20=-EFG20*EFGFACT
              WRITE(6,7207)  JATOM,EFG20,R(JRI(JATOM)-I1)
              IF (I1.EQ.0) WRITE(21,7720) jatom,JATOM,EFG20,R(JRI(JATOM)-I1)
           ENDIF
        ENDDO
     
     
7371 CONTINUE                                     
     
        REWIND fh_vec
        
        CALL cputim(time3) 
        call walltim(time3_w)
        time_writeefg=time_writeefg+time3-time4
        time_writeefg_w=time_writeefg_w+time3_w-time4_w
        
        
        CALL cputim(time2)
        CALL walltim(time2_w)
        time_writescf=time_writescf+time2-time3
        time_write=time_write+time2-time1
        time_writescf_w=time_writescf_w+time2_w-time3_w
        time_write_w=time_write_w+time2_w-time1_w
     ENDDO
     if (Qforce) then
        deallocate( Vlm, Vr )
        CALL close_V_vns()
        CALL close_V_vsp()
     endif
     
     write(6,*)
     write(6,'(a,2f10.1)') '   atpar              (cpu,wall):',time_atpar_c,time_atpar_w
     write(6,'(a,2f10.1)') '   readvec            (cpu,wall):',time_rd_c,time_rd_w
     write(6,'(a,2f10.1)') '   readvec, read only (cpu,wall):',time_r_c,time_r_w
     WRITE(6,'(a,2f10.1)') '   cmp_MT_density     (cpu,wall):',time_mt,time_mt_w
     WRITE(6,'(a,2f10.1)') '   alm_blm_clm             (cpu):',time_reduc,time_reduc_w
     WRITE(6,'(a,2f10.1)') '   blocked loop            (cpu):',time_bl,time_bl_w
     WRITE(6,'(a,2f10.1)') '   ilm-tot loop            (cpu):',time_ilm,time_ilm_w
     WRITE(6,'(a,2f10.1)') '   cmm_MT_force            (cpu):',time_force,time_force_w
     WRITE(6,'(a,2f10.1)') '   local orbit loop        (cpu):',time_lo,time_lo_w
     WRITE(6,'(a,2f10.1)') '   ilm-1 loop              (cpu):',time_radprod,time_radprod_w
     WRITE(6,'(a,2f10.1)') '   ilm-2 loop              (cpu):',time_m,time_m_w
     WRITE(6,'(a,2f10.1)') '   ilm-3 loop              (cpu):',time_rad,time_rad_w
     WRITE(6,'(a,2f10.1)') '   cmp_dmft_weights        (cpu):',time_dmfwgh, time_dmfwgh_w
     WRITE(6,'(a,2f10.1)') '   write-clm          (cpu,wall):',time_writeclm,time_writeclm_w
     WRITE(6,'(a,2f10.1)') '   write-scf          (cpu,wall):',time_writescf,time_writescf_w
     WRITE(6,'(a,2f10.1)') '   write-efg          (cpu,wakk):',time_writeefg,time_writeefg_w
     WRITE(6,'(a,2f10.1)') '   write              (cpu,wall):',time_write,time_write_w
     WRITE(6,'(a,2f10.1)') '   dmft renorm-overlap(cpu,wall):',time_dmf0, time_dmf0w
     WRITE(6,'(a,2f10.1)') '   dmft chemical pot  (cpu,wall):',time_dmf1, time_dmf1w
     WRITE(6,'(a,2f10.1)') '   dmft transformation(cpu,wall):',time_dmf, time_dmfw
     WRITE(6,'(a,2f10.1)') '   dmft exact-diag    (cpu,wall):',time_ed, time_edw
     WRITE(6,'(a,2f10.1)') '   interstitial charge(cpu,wall):',time_int, time_intw
     WRITE(6,'(a,2f10.1)') '   cmp_log_gdloc      (cpu,wall):',time_lnG,time_lnGw
     
     WRITE(6,881)  XWT
     WRITE(21,881) XWT

     dEtot=0
     if (recomputeEF.EQ.0) then
        dEtot = DM_EF*(elecn-xwt)  ! dE = mu*(N0-N)  ; If chemical potential is fixed, we need to correct Energy  
     endif
     CALL PrintLocalDensityMatrix(DM, 6, s_oo, cixdim, .True.)
     CALL PrintLocalDensityMatrix(DM, 21, s_oo, cixdim, .False.)

     WRITE(6,*)
     WRITE(6,'(A,f15.12,1x,A,f15.10)') 'Ratio to renormalize=', elecn/xwt, 'rho-rho_expected=', xwt-elecn
     !WRITE(21,'(A,f15.12,1x,A,f15.10)') '# Ratio to renormalize=', elecn/xwt, 'rho-rho_expected=', xwt-elecn
     WRITE(21,'(A)') '  rho-rho_expected in dmft2:'
     WRITE(21,'(A,f20.12,1x,A,f15.10)') ':DRHO  ', xwt-elecn
     if (abs(xwt-elecn).gt.0.1 ) then
        WRITE(21,'(A)') 'WARNING : Electron charge density is very different than expected. You should change recomputeEF and set it to 1 (switch in on)'
     endif

     WRITE(21,725) ETOT2+dEtot
     WRITE(21,'(A,1x,A,F20.9,3x,A)') ':XSUM :', 'SUM OF EIGENVALUES =  ', FTOT2+DM_EF*elecn-logGloc+eimp_nd+DeltaG, '# Tr(log(G))-Tr(log(G_loc))+Tr((eimp+Vdc)G)+Tr(G(D-w*dD))'
     WRITE(21,'(A,1x,A,F20.9,3x,A)') ':YSUM :', 'SUM OF EIGENVALUES =  ', FTOT2+DM_EF*elecn-TrGSigVdc,       '# Tr(log(G))+mu*N-Tr((Sig-Vdc)*G)'
     WRITE(21,'(A,1x,A,F20.9,3x,A)') ':ZSUM :', 'SUM OF EIGENVALUES =  ', FTOT2+DM_EF*elecn-logGloc+eimp_nd, '# Tr(log(G))-Tr(log(G_loc))+Tr((eimp+Vdc)G)'
     WRITE(21,'(A,1x,A,F20.9,3x,A)') ':WSUM :', 'SUM OF EIGENVALUES =  ', FTOT2+DM_EF*elecn-logGloc+DeltaG,  '# Tr(log(G))-Tr(log(G_loc))+Tr(G(D-w*dD))'
     WRITE(21,'(A,1x,A,F20.9,3x,A)') ':QSUM :', 'SUM OF EIGENVALUES =  ', FTOT2+DM_EF*elecn-logGloc+eimp_nd2,'# Tr(log(G))-Tr(log(G_loc))+Tr((eimp+Vdc+s_oo)G)'
     WRITE(21,'(A,1x,A,F20.9,3x,A)') ':VSUM :', 'SUM OF EIGENVALUES =  ', FTOT2+DM_EF*elecn, '# Tr(log(G))+mu*N'
     WRITE(6,'(A,1x,F15.6,2x,A,1x,F15.6,2x,A,1x,F15.6)') 'TrLog(-G)=', FTOT2+DM_EF*elecn, 'TrLog(-Gloc)=', logGloc
     WRITE(6,882)  ETOT2
     WRITE(6,883)  FTOT2+DM_EF*elecn-logGloc

     !if (Qforce) CALL PrintForces(forcea,fsph,fsph2,fnsp,fvdrho,forb) ! former FOMAI3
  endif

  deallocate( iorbital )
  CALL fini_charp
  CALL fini_chard
  CALL fini_charf

  REWIND fh_vec

6170 CONTINUE  

  DEALLOCATE( sigma, omega, s_oo )
  CALL fini_muzero
  DEALLOCATE( Rspin )
  deallocate (e_store,elo_store)
  DEALLOCATE( nindo, cix_orb, cixdim, iSx )
  deallocate( cfX, noccur )
  CALL fini_lohelp
  CALL fini_xa
  CALL fini_xa3
  CALL w_deallocate1

  if (Qrenormalize .and. abs(projector).le.5) then
     deallocate(Olapm0, SOlapm )
  endif

  
  CALL cputim(time2)
  TCLM = time2
  CALL walltim(time2_w)
  TCLM_w = time2_w

  if (myrank.EQ.master .and. mode.NE.'m' .and. mode.NE.'e') then
     CALL cputim(time1)
     CALL walltim(time1_w)

     ! ---- Interstitials ---
     ! Back transform cumulative grids
     ! computes rho(K)=FFT(rho(r))
     call fft_init(iff1, iff2, iff3, 1)
     fft(:,:,:)=sumfft(:,:,:)
     call fft_run()
     sumfft(:,:,:) = fft(:,:,:)/dble(iff1*iff2*iff3)
     IF(Qforce) THEN
        fft(:,:,:) = tloc(:,:,:)
        call fft_run()
        tloc(:,:,:) = fft(:,:,:)/dble(iff1*iff2*iff3)
     ENDIF
     call fft_fini()

     ALLOCATE( rho1(nwave*nsym) )
     rho1(:)=0.0 
     call fftget(nwave,iff1,iff2,iff3,rho1,sumfft,kmax,nsym,kzz,iz,tau,Qcomplex)
     
     ALLOCATE( rhok(nwave) )
     rhok(:)=0.0
     
     if(Qforce) then
        allocate( ekin1(nwave*nsym) )
        ekin1(:)=0.0
        call fftget(nwave,iff1,iff2,iff3,ekin1,tloc,kmax,nsym,kzz,iz,tau,Qcomplex)
        allocate( ekink(nwave) )
        ekink(:)=0.0;
     endif
     CALL cputim(time2)
     CALL walltim(time2_w)
     
     time_int  = time_int +time2-time1
     time_intw = time_intw +time2_w-time1_w
     !.....SUM OVER ALL RECPR. LATTIC VECTORS OF ONE STAR
     ia1=1
     volin = 1.0d0/vol
     DO j=1,nwave
        ia2=ia1+inst(j)-1
        rhok(j) = sum(  rho1(ia1:ia2)*dconjg(tauk(ia1:ia2)) )*volin
        if (Qforce) ekink(j)= sum( ekin1(ia1:ia2)*dconjg(tauk(ia1:ia2)) )*volin
        ia1=ia2+1
     ENDDO
     
     if (Rho_Renormalize) then
        rhok = rhok * (elecn/xwt)
        if (Qforce) ekink= ekink* (elecn/xwt)
     endif
     
     WRITE(8,*) '   VALENCE CHARGE DENSITY IN INTERSTITIAL '
     WRITE(8,*)
     WRITE(8,2061)  NWAVE
     WRITE(8,2071)  ( (KZZ(JX,J),JX=1,3),RHOK(J), J=1,NWAVE)
     WRITE(6,204) NWAVE,kmax
     close(8)

     if (Qforce2 .and. Qforce) then
        call cputim(time1)
        
        CALL init_V_vns(filename_V_vns, 9919)  ! Start reading non-spherical potential
        CALL init_V_vsp(filename_V_sph, 18, ISCF)    ! Start reading spherical potential
        allocate( Vlm(nrad,0:lm_max-1) )
        allocate( Vr(nrad) )
        lfirst=1
        DO jatom=1,nat
           lm(1:2,1:NCOM) = w_lm(1:2,1:NCOM,jatom)
           lmmax = w_lmmax(jatom)
           imax = jri(jatom)
           CALL read_V_vns(Vlm(:,0:lmmax-1),jatom,nrad,imax,lmmax)
           CALL read_V_vsp(Vr(1:imax),jatom,imax)
           
           R(:)=0.d0
           DO i=1,imax
              R(i)=r0(jatom)*Exp((i-1.)*dx(jatom)) ! Radial mesh
           ENDDO
           Vlm(1:imax,0) = Vr(1:imax)/R(1:imax)*(sqrt(4.d0*pi)) ! Spherically symmetric part. No need for conversion back to Rydbergs
           
           fsur2 = Force_surface_extra(rhok,ekink,Vlm(imax,0:lmmax-1),rh3_tmp(1:lmmax,jatom),kzz,pos(:,lfirst),rotij(:,:,lfirst),Rmt(jatom),BR1,rotloc(:,:,jatom),iz,tau,lm,iord,nwave,lmmax,ncom,jatom)
           fsur(:,jatom) = fsur2(1:3)
           fextra(:,jatom) = fsur2(4:6)
           fkdc(:) = fsur2(7:9)
           write(6,'(A,I3,A,3F15.7)') 'Fext_fromdrho[', jatom,']=', (fvdrho(:,jatom) + fdvrho(:,jatom))*1000
           write(6,'(A,I3,A,3F15.7)') 'Fext_differ  [', jatom,']=', fextra(:,jatom)*1000.d0
           write(6,'(A,I3,A,3F15.7)') 'Fsph2_interst[', jatom,']=', fkdc(:)*1000.d0
           write(6,'(A,I3,A,3F15.7)') 'Tr(dV/dr*rho)[', jatom,']=', fdvrho(:,jatom)*1000

           !write(6,'(A,F7.4,A)') 'Radial Derivative of the Density at the MT-boundary r=', Rmt(jatom), ' from augmented PW:'
           !DO lm1=1,lmmax
           !   write(6,'(I3,1x,I3,1x,F20.10)') lm(1,lm1), lm(2,lm1), drh3_tmp(lm1,jatom)
           !END DO
           
           lfirst = lfirst+mult(jatom)
        ENDDO

        CALL PrintForces2(forcea,fsph,fsph2,fnsp,fvdrho,forb,fextra) ! former FOMAI3
        
        DO jatom=1,nat
           fsur_norm = sqrt(sum(fsur(:,jatom)**2))
           WRITE(6, '(2i3,a7,4f17.9)') jatom,1,'+ SUR', fsur_norm*1000, fsur(:,jatom)*1000
           WRITE(21,'(A,i3.3,A,1x,i3,A,4F17.9) ') ':FSU', jatom, ':', jatom, '.ATOM', fsur_norm*1000, fsur(:,jatom)*1000
        ENDDO
        
        deallocate( Vlm, Vr )
        CALL close_V_vns()
        CALL close_V_vsp()
        call cputim(time2)
        write(6,*) 'time of fsumai0', time2-time1
        deallocate( rh3_tmp )
        !deallocate( drh3_tmp )
     endif
     
     CALL w_deallocate2
     DEALLOCATE(rhok)
     DEALLOCATE(rho1)

     if(Qforce .and. .not. Qforce2) then
        CALL PrintForces(forcea,fsph,fsph2,fnsp,fvdrho,forb) ! former FOMAI3
        call cputim(time1)
        call Force_surface(ekink,kzz,pos,rotij,mult,Rmt,BR1,rotloc,iz,tau,iord,nwave,nat,natm,Qcomplex)
        call cputim(time2)
        write(6,*) 'time of fsumai1', time2-time1 
        DEALLOCATE(ekin1, ekink)
     endif
  endif

  call fini_forces()
  
  if (mode.ne.'m' .and. mode.ne.'e') then
     DEALLOCATE ( sumfft, tloc )
  endif
  

  DEALLOCATE( pr_procs )
  ! --- Interstitials ---

  CALL cputim(time2)
  TFOUR = time2
  CALL walltim(time2_w)
  TFOUR_w = time2_w

  RETURN
  

207 FORMAT(1H0,' TOTAL CHARGE INSIDE SPHERE',I5,1H:,F12.6)            
250 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PZ,PXY',/,':QTL',i3.3,':',12F7.4)               
251 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,      ','D-EG,D-T2G ',/,':QTL',i3.3,':',12F7.4)  
252 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,      ','D-Z2,D-XY,X2Y2,D-XZ,YZ ',/,':QTL',i3.3,':',12F7.4)
253 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PZ,PXY,,','D-Z2,D-XY,D-X2Y2,D-XZ,YZ ',/,':QTL',i3.3,':',12F7.4)                          
254 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PZ,PXY, ','D-Z2,D-XY,X2Y2,D-XZ,YZ ',/,':QTL',i3.3,':',12F7.4)                            
255 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,      ','D-Z2,D-X2Y2,D-XY,D-XZ,D-YZ ',/,':QTL',i3.3,':',12F7.4)                        
256 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PX,PY,PZ,',' ',/,':QTL',i3.3,':',12F7.4)
257 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PX,PY,PZ,','D-Z2,D-X2Y2,D-XY,D-XZ,D-YZ ',/,':QTL',i3.3,':',12F7.4)                        
258 FORMAT(':PCS',i3.3,':',1X,'PARTIAL CHARGES SPHERE =',I3,' S,P,D,F,PX,PY,PZ,D-Z2,D-X2Y2,D-XY,D-XZ,D-YZ,F00,F11,F22',',F33,','F1M,F2M,F3M ',/,':QTL',i3.3,':',19F7.4)
260 FORMAT(8x,'Q-s-low E-s-low   Q-p-low E-p-low   Q-d-low E-d-low','   Q-f-low E-f-low')
261 FORMAT(':EPL',i3.3,':',4(2F8.4,2x))
262 FORMAT(8x,'Q-s-hi  E-s-hi    Q-p-hi  E-p-hi    Q-d-hi  E-d-hi ','   Q-f-hi  E-f-hi ')
263 FORMAT(':EPH',i3.3,':',4(2F8.4,2x))
720 FORMAT(/,':CHA',i3.3,':',1X,'TOTAL CHARGE INSIDE SPHERE ',I3,' = ',F12.6)         
725 FORMAT(/,':SUM  :',1X,'SUM OF EIGENVALUES =  ',F20.9,/)            
!726 FORMAT(':XSUM :',1X,'SUM OF EIGENVALUES =  ',F20.9,/)            
787 FORMAT('    VALENCE CHARGE DENSITY   IN  MT SPHERES',5X,I3,' ITERATION')
800 FORMAT(/,10X,68(1H-),/,11X,'W A V E F U N C T I O N S','   AND   C H A R G E S   IN   S P H E R E S',/,10X,68(1H-))                                         
881 FORMAT(/,':CHA  :',' TOTAL CHARGE INSIDE UNIT CELL =',F15.6) 
882 FORMAT(//,3X,'SUM OF EIGENVALUES:',9X,F15.6)                    
883 FORMAT(3X,'Tr(log(G))-Tr(log(Gloc)):',3X,F15.6/)                    
1003 FORMAT(251(I3,I2))                                                 
1005 FORMAT(/,7X,'LMMAX',I3,/,7x,'LM= ',17(i3,I2),(/,7x,18(i3,i2)))          
1990 FORMAT(3X,'ATOMNUMBER   ',I3,5X,10A4)                             
2001 FORMAT(3X,'NUMBER OF LM',I3//)                                   
2011 FORMAT(3X,'CLM(R) FOR L',I3,3X,'M=',I2/)                         
2022 FORMAT(3X,4E19.12)                                                 
2030 FORMAT(///)                                                       
2031 FORMAT(/)                                                         
!2032 FORMAT(49X,I3,//)  
3017 FORMAT(/,22X,'VXX',9X,'VYY',9X,'VZZ',7X,'UP TO R',/)
3019 FORMAT(/,22X,'QXX',9X,'QXY',9X,'QYY',9X,'QZZ',7X,'UP TO R',/)
3121 FORMAT(':VZZ',i3.3,':',8x,4F12.5,F12.3)
3021 FORMAT(':VZZ',i3.3,':',8x,3F12.5,F12.3)
13   FORMAT(////,':POS',i3.3,':',1x,'AT.NR.',I4,1X,'POSITION =',3F8.5,2X,'MULTIPLICITY =',I3)
7207 FORMAT(1H ,' EFG INSIDE SPHERE',I5,1H:,F12.6,5X,' UP TO R = ',F10.7)                                                          
7720 FORMAT(':VZZ',i3.3,':',1X,'EFG INSIDE SPHERE ',I3,' = ',F12.6,5X,' UP TO R =',F10.5)


204 FORMAT(I10,' FOURIER COEFFICIENTS CALCULATED',/' LARGEST COMPONENTS:',3I6)
2061 FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')
2071 FORMAT(3X,3I5,2E19.12)
1060 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'=',F8.3,//   ,':FER  :',1X,'F E R M I - ENERGY',11X,'= ',F9.5)            

END SUBROUTINE L2MAIN


SUBROUTINE cmp_interstitial(sumfft, tloc, Aweight, AEweight, zwe, nedim, iff1, iff2, iff3, ift1, ift2, ift3, n0, nnlo, nbands, nbandsx, nemin, DM_nemin, DM_nemaxx, Qforce)
  ! INTERSTITIAL CHARGE
  ! FFT eigenvectors and accumulate weighted square in sumfft
  USE param, ONLY: nmat
  USE xa3, ONLY: k3, bk3, As
  USE dmfts, ONLY: iso
  USE xa,  ONLY: weight, E
  USE structure,ONLY: BR1
  USE fftw3_omp, ONLY: fft_init, fft_fini, fft_run, fft
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: sumfft(iff1,iff2,iff3), tloc(ift1,ift2,ift3)
  COMPLEX*16, intent(in)   :: Aweight(nbands,nbands), AEweight(nedim,nedim)
  REAL*8, intent(in)       :: zwe(nedim)
  INTEGER, intent(in)      :: iff1, iff2, iff3, ift1, ift2, ift3, n0, nnlo, nbands, nbandsx, nemin, DM_nemin, DM_nemaxx, nedim
  LOGICAL, intent(in)      :: Qforce
  ! locals
  INTEGER :: is, num, ii !, isig, ierr, csize, dsize
  COMPLEX*16, ALLOCATABLE :: Asc(:,:), Asq(:,:)
  complex*16 :: cone, czero
  cone  = cmplx(1.d0, 0.d0, 8)
  czero = cmplx(0.d0, 0.d0, 8)
  
  !isig=-1
  !csize = iff1+iff2+iff3
  !dsize = 4*(iff1+iff2+iff3)+15
  !allocate( fft(iff1,iff2,iff3), dwork(dsize),cwork(csize) )
  call fft_init(iff1, iff2, iff3, -1)
  
  ALLOCATE( Asc(n0-nnlo,1:DM_nemaxx) )
  if (Qforce) ALLOCATE( Asq(n0-nnlo,1:DM_nemaxx) )
  
  do is=1,iso
     Asc(:,:) = As(:n0-nnlo,:DM_nemaxx,is)
     CALL zgemm('N','N',n0-nnlo,nbandsx,nbands,cone,As(1,DM_nemin,is),nmat,Aweight,nbands,czero,Asc(1,DM_nemin),n0-nnlo)

     if (Qforce) then
        Asq(:,:) = As(:n0-nnlo,1:DM_nemaxx,is)
        CALL zgemm('N','N',n0-nnlo,nbandsx,nbands,cone,As(1,DM_nemin,is),nmat,AEweight,nbands,czero,Asq(1,DM_nemin),n0-nnlo)
     endif
     
     DO num=nemin,DM_nemaxx
        !  puts eigenvector into fft array for FFT
        CALL fft1set(n0-nnlo,iff1,iff2,iff3,Asc(:,num),fft,k3,nmat)   ! fft becomes equivalent to Asc, which is the eigenvector, transformed to the basis with diagonal DMFT density matrix.
        !  computes a(r) = FFT(a(K))
        call fft_run()
        !  |a(r)|^2 w(k)
        CALL fftsumup(iff1,iff2,iff3,weight(num)/iso,fft, sumfft) ! here only sumfft is changed (updated)
        IF (Qforce) THEN
           if (num.LT.DM_nemin) then
              CALL fftsumup(iff1,iff2,iff3,-weight(num)*E(num)/iso,fft,tloc)                !  tloc <--  -E(i)*|a(r)|^2 w(i)
           else
              CALL fft1set(n0-nnlo,iff1,iff2,iff3,Asq(:,num),fft,k3,nmat)                   ! fft becomes equivalent Asq, which is transformed to the basis which diagonalizes (rho*E) type of density matrix.
              call fft_run()
              CALL fftsumup(iff1,iff2,iff3, -zwe(num-DM_nemin+1)/iso,fft, tloc)                        ! zwe is like E(num)*weight(num)    
           endif
           DO ii=1,3
              CALL fft2set(n0-nnlo,iff1,iff2,iff3,Asc(:,num),fft,ii,k3,bk3,BR1,nmat)     ! Puts into fft <- (k+K)*a(K,r)
              call fft_run()
              CALL fftsumup(iff1,iff2,iff3, weight(num)/iso, fft, tloc)                  ! |ka(r)|^2 *w(i)
           ENDDO
        END IF
     ENDDO
  enddo
  call fft_fini()
  DEALLOCATE( Asc )
  if (Qforce) DEALLOCATE( Asq )
END SUBROUTINE cmp_interstitial


SUBROUTINE Debug_Print_Projector(DMFTU,nindo,norbitals,nbands,maxdim2,nipc)
  IMPLICIT NONE
  complex*16, intent(in)  :: DMFTU(nbands,maxdim2,norbitals,nipc)
  INTEGER, intent(in)     :: nindo(norbitals), norbitals, nbands, maxdim2, nipc
  ! locals
  INTEGER :: iorb, ind, iband
  open(988,FILE='U_2.dat',form='formatted',status='unknown',access='append')
  do iorb=1,norbitals
     do ind=1,nindo(iorb)
        WRITE(988,*) 'iorb=', iorb, 'ind=', ind
        do iband=1,nbands
           WRITE(988,'(2F20.14,1x)',advance='no') DMFTU(iband,ind,iorb,1)
        enddo
        WRITE(988,*)
     enddo
  enddo
  close(988)
END SUBROUTINE Debug_Print_Projector

     


        
