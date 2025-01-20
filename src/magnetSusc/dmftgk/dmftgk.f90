! @Copyright 2007 Kristjan Haule
program print_gk
  USE com_mpi, ONLY: myrank, master, nprocs, start_MPI, stop_MPI, nargs, argv, Qprint, cpuID, fUdmft
  IMPLICIT NONE
  ! functions
  INTEGER :: CountSelfenergy
  ! variables
  CHARACTER*100 :: case, fenergy, fbasicArrays, fklist, frotlm
  CHARACTER*2   :: updn, so
  CHARACTER*1   :: mode
  INTEGER :: fhb, fhi
  CHARACTER*100 :: STR, feigenvals, fUR,fUL, filename
  INTEGER :: nom, nk, iso, nat, norbitals, ncix, natom, nkpt, nmat, nume, lmax2, maxdim2, maxdim, maxsize, icix
  LOGICAL :: Qcomplex, matsubara
  REAL*8  :: gamma, gammac, BRX(3,3), emin, emax
  INTEGER :: i, j, l, imatsubara, slen, ierr!, nargs
  REAL*8, PARAMETER      :: Ry2eV = 13.60569193
  CHARACTER*100, allocatable :: fsigname(:), fgkout(:), fglocout(:), fgloc(:)
  ! Functions
  CHARACTER*100 :: ReadWord
  CALL start_MPI()
  if (nargs.LT.1) then
     CALL stop_MPI()
     if (Qprint) then
        print *, 'Missing the input file to proceed. Create input file with the following lines'
        print *, 'g/e              # mode'
        !print *, 'BasicArrays.dat  # filename for projector'
        print *, '0/1              # mastubara or not'
        print *, 'case.energy      # wien2k energy file case.energy[so]'
        print *, 'case.klist[_band]# wien2k klist compatible with Udmft.proc'
        print *, 'case.rotlm       # wien2k rotlm file'
        print *, 'Udmft.0         # dmft projector'
        print *, 'float            # broadening for non-correlated states from case.indmfl'
        print *, 'float            # broadening for correlated states from case.indmfl'
        print *, 'sig.inp1         # input self-energies, could be several'
        print *, '------- for mode==e ----------------'
        print *, 'eigenvalues.dat  # filename with output frequency dependent eigenvalues'
        print *, 'UR.dat           # right eigenvector of non-hermitian DMFT problem'
        print *, 'UL.dat           # left eigenvector of non-hermitian DMFT problem'
        print *, '-float           # lower cutoff for output energies measured from EF'
        print *, ' float           # upper cutoff for output energies measured from EF'
        print *, '------- for mode==g ----------------'
        print *, 'G_k1             # output greens function'
        print *, 'G_local1         # output local greens function'
        print *, 'g_local1         # output diagonal local greens function'
     endif
     STOP 'ERROR dmftgk: missing arguments'
  endif

!!! Eventually I should read this on master only, and Bcast
  fhi = 995
  open(fhi, file=argv(1), status='old', ERR=900, form='formatted')
  READ(fhi,*) mode
  WRITE(*,'(A,A)') 'mode=', mode
  fbasicArrays = 'BasicArrays.dat'
  !READ(fhi, *) fbasicArrays
  !WRITE(*,'(A,A)') 'fbasicArrays=', fbasicArrays
  READ(fhi, *, iostat=ierr) imatsubara
  if (ierr /= 0) then
     print *, 'ERROR with input file. We eliminated BasicArrays.dat line and this line should contain matsubara'
     STOP 'ERROR dmftgk'
  endif
  WRITE(*,'(A,I2)') 'imatsubara=', imatsubara
  READ(fhi, *) fenergy
  WRITE(*,'(A,A)') 'fenergy=', fenergy
  READ(fhi, *) fklist
  WRITE(*,'(A,A)') 'fklist=', fklist
  READ(fhi, *) frotlm
  WRITE(*,'(A,A)') 'frotlm=', frotlm
  READ(fhi, *) fUdmft
  slen = len_trim(fUdmft)
  if (fUdmft(slen:slen)=='0') then
     fUdmft = fUdmft(:slen-1)
  endif
  WRITE(*, '(A,A)') 'fUdmft=', fUdmft
  READ(fhi, *) gamma
  WRITE(*, '(A,F10.4)') 'gamma= ', gamma
  READ(fhi, *) gammac
  WRITE(*, '(A,F10.4)') 'gammac=', gammac
!!! Eventually I should read this on master only, and Bcast
!!! Reading the first part of the header file of BasicArrays.dat to get dimensions
  fhb=996
  open(fhb,file=fbasicArrays,status='old', ERR=901, form='formatted')
  READ(fhb, *) STR ! 'nat, iso, norbitals, ncix, natom'
  READ(fhb, *) nat, iso, norbitals, ncix, natom
  READ(fhb, *) STR ! nkpt, nmat, nume
  READ(fhb, *) nkpt, nmat, nume
  READ(fhb, *) STR ! Qcomplex
  READ(fhb, *) Qcomplex
  READ(fhb, *) STR ! 'lmax2, maxdim2, maxdim, maxsize'
  READ(fhb, *) lmax2, maxdim2, maxdim, maxsize

  allocate( fsigname(ncix), fgkout(ncix), fglocout(ncix), fgloc(ncix) )
  WRITE(*,'(A,I3)') 'ncix=', ncix
  do icix=1,ncix
     fsigname(icix) = ReadWord(fhi)
     WRITE(*, '(A,A,1x)', advance='no') 'fsignamex=', TRIM(fsigname(icix))
  enddo
  READ(fhi,*)
  WRITE(*,*)

  if (mode.EQ.'g') then
     do icix=1,ncix
        fgkout(icix) = ReadWord(fhi)
        WRITE(*, '(A,A,1x)', advance='no') 'fgkoutx=', TRIM(fgkout(icix))
     enddo
     READ(fhi,*)
     WRITE(*,*)
     do icix=1,ncix
        fglocout(icix) = ReadWord(fhi)
        WRITE(*, '(A,A,1x)', advance='no') 'fglocoutx=', TRIM(fglocout(icix))
     enddo
     READ(fhi,*)
     WRITE(*,*)
     do icix=1,ncix
        fgloc(icix) = ReadWord(fhi)
        WRITE(*, '(A,A,1x)', advance='no') 'fglocx=', TRIM(fgloc(icix))
     enddo
     READ(fhi,*)
     WRITE(*,*)
  else if (mode.EQ.'e') then
     READ(fhi,*) feigenvals
     WRITE(*,'(A,A)') 'feigenvals=', feigenvals
     READ(fhi,*) fUR
     WRITE(*,'(A,A)') 'fUR=', fUR
     READ(fhi,*) fUL
     WRITE(*,'(A,A)') 'fUL=', fUL
     READ(fhi,*) emin
     WRITE(*,'(A,F10.3)') 'emin=', emin
     READ(fhi,*) emax
     WRITE(*,'(A,F10.3)') 'emax=', emax
  else
     print *, 'ERROR: Mode is not one of "e" or "g". Bailling out!  mode=', mode
     call stop_MPI()
     STOP 'ERROR dmftgk: mode not recognized. Need to be e or g'
  endif

  if (imatsubara.EQ.0) then
     matsubara = .FALSE.
  else
     matsubara = .TRUE.
  endif
  gamma = gamma/Ry2eV   ! In dmft1 we use so small broadening !

!!! Counting the number of frequency points in self-energy
  !print *, 'fsigname=', fsigname(1)
  open(81, file=fsigname(1), status='old', ERR=902, form='formatted')
  nom = CountSelfenergy(81, 1) !--- How many frequency point exists in the input file? ---!
  close(81)
  if (Qprint) WRITE(6,'(A,I4)') 'nom=', nom

  CALL ReadBR(frotlm, BRX)

  CALL PrintGk(fenergy, fsigname, fgkout, fglocout, fgloc, fklist, fhb, feigenvals, fUR, fUL, mode, emin, emax, &
       gammac, gamma, nom, nat, iso, norbitals, ncix, natom, nkpt, nmat, nume, Qcomplex, matsubara, lmax2, maxdim2, maxdim, maxsize, BRX)

  deallocate( fsigname, fgkout, fglocout, fgloc )
  !deallocate( argv )
  CALL stop_MPI()
  if (myrank.eq.master) then
     STOP 'dmftgk'
  else
     close(6)
     filename = 'nohup.'//trim(ADJUSTL(cpuID))
     call execute_command_line('rm -f '//trim(filename))
     STOP
  endif
900 print *, 'ERROR opening ',TRIM(argv(1)),' file'
901 print *, 'ERROR opening ',TRIM(fbasicArrays), ' file'
902 print *, 'ERROR opening ',TRIM(fsigname(1)), ' file'
end program print_gk


SUBROUTINE PrintGk(fenergy, fsigname, fgkout, fglocout, fgloc, fklist, fhb, feigvals, fUR, fUL, mode, emin, emax, gammac, gamma, nom, nat, iso, norbitals, ncix, natom, nkpt, nmat, nume, Qcomplex, matsubara, lmax2, maxdim2, maxdim, maxsize, BRX)
  USE com_mpi, ONLY: myrank, master, nprocs, Qprint, Reduce_MPI, cpuID, fUdmft, stop_MPI, Barrier
  IMPLICIT NONE
  CHARACTER*100, intent(in) :: fenergy, fklist, feigvals, fUR, fUL
  CHARACTER*100, intent(in) :: fsigname(ncix), fgkout(ncix), fglocout(ncix), fgloc(ncix)
  CHARACTER*1, intent(in)   :: mode
  REAL*8, intent(in)  :: emin, emax, gammac, gamma, BRX(3,3)
  INTEGER, intent(in) :: fhb, nom, nat, iso, norbitals, ncix, natom, nkpt, nmat, nume
  LOGICAL, intent(in) :: Qcomplex, matsubara
  INTEGER, intent(in) :: lmax2, maxdim2, maxdim, maxsize
  ! locals
  CHARACTER*10:: skii
  CHARACTER*100:: filename
  LOGICAL :: pform
  INTEGER :: ios
  INTEGER :: itape, i, j, k, l, is, N, nemin, nemax, ikp, ip, iq, NE, nb_min, nb_max, dir, n_min, n_max
  REAL*8  :: EF, VOL, EMIST, renorm_wgh, Ek(nume), wgh(nkpt), wg, womega, tweight, kp(3), posc(3,norbitals), phase
  INTEGER :: kvec(4,nkpt)
  integer :: fh_p, fh_gc, fh_gl, fh_gk, fh_ene, fh_ee
  integer :: numk, nsymop, nbands, isym, ddi, iw, istart, ix, iw_p, iw_m, Nd, ibnd, ibnd1, ibnd2, ind1, nxmin, nxmax
  integer :: iorb, wkii, wkis, iom, iband, nord, icix, tnorbitals, slave, ierr
  integer :: nindo(norbitals), csize(ncix), tnindo(norbitals), cixdim(ncix), nl(natom), ll(natom,4), cixdm, cixdms
  INTEGER :: cix(natom,4), iorbital(natom,4), nind(natom,4), iSx(maxdim2,norbitals), Sigind(maxdim,maxdim,ncix), iorb1, iorb2, nind1, nind2
  INTEGER :: maxdim_, norbitals_, iikp, ik_before, pr_proc, pr_procr
  COMPLEX*16 :: xomega, ww, csum, ex
  CHARACTER*3:: cpuID_
  CHARACTER*1000 :: single_line
  REAL*8,     ALLOCATABLE :: omega(:)
  COMPLEX*16, ALLOCATABLE :: STrans(:,:,:,:), DMFTU(:,:,:), gc(:)
  COMPLEX*16, ALLOCATABLE :: sigma(:,:,:), gij(:,:)
  COMPLEX*16, ALLOCATABLE :: gmk(:,:,:), gmloc(:,:,:,:), Glc(:,:,:)
  !COMPLEX*16, ALLOCATABLE :: olp(:,:,:,:)
  INTEGER,    ALLOCATABLE  :: cini(:), Sigini(:,:)
  COMPLEX*16, ALLOCATABLE :: Al(:,:), Ar(:,:), zek(:), UAl(:,:,:), UAr(:,:,:)
  ! constants
  REAL*8, PARAMETER      :: Ry2eV = 13.60569193
  REAL*8, PARAMETER      :: PI = 3.14159265358979
  COMPLEX*16, PARAMETER  :: IMAG = (0.0D0,1.0D0)
  character(len=100), dimension(3) :: fnames

  fh_p = 399   ! file for DMFT transformation UDMFT
  fh_gc  = 180 ! file for local green's function
  fh_gl = 280  ! vector form of the local Greens function
  fh_gk = 1000 ! k-dependent Green's function is fh_hk+icix
  fh_ene = 59  ! eigenvalues from case.energy[so] file
  fh_ee = 480  ! k-dependent eigenvalues, 481 left and 482 right eigenvectors

  CALL ReadKlist(fklist, nkpt, kvec, wgh)
  tweight = sum(wgh)

  ! Many index arrays are imported from DMFT
  CALL Read_The_Rest_Basic_Arrays(fhb, nindo, cixdim, nl, ll, iorbital, cix, nind, csize, iSx, Sigind, EF, VOL, posc, norbitals, ncix, natom, maxdim, maxdim2)
  
  if (Qprint) then
     WRITE(6,'(A)', advance='no') 'csize=['
     do icix=1,ncix
        WRITE(6,'(I3,A)',advance='no') csize(icix),','
     enddo
     WRITE(6,'(A)') ']'
  endif
  
  EF = EF * Ry2eV

  ! Reading the self-energy for all orbitals
  ALLOCATE( sigma(nom,maxsize,ncix), omega(nom)  )
  do icix=1,ncix
     itape = 80+icix
     open(itape, file=TRIM(fsigname(icix)), status='old', ERR=903, form='formatted')
     CALL ReadSelfenergy(itape, sigma(:,:,icix), omega, gammac, csize(icix), nom, maxsize)
  enddo

  ! Each processor reads its own Udmft.proc file
  pform = .false.  ! Should be compatible with dmft1/dmftu code. For debugging we might want to swtich pform=.true.
  filename = TRIM(fUdmft)//ADJUSTL(TRIM(cpuID))  ! fUdmft comes from def file
  if (pform) then
     open(fh_p, file=TRIM(filename), status='old', ERR=904, form='formatted')
  else
     open(fh_p, file=TRIM(filename), status='old', ERR=904, form='unformatted')
  endif
  
  if (pform) then
     READ(fh_p,*) pr_procr, nsymop, tnorbitals     ! real number of k-points in this processor
     READ(fh_p,*) (tnindo(iorb), iorb=1,norbitals)
  else
     READ(fh_p) pr_procr, nsymop, tnorbitals
     READ(fh_p) (tnindo(iorb), iorb=1,norbitals)
  endif
  
  WRITE(6,'(A,I6,2x,A,I4)') 'tot-k=', nkpt, 'pr_procr=', pr_procr

  wkii=0
  if (mode.EQ.'g') then
     allocate( gmk(maxdim,maxdim,ncix) )
     allocate( gmloc(maxdim,maxdim,ncix,nom) )

     if (myrank.eq.master) then
        DO icix=1,ncix
           cixdm = cixdim(icix)
           allocate( Sigini(cixdm,cixdm), cini(cixdm) )
           CALL GetSiginiCini(Sigini, cini, cixdms, Sigind(:,:,icix), cixdm, maxdim )
           deallocate( Sigini, cini)
           open(fh_gk+icix, file=TRIM(fgkout(icix)), status='unknown')
           !
           WRITE(fh_gk+icix, '(A,I6,1x,I3,1x,I6,1x,I3,1x,I3,1x,A)') '#', nkpt, nsymop, nom, cixdms, norbitals, ' # nkpt, nsymop, nom, cixdms norbitals'
           WRITE(fh_gk+icix, '(A)',advance='no') '#'
           do iorb=1,norbitals
              WRITE(fh_gk+icix, '(3F9.5,2x)',advance='no') (POSC(i,iorb),i=1,3)
           enddo
           WRITE(fh_gk+icix,'(A)') '  # actual position of correlated atoms in the unit cell'
        ENDDO
     else
        DO icix=1,ncix
           open(fh_gk+icix, file=TRIM(fgkout(icix))//'.'//ADJUSTL(TRIM(cpuID)) , status='unknown')
        ENDDO
     endif
  else if (mode.EQ.'e') then
     if (myrank.eq.master) then
        open(fh_ee,   file=TRIM(feigvals), status='unknown')
        open(fh_ee+1, file=TRIM(fUL), status='unknown')
        open(fh_ee+2, file=TRIM(fUR), status='unknown')
        do j=1,2
           WRITE(fh_ee+j, '(A,I7,1x,I3,1x,I6,1x,I3,3x)',advance='no') '#', nkpt, nsymop, nom, norbitals
           DO iorb1=1,norbitals
              WRITE(fh_ee+j, '(I3,1x)',advance='no') nindo(iorb1)
           ENDDO
           WRITE(fh_ee+j, '(A)') ' # nkpt, nsymop, nom, norbitals, size(iorb1),...'
           WRITE(fh_ee+j, '(A)',advance='no') '#'
           do iorb=1,norbitals
              WRITE(fh_ee+j, '(3F9.5,2x)',advance='no') (POSC(i,iorb),i=1,3)
           enddo
           WRITE(fh_ee+j,'(A)') '  # actual position of correlated atoms in the unit cell'
        enddo
     else
        open(fh_ee,   file=TRIM(feigvals)//'.'//ADJUSTL(TRIM(cpuID)), status='unknown')
        open(fh_ee+1, file=TRIM(fUL)//'.'//ADJUSTL(TRIM(cpuID)), status='unknown')
        open(fh_ee+2, file=TRIM(fUR)//'.'//ADJUSTL(TRIM(cpuID)), status='unknown')
     endif
  endif

  ! Kohn-Sham eigenvalues from the file
  open (fh_ene, file=fenergy, status='old', ERR=905, form='formatted')
  !! Reading linearization energies in Energy file
  do i=1,nat
     READ(fh_ene, fmt='(f9.5)') EMIST !---- At the beginninge we have linearization energies --!
     READ(fh_ene, fmt='(f9.5)') EMIST !---- We just skip them ---------------------------------!
  enddo

  !allocate( olp(maxdim2,maxdim2,norbitals,norbitals) )
  !olp=0.0

  pr_proc  = floor(nkpt/DBLE(nprocs)+0.999)  ! The maximum number of points calculated per processor     
  ! correcting pr_proc is the number of k-points is not dividable
  !if (Qprint)
  !WRITE(6,'(A,1x,A,I3,1x,A,I3,1x,A,I3)') cpuID,'pr_proc=', pr_proc, 'pr_procr=', pr_procr, 'tot-k=', nkpt

  ! Reading band energies from energy fileneed for those k-points that are skept on this processor
  do ik_before=1,myrank*pr_proc
     CALL ReadEnergiesK(fh_ene, Ek, kp, wg, NE, nume)
  enddo
  do iikp=1,pr_procr  ! only over k-points that are computed on this processor
     ikp = iikp + myrank*pr_proc  ! should be real index of k-point
     ! Reading band energies from energy file
     CALL ReadEnergiesK(fh_ene, Ek, kp, wg, NE, nume) ! k-points that is calculated here.

     if (abs(wgh(ikp)- wg).gt.1e-5) print *, 'ERROR: weight in case.klist and case.energy are not equal', wg, wgh(ikp)
     ! Ek was in Ry -> transform to eV
     Ek = Ek * Ry2eV

     ! Start reading DMFT transformation, which needs to be in separate files for each processor Udmft.proc
     CALL ReadDMFT_TransK_Outside(nbands, nemin, maxdim_, norbitals_, wkii, fh_p, pform, ikp)
     if (maxdim_.ne.maxdim2) print*, 'ERROR maxim=', maxdim_, 'maxdim2=', maxdim2
     if (norbitals.ne.norbitals_) print*, 'ERROR norbitals_=', norbitals_, 'norbitals=', norbitals

     nemax = nemin+nbands-1

     ALLOCATE( DMFTU(nbands,maxdim_,norbitals_) )
     ALLOCATE( STrans(maxsize,ncix,nbands,nbands) )
     allocate( gij(nbands,nbands) )

     if (mode.EQ.'e') then
        ALLOCATE( Al(nbands,nbands), Ar(nbands,nbands), zek(nbands) )
        allocate( UAl(nbands,maxdim_,norbitals_), UAr(maxdim_,nbands,norbitals_) )
     endif

     DO isym=1, nsymop
        !if (Qprint)
        WRITE(6,'(A,1x,A,I3,1x,A,I2)') cpuID,'ikp=', ikp, 'isym=', isym
        call flush(6)
        ! Continuing reading DMFT transformation
        CALL ReadDMFT_TransK_Inside(DMFTU, fh_p, pform, isym, nindo, nbands, norbitals_, maxdim_)
        ! For more efficient transformation of self-energy, we create special array of transformation STrans
        CALL CompressSigmaTransformation2(STrans, DMFTU, Sigind, iSx, cix, csize, iorbital, ll, nl, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize)

        if (mode.EQ.'g') then
           do iom=1,nom
              gij=0
              if (matsubara) then
                 xomega = omega(iom)*IMAG
              else
                 xomega = omega(iom)
              endif
              DO i=1,nbands
                 gij(i,i) = xomega+EF-Ek(i+nemin-1)  !-----  g^-1 of the LDA part in band representation ----!
              ENDDO
              CALL AddSigma_optimized2(gij, sigma(iom,:,:), STrans, csize, -1, nbands, ncix, maxsize)           

              DO i=1,nbands               !-------- adding minimum broadening for all bands -------------------!
                 gij(i,i) = gij(i,i) + (0.d0, 1.d0)*gamma
              ENDDO

              CALL zinv(gij,nbands)    !-------- inversion of matrix to get g -------------------------------!

              CALL CmpGkc2(gmk, gij, DMFTU, iSx, iorbital, ll, nl, cix, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize)
              gmloc(:,:,:,iom) = gmloc(:,:,:,iom) + gmk(:,:,:)*(wgh(ikp)/(nsymop*tweight))

              DO icix=1,ncix
                 cixdm = cixdim(icix)
                 allocate( Sigini(cixdm,cixdm), cini(cixdm) )
                 CALL GetSiginiCini(Sigini, cini, cixdms, Sigind(:,:,icix), cixdm, maxdim )
                 !
                 WRITE(fh_gk+icix,'(f14.8,2x)',advance='no') omega(iom)
                 DO ip=1,cixdms
                    do iq=1,cixdms
                       WRITE(fh_gk+icix, '(f14.8,1x,f14.8,2x)',advance='no') gmk(cini(ip),cini(iq),icix)
                    enddo
                 ENDDO
                 WRITE(fh_gk+icix,*)
                 deallocate( Sigini, cini )
              ENDDO

           enddo
        else if (mode.EQ.'e') then
           WRITE(fh_ee,'(A,5I5,2x,F23.20,2x,A)') '!', ikp, isym, nbands, nemin, nom, wgh(ikp)/(nsymop*tweight), ': ikp, isym, nbands nemin, nomega, kweight'
           do iom=1,nom
              gij=0
              DO i=1,nbands
                 gij(i,i) = Ek(i+nemin-1)  !-----  g^-1 of the LDA part in band representation ----!
              ENDDO
              CALL AddSigma_optimized2(gij, sigma(iom,:,:), STrans, csize, 1, nbands, ncix, maxsize)
              CALL eigsys(gij, zek, Al, Ar, nbands)
              UAl(:,:,:)=0
              UAr(:,:,:)=0
              DO iorb1=1,norbitals
                 nind1 = nindo(iorb1)
                 call zgemm('N','N', nbands,  nind1,  nbands, (1.d0,0.d0), Al, nbands, DMFTU(:,:,iorb1), nbands, (0.d0,0.d0), UAl(:,:,iorb1), nbands)
                 call zgemm('C','N', nind1,  nbands,  nbands, (1.d0,0.d0), DMFTU(:,:,iorb1), nbands, Ar, nbands, (0.d0,0.d0), UAr(:,:,iorb1), maxdim2)
              ENDDO

              nxmin=1
              nxmax=1
              do ibnd=1,nbands
                 if ( dreal(zek(ibnd))-EF < emin ) nxmin=ibnd
                 if ( dreal(zek(ibnd))-EF <= emax ) nxmax=ibnd
              enddo
              if (nxmin<nbands .and. (dreal(zek(ibnd))-EF < emin)) nxmin=nxmin+1
              !
              WRITE(fh_ee+0,'(F19.14,2x,2I5,1x)',advance='no') omega(iom), nxmax-nxmin+1, nxmin+nemin-1
              WRITE(fh_ee+1,'(F19.14,2x,I5,1x)') omega(iom), nxmax-nxmin+1
              WRITE(fh_ee+2,'(F19.14,2x,I5,1x)') omega(iom), nxmax-nxmin+1
              do ibnd=nxmin,nxmax
                 WRITE(fh_ee+0,'(2E24.16,1x)',advance='no') zek(ibnd)!-EF
                 WRITE(fh_ee+1,'(I4,1x)',advance='no') ibnd+nemin-1
                 WRITE(fh_ee+2,'(I4,1x)',advance='no') ibnd+nemin-1
                 DO iorb1=1,norbitals
                    nind1 = nindo(iorb1)
                    do ind1=1,nind1
                       WRITE(fh_ee+1,'(F12.5,1x,F12.5,3x)',advance='no') UAl(ibnd,ind1,iorb1)
                       WRITE(fh_ee+2,'(F12.5,1x,F12.5,3x)',advance='no') UAr(ind1,ibnd,iorb1)
                    enddo
                 ENDDO
                 WRITE(fh_ee+1,*)
                 WRITE(fh_ee+2,*)
              enddo
              WRITE(fh_ee+0,*)
           enddo
        endif

        !DO iorb1=1,norbitals
        !   nind1 = nindo(iorb1)
        !   DO iorb2=1,norbitals
        !      nind2 = nindo(iorb2)
        !      ! kvec(:3,ikp)/REAL(kvec(4,ikp))
        !      phase = dot_product(kp,POSC(:,iorb1)-POSC(:,iorb2))
        !      ex = exp(-2*pi*IMAG*phase)
        !      print *, iorb1, iorb2, 'phase=', phase
        !      ww = (wgh(ikp)/(nsymop*tweight))*ex
        !      !call zgemm('C','N', nind1, nind2, nbands, ww, DMFTU(:,:,iorb1), nbands, DMFTU(:,:,iorb2), nbands, (1.d0,0.d0), olp(:,:,iorb1,iorb2),maxdim2)
        !      
        !      do ind1=1,nind1
        !         csum=0
        !         do ibnd=1,nbands
        !            if (abs(DMFTU(ibnd,ind1,iorb1)).gt.1e-3 .and. abs(DMFTU(ibnd,ind1,iorb2)).gt.1e-3 ) then
        !               csum = csum + DMFTU(ibnd,ind1,iorb1)/DMFTU(ibnd,ind1,iorb2)*ww
        !               if (iorb1.ne.iorb2) print *, iorb1, iorb2, ind1, ibnd, dreal(DMFTU(ibnd,ind1,iorb1)/DMFTU(ibnd,ind1,iorb2)*ex)
        !            endif
        !         enddo
        !         olp(ind1,ind1,iorb1,iorb2) = olp(ind1,ind1,iorb1,iorb2) + csum
        !      enddo
        !      
        !   ENDDO
        !ENDDO

     ENDDO
     DEALLOCATE( DMFTU )
     DEALLOCATE( gij )
     DEALLOCATE( STrans )
     if (mode.EQ.'e') then
        DEALLOCATE( Al, Ar, zek)
        deallocate(  UAl, UAr )
     endif
  ENDDO
  close(fh_ene)
  close(fh_p)


  if (mode.EQ.'g') then
     call Reduce_MPI(gmloc, norbitals, nom, maxdim, ncix)

     do icix=1,ncix
        close(fh_gk+icix)
     enddo
     deallocate( gmk )

     ! Writing out gloc to fh_gc+icix and its vector form to fh_gl+icix
     if (myrank.eq.master) then
        do icix=1,ncix
           cixdm = cixdim(icix)
           allocate( Sigini(cixdm,cixdm), cini(cixdm) )
           CALL GetSiginiCini(Sigini, cini, cixdms, Sigind(:,:,icix), cixdm, maxdim )
           !
           open(fh_gc+icix,FILE=fglocout(icix),STATUS='unknown')
           do iom=1,nom
              WRITE(fh_gc+icix,'(f14.8,2x)',advance='no') omega(iom)
              DO ip=1,cixdms
                 do iq=1,cixdms
                    WRITE(fh_gc+icix, '(f14.8,1x,f14.8,2x)',advance='no') gmloc(cini(ip),cini(iq),icix,iom)
                 enddo
              ENDDO
              WRITE(fh_gc+icix,*)
           enddo
           close(fh_gc+icix)
           deallocate( Sigini, cini )
        enddo

        ALLOCATE( Glc(maxsize, ncix, nom) )
        CALL GtoVectorForm(Glc, gmloc, Sigind, csize, cixdim, ncix, maxsize, maxdim, nom )

        do icix=1,ncix
           open(fh_gl+icix,FILE=fgloc(icix),STATUS='unknown')
        enddo
        CALL PrintGloc(fh_gl, Glc, omega, csize, ncix, nom, maxsize)
        DEALLOCATE( Glc )
        do icix=1,ncix
           close(fh_gl+icix)
        enddo
     endif
     deallocate( gmloc )

     call Barrier() ! some processors might not close the file above, and copying would not have all the data.
     ! Here we combine all g_k(omega) into single file per cix.
     ! Before each processor printed into its own file. But now we need to create one output file.
     if (myrank.eq.master) then
        DO icix=1,ncix
           open(fh_gk+icix, file=TRIM(fgkout(icix)), status='old', position='append', action='write', iostat=ierr)
           if (ierr /= 0) then
              print *, 'Error opening ', TRIM(fgkout(icix)), ' in append mode.'
              stop 'ERROR in dmftgk'
           endif
           do slave=1,nprocs-1
              write(cpuID_,'(I3)') slave
              filename = TRIM(fgkout(icix))//'.'//ADJUSTL(TRIM(cpuID_))
              open(fh_ene+icix, file=filename, status='old',iostat=ierr)
              if (ierr /= 0) then
                 print *, 'Error opening ', filename, ' in read mode.'
                 stop 'ERROR in dmftgk'
              endif
              do
                 read(fh_ene+icix, '(A)', iostat=ierr) single_line
                 if (ierr /= 0) exit  ! end of file or error
                 !if (single_line(1:1).eq.'!') then
                 !write(fh_gk+icix, '(A)') trim(single_line)
                 !else
                 write(fh_gk+icix, '(A)') trim(single_line)//' ' ! It turns out fortran adds a space when printing numbers
                 !endif
              end do
              close(fh_ene+icix)
              call execute_command_line('rm -f '//trim(filename))
           enddo
           close(fh_gk+icix)
        ENDDO
     endif
  else if (mode.EQ.'e') then
     close( fh_ee )
     close( fh_ee+1 )
     close( fh_ee+2 )
  endif

  !DO iorb1=1,norbitals
  !   nind1 = nindo(iorb1)
  !   DO iorb2=1,norbitals
  !      nind2 = nindo(iorb2)
  !      WRITE(123,'(I4,1x,I4)') iorb1,iorb2
  !      do i=1,nind1
  !         do j=1,nind2
  !            WRITE(123,'(f14.8,1x,f14.8,2x)',advance='no') olp(i,j,iorb1,iorb2)
  !         enddo
  !         WRITE(123,*)
  !      enddo
  !      WRITE(123,*)
  !   ENDDO
  !   WRITE(123,*)
  !ENDDO
  !deallocate( olp )

  DEALLOCATE( sigma, omega )


  if (mode.eq.'e' .and. nprocs>1) then
     ! Here we combine all eigenvalues.dat, UL.dat, UR.dat together.
     ! Before each processor printed into its own file. But now we need to create one output file.
     call Barrier() ! some processors might not close the file above, and copying would not have all the data.
     ! Now combine all files into one output file
     fnames(1) = feigvals
     fnames(2) = fUL
     fnames(3) = fUR
     if (myrank.eq.master) then
        do ip=1,3
           open(fh_ee+ip, file=TRIM(fnames(ip)), status='old', position='append', action='write', iostat=ierr)
           if (ierr /= 0) then
              print *, 'Error opening ', TRIM(fnames(ip)), ' in append mode.'
              stop 'ERROR in dmftgk'
           endif
           do slave=1,nprocs-1
              write(cpuID_,'(I3)') slave
              filename = TRIM(fnames(ip))//'.'//ADJUSTL(TRIM(cpuID_))
              open(fh_ene+ip, file=filename, status='old',iostat=ierr)
              if (ierr /= 0) then
                 print *, 'Error opening ', filename, ' in read mode.'
                 stop 'ERROR in dmftgk'
              endif
              do
                 read(fh_ene+ip, '(A)', iostat=ierr) single_line
                 if (ierr /= 0) exit  ! end of file or error
                 if (single_line(1:1).eq.'!') then
                    write(fh_ee+ip, '(A)') trim(single_line)
                 else
                    write(fh_ee+ip, '(A)') trim(single_line)//' ' ! It turns out fortran adds a space when printing numbers
                 endif
              end do
              close(fh_ene+ip)
              call execute_command_line('rm -f '//trim(filename))
           enddo
           close(fh_ee+ip)
        enddo
     endif
  endif



  return
9010 FORMAT(/,2X,' KP:',I6,' NEMIN NEMAX : ',2I5, ' dE:',2f5.2,' K:',a10 /)
9040 FORMAT(3X,2I4,6E13.6,F13.8)
903 print *, 'ERROR opening ', TRIM(fsigname(icix)), ' file'
904 print *, 'ERROR opening ', TRIM(filename), ' file'
905 print *, 'ERROR opening ', TRIM(fenergy), ' file'
END SUBROUTINE PrintGk


SUBROUTINE ReadKlist(fklist, nkpt, kvec, wgh)
  IMPLICIT NONE
  CHARACTER*100, intent(in) :: fklist
  INTEGER, intent(in)       :: nkpt
  REAL*8, intent(out)       :: wgh(nkpt)
  INTEGER, intent(out)      :: kvec(4,nkpt)
  ! locals
  logical      :: newform
  INTEGER      :: numkpt, ik, ios
  CHARACTER*10 :: KNAME
  CHARACTER*161:: line
  open(14,file=fklist,status='old',form='formatted')
  newform = .TRUE.
  read(14, '(A20)', IOSTAT=ios) line
  IF (line(15:16) .NE. ' ') THEN
     newform=.FALSE.
  ENDIF
  close(14)
  ! reopen
  open(14,file=fklist,status='old',form='formatted')
  DO ik=1,nkpt
     !READ (14,'(A10,4I5,3F5.2,A3)',IOSTAT=ios) KNAME, kvec(1,ik), kvec(2,ik), kvec(3,ik), kvec(4,ik), wgh(ik)
     IF(newform) THEN
        READ (14,'(A10,4I10,F5.2)',IOSTAT=ios) KNAME, kvec(1,ik), kvec(2,ik), kvec(3,ik), kvec(4,ik), wgh(ik) 
     ELSE
        READ (14,'(A10,4I5,F5.2)',IOSTAT=ios) KNAME, kvec(1,ik), kvec(2,ik), kvec(3,ik), kvec(4,ik), wgh(ik)
     ENDIF
     !print *, 'KNAME=', KNAME, 'kv=', kvec(1,ik), kvec(2,ik), kvec(3,ik), kvec(4,ik), wgh(ik)
     IF (KNAME .EQ. 'END       ' .OR. ios.ne.0) then
        print *, 'ERROR in reading ', TRIM(fklist)
        EXIT
     ENDIF
  ENDDO
  close(14)
END SUBROUTINE ReadKlist


SUBROUTINE ReadBR(frotlm, BRX)
  IMPLICIT NONE
  CHARACTER*100, intent(in):: frotlm
  REAL*8, intent(out)      :: BRX(3,3)
  ! locals
  INTEGER :: i
  REAL*8  :: BR1(3,3), BR2(3,3)
  !
  open(22,file=frotlm,status='old', form='formatted')
  READ(22,*) ! BR1
  DO i=1,3
     READ(22, '(3F10.5)') BR1(i,1), BR1(i,2), BR1(i,3)
  ENDDO
  READ(22,*) ! BR2
  DO i=1, 3
     READ(22, '(3F10.5)') BR2(i,1), BR2(i,2), BR2(i,3)
  ENDDO
  close(22)

  CALL dinv(BR1,3)
  BRX = matmul(BR2,BR1)
  
END SUBROUTINE ReadBR



CHARACTER*100 FUNCTION ReadWord(unit)
  IMPLICIT NONE
  INTEGER, intent(in) :: unit
  CHARACTER      :: buffer
  CHARACTER*100  :: word
  INTEGER :: size, i
  
  !print *, 'Starting ReadWord'
  
  DO i=1,100 ! Reads spaces
     READ(unit, "(A1)", ADVANCE='NO', SIZE=size, EOR=20, END=20) buffer
     if (buffer.NE.' ') EXIT
  ENDDO
  word=buffer
  DO i=1,100
     READ(unit, "(A1)", ADVANCE='NO', SIZE=size, EOR=20, END=20) buffer
     if (buffer.EQ.' ') EXIT
     word = TRIM(word)//buffer
     !print *, i, buffer, word
  ENDDO
20 CONTINUE
  ReadWord = word
  RETURN 
END FUNCTION ReadWord

SUBROUTINE ReadEnergiesK(fh_ene, Ek, kp, wg, NE, nume)
  IMPLICIT NONE
!!!! Reading band energies from energy file
  !REAL*8, intent(out) :: wgh(nkpt), Ek(nume)
  INTEGER, intent(in) :: fh_ene
  REAL*8, intent(out) :: kp(3), wg, Ek(nume)
  INTEGER, intent(out):: NE
  INTEGER, intent(in) :: nume
  ! locals
  INTEGER :: is, ios, itape, ii, NUM
  REAL*8  :: S, T, Z, E1
  CHARACTER*10 :: KNAME
  INTEGER :: N
  Ek(:)=0.0
  READ(fh_ene,'(3e19.12,a10,2i6,e19.12)',IOSTAT=ios) kp(1),kp(2),kp(3),KNAME,N,NE,WG
  DO ii=1,NE
     READ(fh_ene,*) NUM, E1
     Ek(NUM)=E1
  ENDDO
END SUBROUTINE ReadEnergiesK
  

SUBROUTINE ReadDMFT_TransK_Outside(nbands, nemin, maxdim2, norbitals, wkii, fh_p, pform, ikp)
  USE com_mpi, ONLY: cpuID, stop_MPI, fUdmft
!!! Start reading DMFT transformation
  IMPLICIT NONE
  INTEGER, intent(out)   :: nbands, nemin, maxdim2, norbitals
  INTEGER, intent(inout) :: wkii
  INTEGER, intent(in)    :: fh_p, ikp
  LOGICAL, intent(in)    :: pform
  ! locals
  CHARACTER*10 :: skii
  CHARACTER*100:: filename
  INTEGER :: iikp
  if (pform) then
     READ(fh_p,*) iikp, nbands, maxdim2, norbitals, nemin
  else
     READ(fh_p) iikp, nbands, maxdim2, norbitals, nemin
  endif
  if (iikp.NE.ikp) then
     filename = TRIM(fUdmft)//ADJUSTL(TRIM(cpuID))
     WRITE(6,'(A,I3,A,A,A,I3)') 'ERROR ikp=', ikp, 'but ikp from ',filename,' is', iikp
     call flush(6)
     call stop_MPI()
     STOP 'ERROR dmftgk'
  endif
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
           READ(fh_p,*) (DMFTU(i,ind,iorb),i=1,nbands)
        else
           READ(fh_p) (DMFTU(i,ind,iorb),i=1,nbands)
        endif
     ENDDO
  ENDDO
END SUBROUTINE ReadDMFT_TransK_Inside

SUBROUTINE Read_The_Rest_Basic_Arrays(fhb, nindo, cixdim, nl, ll, iorbital, cix, nind, csize, iSx, Sigind, EF, VOL, posc, norbitals, ncix, natom, maxdim, maxdim2)
  IMPLICIT NONE
  INTEGER, intent(in)  :: fhb
  INTEGER, intent(out) :: nindo(norbitals), cixdim(ncix), nl(natom)
  INTEGER, intent(out) :: ll(natom,4), iorbital(natom,4), cix(natom,4), nind(natom,4)
  INTEGER, intent(out) :: csize(ncix), iSx(maxdim2,norbitals), Sigind(maxdim,maxdim,ncix)
  REAL*8, intent(out)  :: EF, VOL
  REAL*8, intent(out)  :: posc(3,norbitals)
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

  DO iorb=1,norbitals
     READ(fhb,*) (posc(ip,iorb),ip=1,3)
  enddo
END SUBROUTINE Read_The_Rest_Basic_Arrays




SUBROUTINE CmpGkc2(gmk, gij, DMFTU, iSx, iorbital, ll, nl, cix, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: gmk(maxdim,maxdim,ncix)
  COMPLEX*16, intent(in) :: gij(nbands,nbands)
  COMPLEX*16, intent(in) :: DMFTU(nbands,maxdim2,norbitals)
  INTEGER, intent(in)     :: iSx(maxdim2, norbitals)
  INTEGER, intent(in)     :: iorbital(natom,lmax2+1), ll(natom,4), nl(natom), cix(natom,4)
  INTEGER, intent(in)     :: natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize
  ! locals
  COMPLEX*16, allocatable :: tgk(:,:), tmp(:,:)
  INTEGER    :: it, i, j, icase, jcase, l1case, l2case, icix, l1, l2, nind1, nind2, iorb1, iorb2, ind1, ind2
  gmk=0
  allocate( tmp(maxdim2,nbands), tgk(maxdim2,maxdim2) )
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
              !print *, 'iorb1,iorb2=', iorb1, iorb2, maxdim2, nind1, nind2, nbands, maxdim2
              
              call zgemm('C','N', nind1, nbands, nbands, (1.d0,0.d0), DMFTU(:,:,iorb1), nbands, gij(:,:),nbands, (0.d0,0.d0), tmp(:,:),maxdim2)
              call zgemm('N','N', nind1, nind2, nbands, (1.d0,0.d0), tmp,maxdim2, DMFTU(:,:,iorb2),nbands, (0.d0,0.d0), tgk,maxdim2)
              
              do ind1=1,nind1
                 do ind2=1,nind2
                    gmk( iSx(ind1,iorb1), iSx(ind2,iorb2), icix ) = tgk(ind1,ind2)
                 enddo
              enddo
           enddo
        ENDDO
     enddo
  ENDDO

  deallocate( tmp, tgk )
END SUBROUTINE CmpGkc2

!SUBROUTINE FindEigsys(zek, Al, Ar, ix, wsigma, gij, womega, omega, STrans, Ek, nom, ncix, nbands, nemin, maxsize, nume, csize, sigma)
!  IMPLICIT NONE
!  COMPLEX*16, intent(out)   :: zek(nbands), Al(nbands,nbands), Ar(nbands, nbands)
!  INTEGER, intent(inout)    :: ix
!  COMPLEX*16, intent(inout) :: wsigma(maxsize, ncix)
!  COMPLEX*16, intent(inout) :: gij(nbands,nbands)
!  REAL*8, intent(in)     :: womega
!  REAL*8, intent(in)     :: omega(nom)
!  COMPLEX*16, intent(in) :: STrans(maxsize,ncix,nbands,nbands)
!  REAL*8, intent(in)     :: Ek(nume)
!  INTEGER, intent(in)    :: nom, ncix, nbands, nemin, maxsize, nume, csize(ncix)
!  COMPLEX*16, intent(in) :: sigma(nom, maxsize, ncix)
!  ! locals
!  INTEGER :: icix, ind, i
!  
!  gij=0
!  DO i=1,nbands
!     gij(i,i) = Ek(i+nemin-1)  !-----  g^-1 of the LDA part in band representation ----!
!  ENDDO
!  CALL findNext(womega, ix, omega, nom)              
!  do icix=1,ncix
!     do ind=1,csize(icix)
!        CALL interp(wsigma(ind,icix), sigma(:,ind,icix), womega, ix, omega, nom)
!     enddo
!  enddo
!  CALL AddSigma_optimized2(gij, wsigma, STrans, csize, 1, nbands, ncix, maxsize)           
!  CALL eigsys(gij, zek, Al, Ar, nbands)
!END SUBROUTINE FindEigsys

SUBROUTINE PrintGloc(fh_gl, Glc, omega, csize, ncix, nomega, maxsize)
  !---Currently we use          :  fh_dos = 500, fh_gl = 180, fh_dl = 280
  IMPLICIT NONE
  INTEGER, intent(in)    :: fh_gl
  COMPLEX*16, intent(in) :: Glc(maxsize,ncix,nomega)
  REAL*8, intent(in)     :: omega(nomega)
  INTEGER, intent(in)    :: csize(ncix)
  INTEGER, intent(in) :: nomega, ncix, maxsize
  ! local
  INTEGER :: L, wndim, iom, lcase, i, j, icix, itape, jtape, icase, iorb, wmaxsize, it
  COMPLEX*16 :: csum
  REAL*8     :: pi

  pi=ACOS(-1.0D0)
  
  do iom=1,nomega
     do icix=1,ncix
        ! Header
        itape = fh_gl+icix
        write(itape,'(f14.8,1x)',advance='no') omega(iom)
        do i=1,csize(icix)
           write(itape,'(2f14.6)',advance='no') Glc(i,icix,iom)
        enddo
        write(itape,*)
     enddo
  enddo
  return
END SUBROUTINE PrintGloc



SUBROUTINE GtoVectorForm(Glc, gmloc, Sigind, csize, cixdim, ncix, maxsize, maxdim, nomega )
  IMPLICIT NONE
  COMPLEX*16, intent(out):: Glc(maxsize,ncix,nomega)
  COMPLEX*16, intent(in) :: gmloc(maxdim,maxdim,ncix,nomega)
  INTEGER, intent(in)    :: Sigind(maxdim,maxdim,ncix), csize(ncix), cixdim(ncix)
  INTEGER, intent(in)    :: ncix, maxsize, maxdim, nomega
  !----- locals
  INTEGER    :: icix, ip, iq, it, iom, cixdm, cixdms
  REAL*8 :: PI, beta, omw, ywr, ywi, ypwr, ypwi, yppwr, yppwi
  !REAL*8, allocatable  :: x(:), yr(:,:), yi(:,:), yppr(:,:), yppi(:,:)
  INTEGER, ALLOCATABLE :: cind(:), cini(:), Sigini(:,:)
  INTEGER :: noccur(maxsize,ncix)

  noccur=0
  DO icix=1,ncix
     DO ip=1,cixdim(icix)
        DO iq=1,cixdim(icix)
           it = Sigind(ip,iq,icix)
           if (it.gt.0)  noccur(it,icix) = noccur(it,icix) + 1
        ENDDO
     ENDDO
  ENDDO
  
  Glc=0
  DO icix=1,ncix

     cixdm = cixdim(icix)
     allocate( Sigini(cixdm,cixdm), cini(cixdm) )
     CALL GetSiginiCini(Sigini, cini, cixdms, Sigind(:,:,icix), cixdm, maxdim )


     !cixdm = cixdim(icix)
     !allocate( cind(cixdm) )
     !! If Sigind(i,i)=0, we eliminate the i-th column and rown, because such orbital should be treated as non-correlated
     !cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated
     !cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block
     !DO ip=1,cixdm
     !   it = Sigind(ip,ip,icix)
     !   if (it.gt.0) then
     !      cixdms = cixdms + 1
     !      cind(ip) = cixdms
     !   endif
     !ENDDO
     !allocate( cini(cixdms), Sigini(cixdms,cixdms) )
     !do ip=1,cixdm
     !   if (cind(ip).gt.0) cini(cind(ip))=ip
     !enddo
     !DO ip=1,cixdms
     !   do iq=1,cixdms
     !      Sigini(ip,iq) = Sigind(cini(ip),cini(iq),icix)
     !   enddo
     !ENDDO
     !deallocate( cind )
     !print *, 'cini=', cini
     !print *, 'Sigini=', Sigini
     do iom=1,nomega !------- for all frequencies --!
        !---- packing delta to vector form
        DO ip=1,cixdms
           do iq=1,cixdms
              it = Sigini(ip,iq)
              if (it.gt.0) then
                 Glc( it, icix, iom ) =  Glc( it, icix, iom ) + gmloc(cini(ip),cini(iq),icix,iom)
              endif
           enddo
        ENDDO
        DO it=1,csize(icix)
           Glc(it, icix, iom) = Glc(it, icix, iom)/noccur(it,icix)
        ENDDO
     enddo
     deallocate( cini, Sigini )
  ENDDO
END SUBROUTINE GtoVectorForm



SUBROUTINE GetSiginiCini(Sigini, cini, cixdms, Sigind, cixdim, maxdim )
  IMPLICIT NONE
  INTEGER, intent(in)    :: Sigind(maxdim,maxdim), cixdim
  INTEGER, intent(in)    :: maxdim
  INTEGER, intent(out)   :: cini(cixdim), Sigini(cixdim,cixdim), cixdms
  !----- locals
  INTEGER :: cind(cixdim), ip, iq, it

  ! If Sigind(i,i)=0, we eliminate the i-th column and rown, because such orbital should be treated as non-correlated
  cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated
  cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block
  DO ip=1,cixdim
     it = Sigind(ip,ip)
     if (it.gt.0) then
        cixdms = cixdms + 1
        cind(ip) = cixdms
     endif
  ENDDO

  !allocate( cini(cixdms), Sigini(cixdms,cixdms) )
  do ip=1,cixdim
     if (cind(ip).gt.0) cini(cind(ip))=ip
  enddo
  DO ip=1,cixdms
     do iq=1,cixdms
        Sigini(ip,iq) = Sigind(cini(ip),cini(iq))
     enddo
  ENDDO
END SUBROUTINE GetSiginiCini
