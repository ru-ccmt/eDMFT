! @Copyright 2007 Kristjan Haule
! 

PROGRAM DMFTMAIN  ! Program DMFT calculates 
  USE param
  USE structure, ONLY: ReadStructure, WriteInfoStructure, DeallocateStructure, mult, ndif, iz, tau, rotloc, pos, rotij, tauij, BR1, BR2, LATGEN, ROTDEF, iord, nat, rot_spin_quantization
  USE case,      ONLY: nl, ll, qsplit, cix, iatom, crotloc, isort, ifirst, shft, CF, Sigind, legend, csize, maxdim, maxsize, natom, ncix
  USE sym2,      ONLY: idet, phase, init_sym2
  USE com,       ONLY: EF, emin, emax, iso, iso2, nspin1
  USE kpts,      ONLY: numkpt, MKNAME, MIPGR, mweight, tweight, MSX, MSY, MSZ, allocate_kpts, deallocate_kpts
  USE com_mpi,   ONLY: VECFN, vectors, fvectors, myrank, master, nprocs, nvector, vector_para, filenamempi, filenamempi2, start_mpi, stop_mpi, Scatter_Vector_data, fastFilesystem, Qprint
  use abc,       ONLY: init_abc, deallocate_abc
  use error,     ONLY: errflg, errclr, outerr
  IMPLICIT NONE
  complex*16   :: imag
  real*8       :: xx(3), xz(3), check, ddd, E1, E2, ident, phix, phiz, thetax, thetaz, ttime, time1, time0
  integer      :: i, j, j1, jc, jr, jatom, lateq, icase,  ISX, ISY, ISZ, IDV, index, index1, info, ios, irecl, itape, iunit, loro, m1, m2
  CHARACTER*4  :: adum
  CHARACTER*5  :: CHAR
  CHARACTER*10 :: KNAME
  CHARACTER*11 :: STATUS,FORM                                      
  CHARACTER*200:: ERRMSG
  CHARACTER*80 :: DEFFN, ERRFN
  CHARACTER*200:: FNAME, fUdmft, vec_up, vec_dn, ene_up, ene_dn
  CHARACTER*161:: TEXT
  LOGICAL      :: SO, Qcomplex, newform, Qident, Qrenormalize
  CHARACTER*1  :: mode
  INTEGER      :: kindex, nsymop
  INTEGER      :: wndim, icix, wicix, size, imatsubara, irenormalize, icmp_partial_dos
  REAL*8       :: EF_LDA, Ry2eV, cost, cosf, sint, sinf, theta, fi, wdet
  REAL*8       :: POST(3), BKRLOC(3,3), tmp(3,3), BR1inv(3,3), S(3)
  INTEGER      :: projector
  INTEGER      :: strlen, locrot, shift, latom, wfirst, iat, im, wat, minsigind_p, maxsigind_p, minsigind_m, lngth
  INTEGER      :: ivector, nkpt, inkp, ik_start, il
  INTEGER      :: identity3(3,3), diff, ig
  interface
     REAL*8 FUNCTION detx(a)
       IMPLICIT NONE
       REAL*8, intent(in) :: a(3,3)
     end FUNCTION detx
  end interface

  REAL*8,ALLOCATABLE:: wtmp(:)
  DATA SO /.false./, Qcomplex /.false./, IMAG/(0.0D0,1.0D0)/, Ry2eV/13.60569193/
  DATA strlen/200/
  !------------------------------------------------------------------     

  fUdmft = 'Udmft.'

  CALL start_MPI()
  
  call cputim(time0)
  
  CALL GTFNAM(DEFFN,ERRFN)
  
  if (myrank.EQ.master) CALL ERRFLG(ERRFN,'Error in dmft1')
  
  nspin1=1
  !WRITE(*,*) myrank, ('Definition file is ' // DEFFN)
  OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910,IOSTAT=ios) ! opens dmftx.def file
  
  vector_para = .False.
  iso=1
  iso2=1
  ! The following few lines read dmft1.def file containing all filename definitions
  !
  DO 
     READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL  ! can jump to 20
     lngth = Len_Trim(FNAME)
     if ( FNAME(lngth-1:lngth).EQ.'_x' ) then
        FNAME = FNAME(:lngth-2)
        if(iunit.eq.9 .or. iunit.eq.10) vector_para = .True.
     else if (iunit.eq.59 .or.iunit.eq.60) then  ! should this be 59 & 60?
        !
     else if (iunit.EQ.6) then
        if (myrank.eq.master) then
           print *, 'Opening', FNAME
           OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        else if (fastFilesystem) then
           CALL FilenameMPI(FNAME)    ! Each processor will have different information file
           OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        endif
        cycle
     else if (iunit.eq.21) then
        if (myrank.eq.master) then
           OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        endif
        cycle
     else if ( iunit.lt.180  .or. myrank.eq.master) then
        ! output files (iunit>=100) should be open on master node only
        OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
     endif
     if(iunit.eq.9)  VECFN(1)=FNAME
     if(iunit.eq.10) VECFN(2)=FNAME
     if(iunit.eq.59) VECFN(3)=FNAME
     if(iunit.eq.60) VECFN(4)=FNAME
     if(iunit.eq.18) then
        do i=strlen,5,-1
           if(fname(i:i).ne.' ') then
              if(fname(i-2:i).eq.'pup') nspin1=2
              if(fname(i-2:i).eq.'pdn') nspin1=2  ! added 10.07.08 PN
           endif
        enddo
     endif
     if(iunit.eq.7) then
        read(iunit,123,end=12)CHAR
        Qcomplex=.true. ! switches Qcomplex  on if case.in1c exists and is nonzero
12      continue
     end if
     if(iunit.eq.501) then
        do i=strlen,2,-1
           if(fname(i:i).ne.' ') then
              fUdmft = fname(:i-1)
              exit
           endif
        enddo
        if (myrank.eq.master) close(501) ! Will reopen later
     endif
     if(iunit.eq.9) then
        iloop1: do i=strlen,2,-1
           if(fname(i:i).ne.' ') then
              if(fname(i-1:i).eq.'so') then
                 SO=.true.
                 iso=2
                 iso2=1
                 Qcomplex=.true.
                 exit iloop1
              endif
           endif
        enddo iloop1
        iloop2: do i=strlen,4,-1
           if(fname(i:i).ne.' ') then
              if(fname(i-3:i).eq.'soup') then
                 SO=.true.
                 iso=2
                 iso2=2
                 Qcomplex=.true.
                 exit iloop2
              endif
           endif
        enddo iloop2
     endif
  ENDDO ! 10
20 CONTINUE
  CLOSE (1)
  ! finished dmft1.def

  !.....READING STRUCT 
  ! It reads case.struct file and saves data into module named "param". These variables are available in dmftmain.f.
  !CALL init_struct
  call ReadStructure(20)
  if (Qprint) then
     write(6,'(A)')'******Program DMFT: Computes Gloc and Delta for Dynamical Mean Field Theory calculation ******'
     if (Qcomplex) write(6,*)'CALCULATION WITH COMPLEX VECTORS'
     if (so)then
        if(iso2.eq.2)write(6,*)'CALCULATION WITH SOC-UP AND DN VECTORS REQUIRED'
        if(iso2.eq.1)write(6,*)'CALCULATION WITH SOC VECTORS REQUIRED'
     endif
  endif
  
  if (nprocs.EQ.1) then
     write(6,*)'Running DMFT in single processor mode'
     !write(6,*)' '
  else
     if (Qprint) write(6,*)'Running DMFT in mpi mode'
     !write(6,*)' '
  endif

  
  ! find Fermi energy using scf2 file
  EF = 0.0
  EF_LDA =0.0
  DO 
     read(8,666,end=14) adum
     if(adum.eq.':NOE') then
        read(8,*)
        read(8,667) EF_LDA
        exit
     endif
  ENDDO
  EF = EF_LDA
14 continue  ! Fermi energy end


  ! read klist and associated weights from file 66==case.klist

  ! first read how many k-points exists and the old/new file type
  ITAPE = 14
  numkpt = 0
  newform = .TRUE.
  DO
     READ (ITAPE, '(A20)', IOSTAT=ios) TEXT
     IF (numkpt==0 .AND. TEXT(15:16) .NE. ' ') THEN
        newform=.FALSE.
     ENDIF
     KNAME = TEXT(1:10)
     IF (KNAME .EQ. 'END       ' .OR. ios.ne.0) EXIT
     numkpt = numkpt + 1
  ENDDO
  REWIND(ITAPE)

  if (Qprint) WRITE(6,*) 'numkpt=', numkpt
  
  CALL allocate_kpts(numkpt)
  MWEIGHT=0
  DO kindex=1,numkpt
     IF(newform) THEN
        READ (ITAPE,5101) MKNAME(KINDEX), ISX, ISY, ISZ, IDV, MWEIGHT(KINDEX), E1, E2, MIPGR(KINDEX)
     ELSE
        READ (ITAPE,5100) MKNAME(KINDEX), ISX, ISY, ISZ, IDV, MWEIGHT(KINDEX), E1, E2, MIPGR(KINDEX)
     ENDIF
     MSX(KINDEX) = DBLE(ISX)/DBLE(IDV)
     MSY(KINDEX) = DBLE(ISY)/DBLE(IDV)
     MSZ(KINDEX) = DBLE(ISZ)/DBLE(IDV)
     if (Qprint) WRITE(6,'(4f10.5,1x)') MSX(kindex), MSY(kindex), MSZ(kindex), MWEIGHT(kindex)
  END DO
  TWEIGHT = sum(MWEIGHT)
  !print *, 'twgh=', myrank, tweight

  ! iord : number of symmetry operations from case.struct
  CALL init_sym2(iord) ! allocates memory for: idet(iord),iz_cartesian(3,3,iord),phase(iord),tmat(3,3,iord)

  identity3(:,:)=0
  identity3(1,1)=1
  identity3(2,2)=1
  identity3(3,3)=1
  do ig=1,iord
     diff = sum(abs(iz(:,:,ig)-identity3(:,:)))
     if (diff.eq.0) exit
  enddo
  !print *, 'Found identity at ig=', ig
  if (ig.ne.1) then
     if (ig.eq.iord+1) then
        WRITE(6,*) 'ERROR : Could not find indenity among symmetry operations'
     endif
     CALL swapGroup(iz(:,:,ig),tau(:,ig),iz(:,:,1),tau(:,1))
  endif
  
  if (Qprint) then
     call WriteInfoStructure(6)
  endif

  if (vector_para) then
     FNAME = '_processes_'
     !WRITE(6,*) 'Trying to open ', FNAME
     open(999,FILE=FNAME,STATUS='OLD',ERR=77,IOSTAT=ios) ! opens _processes_, hence this parallel dmft_lapw1 run, and not parallel w2k_lapw1 run
     !WRITE(6,*) 'Found ', FNAME
     call Scatter_Vector_data()
     goto 88
77   CONTINUE
     CALL FilenameMPI2(FNAME)
     !WRITE(6,*) 'Trying old method with ', FNAME
     open(999,FILE=FNAME,STATUS='OLD',ERR=88,IOSTAT=ios) ! opens _processes_x                                                                                  
     nvector=0
     vectors=0
     DO 
        READ (999,*,END=88,ERR=970) ivector,inkp,ik_start,vec_up,vec_dn,ene_dn,ene_up  ! can jump to 20
        nvector = nvector+1
        vectors(nvector,1) = ivector    ! the successive number of vector file
        vectors(nvector,2) = inkp       ! number of k-points
        vectors(nvector,3) = ik_start   ! what is the first k-point in this vector file
        fvectors(nvector,1) = vec_up    ! filenames....                                 
        fvectors(nvector,2) = vec_dn
        fvectors(nvector,3) = ene_dn
        fvectors(nvector,4) = ene_up
        if (Qprint) WRITE(6,*) 'FNAMES=', fvectors(nvector, :)
     ENDDO
88   CONTINUE
     if (Qprint) print *, 'Vector files are read in parallel mode'
     close(999)
  else
     nvector=1
  endif
  
  ! LATGEN calculates BR1 and BR2
  CALL LATGEN(myrank.eq.master)
  ! ROTDEF calculates tauij, rotij for non-symorphic lattices
  CALL ROTDEF(myrank.EQ.master)

! Reads current input file case.indmf1
! 5 is connected to current input file case.indmfl
  !print *, 'Start reading 5'
  READ(5,*) EMIN,EMAX,irenormalize,projector
  if (abs(projector).lt.4) then
     EMIN = EMIN/Ry2eV
     EMAX = EMAX/Ry2eV
  else
     nemin0=int(EMIN)
     nemax0=int(EMAX)
  endif
  !print *, 'Emin,Emax=', Emin, Emax

  icmp_partial_dos = 1
  READ(5,*,iostat=ios) imatsubara, gammac, gamma, nom_default, aom_default, bom_default, icmp_partial_dos  ! We added a new flag, cmp_partial_dos
  if (ios.eq.0) then ! icmp_partial_dos is specified
     matsubara = .FALSE.
     if (imatsubara.EQ.1) matsubara = .TRUE.
     cmp_partial_dos = .FALSE.
     if (icmp_partial_dos.EQ.1) cmp_partial_dos = .TRUE.
  else
     BACKSPACE(5)
     READ(5,*, iostat=ios) imatsubara, gammac, gamma, nom_default, aom_default, bom_default                ! We might not specify that, to be compatible with other parts of the code
     if (ios .ne. 0) then
        print *, 'Wrong format for case.indmfl'
        call stop_MPI
        STOP 'ERROR dmft1'
     endif
     if (imatsubara.EQ.1) then
        matsubara = .TRUE.        
        cmp_partial_dos = .FALSE. ! By default, on imaginary axis we do not compute partial dos
     else
        matsubara = .FALSE.
        cmp_partial_dos = .TRUE.     ! By default, we only compute partial dos on the real axis
     endif
  endif
  
  read(5,*) natom
  if (irenormalize.EQ.0) then
      Qrenormalize = .FALSE.
  else
      Qrenormalize = .TRUE.
   endif
   
  gammac = gammac/Ry2eV !--- broadening for correlated orbitals to Ry --!
  gamma  = gamma/Ry2eV  !--- broadening for all orbitals to Ry ---------!
  

  if (Qprint) then
     write(6, 877) Emin,Emax  ! Energy window
     write(6, 887) natom      ! Projected density of states calculation for xxx atoms
     write(6,*) 'cmp_partial_dos=', cmp_partial_dos
  endif
  !write(21,877) Emin,Emax
  !write(21,887) natom
  
  ALLOCATE(nl(natom))
  allocate(ll(natom,4), qsplit(natom,4), cix(natom,4), iatom(natom))
  
  call FindMax_nl(natom, max_nl)
  ALLOCATE(crotloc(3,3,max_nl,natom))
  crotloc=0

  !natm = sum(mult) = ndif
  ALLOCATE( isort(ndif), ifirst(ndif) )
  wfirst  = 1          ! first atom of this particular sort
  do iat=1,nat         ! over all sorts
     do im=1,MULT(iat) ! over all atoms of this sort
        wat = wfirst + im-1  ! atom number
        isort(wat)=iat       ! sort of each atom
        ifirst(wat)=wfirst
     enddo
     wfirst = wfirst + MULT(iat)
  enddo

  
  ALLOCATE(shft(3,ndif))
  shft = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i=1,natom
     READ(5,*) latom,nL(i),locrot ! read from case.indmfl
     
     iatom(i) = latom       ! The succesive number of atom (all atoms counted)
     jatom = isort(latom)   ! The sort of the atom

     loro = MOD(locrot, 3)
     shift = locrot/3
     
     if (Qprint) WRITE(6,'(A,I2,A)') '--------- atom ', i, ' -------------'

     ! loro can be (0 -original coordinate system), (1 - new z-axis), (2 - new z-axis, x-axis)
     do j=1,nL(i)
        read(5,*) LL(i,j), qsplit(i,j), cix(i,j) ! LL contains all L-quantum numbers which need to be computed.
     enddo

     if (Qprint) then
        do j=1,nL(i)
           write(6, '(A,I3,2x,A,I3)') 'l=', LL(i,j), 'qsplit=', qsplit(i,j)!, 'nrfmax=', nrfmax(i,j)
        enddo
        write(6,*)'Symmetrization over eq. k-points is performed'
     endif

     !crotloc(:,:,i) = rotloc(:,:,jatom)  ! Rotation for this type from the struct file
     do j=1,nL(i)
        crotloc(:,:,j,i) = rotloc(:,:,jatom)  ! Rotation for this type from the struct file
     enddo
     
     ! Rotation matrix in case the coordinate system is changed!
     if (loro.eq.0) then
        ! nothing to do
     elseif (loro.gt.0)then        ! change of local rotation matrix
        print *, 'positive locrot not supported anymore. Please use locrot=-1,-2,-3,...'
        call stop_MPI
        STOP 'ERROR dmft1'
        !read(5,*)(xz(j),j=1,3)  ! new z-axis expressed in unit cell vectors
        !if (Qprint) write(6,120)(xz(j),j=1,3)
        !call angle(xz,thetaz,phiz)
        !crotloc(1,3,i)=sin(thetaz)*cos(phiz)
        !crotloc(2,3,i)=sin(thetaz)*sin(phiz)
        !crotloc(3,3,i)=cos(thetaz)          
        !if(loro.eq.1)then       ! only z-axis fixed
        !   crotloc(3,1,i)= 0.     ! new x perpendicular to new z       
        !   ddd=abs(crotloc(3,3,i)**2-1.)
        !   if(ddd.gt.0.000001)then
        !      crotloc(1,1,i)=-crotloc(2,3,i)
        !      crotloc(2,1,i)= crotloc(1,3,i)
        !   else
        !      crotloc(1,1,i)=1.
        !      crotloc(2,1,i)=0.
        !   endif
        !   ddd=0.                 ! normalize new x
        !   do j=1,3
        !      ddd=ddd+crotloc(j,1,i)**2
        !   enddo
        !   do j=1,3
        !      crotloc(j,1,i)=crotloc(j,1,i)/sqrt(ddd)
        !   enddo
        !elseif(loro.eq.2)then   ! also new x-axis fixed
        !   read(5,*)(xx(j),j=1,3)
        !   if (Qprint) write(6,121)(xx(j),j=1,3)
        !   !write(21,121)(xx(j),j=1,3)
        !   call angle(xx,thetax,phix)
        !   crotloc(1,1,i)=sin(thetax)*cos(phix)
        !   crotloc(2,1,i)=sin(thetax)*sin(phix)
        !   crotloc(3,1,i)=cos(thetax)          
        !   !  check orthogonality of new x and z axes
        !   check=0.
        !   do j=1,3
        !      check=check+crotloc(j,1,i)*crotloc(j,3,i)
        !   enddo
        !   if(abs(check).gt.0.00001)then
        !      write(6,*)' new x and z axes are not orthogonal'
        !      print *,' new x and z axes are not orthogonal'
        !      stop
        !   endif
        !endif
        !crotloc(1,2,i)=crotloc(2,3,i)*crotloc(3,1,i)-crotloc(3,3,i)*crotloc(2,1,i)
        !crotloc(2,2,i)=crotloc(3,3,i)*crotloc(1,1,i)-crotloc(1,3,i)*crotloc(3,1,i)
        !crotloc(3,2,i)=crotloc(1,3,i)*crotloc(2,1,i)-crotloc(2,3,i)*crotloc(1,1,i)
        !if (Qprint) then
        !   write(6,*)' New local rotation matrix in global orthogonal system'
        !   write(6,*)'                      new x     new y     new z'
        !   write(6,1013)((crotloc(j,j1,i),j1=1,3),j=1,3)  !written as in .struct
        !endif
     elseif (loro.eq.-1) then               ! in cartesian, and equal for all l cases
        read(5,*)( crotloc(1,j,1,i),j=1,3)  ! new x-axis expressed in unit cell vectors
        read(5,*)( crotloc(2,j,1,i),j=1,3)  ! new y-axis expressed in unit cell vectors
        read(5,*)( crotloc(3,j,1,i),j=1,3)  ! new z-axis expressed in unit cell vectors
        do il=2,nL(i)
           crotloc(:,:,il,i) = crotloc(:,:,1,i)
        enddo
        if (Qprint) then
           write(6,*)' New local rotation matrix in global orthogonal system'
           write(6,*)'                      new x     new y     new z'
           write(6,1013)((crotloc(j,j1,1,i),j1=1,3),j=1,3)  !written as in .struct
        endif
     elseif(loro.eq.-2) then
        do il=1,nL(i)
           read(5,*)( crotloc(1,j,il,i),j=1,3)  ! new x-axis expressed in unit cell vectors
           read(5,*)( crotloc(2,j,il,i),j=1,3)  ! new y-axis expressed in unit cell vectors
           read(5,*)( crotloc(3,j,il,i),j=1,3)  ! new z-axis expressed in unit cell vectors
        enddo
        if (Qprint) then
           do il=1,nL(i)
              write(6, '(A,I3,2x,A,I3)') 'T2C for l=', LL(i,il), 'qsplit=', qsplit(i,il)
              write(6,*)' New local rotation matrix in global orthogonal system'
              write(6,*)'                      new x     new y     new z'
              write(6,1013)((crotloc(j,j1,il,i),j1=1,3),j=1,3)  !written as in .struct
           enddo
        endif
     else
        print *, 'locrot=', locrot, 'with loro=', loro, 'and shift=',shift, 'not implemented'
        WRITE(ERRMSG,*) 'locrot=', locrot, 'with loro=', loro, 'and shift=',shift, 'not implemented'
        call OUTERR('dmft1', ERRMSG)
     endif

        
     ! end Rotation matrix loro

     if (shift.ne.0) then
        read(5,*)(shft(j,latom),j=1,3)
        if (Qprint) then
           write(6,'(A,1x,I2,A,2x,F6.1,F6.1,F6.1)') '** atom', i, ' has nonzero shift: ', shft(1,latom), shft(2,latom), shft(3,latom)
           write(6,*) shft(:,latom)
        endif
     endif
        
  enddo

  READ(2,*) ! comment: mode and the current chemical potential 
  READ(2,*) mode, EF, nsymop
  ! Check if we want to force Hermitian H
  Hrm = '  '
  READ(2,*,iostat=ios)  Hrm
  if (ios .eq. 0 .and. Hrm(1:1).eq.'H') then
     if (myrank.EQ.master) WRITE(6,*) 'Hermitian mode'
  else
     if (myrank.EQ.master) WRITE(6,*) 'Non hermitian mode'
  endif
  
  EF = EF/Ry2eV

  READ(5,*) ! comment: Next few lines contain instructions (transformation,index) for all correlated orbitals
  READ(5,*) ncix, maxdim, maxsize ! ncix -- number of cix blocks, maxdim -- largest dimension of any cix block, maxsize -- maximum number of independent components to store (can be between 1 and maxdim**2)
  if (Qprint) WRITE(6,*) '********** Start Reading Cix file *************'
  ALLOCATE(CF(maxdim,maxdim,ncix))      ! transformation matrix
  ALLOCATE(Sigind(maxdim,maxdim,ncix))  ! correlated index
  ALLOCATE(legend(maxsize, ncix))       ! names of correlated orbitals
  ALLOCATE(csize(ncix))                 ! independent components
  ALLOCATE(wtmp(2*maxdim))              ! temp
  
  CF=0
  Sigind=0  
  Qident=.True.
  do icix=1,ncix
     READ(5,*) wicix, wndim, size  ! icix, size-of-matrix, L, number of independent components of the matrix

     csize(icix)=size
     if (wicix.ne.icix) then
        print *, 'Something wrong reading case.indmf0 file. Boilig out...'
        goto 999
     endif
     READ(5,*) ! Comment: Independent components are
     READ(5,*,ERR=1999) (legend(i,icix),i=1,size)
1999  CONTINUE
     READ(5,*) ! Comment: Sigind follows
     do i=1,wndim
        READ(5,*) (Sigind(i,j,icix),j=1,wndim)
     enddo
     ! modify Sigind so that indices always begin with 1
     ! Example:
     ! 3 0 0 0 0        1 0 0 0 0
     ! 0 3 0 0 0        0 1 0 0 0
     ! 0 0 4 0 0  --->  0 0 2 0 0
     ! 0 0 0 4 0        0 0 0 2 0
     ! 0 0 0 0 4        0 0 0 0 2
     minsigind_p=10000
     maxsigind_p=0
     minsigind_m=10000
     do i=1,wndim
        do j=1,wndim
           if (Sigind(i,j,icix).gt.0) then
              minsigind_p = min(minsigind_p,Sigind(i,j,icix))
              maxsigind_p = max(maxsigind_p,Sigind(i,j,icix))
           endif
           if (Sigind(i,j,icix).lt.0) minsigind_m = min(minsigind_m,-Sigind(i,j,icix))
        enddo
     enddo
     maxsigind_p = maxsigind_p-minsigind_p+1
     if (maxsigind_p.lt.0) maxsigind_p = 0
     if (minsigind_p.ne.1 .or. minsigind_m.ne.10000) then
        do i=1,wndim
           do j=1,wndim
              if (Sigind(i,j,icix).gt.0) then
                 Sigind(i,j,icix) = Sigind(i,j,icix)-minsigind_p+1
              endif
              if (Sigind(i,j,icix).lt.0) then
                 Sigind(i,j,icix) = -(abs(Sigind(i,j,icix))-minsigind_m+maxsigind_p+1)
              endif
           enddo
        enddo
        
        if (Qprint) then
           WRITE(6,*) 'Sigind corrected to'
           do m1=1,wndim
              do m2=1,wndim
                 write(6,'(I3)',advance='no') Sigind(m1,m2,icix)
              enddo
              write(6,*)
           enddo
        endif
     endif

     READ(5,*) ! Comment: Transformation matrix follows
     do i=1,wndim
        READ(5,*) (wtmp(j),j=1,2*wndim)
        do j=1,wndim
           CF(i,j,icix) = dcmplx(wtmp(2*j-1),wtmp(2*j))
        enddo
     enddo

     do i=1,wndim
        if (.not.Qident) exit
        do j=1,wndim
           if (i.eq.j) then
              ident=1.d0
           else 
              ident=0.d0
           endif
           if (abs(CF(i,j,icix)-ident)>1e-5) then
              Qident=.False.
              exit
           endif
        enddo
     enddo

     if (Qprint) then
        write(6,*)' Correlated block number', icix
        write(6,*)' Correlated index='
        do m1=1,wndim
           do m2=1,wndim
              write(6,'(I3)',advance='no') Sigind(m1,m2,icix)
           enddo
           write(6,*)
        enddo
        write(6,*)' Real part of unitary matrix='
        do m1=1,wndim
           do m2=1,wndim
              write(6,'(F8.4)',advance='no') dble(cf(m1,m2,icix))
           enddo
           write(6,*)
        enddo
        write(6,*)' Imaginary part of unitary matrix='
        do m1=1,wndim
           do m2=1,wndim
              write(6,'(F8.4)',advance='no') aimag(cf(m1,m2,icix))
           enddo
           write(6,*)
        enddo
     endif
  enddo
  DEALLOCATE(wtmp)

  if (Qprint) then
     if (Qident) then
        write(6,*)' All transformation matrices are Identity'
     else
        write(6,*)' At least one transformation matrix is not Identity'
     endif
  endif


  ! For spin-orbit coupling crates transformation "rot_spin_quantization", which rotates from spin quantization axis to global cartesian coordinate system
  !....reading *.inso
  if (so) then
     do i=1,3
        read(4,*)
     end do
     read(4,*)s(1),s(2),s(3)
     call angle(s,theta,fi)
     if (Qprint) write(6,126) theta,fi
     !.... test of collinearity of the z-axis with spin quantization axis
     cost=dcos(theta)
     cosf=dcos(fi)
     sint=dsin(theta)
     sinf=dsin(fi)
     rot_spin_quantization(1,1)=cosf*cost
     rot_spin_quantization(1,2)=-sinf
     rot_spin_quantization(1,3)=cosf*sint
     rot_spin_quantization(2,1)=sinf*cost
     rot_spin_quantization(2,2)=cosf
     rot_spin_quantization(2,3)=sinf*sint
     rot_spin_quantization(3,1)=-sint
     rot_spin_quantization(3,2)=0
     rot_spin_quantization(3,3)=cost
  else
     rot_spin_quantization(:,:)=0.d0
     rot_spin_quantization(1,1)=1.d0
     rot_spin_quantization(2,2)=1.d0
     rot_spin_quantization(3,3)=1.d0
  end if  ! spin-orbit
  
  !---------------------------------------------------------------------------------!
  !--- WRITE Output for the BR1 matrix and the ROTIJ Matrix in file case.rotlm:  ---!
  !---  they are read by LectureROTIJ in lecture.f in SRC_elnes / or SRC_mdff/.  ---!
  !---  These lines are in SRC_lapw2lm/l2main.frc: surch *PH*                    ---!
  !---  The number of the output file is 22. his extention must be .rotlm        ---!
  !---------------------------------------------------------------------------------!
  if (myrank.EQ.master) then

     WRITE(22,*) 'BR1'
     DO JR=1, 3
        WRITE(22, 9101) BR1(1,JR), BR1(2,JR), BR1(3,JR) !--- Writting BR1 ----!
     ENDDO
     WRITE(22,*) 'BR2'
     DO JR=1, 3
        WRITE(22, 9101) BR2(1,JR), BR2(2,JR), BR2(3,JR) !--- Writting BR1 ----!
     ENDDO
 
     INDEX = 0
     DO LATOM=1, NAT
        INDEX1 = INDEX+1
        DO LATEQ=1,MULT(LATOM)
           INDEX = INDEX + 1
           WRITE(22, 9100) LATOM, LATEQ, INDEX !---- inequivalent atomnumber  ---!
           DO JR=1, 3
              WRITE(22, 9101) (ROTIJ(JC,JR,INDEX),JC=1,3) !--- Writting ROTIJ ---!
           ENDDO
           ! Rotates first atom by rotij, the transformation which makes the first and current atom equivalent.
           POST = matmul(POS(:,index1), ROTIJ(:,:,index))
           POST = POST + TAUIJ(:,index) + shft(:,INDEX)
           write(6,'(A,I2,A)',advance='no') 'Actual position of atom ', INDEX, ' is:'
           write(6,'(f10.6,2x,f10.6,2x,f10.6)') POST
        ENDDO
     ENDDO
     write(6,*)
     write(6,'(A)') 'Combined transformation (acting on k-point) for correlated atoms'
     write(6,'(A)') 'including users rotation and internal local rotation'
     call inv_3x3(BR1,BR1inv)
     DO icase=1,natom
        latom = iatom(icase)
        write(6, '(A,I3,1x,A,1x,I3)') 'catom', icase, 'atom', latom
        BKRLOC = matmul(crotloc(:,:,1,icase),matmul(BR1,rotij(:,:,latom)))
        tmp = matmul(BKRLOC,BR1inv)
        DO JR=1,3
           WRITE(6, '(3F10.5)') (tmp(JR,JC),JC=1,3)  ! (BKRLOC(JR,JC),JC=1,3) 
        ENDDO
        if (so) then
           WRITE(6,'(A,I3,A)', advance='no') '  spin-z-quantization axis for corr-atom', icase, ':'
           wdet = detx(rotij(:,:,latom)) ! note that when det<0 we have inversion, which does not change spin
           WRITE(6, '(3F10.5)') (tmp(3,JC)*wdet,JC=1,3)
        endif
        WRITE(6,*)
     ENDDO
     write(6,*)
  endif
  


  ! If SO-coupling and symmetrization over all k-points
  if (so) then
     CALL symoper
     idet(:)=1
     phase(:)=0.d0
     !CALL SYM(THETA,FI) ! This should be corrected!
  endif
  
  CALL INIT_ENERGY_FILE(nkpt, nmat, nume, maxbands,projector,Qcomplex)
  
  if (nemin0.LT.1) nemin0=1
  if (nemax0.GT.nume) nemax0=nume

  !print *, 'Maximum number of bands on rank', myrank,'=',maxbands
  CALL init_abc(nume,nmat,lmax2,ndim2,nrf) ! Allocates variables:  

  !!!!!!
  Cohfacts = .False.
  IF (mode.EQ.'p') THEN
     ! Check if we want to use coherence factors for plotting bands
     OPEN(1001,FILE='cohfactors.dat',STATUS='OLD',ERR=992)
     READ(1001,*,IOSTAT=ios) TEXT
     if (ios /= 0) goto 992
     Cohfacts = .True.
  ENDIF
992 CONTINUE
  !!!!!!
  IF (mode.EQ.'p') mode='e' ! for plotting, we just print eigenvalues

  ! If 'EF.dat' exists, we should use it instead of case.indm2
  OPEN(1000,FILE='EF.dat',STATUS='OLD',ERR=991)
  READ(1000,*) EF
  EF = EF/Ry2eV
  if (Qprint) WRITE(6,*) 'Found EF.dat and EF is ', EF
991 CONTINUE
  if (Qprint) WRITE(6,*) 'EF finally set to ', EF*Ry2eV
  if (Qprint) WRITE(*,'(A,A1,3x,A,F10.7)') 'mode=', mode, 'EF=', EF*Ry2eV
  if (nsymop.EQ.0) nsymop = iord

  !call init_xa(NRAD)

  call l2main(Qcomplex,nsymop,mode,projector,Qrenormalize,fUdmft)  !!!! This calls the main subroutine

  !call fini_xa()
  DEALLOCATE(CF, Sigind, csize, legend)
  DEALLOCATE(shft)
  DEALLOCATE(isort, ifirst)
  DEALLOCATE(crotloc)
  CALL deallocate_kpts
  CALL deallocate_abc()
  call DeallocateStructure()
  call cputim(time1)
  ttime = time1-time0
  
!.....CALCULATE CPUTIME REQUIRED                                        
  !if (Qprint) WRITE(6,2000)                                                     
  if (Qprint) WRITE(6,2010) ttime
  if (myrank.EQ.master) CALL ERRCLR
  
  CALL stop_MPI()

  STOP !'DMFT1 END'
  
  !   error handling
  !
910 INFO = 1
  !
  !   'DMFTx.def' couldn't be opened
  !910
  print *, myrank, 'Definition file error: ', TRIM(DEFFN), ' could not be opened!'
  print *, 'ios=', ios
  WRITE (ERRMSG,9000) DEFFN
  CALL OUTERR('DMFT',ERRMSG)
  GOTO 999

920 INFO = 2
!
!  file FNAME couldn't be opened
!
  WRITE (ERRMSG,9010) IUNIT
  CALL OUTERR('DMFT',ERRMSG)
  WRITE (ERRMSG,9020) FNAME
  CALL OUTERR('DMFT',ERRMSG)
  WRITE (ERRMSG,9030) STATUS, FORM
  CALL OUTERR('DMFT',ERRMSG)
  GOTO 999
960 INFO = 7
!
!        Error reading file 'dmft1.def'
!
  WRITE (ERRMSG,9040) FNAME
  CALL OUTERR('DMFT',ERRMSG)
  GOTO 999
970 INFO = 8
!
! Error in parallel: could not find .processes.x
!
  WRITE (ERRMSG,'(A,A)')  'file open error:', trim(ADJUSTL(FNAME))
  CALL OUTERR('DMFT',ERRMSG)
  GOTO 999
999 STOP 'DMFT - Error'
!                                                                       
!                                                                       

9100 FORMAT('inequivalent atomnumber ',I3,' number ',I2,' total ',I4)
9101 FORMAT(3F10.5)
123 FORMAT(A5)
666 FORMAT(a4)
667 FORMAT(38x,f10.5)
887 FORMAT(' Projected density of states calculation for',i3,' atoms')
877 FORMAT(' Energy window:',f10.4,' < E < ',f10.3)
120 FORMAT(' New z axis || ',3f9.4)
121 FORMAT(' New x axis || ',3f9.4)
1013 FORMAT('users (crotloc) matrix: ',3f10.7,/,24x,3f10.7,/,24x,3f10.7)
2000 FORMAT(//,3X,'=====>>> CPU TIME SUMMARY',/)
2010 FORMAT(12X,'TOTAL    :',F10.2)
9000 FORMAT('can''t open definition file ',A40)
9010 FORMAT('can''t open unit: ',I2)
9020 FORMAT('       filename: ',A50)
9030 FORMAT('         status: ',A,'  form: ',A)
9040 FORMAT('Error reading file: ',A47)
5100 FORMAT(A10,4I5,3F5.2,A3)
5101 FORMAT(A10,4I10,3F5.2,A3)
126 FORMAT(' Global spin quantization axis is; M theta, phi:',2f8.4) 
END PROGRAM DMFTMAIN


SUBROUTINE INIT_ENERGY_FILE(nkpt, nmat, nume, maxbands,projector,Qcomplex)
  !*****  finds nkpt, nmat and nume in energy file *****
  ! 59 -> case.energydn
  ! 60 -> case.energy
  ! using case.energy, we will find out what is:
  !  nkpt -- number of all k-points
  !  nmat -- maximum number of reciprocal vectors
  !  nume -- maximum number of bands kept
  USE com_mpi,only: nvector, vector_para, fvectors, myrank, master, VECFN, Qprint
  USE structure,only : nat
  USE param, only: nemin0, nemax0
  USE com, only: EMIN,EMAX
  IMPLICIT NONE
  !
  INTEGER, intent(out):: nkpt, nmat, nume, maxbands
  INTEGER, intent(in) :: projector
  LOGICAL, intent(in) :: Qcomplex
  !
  REAL*8        :: sumw0
  CHARACTER*200 :: energy_filename, vector_filename
  LOGICAL       :: energy_file_exists
  CHARACTER*10  :: KNAME
  INTEGER       :: k, ivector, itape, ios
  INTEGER       :: N, NEn, NUM, I, ii, nemin, nemax
  REAL*8        :: EMIST, SS, TT, ZZ, wgh, E1
  INTEGER       :: wkx, wky, wkz, nbands
  REAL*8        :: wAsr
  COMPLEX*16    :: wAsc
  !--------------------------------------------------------------------- 
  !find nkpt, nmat and nume in energy file
  k=0
  nmat=0
  nume=0
  maxbands=0  ! output
  sumw0=0
  DO ivector=1,nvector
     if (vector_para) then
        energy_filename = fvectors(ivector, 3)
        vector_filename = fvectors(ivector, 1)
     else
        energy_filename = VECFN(4)
        vector_filename = VECFN(1)
     endif
     !INQUIRE(FILE=energy_filename, EXIST=energy_file_exists)
     energy_file_exists=.False.
     if (energy_file_exists) then
        itape=60
        open(itape,FILE=energy_filename,STATUS='old',FORM='formatted')
        DO I=1,NAT
           READ(itape,'(f9.5)') EMIST
           READ(itape,'(f9.5)') EMIST
        ENDDO
        ios=0
        DO WHILE (ios == 0)
           READ(itape,'(3e19.12,a10,2i6,F5.1)',IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
           IF (ios /= 0) CYCLE
           k=k+1
           nmat=MAX(n,nmat)
           nume=MAX(nen,nume)
           NEMIN=1
           NEMAX=0
           DO ii=1,NEn
              READ(itape,*) NUM,E1
              if (abs(projector).LT.4) then
                 IF(E1.LT.EMIN) NEMIN=NEMIN+1
                 IF(E1.LT.EMAX) NEMAX=NEMAX+1
              else
                 nemin = nemin0
                 nemax = nemax0
              endif
           ENDDO
           nbands = NEMAX-NEMIN+1
           maxbands = max(maxbands,nbands)
           sumw0=sumw0+wgh
        ENDDO
        close(itape)
     else
        itape=9
        if (vector_para) then
           open(itape,FILE=vector_filename,STATUS='old',FORM='unformatted')
        else
           rewind(itape)
        endif
        DO i=1,nat
           READ(itape) emist
           READ(itape) emist
        enddo
        ios=0
        DO WHILE (ios == 0)
           READ(itape,IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
           !print *, 'GOT', SS, TT, ZZ, KNAME, ios
           IF (ios /= 0) CYCLE
           k=k+1
           nmat=MAX(n,nmat)
           nume=MAX(nen,nume)
           READ(itape) (wkx,wky,wkz,i=1,N)
           NEMIN=1   ! bug 2021
           NEMAX=0   ! bug 2021
           DO ii=1,NEn 
              READ(itape) NUM,E1
              if (abs(projector).LT.4) then
                 IF(E1.LT.EMIN) NEMIN=NEMIN+1
                 IF(E1.LT.EMAX) NEMAX=NEMAX+1
              else
                 nemin = nemin0
                 nemax = nemax0
              endif
              if (Qcomplex) then
                 READ(itape) (wAsc,i=1,N)
              else
                 READ(itape) (wAsr,i=1,N)
              endif
           ENDDO
           nbands = NEMAX-NEMIN+1
           maxbands = max(maxbands,nbands)
           sumw0=sumw0+wgh
        ENDDO
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     endif
  END DO
  nkpt=k
  if (Qprint) WRITE(6,*) 'nume=', nume, 'nkpt=', nkpt
end SUBROUTINE INIT_ENERGY_FILE

subroutine swapGroup(iz1,tau1,iz2,tau2)
  IMPLICIT NONE
  INTEGER, intent(inout) :: iz1(3,3), iz2(3,3)
  REAL*8,  intent(inout) :: tau1(3), tau2(3)
  ! locals
  INTEGER :: izx(3,3)
  REAL*8  :: taux(3)
  izx(:,:) = iz1(:,:)
  iz1(:,:) = iz2(:,:)
  iz2(:,:) = izx(:,:)
  taux(:) = tau1(:)
  tau1(:) = tau2(:)
  tau2(:) = taux(:)
end subroutine swapGroup

subroutine FindMax_nl(natom, max_nl)
  use com_mpi, only: stop_MPI
  use error, only: OUTERR
  IMPLICIT NONE
  INTEGER, intent(out):: max_nl
  INTEGER, intent(in) :: natom
  ! locals
  CHARACTER*200 :: ERRMSG
  INTEGER :: i, j, il, latom,nL,locrot, loro, shift
  max_nl=0
  DO i=1,natom
     READ(5,*) latom,nL,locrot ! read from case.indmfl
     if (nL > max_nl) max_nl = nL
     loro = MOD(locrot, 3)
     shift = locrot/3
     ! loro can be (0 -original coordinate system), (1 - new z-axis), (2 - new z-axis, x-axis)
     do j=1,nL
        read(5,*) !LL(i,j), qsplit(i,j), cix(i,j) ! LL contains all L-quantum numbers which need to be computed.
     enddo
     if (loro.eq.0) then
        ! nothing to do
     elseif(loro.gt.0)then        ! change of local rotation matrix
        print *, 'positive locrot not supported anymore. Please use locrot=-1,-2,-3,...'
        call stop_MPI
        STOP 'ERROR dmft1'
     elseif (loro.eq.-1) then               ! in cartesian, and equal for all l cases
        read(5,*)  !( crotloc(1,j,1,i),j=1,3)  ! new x-axis expressed in unit cell vectors
        read(5,*)  !( crotloc(2,j,1,i),j=1,3)  ! new y-axis expressed in unit cell vectors
        read(5,*)  !( crotloc(3,j,1,i),j=1,3)  ! new z-axis expressed in unit cell vectors
     elseif(loro.eq.-2) then
        do il=1,nL
           read(5,*) !( crotloc(1,j,il,i),j=1,3)  ! new x-axis expressed in unit cell vectors
           read(5,*) !( crotloc(2,j,il,i),j=1,3)  ! new y-axis expressed in unit cell vectors
           read(5,*) !( crotloc(3,j,il,i),j=1,3)  ! new z-axis expressed in unit cell vectors
        enddo
     else
        print *, 'locrot=', locrot, 'with loro=', loro, 'and shift=',shift, 'not implemented'
        WRITE(ERRMSG,*) 'locrot=', locrot, 'with loro=', loro, 'and shift=',shift, 'not implemented'
        call OUTERR('dmft1', ERRMSG)
     endif
     ! end Rotation matrix loro
     if (shift.ne.0) then
        read(5,*) !(shft(j,latom),j=1,3)
     endif
  enddo
  rewind(5)
  READ(5,*) !EMIN,EMAX,irenormalize,projector
  READ(5,*) !imatsubara, gammac, gamma, nom_default, aom_default, bom_default
  read(5,*) !natom
end subroutine FindMax_nl
