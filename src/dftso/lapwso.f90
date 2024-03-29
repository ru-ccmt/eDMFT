! TODO: reading Vns should be done correctly with atpar library. 
PROGRAM lapwso
  USE abcd, ONLY: allocate_abcd
  USE loabc, ONLY: elo, allocate_loabc
  USE loabcr, ONLY: elor, allocate_loabcr
  USE lolog, ONLY: lapw, allocate_lolog, lapw, ilo, nlov, nlon, nlo, loor
  USE rlolog, ONLY: nnrlo, allocate_rlolog, loorext, nrlov, nnrlo, nrlo
  USE rpars, ONLY: allocate_rpars
  USE orb, ONLY: iorbpot, allocate_orb
  USE structure, ONLY: ReadStructure, ReadStructureGroup, ndf, WriteInfoStructure, jri, rotij, tauij, alat, alpha, BR1, BR2, imat, iord, lattic, ortho, pos, vol, tau, mult
  USE rotmat, ONLY: allocate_rotmat, det
  USE vns, ONLY: allocate_vns
  USE ams, ONLY: init_ams
  USE hmsout, ONLY: deallocate_hmsout, allocate_hmsout
  USE meigve_mod, ONLY: meigve
  USE param, ONLY: hblock, labc, labc2, lmax, lomax, nato, ndif, nloat, nmat, nrad, num2, nume, nume2, filename_V_sp !,  , filename_V_sp, filename_V_vns
  USE mpi, ONLY: Qprint, FilenameMPI, FilenameMPI2, start_MPI, stop_MPI, myrank, master, vector_para, nprocs, filename_vector, filename_energy, filename_vectorso, filename_energyso, filename_norm, ikps, Scatter_Vector_data, mpi_bcast_V_vsp, pr_proc, WriteProcesses
  use readPotential, ONLY: init_V_vsp, read_V_vsp
  USE couples, ONLY: allocate_couplo
  USE store_vec, ONLY: vec
  IMPLICIT NONE
  interface
     subroutine Read_Next_Kpoint(SS, bname, NV, ne, weight, KV, ee, meigve, should_store, jspin, Qcomplex)
       integer, intent(in) :: jspin
       logical, intent(in) :: Qcomplex, should_store
       real*8, intent(out) :: SS(:,:), weight(:) ! SS(3,2), weight(2)
       character*10, intent(out) :: bname
       integer, intent(out):: nv(:), ne(:)     !nv(2), ne(2)
       integer, intent(out):: KV(:,:,:)        !KV(3,nmat,2)
       real*8, intent(out) :: ee(:,:)          !ee(nume,2)
       complex*16, intent(out):: meigve(:,:,:) ! meigve(nmat,nume,jspin)
     end subroutine Read_Next_Kpoint
  end interface
  real*8, parameter :: pi = 3.141592653589793238462643d0
  character*80   :: deffn,errfn!,fname
  CHARACTER*67   :: ERRMSG
  CHARACTER*180  :: FNAME
  character*10   :: bname
  logical        :: fl,Qcomplex, rel, must_compute, vector_para_so_only
  INTEGER :: nv(2),ne(2)
  REAL*8 :: emm(2)  ! emin, emax
  REAL*8,ALLOCATABLE     :: Vru(:,:,:),e(:,:,:)
  REAL*8,ALLOCATABLE     :: p(:,:,:),dp(:,:,:),pe(:,:,:),dpe(:,:,:)
  real*8,ALLOCATABLE     :: ri_mat(:,:,:,:,:,:),ee(:,:)
  real*8,ALLOCATABLE     :: ri_orb(:,:,:,:,:,:)
  INTEGER,ALLOCATABLE    :: kv(:,:,:)
  REAL*8 ::  ss(3),cp(5)
  integer :: nkp, iscf, index, ierr
  real*8 :: WEIGHT_(2)
  real*8 :: SS_(3,2), xms(3)
  real*8 :: b2Mb, dt0, dtime0, dtime1, dt4
  real*8 :: phi, theta, weight, restsize
  integer :: ii, i, j, ios, ipr, irlotot, isi, ispin, jatom, jspin, ikp, iikp, kpot, l, n_scr, nban2, nn, num, info, lxdos
  integer, allocatable :: kpoints(:)
  REAL*8, allocatable  ::  APA(:)
  !
  vector_para_so_only=.False.
  !vector_para_so_only=.True.
  CALL start_MPI()
  Qcomplex=.false.
  call gtfnam(deffn,errfn)
  CALL ERRFLG(ERRFN,'Error in LAPWSO')

  if (nprocs.EQ.1) then
     write(6,*)'Running lapwso in single processor mode'
     !write(6,*)' '
  else if (Qprint) then
     write(6,*)'Running lapwso in mpi mode'
     !write(6,*)' '
  endif

  ! Reading lapwso.def  file *
  CALL read_def_file(deffn, Qcomplex, iorbpot,ERRMSG,info)
  if (info.ne.0) then
     CALL OUTERR('lapwso.f90',ERRMSG)
  endif

  if (vector_para) then
     call Scatter_Vector_data(ierr)

     if (vector_para_so_only) then
        ! In this mode, we read from a single vector file case.vector, but we print to
        ! multiple vectorso files, so that each process produces its own vectorso_{myrank} file.
        !print *, 'filename before=', filename_energy(1)
        call Remove_xend(filename_vector(1))
        call Remove_xend(filename_vector(2))
        call Remove_xend(filename_energy(1))
        call Remove_xend(filename_energy(2))
        !print *, 'filename after=', filename_energy(1)
     endif
     
     if (ierr.ne.0) then
        ! if the file _processes_ does not exist, we just have a single vector file
        vector_para = .False.
        ! correct filenames if they were wrong in def file
        call add_xend('', filename_vector(1))    !'case.vector_x'
        call add_xend('', filename_vector(2))    !'case.vectorup_x'
        call add_xend('', filename_energy(1))    !'case.energy_x'
        call add_xend('', filename_energy(2))    !'case.energyup_x'
        call add_xend('', filename_vectorso(1))  !'case.vectorsodn'
        call add_xend('', filename_vectorso(2))  !'case.vectorso'
        call add_xend('', filename_energyso(1))  !'case.energysodn'
        call add_xend('', filename_energyso(2))  !'case.energyso'
        call add_xend('', filename_energyso(3))  !'case.energydum'
        call add_xend('', filename_norm(1))      !'case.normsodn'
        call add_xend('', filename_norm(2))      !'case.normsoup'
     endif
     ! Opens all necessary files
     ! files to read
     FNAME=filename_vector(1)
     OPEN (9, FILE=filename_vector(1),STATUS='old',FORM='unformatted',ERR=920)       ! vector
     FNAME=filename_vector(2)
     OPEN (10,FILE=filename_vector(2),STATUS='unknown',FORM='unformatted',ERR=920)   ! vectorup
     FNAME=filename_energy(1)
     OPEN (54,FILE=filename_energy(1),STATUS='old',FORM='formatted',ERR=920)         ! energy
     FNAME=filename_energy(2)
     OPEN (55,FILE=filename_energy(2),STATUS='unknown',FORM='formatted',ERR=920)     ! energyup
     ! files to write
     FNAME=filename_vectorso(1)
     OPEN (41,FILE=filename_vectorso(1),STATUS='unknown',FORM='unformatted',ERR=920) ! vectorsodn
     OPEN (42,FILE=filename_vectorso(2),STATUS='unknown',FORM='unformatted',ERR=920) ! vectorso
     OPEN (51,FILE=filename_energyso(1),STATUS='unknown',FORM='formatted',ERR=920)   ! energysodn
     OPEN (52,FILE=filename_energyso(2),STATUS='unknown',FORM='formatted',ERR=920)   ! energyso
     OPEN (53,FILE=filename_energyso(3),STATUS='unknown',FORM='formatted',ERR=920)   ! energydum
     OPEN (45,FILE=filename_norm(1),STATUS='unknown',FORM='formatted',ERR=920)       ! normsodn
     OPEN (46,FILE=filename_norm(2),STATUS='unknown',FORM='formatted',ERR=920)       ! normsoup
  endif
  
  call ReadStructure(20,nato,rel,lxdos,ERRMSG,info)
  ndif = ndf
  if (info.ne.0) then
     CALL OUTERR('structure.f90',ERRMSG)
  endif
  if (Qprint) CALL WriteInfoStructure(6, nato)
  CALL ReadStructureGroup(20, nato)
  
  call get_new_nloat(4,nloat)
  nloat = nloat + 1 ! Here we have different definition on nloat. Some inconsistency with the rest of the code.
  !call get_nloat(lomax,nato,nloat)

  allocate( elo(0:lomax,nloat,nato,2) )
  ALLOCATE( e(0:lmax,nato,2) )
  if (vector_para_so_only) vector_para=.False.
  CALL getMaxDim(nmat,nume,nkp,jspin,e,elo, Qcomplex)

  if (vector_para_so_only) then
     vector_para=.True.
     pr_proc = floor(nkp/DBLE(nprocs)+0.999)
     ikps(2)=pr_proc
     ikps(1)=myrank*pr_proc
  endif
  if (vector_para) then
     allocate( kpoints(0:pr_proc) )
     kpoints(0)=0     ! this is just for convenience in the below k-point loop
     do i=0,pr_proc-1 ! k-points are distributed in the simplest way: we know the first point (ikps(1)) and the number of points (ikps(2))
                      ! BUG corrected Nov 30, 2016  (Walber). This loop was going up to pr_proc rather than pr_proc-1, and it would coredump
        kpoints(i+1) = ikps(1)+i+1
     enddo
     !kpoints(1) = ikps(1)+1
     !kpoints(pr_proc) = ikps(1)+pr_proc
  else
     ! We do not have the redistribution of k-points over processors yet. We need to provide it here.
     pr_proc  = floor(nkp/DBLE(nprocs)+0.999)  ! The maximum number of points calculated per processor     
     allocate( kpoints(0:pr_proc) )
     ! correcting number of points if the number of k-points is not divisible
     ikps(2) = pr_proc            ! number of k-points to be computed
     ikps(1) = 1+myrank           ! what is the first k-point on this processor
     ! If we output into a common file, the second core needs to calculate kpoint number 2, and third processor kpoint 3, etc...
     kpoints(0)=0  ! this is just for convenience in the below k-point loop
     do i=0,pr_proc-1
        kpoints(i+1) = 1 + myrank + i*nprocs
     enddo
  endif
  if (Qprint) WRITE(6,'(A,I4,1x,A,I4)') 'pr_proc=', pr_proc, 'tot-k=', nkp
  !WRITE(6,*) myrank, 'pr_proc=', pr_proc
  !WRITE(6,*) myrank, 'k-pnts:', kpoints(:)
  !WRITE(6,*) myrank, 'ikps:', ikps(:)
  !print *, 'rank=', myrank, 'pr_proc=', pr_proc, 'tot-k=', nkp, 'kpoints=', kpoints
  !ikps(1) = 1     ! what is the first k-point in this vector file
  !ikps(2) = nkp   ! number of k-points
  !pr_proc = nkp


  do jatom=1,nato
     if (vector_para .or. myrank.eq.master) then
        ! write also on dummy vector file 
        write(53,'(100(f12.5))')(e(l,jatom,1),l=0,lmax) 
        write(53,'(100(f12.5))')((elo(l,nn,jatom,1),l=0,lomax),nn=1,nloat-1)
     endif
  enddo

  if (Qprint) then
     write(6,*)'*********** spin-orbit calculation **************'
     if(jspin.eq.1)write(6,*)' non-spin polarized case'
     if(jspin.eq.2)write(6,*)' spin polarized case'
  endif
  
  CALL init_ams
  CALL allocate_rotmat(ndif)
  CALL allocate_loabc(lomax,nloat,nato)
  CALL allocate_loabcr(lomax,nato)
  CALL allocate_lolog(nato,lomax,lmax)
  CALL allocate_rlolog(nato,lomax)
  CALL allocate_rpars(nato,lomax)
  CALL allocate_orb(ndif)
  CALL allocate_vns(nato)
  
  call cputim(dt0)

  allocate( Vru(nrad,nato,jspin) )
  ! Opens file with spherically symmetric potential
  if (myrank.eq.master) then
     Vru(:,:,:)=0.d0
     do isi=1,jspin
        CALL init_V_vsp(filename_V_sp(isi), 18+isi-1, ISCF)
        do jatom=1,nato  ! Reads the spherically symmetric potential
           CALL read_V_vsp(Vru(:jri(jatom),jatom,isi),jatom,jri(jatom))  ! This potential is in Rydbergs
        enddo
        close(18+isi-1)
     enddo
     !Vru(:,:) = Vru(:,:)/2.d0    ! Converst from Rydbergs to Hartrees
  endif
  call mpi_bcast_V_vsp(Vru,nrad,nato,jspin)
  
  call read_inso(emm,xms,fl,jspin,kpot,ipr,irlotot) !

  !.....Generate symmop
  call allocate_couplo(ndif,labc)
  CALL generate_lattice(alpha, BR1, BR2, VOL, ORTHO, ALAT(1), ALAT(2), ALAT(3), LATTIC)
  if (Qprint) then
     WRITE(6,*) 'BR1='
     WRITE(6,'(3F12.6)') BR1
     WRITE(6,*) 'BR2='
     WRITE(6,'(3F12.6)') BR2
  endif
  ! define rotation matrices if required
  allocate( rotij(3,3,ndif) )
  allocate( tauij(3,ndif) )
  call get_rotij_tauij(rotij,tauij,pos,alat,imat,tau,iord,nato,ndif,mult,lattic)
  ! direction of the spin quantization axis
  call angle(xms,theta,phi)
  ! group operations with SO
  call symop(theta,phi)

  if (Qprint) then
     index=0
     do i=1,nato
        do j=1,mult(i)
           index=index+1
           if(det(index).gt.0.)then
              write(6,557)i,j
           else
              write(6,558)i,j
           endif
557        format(' for atom type',i3,' eq. atom',i3,' sym. oper. conserves spin')
558        format(' for atom type',i3,' eq. atom',i3,' sym. oper.   inverts spin')
        end do
     end do
     write(8,111) (theta*180/pi),(phi*180/pi)
     write(6,111) (theta*180/pi),(phi*180/pi)
111  format(2f6.1,' angle (M,z), angle (M,x) deg')
  endif
  ! How many and which local orbitals present
  call Find_nlos(lapw, loor, ilo, nlo, nlov, nlon, nrlo, nrlov, nnrlo, loorext, jspin, nato, lmax, lomax, nloat, e, elo)
  ! For LDA+U or similar, we read orbital potential here
  call read_orbital_potential()

  do ispin=1,2
     do jatom=1,nato
        do l=0,lomax
           elor(l,jatom,ispin)=1.e+6
        enddo
     enddo
  enddo
  
  ALLOCATE(p(labc+1,nato,2),dp(labc+1,nato,2),pe(labc+1,nato,2),dpe(labc+1,nato,2))
  ALLOCATE(ri_mat(nloat,nloat,0:labc,nato,2,2))
  ALLOCATE(ri_orb(nloat,nloat,0:labc,nato,2,2))
  
  if (Qprint) write(6,*)nmat,nnrlo
  nmat=nmat+nnrlo
  nume=nume+nnrlo
  nume2=  2*nume
  num2= nume2*(nume2+1)/2+(nume2/hblock+1)*(hblock*(hblock-1)/2)   
  b2Mb=1.d0/(1024.d0**2.d0)
  
  allocate(kv(3,nmat,2),ee(nume,2))
  
  call garadme(e,Vru,p,dp,pe,dpe,ri_mat,jspin,kpot,ipr)
  
  if(iorbpot.ne.0)then
     call garadorb(e,Vru,p,dp,pe,dpe,ri_orb,jspin,kpot,ipr)
  endif
  call gaunt2
  
  call cputim(dtime1)

  !memory usage
  if (Qprint) then
     write(6,*) "-------------memory usage-----------------"
     write(6,"(A,F15.3,A)") "meigve:",16.0d0*nmat*nume* 2.d0*b2Mb," MB"
     write(6,"(A,F15.3,A)") "abclm:",16.0d0*(nloat)*labc2*nume*2.d0*b2Mb," MB"
     write(6,"(A,F15.3,A)") "h_:",16.0d0*nume2*nume2*b2Mb," MB"
     if (nnrlo.eq.0) then
        write(6,"(A,F15.3,A)") "s_:",16.0d0*1.d0*b2Mb," MB"
     else
        write(6,"(A,F15.3,A)") "s_:",16.0d0*nume2*nume2*b2Mb," MB"
     endif
     write(6,"(A,F15.3,A)") "vec:",16.0d0*nume2*nume2*b2Mb," MB"
     write(6,"(A,F15.3,A)") "vect:",16.0d0*nmat*nume2*2d0*b2Mb," MB"
     restsize=(4*(77*nato+9*ndif+3*lomax*nato) + 8*(157*nato+ 14*ndif+22*lomax*nato*nloat+28*lomax*nato+2*nrad*nato+8*labc*nato+3*nume2+8*nloat*labc*nloat*nato) + 16*(441*ndif+16*ndif*labc*labc*labc))*b2Mb
     write(6,"(A,F15.3,A)") "rest:",restsize," MB"
     if (nnrlo.eq.0) then
        write(6,"(A)") "Totmem1=meigve+abclm+h_"
        write(6,"(A,F15.3,A)") ":TOTMEM1",(nmat*nume*2.d0+(nloat)*labc2*nume*2.d0 + nume2*nume2)*16.0d0*b2Mb+restsize," MB"
        write(6,"(A)") "Totmem2=meigve+vec+vect+h_"
        write(6,"(A,F15.3,A)") ":TOTMEM2",(nmat*nume*2.d0+nume2*nume2+nmat*nume2*2d0+nume2*nume2)*16.0d0*b2Mb+restsize," MB"
     else
        write(6,"(A)") "Totmem1=meigve+abclm+h_+s_"
        write(6,"(A,F15.3,A)") ":TOTMEM1",( nmat*nume*2.d0+(nloat)*labc2*nume*2.d0+nume2*nume2+nume2*nume2)*16.0d0*b2Mb+restsize," MB"
        write(6,"(A)") "Totmem2=meigve+vec+vect+h_+s_"
        write(6,"(A,F15.3,A)") ":TOTMEM2",(nmat*nume*2.d0+nume2*nume2+nmat*nume2*2d0+nume2*nume2+nume2*nume2)*16.0d0*b2Mb+restsize," MB"
     endif
     write(6,*) "-------------memory usage-----------------"
  endif
  !memory usage


  !*** Writing linearization energies into vectorso files
  if (vector_para .or. myrank.eq.master) then
     do ii=1,2
        isi=ii
        if (jspin.eq.1) isi=1
        do jatom=1,nato
           do l=0,lmax
              if(.not.lapw(l,jatom)) e(l,jatom,isi)=e(l,jatom,isi)+200.d0 ! Correcting energies back before storing them
           enddo
           write(40+ii)(e(l,jatom,isi),l=0,lmax)                        ! Writing into vectorso file
           write(40+ii)((elo(l,nn,jatom,isi),l=0,lomax),nn=1,nloat-1)   ! Writing into vectorso file
           !
           if (nloat-1.gt.3) then                                      ! writting into energyso file
              write(50+ii,'(100(f12.5))')(e(l,jatom,isi),l=0,lmax)   ! some more precision when more local orbitals
              write(50+ii,'(100(f12.5))')((elo(l,nn,jatom,isi),l=0,lomax),nn=1,nloat-1)
           else
              write(50+ii,'(100(f9.5))')(e(l,jatom,isi),l=0,lmax)
              write(50+ii,'(100(f9.5))')((elo(l,nn,jatom,isi),l=0,lomax),nn=1,nloat-1)
           endif
           do l=0,lmax
              if(.not.lapw(l,jatom)) e(l,jatom,isi)=e(l,jatom,isi)-200.d0 ! Correcting lo-energies for atpar
           enddo
        enddo
     enddo
  endif

  if (.not.Qcomplex) allocate( APA(nmat) )
  
  !print *, myrank, 'nkp=', nkp, 'nmat=', nmat, 'nume=', nume
  DO iikp=1,pr_proc
     ! Here we could just loop over the k-points, which need to be computed on this processor.
     ! However, we loop over the same number of k-points on each processor, even when some processor have nothing to do (if number of k-points is not commensurate with the number of processors)
     ! This is because of the MPI communication, which needs to be initiated on all processors, even on those that have nothing to do.
     ! All processes have to be enter equal number of times to MPI_send/receive, so that we can call MPI_Barrier. The work is done only if needed.
     ikp = kpoints(iikp)
     must_compute = (ikp .le. nkp)
     
     if (.not.vector_para .or. vector_para_so_only) then
        ! We need to read all k-points in-between, which will not be computed on this processor
        do j=kpoints(iikp-1)+1,kpoints(iikp)-1  ! do not store
           if (j > nkp) exit
           call Read_Next_Kpoint(SS_, bname, NV, NE, weight_, KV, ee, meigve, .False., jspin, Qcomplex)
           !print *, 'On proc=',myrank,'reading a point'
        enddo
     endif

     if (must_compute) then
        ! This k-point has to be computed here
        ALLOCATE(meigve(nmat,nume,jspin))
        CALL allocate_abcd(labc2,nloat,nume)
        
        call cputim(dtime0)

        call Read_Next_Kpoint(SS_, bname, NV, NE, weight_, KV, ee, meigve, .True., jspin, Qcomplex)
        !print *, ikp, 'computed on processor', myrank, bname

        call cputim(dtime1)
        cp(1)=cp(1)+dtime1-dtime0
        
        weight=WEIGHT_(1) ! both up&down should have the same k-points
        SS=SS_(:,1)       ! both up&down should have the same k-points
        if (jspin.eq.1) then
           NV(2) = NV(1)
           NE(2) = NE(1)
        endif
        nban2=ne(1)+ne(2)+2*nnrlo
        n_scr=ne(1)+ne(2)
        
        if(.not.fl) then  ! Only eigenvalues but not eigenvectors
           nv(1)=0
           nv(2)=0
        endif
        call cputim(dtime0)
        call hmsec(fl,emm,ne,nv,ee,ri_mat,ri_orb,jspin,ipr,ss,kv,p,dp,pe,dpe)
        call cputim(dtime1)
        cp(3)=cp(3)+dtime1-dtime0
        deallocate(vec)
     else
        CALL allocate_hmsout(nmat,nume2)
     endif
  
     call cputim(dtime0)     
     call kptout(ss,bname,weight,ikp,kv,jspin,nv,ne,must_compute)
     call cputim(dtime1)
     cp(4)=cp(4)+dtime1-dtime0
     
     !if (must_compute) then
     call deallocate_hmsout
     !endif
     
  ENDDO
  call cputim(dt4)
  deallocate( kpoints )
  if (.not.Qcomplex) deallocate( APA )
  if (vector_para) call WriteProcesses()
    
  if (Qprint) then
     write(6,*) 'TOTAL NUMBER OF K-POINTS:',nkp
     write(6,*) 'TOTAL TIME:',dt4-dt0
     write(6,*) 'read:',cp(1)
     write(6,*) 'hmsec:',cp(3)
     write(6,*) 'write:',cp(4)
     write(6,*) ':TOTTIME on proc 0:',dt4-dt0
  endif
  CALL ERRCLR(ERRFN)

  call stop_MPI
  STOP 'LAPWSO END'

920 CONTINUE
  WRITE(6,*) ' ERROR IN OPENING FILE:', TRIM(ADJUSTL(FNAME))
  ERRMSG = 'ERROR READ/WRITE :'//TRIM(ADJUSTL(FNAME))
  CALL OUTERR('lapwso',ERRMSG)
998 CONTINUE
  WRITE(6,*) ' ERROR IN READING VECTOR FILE. likely not enough k-points'
  ERRMSG = ' ERROR IN READING VECTOR FILE. likely not enough k-points'
  CALL OUTERR('lapwso',ERRMSG)
END PROGRAM lapwso


subroutine Read_Next_Kpoint(SS, bname, NV, ne, weight, KV, ee, meigve, should_store, jspin, Qcomplex)
  use param, only: nmat, nume
  use mpi, only: myrank, stop_mpi
  IMPLICIT NONE
  integer, intent(in) :: jspin
  logical, intent(in) :: Qcomplex, should_store
  real*8, intent(out) :: SS(:,:), weight(:) ! SS(3,2), weight(2)
  character*10, intent(out) :: bname
  integer, intent(out):: nv(:), ne(:)     !nv(2), ne(2)
  integer, intent(out):: KV(:,:,:)        !KV(3,nmat,2)
  real*8, intent(out) :: ee(:,:)          !ee(nume,2)
  complex*16, intent(out):: meigve(:,:,:) ! meigve(nmat,nume,jspin)
  ! locals
  INTEGER :: isi, num, i, j, ios
  CHARACTER*67   :: ERRMSG
  REAL*8, allocatable  ::  APA(:)
  complex*16 :: meig_dum
  real*8 :: apa_dum
  !
  if ((.not.Qcomplex) .and. should_store) allocate( APA(nmat) )


  do isi=1,jspin
     ! Reading the vector file, which was obtained by lapw1
     READ(8+isi,IOSTAT=ios,err=997) SS(1,isi),SS(2,isi),SS(3,isi),BNAME,NV(isi),NE(isi),WEIGHT(isi)
     IF(NV(isi).GT.NMAT) write(6,*) 'TOO MANY KJS, NV.GT.NMAT ',nv(isi),nmat
     IF(NE(isi).GT.NUME) write(6,*) 'TOO MANY BANDS, NE.GT.NUME',ne(isi),nume
     IF(NV(isi).GT.NMAT)  STOP 'TOO MANY KJS: '
     IF(NE(isi).GT.NUME)  STOP 'TOO MANY BANDS: '

     READ(8+isi)(KV(1,i,isi),KV(2,i,isi),KV(3,i,isi),i=1,NV(isi))
     
     DO j=1,NE(isi)
        READ(8+isi) NUM,ee(j,isi)
        
        IF(Qcomplex) THEN 
           if (should_store) then
              READ(8+isi) (meigve(i,j,isi),I=1,NV(isi))
           else
              READ(8+isi) (meig_dum,I=1,NV(isi)) ! do not store
           endif
        ELSE
           if (should_store) then
              READ(8+isi)(APA(I),I=1,NV(isi)) ! temporary real storage
              DO i=1,NV(isi)
                 meigve(i,j,isi)=dcmplx(APA(i),0.D0)
              ENDDO
           else
              READ(8+isi)(apa_dum,I=1,NV(isi)) ! do not store
           endif
        ENDIF
     ENDDO
  enddo
        
  if ((.not.Qcomplex) .and. should_store) deallocate( APA )
  return
997 CONTINUE
  WRITE(6,*) ' ERROR IN READING VECTOR FILE. likely not enough k-points', SS(1,1), SS(2,1), SS(3,1), BNAME, NV(1), NE(1), WEIGHT(1)
  ERRMSG = ' ERROR IN READING VECTOR FILE. likely not enough k-points'
  CALL OUTERR('lapwso',ERRMSG)
end subroutine Read_Next_Kpoint


subroutine Remove_xend(filnm)
  IMPLICIT NONE
  CHARACTER*180, intent(inout) :: filnm
  CHARACTER*1 :: c
  INTEGER :: ln, i
  ln=len_trim(filnm)
  do i=ln,1,-1
     c=filnm(i:i)
     if (.not.(c.eq.'0'.or.c.eq.'1'.or.c.eq.'2'.or.c.eq.'3'.or.c.eq.'4'.or.c.eq.'5'.or.c.eq.'6'.or.c.eq.'7'.or.c.eq.'8'.or.c.eq.'9'.or.c.eq.'_')) then
        exit
     endif
  enddo
  filnm = filnm(1:i)
end subroutine Remove_xend
  
