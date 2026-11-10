 MODULE defs
  REAL*8,PARAMETER             :: zero= 0.0d0
  COMPLEX*16,PARAMETER         :: zeroc= (0.0d0,0.0d0)
  REAL*8             :: tlev=2D-9
  REAL*8             :: zlevel=1D-15
  REAL*8             :: mgga_rhocut=1.d-5
  REAL*8             :: mgga_vcut=0.5
END MODULE defs

MODULE parallel
  use :: ISO_C_BINDING, only: C_INTPTR_T

  ! Generic OpenMPI parameters
  INTEGER             :: ierr,npe,myid,npe_coul,myid_coul
  INTEGER,ALLOCATABLE :: nstart(:),nstop(:)
  INTEGER             :: comm_coul,comm_xc
  INTEGER             :: icoul
  INTEGER,ALLOCATABLE :: statusmpi(:)
  LOGICAL             :: coul,xc
  
  CONTAINS
  
  SUBROUTINE init_parallel
#ifdef Parallel
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    ALLOCATE(statusmpi(MPI_STATUS_SIZE))
    CALL MPI_INIT(ierr)
    IF(ierr.NE.MPI_SUCCESS) THEN
       PRINT*, 'MPI ini error'
       call abort_parallel
       STOP
    ENDIF
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, npe, ierr)
    IF(ierr.NE.MPI_SUCCESS) THEN
       PRINT*, 'MPI_comm_size error'
       call abort_parallel
       STOP
    ENDIF
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
#else
    myid=0
    npe=1
    myid_coul=0
    npe_coul=1
#endif
    ALLOCATE(nstart(0:npe-1),nstop(0:npe-1))
  END SUBROUTINE init_parallel

        SUBROUTINE abort_parallel
#ifdef Parallel
        IMPLICIT NONE
        include 'mpif.h'
        INTEGER :: IERR
!       Flush outputs
        call flush(6)
        close(unit=6)
        if(myid .eq. 0)then
                call flush(99)
                close(unit=99)
        endif
        CALL MPI_ABORT(MPI_COMM_WORLD, IERR)
#endif
        END SUBROUTINE abort_parallel

  SUBROUTINE begend(nat)
#ifdef Parallel
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER nat,modu,i,icolor1
    icolor1=MPI_UNDEFINED
    coul=.FALSE.
    xc=.true.
    modu=(nat-1)/npe+1
    icoul=0
    icolor1=MPI_UNDEFINED
    DO i=0,npe-1
       nstart(i)=i*modu+1
       nstop(i)=MIN(nat,(i+1)*modu)
       IF(nstart(i).LE.nstop(i)) icoul=i
    ENDDO
    IF(nstart(myid).LE.nstop(myid)) THEN
       icolor1=1
       coul=.true.
    ENDIF

    call MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor1, myid, comm_coul,ierr)

    if (nstart(myid).LE.nstop(myid)) then
       CALL MPI_COMM_SIZE( comm_coul, npe_coul, ierr)
       CALL MPI_COMM_RANK( comm_coul, myid_coul, ierr)
    else
       npe_coul=0
       myid_coul=myid
    endif

!!$    icolor2=MPI_UNDEFINED
!!$    IF(myid.LT.nat) THEN
!!$       DO i=
!!$       nstart=myid*modu+1
!!$       nstop=MIN(nat,(myid+1)*modu)
!!$       icolor1=1
!!$       IF(nstop.EQ.nat) log_write=.TRUE.
!!$       IF(nstart.gt.nstop) THEN
!!$          icolor1=MPI_UNDEFINED
!!$          coul=.false.
!!$       ENDIF
!!$!       icolor2=1
!!$    endif
!!$    elseif(npe.ge.(2*nat)) then
!!$       icoul=nat-1
!!$       ixc=2*nat-1
!!$!       write(6,*) 'icoul,ixc',icoul,ixc
!!$       if(myid.lt.nat) then
!!$          modu=(nat-1)/(npe-nat)+1
!!$          nstart=myid*modu+1
!!$          nstop=MIN(nat,(myid+1)*modu)
!!$          icolor1=1
!!$       elseif(myid.ge.nat.and.myid.lt.(2*nat)) then
!!$          modu=(nat-1)/(npe-nat)+1
!!$          nstart=(myid-nat)*modu+1
!!$          nstop=MIN(nat,((myid-nat)+1)*modu)
!!$          icolor2=1
!!$      endif
!!$    endif
!    call MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor2, myid, comm_xc,ierr)
#else
    !BUG: this was missing: the argument should be declared even if Parallel is not defined!
    integer nat
#endif
  END SUBROUTINE begend
END MODULE parallel

MODULE ams
  REAL*8          :: atom_mass(118)

 CONTAINS
  SUBROUTINE init_ams
  REAL*8          :: a_m(118)
     DATA a_m /1.,4.,6.9,9.,10.8,12.,14.,16.,19.,20.2, &
          23.,24.3,27.,28.1,31.,32.,35.4,40.,39.1,40.,45., &
          47.9,50.9,52.,54.9,55.8,58.9,58.7,63.5,65.4,69.7, &
          72.6,74.9,79.,79.9,83.8,85.5,87.6,88.9,91.2,92.9, &
          95.9,98.,101.1,102.9,106.4,107.9,112.4,114.8, &
          118.7,121.8,127.6,126.9,131.3,132.9,137.3,138.9,140.1, &
          140.9,144.2,145.,150.4,152.,157.3,158.9,162.5, &
          164.9,167.3,168.9,173.,175.,178.5,180.9,183.8,186.2, &
          190.2,192.2,195.1,197.,200.6,204.4,207.2,209.,209., &
          210.,222.,223.,226.,227.,232.,231.,238.,237.,244.,243., &
          247.,247.,251.,252.,257.,258.,259.,262.,261.,262.,266.,&
          264.,269.,278.,281.,280.,285.,286.,289.,289.,293.,294.,294./     
     atom_mass(1:118) = a_m(1:118) 
   END SUBROUTINE init_ams
 END MODULE ams

module efeld
!!!  real*8  :: refeld
!!!  integer :: iefeld,jefeld(100)
  real*8  :: refeld, wefeld,fefeldpw(800)
  integer :: iefeld,jefeld(2000),kefeld,kebig
end module efeld

MODULE struct
  LOGICAL     :: ortho,inversion
  INTEGER     :: nat,nsym,ndif
  CHARACTER*4  :: lattic,irel,cform
  CHARACTER*80 :: title
  CHARACTER*10,ALLOCATABLE :: aname(:)
  REAL*8      :: br1(3,3),aa,bb,cc,br2(3,3),avec(3,3)
  REAL*8      :: a(3),alpha(3),pia(3),vol
  INTEGER,ALLOCATABLE :: iatnr(:),mult(:),jri(:),isplit(:)      !NATO
  REAL*8,ALLOCATABLE  :: r0(:),dx(:)                            !NATO   
  REAL*8,ALLOCATABLE  :: rotloc(:,:,:),rmt(:),zz(:),vsph(:)     !NATO   
  REAL*8,ALLOCATABLE  :: rotij(:,:,:),tauij(:,:)                !NDIF
  REAL*8,POINTER      :: pos(:,:)
END MODULE struct

MODULE symetr
  INTEGER   :: iord
  INTEGER,ALLOCATABLE :: iz(:,:,:),inum(:),kkk(:,:)                 !NSYM
  INTEGER,ALLOCATABLE :: kzz(:,:),kzzr2v(:,:),kzznonloc(:,:),inst(:)!NVAW
  REAL*8,ALLOCATABLE  :: tau(:,:)                                   !NSYM
END MODULE symetr

MODULE densit
  integer :: lmmaxx
  REAL*8,ALLOCATABLE :: clm(:,:,:),clmsp(:,:,:,:)
  REAL*8,ALLOCATABLE :: rhopw0(:,:,:),rhopw0sp(:,:,:,:),rhopw0r(:,:)
  COMPLEX*16,ALLOCATABLE :: rhok(:),rhoksp(:,:)
END MODULE densit

!!$module work
!!$  complex*16,allocatable :: YKA(:)
!!$  real*8,allocatable     :: V(:,:,:)
!!$  real*8,allocatable     :: r(:),vcoul2(:,:),value(:)
!!$end module work

module work
  complex*16,allocatable :: YKA(:)
  real*8,allocatable     :: V(:,:,:),exc(:,:,:,:),vxc(:,:,:,:),vxcold(:,:,:,:),vxcnonloc(:,:,:,:)
  real*8,allocatable     :: vtau(:,:,:,:),vlapl(:,:,:,:) ! non-multiplicative MGGA potential part
  real*8,allocatable     :: r(:),vcoul2(:,:),value(:)
! slater_tau_lm(r,lm,jatom,spin)
! slater form of kinetic energy density: |grad psi|**2
  real*8,allocatable     :: slater_tau_lm(:,:,:,:)
end module work

module vresp1 
  real*8,allocatable     :: vclmsp(:,:,:,:)
end module vresp1

module xcparam
  character*50 EX_SWITCH,EC_SWITCH,VX_SWITCH,VC_SWITCH,VDISP_SWITCH,XC_SWITCH
  character*50 X_SWITCH_LIBXC, C_SWITCH_LIBXC
  integer :: id_func_x, id_func_c
  real*8 :: xc1, xc2, xc3, xc4
!$omp threadprivate(xc1,xc2,xc3,xc4)
  real*8 :: grr, xca, xcb, xcp, xcconst, ycconst, vxmix, r_th
!$omp threadprivate(xcconst)
  real*8 :: alphax, slamb, rhocutvdw
  real*8, allocatable :: stagfield(:)
  logical :: screening,VRESP_SWITCH,GGA_SWITCH,MGGA_SWITCH,MGGA_LAPL_SWITCH,&
       MGGA_TAU_SWITCH,AEQ_SWITCH,GRR_SWITCH,read_nlvdw,READ_TAU_SWITCH,lmbj_switch
  logical :: read_vhalf, dostagfield
  character*20 VHALF_SWITCH
! various print counters
  integer(kind=1) :: icount = 0
  integer(kind=1) :: iint = 0
  integer(kind=1) :: isphere = 0
  integer :: windex = 1
! BUG: should be removed (index_vxclm2 was never declared!)
!!!$omp threadprivate(index_vxclm2)
end module xcparam

module norm_kub
  
  real*8   ::  c_kub(0:10,0:10)

contains
  subroutine init_c_kub
    
    c_kub=0.0d0
    c_kub(0,0)=1.d0
    c_kub(3,2)=1.d0
    c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
    c_kub(4,4)=.5d0*SQRT(5.d0/3.d0)
    c_kub(6,0)=.5d0*SQRT(.5d0)
    c_kub(6,2)=.25d0*SQRT(11.d0)
    c_kub(6,4)=-.5d0*SQRT(7.d0/2.d0)
    c_kub(6,6)=-.25d0*SQRT(5.d0)
    c_kub(7,2)=.5d0*SQRT(13.d0/6.d0)
    c_kub(7,6)=.5d0*SQRT(11.d0/6.d0)
    c_kub(8,0)=.125d0*SQRT(33.d0)
    c_kub(8,4)=.25d0*SQRT(7.d0/3.d0)
    c_kub(8,8)=.125d0*SQRT(65.d0/3.d0)
    c_kub(9,2)=.25d0*SQRT(3.d0)
    c_kub(9,4)=.5d0*SQRT(17.d0/6.d0)
    c_kub(9,6)=-.25d0*SQRT(13.d0)
    c_kub(9,8)=-.5d0*SQRT(7.d0/6.d0)
    c_kub(10,0)=.125d0*SQRT(65.D0/6.D0)
    c_kub(10,2)=.125d0*SQRT(247.D0/6.D0)
    c_kub(10,4)=-.25d0*SQRT(11.D0/2.D0)
    c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
    c_kub(10,8)=-.125d0*SQRT(187.D0/6.D0)
    c_kub(10,10)=-.0625d0*SQRT(85.d0)

  end subroutine init_c_kub
end module norm_kub


module str_switch
  CHARACTER*5  :: SW_STRESS
  logical  :: SW_TOT_STR
  logical  :: STRESS  
end module str_switch

module stress_time
   real*8   :: time1_st, time2_st
   real*8   :: time_electro
end module stress_time

MODULE PW_mult
  COMPLEX*16,ALLOCATABLE     :: DECOUL_POT_PW(:,:), CFFT_STR(:,:,:,:), CFFT8(:,:,:)
  COMPLEX*16                 :: cabt(3,-1:1)
  REAL*8,ALLOCATABLE         :: SGGP(:,:), SGGP_GGA(:,:)
  COMPLEX*16, ALLOCATABLE    :: CFFT_STR_GGA(:,:,:,:)
  COMPLEX*16,ALLOCATABLE     ::  meta_fft(:,:,:,:)

END MODULE PW_mult

module gga_str_kb
!  real*8,allocatable      :: sigma_str(:,:,:,:)
  complex*16, allocatable :: gga_int(:,:)
!  real*8                  :: sig_str_u, sig_str_d, sig_str_ud  
  real*8                  :: GGA_INTER(1:9,1:2) 
  real*8                  :: gga_fin_str(1:9)
end module gga_str_kb
