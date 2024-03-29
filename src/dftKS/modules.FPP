module param
  INTEGER            LMAX,  LMMX, LOMAX, restrict_output
  INTEGER            NATO, NDIF,  NGAU, NKPTSTART, NMATMAX,  NRAD
  INTEGER            NSLMAX, NSYM, NUME,NVEC1, NVEC2, NVEC3
  INTEGER            NWAV,NMATIT,NUMEIT,HB,NMATBL,nloat
  INTEGER            BLOCKSIZE, nume2
  PARAMETER          (LMAX=   13)
  PARAMETER          (LMMX=  120)
  PARAMETER          (LOMAX=   3) 
  PARAMETER          (NGAU= 2350)
  PARAMETER          (NKPTSTART= 100)
  PARAMETER          (NMATMAX=   13000)   ! largest possible matrix size ! should change !
  PARAMETER          (NRAD=  881)
  PARAMETER          (NSLMAX=  5)
  PARAMETER          (NSYM=   48)
  PARAMETER          (NUME=   3000)   
  PARAMETER          (NVEC1=  35)
  PARAMETER          (NVEC2=  35)
  PARAMETER          (NVEC3=   95)
  !      PARAMETER          (nloat= 60)      !  should be 3 for only one LO
  PARAMETER          (restrict_output= 9999) ! 1 for mpi with less output-files
  PARAMETER          (BLOCKSIZE=128)
  CHARACTER*180     :: filename_V_sph
  CHARACTER*180     :: filename_V_nsh
  CHARACTER*180     :: filename_V_vns
  CHARACTER*180     :: filename_vector, filename_energy
end module param

module lapw_timer
  !  64 Timers are available in ScaLapack
  double precision :: cpu_time(64)
  double precision :: old_cpu_time(64)
  double precision :: wall_time(64)
  double precision :: old_wall_time(64)
  integer, parameter :: time_total=1, time_hamilt=10, time_hns=11, time_diag=12, time_horb=13,  &
       time_coor=15, time_setwar=17, time_pkp=18, time_cholesky=20, time_transform=21, time_sep=22, time_backtransform=23,&
       time_pwpw=25, time_pwlo=26,  time_lolo=27,  time_hns30=28, time_iter=30,time_H1=31, time_HxF=32, time_F=33,&
       time_Ffac=34, time_Fleq=35, time_SxF=36, time_projH=37, time_projH_seq=38, time_backproject=39, time_albl=40,&
       time_loop260=41, time_phase=42, time_legendre=43, time_step_function=44, time_h=45,time_us=46, time_loop210=47,&
       time_loop230=48, time_loop240=49, time_overlap=50, time_lo=51, time_g1=55, time_atpar=56, time_distiouter=57, time_horb_setup=58, &
       time_debug1=2, time_debug2=3
  logical    usecputime, usewalltime
contains
  subroutine init_all_timer
    cpu_time(:)=0.0D0
    wall_time(:)=0.0D0
  end subroutine init_all_timer
  
  subroutine init_timer(n)
    IMPLICIT NONE
    integer :: n
    cpu_time(n) = 0.0D0
    wall_time(n) = 0.0D0
  end subroutine init_timer
  
  subroutine start_timer(n)
    integer :: n
    call CPUTIM(old_cpu_time(n))
    cpu_time(n)=0.0D0
    call WALLTIM(old_wall_time(n))
    wall_time(n)=0.0D0
  end subroutine start_timer

  subroutine stop_timer(n)
    integer :: n
    double precision :: current_time
    call CPUTIM(current_time)
    cpu_time(n) = cpu_time(n) + current_time - old_cpu_time(n)
    call WALLTIM(current_time)
    wall_time(n) = wall_time(n) + current_time - old_wall_time(n)
  end subroutine stop_timer
         
  double precision function read_cpu_time(n)
    integer :: n
    read_cpu_time=cpu_time(n)
    return
  end function read_cpu_time

  double precision function read_wall_time(n)
    integer :: n
    read_wall_time=wall_time(n)
  end function read_wall_time
end module lapw_timer
          
!--------------------------------------------------------------

module reallocate
  interface doreallocate
     module procedure doreallocate_r8_d1
     module procedure doreallocate_r8_d2
     module procedure doreallocate_r8_d3
     module procedure doreallocate_i4_d1
     module procedure doreallocate_i4_d2
     module procedure doreallocate_warp_c16_d2x3
     module procedure doreallocate_c10_d1
     module procedure hugo     !   ;)
  end interface doreallocate
contains
  
  !     leider sind mehrere subroutines notwendig fuer verschiedene Typen
  subroutine doreallocate_r8_d1(tf, newdimension)
    implicit none
    real*8, pointer :: hilfsfeld(:), tf(:)
    integer min1, newdimension
    allocate(hilfsfeld(newdimension))
    hilfsfeld=0.0D0
    !     nur 1 mal kopieren reicht
    !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
    min1=min(newdimension,size(tf,1))
    hilfsfeld(1:min1)=tf(1:min1)
    deallocate(tf)
    !     der Zeiger wird nur auf das neue Feld umgebogen, nicht neu alloziert
    tf=>hilfsfeld
  end subroutine doreallocate_r8_d1
  
  subroutine doreallocate_r8_d2(tf, newdimension1, newdimension2)
    implicit none
    real*8, pointer :: hilfsfeld(:,:), tf(:,:)
    integer min1, min2, newdimension1, newdimension2
    allocate(hilfsfeld(newdimension1,newdimension2))
    hilfsfeld=0.0D0
    min1=min(newdimension1,size(tf,1))
    min2=min(newdimension2,size(tf,2))
    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
    deallocate(tf)
    tf=>hilfsfeld
  end subroutine doreallocate_r8_d2
  
  subroutine doreallocate_r8_d3(tf, newdimension1, newdimension2, newdimension3)
    implicit none
    real*8, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
    integer min1, min2, min3, newdimension1, newdimension2, newdimension3
    allocate(hilfsfeld(newdimension1,newdimension2,newdimension3))
    hilfsfeld=0.0D0
    min1=min(newdimension1,size(tf,1))
    min2=min(newdimension2,size(tf,2))
    min3=min(newdimension3,size(tf,3))
    hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
    deallocate(tf)
    tf=>hilfsfeld
  end subroutine doreallocate_r8_d3
  
  subroutine doreallocate_i4_d1(tf, newdimension)
    implicit none
    integer*4, pointer :: hilfsfeld(:), tf(:)
    integer min1, newdimension
    allocate(hilfsfeld(newdimension))
    hilfsfeld=0
    min1=min(newdimension,size(tf,1))
    hilfsfeld(1:min1)=tf(1:min1)
    deallocate(tf)
    tf=>hilfsfeld
  end subroutine doreallocate_i4_d1
  
  subroutine doreallocate_i4_d2(tf, newdimension1, newdimension2)
    implicit none
    integer*4, pointer :: hilfsfeld(:,:), tf(:,:)
    integer min1, min2, newdimension1, newdimension2
    allocate(hilfsfeld(newdimension1,newdimension2))
    hilfsfeld=0
    min1=min(newdimension1,size(tf,1))
    min2=min(newdimension2,size(tf,2))
    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
    deallocate(tf)
    tf=>hilfsfeld
  end subroutine doreallocate_i4_d2
  
  subroutine doreallocate_warp_c16_d2x3(tf, n1a, n1b, n2a, n2b, n3a, n3b)
    implicit none
    complex*16, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
    integer n1a, n1b, n2a, n2b, n3a, n3b
    allocate( hilfsfeld(-n1b:n1b, -n2b:n2b, -n3b:n3b) )
    hilfsfeld(-n1b:n1b, -n2b:n2b, -n3b:n3b)=(0.0D0,0.0D0)
    hilfsfeld(-n1a:n1a,-n2a:n2a,-n3a:n3a)=tf(-n1a:n1a,-n2a:n2a,-n3a:n3a)
    deallocate(tf)
    tf=>hilfsfeld
  end subroutine doreallocate_warp_c16_d2x3
  
  subroutine doreallocate_c10_d1(tf, newdimension1)
    implicit none
    character*10, pointer :: hilfsfeld(:), tf(:)
    integer min1, newdimension1
    allocate(hilfsfeld(newdimension1))
    min1=min(newdimension1,size(tf,1))
    hilfsfeld(1:min1)=tf(1:min1)
    deallocate(tf)
    tf=>hilfsfeld
  end subroutine doreallocate_c10_d1
  
  !     Es gibt auch Methoden, um das Programm unleserlich zu machen :-)
  !     das sollten wir besser vermeiden ;-)
  subroutine hugo
    write(*,*) " Hier ist Hugo"
  end subroutine hugo
end module reallocate

!--------------------------------------------------------------
module matrices
  use param
  integer hsrows
  !        HS - Hamilton matrix (hermitian matrix, stored below diagonal in conjugated form)
  !           - Overlap matrix  (hermitian matrix, stored above diagonal)
  !        HSDIAG - Diagonal elements of Overlap matrix
  !               - several times swapped with 
  !                 diagonal elements of Hamilton matrix
#ifdef _REAL  
  DOUBLE PRECISION, allocatable ::  HS(:, :)
  DOUBLE PRECISION, allocatable ::  HSDIAG(:)
  DOUBLE PRECISION, allocatable ::  Z(:,:)
#else  
  DOUBLE COMPLEX, allocatable :: HS(:, :)
  DOUBLE COMPLEX, allocatable :: HSDIAG(:)
  DOUBLE COMPLEX, allocatable :: Z(:, :)
#endif  
  DOUBLE PRECISION, allocatable ::  EIGVAL(:)
  DOUBLE PRECISION, pointer ::  Kn(:,:)
  DOUBLE PRECISION, pointer ::  RK(:)
  INTEGER, pointer :: KZZ(:,:)
CONTAINS
  SUBROUTINE INIT_MATRICES(N, NUME2_,Qprint)
    use reallocate, only: doreallocate
    use param, only : nume2
    IMPLICIT NONE        
    INTEGER, intent(in)  :: N, NUME2_
    LOGICAL, intent(in)  :: Qprint
!    nmatmax1=nmatmax
!#ifndef _REAL    
!    nmatmax1=int(nmatmax1/sqrt(2.d0))
!#endif    
    if( n.eq.0 .and. nume2_.eq.0) then
       !allocate( RK(NMATMAX1+1) )
       !allocate( Kn(3,NMATMAX1+1) )
       !allocate( KZZ(3,NMATMAX1+1) )
       !RK=0.0D0
       !Kn(:,:)=0.d0
       !KZZ=0
       WRITE(6,*) 'ERROR : We should not come here !'
       STOP 'ERROR LAPW1'
    else
       if (Qprint) write(6,*) 'Matrix size ', N
       HSROWS = N
       nume2 = nume2_
       call doreallocate( RK, N+1 )
       call doreallocate( Kn, 3,N+1 )
       call doreallocate( KZZ,3,N+1 )
       allocate( HS(N,N) )
       allocate( HSDIAG(N) )
       allocate( Z(N,NUME2) )
       if (Qprint) then
#ifdef _REAL       
          write(6,'(a20,f12.1,a)') 'allocate HS',HSROWS**2*8.0/(1024.0*1024.0),' MB '
          write(6,'(a20,f12.1,a)') 'allocate Z ',HSROWS*NUME2*8.0/(1024.0*1024.0),' MB '
#else       
          write(6,'(a20,f12.1,a)') 'allocate HS',HSROWS**2*16.0/(1024.0*1024.0),' MB '
          write(6,'(a20,f12.1,a)') 'allocate Z ',HSROWS*NUME2*16.0/(1024.0*1024.0),' MB '
#endif
       endif
       HS=0.0D0
       HSDIAG=0.0D0
       Z=0.0D0
       allocate( EIGVAL(N) )
       EIGVAL=0.0D0
    end if

  END SUBROUTINE INIT_MATRICES
        
  SUBROUTINE END_MATRICES
    deallocate( HS, HSDIAG )
    deallocate( KZZ )
    deallocate( Z )
    deallocate( EIGVAL )
  END SUBROUTINE END_MATRICES

  SUBROUTINE END_MATRICES_PART
    deallocate( Rk, Kn)
  END SUBROUTINE END_MATRICES_PART

end module matrices

module albl
  !  AL(n,l,i) - coefficient al(kn) for atom i (including equivalent atoms)
  !  BL(n,l,i) - coefficient bl(kn) for atom i (including equivalent atoms)
  DOUBLE PRECISION, allocatable ::   AL_r(:,:)
  DOUBLE PRECISION, allocatable ::   BL_r(:,:)
  DOUBLE PRECISION, allocatable ::   AL_c(:,:)
  DOUBLE PRECISION, allocatable ::   BL_c(:,:)
  DOUBLE PRECISION, allocatable ::   AL(:,:,:)
  DOUBLE PRECISION, allocatable ::   BL(:,:,:)
contains
  subroutine init_albl(lmax,SUB)
    use matrices, only : HSrows
    integer  lmax
    character*3 sub
    if (sub.eq.'HAM') then
       allocate( AL_r(HSROWS, 0:LMAX) )
       allocate( BL_r(HSROWS, 0:LMAX) )
       allocate( AL_c(HSROWS, 0:LMAX) )
       allocate( BL_c(HSROWS, 0:LMAX) )
    else
       allocate( AL(HSROWS, 0:LMAX,1) )
       allocate( BL(HSROWS, 0:LMAX,1) )
       AL=0.0D0
       BL=0.0D0
    endif
  end subroutine init_albl

  subroutine end_albl(sub)
    character*3 sub
    if (sub.eq.'HAM') then
       deallocate( AL_r, BL_r)
       deallocate( AL_c, BL_c)
    else
       deallocate( AL, BL)
    endif
  end subroutine end_albl
  
end module albl
!--------------------------------------------------------------

module loabc
  !
  !        ALO(l,i)    - local orbital coefficient alm for atom i
  !        BLO(l,i)    - local orbital coefficient blm for atom i
  !        CLO(l,i)    - local orbital coefficient clm for atom i
  !        DPELO(l,i)  - derivative of the energy derivative ul_dot(r,E)
  !                      of the local orbital radial-function ul(r,E) with
  !                      respect to r
  !                      (evaluated at r = muffin tin radius of atom i)
  !        DPLO(l,i)   - derivative of the local orbital radial-function
  !                      ul(r,E) with respect to radius r
  !                      (evaluated at r = muffin tin radius of atom i)
  !        ELO(l,i)    - local orbital expansion energy El for atom i
  !        PEILO(l,i)  - the norm of the energy derivative ul_dot(r,E)
  !                      of the local orbital radial function ul(r,E1)
  !                      integrated over the muffin tin sphere
  !        PELO(l,i)   - ul_dot(r,E)
  !                      (evaluated at r = muffin tin radius of atom i)
  !        PLO(l,i)    - local orbital radial-function ul(r,E)
  !                      (evaluated at r = muffin tin radius of atom i)
  !        PI12LO(l,i) - the product of the local orbital radial functions
  !                      ul(r,E1) and ul(r,E2) integrated over the muffin
  !                      tin sphere (for atom i)
  !        PE12LO(l,i) - the product of the derivative of ul(r,E1) with
  !                      respect to r and ul(r,E2) integrated over the
  !                      muffin tin sphere (for atom i)
  integer nloat_new
  DOUBLE PRECISION, allocatable ::   ALO(:,:,:)
  DOUBLE PRECISION, allocatable ::   BLO(:,:,:)
  DOUBLE PRECISION, allocatable ::   CLO(:,:,:)
  DOUBLE PRECISION, allocatable ::   ELO(:,:,:)
  DOUBLE PRECISION, allocatable ::   PLO(:,:,:)
  DOUBLE PRECISION, allocatable ::   DPLO(:,:,:)
  DOUBLE PRECISION, allocatable ::   PELO(:,:,:)
  DOUBLE PRECISION, allocatable ::   DPELO(:,:,:)
  DOUBLE PRECISION, allocatable ::   PEILO(:,:,:)
  DOUBLE PRECISION, allocatable ::   PI12LO(:,:,:)
  DOUBLE PRECISION, allocatable ::   PE12LO(:,:,:)
  DOUBLE PRECISION, allocatable ::   PILOLO(:,:,:,:)
contains
  SUBROUTINE init_loabc(LOMAX, NLOAT, NATO)
    IMPLICIT NONE
    integer LOMAX, NLOAT, NATO
    allocate( ALO(0:LOMAX,nloat,NATO) )
    allocate( BLO(0:LOMAX,nloat,NATO) )
    allocate( CLO(0:LOMAX,nloat,NATO) )
    allocate( ELO(0:LOMAX,nloat,NATO) )
    allocate( PLO(0:LOMAX,nloat,NATO) )
    allocate( DPLO(0:LOMAX,nloat,NATO) )
    allocate( PELO(0:LOMAX,nloat,NATO) )
    allocate( DPELO(0:LOMAX,nloat,NATO) )
    allocate( PEILO(0:LOMAX,nloat,NATO) )
    allocate( PI12LO(0:LOMAX,nloat,NATO) )
    allocate( PE12LO(0:LOMAX,nloat,NATO) )
    allocate( PILOLO(0:LOMAX,nloat,nloat,NATO) )
    ALO=0.0D0
    BLO=0.0D0
    CLO=0.0D0
    ELO=0.0D0
    PLO=0.0D0
    DPLO=0.0D0
    PELO=0.0D0
    DPELO=0.0D0
    PEILO=0.0D0
    PI12LO=0.0D0
    PE12LO=0.0D0
    PILOLO=0.0D0
  END SUBROUTINE init_loabc
end module loabc

!--------------------------------------------------------------

module radial
  DOUBLE PRECISION, allocatable ::   A(:), A1(:,:)
  DOUBLE PRECISION, allocatable ::   A1LO(:,:,:)
  DOUBLE PRECISION, allocatable ::   AE(:), AE1(:,:)
  DOUBLE PRECISION, allocatable ::   AE1LO(:,:,:)
  DOUBLE PRECISION, allocatable ::   AP(:), B(:), B1(:,:)
  DOUBLE PRECISION, allocatable ::   B1LO(:,:,:)
  DOUBLE PRECISION, allocatable ::   BE(:), BE1(:,:)
  DOUBLE PRECISION, allocatable ::   BE1LO(:,:,:)
  DOUBLE PRECISION, allocatable ::   BP(:), VLM(:,:)
contains
  SUBROUTINE init_radial(NRAD, LOMAX, NSLMAX, NLOAT, LMMX)
    IMPLICIT NONE
    integer NRAD, NSLMAX, LOMAX, NLOAT, LMMX
    allocate( A(NRAD) )
    allocate( A1(NRAD,NSLMAX) )
    allocate( A1LO(NRAD,nloat,0:LOMAX) )
    allocate( AE(NRAD) )
    allocate( AE1(NRAD,NSLMAX) )
    allocate( AE1LO(NRAD,nloat,0:LOMAX) )
    allocate( AP(NRAD) )
    allocate( B(NRAD) )
    allocate( B1(NRAD,NSLMAX) )
    allocate( B1LO(NRAD,nloat,0:LOMAX) )
    allocate( BE(NRAD) )
    allocate( BE1(NRAD,NSLMAX) )
    allocate( BE1LO(NRAD,nloat,0:LOMAX) )
    allocate( BP(NRAD) )
    allocate( VLM(NRAD,LMMX) )
    A=0.0D0
    A1=0.0D0
    A1LO=0.0D0
    AE=0.0D0
    AE1=0.0D0
    AP=0.0D0
    B=0.0D0
    B1=0.0D0
    B1LO=0.0D0
    BE=0.0D0
    BE1=0.0D0
    BE1LO=0.0D0
    BP=0.0D0
    VLM=0.0D0
  END SUBROUTINE init_radial
  SUBROUTINE end_radial
    deallocate(A, A1, A1LO, AE, AE1, AE1LO, AP, B, B1, B1LO, BE, BE1, BE1LO, BP, VLM)
  END SUBROUTINE end_radial
end module radial

!--------------------------------------------------------------

module orb
  !   NMOD=1 LDA+U, NMOD=2 OP, NMOD=3 Bext calculation
  !   NSP=1 spin up, NSP=-1 spin down
  !   NATORB number of atoms for which orbital potential is included
  !   Iatom are the indexes of atom types selected for orb.pot.
  !   NLORB number of L's of given atom for which OP is considered
  !   LORB L's of given atom
  !   Bext is energy of ext. magn field muB*Bext in Hartree
  !   VORB is orbital potential matrix created by ORB package
  !       COMPLEX*16 VORB(NATO,0:3,-3:3,-3:3)
  !       DOUBLE PRECISION BEXT
  !       INTEGER NMOD, NSP, NATORB, IATOM(NATO), NLORB(NATO), LORB(3,NATO)
  !       COMMON /ORB/  Bext,VORB, NMOD,NSP,Natorb,IATOM,NLORB,LORB
  !       SAVE   /ORB/
  DOUBLE PRECISION            :: BEXT
  DOUBLE COMPLEX, allocatable :: VORB(:,:,:,:)
  INTEGER                     :: NMOD, NSP, NATORB
  INTEGER, allocatable        :: IATOM(:), NLORB(:), LORB(:,:)
contains
  SUBROUTINE init_orb(NAT)
    IMPLICIT NONE
    integer NAT
    allocate( VORB(NAT,0:3,-3:3,-3:3) )
    allocate( IATOM(NAT) )
    allocate( NLORB(NAT) )
    allocate( LORB(3,NAT) )
    VORB=(0.0D0,0.0D0)
    IATOM=0
    NLORB=0
    LORB=0
  END SUBROUTINE init_orb
end module orb
!--------------------------------------------------------------

module lolog
  !        LOOR(l,i) - selects which local orbitals should be used for
  !                    atom i
  !        NLO       - total number of used local orbitals
  !
  !      INTEGER            NLO
  !      LOGICAL            LOOR(0:LOMAX,NATO),lapw(0:lmax-1,nato)
  !      COMMON  /LOLOG/    NLO, LOOR, lapw, ilo(0:lomax,nato)
  !      SAVE    /LOLOG/
  INTEGER              :: NLO
  LOGICAL, allocatable :: LOOR(:,:), LAPW(:,:)
  INTEGER, allocatable :: ILO(:,:), nlov(:)
contains
  SUBROUTINE init_lolog(LOMAX, NATO, LMAX)
    IMPLICIT NONE
    integer LOMAX, NATO, LMAX
    allocate( LOOR(0:LOMAX, NATO) )
    allocate( LAPW(0:LMAX-1, NATO) )
    allocate( ILO(0:LOMAX, NATO) )
    LOOR=.true.
    LAPW=.true.
    ILO=0
  END SUBROUTINE init_lolog
end module lolog

!--------------------------------------------------------------
module atspdt
  !        DP(l,i)      - derivative of the radial-function ul(r,E) with
  !                       respect to radius r
  !                       (evaluated at r = muffin tin radius of atom i)
  !        DPE(l,i)     - derivative of the energy derivative ul_dot(r,E)
  !                       of the radial-function ul(r,E) with respect to
  !                       radius r
  !                       (evaluated at r = muffin tin radius of atom i)
  !        GFAC(:,i)    - gaunt-factors for atom i
  !        LM(:,j,i)    - L,M combination for lattice harmonics j for
  !                       atom i
  !                       LM(1,j,i) ... L for j-th gaunt factor
  !                       LM(2,j,i) ... M for j-th gaunt factor
  !        LMMAX(i)     - number of L,M combinations for the
  !                       nonspherical contribution of atom i
  !        LQIND(i)     - number of gaunt factors for atom i
  !        LQNS(:,j,i)  - l,m,LM,l',m' combination corresponding to
  !                       the gaunt-factor GFAC(j,i) for atom i
  !                       LQNS(1,1:LQIND(i),i) ... l'+1 for atom i
  !                       LQNS(2,1:LQIND(i),i) ... L
  !                       LQNS(3,1:LQIND(i),i) ... l+1 for atom i
  !                       LQNS(4,1:LQIND(i),i) ... index of LM for atom i
  !                       LQNS(5,1:LQIND(i),i) ... m' for atom i
  !                       LQNS(6,1:LQIND(i),i) ... m for atom i
  !        E(l,i)       - expansion energy El for atom i
  !        P(l,i)       - radial-function ul(r,E)
  !                       (evaluated at r = muffin tin radius of atom i)
  !        PE(l,i)      - derivative of the radial-function ul(r,E) with
  !                       respect to energy E
  !                       (evaluated at r = muffin tin radius of atom i)
  !        PEI(l,i)     - norm of ul_dot(r,E) integrated over the muffin
  !                       tin sphere
  !        VNS1(l',L,l,i) - non-spherical matrix elements for atom i
  !                         integral(u_{l'}*V_{L}*u_{l} dr)
  !        VNS2(l',L,l,i) - non-spherical matrix elements for atom i
  !                         integral(u_{l'}_dot*V_{L}*u_{l}_dot dr)
  !        VNS3(l',L,l,i) - non-spherical matrix elements for atom i
  !                         integral(u_{l'}*V_{L}*u_{l}_dot dr)
  !                         VNS3(l,L,l') ...
  !                            ... integral(u_{l'}_dot*V_{L}*u_{l} dr)
  !
  !      INTEGER            LMMAX(NATO), LQIND(NATO)
  !      INTEGER            LM(2,LMMX,NATO), LQNS(6,NGAU,NATO)
  !      DOUBLE PRECISION   DP(LMAX,NATO), DPE(LMAX,NATO), E(LMAX,NATO)
  !      DOUBLE PRECISION   P(LMAX,NATO), PE(LMAX,NATO), PEI(LMAX,NATO)
  !      DOUBLE PRECISION   VNS1(NSLMAX,LMMX,NSLMAX,NATO)
  !      DOUBLE PRECISION   VNS2(NSLMAX,LMMX,NSLMAX,NATO)
  !      DOUBLE PRECISION   VNS3(NSLMAX,LMMX,NSLMAX,NATO)
  !      COMPLEX*16         GFAC(NGAU,NATO)
  INTEGER                       :: NL
  INTEGER, allocatable          :: LMMAX(:), LQIND(:), LM(:,:,:), LQNS(:,:,:)
  DOUBLE PRECISION, allocatable :: DP(:,:), DPE(:,:), E(:,:), P(:,:), PE(:,:), PEI(:,:)
  DOUBLE PRECISION, allocatable :: VNS1(:,:,:,:), VNS2(:,:,:,:), VNS3(:,:,:,:)
  DOUBLE COMPLEX, allocatable   :: GFAC(:,:)
contains
  SUBROUTINE init_atspdt(NATO, LMMX, NGAU, LMAX, NSLMAX)
    IMPLICIT NONE
    integer NATO, LMMX, NGAU, LMAX, NSLMAX
    allocate( LMMAX(NATO) )
    allocate( LQIND(NATO) )
    allocate( LM(2,LMMX,NATO) )
    allocate( LQNS(6,NGAU,NATO) )
    allocate( DP(LMAX,NATO) )
    allocate( DPE(LMAX,NATO) )
    allocate( E(LMAX,NATO) )
    allocate( P(LMAX,NATO) )
    allocate( PE(LMAX,NATO) )
    allocate( PEI(LMAX,NATO) )
    allocate( VNS1(NSLMAX,LMMX,NSLMAX,NATO) )
    allocate( VNS2(NSLMAX,LMMX,NSLMAX,NATO) )
    allocate( VNS3(NSLMAX,LMMX,NSLMAX,NATO) )
    allocate( GFAC(NGAU,NATO) )
    LMMAX=0
    LQIND=0
    LM=0
    LQNS=0
    DP=0.0D0
    DPE=0.0D0
    E=0.0D0
    P=0.0D0
    PE=0.0D0
    PEI=0.0D0
    VNS1=0.0D0
    VNS2=0.0D0
    VNS3=0.0D0
    GFAC=(0.0D0,0.0D0)
  END SUBROUTINE init_atspdt
end module atspdt

!--------------------------------------------------------------

module comc
  !        IPGR(j)  - point group of j-th k-point (only used for output
  !                   (documentation) purposes)
  !        KNAME(j) - name of j-th k-point (optional)
  CHARACTER*3, pointer :: IPGR(:)
  CHARACTER*10, pointer :: KNAME(:)
end module comc

!--------------------------------------------------------------

module comi
  !        LNSMAX - maximum considered azimuthal quantum number l of ul
  !                 for non-spherical hamilton matrix contributions
  !        NAT    - number of inequivalent atoms
  !        NBELW  - number of Eigenvalues below given energy window
  !        NE     - number of calculated Eigenvalues
  !        NT     - maximum considered azimuthal quantum number + 1
  !                 (l+1 of ul) for the spherical contribution
  !        NVAA   - number of plane waves
  INTEGER :: NVAA, NE, NBELW, NAT, NT, LNSMAX,nemax_it
end module comi

!--------------------------------------------------------------

module coml
  !        PRNTWF - if .TRUE. print wave functions to output-file
  !        SPRSWF - if .TRUE. suppress wave function calculation
  !        WFTAPE - if .TRUE. write wave functions to vector-file
  !        REL    - if .TRUE. perform relativistic calculations
  !        ITER   - if .TRUE. perform iterative diagonalization
  !        NOHNS  - if .TRUE. supress nonspherical matrix elements
  LOGICAL :: SPRSWF,  PRNTWF,  WFTAPE,  REL, ITER,NOHNS
end module coml

!--------------------------------------------------------------

module comr
  !        EL        - lower boundary of energy window for Eigenvalues
  !        EU        - upper boundary of energy window for Eigenvalues
  !        EIGVAL(j) - calculated Eigenvalues
  !        RKM       - Rmt*Kmax determines matrix size (convergence)
  !                    Rmt .... muffin tin radius
  !                    Kmax ... plane wave cut-off
  !        WEIGHT(j) - weight of the j-th k-vector (order of group of k)
  !
  DOUBLE PRECISION :: RKM, ELow, EUp,sproj_limit, eshift_iter
  DOUBLE PRECISION, pointer :: WEIGHT(:)
  logical readHinv
end module comr

!--------------------------------------------------------------

module kpts
  ! k-points in the first Brillouine zone
  REAL*8, allocatable :: K3(:,:)
end module kpts

!--------------------------------------------------------------

module loint
  !        VNS1LO(l',L,l) - local orbital non-spherical matrix elements
  !                         integral(u_{l'}*V_{L}*u_{l} dr)
  !        VNS2LO(l',L,l) - local orbital non-spherical matrix elements
  !                         integral(u_{l'}_dot*V_{L}*u_{l}_dot dr)
  !        VNS3LO(l',L,l) - local orbital non-spherical matrix elements
  !                         integral(u_{l'}*V_{L}*u_{l}_dot dr)
  !                         VNS3(l,L,l') ...
  !                           ... integral(u_{l'}_dot*V_{L}*u_{l} dr)
  !
  !      DOUBLE PRECISION   VNS1LO(LOMAX+1,LMMX,NSLMAX,NATO)
  !      DOUBLE PRECISION   VNS2LO(LOMAX+1,LMMX,NSLMAX,NATO)
  !      DOUBLE PRECISION   VNS3LO(LOMAX+1,LMMX,NSLMAX,NATO)
  !      COMMON  /LOINT/    VNS1LO, VNS2LO, VNS3LO
  !      SAVE    /LOINT/
  DOUBLE PRECISION, allocatable :: VNS1LO(:,:,:,:,:), VNS2LO(:,:,:,:,:), VNS3LO(:,:,:,:,:,:)
contains
  SUBROUTINE init_loint(LOMAX, LMMX, NSLMAX, NATO,nloat)
    IMPLICIT NONE
    integer LOMAX, LMMX, NSLMAX, NATO,nloat
    allocate( VNS1LO(LOMAX+1, LMMX, NSLMAX, nloat, NATO) )
    allocate( VNS2LO(LOMAX+1, LMMX, NSLMAX, nloat, NATO) )
    allocate( VNS3LO(LOMAX+1, LMMX, NSLMAX, nloat, nloat,NATO) )
    VNS1LO=0.0D0
    VNS2LO=0.0D0
    VNS3LO=0.0D0
  END SUBROUTINE init_loint
end module loint

!--------------------------------------------------------------

module lstapw
  !        NV    - total number of considered reciprocal lattice vectors
  !                (number of plane waves)
  !        RKM       - Rmt*Kmax determines matrix size (input)
  DOUBLE PRECISION :: RKMT
  INTEGER          :: NV
end module lstapw

!--------------------------------------------------------------
module out
  INTEGER     :: NKK
#ifdef _REAL  
  REAL*8, allocatable     :: WARP(:,:,:)
#else  
  COMPLEX*16, allocatable :: WARP(:,:,:)
#endif  
  INTEGER                 :: kmax(3)
contains
  SUBROUTINE end_out
    if( allocated(WARP) ) then
       deallocate( WARP )
    end if
  END SUBROUTINE end_out
end module out

!--------------------------------------------------------------
module potnlc
  !        VR(j,i) - spherical part (l=0, m=0) of the total potential r*V
  !                  at mesh point j for atom i
  DOUBLE PRECISION, allocatable :: VR(:,:)
end module potnlc

!--------------------------------------------------------------
module cut1
  LOGICAL :: CUT
end module cut1

!--------------------------------------------------------------
module totpot
  INTEGER, allocatable        :: KPxyz(:,:)
#ifdef _REAL
  REAL*8,    allocatable :: POTK(:)
#else  
  COMPLEX*16, allocatable :: POTK(:)
#endif  
contains
end module totpot

!--------------------------------------------------------------
!module kpthold
!  integer number_kpt
!end module kpthold
      
!--------------------------------------------------------------
module nmr
  logical  writeham,  dokorig
  real*8    korig(3)
end module nmr
