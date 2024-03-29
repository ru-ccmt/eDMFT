MODULE param
  INTEGER            LMAX2, LOMAX
  INTEGER            NRAD, nloat, nrf
  integer            iblck
  real*8             clight
  parameter (IBLCK=   64)
  PARAMETER (IBLOCK= 128)
!.....Optimize IBLOCK and IBLCK for your hardware (32-255)
  PARAMETER (LMAX2=   3)  
  PARAMETER (lmaxr=  13)
  PARAMETER (LOMAX=   3)
  PARAMETER (NDIM=  2*lmax2+1)
  PARAMETER  (NDIM2= 2*NDIM)                                              
  INTEGER          :: NATO=   0
! for ncom parameter check format 1003 in l2main.frc
  !INTEGER          :: NDIF=   0
  INTEGER          :: NMAT=0
  PARAMETER (NRAD=  881)                                              
  INTEGER          :: NSYM=   0 
  INTEGER          :: NUME=  0 ! 1000
  INTEGER          :: maxbands = 0
  REAL*8           :: gamma, gammac, aom_default, bom_default
  INTEGER          :: nom_default
  LOGICAL          :: matsubara, Cohfacts
  LOGICAL          :: ComputeLogGloc
  INTEGER          :: nemin0, nemax0
  INTEGER          :: max_nl
  LOGICAL          :: cmp_partial_dos
  CHARACTER*2      :: Hrm
! for x-dos set lxdos to 3
  parameter (lxdos= 3)
  parameter (nloat= 3)
  parameter (nrf=4)
  PARAMETER (CLIGHT=137.0359895D0)
END MODULE param

MODULE xa
  USE param
  REAL*8 :: R(nrad)!,RHOLM(NRAD,NCOM),AVEC(3),BK(3),BKROT(3),BKRLOC(3)!,XWT1(0:21)
  !INTEGER :: LM(2,NCOM)
  !REAL*8,ALLOCATABLE     :: fj(:,:,:),dfj(:,:,:)
  !REAL*8,ALLOCATABLE     :: E(:)
  !REAL*8,ALLOCATABLE     :: WEIGHT(:)
CONTAINS
  SUBROUTINE init_xa()
    IMPLICIT NONE
    !ALLOCATE(FJ(0:LMAX2,NMAT,nat),DFJ(0:LMAX2,NMAT,nat))
    !ALLOCATE(E(NUME),WEIGHT(NUME))
  END SUBROUTINE init_xa

  SUBROUTINE fini_xa
    !DEALLOCATE(FJ,DFJ)
    !DEALLOCATE(E,WEIGHT)
    !DEALLOCATE(TC100,TCA100,TCB100)
    !DEALLOCATE(SUMA,SUMB,SUMAB,SUMBA)
  END SUBROUTINE fini_xa
END MODULE xa

!MODULE struct
!  USE param
!  LOGICAL                  :: rel
!  REAL*8                   :: AA,BB,CC,alpha(3),pia(3),VOL
!  REAL*8,ALLOCATABLE       :: R0(:),DX(:),RMT(:),zz(:),rotloc(:,:,:),v(:)
!  REAL*8,ALLOCATABLE       :: tau(:,:)
!  REAL*8,POINTER           :: pos(:,:)
!  CHARACTER*4              :: lattic,irel,cform
!  CHARACTER*80             :: title
!  CHARACTER*10,ALLOCATABLE :: aname(:)
!  INTEGER                  :: nat,iord
!  INTEGER,ALLOCATABLE      :: mult(:),jrj(:),iatnr(:),isplit(:)
!  INTEGER,ALLOCATABLE      :: iz(:,:,:),inum(:)
!  REAL*8,ALLOCATABLE       :: rotij(:,:,:),tauij(:,:) ! ndif
!
! CONTAINS
!  SUBROUTINE init_struct
!    USE reallocate
!    IMPLICIT NONE
!
!    INTEGER                :: ios
!    REAL*8                 :: test,ninety
!!loop indexs
!    INTEGER                :: index,i,j,j1,j2,m,jatom
!
!    test=1.D-5
!    ninety=90.0D0
!
!    read (20,1000) title
!    read (20,1010) lattic,nat,cform,irel
!    nato=nat
!    REL=.TRUE.                                     
!    IF(IREL.EQ.'NREL') REL=.FALSE.                                    
!    ALLOCATE(aname(nato),mult(nato),jrj(nato),r0(nato),dx(nato),rmt(nato))
!    allocate(zz(nato),rotloc(3,3,nato),iatnr(nato),isplit(nato),v(nato))
!    v=0.0d0
!    ALLOCATE (pos(3,48*nat))
!    read (20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
!    IF(ABS(ALPHA(1)).LT.test) ALPHA(1)=ninety
!    IF(ABS(ALPHA(2)).LT.test) ALPHA(2)=ninety
!    IF(ABS(ALPHA(3)).LT.test) ALPHA(3)=ninety
!    INDEX=0
!    
!    DO jatom=1,NAT
!       INDEX=INDEX+1
!       READ(20,1030,iostat=ios) iatnr(jatom),( pos(j,index),j=1,3 ), mult(jatom),isplit(jatom) 
!       IF(ios /= 0 ) THEN
!          WRITE(6,*) iatnr(jatom),( pos(j,index),j=1,3 ), mult(jatom),isplit(jatom) 
!          WRITE(6,*) 'ERROR IN STRUCT FILE READ'
!          STOP
!       ENDIF
!       IF (mult(jatom) .EQ. 0) THEN
!          WRITE (6,6000) jatom, index, mult(jatom)
!          STOP
!       ENDIF
!       DO m=1,mult(jatom)-1                                     
!          index=index+1                                            
!          READ(20,1031) iatnr(jatom),( pos(j,index),j=1,3)   ! pos -- position inside the unit cell read from case.struct
!       ENDDO
!       READ(20,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom),zz(jatom) ! zz-nuclear charge, jrj--number of radial data points
!       dx(jatom)=LOG(rmt(jatom)/r0(jatom)) / (jrj(jatom)-1)           
!       rmt(jatom)=r0(jatom)*EXP( dx(jatom)*(jrj(jatom)-1) )           
!       READ(20,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)                
!    ENDDO
!    ndif=index
!    CALL doreallocate(pos, 3, ndif)
!    READ(20,1151) iord
!    nsym=iord
!    ALLOCATE(iz(3,3,nsym),tau(3,nsym),inum(nsym))
!    DO j=1,iord  ! iz(:,:,iord) - all symmetry transformations
!                 ! tau(:,iord)  - translations
!       READ(20,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
!    ENDDO
!
!1000 FORMAT(A80)                                                       
!1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
!1020 FORMAT(6F10.7,10X,F10.7)                                          
!1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
!1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
!1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
!1051 FORMAT(20X,3F10.8)
!1101 FORMAT(3(3I2,F10.8/),I8)
!1151 FORMAT(I4)
!6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
!          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
!  END SUBROUTINE init_struct
!
!END MODULE struct

MODULE case
  COMPLEX*16,ALLOCATABLE :: cf(:,:,:)
  INTEGER,ALLOCATABLE    :: nl(:),ll(:,:),iatom(:),qsplit(:,:),cix(:,:)
  INTEGER,ALLOCATABLE    :: Sigind(:,:,:), csize(:)
  CHARACTER*30,ALLOCATABLE::legend(:,:)
  REAL*8, ALLOCATABLE    :: shft(:,:)
  REAL*8, ALLOCATABLE    :: crotloc(:,:,:,:)
  INTEGER                :: natom, ncix, maxsize, maxdim
  INTEGER, ALLOCATABLE   :: isort(:), ifirst(:)
END MODULE case

MODULE sym2
   INTEGER,ALLOCATABLE    :: idet(:)
   REAL*8,ALLOCATABLE     :: iz_cartesian(:,:,:),phase(:),tmat(:,:,:)
 CONTAINS
  SUBROUTINE init_sym2(iord)
   ALLOCATE(idet(iord),iz_cartesian(3,3,iord),phase(iord),tmat(3,3,iord))
   idet=1
   phase=0
  END SUBROUTINE init_sym2
END MODULE sym2

MODULE com
  INTEGER                :: nband,iso,iso2,minwav,maxwav,nspin1
  INTEGER                :: klmax= 0          !pb
  REAL*8                 :: emin,emax,ef
END MODULE com

MODULE abc
   INTEGER,ALLOCATABLE    ::  kx(:),ky(:),kz(:)
   !REAL*8,ALLOCATABLE    ::  bkx(:),bky(:),bkz(:)
   REAL*8,ALLOCATABLE     ::  bk3(:,:), aK(:)
   REAL*8,ALLOCATABLE     ::  e(:),xden(:,:,:)
   COMPLEX*16,ALLOCATABLE ::  a(:,:,:),alm(:,:,:,:),alml(:,:,:,:,:)
   COMPLEX*16,ALLOCATABLE ::  dmat(:,:,:,:)
 CONTAINS
  SUBROUTINE init_abc(nume,nmat,lmax2,ndim2,nrf)
    ALLOCATE(kx(nmat),ky(nmat),kz(nmat) )
    ALLOCATE( bk3(3,nmat), aK(nmat) )
    ALLOCATE(a(nmat,nume,2),alm(2*lmax2+1,nume,nrf,2))
    ALLOCATE(alml(0:lmax2,2*lmax2+1,nume,nrf,2))
    ALLOCATE(e(nume),xden(0:lmax2,ndim2,nume)) 
    ALLOCATE(dmat(0:lmax2,0:lmax2,ndim2,ndim2))
  END SUBROUTINE init_abc
  SUBROUTINE deallocate_abc()
    DEALLOCATE( kx,ky,kz )
    DEALLOCATE( bk3, aK )
    DEALLOCATE( a, alm )
    DEALLOCATE( alml )
    DEALLOCATE( e, xden ) 
    DEALLOCATE( dmat )
  END SUBROUTINE deallocate_abc
END MODULE abc
  

! JR: from lapw2 - for rotdef.f functionality
MODULE defs
  REAL*8,PARAMETER       :: CLIGHT= 137.0359895d0
  REAL*8,PARAMETER       :: PI=     3.1415926535897932d0
  REAL*8,PARAMETER       :: TEST=   1.D-12
  REAL*8,PARAMETER       :: ZERO=   0.0d0
  REAL*8,PARAMETER       :: TWO=    2.0d0
  REAL*8,PARAMETER       :: NINETY=  90.0d0
  COMPLEX*16,PARAMETER   :: ZEROC=  (0.0d0,0.0d0)
  COMPLEX*16,PARAMETER   :: IMAG=   (0.0D0,1.0D0)
END MODULE defs


! from lapw1 for k-point information
module kpts
! MSX(j) - x-component of K-point j
! MSY(j) - y-component of K-point j
! MSZ(j) - z-component of K-point j
  DOUBLE PRECISION, allocatable :: MSX(:),  MSY(:),  MSZ(:), MWEIGHT(:)
  CHARACTER*10, allocatable :: MKNAME(:)
  CHARACTER*3, allocatable :: MIPGR(:)
  INTEGER                  :: numkpt
  REAL*8                   :: TWEIGHT
contains
  SUBROUTINE allocate_kpts(numkpt)
    IMPLICIT NONE
    integer numkpt
    allocate( MSX(numkpt), MSY(numkpt), MSZ(numkpt), MWEIGHT(numkpt) , MKNAME(numkpt), MIPGR(numkpt))
  END SUBROUTINE allocate_kpts
  SUBROUTINE deallocate_kpts()
    IMPLICIT NONE
    deallocate( MSX, MSY, MSZ, MWEIGHT, MKNAME, MIPGR )
  END SUBROUTINE deallocate_kpts
end module kpts

