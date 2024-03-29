MODULE param
  INTEGER,PARAMETER :: LMAX=   12
  INTEGER           :: nato=    0
  INTEGER           :: ndif=    0
  INTEGER           :: nume=    0
  INTEGER           :: nmat=    0
  INTEGER,PARAMETER :: NRAD=  881
  ! LOMAX must be = LOMAX in LAPW1 otherwise conflict in INIT
  INTEGER,PARAMETER :: LOMAX=   3
  INTEGER,PARAMETER :: FLMAX=   3 
  INTEGER,PARAMETER :: LMAX2=  14
  ! LMX (not LMAX!) must be = LMAX in LAPW1 otherwise conflict in INIT
  INTEGER           :: LABC =    0
  INTEGER           :: LABC2=    0
  INTEGER,PARAMETER :: LMX =  LMAX+1
  INTEGER,PARAMETER :: LMX2=  LMX*LMX
  INTEGER,PARAMETER :: MMAX= 2*LMAX+1
  INTEGER           :: nume2=    0
  INTEGER           :: num2=     0
  REAL*8,PARAMETER  :: CLIGHT=137.0359895D0
  INTEGER           :: nloat=0
  INTEGER           :: nsym=0
  !rschmid
  !  Extension of relativistic local orbitals.
  !rschmid
  INTEGER,PARAMETER :: HBLOCK = 32
  INTEGER,PARAMETER :: NSLMAX=  5
  CHARACTER*180     :: filename_V_sp(2), filename_V_vns(2)
END MODULE param

MODULE abcd
  COMPLEX*16,ALLOCATABLE  :: abcdlm(:,:,:,:)
CONTAINS 
  SUBROUTINE allocate_abcd(labc2,nloat,nume)
    ALLOCATE(abcdlm(nloat,labc2,nume,2))
    abcdlm=(0.0d0,0.0d0)
  END SUBROUTINE allocate_abcd
  SUBROUTINE deallocate_abcd
    DEALLOCATE(abcdlm)
  END SUBROUTINE deallocate_abcd
END MODULE abcd

MODULE ams
  REAL*8,ALLOCATABLE          :: atom_mass(:)
CONTAINS
  SUBROUTINE init_ams
    REAL*8          :: a_m(103)
    ALLOCATE(atom_mass(103))
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
         247.,247.,251.,252.,257.,258.,259.,262./     
    atom_mass(1:103) = a_m(1:103) 
  END SUBROUTINE init_ams
END MODULE ams

MODULE couples
  COMPLEX*16,ALLOCATABLE  ::couplo(:,:,:,:,:,:)
CONTAINS
  SUBROUTINE allocate_couplo(ndif,labc)
    ALLOCATE(couplo(ndif,labc,-labc:labc,-labc:labc,2,2))
  END SUBROUTINE allocate_couplo
END MODULE couples


MODULE hmsout
  INTEGER                :: neig
  REAL*8,ALLOCATABLE     :: en(:),vnorm(:,:)
  COMPLEX*16,ALLOCATABLE :: vect(:,:,:)
CONTAINS
  SUBROUTINE allocate_hmsout(nmat,nume2)
    ALLOCATE(en(nume2),vnorm(nume2,2))
    ALLOCATE(vect(nmat,nume2,2))
    en=0.0d0; vnorm=0.0d0
    vect=(0.0d0,0.0d0)
  END SUBROUTINE allocate_hmsout
  SUBROUTINE deallocate_hmsout
    DEALLOCATE(vect)
    DEALLOCATE(en)
    DEALLOCATE(vnorm)
  END SUBROUTINE deallocate_hmsout
END MODULE hmsout

MODULE orb
  INTEGER                :: nmod,nsp,natorb,nonso,iorbpot
  INTEGER,ALLOCATABLE    :: iat(:),iiat(:),ll(:)
  REAL*8                 :: Bext
  COMPLEX*16,ALLOCATABLE :: vv(:,:,:,:)
CONTAINS
  SUBROUTINE allocate_orb(ndif)
    ALLOCATE(iat(3*ndif),iiat(3*ndif),ll(3*ndif))
    ALLOCATE(vv(3*ndif,-3:3,-3:3,3))
  END SUBROUTINE allocate_orb
END MODULE orb 

MODULE loabc
  REAL*8,ALLOCATABLE  :: alo(:,:,:,:),blo(:,:,:,:),clo(:,:,:,:)
  REAL*8,ALLOCATABLE  :: dplo(:,:,:,:),dpelo(:,:,:,:),elo(:,:,:,:)
  REAL*8,ALLOCATABLE  :: peilo(:,:,:,:),pelo(:,:,:,:)
  REAL*8,ALLOCATABLE  :: plo(:,:,:,:),pi12lo(:,:,:,:),pe12lo(:,:,:,:)
  
CONTAINS
  SUBROUTINE allocate_loabc(lomax,nloat,nato)
    ALLOCATE(alo(0:lomax,nloat,nato,2),blo(0:lomax,nloat,nato,2),clo(0:lomax,nloat,nato,2))
    ALLOCATE(dplo(0:lomax,nato,2,nloat),dpelo(0:lomax,nato,2,nloat) ) !elo(0:lomax,nloat,nato,2))
    ALLOCATE(peilo(0:lomax,nato,2,nloat),pelo(0:lomax,nato,2,nloat))
    ALLOCATE(plo(0:lomax,nato,2,nloat),pi12lo(0:lomax,nato,2,nloat),pe12lo(0:lomax,nato,2,nloat))
  END SUBROUTINE allocate_loabc
END MODULE loabc

MODULE loabcr
  real*8,ALLOCATABLE  :: alor(:,:,:), blor(:,:,:)
  real*8,ALLOCATABLE  :: clor(:,:,:), dplor(:,:,:)
  real*8,ALLOCATABLE  :: dpelor(:,:,:), elor(:,:,:)
  real*8,ALLOCATABLE  :: peilor(:,:,:), pelor(:,:,:)
  real*8,ALLOCATABLE  :: plor(:,:,:)
  real*8,ALLOCATABLE  :: pi2lor(:,:,:),pe2lor(:,:,:)
  real*8,ALLOCATABLE  :: elor2(:,:,:)

CONTAINS
  SUBROUTINE allocate_loabcr(lomax,nato)
    ALLOCATE(alor(0:lomax,nato,2), blor(0:lomax,nato,2))
    ALLOCATE(clor(0:lomax,nato,2), dplor(0:lomax,nato,2))
    ALLOCATE(dpelor(0:lomax,nato,2), elor(0:lomax,nato,2))
    ALLOCATE(peilor(0:lomax,nato,2), pelor(0:lomax,nato,2))
    ALLOCATE(plor(0:lomax,nato,2))
    ALLOCATE(pi2lor(0:lomax,nato,2),pe2lor(0:lomax,nato,2))
    ALLOCATE(elor2(0:lomax,nato,2))
  END SUBROUTINE allocate_loabcr
END MODULE loabcr

MODULE lolog
  INTEGER,ALLOCATABLE :: nlo(:),nlov(:),nlon(:),ilo(:,:),mrf(:,:)
  LOGICAL,ALLOCATABLE :: loor(:,:),lapw(:,:),lso(:)

CONTAINS
  SUBROUTINE allocate_lolog(nato,lomax,lmax)
    ALLOCATE(nlo(nato),nlov(nato),nlon(nato),ilo(0:lomax,nato),mrf(0:lmax,nato))
    ALLOCATE(loor(0:lomax,nato),lapw(0:lmax,nato),lso(nato))
    lso=.true.
  END SUBROUTINE allocate_lolog
END MODULE lolog

MODULE radovlp
  real*8,ALLOCATABLE :: ruu (:,:,:,:,:,:) 

CONTAINS
  SUBROUTINE allocate_radovlp(lmax,nato,nrf)
    ALLOCATE(ruu  (0:lmax,nato,2,2,nrf,nrf))
  END SUBROUTINE allocate_radovlp
END MODULE radovlp

MODULE rlolog
  INTEGER              :: nnrlo
  integer,ALLOCATABLE  :: nrlo(:),nrlov(:),nrlon(:)
  logical,ALLOCATABLE  :: loorext(:,:)

CONTAINS
  SUBROUTINE allocate_rlolog(nato,lomax)
    ALLOCATE(nrlo(nato),nrlov(nato),nrlon(nato))
    ALLOCATE(loorext(0:lomax,nato))
  END SUBROUTINE allocate_rlolog
END MODULE rlolog

MODULE rpars
  INTEGER,allocatable :: extl(:,:),nlr(:)
  REAL*8,ALLOCATABLE  :: extei(:,:),extde(:,:)
  CHARACTER*4,allocatable  :: extem(:,:)

CONTAINS
  SUBROUTINE allocate_rpars(nato,lomax)
    ALLOCATE(extl(nato, lomax),nlr(nato))
    ALLOCATE(extei(nato, lomax),extde(nato, lomax))
    ALLOCATE(extem(nato, lomax))
  END SUBROUTINE allocate_rpars
END MODULE rpars

MODULE rotmat
  REAL*8,ALLOCATABLE       :: det(:),phase(:)
  CONTAINS
    SUBROUTINE allocate_rotmat(ndif)
      ALLOCATE(det(ndif),phase(ndif))
    END SUBROUTINE allocate_rotmat
END MODULE rotmat

MODULE vns
  LOGICAL,ALLOCATABLE   :: lvns(:)
  INTEGER,ALLOCATABLE   :: nvns(:),mvns(:,:)
  REAL*8,ALLOCATABLE    :: vaa(:,:,:),vab(:,:,:),vad(:,:,:),vbb(:,:,:),vbd(:,:,:),vdd(:,:,:)
CONTAINS
  SUBROUTINE allocate_vns(nato)
    ALLOCATE(lvns(nato),nvns(nato),mvns(5,nato))
    ALLOCATE(vaa(nato,5,2),vab(nato,5,2),vad(nato,5,2),vbb(nato,5,2),vbd(nato,5,2),vdd(nato,5,2))
    lvns=.FALSE.
  END SUBROUTINE allocate_vns
END MODULE vns

MODULE peic
  REAL*8,ALLOCATABLE    :: pei(:,:,:)
CONTAINS
  SUBROUTINE allocate_peic(lmax,nato)
    ALLOCATE(pei(1:lmax+1,nato,2))
  END SUBROUTINE allocate_peic
END MODULE peic

MODULE hexpt
  REAL*8,ALLOCATABLE    :: hexl(:,:,:,:)
CONTAINS
  SUBROUTINE allocate_hexpt(lomax,nato)
    ALLOCATE(hexl(0:lomax,nato,2,9))
    hexl = 0.d0
  END SUBROUTINE allocate_hexpt
END MODULE hexpt

!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------
module meigve_mod
  complex*16,allocatable :: meigve(:,:,:)
end module meigve_mod

!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------

module store_vec
  complex*16,allocatable ::  vec(:,:)
end module store_vec
