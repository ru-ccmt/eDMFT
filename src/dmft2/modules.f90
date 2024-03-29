MODULE param
  INTEGER,PARAMETER :: CHUNK= 1
  INTEGER,PARAMETER :: IBLOCK= 128     ! 256
  !c.....Optimize IBLOCK for your hardware (32-255)
  INTEGER,PARAMETER :: LMAX2=    8     ! 10
  INTEGER,PARAMETER :: LOMAX=    3
  INTEGER,PARAMETER :: NCOM=   121
  ! for ncom parameter check format 1003 in l2main.frc
  INTEGER,PARAMETER :: NLOAT=    3
  INTEGER,PARAMETER :: NRAD=   881
  INTEGER           :: NGAU=  23500     ! 3200 CAREFUL: This was incresed for factor of 10 by email Rishi Rao. Check
  INTEGER           :: LXDOS=    1

  INTEGER           :: NMAT=     0
  INTEGER           :: NUME=     0
  INTEGER           :: NSYM=     0 
  INTEGER           :: NKPT=     0
  INTEGER           :: NUMKPT=   0
  INTEGER, PARAMETER:: fh_vec  = 9
  INTEGER, PARAMETER:: fh_vecdn=10
  INTEGER           :: nemin0, nemax0
  CHARACTER*180     :: filename_V_sph
  CHARACTER*180     :: filename_V_nsh
  CHARACTER*180     :: filename_V_vns
  INTEGER           :: nwave, iff1, iff2, iff3, kmax(3)  ! for fft
END MODULE param

MODULE defs
  REAL*8,PARAMETER       :: CLIGHT= 137.0359895d0
  REAL*8,PARAMETER       :: PI=     3.1415926535897932d0
  REAL*8,PARAMETER       :: TEST=   1.D-12
  REAL*8,PARAMETER       :: ZERO=   0.0d0
  REAL*8,PARAMETER       :: TWO=    2.0d0
  REAL*8,PARAMETER       :: NINETY=  90.0d0
  COMPLEX*16,PARAMETER   :: ZEROC=  (0.0d0,0.0d0)
  COMPLEX*16,PARAMETER   :: IMAG=   (0.0D0,1.0D0)
  REAL*8,PARAMETER       :: Ry2eV= 13.60569253d0 
END MODULE defs

module fact
  REAL*8  :: FCT(100)
end module fact

MODULE ams
  REAL*8          :: atom_mass(103)

 CONTAINS
  SUBROUTINE init_ams
  REAL*8          :: a_m(103)
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

MODULE atspdt
  USE param
  REAL*8                :: EL(0:LMAX2),P(0:LMAX2),DP(0:LMAX2),PE(0:LMAX2)
  REAL*8                :: DPE(0:LMAX2),PEI(0:LMAX2)
  real*8, allocatable   :: e_store(:,:)
END MODULE atspdt

!MODULE normal
!  USE param
!  REAL*8,allocatable  ::    pu1u1(:,:),pu1ue(:,:),pu1u2(:,:,:),pueue(:,:),pueu2(:,:,:),pu2u2(:,:,:,:)
!END MODULE normal

MODULE bandm
  REAL*8                :: eseper
END MODULE bandm

MODULE char
  CHARACTER*5           :: modus,efmod,modus1
END MODULE char

MODULE charp
  ! VARIABLES FOR TOTAL P SPLIT
  REAL*8,ALLOCATABLE     :: apx(:),bpx(:),cpx(:)
  REAL*8,ALLOCATABLE     :: capx(:),cbpx(:),acpx(:),bcpx(:)
  REAL*8,ALLOCATABLE     :: apy(:),bpy(:),cpy(:)
  REAL*8,ALLOCATABLE     :: capy(:),cbpy(:),acpy(:),bcpy(:)
  REAL*8,ALLOCATABLE     :: apz(:),bpz(:),cpz(:)
  REAL*8,ALLOCATABLE     :: capz(:),cbpz(:),acpz(:),bcpz(:)

 CONTAINS
  SUBROUTINE init_charp(nume)
    INTEGER nume
    ALLOCATE(apx(nume),bpx(nume),cpx(nume))
    ALLOCATE(capx(nume),cbpx(nume),acpx(nume),bcpx(nume))
    ALLOCATE(apy(nume),bpy(nume),cpy(nume))
    ALLOCATE(capy(nume),cbpy(nume),acpy(nume),bcpy(nume))
    ALLOCATE(apz(nume),bpz(nume),cpz(nume))
    ALLOCATE(capz(nume),cbpz(nume),acpz(nume),bcpz(nume))
  END SUBROUTINE init_charp

  SUBROUTINE fini_charp
    DEALLOCATE(apx,bpx,cpx)
    DEALLOCATE(capx,cbpx,acpx,bcpx)
    DEALLOCATE(apy,bpy,cpy)
    DEALLOCATE(capy,cbpy,acpy,bcpy)
    DEALLOCATE(apz,bpz,cpz)
    DEALLOCATE(capz,cbpz,acpz,bcpz)
  END SUBROUTINE fini_charp

  SUBROUTINE zero_charp(nemin,nemax)
    USE defs
    INTEGER nemin,nemax
    apx(nemin:nemax)=zero; bpx(nemin:nemax)=zero; cpx(nemin:nemax)=zero
    capx(nemin:nemax)=zero; cbpx(nemin:nemax)=zero; acpx(nemin:nemax)=zero
    bcpx(nemin:nemax)=zero
    apy(nemin:nemax)=zero; bpy(nemin:nemax)=zero; cpy(nemin:nemax)=zero
    capy(nemin:nemax)=zero; cbpy(nemin:nemax)=zero; acpy(nemin:nemax)=zero
    bcpy(nemin:nemax)=zero
    apz(nemin:nemax)=zero; bpz(nemin:nemax)=zero; cpz(nemin:nemax)=zero
    capz(nemin:nemax)=zero; cbpz(nemin:nemax)=zero; acpz(nemin:nemax)=zero
    bcpz(nemin:nemax)=zero
  END SUBROUTINE zero_charp
END MODULE charp

MODULE chard
  ! VARIABLES FOR TOTAL D SPLIT
  REAL*8,ALLOCATABLE    :: adz2(:),bdz2(:),cdz2(:)
  REAL*8,ALLOCATABLE    :: cadz2(:),cbdz2(:),acdz2(:),bcdz2(:)
  REAL*8,ALLOCATABLE    :: adx2y2(:),bdx2y2(:),cdx2y2(:)
  REAL*8,ALLOCATABLE    :: cadx2y2(:),cbdx2y2(:),acdx2y2(:),bcdx2y2(:)
  REAL*8,ALLOCATABLE    :: adxy(:),bdxy(:),cdxy(:)
  REAL*8,ALLOCATABLE    :: cadxy(:),cbdxy(:),acdxy(:),bcdxy(:)
  REAL*8,ALLOCATABLE    :: adxz(:),bdxz(:),cdxz(:)
  REAL*8,ALLOCATABLE    :: cadxz(:),cbdxz(:),acdxz(:),bcdxz(:)
  REAL*8,ALLOCATABLE    :: adyz(:),bdyz(:),cdyz(:)
  REAL*8,ALLOCATABLE    :: cadyz(:),cbdyz(:),acdyz(:),bcdyz(:)

 CONTAINS
  SUBROUTINE init_chard(nume)
    INTEGER :: nume
    ALLOCATE(adz2(nume),bdz2(nume),cdz2(nume))
    ALLOCATE(cadz2(nume),cbdz2(nume),acdz2(nume),bcdz2(nume))
    ALLOCATE(adx2y2(nume),bdx2y2(nume),cdx2y2(nume))
    ALLOCATE(cadx2y2(nume),cbdx2y2(nume),acdx2y2(nume),bcdx2y2(nume))
    ALLOCATE(adxy(nume),bdxy(nume),cdxy(nume))
    ALLOCATE(cadxy(nume),cbdxy(nume),acdxy(nume),bcdxy(nume))
    ALLOCATE(adxz(nume),bdxz(nume),cdxz(nume))
    ALLOCATE(cadxz(nume),cbdxz(nume),acdxz(nume),bcdxz(nume))
    ALLOCATE(adyz(nume),bdyz(nume),cdyz(nume))
    ALLOCATE(cadyz(nume),cbdyz(nume),acdyz(nume),bcdyz(nume))
  END SUBROUTINE init_chard

  SUBROUTINE fini_chard
    DEALLOCATE(adz2,bdz2,cdz2)
    DEALLOCATE(cadz2,cbdz2,acdz2,bcdz2)
    DEALLOCATE(adx2y2,bdx2y2,cdx2y2)
    DEALLOCATE(cadx2y2,cbdx2y2,acdx2y2,bcdx2y2)
    DEALLOCATE(adxy,bdxy,cdxy)
    DEALLOCATE(cadxy,cbdxy,acdxy,bcdxy)
    DEALLOCATE(adxz,bdxz,cdxz)
    DEALLOCATE(cadxz,cbdxz,acdxz,bcdxz)
    DEALLOCATE(adyz,bdyz,cdyz)
    DEALLOCATE(cadyz,cbdyz,acdyz,bcdyz)
  END SUBROUTINE fini_chard

  SUBROUTINE zero_chard(nemin,nemax)
    USE defs
    INTEGER :: nemin,nemax
    adz2(nemin:nemax)=zero;bdz2(nemin:nemax)=zero;cdz2(nemin:nemax)=zero
    cadz2(nemin:nemax)=zero;cbdz2(nemin:nemax)=zero;acdz2(nemin:nemax)=zero;bcdz2(nemin:nemax)=zero
    adx2y2(nemin:nemax)=zero;bdx2y2(nemin:nemax)=zero;cdx2y2(nemin:nemax)=zero
    cadx2y2(nemin:nemax)=zero;cbdx2y2(nemin:nemax)=zero;acdx2y2(nemin:nemax)=zero;bcdx2y2(nemin:nemax)=zero
    adxy(nemin:nemax)=zero;bdxy(nemin:nemax)=zero;cdxy(nemin:nemax)=zero
    cadxy(nemin:nemax)=zero;cbdxy(nemin:nemax)=zero;acdxy(nemin:nemax)=zero;bcdxy(nemin:nemax)=zero
    adxz(nemin:nemax)=zero;bdxz(nemin:nemax)=zero;cdxz(nemin:nemax)=zero
    cadxz(nemin:nemax)=zero;cbdxz(nemin:nemax)=zero;acdxz(nemin:nemax)=zero;bcdxz(nemin:nemax)=zero
    adyz(nemin:nemax)=zero;bdyz(nemin:nemax)=zero;cdyz(nemin:nemax)=zero
    cadyz(nemin:nemax)=zero;cbdyz(nemin:nemax)=zero;acdyz(nemin:nemax)=zero;bcdyz(nemin:nemax)=zero
  END SUBROUTINE zero_chard
END MODULE chard

MODULE charf
  ! VARIABLES FOR TOTAL F SPLIT
  REAL*8,ALLOCATABLE    :: af00(:),bf00(:),af11(:)
  REAL*8,ALLOCATABLE    :: bf11(:),af22(:),bf22(:),af33(:)
  REAL*8,ALLOCATABLE    :: bf33(:),af1m(:),bf1m(:),af2m(:)
  REAL*8,ALLOCATABLE    :: bf2m(:),af3m(:),bf3m(:)

 CONTAINS
  SUBROUTINE init_charf(nume)
    INTEGER :: nume
    ALLOCATE(af00(nume),bf00(nume),af11(nume))
    ALLOCATE(bf11(nume),af22(nume),bf22(nume),af33(nume))
    ALLOCATE(bf33(nume),af1m(nume),bf1m(nume),af2m(nume))
    ALLOCATE(bf2m(nume),af3m(nume),bf3m(nume))
  END SUBROUTINE init_charf

  SUBROUTINE fini_charf
    DEALLOCATE(af00,bf00,af11)
    DEALLOCATE(bf11,af22,bf22,af33)
    DEALLOCATE(bf33,af1m,bf1m,af2m)
    DEALLOCATE(bf2m,af3m,bf3m)
  END SUBROUTINE fini_charf

  SUBROUTINE zero_charf(nemin,nemax)
    USE defs
    INTEGER :: nemin,nemax
    af00(nemin:nemax)=zero;bf00(nemin:nemax)=zero;af11(nemin:nemax)=zero
    bf11(nemin:nemax)=zero;af22(nemin:nemax)=zero;bf22(nemin:nemax)=zero;af33(nemin:nemax)=zero
    bf33(nemin:nemax)=zero;af1m(nemin:nemax)=zero;bf1m(nemin:nemax)=zero;af2m(nemin:nemax)=zero
    bf2m(nemin:nemax)=zero;af3m(nemin:nemax)=zero;bf3m(nemin:nemax)=zero
  END SUBROUTINE zero_charf
END MODULE charf

MODULE com
  USE param
  LOGICAL               :: rel
  INTEGER               :: NSPIN,NAT,NBAND,NK,MINWAV,MAXWAV
  REAL*8                :: EMIN, EMAX, EF,ELECN,XWT
END MODULE com

MODULE lo
  USE param
  LOGICAL    ::   loor(nloat,0:lomax),rlo(nloat,0:lomax),lapw(0:lmax2)
  ! nlo #LO on this atom, nlov #LO up til now, nlon #LO left
  INTEGER    ::   nlo=0,nlov=0,nlon=0,ilo(0:lmax2)
  REAL*8     ::   alo(0:lomax,nloat),blo(0:lomax,nloat)
  REAL*8     ::   clo(0:lomax,nloat)
  REAL*8     ::   elo(0:lomax,nloat),plo(nloat,0:lomax)
  REAL*8     ::   dplo(nloat,0:lomax)
  REAL*8     ::   pi12lo(nloat,0:lomax),pe12lo(nloat,0:lomax),pr12lo(nloat,nloat,0:lomax)
  REAL*8     ::   a1lo(nrad,nloat,0:lomax),b1lo(nrad,nloat,0:lomax)
  real*8 , allocatable :: elo_store(:,:,:)
END MODULE lo

MODULE lohelp
  USE param
  REAL*8  :: u21(nrad,nloat),ue21(nrad,nloat),u12(nrad,nloat)
  REAL*8  :: ue12(nrad,nloat),u22(nrad,nloat,nloat)

  COMPLEX*16,ALLOCATABLE  :: sum12(:,:),sum21(:,:),sum22(:,:,:)
  COMPLEX*16,ALLOCATABLE  :: sume21(:,:),sume12(:,:)

  REAL*8,ALLOCATABLE  :: tc12(:,:),tc21(:,:)
  REAL*8,ALLOCATABLE  :: tce12(:,:),tce21(:,:)
  REAL*8,ALLOCATABLE  :: tc22(:,:)
  
 CONTAINS
  SUBROUTINE init_lohelp
    ALLOCATE(sum12(nume,nloat),sum21(nume,nloat),sum22(nume,nloat,nloat))
    ALLOCATE(sume21(nume,nloat),sume12(nume,nloat))
    ALLOCATE(tc12(0:lmax2,nume),tc21(0:lmax2,nume))
    ALLOCATE(tce12(0:lmax2,nume),tce21(0:lmax2,nume))
    ALLOCATE(tc22(0:lmax2,nume))
  END SUBROUTINE init_lohelp

  SUBROUTINE zerotc_lohelp(nemin,nemax)
    USE defs
    tc12(0:lmax2,nemin:nemax)=zero;tc21(0:lmax2,nemin:nemax)=zero
    tce12(0:lmax2,nemin:nemax)=zero; tce21(0:lmax2,nemin:nemax)=zero
    tc22(0:lmax2,nemin:nemax)=zero
  END SUBROUTINE zerotc_lohelp

  SUBROUTINE fini_lohelp
    DEALLOCATE(sum12,sum21,sum22)
    DEALLOCATE(sume21,sume12)
    DEALLOCATE(tc12,tc21)
    DEALLOCATE(tce12,tce21)
    DEALLOCATE(tc22)
  END SUBROUTINE fini_lohelp
END MODULE lohelp

MODULE ORB
  REAL*8                   :: BEXT
  COMPLEX*16, ALLOCATABLE  :: VORB(:,:,:,:)
  INTEGER                  :: NMOD, NSP, NATORB
  INTEGER, ALLOCATABLE     :: IATOM(:), NLORB(:), LORB(:,:)
  logical                  :: forb_log

  CONTAINS

    SUBROUTINE INIT_ORB(NAT)
      IMPLICIT NONE
      INTEGER NAT
      ALLOCATE(VORB(NAT,0:3,-3:3,-3:3))
      ALLOCATE(IATOM(NAT))
      ALLOCATE(NLORB(NAT))
      ALLOCATE(LORB(4,NAT))
      VORB = (0D0,0D0)
      IATOM = 0
      NLORB = 0
      LORB = 0
      NATORB=0
    END SUBROUTINE INIT_ORB
END MODULE ORB

MODULE reclat
  INTEGER,POINTER        :: kzz(:,:) !vectors for dens ekspansion
  INTEGER,ALLOCATABLE    :: inst(:)  !size of star
  COMPLEX*16,ALLOCATABLE :: tauk(:)
END MODULE reclat

MODULE sym2
  INTEGER               :: iord
  INTEGER,ALLOCATABLE   :: iz(:,:,:), idet(:)
  REAL*8,ALLOCATABLE    :: tau(:,:), phase(:)

CONTAINS
  SUBROUTINE init_sym2(nsym)
    INTEGER nsym
    ALLOCATE(tau(3,nsym),iz(3,3,nsym),idet(nsym),phase(nsym))
    idet(:)=1
    phase(:)=0.d0
  END SUBROUTINE init_sym2
END MODULE sym2

MODULE xa
  USE param
  INTEGER :: LM(2,NCOM)
  REAL*8  :: R(NRAD),RHOLM(NRAD,NCOM),AVEC(3),BK(3),BKROT(3),BKRLOC(3)!,XWT1(0:21)
  REAL*8,ALLOCATABLE     :: fj(:,:,:),dfj(:,:,:)
  REAL*8,ALLOCATABLE     :: E(:)
  REAL*8,ALLOCATABLE     :: WEIGHT(:)
  REAL*8,ALLOCATABLE     :: TC100(:,:),TCA100(:,:),TCB100(:,:)
  COMPLEX*16,ALLOCATABLE :: SUMA(:),SUMB(:),SUMAB(:),SUMBA(:)
  !COMPLEX*16,ALLOCATABLE :: phs(:)
CONTAINS
  SUBROUTINE init_xa(nat)
    IMPLICIT NONE
    INTEGER :: nat
    ALLOCATE(FJ(0:LMAX2,NMAT,nat),DFJ(0:LMAX2,NMAT,nat))
    ALLOCATE(E(NUME),WEIGHT(NUME))
    ALLOCATE(TC100(0:LMAX2,NUME),TCA100(0:LMAX2,NUME),TCB100(0:LMAX2,NUME))
    ALLOCATE(SUMA(NUME),SUMB(NUME),SUMAB(NUME),SUMBA(NUME))
  END SUBROUTINE init_xa

  SUBROUTINE fini_xa
    DEALLOCATE(FJ,DFJ)
    DEALLOCATE(E,WEIGHT)
    DEALLOCATE(TC100,TCA100,TCB100)
    DEALLOCATE(SUMA,SUMB,SUMAB,SUMBA)
  END SUBROUTINE fini_xa
END MODULE xa

MODULE xa2
  !REAL*8,allocatable   :: WEIGHT(:),E(:)
  !INTEGER,allocatable  :: NE(:)
  REAL*8, allocatable   :: lda_weight(:,:)
  INTEGER :: lda_nbmax
 CONTAINS
  SUBROUTINE init_xa2(nume,nkpt)
    IMPLICIT NONE
    INTEGER nume,nkpt
    !ALLOCATE(WEIGHT(nume*2*NKPT),E(2*NKPT*NUME))
    !ALLOCATE(NE(2*NKPT))
    !ne=0; weight=0.0d0; e=0.0d0
  END SUBROUTINE init_xa2
END MODULE xa2

MODULE xa3
  COMPLEX*16,ALLOCATABLE ::  As_lo(:,:,:)
  COMPLEX*16,ALLOCATABLE ::  As(:,:,:)
  REAL*8,ALLOCATABLE   :: bk3(:,:), aK(:)
  INTEGER,ALLOCATABLE  :: k3(:,:)
  REAL*8,ALLOCATABLE    :: bk3lo(:,:)
  INTEGER,ALLOCATABLE   :: k3lo(:,:)
CONTAINS
  SUBROUTINE init_xa3(nmat,nnlo,nume,iso)
    INTEGER nmat,nnlo,nume
    if (.not.allocated(k3)) then
       ALLOCATE(k3(3,nmat))
       ALLOCATE(bk3(3,nmat))
       ALLOCATE(aK(nmat))
       ALLOCATE(k3lo(3,nnlo))
       ALLOCATE(bk3lo(3,nnlo))
       ALLOCATE(As(nmat,nume,iso),As_lo(nnlo,nume,iso))
    endif
  END SUBROUTINE init_xa3

  SUBROUTINE fini_xa3
    if (allocated(k3)) then
       DEALLOCATE(k3)
       DEALLOCATE(bk3)
       DEALLOCATE(aK)
       DEALLOCATE(k3lo)
       DEALLOCATE(bk3lo)
       DEALLOCATE(As, As_lo)
    endif
  END SUBROUTINE fini_xa3
END MODULE xa3


!--------------------------------------------------------
! For determining the chemical potential
!--------------------------------------------------------
MODULE muzero
  REAL*8                :: w_sum, w_gamma, w_norm1, w_beta
  INTEGER               :: w_nomega, w_npomega
  INTEGER               :: n0_om, nkp, max_nbands
  LOGICAL               :: wprint, wmatsubara, wprint1, Qcheckbands
  REAL*8,    ALLOCATABLE:: Ek(:,:), wgh(:)
  COMPLEX*16,ALLOCATABLE:: zEk(:,:,:)
  INTEGER,   ALLOCATABLE:: nemm(:,:)
  REAL*8,   ALLOCATABLE :: abom(:,:), w_omega(:)
  INTEGER,  ALLOCATABLE :: nomq(:), iomq(:,:), jomq(:)
  REAL*8,   ALLOCATABLE :: womq(:,:)
  !
CONTAINS

  SUBROUTINE init_muzero(nomega, qmax)
    INTEGER, intent(in) :: nomega, qmax
    ALLOCATE( nomq(nomega), jomq(nomega), iomq(nomega,qmax), womq(nomega,qmax), w_omega(nomega) )
    ALLOCATE( abom(2,nomega) )
    wprint1=.FALSE.
  END SUBROUTINE init_muzero
  
  SUBROUTINE fini_muzero()
    DEALLOCATE( abom )
    DEALLOCATE( nomq, jomq, iomq, womq, w_omega )
  END SUBROUTINE fini_muzero
  SUBROUTINE print_muzero(fh,elecn,DM_EF)
    IMPLICIT NONE
    INTEGER, intent(in) :: fh
    REAL*8,  intent(in) :: elecn, DM_EF
    !
    INTEGER :: ikp, iom, j
    REAL*8  :: wkp
    WRITE(fh) elecn, DM_EF, nkp, size(zEk,1)
    WRITE(fh) w_nomega, n0_om, w_beta, size(iomq,2), wmatsubara
    WRITE(fh) w_npomega, w_gamma
    WRITE(fh) wprint,wprint1,Qcheckbands
    do iom=1,w_nomega
       WRITE(fh) w_omega(iom), abom(1,iom), abom(2,iom), nomq(iom), jomq(iom)
       do j=1,nomq(iom)
          WRITE(fh) womq(iom,j), iomq(iom,j)
       enddo
    enddo
    do ikp=1,nkp
       wkp = (wgh(ikp)/w_sum*2.0)*w_norm1
       WRITE(fh) ikp, wkp, nemm(1,ikp), nemm(2,ikp), nemm(3,ikp)
       do iom=1,w_nomega
          do j = 1,nemm(3,ikp)-nemm(2,ikp)+1
             WRITE(fh) zEk(j,iom,ikp)
          enddo
       enddo
    enddo
  END SUBROUTINE print_muzero
END MODULE muzero

MODULE mVks
  INTEGER    :: nwavev
  INTEGER,    ALLOCATABLE :: Gs(:,:)
  COMPLEX*16, ALLOCATABLE :: Vint(:)
END MODULE mVks
