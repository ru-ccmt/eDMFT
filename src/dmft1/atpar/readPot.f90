module readPotential
  INTEGER :: ifhs ! case.vsp (old 18) : spherical potential on real space grid
  INTEGER :: ifhr ! case.vns (old 19) : non-spherical potential on real space grid
  INTEGER :: ifhn ! case.nsh (old  2) : matrix elements of the non-spherical potential in terms of u,udot,ulo
contains

!!!************** SPHERICAL POTENTIAL ON REAL SPACE GRID ***************************
  SUBROUTINE init_V_vsp(fname, ifhs_, ISCF) ! ifhs is normally 18
    IMPLICIT NONE
    CHARACTER*180, intent(in) :: fname
    INTEGER, intent(in)       :: ifhs_
    INTEGER, intent(out)      :: ISCF
    ! local
    ifhs = ifhs_
    open(ifhs, file=fname, status='old', err=87)
    read(ifhs,'(49X,I3,//)',err=87) ISCF
    !print*, ISCF
    return
87  continue
    WRITE(6,*) 'ERROR opening file case.vsp'
  END SUBROUTINE init_V_vsp
  
  SUBROUTINE close_V_vsp() 
    close(ifhs)
  END SUBROUTINE close_V_vsp

  SUBROUTINE read_V_vsp(VR,jatom,nr)
    IMPLICIT NONE
    REAL*8, intent(out):: VR(nr)
    INTEGER,intent(in) :: nr, jatom
    ! locals
    INTEGER :: jatom_read, ir
    !-----------------------------------------------------------------------------------------------
    ! READ TOTAL SPHERICAL POTENTIAL V(0,0) OF TAPE18=VSP
    ! NORM OF V(0,0)=V00(R)*R/SQRT(4.D0*PI)
    ! Notice that V on file (V_file) is in special form:  V_file = V_{Rydberg-units}(r)*sqrt(4*pi)/r
    !
    VR(:)=0.0d0
    READ(ifhs,'(15X,i3)') jatom_read
    if (jatom_read.NE.jatom) print*, 'ERROR in read_V_vsp reading case.vsp!', jatom_read, jatom
    READ(ifhs,'(16X,//)')
    READ(ifhs,'(/)')
    READ(ifhs,'(3X,4E19.12)') ( VR(ir),ir=1,nr )
    READ(ifhs,'(/)')
    READ(ifhs,'(///)')
  END SUBROUTINE read_V_vsp
  
!!!************** NON SPHERICAL POTENTIAL ON REAL SPACE GRID ***************************
  SUBROUTINE init_V_vns(fname, ifhr_) ! ifhr is normally 9919 (old was 19)
    ! Opens case.vns :  non-spherical potential on real space grid
    IMPLICIT NONE
    CHARACTER*180, intent(in) :: fname
    INTEGER, intent(in)       :: ifhr_
    ifhr = ifhr_
    open(ifhr, file=fname, status='old', err=88)
    read(ifhr,'(//)',err=88)
    return
88  continue
    WRITE(6,*) 'ERROR opening file case.vns'
  END SUBROUTINE init_V_vns
  
  SUBROUTINE close_V_vns() 
    close(ifhr)
  END SUBROUTINE close_V_vns

  !SUBROUTINE add_V00_to_vns(Vlm, Vr, R, Nrad, nr, lmmax)
  !  REAL*8, intent(in)   :: Vr(Nrad), R(Nrad)
  !  INTEGER, intent(in)  :: Nrad, nr, lmmax
  !  REAL*8, intent(inout):: Vlm(nr,0:lmmax-1)
  !  ! Spherically symmetric part was read before and is stored in Vr. Just add it to Vlm.
  !  ! Note that Vr in the file was kept in the following form:     V_file = V_{Rydberg} * R/sqrt(4*pi),
  !  !  because the quantity (V*R-2*Z)/R enters Schroedinger equation
  !  ! Note that Vr was changed to Hartree units in atpar, but here it is changed back to Rydbergs
  !  Vlm(1:nr,0) = Vr(1:nr)/R(1:nr) * (sqrt(4.d0*pi)) ! No need for conversion back to Rydbergs
  !END SUBROUTINE add_V00_to_vns

  SUBROUTINE read_V_vns(Vlm,jatom,nrad,nr,lmmax)
    IMPLICIT NONE
    REAL*8, intent(out):: Vlm(nrad,0:lmmax-1)
    INTEGER,intent(in) :: nr, nrad, lmmax, jatom
    ! locals
    INTEGER :: jatom_read, lmp, ilm, llp, mp, ir
    REAL*8, allocatable :: dummy(:)
    Vlm(:,:)=0.d0
    ! Reading case.vns, which contains non-spherical potential on radial mesh grid
    !     loop over potential
    read(ifhr,'(15x,i3)') jatom_read
    if (jatom_read.NE.jatom) print*, 'ERROR in read_V_vns reading case.vns!'
    read(ifhr,'(15x,i3//)') lmp
    DO ilm=1,lmp
       !  normalize potential with vlm=V(l,m)*r*r
       read(ifhr,'(15x,i3,5x,i2/)') llp,mp
       if (ilm.le.lmmax-1) then
          read(ifhr,'(3x,4e19.12)') ( Vlm(ir,ilm), ir=1,nr )
       else ! It should not happen, but just in case the file contains too many lm components
          allocate( dummy(nr) )
          read(ifhr,'(3x,4e19.12)') ( dummy(ir), ir=1,nr )
          deallocate( dummy )
       endif
       read(ifhr,'(/)')
    ENDDO
    read(ifhr,'(///)')
  END SUBROUTINE read_V_vns
  
  !!!************** NON SPHERICAL POTENTIAL MATRIX ELEMENTS ***************************
  SUBROUTINE init_V_nsh(fname, ifhn_) ! ifhn is normally 9902 (old was 2)
    ! Opens case.nsh :  matrix elements of the non-spherical potential in terms of u,udot,ulo    
    IMPLICIT NONE
    CHARACTER*180, intent(in) :: fname
    INTEGER, intent(in)       :: ifhn_
    ifhn = ifhn_
    open(ifhn, file=fname, status='old', err=89)
    read(ifhn,'(//)',err=89)
    return
89  continue
    WRITE(6,*) 'ERROR opening file case.nsh'
  END SUBROUTINE init_V_nsh
  
  SUBROUTINE close_V_nsh() 
    close(ifhn)
  END SUBROUTINE close_V_nsh

  SUBROUTINE read_V_nsh(Vts,lv,lpv,mv,mpv,ihmx,jatom,lomax,ngau)
    IMPLICIT NONE
    COMPLEX*16, intent(out):: Vts(3,3,ngau)       ! non-spherical potential matrix elements
    INTEGER, intent(out)   :: lpv(ngau), lv(ngau), mpv(ngau), mv(ngau)
    INTEGER, intent(out)   :: ihmx                ! How many matrix elements for this jatom
    INTEGER, intent(in)    :: jatom, ngau, lomax
    ! locals
    INTEGER :: jatom_read, ih, ll, mm
    ! reading in thir order  : tuu,tdd,tud,tdu,tuu21,tuu12,tuu22,tud21,tdu12
    ! which is equivalent to : t11,t22,t12,t21, t31,  t13,  t33,  t32,  t23
    ! The potential is stored in this order
    !Vts[:,:]=[[Vtuu,   Vtud,  Vtuu12],
    !          [Vtdu,   Vtdd,  Vtdu12],
    !          [Vtuu21, Vtud21,Vtuu22]]
    !
    ! Vtuu,   Vtud,   Vtuu12
    ! Vtdu,   Vtdd,   Vtdu12
    ! Vtuu21, Vtud21, Vtuu22
    !
    Vts(:,:,:)=0.d0
    read(ifhn,*) jatom_read
    if (jatom_read.ne.jatom) then
       WRITE(6,*) "ERROR: reading case.nsh, jatom_read != jatom"
       WRITE(6,*) 'jatom_read=', jatom_read, 'jatom=', jatom
    endif
    DO ih=1,ngau
       ! first 4 complex numbers are tuu,tdd,tud,tdu
       read(ifhn,'(6i3,4e19.12,/,6(3x),4e19.12)',end=669) lpv(ih),ll,lv(ih),mpv(ih),mm,mv(ih),Vts(1,1,ih),Vts(2,2,ih),Vts(1,2,ih),Vts(2,1,ih)
       if (ll.eq.0) EXIT
       ! next 5 complex numbers are tuu21, tuu12, tuu22, tud21, tdu12
       IF ((LPV(IH).LE.LOMAX).OR.(LV(IH).LE.LOMAX)) read(ifhn,'(6(3x),6e19.12,/,6(3x),4e19.12)',end=669) Vts(3,1,ih),Vts(1,3,ih),Vts(3,3,ih),Vts(3,2,ih),Vts(2,3,ih)
    ENDDO
669 continue
    ihmx=ih-1
  END SUBROUTINE read_V_nsh

  Function CheckFilesExist(filename_V_nsh, filename_V_vns) result(exist)
    CHARACTER*180, intent(in)  :: filename_V_nsh
    CHARACTER*180, intent(in)  :: filename_V_vns
    LOGICAL :: exist
    ! locals
    INTEGER :: ios
    open(9902, file=filename_V_nsh, status='old', err=88)
    read(9902,'(//)',IOSTAT=ios)
    close(9902)
    if (ios.ne.0) goto 88
    open(9919, file=filename_V_vns, status='old', err=88)
    read(9919,'(//)',IOSTAT=ios)
    close(9919)
    if (ios.ne.0) goto 88
    exist=.True.
    return
88  continue
    exist=.False.
    return 
  END Function CheckFilesExist
  
End module readPotential



  !if(mode.eq.'FOR ')   force=.true.
!!!! NOTES: NCOM=121
  ! nr=jri(jatom)
  ! lmmax
  !use param, ONLY: ncom, nrad, lmax2
  !use defs,  ONLY: PI
  !use struk, ONLY: jri, dx, ndif
  !use xa,    ONLY: lm, R, rholm
  !integer, intent(in) :: jatom, lmmax
  !real*8, intent(out) :: fvdrho(0:3,ndif),Vlm(nrad,ncom+1)
  !
  !COMMON /POTNLC/ VR(NRAD)
  ! nr=jri(jatom)

