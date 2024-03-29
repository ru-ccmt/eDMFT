MODULE structure
  CHARACTER                :: title*80,lattic*4  ! structure title and lattice type
  REAL*8                   :: alat(3),alpha(3),Vol ! lattice constants, angles, vol
  LOGICAL                  :: ortho        ! orthogonal unit cell
  REAL*8                   :: BR1(3,3)     ! conventional reciprocal
  REAL*8                   :: BR2(3,3)     ! primitive reciprocal vector b_i
  INTEGER                  :: ndf          ! number of all atoms in the unit cell
  INTEGER,ALLOCATABLE      :: iatnr(:)
  INTEGER,ALLOCATABLE      :: mult(:)      ! number of equivalent atoms
  INTEGER,ALLOCATABLE      :: isplit(:)
  INTEGER,ALLOCATABLE      :: jri(:)       ! jri: number of radial (logarithmic) mesh points
  REAL*8,ALLOCATABLE       :: pos(:,:)     ! positions of all atoms in the unit cell
  REAL*8,ALLOCATABLE       :: Rmt(:)       ! mufin-thin sphere radius
  REAL*8,ALLOCATABLE       :: v(:)         ! relative volume of this sphere Vol_Sphere/Vol_Unit_cel
  REAL*8,ALLOCATABLE       :: rotloc(:,:,:)! local rotation matrix from global to local coordinate system
  REAL*8,ALLOCATABLE       :: r0(:),dx(:)  ! smallest point in logarithmic difference for the mesh
  CHARACTER*10,ALLOCATABLE :: aname(:)     ! name of the atom
  REAL*8, allocatable      :: ZZ(:)        ! nucleous charge
  CHARACTER*4              :: irel, cform  ! relativistic or not.
  REAL*8                   :: rot_spin_quantization(3,3)
  REAL*8,ALLOCATABLE       :: rotij(:,:,:) ! rotate atom (including equivalent ones) according to the symmetry. For atoms with only one equiv. atom this is unity.
  REAL*8,ALLOCATABLE       :: tauij(:,:)   ! shift for nonsymorfic groups
  ! Space group information
  INTEGER                  :: IORD         ! number of symmetry operations of space group
  INTEGER, allocatable     :: IMAT(:,:,:)  ! matrix representation of (space group) symmetry operation 
  REAL*8,  allocatable     :: TAU(:,:)     ! non-primitive translation vector for symmetry operation
contains

  SUBROUTINE ReadStructure(fh_str,nat,rel,lxdos,ERRMSG,info)
    IMPLICIT NONE
    INTEGER, intent(in)      :: fh_str ! should be 20
    INTEGER, intent(out)     :: nat
    LOGICAL, intent(out)     :: rel
    INTEGER, intent(inout)   :: lxdos
    INTEGER, intent(out)     :: info
    CHARACTER*67, intent(out):: ERRMSG
    !
    REAL*8,PARAMETER          :: test = 1.D-12
    REAL*8, allocatable       :: tpos(:,:)
    INTEGER :: index, jatom, i, j, m
    REAL*8  :: PI
    !
    READ(fh_str,1000) title
    READ(fh_str,1010) lattic,nat,cform,irel
    
    allocate( rmt(nat),v(nat),iatnr(nat),mult(nat),isplit(nat) )
    allocate( rotloc(3,3,nat) )
    allocate( r0(nat), dx(nat), jri(nat) )
    allocate( aname(nat) )
    allocate( ZZ(nat) )
    allocate ( tpos(3,48*nat) ) ! temporary allocate large array

    !.....READ IN LATTICE CONSTANTS                                         
    READ(fh_str,1020) alat(1),alat(2),alat(3),alpha(1),alpha(2),alpha(3)
    if (abs(alpha(1)).LT.test) alpha(1)=90.0d0
    if (abs(alpha(2)).LT.test) alpha(2)=90.0d0
    if (abs(alpha(3)).LT.test) alpha(3)=90.0d0

    PI=ACOS(-1.0D0)
    alpha(1) = alpha(1)*PI/180.0D0
    alpha(2) = alpha(2)*PI/180.0D0
    alpha(3) = alpha(3)*PI/180.0D0

    rel=.FALSE.
    IF(IREL.EQ.'RELA') rel=.TRUE.
    !IF(IREL.EQ.'NREL') REL=.FALSE.
    !
    !  read crystal-structure (atompositions, symmetry-parameters,muffin-tin radius, ...)
    !  'INDEX' counts all atoms in the unit cell,
    !  'JATOM' counts only the non-equivalent atoms
    index = 0                                                          
    DO jatom = 1,NAT                                               
       index = index+1
       READ(fh_str,1030) iatnr(jatom),( tpos(j,index),j=1,3), mult(jatom),isplit(jatom)
       if(isplit(jatom).eq.99)    lxdos=   3
       if(isplit(jatom).eq.88)    lxdos=   3
       if ( mult(jatom) .EQ. 0 .OR. mult(jatom).gt.(48*nat)) THEN !...illegal number of equivalent atoms
          WRITE (6,6000) jatom, index, mult(jatom)
          info = 3
          ERRMSG = 'In structure.f90 MULT (number of eq. atoms) has illegal value'
          return
       ENDIF
       DO m=1,mult(jatom)-1
          index = index+1                                            
          READ(fh_str,1031) iatnr(jatom),( tpos(j,index),j=1,3)
       ENDDO
       READ(fh_str,1050) aname(jatom), jri(jatom), r0(jatom), rmt(jatom), ZZ(jatom)
       dx(jatom) = log(rmt(jatom)/r0(jatom)) / (jri(jatom)-1)
       rmt(jatom) = r0(jatom)*exp(dx(jatom)*(jri(jatom)-1) )
       READ(fh_str,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)
    ENDDO
    ndf = sum(mult)
    allocate( pos(3,ndf) )
    pos(:,:) = tpos(:,:ndf)
    deallocate( tpos )
    
    call SymmRot(rotloc,nat)   !     Correct rotloc
    info=0
    ERRMSG=''
    !print *, 'ERRMSG=', ERRMSG
    !print *, 'info=', info
1000 FORMAT(A80)
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
!1040 FORMAT(///,3X,'ERROR IN DMFT2 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
1051 FORMAT(20X,3F10.8)
6000 FORMAT(///,3X,'ERROR IN DMFT2 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE ReadStructure


  SUBROUTINE ReadStructureGroup(fh_str, nat)
    IMPLICIT NONE
    INTEGER, intent(in)      :: fh_str, nat
    ! locals
    CHARACTER*80  ::  WSPACE
    LOGICAL       :: there
    integer       :: i, j, j1, j2, jj, jatom
    REAL*8        :: postest(3)
    !
    READ (fh_str,'(I4)') IORD
    call spaceGroupInit(IORD)
    DO  j=1,IORD
       READ (20,50101) ((IMAT(j1,j2,j),j1=1,3),TAU(j2,j),j2=1,3)
    ENDDO
50101 FORMAT (3 (3I2,F11.8,/))
    !
    !     Are high-precision atoms available ?
    read(fh_str,'(A)',err=1999,end=1999)WSPACE
    there=.false.
    if(WSPACE(1:17).eq.'Precise positions') there=.true.

    if(there)then
       do jj=1,sum(mult)
          read(fh_str,*)postest(1:3)
          do i=1,3
             if(abs(modulo(postest(i),1.d0)-modulo(pos(i,jj),1.d0)).gt.1.d-8) then
                print*, 'Precise positions overwritten',postest(i),pos(i,jj) 
                !goto 1999   !KH: This jump looks really strange!!!!
             endif
          enddo
          pos(1:3,jj)=postest(1:3)
       enddo
    endif
1999 continue
  END SUBROUTINE ReadStructureGroup

  
  SUBROUTINE WriteInfoStructure(fh_stdout, nat)
    INTEGER, intent(in) :: fh_stdout  ! should be 6
    INTEGER, intent(in) :: nat
    WRITE(fh_stdout,800)                                                      
    WRITE(fh_stdout,805)  title
    WRITE(fh_stdout,810)  lattic
    WRITE(fh_stdout,820)  alat(1),alat(2),alat(3)                                            
    WRITE(fh_stdout,840)  nat
    WRITE(fh_stdout,850)  irel
800 FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ','I N F O R M A T I O N',/,30X,50(1H-),//)                  
805 FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)                             
810 FORMAT(3X,'LATTICE',22X,'= ',A4)                                  
820 FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)                 
840 FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)                   
850 FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)                    
  END SUBROUTINE WriteInfoStructure

  SUBROUTINE DeallocateStructure()
    IMPLICIT NONE
    deallocate( rmt, v, iatnr, mult, isplit )
    deallocate( rotloc )
    deallocate( r0, dx, jri )
    deallocate( aname )
    if (allocated(ZZ)) deallocate( ZZ )
    deallocate( pos )
    if (allocated(imat)) deallocate( imat )
    if (allocated(tau)) deallocate( tau )
  END SUBROUTINE DeallocateStructure

  subroutine spaceGroupInit(nsym)
    IMPLICIT NONE
    integer, intent(in) :: nsym
    allocate( IMAT(3,3,NSYM) )
    allocate( TAU(3,NSYM) )
    IMAT(1:3,1:3,1:nsym)=0
    TAU(1:3,1:nsym)=0.0D0
  end subroutine spaceGroupInit
  
  SUBROUTINE SymmRot(rotloc,NAT)
    IMPLICIT NONE
    !
    !---------------------------------------------------------------------  
    !  Patch up symmetry issues of rotation matrices
    INTEGER, intent(in) :: nat
    DOUBLE PRECISION, intent(inout) :: rotloc(3,3,NAT)
    ! locals
    DOUBLE PRECISION :: rold(3,3,NAT), RotTest, test1
    INTEGER :: IRtest, jatom, j1, j2
    rold=rotloc
    !
    !     Test various integer combinations -- this will catch most things
    DO IRTest = 2,6
       RotTest = 1.D0/sqrt(DBLE(IRTEST))
       do jatom=1,NAT
          do j1=1,3
             do j2=1,3
                test1=abs(rotloc(j2,j1,jatom))-RotTest
                if(abs(test1).lt.1d-6)rotloc(j2,j1,jatom)=sign(RotTest,rotloc(j2,j1,jatom))
             enddo
          enddo
       enddo
    ENDDO
    RotTest = sqrt(0.75D0)
    do jatom=1,NAT
       do j1=1,3
          do j2=1,3
             test1=abs(rotloc(j2,j1,jatom))-RotTest
             if(abs(test1).lt.1d-6)rotloc(j2,j1,jatom)=sign(RotTest,rotloc(j2,j1,jatom))
          enddo
       enddo
    enddo

    RotTest = sqrt(2.D0/3.D0)
    do jatom=1,NAT
       do j1=1,3
          do j2=1,3
             test1=abs(rotloc(j2,j1,jatom))-RotTest
             if(abs(test1).lt.1d-6)rotloc(j2,j1,jatom)=sign(RotTest,rotloc(j2,j1,jatom))
          enddo
       enddo
    enddo
    !     do Jatom=1,NAT
    !       write(76,*)rotloc(1:3,1:3,jatom)
    !       write(75,*)rold  (1:3,1:3,jatom)
    !     enddo
    !
    return
  END SUBROUTINE SymmRot
  
END MODULE structure
