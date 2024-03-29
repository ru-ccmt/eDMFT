! @Copyright 2007 Kristjan Haule
module Qlist
  IMPLICIT NONE
  INTEGER :: nQ, SCALE
  INTEGER, allocatable :: Qp(:,:)
  REAL*8, allocatable :: wghq(:)
CONTAINS
  SUBROUTINE Qlist__Init__(filename)
    IMPLICIT NONE
    CHARACTER*200, intent(in) :: filename
    !
    INTEGER, PARAMETER :: fh_Qlist=92
    INTEGER :: ios, i
    CHARACTER*10 :: KNAME
    
    open(fh_Qlist, FILE=filename, STATUS='old', IOSTAT=ios)
    
    !READ(fh_Qlist,*) nQ, SCALE
    DO i=1,1000000
       READ (fh_Qlist,'(A10)',IOSTAT=ios) KNAME
       IF (KNAME .EQ. 'END       ') then
          nQ=i-1
          EXIT
       ENDIF
       IF (ios.ne.0) then
          print *, 'ERROR in reading ', TRIM(filename)
          EXIT
       ENDIF
    ENDDO
    ALLOCATE( Qp(4,nQ), wghq(nQ) )
    close(fh_Qlist)
    
    !print *, 'Qlist size=', nQ
    
    !do i=1,nQ
    !   READ(fh_Qlist,*) Qp(1,i), Qp(2,i), Qp(3,i)
    !enddo
    open(fh_Qlist, FILE=filename, STATUS='old', IOSTAT=ios)
    DO i=1,nQ
       READ (fh_Qlist,'(A10,4I5,3F5.2,A3)',IOSTAT=ios) KNAME, Qp(1,i), Qp(2,i), Qp(3,i), Qp(4,i), wghq(i)
       !print *, i, Qp(:,i)
    ENDDO

    
    close(fh_Qlist)
  END SUBROUTINE Qlist__Init__

  SUBROUTINE Qlist__Destruct__()
    DEALLOCATE( Qp, wghq )
  END SUBROUTINE Qlist__Destruct__
end module Qlist

module Klist
  IMPLICIT NONE
  INTEGER :: nkp, SCALE
  INTEGER, allocatable :: kvec(:,:)
  REAL*8, allocatable  :: wgh(:)
  REAL*8 :: BSI(3,3)
  INTEGER, allocatable :: k_m_q(:,:)
  INTEGER, allocatable :: ind1(:,:,:)
CONTAINS

  SUBROUTINE Klist__Init__(filename)
    IMPLICIT NONE
    CHARACTER*200, intent(in) :: filename
    !
    INTEGER, PARAMETER :: fh_klist=93
    INTEGER :: ios, ik1, ik2, ik, iq, i, kin(3), kq(3)
    CHARACTER*10 :: KNAME
    
    open(fh_klist, FILE=filename, STATUS='old', IOSTAT=ios)
    DO ik=1,1000000
       READ (fh_klist,'(A10)',IOSTAT=ios) KNAME
       IF (KNAME .EQ. 'END       ') then
          nkp=ik-1
          EXIT
       ENDIF
       IF (ios.ne.0) then
          print *, 'ERROR in reading ', TRIM(filename)
          EXIT
       ENDIF
    ENDDO
    close(fh_klist)

    ALLOCATE(kvec(4,nkp), wgh(nkp) )
    
    open(fh_klist, FILE=filename, STATUS='old', IOSTAT=ios)
    DO ik=1,nkp
       READ (fh_klist,'(A10,4I5,3F5.2,A3)',IOSTAT=ios) KNAME, kvec(1,ik), kvec(2,ik), kvec(3,ik), kvec(4,ik), wgh(ik)
       !print *, 'KNAME=', KNAME, 'kv=', kvec(1,ik), kvec(2,ik), kvec(3,ik), kvec(4,ik)
    ENDDO

    READ (fh_klist,'(A10)',IOSTAT=ios) KNAME
    READ (fh_klist, *, IOSTAT=ios) (BSI(i,1),i=1,3)
    READ (fh_klist, *, IOSTAT=ios) (BSI(i,2),i=1,3)
    READ (fh_klist, *, IOSTAT=ios) (BSI(i,3),i=1,3)
    
    close(fh_klist)
    
    IF (ios.ne.0) then
       print *, 'ERROR in reading BS at the end of ', TRIM(filename)
    ENDIF

    CALL dinv(BSI,3)
    
    ALLOCATE( k_m_q(nkp, nkp) )

    SCALE = int(kvec(4,1))

    ALLOCATE( ind1(SCALE+1,SCALE+1,SCALE+1) )
    do ik=1,nkp
       CALL k_index__(kin, kvec(:3,ik))
       ind1(kin(1),kin(2),kin(3)) = ik;
    end do
    
    do ik=1,nkp
       do iq=1,nkp
          kq = kvec(:3,ik)-kvec(:3,iq)
          k_m_q(ik,iq) = k_index(kq)
       end do
    end do
    
  END SUBROUTINE Klist__Init__

  SUBROUTINE Klist__Destruct__()
    DEALLOCATE(kvec, wgh)
    DEALLOCATE( k_m_q )
    DEALLOCATE( ind1 )
  END SUBROUTINE Klist__Destruct__
    
  SUBROUTINE k_index__(kin, kq)
    IMPLICIT NONE
    INTEGER, intent(out) :: kin(3)
    INTEGER, intent(in) :: kq(3)
    !
    REAL*8 :: ki(3), kq_(3)
    INTEGER:: i
    !
    kq_(:) = kq(:)
    ki = matmul(BSI,kq_)
    do i=1,3
       kin(i) = MODULO( int(ki(i)+SCALE), SCALE )+1
    enddo
    
  END SUBROUTINE k_index__
  
  INTEGER FUNCTION k_index(kq)
    IMPLICIT NONE
    INTEGER, intent(in) :: kq(3)
    !
    INTEGER :: kin(3)
    CALL k_index__(kin, kq)
    k_index=ind1(kin(1),kin(2),kin(3))
    return
  END FUNCTION k_index
  
end module Klist
