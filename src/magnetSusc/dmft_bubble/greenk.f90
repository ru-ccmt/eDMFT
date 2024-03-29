! @Copyright 2007 Kristjan Haule
module greenk
  IMPLICIT NONE
  REAL*8, allocatable :: oml(:)
  COMPLEX*16, allocatable :: gk(:,:,:,:)
  INTEGER :: noml, nkpl, cixdm
CONTAINS
  SUBROUTINE greenk__Init__(filename)
    IMPLICIT NONE
    INTEGER, PARAMETER :: fh_gk=94
    CHARACTER*200 :: filename
    !
    INTEGER :: nsymop, norbitals
    INTEGER :: ios, i, j, ik, isym, iom
    CHARACTER*1 :: ch
    REAL*8 :: x, re, im
    REAL*8, allocatable :: gg2(:)
    COMPLEX*16, PARAMETER  :: IMAG = (0.0D0,1.0D0)
    
    open(fh_gk, FILE=filename, STATUS='old', IOSTAT=ios)
    
    READ(fh_gk,*) ch, nkpl,nsymop,noml,cixdm,norbitals
    READ(fh_gk,*) ! actual positions of atoms

    allocate( gk(noml, cixdm, cixdm, nkpl), oml(noml) )
    
    allocate( gg2(cixdm*cixdm*2) )
    do ik=1,nkpl
       do isym=1,nsymop
          do iom=1,noml
             READ(fh_gk,*) x, (gg2(i),i=1,cixdm*cixdm*2)
             if (ik.eq.1 .and. isym.eq.1) then 
                oml(iom)=x
             endif
             do i=1,cixdm
                do j=1,cixdm
                   gk(iom,i,j,ik) = gk(iom,i,j,ik) + gg2((i-1)*2*cixdm+j*2-1) + gg2((i-1)*2*cixdm+j*2)*IMAG
                enddo
             enddo
          enddo
       enddo
       if (MODULO(ik,10).EQ.0) WRITE(*,'(A,I5)') 'Finished reading k-point=', ik
    enddo
    deallocate( gg2 )
    
    if (nsymop.GT.1) gk = gk * (1./nsymop)
    close(fh_gk)
  END SUBROUTINE greenk__Init__

  SUBROUTINE greenk__Destruct__()
    deallocate( gk, oml )
  END SUBROUTINE greenk__Destruct__
  
end module greenk
