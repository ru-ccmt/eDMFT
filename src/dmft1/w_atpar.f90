module w_atpar
  IMPLICIT NONE
  !*******************************************
  ! atpar changes the following common blocks:
  REAL*8, allocatable :: w_alo(:,:,:,:,:)
  INTEGER, allocatable:: w_nlo(:), w_nlov(:), w_nlon(:), w_ilo(:,:)
  LOGICAL, allocatable:: w_lapw(:,:)!, w_loor(:,:)
  REAL*8, allocatable :: w_ri_mat(:,:,:,:,:)
  REAL*8, allocatable :: w_P(:,:,:,:), w_DP(:,:,:,:)
  INTEGER, allocatable:: w_jatom(:)
  REAL*8, allocatable :: w_FJ(:,:,:)
  REAL*8, allocatable :: w_DFJ(:,:,:)
  ! Wave functions for projector and current solution within MT-sphere
  REAL*8, allocatable :: w_rfk(:,:,:,:)
  REAL*8, allocatable :: el_store(:,:,:), elo_store(:,:,:,:)
  INTEGER :: nnlo
  CONTAINS

  SUBROUTINE w_allocate(maxucase,nat,iso2)
    use param
    use case
    IMPLICIT NONE
    integer :: maxucase, nat, iso2
    allocate (el_store(0:lmax2,nat,iso2),elo_store(0:lomax,1:nloat,nat,iso2))
    allocate( w_alo(0:lomax,nloat,nrf,2,maxucase) )
    allocate( w_nlo(maxucase), w_nlov(maxucase), w_nlon(maxucase), w_ilo(0:lomax,maxucase) )
    !allocate( w_loor(0:lomax,maxucase) )
    allocate( w_lapw(0:lmax2,maxucase) )
    allocate( w_ri_mat(nrf,nrf,0:lmax2,2,maxucase) )
    !allocate( w_rf1(NRAD,0:LMAX2,2,nrf,maxucase), w_rf2(NRAD,0:LMAX2,2,nrf,maxucase) ) ! It seems to me that this one does not need to be saved!
    allocate( w_P(0:LMAX2,2,nrf,maxucase), w_DP(0:LMAX2,2,nrf,maxucase) )              ! Needs to be saved!
    !allocate( w_lfirst(maxucase) )
    allocate( w_jatom(maxucase))    
    allocate( w_FJ(0:lmax2,nmat,maxucase) )
    allocate( w_DFJ(0:lmax2,nmat,maxucase) )

  END SUBROUTINE w_allocate
  
  SUBROUTINE w_allocate_rfk(n_al_ucase)
    use param
    use case
    IMPLICIT NONE
    integer :: n_al_ucase
    allocate( w_rfk(nrad,2,nrf,n_al_ucase) )
    w_rfk=0.0
  END SUBROUTINE w_allocate_rfk
  
  SUBROUTINE w_deallocate
    IMPLICIT NONE
    deallocate( w_alo )
    deallocate( w_nlo, w_nlov, w_nlon, w_ilo ) 
    !deallocate( w_loor )
    deallocate( w_lapw )
    deallocate( w_ri_mat )
    !deallocate( w_rf1, w_rf2 )
    deallocate( w_P, w_DP )
    !deallocate( w_lfirst )
    deallocate( w_jatom )
    deallocate( w_FJ )
    deallocate( w_DFJ )
    deallocate( el_store, elo_store )
  END SUBROUTINE w_deallocate
  
  SUBROUTINE w_deallocate_rfk()
    IMPLICIT NONE
    deallocate( w_rfk )
  END SUBROUTINE w_deallocate_rfk

end module w_atpar
