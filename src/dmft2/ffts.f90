module fftw3_omp
  ! This is interface to FFTW3. It works for OMP and non-OMP case.
  use, intrinsic :: iso_c_binding
  !$ use omp_lib
  include 'fftw3.f03'
  complex(C_DOUBLE_COMPLEX), pointer :: fft(:,:,:)
  type(C_PTR) :: pfft
  type(C_PTR) :: plan
CONTAINS

  SUBROUTINE fft_init_step1(M1,M2,M3)
    IMPLICIT NONE
    INTEGER, intent(in) :: M1, M2, M3
    INTEGER :: void, nthreads, id
    pfft  = fftw_alloc_complex(int(M1 * M2 * M3, C_SIZE_T))
    call c_f_pointer(pfft, fft, [M1,M2,M3])
  END SUBROUTINE fft_init_step1
  SUBROUTINE fft_init_step2(M1, M2, M3, isig)
    IMPLICIT NONE
    INTEGER, intent(in) :: M1, M2, M3, isig
    INTEGER :: void, nthreads, id
    !$omp parallel private(id)
    !$ id = omp_get_thread_num()
    !$omp barrier
    !$ if ( id == 0 ) then
    !$   nthreads = omp_get_num_threads()
    !   !$   write (*,*) 'There are', nthreads, 'threads'
    !$ end if
    !$omp end parallel
    !$ void = fftw_init_threads()
    !$ call fftw_plan_with_nthreads(nthreads);
    if (isig.eq.1) then
       plan = fftw_plan_dft_3d(M3,M2,M1, fft, fft, FFTW_FORWARD,FFTW_ESTIMATE)
    else
       plan = fftw_plan_dft_3d(M3,M2,M1, fft, fft, FFTW_BACKWARD,FFTW_ESTIMATE)
    endif
  END SUBROUTINE fft_init_step2
  
  SUBROUTINE fft_init(M1, M2, M3, isig)
    IMPLICIT NONE
    INTEGER, intent(in) :: M1, M2, M3, isig
    CALL fft_init_step1(M1,M2,M3)
    CALL fft_init_step2(M1, M2, M3, isig)
  END SUBROUTINE fft_init

  SUBROUTINE fft_run()
    IMPLICIT NONE
    call fftw_execute_dft(plan, fft, fft )
  END SUBROUTINE fft_run
  
  SUBROUTINE fft_fini_step1()
    IMPLICIT NONE
    call fftw_free(pfft)
  END SUBROUTINE fft_fini_step1
  SUBROUTINE fft_fini_step2()
    IMPLICIT NONE
    call fftw_destroy_plan(plan)
  END SUBROUTINE fft_fini_step2
  SUBROUTINE fft_fini()
    IMPLICIT NONE
    CALL fft_fini_step1()
    CALL fft_fini_step2()
  END SUBROUTINE fft_fini
end module fftw3_omp

module fftold
  ! This is used for old FFT which does not use external library
  COMPLEX*16, allocatable :: fft(:,:,:)
  COMPLEX*16, allocatable :: cwork(:)
  COMPLEX*16, allocatable :: dwork(:)
  INTEGER :: isig, M1, M2, M3, ierr, csize_fft, dsize_fft
CONTAINS
  SUBROUTINE fft_init(M1_, M2_, M3_, isig_)
    IMPLICIT NONE
    INTEGER, intent(in) :: M1_, M2_, M3_, isig_
    ! locals
    isig = isig_
    M1 = M1_
    M2 = M2_
    M3 = M3_
    dsize_fft = M1+M2+M3
    csize_fft = 4*(M1+M2+M3)+15
    allocate( cwork(csize_fft), dwork(dsize_fft) )
    allocate( fft(M1,M2,M3) )
  END SUBROUTINE fft_init
  
  SUBROUTINE fft_run()
    IMPLICIT NONE
    call C3FFT(M1,M2,M3,fft,M1,M2,isig,cwork,csize_fft,dwork,dsize_fft,ierr)
  END SUBROUTINE fft_run
  
  SUBROUTINE fft_fini()
    IMPLICIT NONE
    deallocate( fft )
    deallocate( cwork, dwork )
  END SUBROUTINE fft_fini
end module fftold

SUBROUTINE fft1set(nk,N1,N2,N3,A,fft,k3,nmat)
  !
  ! Setup fft-fields for the eigenvectors A(K)
  !
  IMPLICIT NONE
  INTEGER, intent(in)     :: nk,N1,N2,N3
  COMPLEX*16, intent(in)  :: A(nk)
  INTEGER, intent(in)     :: nmat
  INTEGER, intent(in)     :: k3(3,nmat)
  COMPLEX*16, intent(out) :: FFT(N1,N2,N3)
  ! locals
  INTEGER                 :: i, i1, i2, i3
  !
  fft(:,:,:)=(0.0d0,0.0d0)
  DO i=1,nk                                                      
     i1 = k3(1,i)
     i2 = k3(2,i)
     i3 = k3(3,i)
     IF(I1.LT.0) i1 = i1+N1
     IF(I2.LT.0) i2 = i2+N2
     IF(I3.LT.0) i3 = i3+N3
     fft(i1+1,i2+1,i3+1)=A(i)
  ENDDO
END SUBROUTINE fft1set

SUBROUTINE fft2set(nk,N1,N2,N3,A,fft,ii,k3,bk3,BR1,nmat)
  !
  ! Sets the fft-fields for (k+K)*A(K), where A(K) is the eigenvector
  !
  IMPLICIT NONE
  INTEGER, intent(in)    :: nk,N1,N2,N3
  COMPLEX*16, intent(in) :: A(nk)
  INTEGER, intent(in)    :: ii, nmat
  INTEGER, intent(in)    :: k3(3,nmat)
  REAL*8, intent(in)     :: bk3(3,nmat)
  REAL*8, intent(in)     :: BR1(3,3)
  COMPLEX*16, intent(out):: FFT(N1,N2,N3)
  ! locals
  INTEGER                :: i, i1, i2, i3
  REAL*8                 :: kpK_xyz
  !
  fft(:,:,:)=(0.0d0,0.0d0)
  DO i=1,nk
     i1 = k3(1,i)
     i2 = k3(2,i)
     i3 = k3(3,i)
     IF(i1.LT.0) i1=i1+N1
     IF(i2.LT.0) i2=i2+N2
     IF(i3.LT.0) i3=i3+N3
     kpK_xyz = dot_product(BR1(ii,:),bk3(:,i))
     fft(i1+1,i2+1,i3+1) = kpK_xyz*A(i)
  ENDDO
END SUBROUTINE fft2set

SUBROUTINE fftsumup(N1,N2,N3,w,fft,sumfft)
  !
  ! Accumulates fft fields
  !
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: sumfft(N1,N2,N3)
  INTEGER, intent(in)      :: N1,N2,N3
  REAL*8, intent(in)       :: w
  COMPLEX*16, intent(in)   :: fft(N1,N2,N3)
  ! locals
  INTEGER    :: i1,i2,i3
  do i3=1,N3
     do i2=1,N2
        do i1=1,N1
           sumfft(i1,i2,i3)=sumfft(i1,i2,i3)+ w*DBLE(FFT(i1,i2,i3)*conjg(FFT(i1,i2,i3)))
        enddo
     enddo
  enddo
  return
end SUBROUTINE fftsumup

SUBROUTINE fftget(nwave,N1,N2,N3,rho1,fft,kmax,nsym,kzz,iz,tau,Qcomplex)
  !
  !..... copies result from fft data to rho(k) output
  !
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: rho1(*)
  INTEGER, intent(in)     :: nwave, nsym
  INTEGER, intent(in)     :: N1,N2,N3
  COMPLEX*16, intent(in)  :: fft(N1,N2,N3)
  INTEGER, intent(in)     :: kmax(3)
  INTEGER, intent(in)     :: kzz(3,nwave)
  INTEGER, intent(in)     :: iz(3,3,nsym)  ! iord==nsym
  REAL*8, intent(in)      :: tau(3,nsym)
  LOGICAL, intent(in)     :: Qcomplex
  ! locals
  INTEGER     :: NST,STG(3,NSYM),IND(NSYM)
  INTEGER     :: j,jj,index,i1,i2,i3
  COMPLEX*16  :: TAUP(NSYM)
  taup=0.0
  index=0
  do j=1,nwave
     call STERN(kzz(:,j),NST,STG,TAUP, iz, tau, nsym, Qcomplex) ! calculates the star of momentum kzz
     do jj=1,NST               ! over all star members
        index=index+1
        i1=stg(1,jj)
        i2=stg(2,jj)
        i3=stg(3,jj)
        IF(IABS(i1).GT.2*kmax(1)) EXIT 
        IF(IABS(i2).GT.2*kmax(2)) EXIT 
        IF(IABS(i3).GT.2*kmax(3)) EXIT 
        IF(i1.LT.0) i1=i1+N1
        IF(i2.LT.0) i2=i2+N2
        IF(i3.LT.0) i3=i3+N3
        rho1(index)=FFT(i1+1,i2+1,i3+1)
     enddo
  enddo
  return
END SUBROUTINE fftget
