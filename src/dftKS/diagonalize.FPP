#ifdef _REAL
#define tpcast(X) dble(X)  !! typecasting in fortran
#define XSWAP DSWAP
#define Ttype REAL*8
#define tsize 8
#define routine_name 'DSYGVX'
#else
#define tpcast(X) X
#define XSWAP ZSWAP
#define Ttype COMPLEX*16
#define tsize 16
#define routine_name 'ZHEGVX'
#endif
SUBROUTINE Diagonalize(Elow, Eup, NE, HS, HSDIAG, Z, EIGVAL, HSROWS, LCOLZ, Qprint)
  use omp_lib
  IMPLICIT NONE
  !..................................................................
  ! The interface between matrix set up in LAPW1 and LAPACK (BLAS) matrix subroutines
  !
  !   Solve the secular equations H*Z = EIGVAL*S*Z
  !   H and S are stored in the array HS,
  !   H is stored below the diagonal, whereas S is stored above
  !   the diagonal. The diagonal of H is stored in the 
  !   one dimensional array HSDIAG.
  !   HSROWS is the dimension of the secular matrices
  !   NE is returned as the number of Eigenvalues between Elow and Eup
  !..................................................................
  REAL*8,  intent(in) :: Elow, Eup
  INTEGER, intent(in) :: HSROWS, LCOLZ
  INTEGER, intent(out):: NE
  Ttype, intent(inout):: HS(HSROWS,HSROWS)
  Ttype, intent(inout):: HSDIAG(HSROWS)
  Ttype,   intent(out):: Z(HSROWS,LCOLZ)
  REAL*8,  intent(out):: eigval(HSROWS)
  LOGICAL, intent(in) :: Qprint
  !   Locals
  Ttype,   allocatable :: Olap(:,:)
  Ttype,   allocatable :: work(:)
  INTEGER, allocatable :: ifail(:), iwork(:)
#ifndef _REAL  
  REAL*8,  allocatable :: rwork(:)
#endif  
  REAL*8             :: t1, t2, t1w, t2w, abstol, Flop
  INTEGER            :: i, j, info, lwork, ierr
  CHARACTER*67       :: ERRMSG
  LOGICAL            :: out, out2
  
  out =.false.
  out2=.false.
  !out=.true.
  if (out) then
     open(86,FILE='_eigvals_',status='unknown')
     open(89,FILE='_eigvect_',status='unknown')
  endif
  if (out2) then
     open(87,FILE='_h_diag_',status='unknown')
     open(88,FILE='_h_and_o_',status='unknown')
     !
     write(87,*) HSrows!,Nvaa,Nlo
     write(88,'(4e23.15)') HS
     write(87,'(4e23.15)') HSdiag
  endif

#ifndef _REAL  
  !  conjugate entire matrix. When we computed H,O we stored it transposed. Now we conjugate, to get correct O & H.
  HS(:,:) = conjg(HS(:,:))
  HSDIAG(:) = conjg(HSDIAG)
#endif

  
  allocate( Olap(HSROWS, HSROWS), STAT=IERR )
  if(IERR.ne.0) THEN 
     write(6,*) 'not enough memory: ', IERR
     stop 'not enough memory' 
  end if
  Olap(:,:)=0.d0
  do i=1,HSROWS
     Olap(i,1:i) = HS(1:i,i)  ! By transposing here, we effectively conjugated overlap
     HS(i,i) = HSDIAG(i)      ! Diagonal of Ham is set to HS
  enddo  ! Now Olap and HS are lower-diagonal matrices of Overlap and Hamiltonian
  
  allocate( iwork(5*HSROWS), ifail(HSROWS), work(2), STAT=IERR )
  if(IERR.ne.0) THEN 
     write(6,*) 'not enough memory: ', IERR
     stop 'not enough memory' 
  end if
  
  abstol=-1.d0 ! set best tolerance
  ! First we do a query for the best size of lwork
  lwork=-1  ! Set lwork to -1 for query
#ifdef _REAL  
  call DSYGVX(1, 'V', 'V', 'L', HSROWS, HS, HSROWS, Olap, HSROWS, Elow, Eup, 0, 0, abstol, NE, EIGVAL, Z, HSROWS, WORK, lwork, IWORK, ifail, info)
#else
  allocate( rwork(7*HSROWS) )
  call ZHEGVX(1, 'V', 'V', 'L', HSROWS, HS, HSROWS, Olap, HSROWS, Elow, Eup, 0, 0, abstol, NE, EIGVAL, Z, HSROWS, WORK, lwork, RWORK, IWORK, ifail, info)
#endif
  lwork = int(work(1)+0.1) ! This is the best lwork
  !print *, 'lwork=', work(1), lwork
  deallocate( work )
  allocate(work(lwork), STAT=IERR )
  if(IERR.ne.0) THEN 
     write(6,*) 'not enough memory: ', IERR
     stop 'not enough memory' 
  end if
  ! Real diagonalization
  t1w = omp_get_wtime ( )
  call cpu_time ( t1 )
#ifdef _REAL  
  call DSYGVX(1, 'V', 'V', 'L', HSROWS, HS, HSROWS, Olap, HSROWS, Elow, Eup, 0, 0, abstol, NE, EIGVAL, Z, HSROWS, WORK, lwork, IWORK, ifail, info)
#else  
  call ZHEGVX(1, 'V', 'V', 'L', HSROWS, HS, HSROWS, Olap, HSROWS, Elow, Eup, 0, 0, abstol, NE, EIGVAL, Z, HSROWS, WORK, lwork, RWORK, IWORK, ifail, info)
  deallocate( rwork )
#endif
  t2w = omp_get_wtime ( )
  call cpu_time ( t2 )
  
  deallocate(iwork, work)
  deallocate( Olap )
  
  if (INFO.NE.0) THEN
     write(6,*) routine_name, 'INFO = ', INFO
     write(99,*) routine_name, 'INFO = ', INFO
     CALL OUTERR(routine_name,'diagonalize(LAPACK) failed.')
     STOP 'diagonalize - Error'
  ENDIF
  IF (NE>LCOLZ) THEN
     CALL OUTERR('diagonalize','More than NUME eigenvalues.')
     WRITE (ERRMSG,'(A,i4,A,i4,A)') 'Increase NUME from',LCOLZ,'to',NE,'   OR'
     CALL OUTERR('diagonalize',ERRMSG)
     WRITE(ERRMSG, '(A)') 'reduce E2 in *.klist.(see *.output1)'
     CALL OUTERR('diagonalize',ERRMSG)
     STOP 'diagonalize - Error'
  ENDIF
  
  if (Qprint) then
     write(6,*) NE, ' Eigenvalues computed '
     write(6,*)  'info=', info
  endif
  deallocate( ifail )
  
  Flop = 4.d0*dble(HSROWS)**3* (tsize/8.d0)**2/3.d0
  if (Qprint) WRITE (6,'(1X,A,A,1x,f10.3,2x,f10.3,2x,f10.2,A)') routine_name,'(Compute eigenvalues (CPU)) :', (t2-t1), (t2w-t1w), (1E-9*Flop/(t2-t1)),' Gflops'

  ! Below is printing for debuggig
  if(out) then
     do i=1,ne
        write(86,*) eigval(i)
     enddo
     
     write(89,'(A)',advance='no') '# basis   '
     do i=1,ne
#ifdef _REAL        
        write(89,'(A5,I4,14X)',advance='no') 'ibnd=',i
#else
        write(89,'(A7,I4,12X)',advance='no') 'RE:ibnd=',i
        write(89,'(A7,I4,12X)',advance='no') 'Im:ibnd=',i
#endif        
     enddo
     write(89,*)
     do j=1,hsrows
        write(89,'(I4)',advance='no') j
        do i=1,ne
#ifdef _REAL           
           write(89,'(f23.15)',advance='no') Z(j,i)
#else
           write(89,'(2f23.15)',advance='no') Z(j,i)
#endif           
        enddo
        write(89,*)
     enddo
     !
     close(86)
     close(89)
  endif
  if (out2) then
     close(87)
     close(88)
  endif
END SUBROUTINE Diagonalize
