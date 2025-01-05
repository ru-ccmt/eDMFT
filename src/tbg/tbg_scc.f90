
subroutine inv2(B, Am)
  IMPLICIT NONE
  !! Performs a direct calculation of the inverse of a 2×2 matrix.
  complex*16, intent(in) :: Am(2,2)   !! Matrix
  complex*16, intent(out):: B(2,2)   !! Inverse matrix
  complex*16             :: detinv
  ! Calculate the inverse determinant of the matrix
  detinv = 1/(Am(1,1)*Am(2,2) - Am(1,2)*Am(2,1))
  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * Am(2,2)
  B(2,1) = -detinv * Am(2,1)
  B(1,2) = -detinv * Am(1,2)
  B(2,2) = +detinv * Am(1,1)
end subroutine inv2

subroutine inv3(B, A)
  IMPLICIT NONE
  !! Performs a direct calculation of the inverse of a 3×3 matrix.
  complex*16, intent(in) :: A(3,3)   !! Matrix
  complex*16, intent(out):: B(3,3)   !! Inverse matrix
  complex*16             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
       - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
       + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end subroutine inv3

subroutine inv4(B, A)
  IMPLICIT NONE
  !! Performs a direct calculation of the inverse of a 4×4 matrix.
  complex*16, intent(in) :: A(4,4)   !! Matrix
  complex*16, intent(out):: B(4,4)   !! Inverse matrix
  complex*16             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = &
       1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

  ! Calculate the inverse of the matrix
  B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
  B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
  B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
  B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
  B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
  B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
  B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
  B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
  B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
  B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
  B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
  B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
  B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
  B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
  B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
  B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
end subroutine inv4

subroutine inv(Ainv, Am, n)
  IMPLICIT NONE
  complex*16, intent(in)  :: Am(n,n)
  complex*16, intent(out) :: Ainv(n,n)
  INTEGER, intent(in)     :: n
  !
  complex*16, dimension(n) :: work  ! work array for LAPACK
  integer,    dimension(n) :: ipiv   ! pivot indices
  integer :: info
  ! External procedures defined in LAPACK
  external ZGETRF
  external ZGETRI
  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = Am
  !n = size(Am,1)
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call ZGETRF(n, n, Ainv, n, ipiv, info)
  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if
  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call ZGETRI(n, Ainv, n, ipiv, work, n, info)
  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end subroutine inv


subroutine Cmp_local_G(gc, dos, Nd0, UAA, Sig, wn, eks, wki, beta, imag_axis, Nw, Nband, Nfunc, Nsite)
  IMPLICIT NONE
  INTEGER,    intent(in) :: Nw, Nband, Nfunc, Nsite
  COMPLEX*16, intent(in) :: wn(Nw)
  REAL*8,     intent(in) :: eks(Nband), wki, beta
  COMPLEX*16, intent(in) :: UAA(Nfunc,Nband,Nsite)
  COMPLEX*16, intent(in) :: Sig(Nfunc,Nw,Nsite)
  COMPLEX*16, intent(out):: gc(Nfunc,Nfunc,Nw,Nsite)
  REAL*8,     intent(out):: dos(Nw), Nd0
  LOGICAL,    intent(in) :: imag_axis
  !f2py integer intent(hide), depend(wn)   :: Nw    = len(wn)
  !f2py integer intent(hide), depend(eks)  :: Nband = len(eks)
  !f2py integer intent(hide), depend(UAA)  :: Nfunc = shape(UAA,0)
  !f2py integer intent(hide), depend(UAA)  :: Nsite = shape(UAA,2)
  !
  COMPLEX*16, allocatable :: tmp(:,:), USigU(:,:), gk_1(:,:), gk(:,:)
  REAL*8, allocatable :: eps(:), tdos(:)
  REAL*8  :: dsum, mu, x
  INTEGER :: iw, ifc, in, iAA
  REAL*8, parameter :: pi = 3.1415926535897932D0
  allocate( tmp(Nfunc,Nband), USigU(Nband,Nband), gk_1(Nband,Nband), gk(Nband,Nband), eps(Nband), tdos(Nw) )
  !$OMP  PARALLEL DO SHARED(gc,dos,Sig,UAA,wn,eks,Nw,Nfunc,Nband,eps)&
  !$OMP& PRIVATE(iw,iAA,tmp,ifc,USigU,gk_1,gk,dsum,in)&
  !$OMP& SCHEDULE(STATIC)
  do iw=1,Nw
     USigU=0.0
     do iAA=1,Nsite
        do ifc=1,Nfunc
           tmp(ifc,:) = Sig(ifc,iw,iAA) * UAA(ifc,:,iAA)
        enddo
        ! USigU = UAA.H * Sig * UAA = UAA.H(ib1,ifc) * Sig_AUU(ifc,ib2)
        ! USigU(:,:) += matmul(transpose(conjg(UAA)), tmp)
        call zgemm('C','N', Nband, Nband, Nfunc, (1.d0,0.d0), UAA(:,:,iAA), Nfunc, tmp, Nfunc, (1.d0,0.d0), USigU, Nband)
     enddo
     
     gk_1 =  0.0
     do in=1,Nband
        gk_1(in,in) = wn(iw)-eks(in)
     enddo

     gk_1(:,:) = gk_1(:,:) - USigU(:,:)
     gk=0.0
     call inv(gk, gk_1, Nband)

     if (imag_axis) then
        dsum = 0.0
        do in=1,Nband
           dsum = dsum + dble(gk(in,in))
        enddo
        dos(iw) = dsum
        if (iw .eq. Nw) then
           do in=1,Nband
              eps(in) = eks(in)+dble(USigU(in,in))  ! the large omega static limit
           enddo
        endif
     else
        dsum = 0.0
        do in=1,Nband
           dsum = dsum - aimag(gk(in,in))/pi
        enddo
        dos(iw) = dsum
     endif
     
     do iAA=1,Nsite
        ! UAA * gk * UAA.H
        ! tmp = matmul(UAA,gk)
        call zgemm('N','N', Nfunc, Nband, Nband, (1.d0,0.d0), UAA(:,:,iAA), Nfunc, gk, Nband, (0.d0,0.d0), tmp, Nfunc)
        !gc(:,:,iw) = matmul(tmp, transpose(conjg(UAA)))
        call zgemm('N','C', Nfunc, Nfunc, Nband, (1.d0,0.d0), tmp, Nfunc, UAA(:,:,iAA), Nfunc, (0.d0,0.d0), gc(:,:,iw,iAA), Nfunc)
     enddo
  enddo
  !$OMP END PARALLEL DO
  gc = gc * wki
  
  deallocate( tmp, USigU, gk_1, gk )      

  if (imag_axis) then
     do iw=1,Nw
        dsum = 0.0
        do in=1,Nband
           dsum = dsum + dble(1.0/(wn(iw)-eps(in)))
        enddo
        tdos(iw) = dsum
     enddo
     dos(:) = dos(:) - tdos(:)
     mu = dble(wn(1))
     Nd0=0.0
     do in=1,Nband
        ! computes the fermi function 
        x = (eps(in)-mu)*beta
        if (x.lt.-100) then
           Nd0 = Nd0 + 1.0
        else if (x .lt. 100) then
           Nd0 = Nd0 + 1.0/(1.0 + exp(x))
        endif
     enddo
     Nd0 = Nd0 * wki
  endif
  
  dos = dos * wki
  
  deallocate( eps, tdos )
end subroutine Cmp_local_G

subroutine Cmp_eigenvals(zek, UAA, Sig, eks,  Nw, Nband, Nfunc, Nsite)
  IMPLICIT NONE
  INTEGER,    intent(in) :: Nw, Nband, Nfunc, Nsite
  REAL*8,     intent(in) :: eks(Nband)
  COMPLEX*16, intent(in) :: UAA(Nfunc,Nband,Nsite)
  COMPLEX*16, intent(in) :: Sig(Nfunc,Nw,Nsite)
  COMPLEX*16, intent(out):: zek(Nband,Nw)
  !f2py integer intent(hide), depend(Sig)  :: Nw    = shape(Sig,1)
  !f2py integer intent(hide), depend(eks)  :: Nband = len(eks)
  !f2py integer intent(hide), depend(UAA)  :: Nfunc = shape(UAA,0)
  !f2py integer intent(hide), depend(UAA)  :: Nsite = shape(UAA,2)
  !
  COMPLEX*16, allocatable :: tmp(:,:), USigU(:,:), ham(:,:)
  INTEGER :: iw, ifc, in
  !REAL*8, parameter :: pi = 3.1415926535897932D0
  COMPLEX*16 evl(1,1),evr(1,1)
  COMPLEX*16 :: cworkvec(4*Nband) ! Work array for zgeev                                                                                                           
  REAL*8     :: rworkvec(4*Nband) ! Work array for zgeev                                                                                                           
  INTEGER    :: ierr, iAA
  INTEGER    :: idxarr(Nband)
  
  allocate( tmp(Nfunc,Nband), USigU(Nband,Nband), ham(Nband,Nband) )

  !$OMP  PARALLEL DO SHARED(zek,Sig,UAA,eks,Nw,Nfunc,Nband)&
  !$OMP& PRIVATE(iw,tmp,ifc,USigU,ham,evl,evr,cworkvec,rworkvec,in,iAA,ierr,idxarr)&
  !$OMP& SCHEDULE(STATIC)
  do iw=1,Nw
     USigU=0.0
     do iAA=1,Nsite
        do ifc=1,Nfunc
           tmp(ifc,:) = Sig(ifc,iw,iAA) * UAA(ifc,:,iAA)
        enddo
        ! USigU = UAA.H * Sig * UAA = UAA.H(ib1,ifc) * Sig_AUU(ifc,ib2)
        ! USigU(:,:) = matmul(transpose(conjg(UAA)), tmp)
        call zgemm('C','N', Nband, Nband, Nfunc, (1.d0,0.d0), UAA(:,:,iAA), Nfunc, tmp, Nfunc, (1.d0,0.d0), USigU, Nband)
     enddo
     
     ham =  0.0
     do in=1,Nband
        ham(in,in) = eks(in)
     enddo
     ham(:,:) = ham(:,:) + USigU(:,:)
     
     CALL ZGEEV('N','N',Nband,ham,Nband,zek(:,iw),evl,Nband,evr,Nband,cworkvec,4*Nband,rworkvec,ierr)
     IF(ierr.NE.0) THEN
        WRITE(6,*) 'Error code of zgeev ', ierr
     ENDIF
     
     !========= Sorting acording to real parts
     CALL eig_order_real_part(zek(:,iw), idxarr, Nband)
     CALL permute_eigensystem(idxarr, zek(:,iw), Nband)
     
  enddo
  !$OMP END PARALLEL DO
  deallocate( tmp, USigU, ham )      
end subroutine Cmp_eigenvals

subroutine Cmp_dens(Nd, Nd0, mu, wn, zek, eks, wk, beta, Nw, Nband, Nik)
  IMPLICIT NONE
  INTEGER,    intent(in) :: Nw, Nband, Nik
  REAL*8,     intent(in) :: mu, wk(Nik), eks(Nband,Nik), beta
  REAL*8,     intent(in) :: wn(Nw)
  COMPLEX*16, intent(in) :: zek(Nband,Nw,Nik)
  REAL*8,     intent(out):: Nd(Nw), Nd0
  !f2py integer intent(hide), depend(wn)   :: Nw    = len(wn)
  !f2py integer intent(hide), depend(zek)  :: Nband = shape(zek,0)
  !f2py integer intent(hide), depend(zek)  :: Nik   = shape(zek,2)
  !
  COMPLEX*16 :: zsum, wsum, om
  REAL*8     :: wghk, sm_wk, dsum, nn, x
  INTEGER    :: iw, ik, n
  COMPLEX*16 :: i
  i = (0.0D+00,1.0D+00)
  sm_wk = sum(wk)
  !print *, 'sm_wk=', sm_wk, 'wghk=', wk
  Nd=0
  do iw=1,Nw
     zsum = 0.0
     wsum = 0.0
     om = i*wn(iw)
     do ik=1,Nik
        wghk = wk(ik)
        do n=1,Nband
           zsum = zsum + wghk/(om + mu - zek(n,iw,ik))
           wsum = wsum + wghk/(om + mu - dble(zek(n,Nw,ik)) )
        enddo
     enddo
     Nd(iw) = dble((zsum-wsum)/sm_wk)
  enddo
  
  Nd0=0 ! sum of fermi functions
  do ik=1,Nik
     wghk = wk(ik)
     nn=0
     do n=1,Nband
        !x = (eks(n,ik)-mu)*beta
        x = (dble(zek(n,Nw,ik))-mu)*beta
        if (x.lt.-100) then
           Nd0 = Nd0 + wghk
        else if (x .lt. 100) then
           Nd0 = Nd0 + wghk/(1.0 + exp(x))
        endif
     enddo
     !print *, ik, nn, eks(:,ik)
  enddo
  Nd0 = Nd0/sm_wk
end subroutine Cmp_dens

SUBROUTINE eig_order_real_part(ev, idxarr, ndim)
  IMPLICIT NONE 
!!!-----------------------------------------------------------------!!!
!!! This routine sorts complex eigenvalues of a matrix according to !!!
!!! its real parts with the smallest in the first slot and reorders !!!
!!! the matrices of left (row) and right (column) eigenvectors in a !!!
!!! corresponding manner.                                           !!!
!!!-----------------------------------------------------------------!!!
  !---------- Passed variables ----------
  COMPLEX*16, intent(in) :: ev(ndim)         ! Array of eigenvalues
  INTEGER, intent(out)   :: idxarr(ndim)     ! Index array which gives proper order
  INTEGER :: ndim                            ! Dimension of matrices 
  !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
  !---------- Parameters ----------
  REAL*8, PARAMETER :: maxval = 1000.d0
  !---------- Local variables ----------
  LOGICAL, ALLOCATABLE :: sorted(:)
  REAL*8,  ALLOCATABLE :: sortonr(:)
  INTEGER :: p
  INTEGER :: q
  INTEGER :: idx
  REAL*8  :: min
  !---------- Allocate dynamic memory storage ----------
  ALLOCATE(sortonr(ndim), sorted(ndim))
  !---------- Initialize arrays ----------
  idxarr = 0
  sortonr = DBLE(ev)
  sorted = .FALSE.
  !---------- Create index array for real value ----------
  sorted = .FALSE.
  DO p = 1,ndim
     min = maxval
     DO q = 1,ndim
        IF(.NOT.sorted(q).AND.min.GT.sortonr(q)) THEN
           min = sortonr(q)
           idx = q
        ENDIF
     ENDDO
     idxarr(p) = idx
     sorted(idx) = .TRUE.
  ENDDO
  DEALLOCATE(sortonr, sorted)
  RETURN
END SUBROUTINE eig_order_real_part

SUBROUTINE permute_eigensystem(idxarr, ev, ndim)
  IMPLICIT NONE
  !---------- Passed variables ----------
  INTEGER, intent(in)       :: idxarr(ndim)     ! Index array which gives proper order
  COMPLEX*16, intent(inout) :: ev(ndim)         ! Array of eigenvalues
  INTEGER :: ndim                               ! Dimension of matrices 
  !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
  !---------- Local variables ------------------
  INTEGER :: p
  COMPLEX*16, ALLOCATABLE :: eval(:)
  ALLOCATE(eval(ndim) ) 
  !---------- Permute the eigenvalues ----------
  DO p = 1,ndim
     eval(p) = ev(idxarr(p))
  ENDDO
  ev = eval
  !---------- Deallocate dynamic memory storage ----------
  DEALLOCATE(eval) 
END SUBROUTINE permute_eigensystem
