subroutine prod_zher (uplo,transa,transb,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
!
!  Purpose
!  =======
!
!  zher2m performs one of the matrix-matrix operations
!
!     op1(C) := alpha*op2( A )*op2( B ) + beta*op1(C),
!
!  where  op1( X ) means that only the UPPER or only the LOWER
!  part of X is accessed
!
!  and  op2( X ) is one of
!
!     op2( X ) = X   or   op2( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with 
!  op2( A ) an m by n matrix, 
!  op2( B ) an n by m matrix, and 
!  C an m by m matrix.
!
!     uplo: 'U'pper or 'L'ower
!     transa,transb: op2( X ) = X   or   op2( X ) = X', with X=A,B
!     m: number of rows of op2(A) and op2(B) and of
!        rows and columns of C
!     n: number of columns of op2(A) and op2(B)
!     alpha
!     op2(A): mXn matrix
!     lda: leading dimension of a
!     op2(B): mXn matrix
!     ldb: leading dimension of b
!     beta
!     c: symmetric mXm matrix, 
!        stored either above (uplo='U') or below (uplo='L') the
!        diagonal
!     ldc: leading dimension of c
!
  implicit none
  character*1 uplo,transa,transb
  integer m,n,lda,ldb,ldc
  complex*16 alpha,beta
  complex*16 a(lda,*), b(ldb,*), c(ldc,m)
  
  integer nblock
  parameter (nblock=128)
  complex*16 one,zero
  parameter (one=(1.0d0,0.0d0), zero=(0.0d0,0.0d0))
  
  integer j,ii,jj,imax,jmax
  logical upper,atransC,btransC,atransT,btransT,atransN,btransN
  
  logical lsame
  external lsame
  
  upper=lsame(uplo, 'U')
  atransC=lsame(transa, 'C')
  btransC=lsame(transb, 'C')
  atransT=lsame(transa, 'T')
  btransT=lsame(transb, 'T')
  atransN=lsame(transa, 'N')
  btransN=lsame(transb, 'N')
  
  if(.not.upper .and. atransN .and. (btransC.or.btransT)) then
     do ii=1,m,nblock
        imax=min(m,ii+nblock-1)
        do j=ii,imax
           call zgemv('n', imax-j+1, n, alpha, a(j,1), lda,b(j,1), ldb, beta, c(j,j), 1)
        end do
        if(imax.lt.m) call zgemm('n',transb,m-imax,imax-ii+1,n,alpha,a(imax+1,1),lda,b(ii,1),ldb,beta,c(imax+1,ii),ldc)
     end do
  else if(.not.upper .and. (atransC.or.atransT) .and. btransN) then
     do ii=1,m,nblock
        imax=min(m,ii+nblock-1)
        do j=ii,imax
           call zgemv(transa, n, imax-j+1, alpha, a(1,j), lda,b(1,j), 1, beta, c(j,j), 1)
        end do
        if(imax.lt.m) call zgemm(transa,'n',m-imax,imax-ii+1,n,alpha,a(1,imax+1),lda,b(1,ii),ldb,beta,c(imax+1,ii),ldc)
     end do
  else if(upper .and. (atransC.or.atransT) .and. btransN) then
     do jj=1,m,nblock
        jmax=min(m,jj+nblock-1)
        do j=jj,jmax
           call zgemv(transa, n, j-jj+1, alpha, a(1,jj), lda,b(1,j), 1, beta, c(jj,j), 1)
        end do
        if(jj.ne.1) call zgemm(transa,'n',jj-1,jmax-jj+1,n,alpha,a(1,1),lda,b(1,jj),ldb,beta,c(1,jj),ldc)
     end do
  else
     stop 'wrong parameter in zher2m'
  end if
end subroutine prod_zher
