subroutine dsyr2m (uplo,transa,transb,m,n, alpha,a,lda,b,ldb,beta,c,ldc)
  !
  !  Purpose
  !  =======
  !
  !  dsyr2m performs one of the matrix-matrix operations
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
  double precision alpha,beta
  double precision a(lda,*), b(ldb,*), c(ldc,m)
  
  integer nblock
  parameter (nblock=128)
  double precision one,zero
  parameter (one=1.0d0, zero=0.0d0)
  
  integer j,ii,jj,imax,jmax
  logical upper,atrans,btrans
  
  logical lsame
  external lsame
  
  upper=lsame(uplo, 'U')
  atrans=lsame(transa, 'T')
  btrans=lsame(transb, 'T')
  
  if(.not.upper .and. .not.atrans .and. btrans) then
     do ii=1,m,nblock
        imax=min(m,ii+nblock-1)
        do j=ii,imax
           call dgemv('n', imax-j+1, n, alpha, a(j,1), lda, b(j,1), ldb, beta, c(j,j), 1)
        end do
        if(imax.lt.m) call dgemm('n','t',m-imax,imax-ii+1,n, alpha,a(imax+1,1),lda,b(ii,1),ldb, beta,c(imax+1,ii),ldc)
     end do
  else if(.not.upper .and. atrans .and. .not.btrans) then
     do ii=1,m,nblock
        imax=min(m,ii+nblock-1)
        do j=ii,imax
           call dgemv('t', n, imax-j+1, alpha, a(1,j), lda, b(1,j), 1, beta, c(j,j), 1)
        end do
        if(imax.lt.m) &
             call dgemm('t','n',m-imax,imax-ii+1,n,alpha,a(1,imax+1),lda,b(1,ii),ldb,beta,c(imax+1,ii),ldc)
     end do
  else if(upper .and. atrans .and. .not.btrans) then
     do jj=1,m,nblock
        jmax=min(m,jj+nblock-1)
        do j=jj,jmax
           call dgemv('t', n, j-jj+1, alpha, a(1,jj), lda,b(1,j), 1, beta, c(jj,j), 1)
        end do
        if(jj.ne.1) call dgemm('t','n',jj-1,jmax-jj+1,n,alpha,a(1,1),lda,b(1,jj),ldb,beta,c(1,jj),ldc)
     end do
  else
     stop 'wrong parameter in dsyr2m'
  end if
end subroutine dsyr2m

