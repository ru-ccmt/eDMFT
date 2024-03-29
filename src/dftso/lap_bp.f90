SUBROUTINE DBP2P( UPLO, N, AP, HB )
  IMPLICIT NONE
  !
  !  DBP2P transforms a block packed stored matrix AP into a packed matrix.
  !
  !     .. Scalar Arguments ..
  CHARACTER*1        UPLO
  INTEGER            N, HB
  REAL*8             AP(*)
  !     
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix AP is stored as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of AP
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of AP
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A. 
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - REAL*8
  !           On entry, AP stores the matrix in packed storage format.
  !           In order to store the (not referenced elements in the 
  !           block diagonal), AP must be at least of size 
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !           On exit, AP stores the matrix in block packed storage.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            IJ2K
  EXTERNAL           LSAME, IJ2K
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            INFO, I, I1, I2, I3, J
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  UPPER = LSAME( UPLO, 'U' )
  !
  INFO = 0
  IF(( .NOT.UPPER).AND. &
       ( .NOT.LSAME( UPLO , 'L' ) ) )THEN
     INFO = 1
  ELSE IF( N  .LT.0 )THEN
     INFO = 2
  ELSE IF( HB .LT.0 ) THEN
     INFO = 4
  END IF
  !
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'DBP2P ', INFO )
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF ( UPPER ) THEN
     I1 = 1
     I2 = 1
     I3 = 0
     DO J = 1,N        ! 20
        DO I = 1,J     ! 10
           AP(I1) = AP(I2)
           I1 = I1 +1
           I2 = I2 +1
        ENDDO          ! 10         CONTINUE
        
        IF (J .GT. ((N/HB)*HB)) THEN
           I3 = N-J
        ELSE
           I3 = MOD(I3+HB-1 , HB)
        END IF
        I2 = I2 + I3
     ENDDO            ! 20         CONTINUE
  ELSE
     I1 = 1
     I2 = 1
     I3 = 0
     DO J = 1,N      ! 40
        DO I = J,N   ! 30
           AP(I1) = AP(I2)
           I1 = I1 +1
           I2 = I2 +1
           
        ENDDO  ! 30               CONTINUE
        I3 = MOD(I3+1 , HB)
        I2 = I2 + I3
     ENDDO ! 40               CONTINUE
  END IF
  !
  !     End of DBP2P
  !
END SUBROUTINE DBP2P

SUBROUTINE DLATD4( UPLO, N, NB, ABP, E, TAU, W, LDW, HB )
  IMPLICIT NONE
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER          UPLO
  INTEGER            LDW, N, NB, HB
  !     ..
  !     .. Array Arguments ..
  DOUBLE PRECISION   ABP( * ), E( * ), TAU( * ), W( LDW, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DLATRD reduces NB rows and columns of a real symmetric matrix A to
  !  symmetric tridiagonal form by an orthogonal similarity
  !  transformation Q' * A * Q, and returns the matrices V and W which are
  !  needed to apply the transformation to the unreduced part of A.
  !
  !  If UPLO = 'U', DLATRD reduces the last NB rows and columns of a
  !  matrix, of which the upper triangle is supplied;
  !
  !  This is an auxiliary routine called by DSYTRD.
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER
  !          Specifies whether the upper or lower triangular part of the
  !          symmetric matrix A is stored:
  !          = 'U': Upper triangular
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.
  !
  !  NB      (input) INTEGER
  !          The number of rows and columns to be reduced.
  !
  !  ABP     (input/output) DOUBLE PRECISION array, 
  !          dimension (1:n*(n+1)/2+(n/hb+1)*hb*(hb-1)/2)
  !          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !          n-by-n upper triangular part of A contains the upper
  !          triangular part of the matrix A, and the strictly lower
  !          triangular part of A is not referenced.  
  !          On exit:
  !          if UPLO = 'U', the last NB columns have been reduced to
  !            tridiagonal form, with the diagonal elements overwriting
  !            the diagonal elements of A; the elements above the diagonal
  !            with the array TAU, represent the orthogonal matrix Q as a
  !            product of elementary reflectors.
  !          See Further Details.
  !
  !  E       (output) DOUBLE PRECISION array, dimension (N-1)
  !          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
  !          elements of the last NB columns of the reduced matrix;
  !
  !  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
  !          The scalar factors of the elementary reflectors, stored in
  !          TAU(n-nb:n-1) if UPLO = 'U'.
  !          See Further Details.
  !
  !  W       (output) DOUBLE PRECISION array, dimension (LDW,NB)
  !          The n-by-nb matrix W required to update the unreduced part
  !          of A.
  !
  !  LDW     (input) INTEGER
  !          The leading dimension of the array W. LDW >= max(1,N).
  !
  !  Further Details
  !  ===============
  !
  !  If UPLO = 'U', the matrix Q is represented as a product of elementary
  !  reflectors
  !
  !     Q = H(n) H(n-1) . . . H(n-nb+1).
  !
  !  Each H(i) has the form
  !
  !     H(i) = I - tau * v * v'
  !
  !  where tau is a real scalar, and v is a real vector with
  !  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
  !  and tau in TAU(i-1).
  !
  !  The elements of the vectors v together form the n-by-nb matrix V
  !  which is needed, with W, to apply the transformation to the unreduced
  !  part of the matrix, using a symmetric rank-2k update of the form:
  !  A := A - V*W' - W*V'.
  !
  !  The contents of A on exit are illustrated by the following example
  !  with n = 5 and nb = 2:
  !
  !    (  a   a   a   v4  v5 )
  !    (      a   a   v4  v5 )  
  !    (          a   1   v5 ) 
  !    (              d   1  )
  !    (                  d  )
  !
  !  where d denotes a diagonal element of the reduced matrix, a denotes
  !  an element of the original matrix that is unchanged, and vi denotes
  !  an element of the vector defining H(i).
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  DOUBLE PRECISION   ZERO, ONE, HALF
  PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
  !     ..
  !     .. Local Scalars ..
  INTEGER            I, IW, LDI, LDJ
  DOUBLE PRECISION   ALPHA
  INTEGER            I1I, I1I1, i1J, IJI, IJJ, II1I, III1, J, JB 
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  DOUBLE PRECISION   DDOT
  EXTERNAL           LSAME, DDOT
  INTEGER            IJ2K
  EXTERNAL           IJ2K
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MIN
  
  !     ..
  !     .. Executable Statements ..
  !
  !     Quick return if possible
  !
  IF( N.LE.0 ) RETURN
  !
  IF( LSAME( UPLO, 'U' ) ) THEN
     !
     !        Reduce last NB columns of upper triangle
     !
     !	 write(*,*) 'DLATD4: N=',N,' NB=',NB
     LDI = N
     DO I = N, N - NB + 1, -1    ! 10
        IW = I - N + NB
        I1I = IJ2K( UPLO, 1, I, N, HB )
        II1I = IJ2K( UPLO, I-1, I, N, HB)
        IF( I.LT.N ) THEN
           I1I1 = IJ2K( UPLO, 1, I+1, N, HB )
           III1 = IJ2K( UPLO, I, I+1, N, HB )
           !
           !              Update A(1:i,i)
           !
           CALL DGEMV( 'No transpose', I, N-I, -ONE, ABP(I1I1),LDI, W( I, IW+1 ), LDW, ONE, ABP(I1I), 1 )
           CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),LDW, ABP(III1), LDI, ONE, ABP(I1I), 1 )
        END IF
        IF( I.GT.1 ) THEN
           !
           !              Generate elementary reflector H(i) to annihilate
           !              A(1:i-2,i)
           !
           CALL DLARFG( I-1, ABP(II1I), ABP(I1I), 1, TAU( I-1 ) )
           E( I-1 ) = ABP(II1I)
           ABP(II1I) = ONE
           !
           !              Compute W(1:i-1,i)
           !
           
           DO J = 1, I-1   ! 4
              W(J, IW) = ZERO

           ENDDO           !  4  CONTINUE
           DO J = 1, I-1, HB  ! 5
              JB = MIN( HB, I-J)
              IF( JB.EQ.HB ) THEN
                 LDJ = J+JB-1
              ELSE
                 LDJ = LDI
              END IF
              I1J = IJ2K( UPLO, 1, J, N, HB )
              IJI = IJ2K( UPLO, J, I, N, HB ) 
              IJJ = IJ2K( UPLO, J, J, N, HB ) 
              IF( J.NE.1) THEN 
                 CALL DGEMV( 'NoTranspose', J-1, JB, ONE, ABP(I1J), LDJ, ABP(IJI), 1, ONE, W(1, IW), 1) 
                 CALL DGEMV( 'Transpose', J-1, JB, ONE, ABP(I1J), LDJ, ABP(I1I), 1, ONE, W(J, IW), 1)
              END IF
              !	         write(*,*) 'HB=',HB,' J=',J,' JB=',JB,' I=',I
              CALL DSYMV( 'U', JB, ONE, ABP(IJJ), LDJ, ABP(IJI), 1, ONE, W( J, IW), 1)
           ENDDO   !5             CONTINUE
           IF( I.LT.N ) THEN
              CALL DGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ),LDW, ABP(I1I), 1, ZERO, W( I+1, IW ), 1 )
              CALL DGEMV( 'No transpose', I-1, N-I, -ONE,ABP(I1I1), LDI, W( I+1, IW ), 1, ONE,W( 1, IW ), 1 )
              CALL DGEMV( 'Transpose', I-1, N-I, ONE, ABP(I1I1),LDI, ABP(I1I), 1, ZERO, W( I+1, IW ), 1 )
              CALL DGEMV( 'No transpose', I-1, N-I, -ONE,W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,W( 1, IW ), 1 )
           END IF
           CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
           ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1, ABP(I1I), 1 )
           CALL DAXPY( I-1, ALPHA, ABP(I1I), 1, W( 1, IW ), 1 )
        END IF
        !
        ! 10      CONTINUE
     ENDDO
  END IF
  !
  RETURN
  !
  !     End of DLATRD
  !
END SUBROUTINE DLATD4

SUBROUTINE DP2BP( UPLO, N, AP, HB )
  IMPLICIT NONE
  !
  !  DP2BP transforms a packed stored matrix AP into a block packed matrix.
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER*1        UPLO
  INTEGER            N, HB
  REAL*8             AP(*)
  !     
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix AP is stored as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of AP
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of AP
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A. 
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - REAL*8
  !           On entry, AP stores the matrix in packed storage format.
  !           In order to store the (not referenced elements in the 
  !           block diagonal), AP must be at least of size 
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !           On exit, AP stores the matrix in block packed storage.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            IJ2K
  EXTERNAL           LSAME, IJ2K
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            INFO, I, I1, I2, I3, J
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  UPPER = LSAME( UPLO, 'U' )
  !
  INFO = 0
  IF(( .NOT.UPPER).AND. &
       ( .NOT.LSAME( UPLO , 'L' ) ) )THEN
     INFO = 1
  ELSE IF( N  .LT.0 )THEN
     INFO = 2
  ELSE IF( HB .LT.0 ) THEN
     INFO = 4
  END IF
  !
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'DP2BP ', INFO )
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF ( UPPER ) THEN
     I1 = N*(N+1)/2
     I2 = IJ2K('U',N,N,N,HB)
     I3 = 0
     DO J = N,1,-1           ! 20
        DO I = J,1,-1        ! 10
           AP(I2) = AP(I1)
           I1 = I1 - 1
           I2 = I2 - 1
        ENDDO  ! 10         CONTINUE
        IF ((J-1) .EQ. ((N/HB)*HB)) I3 = HB-1
        I3 = MOD (I3+1, HB)
        DO I = 0, I3-1     ! 15
           AP(I2-I) = 0.0
        ENDDO              ! 15       CONTINUE
        I2 = I2 - I3
     ENDDO                 ! 20      CONTINUE
  ELSE
     I1 = N*(N+1)/2
     I2 = IJ2K('L',N,N,N,HB)
     I3 = MOD(N-1,HB) 
     DO J = N,1,-1      ! 40
        DO I = N,J,-1   ! 30
           AP(I2) = AP(I1)
           I1 = I1 -1
           I2 = I2 -1
        ENDDO         ! 30             CONTINUE
        DO I=0,I3-1   ! 35
           AP(I2-I)=0.0
        ENDDO         !  35     CONTINUE
        I2 = I2 - I3
        I3 = MOD(I3 - 1+HB, HB)
     ENDDO            !   40    CONTINUE
  END IF
  !
  !     End of DP2BP
  !
END SUBROUTINE DP2BP

SUBROUTINE DPHTR4( N, ABP, HB, INFO )
  IMPLICIT NONE
  !
  !
  !
  !     .. Scalar Arguments ..
  INTEGER            INFO, N, HB
  !     ..
  !     .. Array Arguments ..
  REAL*8             ABP( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DPHTR4 computes an U**T*U factorization of a real symmetric
  !  positive definite matrix A.
  !
  !  This is the hierarchically blocked right-looking (k-variant) of
  !  the algorithm.
  !
  !  Arguments
  !  =========
  !
  !  N       (input) INTEGER
  !          The number of rows and columns of the matrix A.  N >= 0.
  !
  !  ABP     (input/output) REAL array,
  !          dimension (N*(N+1)/2+(N/HB+1)*HB*(HB-1)/2).
  !          Only the block diagonal and the upper triangular part of A
  !          are stored columnwise (block packed storage). All blocks
  !          of the block diagonal are of order HB (the order of the last
  !          one may be less).
  !          On entry, the symmetric matrix A.  Only the upper triangular
  !          part of A is referenced.
  !          On exit, the factor U from the Cholesky factorization.
  !          The matrix U is stored using block packed storage.
  !
  !  HB      (input) INTEGER
  !          The blocksize for the hierarchically blocked algorithm.
  !          HB determines the order of the diagonal blocks. HB > 0.
  !          = 1  : use packed storage format
  !          >= N : use conventional storage format
  !          else : use block packed storage
  !
  !  INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value 
  !          > 0: if INFO = k, the leading minor of order K is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  REAL*8             ONE
  PARAMETER          ( ONE = 1.0D+0 )
  !     ..
  !     .. Local Scalars ..
  INTEGER            K, KK, KKK, KB, KKB, KKKB
  INTEGER            KP, KPEND
  INTEGER            KTMP1, KTMP2, KTMP3
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DSYRK, DTRSM
  !     .. External Functions ..
  INTEGER            IJ2K 
  EXTERNAL           IJ2K
  !     INTEGER            mp_numthreads
  !     EXTERNAL           mp_numthreads
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MIN
  !     ..
  !     .. Executable Statements ..
  !
  INFO = 0
  IF( N.LT.0 ) THEN
     INFO = -1
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DPHTR4', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) RETURN
  IF( HB.EQ.1 ) THEN
     CALL DPPTRF( 'U', N, ABP, INFO )
  ELSEIF ( HB.GE.N ) THEN
     CALL DPOTRF( 'U', N, ABP, N, INFO )
  ELSE
     DO K = 1, N, HB                   ! 10
        KB = MIN( HB, N-K+1 )
        KTMP1 = IJ2K( 'U', K, K, N, HB )
        CALL DPOTRF( 'U', KB, ABP( KTMP1 ), K+KB-1, INFO )
        if ( info.ne.0 ) goto 20
        IF( K+KB.LE.N ) THEN
           KPEND = ((N-K-KB+1)+(HB-1))/HB -1
           DO KP = 0, KPEND             ! 7
              KK = K+KB + KP*HB
              KKB = MIN( HB, N-KK+1 )
              KTMP2 = IJ2K( 'U', K, KK, N, HB )
              CALL DTRSM( 'Left', 'U', 'Transpose', 'Non-unit', KB, KKB, ONE, ABP(KTMP1), K+KB-1, ABP( KTMP2 ), KK+KKB-1 )
           ENDDO  ! 7             CONTINUE
           ! DOACROSS Local (KP,KK,KKB,KKK,KKKB,KTMP1,KTMP2,KTMP3)
           DO KP = 0, KPEND  ! 5
              KK = K+KB + KP*HB
              KKB = MIN( HB, N-KK+1 )
              KTMP1 = IJ2K( 'U', K, KK, N, HB )
              KTMP2 = IJ2K( 'U', KK, KK, N, HB )
              CALL DSYRK( 'U', 'Transpose', KKB, KB, -ONE, ABP( KTMP1 ), KK+KKB-1, ONE, ABP( KTMP2 ), KK+KKB-1  )
              IF( KK+KKB.LE.N ) THEN
                 DO KKK = KK+KKB, N, HB            ! 3
                    KKKB = MIN( HB, N-KKK+1 )
                    KTMP1 = IJ2K( 'U', K, KKK, N, HB )
                    KTMP2 = IJ2K( 'U', K, KK, N, HB )
                    KTMP3 = IJ2K( 'U', KK, KKK, N, HB )
                    CALL DGEMM( 'Transpose', 'No transpose', KKB, KKKB, KB, -ONE, ABP( KTMP2 ), KK+KKB-1, ABP ( KTMP1 ), KKK+KKKB-1, ONE, ABP( KTMP3 ), KKK+KKKB-1 )
                 ENDDO                             ! 3  CONTINUE
              END IF
           ENDDO !    5          CONTINUE
        END IF
        IF( INFO.NE.0 ) GO TO 20
     ENDDO ! 10      CONTINUE
  END IF
  RETURN
20 CONTINUE
  INFO = INFO + K - 1
  RETURN
  !
  !     End of DPHTR4
  !
END SUBROUTINE DPHTR4


SUBROUTINE DPHTRF( UPLO, N, AP, HB, INFO )
  IMPLICIT NONE
  !
  !
  !
  CHARACTER      UPLO
  INTEGER        INFO, N, HB
  DOUBLE         PRECISION AP( * )
  !	
  !       Hierarchically blocked version of DPPTRF
  !
  !	Parameters
  !
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix AP is stored as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of AP
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of AP
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - REAL*8
  !           On entry, AP stores the matrix in packed storage format.
  !           In order to store the (not referenced elements in the
  !           block diagonal), AP must be at least of size
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !           On exit, AP stores the matrix in block packed storage.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !  INFO  - (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order K is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !  =====================================================================
  !           
  Integer            j
  REAL*8             CP(4), FLOP
  logical            UPPER, LSAME
  integer ilaenv
  external ilaenv, LSAME
  !
  !     Test the input parameters.
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( HB.LE.0 ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DPHTRF', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) RETURN
  
  IF( HB.le.1 ) then
     call DPPTRF( UPLO, N, AP, INFO ) 
  ELSE
     
     call CPUTIM(CP(1))
     !
     !      call the subroutines
     !			
     !
     !      convert from LAPACK upper triangular packed storage format
     !                     into a upper triangular block packed matrix
     !
     call dp2bp('U',N,AP,HB)
     call CPUTIM(CP(2))
     !
     !      compute an U**T*U factorization
     !      using the hierarchically blocked k-version of the UTU Algorithm
     !
     call DPHTR4(N,AP,HB,INFO)
     call CPUTIM(CP(3))
     IF (INFO .NE. 0) THEN
        CALL OUTERR('DPHTRF','DPHTR4 aborted unsuccessfully.')
        GOTO 999
     ENDIF
     !
     !      convert back to upper triangular packed storage format
     !
     call dbp2p('U',N,AP,HB)
     call CPUTIM(CP(4))
     !
     !      timing output
     !
     !         DO 30 j = 1, 3
     !           CP(j) = CP(j+1)-CP(j)
     ! 30      CONTINUE
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Packing' , CP(1) , 1E-6*Flop/CP(1)
     !         Flop = (N*DBLE(N+1)*(2*N+1))/6
     !         WRITE (*,1001) 'Cholesky' , CP(2) , 1E-6*Flop/CP(2)
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Unpacking' , CP(3) , 1E-6*Flop/CP(3)
     ! 1001    FORMAT (1X,'DPHTRF(',A,') :',t30,f7.3:t40,f8.3,' Mflops')
     !	 write (*,*) 'INFO = ', INFO
     
  END IF
  RETURN
  !	
  !	Result stored using LAPACK packed format
  !
  !
999 STOP 'LAPW1 - Error'
  !
END SUBROUTINE DPHTRF

SUBROUTINE DSHEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL,HBTRD, NUME, NBVL, INFO )
  !
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER          JOBZ, RANGE, UPLO
  INTEGER            IL, INFO, IU, LDZ, M, N, HBTRD, NUME, NBVL
  DOUBLE PRECISION   ABSTOL, VL, VU
  !     ..
  !     .. Array Arguments ..
  INTEGER            IFAIL( * ), IWORK( * )
  DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DSHEVX computes selected eigenvalues and, optionally, eigenvectors
  !  of a real symmetric matrix A in packed storage, which is converted 
  !  for computation to block packed storage.  
  !  Eigenvalues/vectors
  !  can be selected by specifying either a range of values or a range of
  !  indices for the desired eigenvalues.
  !
  !  Arguments
  !  =========
  !
  !  JOBZ    (input) CHARACTER*1
  !          = 'N':  Compute eigenvalues only;
  !          = 'V':  Compute eigenvalues and eigenvectors.
  !
  !  RANGE   (input) CHARACTER*1
  !          = 'A': all eigenvalues will be found;
  !          = 'V': all eigenvalues in the half-open interval (VL,VU]
  !                 will be found;
  !          = 'I': the IL-th through IU-th eigenvalues will be found.
  !          = 'C': the eigenvalues in the half-open interval (VL,VU]
  !                 will be found, if there are more than M Eigenvalues,
  !                 only the smallest M Eigenvalues are computed;
  !
  !  UPLO    (input) CHARACTER*1
  !          = 'U':  Upper triangle of A is stored;
  !          = 'L':  Lower triangle of A is stored.
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  AP      (input/output) DOUBLE PRECISION array, 
  !          dimension (N*(N+1)/2+(N/HBTRD+1)*HBTRD*(HBTRD-1)/2)
  !          On entry, the upper triangle of the symmetric matrix
  !          A, packed columnwise in a linear array.  The j-th column of A
  !          is stored in the array AP as follows:
  !          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j.
  !
  !          On exit, AP is overwritten by values generated during the
  !          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
  !          and first superdiagonal of the tridiagonal matrix T overwrite
  !          the corresponding elements of A.
  !
  !  VL      (input) DOUBLE PRECISION
  !  VU      (input) DOUBLE PRECISION
  !          If RANGE='V' or RANGE'C', the lower and upper bounds of the 
  !          interval to be searched for eigenvalues. VL < VU.
  !          Not referenced if RANGE = 'A' or 'I'.
  !
  !  IL      (input) INTEGER
  !  IU      (input) INTEGER
  !          If RANGE='I', the indices (in ascending order) of the
  !          smallest and largest eigenvalues to be returned.
  !          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
  !          Not referenced if RANGE = 'A' or 'V'.
  !
  !  ABSTOL  (input) DOUBLE PRECISION
  !          The absolute error tolerance for the eigenvalues.
  !          An approximate eigenvalue is accepted as converged
  !          when it is determined to lie in an interval [a,b]
  !          of width less than or equal to
  !
  !                  ABSTOL + EPS *   max( |a|,|b| ) ,
  !
  !          where EPS is the machine precision.  If ABSTOL is less than
  !          or equal to zero, then  EPS*|T|  will be used in its place,
  !          where |T| is the 1-norm of the tridiagonal matrix obtained
  !          by reducing AP to tridiagonal form.
  !
  !          Eigenvalues will be computed most accurately when ABSTOL is
  !          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
  !          If this routine returns with INFO>0, indicating that some
  !          eigenvectors did not converge, try setting ABSTOL to
  !          2*DLAMCH('S').
  !
  !          See "Computing Small Singular Values of Bidiagonal Matrices
  !          with Guaranteed High Relative Accuracy," by Demmel and
  !          Kahan, LAPACK Working Note #3.
  !
  !  M       (output) INTEGER
  !          The total number of eigenvalues found.  0 <= M <= N.
  !          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1, and
  !          if RANGE = 'C', M <= NUME.
  !
  !  W       (output) DOUBLE PRECISION array, dimension (N)
  !          If INFO = 0, the selected eigenvalues in ascending order.
  !
  !  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,NUME))
  !          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
  !          contain the orthonormal eigenvectors of the matrix A
  !          corresponding to the selected eigenvalues, with the i-th
  !          column of Z holding the eigenvector associated with W(i).
  !          If an eigenvector fails to converge, then that column of Z
  !          contains the latest approximation to the eigenvector, and the
  !          index of the eigenvector is returned in IFAIL.
  !          If JOBZ = 'N', then Z is not referenced.
  !          Note: the user must ensure that at least max(1,M) columns are
  !          supplied in the array Z; if RANGE = 'V', the exact value of M
  !          is not known in advance and an upper bound must be used.
  !
  !  LDZ     (input) INTEGER
  !          The leading dimension of the array Z.  LDZ >= 1, and if
  !          JOBZ = 'V', LDZ >= max(1,N).
  !
  !  WORK    (workspace) DOUBLE PRECISION array, 
  !          dimension (MAX((HBTRD+3),8)*N)
  !
  !  IWORK   (workspace) INTEGER array, dimension (5*N)
  !
  !  IFAIL   (output) INTEGER array, dimension (N)
  !          If JOBZ = 'V', then if INFO = 0, the first M elements of
  !          IFAIL are zero.  If INFO > 0, then IFAIL contains the
  !          indices of the eigenvectors that failed to converge.
  !          If JOBZ = 'N', then IFAIL is not referenced.
  !
  !  HBTRD   (input)  INTEGER.
  !          On entry, HBTRD specifies the block size of the diagonal blocks.
  !          HBTRD must be at least zero.
  !          Unchanged on exit.
  !
  !  NUME    (input) INTEGER
  !          The maximum number of Eigenvectors allowed
  !            (size of array Z)
  !
  !  NBVL    (output) INTEGER
  !          The number of Eigenvalues below the lower bound VL.
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !          > 0:  if INFO = i, then i eigenvectors failed to converge.
  !                Their indices are stored in array IFAIL.
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  DOUBLE PRECISION   ZERO, ONE
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            ALLEIG, INDEIG, VALEIG, WANTZ
  CHARACTER          ORDER
  INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWO, INDTAU, INDWRK, ISCALE, ITMP1, J, JJ, NSPLIT 
  DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
  INTEGER            LLWORK
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  DOUBLE PRECISION   DLAMCH, DLANSP
  EXTERNAL           LSAME, DLAMCH, DLANSP
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DCOPY, DOPGTR, DOPMTR, DSCAL, DSHTRD, DSTEZ2, DSTEIN, DSTEQR, DSTERF, DSWAP, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MIN, SQRT
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  WANTZ = LSAME( JOBZ, 'V' )
  ALLEIG = LSAME( RANGE, 'A' )
  VALEIG = LSAME( RANGE, 'V' ) .OR. LSAME( RANGE, 'C' )
  INDEIG = LSAME( RANGE, 'I' )
  !
  INFO = 0
  IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
     INFO = -1
  ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
     INFO = -2
  ELSE IF( .NOT.( LSAME( UPLO, 'L' ) .OR. LSAME( UPLO, 'U' ) ) ) &
       THEN
     INFO = -3
  ELSE IF( N.LT.0 ) THEN
     INFO = -4
  ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
     INFO = -7
  ELSE IF( INDEIG .AND. IL.LT.1 ) THEN
     INFO = -8
  ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) THEN
     INFO = -9
  ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
     INFO = -14
  END IF
  !
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSHEVX', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  M = 0
  IF( N.EQ.0 ) &
       RETURN
  !
  IF( N.EQ.1 ) THEN
     IF( ALLEIG .OR. INDEIG ) THEN
        M = 1
        W( 1 ) = AP( 1 )
     ELSE
        IF( VL.LT.AP( 1 ) .AND. VU.GE.AP( 1 ) ) THEN
           M = 1
           W( 1 ) = AP( 1 )
        END IF
     END IF
     IF( WANTZ ) Z( 1, 1 ) = ONE
     RETURN
  END IF
  !
  !     Get machine constants.
  !
  SAFMIN = DLAMCH( 'Safe minimum' )
  EPS = DLAMCH( 'Precision' )
  SMLNUM = SAFMIN / EPS
  BIGNUM = ONE / SMLNUM
  RMIN = SQRT( SMLNUM )
  RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
  !
  !     Scale matrix to allowable range, if necessary.
  !
  ISCALE = 0
  ABSTLL = ABSTOL
  IF( VALEIG ) THEN
     VLL = VL
     VUU = VU
  END IF
  ANRM = DLANSP( 'M', UPLO, N, AP, WORK )
  IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
     ISCALE = 1
     SIGMA = RMIN / ANRM
  ELSE IF( ANRM.GT.RMAX ) THEN
     ISCALE = 1
     SIGMA = RMAX / ANRM
  END IF
  IF( ISCALE.EQ.1 ) THEN
     CALL DSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
     IF( ABSTOL.GT.0 ) &
          ABSTLL = ABSTOL*SIGMA
     IF( VALEIG ) THEN
        VLL = VL*SIGMA
        VUU = VU*SIGMA
     END IF
  END IF
  !      write(*,*) 'Scaling (DLANSY + DSCAL)'
  !      call printtime
  !
  !     Call DSHTRD to reduce symmetric packed matrix to tridiagonal form.
  !
  INDTAU = 1
  INDE = INDTAU + N
  INDD = INDE + N
  INDWRK = INDD + N
  LLWORK = N*HBTRD
  !      tarray(1)=dtime(tarray)
  CALL DSHTRD( UPLO, N, AP, WORK( INDD ), WORK( INDE ),WORK( INDTAU ), WORK( INDWRK ), LLWORK, HBTRD,IINFO )
  IF( IINFO.NE.0 ) write(*,*) 'IINFO=',IINFO
  !      write(*,*) 'Reduce to tridiagonal form (DSHTRD)'
  !      call printtime
  !
  !     If all eigenvalues are desired and ABSTOL is less than or equal
  !     to zero, then call DSTERF or DOPGTR and SSTEQR.  If this fails
  !     for some eigenvalue, then try DSTEZ2.
  !
  IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.EQ.N ) ) .AND.( ABSTOL.LE.ZERO ) ) THEN
     CALL DCOPY( N, WORK( INDD ), 1, W, 1 )
     INDEE = INDWRK + 2*N
     IF( .NOT.WANTZ ) THEN
        CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
        CALL DSTERF( N, W, WORK( INDEE ), INFO )
     ELSE
        CALL DOPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO )
        CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
        CALL DSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ,WORK( INDWRK ), INFO )
        IF( INFO.EQ.0 ) THEN
           DO I = 1, N     !10
              IFAIL( I ) = 0
           ENDDO ! 10            CONTINUE
        END IF
     END IF
     IF( INFO.EQ.0 ) THEN
        M = N
        GO TO 20
     END IF
     INFO = 0
  END IF
  !
  !     Otherwise, call DSTEZ2 and, if eigenvectors are desired, SSTEIN.
  !
  IF( WANTZ ) THEN
     ORDER = 'B'
  ELSE
     ORDER = 'E'
  END IF
  INDIBL = 1
  INDISP = INDIBL + N
  INDIWO = INDISP + N
  IF(WANTZ .AND. (RANGE.EQ.'C') ) THEN
     M=NUME
  END IF
  CALL DSTEZ2( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL,WORK( INDD ), WORK( INDE ), M, NSPLIT, W,IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ),IWORK( INDIWO ), NUME, NBVL, INFO )
  !       write(*,*) M, ' Eigenvectors computed'
  !       write(*,*) NBVL, ' Eigenvalues below lower bound'
  !      write(*,*) 'compute eigenvalues (DSTEZ2)'
  !      call printtime
  !
  IF( WANTZ ) THEN
     CALL DSTEIN( N, WORK( INDD ), WORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO )
     !      write(*,*) 'compute eigenvectors (DSTEIN)'
     !      call printtime
     !
     !        Apply orthogonal matrix used in reduction to tridiagonal
     !        form to eigenvectors returned by DSTEIN.
     !
     !         INDWKN = INDE
     !         LLWRKN = LWORK - INDWKN + 1
     !         tarray(1)=dtime(tarray)
     !         CALL DOHMT4( 'L', UPLO, 'N', N, M, A, WORK( INDTAU ), Z,
     !     $                LDZ, WORK( INDWKN ), LLWRKN, IINFO, LDA )
     CALL DOPMTR( 'L', UPLO, 'N', N, M, AP, WORK( INDTAU ), Z, LDZ,WORK( INDWRK ), INFO )
     !      write(*,*) 'backtransform eigenvectors (DOHMT4)'
     !      call printtime
  END IF
  !
  !     If matrix was scaled, then rescale eigenvalues appropriately.
  !
20 CONTINUE
  IF( ISCALE.EQ.1 ) THEN
     IF( INFO.EQ.0 ) THEN
        IMAX = M
     ELSE
        IMAX = INFO - 1
     END IF
     CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
  END IF
  !      write(*,*) 'scaling (DSCAL)'
  !      call printtime
  !
  !     If eigenvalues are not in order, then sort them, along with
  !     eigenvectors.
  !
  IF( WANTZ ) THEN
     DO J = 1, M - 1                 ! 40
        I = 0
        TMP1 = W( J )
        DO JJ = J + 1, M             ! 30
           IF( W( JJ ).LT.TMP1 ) THEN
              I = JJ
              TMP1 = W( JJ )
           END IF
        ENDDO  ! 30         CONTINUE
        !
        IF( I.NE.0 ) THEN
           ITMP1 = IWORK( INDIBL+I-1 )
           W( I ) = W( J )
           IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
           W( J ) = TMP1
           IWORK( INDIBL+J-1 ) = ITMP1
           CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
           IF( INFO.NE.0 ) THEN
              ITMP1 = IFAIL( I )
              IFAIL( I ) = IFAIL( J )
              IFAIL( J ) = ITMP1
           END IF
        END IF
     ENDDO ! 40      CONTINUE
  END IF
  !      write(*,*) 'sorting (DSWAP)'
  !      call printtime
  !
  RETURN
  !
  !     End of DSHEVX
  !
END SUBROUTINE DSHEVX

SUBROUTINE DSHGS4( ITYPE, UPLO, N, ABP, BBP, HB, INFO )
  IMPLICIT NONE
  ! 
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER          UPLO
  INTEGER            INFO, ITYPE, N, HB
  !     ..
  !     .. Array Arguments ..
  DOUBLE PRECISION   ABP( * ), BBP( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DSHGS4 reduces a real symmetric-definite generalized eigenproblem
  !  to standard form. The arrays ABP and BBP are stored in block packed 
  !  mode. All blocks of the diagonal are of order HB (the order of the 
  !  last one may be less).
  !
  !  If ITYPE = 1, the problem is A*x = lambda*B*x,
  !  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
  !
  !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
  !  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
  !
  !  Only ITYPE = 1 is implemented.
  !
  !  B must have been previously factorized as U**T*U or L*L**T by DPHTRF.
  !
  !  Arguments
  !  =========
  !
  !  ITYPE   (input) INTEGER
  !          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
  !          = 2 or 3: compute U*A*U**T or L**T*A*L, not implemented.
  !
  !  UPLO    (input) CHARACTER
  !          = 'U':  Upper triangle of A is stored and B is factored as
  !                  U**T*U;
  !          = 'L':  Lower triangle of A is stored and B is factored as
  !                  L*L**T, not implemented.
  !
  !  N       (input) INTEGER
  !          The order of the matrices A and B.  N >= 0.
  !
  !  ABP     (input/output) DOUBLE PRECISION array.
  !          Dimension (N*(N+1)/2+(N/HB+1)*HB*(HB-1)/2).
  !          On entry, the symmetric matrix A.  
  !
  !          On exit, if INFO = 0, the transformed matrix, stored in the
  !          same format as A.
  !
  !  BBP     (input) DOUBLE PRECISION array.
  !          Dimension (N*(N+1)/2+(N/HB+1)*HB*(HB-1)/2).
  !          The triangular factor from the Cholesky factorization of B,
  !          as returned by DPHTRF.
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  DOUBLE PRECISION   ONE, HALF
  PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            K, KB, KK, KKB, KKK, KKKB
  INTEGER            IKK, IKK2, IK2K2, IKK3, IK3K2
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DSYGS2, DSYMM, DSYR2K, DTRMM, DTRSM, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            ILAENV, IJ2K
  EXTERNAL           LSAME, ILAENV
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( ITYPE.NE.1 ) THEN
     INFO = -1
  ELSE IF( .NOT.UPPER ) THEN
     INFO = -2
  ELSE IF( N.LT.0 ) THEN
     INFO = -3
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSHGS4', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) RETURN
  !
  IF( HB.LE.1 ) THEN
     CALL DSPGST( ITYPE, UPLO, N, ABP, BBP, INFO)
  ELSE IF( HB.GE.N ) THEN
     CALL DSYGST( ITYPE, UPLO, N, ABP, N, BBP, N, INFO )
  ELSE
     !
     !     Determine the block size for this environment.
     !
     !        NB = ILAENV( 1, 'DSYGST', UPLO, N, -1, -1, -1 )
     !
     !        Use hierarchically blocked code
     !
     IF( ITYPE.EQ.1 ) THEN
        IF( UPPER ) THEN
           !
           !              Compute inv(U')*A*inv(U)
           !
           DO K = 1, N, HB      ! 10
              KB = MIN( N-K+1, HB )
              !
              !                 Update the upper triangle of A(k:n,k:n)
              !
              IKK = IJ2K( UPLO, K, K, N, HB )
              !		  write (*,*) 'dsygst: ikk,k,kb=',ikk,k,kb 
              CALL DSYGST( ITYPE, UPLO, KB, ABP( IKK ), K+KB-1, &
                   BBP( IKK ), K+KB-1, INFO )
              IF( INFO.NE.0) GO TO 20
              IF( K+KB.LE.N ) THEN
                 DO KK = K+KB, N, HB   ! 9
                    KKB = MIN( N-KK+1, HB )
                    IKK2 = IJ2K( UPLO, K, KK, N, HB )
                    IK2K2 = IJ2K( UPLO, KK, KK, N, HB )
                    CALL DTRSM( 'Left', UPLO, 'Trans', 'Non-unit', KB, KKB, ONE, BBP( IKK ), K+KB-1,ABP( IKK2 ), KK+KKB-1 )
                    CALL DSYMM( 'Left', UPLO, KB, KKB, -HALF,ABP( IKK ), K+KB-1, BBP( IKK2 ),KK+KKB-1,ONE, ABP( IKK2 ), KK+KKB-1 )
                 ENDDO ! 9                CONTINUE
                 !
                 !                       dsyr2k block-wise
                 !
                 DO KK = K+KB, N, HB                         ! 8
                    KKB = MIN( N-KK+1, HB )
                    IKK2 = IJ2K( UPLO, K, KK, N, HB )
                    IK2K2 = IJ2K( UPLO, KK, KK, N, HB )
                    IF (  KK .GT. K+KB ) THEN
                       DO KKK = K+KB, KK-KKB, HB             ! 2
                          KKKB = MIN( N-KKK+1, HB )
                          IKK3 = IJ2K( UPLO, K, KKK, N, HB )
                          IK3K2 = IJ2K( UPLO, KKK, KK, N, HB )
                          CALL DGEMM( 'Transpose', 'No transpose',HB, KKB, HB, -ONE, ABP( IKK3 ), KKK+KKKB-1,BBP( IKK2 ), KK+KKB-1,ONE, ABP( IK3K2 ), KK+KKB-1)
                          CALL DGEMM( 'Transpose', 'No transpose',HB, KKB, HB, -ONE, BBP( IKK3 ), KKK+KKKB-1,ABP( IKK2 ), KK+KKB-1,ONE, ABP( IK3K2 ), KK+KKB-1)
                       ENDDO                                 ! 2                         CONTINUE
                    ENDIF
                    CALL DSYR2K( UPLO, 'Transpose', KKB, KB, -ONE,ABP( IKK2 ), KK+KKB-1, BBP( IKK2 ),KK+KKB-1,ONE, ABP( IK2K2 ), KK+KKB-1 )
                 ENDDO  !  8                   CONTINUE
                 !123456789012345678901234567890123456789012345678901234567890123456789012
                 DO KK = K+KB, N, HB   ! 5
                    KKB = MIN( N-KK+1, HB )
                    IKK2 = IJ2K( UPLO, K, KK, N, HB )
                    IK2K2 = IJ2K( UPLO, KK, KK, N, HB )
                    CALL DSYMM( 'Left', UPLO, KB, KKB, -HALF,ABP( IKK ), K+KB-1, BBP( IKK2 ),KK+KKB-1,ONE, ABP( IKK2 ), KK+KKB-1 )
                    !
                    !                       Reduce triangular system
                    !
                    IF (  KK .GT. K+KB ) THEN
                       DO KKK = K+KB, KK-KKB, HB     ! 3
                          KKKB = MIN( N-KKK+1, HB )
                          IKK3 = IJ2K( UPLO, K, KKK, N, HB )
                          IK3K2 = IJ2K( UPLO, KKK, KK, N, HB )
                          CALL DGEMM( 'No trans', 'No trans',HB, KKB, HB,-ONE, ABP( IKK3 ), KKK+KKKB-1,BBP( IK3K2 ), KK+KKB-1,ONE, ABP( IKK2 ), KK+KKB-1)
                       ENDDO ! 3                         CONTINUE
                    END IF
                    !
                    !                       compute block-column of reduced system
                    !
                    CALL DTRSM( 'Right', UPLO, 'No transpose','Non-unit', KB, KKB, ONE,BBP( IK2K2 ), KK+KKB-1, ABP( IKK2 ),KK+KKB-1 )
                    !123456789012345678901234567890123456789012345678901234567890123456789012
                 ENDDO ! 5                   CONTINUE
              END IF
           ENDDO ! 10            CONTINUE
        END IF
     END IF
  END IF
  RETURN
  
20 CONTINUE
  INFO = INFO + K - 1
  RETURN
  !
  !     End of DSHGS4
  !
END SUBROUTINE DSHGS4

SUBROUTINE DSHGST( ITYPE, UPLO, N, AP, BP, HB, INFO )
  IMPLICIT NONE
  !
  !
  !
  CHARACTER      UPLO
  INTEGER        INFO, N, ITYPE, HB
  DOUBLE         PRECISION AP( * ), BP( * )
  !
  !  Purpose
  !  =======
  !
  !  DSHGST reduces a real symmetric-definite generalized eigenproblem
  !  to standard form.
  !
  !  If ITYPE = 1, the problem is A*x = lambda*B*x,
  !  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
  !
  !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
  !  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
  !
  !  B must have been previously factorized as U**T*U or L*L**T by DPHTRF.
  !
  !  Arguments
  !  =========
  !
  !  ITYPE   (input) INTEGER
  !          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
  !          = 2 or 3: compute U*A*U**T or L**T*A*L, not implemented.
  !
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix AP is stored as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of AP
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of AP
  !                                  is to be referenced, not implemented.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - REAL*8
  !           On entry, AP stores the matrix in packed storage format.
  !           In order to store the (not referenced elements in the
  !           block diagonal), AP must be at least of size
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !           On exit, AP stores the matrix in block packed storage.
  !
  !  BP     - REAL*8
  !           On entry, BP stores the triangular factor from the 
  !           Cholesky factorization of B, as returned by DPPTRFHB.
  !           In order to store the (not referenced elements in the
  !           block diagonal), BP must be at least of size
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !  INFO  - (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order K is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !  =====================================================================
  !           
  Integer            j, ilaenv 
  REAL*8             CP(4), FLOP
  logical            LSAME, UPPER
  
  external ilaenv, LSAME
  
  !
  !     Test the input parameters.
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
     INFO = -1
  ELSE IF( .NOT.UPPER ) THEN
     INFO = -2
  ELSE IF( N.LT.0 ) THEN
     INFO = -3
  ELSE IF( HB.LT.0 ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSHGST', -INFO )
     RETURN
  END IF
  
  IF(HB.le.1) then
     call DSPGST( ITYPE, UPLO, N, AP, BP, INFO ) 
  ELSE
     
     call CPUTIM(CP(1))
     !
     !      call the subroutines
     !			
     !
     !      convert from LAPACK upper triangular packed storage format
     !                     into a upper triangular block packed matrix
     !
     call dp2bp('U',N,AP,HB)
     call dp2bp('U',N,BP,HB)
     call CPUTIM(CP(2))
     !
     !      reduce the real symmetrix-definite generalized eigenproblem
     !      to standard form
     !
     call DSHGS4(ITYPE, UPLO, N, AP, BP, HB, INFO)
     call CPUTIM(CP(3))
     IF (INFO .NE. 0) THEN
        CALL OUTERR('DSHGST','DSHGS4 aborted unsuccessfully.')
        GOTO 999
     ENDIF
     !
     !      convert back to upper triangular packed storage format
     !
     call dbp2p('U',N,AP,HB)
     call dbp2p('U',N,BP,HB)
     call CPUTIM(CP(4))
     
     !
     !      timing output
     !
     !         DO 30 j = 1, 3
     !           CP(j) = CP(j+1)-CP(j)
     ! 30      CONTINUE
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Packing' , CP(1) , 1E-6*Flop/CP(1)
     !         Flop = 2*N/3
     !         WRITE (*,1001) 'Transformation' , CP(2) , 1E-6*Flop/CP(2)
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Unpacking' , CP(3) , 1E-6*Flop/CP(3)
     ! 1001    FORMAT (1X,'DPPTRFB(',A,') :',t30,f7.3:t40,f8.3,' Mflops')
     !	 write (*,*) 'INFO = ', INFO
     
  END IF
  RETURN
  !	
  !	Result stored using LAPACK packed format
  !
  !
999 STOP 'LAPW1 - Error'
  !
END SUBROUTINE DSHGST

SUBROUTINE DSHTR4( UPLO, N, ABP, D, E, TAU, WORK, LWORK, HB, INFO )
  IMPLICIT NONE
  !
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER          UPLO
  INTEGER            INFO, LWORK, N, HB
  !     ..
  !     .. Array Arguments ..
  DOUBLE PRECISION   ABP( * ), D( * ), E( * ), TAU( * ), WORK( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DSHTR4 reduces a real symmetric matrix A to real symmetric
  !  tridiagonal form T by an orthogonal similarity transformation:
  !  Q**T * A * Q = T.
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER*1
  !          = 'U':  Upper triangle of A is stored;
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  ABP     (input/output) DOUBLE PRECISION array, 
  !          dimension (n*(n+1)/2+(n/hb+1)*hb*(hb-1)/2)
  !          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !          N-by-N upper triangular part of A contains the upper
  !          triangular part of the matrix A, and the strictly lower
  !          triangular part of A is not referenced.  
  !          On exit, if UPLO = 'U', the diagonal and first superdiagonal
  !          of A are overwritten by the corresponding elements of the
  !          tridiagonal matrix T, and the elements above the first
  !          superdiagonal, with the array TAU, represent the orthogonal
  !          matrix Q as a product of elementary reflectors.
  !          See Further Details.
  !
  !  D       (output) DOUBLE PRECISION array, dimension (N)
  !          The diagonal elements of the tridiagonal matrix T:
  !          D(i) = A(i,i).
  !
  !  E       (output) DOUBLE PRECISION array, dimension (N-1)
  !          The off-diagonal elements of the tridiagonal matrix T:
  !          E(i) = A(i,i+1) if UPLO = 'U'.
  !
  !  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
  !          The scalar factors of the elementary reflectors (see Further
  !          Details).
  !
  !  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
  !
  !  LWORK   (input) INTEGER
  !          LWORK = N*HB, where HB is the blocksize.
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !
  !  Further Details
  !  ===============
  !
  !  If UPLO = 'U', the matrix Q is represented as a product of elementary
  !  reflectors
  !
  !     Q = H(n-1) . . . H(2) H(1).
  !
  !  Each H(i) has the form
  !
  !     H(i) = I - tau * v * v'
  !
  !  where tau is a real scalar, and v is a real vector with
  !  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
  !  A(1:i-1,i+1), and tau in TAU(i).
  !
  !  The contents of A on exit are illustrated by the following examples
  !  with n = 5:
  !
  !  if UPLO = 'U':               
  !
  !    (  d   e   v2  v3  v4 )     
  !    (      d   e   v3  v4 )    
  !    (          d   e   v4 )   
  !    (              d   e  )  
  !    (                  d  ) 
  !
  !  where d and e denote diagonal and off-diagonal elements of T, and vi
  !  denotes an element of the vector defining H(i).
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  DOUBLE PRECISION   ONE
  PARAMETER          ( ONE = 1.0D0 )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            I, IB, J, LDWORK, NB
  INTEGER            I1I, I1J, IJI, IJJ, IJ1J, LDI, LDJ
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DLATD4, DSYR2K, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            ILAENV
  real               dtime2
  integer            IJ2K
  EXTERNAL           LSAME, ILAENV, IJ2K
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LWORK.LT.(N*HB) ) THEN
     INFO = -8
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSHTR4', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) THEN
     WORK( 1 ) = 1
     RETURN
  END IF
  !
  !     Determine the block size.
  !
  !
  NB=MIN(HB,N)
  LDWORK = N
  !
  IF( UPPER ) THEN
     !
     !
     !         T1=0
     !	 T2=0
     DO I= NB*((N-1)/NB)+1, 1, -NB   ! 20
        IB = MIN(NB, N-I+1)
        LDI = I+IB-1
        !           write(*,*) 'I=',I,' IB=',IB,' N=',N,' NB=',NB,
        !     $                ' LDWORK',LDWORK
        !
        !           Reduce columns i:i+ib-1 to tridiagonal form and form the
        !           matrix W which is needed to update the unreduced part of
        !           the matrix
        !
        CALL DLATD4( UPLO, LDI, IB, ABP, E, TAU, WORK, LDWORK, NB )
        !	    T1 = T1 + dtime2(lasttime)
        !
        !           Update the unreduced submatrix A(1:i-1,1:i-1), using an
        !           update of the form:  A := A - V*W' - W*V'
        !
        I1I = IJ2K(UPLO, 1, I, N, NB)
        DO J = 1, I-1, NB     ! 5
           LDJ = J+NB-1
           !              write(*,*) 'DSYR2K I=',I,' J=',J
           I1J = IJ2K(UPLO, 1, J, N, NB)
           IJI = IJ2K(UPLO, J, I, N, NB)
           IJJ = IJ2K(UPLO, J, J, N, NB)
           IF( J.NE.1 ) THEN
              CALL DGEMM( 'NoTranspose', 'Transpose', J-1, NB, IB,-ONE, ABP(I1I), LDI, WORK(J), LDWORK, ONE, ABP(I1J), LDJ)
              CALL DGEMM( 'NoTranspose', 'Transpose', J-1, NB, IB,-ONE, WORK, LDWORK, ABP(IJI), LDI, ONE, ABP(I1J), LDJ)
           END IF
           CALL DSYR2K( UPLO, 'NoTranspose', NB, IB,-ONE, ABP(IJI), LDI, WORK(J), LDWORK,ONE, ABP(IJJ), LDJ)
        ENDDO ! 5          CONTINUE
        
        !	    T2 = T2 + dtime2(lasttime)
        !
        !           Copy superdiagonal elements back into A, and diagonal
        !           elements into D
        !
        DO J = I, I + IB - 1    ! 10
           !	       write(*,*) 'Copy: I=',I,' J=',J
           IJJ = IJ2K(UPLO, J, J, N, NB)
           IF( J.GT.1 ) THEN
              IJ1J = IJ2K(UPLO, J-1, J, N, NB)
              ABP(IJ1J) = E( J-1 )
           END IF
           D( J ) = ABP( IJJ )
        ENDDO !10       CONTINUE
     ENDDO ! 20      CONTINUE
  END IF
  !
  !      write(*,*) '           DLATD4: ', T1
  !      write(*,*) '           DSYR2K: ', T2
  RETURN
  !
  !     End of DSHTR4
  !
END SUBROUTINE DSHTR4

SUBROUTINE DSHTRD( UPLO, N, ABP, D, E, TAU, WORK, LWORK, HB, INFO )
  IMPLICIT NONE
  !
  !
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER          UPLO
  INTEGER            INFO, LWORK, N, HB
  !     ..
  !     .. Array Arguments ..
  DOUBLE PRECISION   ABP( * ), D( * ), E( * ), TAU( * ), WORK( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DSHTRD reduces a real symmetric matrix A to real symmetric
  !  tridiagonal form T by an orthogonal similarity transformation:
  !  Q**T * A * Q = T.
  !
  !  The matrix A is stored in packed mode and will be restored in packed 
  !  mode on exit. 
  !  Computation is done using block packed storage mode.
  !
  !  Note: this routine has three parameters which do not exist in the 
  !        call to DSPTRD:
  !        WORK, LWORK, HB
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER*1
  !          = 'U':  Upper triangle of A is stored;
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  ABP     (input/output) DOUBLE PRECISION array, 
  !          dimension (n*(n+1)/2+(n/hb+1)*hb*(hb-1)/2)
  !          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !          N-by-N upper triangular part of A contains the upper
  !          triangular part of the matrix A, and the strictly lower
  !          triangular part of A is not referenced.  
  !          On exit, if UPLO = 'U', the diagonal and first superdiagonal
  !          of A are overwritten by the corresponding elements of the
  !          tridiagonal matrix T, and the elements above the first
  !          superdiagonal, with the array TAU, represent the orthogonal
  !          matrix Q as a product of elementary reflectors.
  !          See Further Details.
  !
  !  D       (output) DOUBLE PRECISION array, dimension (N)
  !          The diagonal elements of the tridiagonal matrix T:
  !          D(i) = A(i,i).
  !
  !  E       (output) DOUBLE PRECISION array, dimension (N-1)
  !          The off-diagonal elements of the tridiagonal matrix T:
  !          E(i) = A(i,i+1) if UPLO = 'U'.
  !
  !  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
  !          The scalar factors of the elementary reflectors (see Further
  !          Details).
  !
  !  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
  !
  !  LWORK   (input) INTEGER
  !          LWORK >= N*HB.
  !
  !  HB      (input) INTEGER
  !          The blocksize used in packed storage scheme.
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !
  !  Further Details
  !  ===============
  !
  !  If UPLO = 'U', the matrix Q is represented as a product of elementary
  !  reflectors
  !
  !     Q = H(n-1) . . . H(2) H(1).
  !
  !  Each H(i) has the form
  !
  !     H(i) = I - tau * v * v'
  !
  !  where tau is a real scalar, and v is a real vector with
  !  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
  !  A(1:i-1,i+1), and tau in TAU(i).
  !
  !  The contents of A on exit are illustrated by the following examples
  !  with n = 5:
  !
  !  if UPLO = 'U':               
  !
  !    (  d   e   v2  v3  v4 )     
  !    (      d   e   v3  v4 )    
  !    (          d   e   v4 )   
  !    (              d   e  )  
  !    (                  d  ) 
  !
  !  where d and e denote diagonal and off-diagonal elements of T, and vi
  !  denotes an element of the vector defining H(i).
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  DOUBLE PRECISION   ONE
  PARAMETER          ( ONE = 1.0D0 )
  INTEGER            MINHB
  PARAMETER          ( MINHB = 10 )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            NB
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DLATD4, DSYR2K, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            ILAENV
  real               dtime2
  integer            IJ2K
  EXTERNAL           LSAME, ILAENV, IJ2K
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LWORK.LT.(N*HB) ) THEN
     INFO = -8
  ELSE IF( HB.LT.1 ) THEN
     INFO = -9
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSHTRD', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) THEN
     WORK( 1 ) = 1
     RETURN
  END IF
  !
  !     Hetermine the block size.
  !
  NB=MIN(HB,N)
  IF( NB.LE.MINHB ) THEN
     CALL DSPTRD( UPLO, N, ABP, D, E, TAU, INFO )
  ELSE IF( UPPER ) THEN
     !
     !      convert from LAPACK upper triangular packed storage format
     !                     into a upper triangular block packed matrix
     !
     !	call CPUTIM(CP(1))
     call dp2bp('U',N,ABP,NB)
     !	call CPUTIM(CP(2))
     !      
     !      compute the tridiagonalization of A
     !
     call DSHTR4( UPLO, N, ABP, D, E, TAU, WORK, LWORK, NB, INFO )
     !      call CPUTIM(CP(3))
     !
     !      convert back to upper triangular packed storage format
     !
     call dbp2p('U',N,ABP,NB)
     !      call CPUTIM(CP(4))
     
     !
     !      timing output
     !
     !         DO 30 j = 1, 3
     !           CP(j) = CP(j+1)-CP(j)
     ! 30      CONTINUE
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Packing' , CP(1) , 1E-6*Flop/CP(1)
     !         Flop = (N*DBLE(N+1)*(2*N+1))/3
     !         WRITE (*,1001) 'Tridiagonalization' , CP(2) , 1E-6*Flop/CP(2)
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Unpacking' , CP(3) , 1E-6*Flop/CP(3)
     ! 1001    FORMAT (1X,'DPHTRF(',A,') :',t30,f7.3:t40,f8.3,' Mflops')
     !	 write (*,*) 'INFO = ', INFO
     !
  END IF
  RETURN
  !	
  !	Result stored using LAPACK packed format
  !
  !
  !
  !     End of DSHTRD
  !
END SUBROUTINE DSHTRD

SUBROUTINE DSTEZ2( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,NUME, NWL, INFO )
  IMPLICIT NONE
  !
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER          ORDER, RANGE
  INTEGER            IL, INFO, IU, M, N, NSPLIT, NUME, NWL
  DOUBLE PRECISION   ABSTOL, VL, VU
  !     ..
  !     .. Array Arguments ..
  INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
  DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DSTEZ2 computes the eigenvalues of a symmetric tridiagonal
  !  matrix T.  The user may ask for all eigenvalues, all eigenvalues
  !  in the half-open interval (VL, VU], or the IL-th through IU-th
  !  eigenvalues.
  !
  !  To avoid overflow, the matrix must be scaled so that its
  !  largest element is no greater than overflow**(1/2) *
  !  underflow**(1/4) in absolute value, and for greatest
  !  accuracy, it should not be much smaller than that.
  !
  !  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
  !  Matrix", Report CS41, Computer Science Dept., Stanford
  !  University, July 21, 1966.
  !
  !  Arguments
  !  =========
  !
  !  RANGE   (input) CHARACTER
  !          = 'A': ("All")   all eigenvalues will be found.
  !          = 'V': ("Value") all eigenvalues in the half-open interval
  !                           (VL, VU] will be found.
  !          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
  !                           entire matrix) will be found.
  !          = 'C': ("Combined") the eigenvalues in the half-open interval
  !                           (VL, VU] will be found, but if there are more
  !                           than M Eigenvektors, only the lowest M will
  !                           be found.
  !
  !  ORDER   (input) CHARACTER
  !          = 'B': ("By Block") the eigenvalues will be grouped by
  !                              split-off block (see IBLOCK, ISPLIT) and
  !                              ordered from smallest to largest within
  !                              the block.
  !          = 'E': ("Entire matrix")
  !                              the eigenvalues for the entire matrix
  !                              will be ordered from smallest to
  !                              largest.
  !
  !  N       (input) INTEGER
  !          The order of the tridiagonal matrix T.  N >= 0.
  !
  !  VL      (input) DOUBLE PRECISION
  !  VU      (input) DOUBLE PRECISION
  !          If RANGE='V', the lower and upper bounds of the interval to
  !          be searched for eigenvalues.  Eigenvalues less than or equal
  !          to VL, or greater than VU, will not be returned.  VL < VU.
  !          Not referenced if RANGE = 'A' or 'I'.
  !
  !  IL      (input) INTEGER
  !  IU      (input) INTEGER
  !          If RANGE='I', the indices (in ascending order) of the
  !          smallest and largest eigenvalues to be returned.
  !          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
  !          Not referenced if RANGE = 'A' or 'V'.
  !
  !  ABSTOL  (input) DOUBLE PRECISION
  !          The absolute tolerance for the eigenvalues.  An eigenvalue
  !          (or cluster) is considered to be located if it has been
  !          determined to lie in an interval whose width is ABSTOL or
  !          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
  !          will be used, where |T| means the 1-norm of T.
  !
  !          Eigenvalues will be computed most accurately when ABSTOL is
  !          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
  !
  !  D       (input) DOUBLE PRECISION array, dimension (N)
  !          The n diagonal elements of the tridiagonal matrix T.
  !
  !  E       (input) DOUBLE PRECISION array, dimension (N-1)
  !          The (n-1) off-diagonal elements of the tridiagonal matrix T.
  !
  !  M       (output) INTEGER
  !          The actual number of eigenvalues found. 0 <= M <= N.
  !          (See also the description of INFO=2,3.)
  !
  !  NSPLIT  (output) INTEGER
  !          The number of diagonal blocks in the matrix T.
  !          1 <= NSPLIT <= N.
  !
  !  W       (output) DOUBLE PRECISION array, dimension (N)
  !          On exit, the first M elements of W will contain the
  !          eigenvalues.  (DSTEZ2 may use the remaining N-M elements as
  !          workspace.)
  !
  !  IBLOCK  (output) INTEGER array, dimension (N)
  !          At each row/column j where E(j) is zero or small, the
  !          matrix T is considered to split into a block diagonal
  !          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
  !          block (from 1 to the number of blocks) the eigenvalue W(i)
  !          belongs.  (DSTEZ2 may use the remaining N-M elements as
  !          workspace.)
  !
  !  ISPLIT  (output) INTEGER array, dimension (N)
  !          The splitting points, at which T breaks up into submatrices.
  !          The first submatrix consists of rows/columns 1 to ISPLIT(1),
  !          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
  !          etc., and the NSPLIT-th consists of rows/columns
  !          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
  !          (Only the first NSPLIT elements will actually be used, but
  !          since the user cannot know a priori what value NSPLIT will
  !          have, N words must be reserved for ISPLIT.)
  !
  !  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
  !
  !  IWORK   (workspace) INTEGER array, dimension (3*N)
  !
  !  NUME    (input) INTEGER
  !          The maximum number of Eigenvalues allowed
  !
  !  NWL     (output) INTEGER
  !          The number of Eigenvalues less or equal VL
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !          > 0:  some or all of the eigenvalues failed to converge or
  !                were not computed:
  !                =1 or 3: Bisection failed to converge for some
  !                        eigenvalues; these eigenvalues are flagged by a
  !                        negative block number.  The effect is that the
  !                        eigenvalues may not be as accurate as the
  !                        absolute and relative tolerances.  This is
  !                        generally caused by unexpectedly inaccurate
  !                        arithmetic.
  !                =2 or 3: RANGE='I' only: Not all of the eigenvalues
  !                        IL:IU were found.
  !                        Effect: M < IU+1-IL
  !                        Cause:  non-monotonic arithmetic, causing the
  !                                Sturm sequence to be non-monotonic.
  !                        Cure:   recalculate, using RANGE='A', and pick
  !                                out eigenvalues IL:IU.  In some cases,
  !                                increasing the PARAMETER "FUDGE" may
  !                                make things work.
  !                = 4:    RANGE='I', and the Gershgorin interval
  !                        initially used was too small.  No eigenvalues
  !                        were computed.
  !                        Probable cause: your machine has sloppy
  !                                        floating-point arithmetic.
  !                        Cure: Increase the PARAMETER "FUDGE",
  !                              recompile, and try again.
  !
  !  Internal Parameters
  !  ===================
  !
  !  RELFAC  DOUBLE PRECISION, default = 2.0e0
  !          The relative tolerance.  An interval (a,b] lies within
  !          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
  !          where "ulp" is the machine precision (distance from 1 to
  !          the next larger floating point number.)
  !
  !  FUDGE   DOUBLE PRECISION, default = 2
  !          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
  !          a value of 1 should work, but on machines with sloppy
  !          arithmetic, this needs to be larger.  The default for
  !          publicly released versions should be large enough to handle
  !          the worst machine around.  Note that this has no effect
  !          on accuracy of the solution.
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  DOUBLE PRECISION   ZERO, ONE, TWO, HALF
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 1.0D0 / TWO )
  DOUBLE PRECISION   FUDGE, RELFAC
  PARAMETER          ( FUDGE = 2.0D0, RELFAC = 2.0D0 )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            NCNVRG, TOOFEW
  INTEGER            IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO, IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX, ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWU
  DOUBLE PRECISION   ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN, TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL
  !     ..
  !     .. Local Arrays ..
  INTEGER            IDUMMA( 1 )
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            ILAENV
  DOUBLE PRECISION   DLAMCH
  EXTERNAL           LSAME, ILAENV, DLAMCH
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DLAEBZ, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
  !     ..
  !     .. Executable Statements ..
  !
  INFO = 0
  !
  !     Decode RANGE
  !
  IF( LSAME( RANGE, 'A' ) ) THEN
     IRANGE = 1
  ELSE IF( LSAME( RANGE, 'V' ) ) THEN
     IRANGE = 2
  ELSE IF( LSAME( RANGE, 'I' ) ) THEN
     IRANGE = 3
  ELSE IF( LSAME( RANGE, 'C' ) ) THEN
     IRANGE = 4
  ELSE
     IRANGE = 0
  END IF
  !
  !     Decode ORDER
  !
  IF( LSAME( ORDER, 'B' ) ) THEN
     IORDER = 2
  ELSE IF( LSAME( ORDER, 'E' ) ) THEN
     IORDER = 1
  ELSE
     IORDER = 0
  END IF
  !
  !     Check for Errors
  !
  IF( IRANGE.LE.0 ) THEN
     INFO = -1
  ELSE IF( IORDER.LE.0 ) THEN
     INFO = -2
  ELSE IF( N.LT.0 ) THEN
     INFO = -3
  ELSE IF( ( ( IRANGE.EQ.2) .OR. (IRANGE.EQ.4) ) .AND. VL.GE.VU ) THEN
     INFO = -5
  ELSE IF( IRANGE.EQ.3 .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) ) THEN
     INFO = -6
  ELSE IF( IRANGE.EQ.3 .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) THEN
     INFO = -7
  END IF
  !
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSTEZ2', -INFO )
     RETURN
  END IF
  !
  !     Initialize error flags
  !
  INFO = 0
  NCNVRG = .FALSE.
  TOOFEW = .FALSE.
  
  !
  !     Quick return if possible
  !
  M = 0
  IF( N.EQ.0 ) &
       RETURN
  !
  !     Simplifications:
  !
  IF( IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.N ) &
       IRANGE = 1
  !
  !     Get machine constants
  !     NB is the minimum vector length for vector bisection, or 0
  !     if only scalar is to be done.
  !
  SAFEMN = DLAMCH( 'S' )
  ULP = DLAMCH( 'P' )
  RTOLI = ULP*RELFAC
  NB = ILAENV( 1, 'DSTEZ2', ' ', N, -1, -1, -1 )
  IF( NB.LE.1 ) NB = 0
  !
  !     Special Case when N=1
  !
  IF( N.EQ.1 ) THEN
     NSPLIT = 1
     ISPLIT( 1 ) = 1
     IF( ( ( IRANGE.EQ.2 ) .OR. ( IRANGE.EQ.4 ) ) .AND.( VL.GE.D( 1 ) .OR. VU.LT.D( 1 ) ) ) THEN
        M = 0
     ELSE
        W( 1 ) = D( 1 )
        IBLOCK( 1 ) = 1
        M = 1
     END IF
     RETURN
  END IF
  !
  !     Compute Splitting Points
  !
  NSPLIT = 1
  WORK( N ) = ZERO
  PIVMIN = ONE
  !
  DO J = 2, N    ! 10
     TMP1 = E( J-1 )**2
     IF( ABS( D( J )*D( J-1 ) )*ULP**2+SAFEMN.GT.TMP1 ) THEN
        ISPLIT( NSPLIT ) = J - 1
        NSPLIT = NSPLIT + 1
        WORK( J-1 ) = ZERO
     ELSE
        WORK( J-1 ) = TMP1
        PIVMIN = MAX( PIVMIN, TMP1 )
     END IF
  ENDDO ! 10   CONTINUE
  ISPLIT( NSPLIT ) = N
  PIVMIN = PIVMIN*SAFEMN
  !
  !     Compute Interval and ATOLI
  !
  IF( IRANGE.EQ.3 ) THEN
     !
     !        RANGE='I': Compute the interval containing eigenvalues
     !                   IL through IU.
     !
     !        Compute Gershgorin interval for entire (split) matrix
     !        and use it as the initial interval
     !
     GU = D( 1 )
     GL = D( 1 )
     TMP1 = ZERO
     !
     DO J = 1, N - 1    ! 20
        TMP2 = SQRT( WORK( J ) )
        GU = MAX( GU, D( J )+TMP1+TMP2 )
        GL = MIN( GL, D( J )-TMP1-TMP2 )
        TMP1 = TMP2
     ENDDO ! 20      CONTINUE
     !
     GU = MAX( GU, D( N )+TMP1 )
     GL = MIN( GL, D( N )-TMP1 )
     TNORM = MAX( ABS( GL ), ABS( GU ) )
     GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
     GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
     !
     !        Compute Iteration parameters
     !
     ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2
     IF( ABSTOL.LE.ZERO ) THEN
        ATOLI = ULP*TNORM
     ELSE
        ATOLI = ABSTOL
     END IF
     !
     WORK( N+1 ) = GL
     WORK( N+2 ) = GL
     WORK( N+3 ) = GU
     WORK( N+4 ) = GU
     WORK( N+5 ) = GL
     WORK( N+6 ) = GU
     IWORK( 1 ) = -1
     IWORK( 2 ) = -1
     IWORK( 3 ) = N + 1
     IWORK( 4 ) = N + 1
     IWORK( 5 ) = IL - 1
     IWORK( 6 ) = IU
     !
     CALL DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E,WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT,IWORK, W, IBLOCK, IINFO )
     !
     IF( IWORK( 6 ).EQ.IU ) THEN
        WL = WORK( N+1 )
        WLU = WORK( N+3 )
        NWL = IWORK( 1 )
        WU = WORK( N+4 )
        WUL = WORK( N+2 )
        NWU = IWORK( 4 )
     ELSE
        WL = WORK( N+2 )
        WLU = WORK( N+4 )
        NWL = IWORK( 2 )
        WU = WORK( N+3 )
        WUL = WORK( N+1 )
        NWU = IWORK( 3 )
     END IF
     !
     IF( NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N ) THEN
        INFO = 4
        RETURN
     END IF
  ELSE
     !
     !        RANGE='A' or 'V' -- Set ATOLI
     !
     TNORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),ABS( D( N ) )+ABS( E( N-1 ) ) )
     !
     DO J = 2, N - 1    ! 30
        TNORM = MAX( TNORM, ABS( D( J ) )+ABS( E( J-1 ) )+ ABS( E( J ) ) )
     ENDDO ! 30      CONTINUE
     !
     IF( ABSTOL.LE.ZERO ) THEN
        ATOLI = ULP*TNORM
     ELSE
        ATOLI = ABSTOL
     END IF
     !
     IF( IRANGE.EQ.2 .OR. IRANGE.EQ.4 ) THEN
        WL = VL
        WU = VU
     END IF
  END IF
  !
  !     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
  !     NWL accumulates the number of eigenvalues .le. WL,
  !     NWU accumulates the number of eigenvalues .le. WU
  !
  M = 0
  IEND = 0
  INFO = 0
  NWL = 0
  NWU = 0
  !
  DO JB = 1, NSPLIT    ! 70
     IOFF = IEND
     IBEGIN = IOFF + 1
     IEND = ISPLIT( JB )
     IN = IEND - IOFF
     !
     IF( IN.EQ.1 ) THEN
        !
        !           Special Case -- IN=1
        !
        IF( IRANGE.EQ.1 .OR. WL.GE.D( IBEGIN )-PIVMIN ) NWL = NWL + 1
        IF( IRANGE.EQ.1 .OR. WU.GE.D( IBEGIN )-PIVMIN ) NWU = NWU + 1
        IF( IRANGE.EQ.1 .OR. ( WL.LT.D( IBEGIN )-PIVMIN .AND. WU.GE.D( IBEGIN )-PIVMIN ) ) THEN
           M = M + 1
           W( M ) = D( IBEGIN )
           IBLOCK( M ) = JB
        END IF
     ELSE
        !
        !           General Case -- IN > 1
        !
        !           Compute Gershgorin Interval
        !           and use it as the initial interval
        !
        GU = D( IBEGIN )
        GL = D( IBEGIN )
        TMP1 = ZERO
        !
        DO J = IBEGIN, IEND - 1    ! 40
           TMP2 = ABS( E( J ) )
           GU = MAX( GU, D( J )+TMP1+TMP2 )
           GL = MIN( GL, D( J )-TMP1-TMP2 )
           TMP1 = TMP2
        ENDDO !40         CONTINUE
        !
        GU = MAX( GU, D( IEND )+TMP1 )
        GL = MIN( GL, D( IEND )-TMP1 )
        BNORM = MAX( ABS( GL ), ABS( GU ) )
        GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
        GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN
        !
        !           Compute ATOLI for the current submatrix
        !
        IF( ABSTOL.LE.ZERO ) THEN
           ATOLI = ULP*MAX( ABS( GL ), ABS( GU ) )
        ELSE
           ATOLI = ABSTOL
        END IF
        !
        IF( IRANGE.GT.1 ) THEN
           IF( GU.LT.WL ) THEN
              NWL = NWL + IN
              NWU = NWU + IN
              CYCLE  !GO TO 70
           END IF
           GL = MAX( GL, WL )
           GU = MIN( GU, WU )
           IF( GL.GE.GU )  CYCLE  !GO TO 70
        END IF
        !
        !           Set Up Initial Interval
        !
        WORK( N+1 ) = GL
        WORK( N+IN+1 ) = GU
        CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM, IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
        !
        NWL = NWL + IWORK( 1 )
        NWU = NWU + IWORK( IN+1 )
        IWOFF = M - IWORK( 1 )
        !
        !           Compute Eigenvalues
        !
        ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2
        CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ), IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT, IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
        !
        !           Copy Eigenvalues Into W and IBLOCK
        !           Use -JB for block number for unconverged eigenvalues.
        !
        DO J = 1, IOUT     ! 60
           TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) )
           !
           !              Flag non-convergence.
           !
           IF( J.GT.IOUT-IINFO ) THEN
              NCNVRG = .TRUE.
              IB = -JB
           ELSE
              IB = JB
           END IF
           DO JE = IWORK( J ) + 1 + IWOFF, IWORK( J+IN ) + IWOFF  ! 50
              W( JE ) = TMP1
              IBLOCK( JE ) = IB
           ENDDO !50            CONTINUE
        ENDDO !   60       CONTINUE
        !
        M = M + IM
     END IF
  ENDDO !70   CONTINUE
  !
  !     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
  !     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
  !
  IF( IRANGE.EQ.3 ) THEN
     IM = 0
     IDISCL = IL - 1 - NWL
     IDISCU = NWU - IU
     !
     IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
        DO JE = 1, M    ! 80
           IF( W( JE ).LE.WLU .AND. IDISCL.GT.0 ) THEN
              IDISCL = IDISCL - 1
           ELSE IF( W( JE ).GE.WUL .AND. IDISCU.GT.0 ) THEN
              IDISCU = IDISCU - 1
           ELSE
              IM = IM + 1
              W( IM ) = W( JE )
              IBLOCK( IM ) = IBLOCK( JE )
           END IF
        ENDDO ! 80         CONTINUE
        M = IM
     END IF
     IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
        !
        !           Code to deal with effects of bad arithmetic:
        !           Some low eigenvalues to be discarded are not in (WL,WLU],
        !           or high eigenvalues to be discarded are not in (WUL,WU]
        !           so just kill off the smallest IDISCL/largest IDISCU
        !           eigenvalues, by simply finding the smallest/largest
        !           eigenvalue(s).
        !
        !           (If N(w) is monotone non-decreasing, this should never
        !               happen.)
        !
        IF( IDISCL.GT.0 ) THEN
           WKILL = WU
           DO JDISC = 1, IDISCL   ! 100
              IW = 0
              DO JE = 1, M  ! 90
                 IF( IBLOCK( JE ).NE.0 .AND.( W( JE ).LT.WKILL .OR. IW.EQ.0 ) ) THEN
                    IW = JE
                    WKILL = W( JE )
                 END IF
              ENDDO ! 90               CONTINUE
              IBLOCK( IW ) = 0
           ENDDO !  100          CONTINUE
        END IF
        IF( IDISCU.GT.0 ) THEN
           !
           WKILL = WL
           DO JDISC = 1, IDISCU    ! 120
              IW = 0
              DO JE = 1, M          ! 110
                 IF( IBLOCK( JE ).NE.0 .AND.( W( JE ).GT.WKILL .OR. IW.EQ.0 ) ) THEN
                    IW = JE
                    WKILL = W( JE )
                 END IF
              ENDDO ! 110              CONTINUE
              IBLOCK( IW ) = 0
           ENDDO !  120          CONTINUE
        END IF
        IM = 0
        DO JE = 1, M   ! 130
           IF( IBLOCK( JE ).NE.0 ) THEN
              IM = IM + 1
              W( IM ) = W( JE )
              IBLOCK( IM ) = IBLOCK( JE )
           END IF
        ENDDO ! 130        CONTINUE
        M = IM
     END IF
     IF( IDISCL.LT.0 .OR. IDISCU.LT.0 ) THEN
        TOOFEW = .TRUE.
     END IF
  END IF
  !
  !     If ORDER='B', do nothing -- the eigenvalues are already sorted
  !        by block.
  !     If ORDER='E', sort the eigenvalues from smallest to largest
  !
  IF( ( IORDER.EQ.1 .AND. NSPLIT.GT.1 ) .OR. ( IORDER.EQ.4 .AND. M.GT.NUME ) ) THEN
     DO JE = 1, M - 1   ! 150
        IE = 0
        TMP1 = W( JE )
        DO J = JE + 1, M    ! 140
           IF( W( J ).LT.TMP1 ) THEN
              IE = J
              TMP1 = W( J )
           END IF
        ENDDO ! 140        CONTINUE
        !
        IF( IE.NE.0 ) THEN
           ITMP1 = IBLOCK( IE )
           W( IE ) = W( JE )
           IBLOCK( IE ) = IBLOCK( JE )
           W( JE ) = TMP1
           IBLOCK( JE ) = ITMP1
        END IF
     ENDDO !  150    CONTINUE
  END IF
  !
  !     Do not return more Eigenvalues than requested
  !
  IF( M.GT.NUME ) THEN
     M = NUME
  END IF
  !
  INFO = 0
  IF( NCNVRG ) &
       INFO = INFO + 1
  IF( TOOFEW ) &
       INFO = INFO + 2
  RETURN
  !
  !     End of DSTEZ2
  !
END SUBROUTINE DSTEZ2

SUBROUTINE DUP2BP( UPLO, N, AP, HB )
  IMPLICIT NONE
  !
  !  DP2BP transforms a unpacked stored matrix AP into a block packed matrix.
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER*1        UPLO
  INTEGER            N, HB
  REAL*8             AP(*)
  !     
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix AP is stored as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of AP
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of AP
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A. 
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - REAL*8
  !           On entry, AP stores the matrix in unpacked storage format.
  !           On exit, AP stores the matrix in block packed storage.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            IJ2K
  EXTERNAL           LSAME, IJ2K
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            INFO, I, J, I1J
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  UPPER = LSAME( UPLO, 'U' )
  !
  INFO = 0
  IF(( .NOT.UPPER)) THEN
     INFO = 1
  ELSE IF( N  .LT.0 )THEN
     INFO = 2
  ELSE IF( HB .LT.0 ) THEN
     INFO = 4
  END IF
  !
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'DUP2BP ', INFO )
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF ( UPPER ) THEN
     DO J = 1,N   ! 20
        I1J = IJ2K( UPLO, 1, J, N, HB)
        DO I = 1, J  ! 10
           AP(I1J+I-1) = AP(I+(J-1)*N)
        ENDDO ! 10         CONTINUE
     ENDDO ! 20   CONTINUE
  END IF
  !
  !     End of DUP2BP
  !
END SUBROUTINE DUP2BP

INTEGER FUNCTION IJ2K( UPLO, I, J, N, HB )
  !
  !  IJ2K computes the index k of the matrix element A(I,J) in the 
  !  corresponding (block) packed array AP.
  !
  !     .. Scalar Arguments ..
  CHARACTER*1        UPLO
  INTEGER            I, J, N, HB
  !     
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix A is referenced as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of A
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of A
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  I      - INTEGER.
  !           On entry, I specifies the row of the matrix A.  
  !           Unchanged on exit.
  !
  !  J      - INTEGER.
  !           On entry, J specifies the column of the matrix A.  
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A. 
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            INFO, ITMP
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  UPPER = LSAME( UPLO, 'U' )
  !
  !      write(*,*) 'IJ2K: i,j,n,hb= ',i,j,n,hb
  INFO = 0
  IF(( .NOT.UPPER).AND. ( .NOT.LSAME( UPLO , 'L' ) ) )THEN
     INFO = 1
  ELSE IF( N  .LT.0 )THEN
     INFO = 4
  ELSE IF(( I  .LT.0 ) .OR. ( I .GT. N ))THEN
     INFO = 2
  ELSE IF(( J  .LT.0 ) .OR. ( J .GT. N ))THEN
     INFO = 3
  ELSE IF (( UPPER ) .AND. ( I .GT. J )) THEN
     INFO = 3
  END IF
  !
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'IJ2K ', INFO )
     RETURN
  END IF
  !
  !     Start the operations.
  !
  ITMP = MOD(J-1,HB)
  IF ( UPPER ) THEN
     IJ2K = I + HB*HB*((J-1)/HB)*((J-1)/HB+1)/2 + ITMP*MIN(((J-1)/HB+1)*HB,N)
  ELSE
     IJ2K = I+(2*N-J)*(J-1)/2 + (HB*(HB-1)/2) * ((J-1)/HB) + ITMP*(ITMP+1)/2
  END IF
  RETURN
  !
  !     End of IJ2K
  !
END FUNCTION IJ2K

SUBROUTINE ZBP2P( UPLO, N, AP, HB )
  IMPLICIT NONE
  !
  !  ZBP2P transforms a block packed stored matrix AP into a packed matrix.
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER*1        UPLO
  INTEGER            N, HB
  COMPLEX*16         AP(*)
  !     
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix AP is stored as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of AP
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of AP
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A. 
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - COMPLEX*16
  !           On entry, AP stores the matrix in packed storage format.
  !           In order to store the (not referenced elements in the 
  !           block diagonal), AP must be at least of size 
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !           On exit, AP stores the matrix in block packed storage.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            IJ2K
  EXTERNAL           LSAME, IJ2K
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            INFO, I, I1, I2, I3, J
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  UPPER = LSAME( UPLO, 'U' )
  !
  INFO = 0
  IF(( .NOT.UPPER).AND.( .NOT.LSAME( UPLO , 'L' ) ) )THEN
     INFO = 1
  ELSE IF( N  .LT.0 )THEN
     INFO = 2
  ELSE IF( HB .LT.0 ) THEN
     INFO = 4
  END IF
  !
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'ZBP2P ', INFO )
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF ( UPPER ) THEN
     I1 = 1
     I2 = 1
     I3 = 0
     DO J = 1,N      ! 20
        DO I = 1,J   ! 10
           AP(I1) = AP(I2)
           I1 = I1 +1
           I2 = I2 +1
        ENDDO !10         CONTINUE
        IF (J .GT. ((N/HB)*HB)) THEN
           I3 = N-J
        ELSE
           I3 = MOD(I3+HB-1 , HB)
        END IF
        I2 = I2 + I3
     ENDDO !   20    CONTINUE
  ELSE
     I1 = 1
     I2 = 1
     I3 = 0
     DO J = 1,N         ! 40
        DO I = J,N      ! 30
           AP(I1) = AP(I2)
           I1 = I1 +1
           I2 = I2 +1
        ENDDO ! 30         CONTINUE
        I3 = MOD(I3+1 , HB)
        I2 = I2 + I3
     ENDDO !40      CONTINUE
  END IF
  !
  !     End of ZBP2P
  !
END SUBROUTINE ZBP2P

SUBROUTINE ZHHEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK,IFAIL, HBTRD, NUME, NBVL, INFO )
  IMPLICIT NONE
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER, intent(in)  :: JOBZ, RANGE, UPLO
  INTEGER,   intent(in)  :: IL, IU, LDZ, N, HBTRD, NUME
  INTEGER,   intent(out) :: INFO, NBVL, M
  DOUBLE PRECISION, intent(in) ::  ABSTOL, VL, VU
  !     ..
  !     .. Array Arguments ..
  DOUBLE PRECISION, intent(out)   :: W( * )
  COMPLEX*16,       intent(out)   :: Z(LDZ, * )
  INTEGER,          intent(out)   :: IFAIL( * )
  INTEGER,          intent(inout) :: IWORK( * )
  DOUBLE PRECISION, intent(inout) :: RWORK( * )
  COMPLEX*16,       intent(inout) :: AP( * ), WORK( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  ZHHEVX computes selected eigenvalues and, optionally, eigenvectors
  !  of a complex Hermitian matrix A in packed storage, which is converted
  !  for computation to block packed storage.
  !  Eigenvalues/vectors can be selected by specifying either a range of
  !  values or a range of indices for the desired eigenvalues.
  !
  !  Arguments
  !  =========
  !
  !  JOBZ    (input) CHARACTER*1
  !          = 'N':  Compute eigenvalues only;
  !          = 'V':  Compute eigenvalues and eigenvectors.
  !
  !  RANGE   (input) CHARACTER*1
  !          = 'A': all eigenvalues will be found;
  !          = 'V': all eigenvalues in the half-open interval (VL,VU]
  !                 will be found;
  !          = 'I': the IL-th through IU-th eigenvalues will be found.
  !          = 'C': the eigenvalues in the half-open interval (VL,VU]
  !                 will be found, if there are more than M Eigenvalues,
  !                 only the smallest M Eigenvalues are computed;
  !
  !  UPLO    (input) CHARACTER*1
  !          = 'U':  Upper triangle of A is stored;
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  AP      (input/output) COMPLEX*16 array,
  !          dimension (N*(N+1)/2+(N/HBTRD+1)*HBTRD*(HBTRD-1)/2)
  !          On entry, the upper triangle of the Hermitian matrix
  !          A, packed columnwise in a linear array.  The j-th column of A
  !          is stored in the array AP as follows:
  !          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
  !
  !          On exit, AP is overwritten by values generated during the
  !          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
  !          and first superdiagonal of the tridiagonal matrix T overwrite
  !          the corresponding elements of A.
  !
  !  VL      (input) DOUBLE PRECISION
  !  VU      (input) DOUBLE PRECISION
  !          If RANGE='V' or RANGE'C', the lower and upper bounds of the
  !          interval to be searched for eigenvalues. VL < VU.
  !          Not referenced if RANGE = 'A' or 'I'.
  !
  !  IL      (input) INTEGER
  !  IU      (input) INTEGER
  !          If RANGE='I', the indices (in ascending order) of the
  !          smallest and largest eigenvalues to be returned.
  !          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
  !          Not referenced if RANGE = 'A' or 'V'.
  !
  !  ABSTOL  (input) DOUBLE PRECISION
  !          The absolute error tolerance for the eigenvalues.
  !          An approximate eigenvalue is accepted as converged
  !          when it is determined to lie in an interval [a,b]
  !          of width less than or equal to
  !
  !                  ABSTOL + EPS *   max( |a|,|b| ) ,
  !
  !          where EPS is the machine precision.  If ABSTOL is less than
  !          or equal to zero, then  EPS*|T|  will be used in its place,
  !          where |T| is the 1-norm of the tridiagonal matrix obtained
  !          by reducing AP to tridiagonal form.
  !
  !          Eigenvalues will be computed most accurately when ABSTOL is
  !          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
  !          If this routine returns with INFO>0, indicating that some
  !          eigenvectors did not converge, try setting ABSTOL to
  !          2*DLAMCH('S').
  !
  !          See "Computing Small Singular Values of Bidiagonal Matrices
  !          with Guaranteed High Relative Accuracy," by Demmel and
  !          Kahan, LAPACK Working Note #3.
  !
  !  M       (output) INTEGER
  !          The total number of eigenvalues found.  0 <= M <= N.
  !          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1, and
  !          if RANGE = 'C', M <= NUME.
  !
  !  W       (output) DOUBLE PRECISION array, dimension (N)
  !          If INFO = 0, the selected eigenvalues in ascending order.
  !
  !  Z       (output) COMPLEX*16 array, dimension (LDZ, max(1,M))
  !          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
  !          contain the orthonormal eigenvectors of the matrix A
  !          corresponding to the selected eigenvalues, with the i-th
  !          column of Z holding the eigenvector associated with W(i).
  !          If an eigenvector fails to converge, then that column of Z
  !          contains the latest approximation to the eigenvector, and
  !          the index of the eigenvector is returned in IFAIL.
  !          If JOBZ = 'N', then Z is not referenced.
  !          Note: the user must ensure that at least max(1,M) columns are
  !          supplied in the array Z; if RANGE = 'V', the exact value of M
  !          is not known in advance and an upper bound must be used.
  !
  !  LDZ     (input) INTEGER
  !          The leading dimension of the array Z.  LDZ >= 1, and if
  !          JOBZ = 'V', LDZ >= max(1,N).
  !
  !  WORK    (workspace) COMPLEX*16 array, dimension ((2+HBTRD)*N)
  !
  !  RWORK   (workspace) DOUBLE PRECISION array, dimension (7*N)
  !
  !  IWORK   (workspace) INTEGER array, dimension (5*N)
  !
  !  IFAIL   (output) INTEGER array, dimension (N)
  !          If JOBZ = 'V', then if INFO = 0, the first M elements of
  !          IFAIL are zero.  If INFO > 0, then IFAIL contains the
  !          indices of the eigenvectors that failed to converge.
  !          If JOBZ = 'N', then IFAIL is not referenced.
  !
  !  HBTRD   (input) INTEGER
  !          storage block size of AP, used in ZHHTRD.
  !
  !  NUME    (input) INTEGER
  !          The maximum number of Eigenvalues allowed
  !          (dimension of Z)
  !
  !  NBVL    (output) INTEGER
  !          The number of Eigenvalues below lower bound VL
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !          > 0:  if INFO = i, then i eigenvectors failed to converge.
  !                Their indices are stored in array IFAIL.
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  DOUBLE PRECISION   ZERO, ONE
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
  COMPLEX*16         CONE
  PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            ALLEIG, INDEIG, VALEIG, WANTZ
  CHARACTER          ORDER
  INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWK, INDRWK, INDTAU, INDWRK, ISCALE, ITMP1, J, JJ, NSPLIT, LLWORK
  DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  DOUBLE PRECISION   DLAMCH, ZLANHP
  EXTERNAL           LSAME, DLAMCH, ZLANHP
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DCOPY, DSCAL, DSTEZ2, DSTERF, XERBLA, ZDSCAL,ZHHTRD, ZSTEIN, ZSTEQR, ZSWAP, ZUPGTR, ZUPMTR
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          DBLE, MIN, SQRT
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  WANTZ = LSAME( JOBZ, 'V' )
  ALLEIG = LSAME( RANGE, 'A' )
  VALEIG = LSAME( RANGE, 'V' ) .OR. LSAME( RANGE, 'C' )
  INDEIG = LSAME( RANGE, 'I' )
  !
  INFO = 0
  IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
     INFO = -1
  ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
     INFO = -2
  ELSE IF( .NOT.LSAME( UPLO, 'U' ) ) THEN
     INFO = -3
  ELSE IF( N.LT.0 ) THEN
     INFO = -4
  ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
     INFO = -7
  ELSE IF( INDEIG .AND. IL.LT.1 ) THEN
     INFO = -8
  ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) THEN
     INFO = -9
  ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
     INFO = -14
  END IF
  !
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZHHEVX', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  M = 0
  IF( N.EQ.0 ) RETURN
  !
  IF( N.EQ.1 ) THEN
     IF( ALLEIG .OR. INDEIG ) THEN
        M = 1
        W( 1 ) = AP( 1 )
     ELSE
        IF( VL.LT.DBLE( AP( 1 ) ) .AND. VU.GE.DBLE( AP( 1 ) ) ) THEN
           M = 1
           W( 1 ) = AP( 1 )
        END IF
     END IF
     IF( WANTZ ) Z( 1, 1 ) = CONE
     RETURN
  END IF
  !
  !     Get machine constants.
  !
  SAFMIN = DLAMCH( 'Safe minimum' )
  EPS = DLAMCH( 'Precision' )
  SMLNUM = SAFMIN / EPS
  BIGNUM = ONE / SMLNUM
  RMIN = SQRT( SMLNUM )
  RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
  !
  !     Scale matrix to allowable range, if necessary.
  !
  ISCALE = 0
  ABSTLL = ABSTOL
  IF( VALEIG ) THEN
     VLL = VL
     VUU = VU
  END IF
  ANRM = ZLANHP( 'M', UPLO, N, AP, RWORK )
  IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
     ISCALE = 1
     SIGMA = RMIN / ANRM
  ELSE IF( ANRM.GT.RMAX ) THEN
     ISCALE = 1
     SIGMA = RMAX / ANRM
  END IF
  IF( ISCALE.EQ.1 ) THEN
     CALL ZDSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
     IF( ABSTOL.GT.0 ) &
          ABSTLL = ABSTOL*SIGMA
     IF( VALEIG ) THEN
        VLL = VL*SIGMA
        VUU = VU*SIGMA
     END IF
  END IF
  !
  !     Call ZHHTRD to reduce Hermitian packed matrix to tridiagonal form.
  !
  INDD = 1
  INDE = INDD + N
  INDRWK = INDE + N
  INDTAU = 1
  INDWRK = INDTAU + N
  LLWORK = N*HBTRD
  CALL ZHHTRD( UPLO, N, AP, RWORK( INDD ), RWORK( INDE ),WORK( INDTAU ), WORK( INDWRK ), LLWORK, HBTRD,IINFO )
  !
  !     If all eigenvalues are desired and ABSTOL is less than or equal
  !     to zero, then call DSTERF or ZUPGTR and ZSTEQR.  If this fails
  !     for some eigenvalue, then try DSTEZ2.
  !
  IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.EQ.N ) ) .AND.( ABSTOL.LE.ZERO ) ) THEN
     CALL DCOPY( N, RWORK( INDD ), 1, W, 1 )
     INDEE = INDRWK + 2*N
     IF( .NOT.WANTZ ) THEN
        CALL DCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
        CALL DSTERF( N, W, RWORK( INDEE ), INFO )
     ELSE
        CALL ZUPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ,WORK( INDWRK ), IINFO )
        CALL DCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
        CALL ZSTEQR( JOBZ, N, W, RWORK( INDEE ), Z, LDZ,RWORK( INDRWK ), INFO )
        IF( INFO.EQ.0 ) THEN
           DO I = 1, N   ! 10
              IFAIL( I ) = 0
           ENDDO         ! 10          CONTINUE
        END IF
     END IF
     IF( INFO.EQ.0 ) THEN
        M = N
        GO TO 20
     END IF
     INFO = 0
  END IF
  !
  !     Otherwise, call DSTEZ2 and, if eigenvectors are desired, ZSTEIN.
  !
  IF( WANTZ ) THEN
     ORDER = 'B'
  ELSE
     ORDER = 'E'
  END IF
  INDIBL = 1
  INDISP = INDIBL + N
  INDIWK = INDISP + N
  IF(WANTZ .AND. (RANGE.EQ.'C') ) THEN
     M=NUME
  END IF
  CALL DSTEZ2( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL,RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W,IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ),IWORK( INDIWK ), NUME, NBVL, INFO )
  !
  IF( WANTZ ) THEN
     CALL ZSTEIN( N, RWORK( INDD ), RWORK( INDE ), M, W,IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO )
     !	 write(*,*) M,' Eigenvectors computed'
     !         write(*,*) NBVL, ' Eigenvalues below lower bound'
     !
     !        Apply unitary matrix used in reduction to tridiagonal
     !        form to eigenvectors returned by ZSTEIN.
     !
     INDWRK = INDTAU + N
     CALL ZUPMTR( 'L', UPLO, 'N', N, M, AP, WORK( INDTAU ), Z, LDZ,WORK( INDWRK ), INFO )
  END IF
  !
  !     If matrix was scaled, then rescale eigenvalues appropriately.
  !
20 CONTINUE
  IF( ISCALE.EQ.1 ) THEN
     IF( INFO.EQ.0 ) THEN
        IMAX = M
     ELSE
        IMAX = INFO - 1
     END IF
     CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
  END IF
  !
  !     If eigenvalues are not in order, then sort them, along with
  !     eigenvectors.
  !
  IF( WANTZ ) THEN
     DO J = 1, M - 1        ! 40
        I = 0
        TMP1 = W( J )
        DO JJ = J + 1, M    ! 30
           IF( W( JJ ).LT.TMP1 ) THEN
              I = JJ
              TMP1 = W( JJ )
           END IF
        ENDDO               ! 30         CONTINUE
        !
        IF( I.NE.0 ) THEN
           ITMP1 = IWORK( INDIBL+I-1 )
           W( I ) = W( J )
           IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
           W( J ) = TMP1
           IWORK( INDIBL+J-1 ) = ITMP1
           CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
           IF( INFO.NE.0 ) THEN
              ITMP1 = IFAIL( I )
              IFAIL( I ) = IFAIL( J )
              IFAIL( J ) = ITMP1
           END IF
        END IF
     ENDDO !  40    CONTINUE
  END IF
  !
  RETURN
  !
  !     End of ZHHEVX
  !
END SUBROUTINE ZHHEVX

SUBROUTINE ZHHGS4( ITYPE, UPLO, N, ABP, BBP, HB, INFO )
  IMPLICIT NONE
  ! 
  !
  !     .. Scalar Arguments ..
  CHARACTER          UPLO
  INTEGER            INFO, ITYPE, N, HB
  !     ..
  !     .. Array Arguments ..
  COMPLEX*16   ABP( * ), BBP( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  ZHHGS4 reduces a real Hermitsch-definite generalized eigenproblem
  !  to standard form. The arrays ABP and BBP are stored in block packed 
  !  mode. All blocks of the diagonal are of order HB (the order of the 
  !  last one may be less).
  !
  !  If ITYPE = 1, the problem is A*x = lambda*B*x,
  !  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
  !
  !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
  !  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
  !
  !  Only ITYPE = 1 is implemented.
  !
  !  B must have been previously factorized as U**T*U or L*L**T by ZPHTRF.
  !
  !  Arguments
  !  =========
  !
  !  ITYPE   (input) INTEGER
  !          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
  !          = 2 or 3: compute U*A*U**T or L**T*A*L, not implemented.
  !
  !  UPLO    (input) CHARACTER
  !          = 'U':  Upper triangle of A is stored and B is factored as
  !                  U**T*U;
  !          = 'L':  Lower triangle of A is stored and B is factored as
  !                  L*L**T, not implemented.
  !
  !  N       (input) INTEGER
  !          The order of the matrices A and B.  N >= 0.
  !
  !  ABP     (input/output) COMPLEX*16 array.
  !          Dimension (N*(N+1)/2+(N/HB+1)*HB*(HB-1)/2).
  !          On entry, the symmetric matrix A.  
  !
  !          On exit, if INFO = 0, the transformed matrix, stored in the
  !          same format as A.
  !
  !  BBP     (input) COMPLEX*16 array.
  !          Dimension (N*(N+1)/2+(N/HB+1)*HB*(HB-1)/2).
  !          The triangular factor from the Cholesky factorization of B,
  !          as returned by ZPPTRFHB.
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  COMPLEX*16   ONE, HALF
  PARAMETER          ( ONE = (1.0D0,0), HALF = (0.5D0,0) )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            K, KB, KK, KKB, KKK, KKKB
  INTEGER            IKK, IKK2, IK2K2, IKK3, IK3K2
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DSYGS2, DSYMM, DSYR2K, DTRMM, DTRSM, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            ILAENV, IJ2K
  EXTERNAL           LSAME, ILAENV
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( ITYPE.NE.1 ) THEN
     INFO = -1
  ELSE IF( .NOT.UPPER ) THEN
     INFO = -2
  ELSE IF( N.LT.0 ) THEN
     INFO = -3
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZHHGS4', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) RETURN
  !
  IF( HB.LE.1 ) THEN
     CALL ZHPGST( ITYPE, UPLO, N, ABP, BBP, INFO)
  ELSE IF( HB.GE.N ) THEN
     CALL ZHEGST( ITYPE, UPLO, N, ABP, N, BBP, N, INFO )
  ELSE
     !
     !     Determine the block size for this environment.
     !
     !        NB = ILAENV( 1, 'ZHEGST', UPLO, N, -1, -1, -1 )
     !        IF( HB.LE.NB ) 
     !     $    BLOCKED = .TRUE.
     !
     !        Use hierarchically blocked code
     !
     IF( ITYPE.EQ.1 ) THEN
        IF( UPPER ) THEN
           !
           !              Compute inv(U')*A*inv(U)
           !
           DO K = 1, N, HB   ! 10
              !		  IF( HB.LE.(N-K+1) ) THEN
              !		    KB = HB
              !		  ELSE
              !		    KB = N-K+1
              !		    BLOCKED = .FALSE.
              !		  END IF
              KB = MIN( N-K+1, HB )
              !
              !                 Update the upper triangle of A(k:n,k:n)
              !
              IKK = IJ2K( UPLO, K, K, N, HB )
              !		  write (*,*) 'zhegst: ikk,k,kb=',ikk,k,kb 
              CALL ZHEGST( ITYPE, UPLO, KB, ABP( IKK ), K+KB-1,BBP( IKK ), K+KB-1, INFO )
              IF( INFO.NE.0) CYCLE !GO TO 10
              IF( K+KB.LE.N ) THEN
                 DO KK = K+KB, N, HB    ! 9
                    KKB = MIN( N-KK+1, HB )
                    IKK2 = IJ2K( UPLO, K, KK, N, HB )
                    IK2K2 = IJ2K( UPLO, KK, KK, N, HB )
                    CALL ZTRSM( 'Left', UPLO, 'Conj', 'Non-unit',KB, KKB, ONE, BBP( IKK ), K+KB-1,ABP( IKK2 ), KK+KKB-1 )
                    CALL ZHEMM( 'Left', UPLO, KB, KKB, -HALF,ABP( IKK ), K+KB-1, BBP( IKK2 ),KK+KKB-1,ONE, ABP( IKK2 ), KK+KKB-1 )
                 ENDDO ! 9                   CONTINUE
                 !
                 !                       dsyr2k block-wise
                 !
                 DO KK = K+KB, N, HB        ! 8
                    KKB = MIN( N-KK+1, HB )
                    IKK2 = IJ2K( UPLO, K, KK, N, HB )
                    IK2K2 = IJ2K( UPLO, KK, KK, N, HB )
                    IF (  KK .GT. K+KB ) THEN
                       DO KKK = K+KB, KK-KKB, HB   ! 2
                          KKKB = MIN( N-KKK+1, HB )
                          IKK3 = IJ2K( UPLO, K, KKK, N, HB )
                          IK3K2 = IJ2K( UPLO, KKK, KK, N, HB )
                          CALL ZGEMM( 'Conjugate', 'No transpose',HB, KKB, HB, -ONE, ABP( IKK3 ), KKK+KKKB-1,BBP( IKK2 ), KK+KKB-1,ONE, ABP( IK3K2 ), KK+KKB-1)
                          CALL ZGEMM( 'Conjugate', 'No transpose',HB, KKB, HB,-ONE, BBP( IKK3 ), KKK+KKKB-1,ABP( IKK2 ), KK+KKB-1,ONE, ABP( IK3K2 ), KK+KKB-1)
                       ENDDO                       !  2                         CONTINUE
                    ENDIF
                    CALL ZHER2K( UPLO, 'Conjugate', KKB, KB, -ONE,ABP( IKK2 ), KK+KKB-1, BBP( IKK2 ),KK+KKB-1,ONE, ABP( IK2K2 ), KK+KKB-1 )
                 ENDDO !  8                 CONTINUE
                 !123456789012345678901234567890123456789012345678901234567890123456789012
                 DO KK = K+KB, N, HB    ! 5
                    KKB = MIN( N-KK+1, HB )
                    IKK2 = IJ2K( UPLO, K, KK, N, HB )
                    IK2K2 = IJ2K( UPLO, KK, KK, N, HB )
                    CALL ZHEMM( 'Left', UPLO, KB, KKB, -HALF,ABP( IKK ), K+KB-1, BBP( IKK2 ),KK+KKB-1,ONE, ABP( IKK2 ), KK+KKB-1 )
                    !
                    !                       Reduce triangular system
                    !
                    IF (  KK .GT. K+KB ) THEN
                       DO KKK = K+KB, KK-KKB, HB  ! 3
                          KKKB = MIN( N-KKK+1, HB )
                          IKK3 = IJ2K( UPLO, K, KKK, N, HB )
                          IK3K2 = IJ2K( UPLO, KKK, KK, N, HB )
                          CALL ZGEMM( 'No trans', 'No trans', HB, KKB, HB,-ONE, ABP( IKK3 ), KKK+KKKB-1,BBP( IK3K2 ), KK+KKB-1,ONE, ABP( IKK2 ), KK+KKB-1)
                       ENDDO !    3                      CONTINUE
                    END IF
                    !
                    !                       compute block-column of reduced system
                    !
                    CALL ZTRSM( 'Right', UPLO, 'No transpose', 'Non-unit', KB, KKB, ONE, BBP( IK2K2 ), KK+KKB-1, ABP( IKK2 ),KK+KKB-1 )
                    !123456789012345678901234567890123456789012345678901234567890123456789012
                 ENDDO ! 5                   CONTINUE
              END IF
           ENDDO   ! 10            CONTINUE
        END IF
     END IF
  END IF
  RETURN
  
20 CONTINUE
  INFO = INFO + K - 1
  RETURN
  !
  !     End of ZHHGS4
  !
END SUBROUTINE ZHHGS4

SUBROUTINE ZHHGST( ITYPE, UPLO, N, AP, BP, HB, INFO )
  IMPLICIT NONE
  !
  !
  CHARACTER      UPLO
  INTEGER        INFO, N, ITYPE, HB
  COMPLEX*16     AP( * ), BP( * )
  !
  !  Purpose
  !  =======
  !
  !  ZHHGST reduces a real symmetric-definite generalized eigenproblem
  !  to standard form.
  !
  !  If ITYPE = 1, the problem is A*x = lambda*B*x,
  !  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
  !
  !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
  !  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
  !
  !  B must have been previously factorized as U**T*U or L*L**T by ZPHTRF.
  !
  !  Arguments
  !  =========
  !
  !  ITYPE   (input) INTEGER
  !          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
  !          = 2 or 3: compute U*A*U**T or L**T*A*L, not implemented.
  !
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix AP is stored as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of AP
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of AP
  !                                  is to be referenced, not implemented.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - REAL*8
  !           On entry, AP stores the matrix in packed storage format.
  !           In order to store the (not referenced elements in the
  !           block diagonal), AP must be at least of size
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !           On exit, AP stores the matrix in block packed storage.
  !
  !  BP     - REAL*8
  !           On entry, BP stores the triangular factor from the 
  !           Cholesky factorization of B, as returned by ZPPTRFHB.
  !           In order to store the (not referenced elements in the
  !           block diagonal), BP must be at least of size
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !  INFO  - (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order K is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !  =====================================================================
  !           
  REAL*8             CP(4), FLOP
  integer j
  logical UPPER
  logical LSAME
  
  integer ilaenv
  external ilaenv, LSAME
  
  !
  !     Test the input parameters.
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
     INFO = -1
  ELSE IF( .NOT.UPPER ) THEN
     INFO = -2
  ELSE IF( N.LT.0 ) THEN
     INFO = -3
  ELSE IF( HB.LT.0 ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZHHGST', -INFO )
     RETURN
  END IF
  
  IF ( HB.le.1 ) THEN
     call ZHPGST( ITYPE, UPLO, N, AP, BP, INFO ) 
  ELSE
     call CPUTIM(CP(1))
     !
     !      call the subroutines
     !			
     !
     !      convert from LAPACK upper triangular packed storage format
     !                     into a upper triangular block packed matrix
     !
     call zp2bp('U',N,AP,HB)
     call zp2bp('U',N,BP,HB)
     call CPUTIM(CP(2))
     !
     !      reduce the real symmetrix-definite generalized eigenproblem
     !      to standard form
     !
     call ZHHGS4(ITYPE, UPLO, N, AP, BP, HB, INFO)
     call CPUTIM(CP(3))
     IF (INFO .NE. 0) THEN
        CALL OUTERR('zhhgst','ZHHGS4 aborted unsuccessfully.')
        GOTO 999
     ENDIF
     !
     !      convert back to upper triangular packed storage format
     !
     call zbp2p('U',N,AP,HB)
     call zbp2p('U',N,BP,HB)
     call CPUTIM(CP(4))
     
     !
     !      timing output
     !
     DO j = 1, 3  ! 30
        CP(j) = CP(j+1)-CP(j)
     ENDDO !30      CONTINUE
     Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Packing' , CP(1) , 1E-6*Flop/CP(1)
     Flop = 0
     !         WRITE (*,1001) 'Transformation' , CP(2) , 1E-6*Flop/CP(2)
     Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Unpacking' , CP(3) , 1E-6*Flop/CP(3)
1001 FORMAT (1X,'ZHBGSTHB(',A,') :',t30,f7.3:t40,f8.3,' Mflops') 
     !	 write (*,*) 'INFO = ', INFO 
  END IF
  RETURN 
  !	
  !	Result stored using LAPACK packed format 
  ! 
  ! 
999 STOP 'LAPW1 - Error'
  !
END SUBROUTINE ZHHGST

SUBROUTINE ZHHTR4( UPLO, N, ABP, D, E, TAU, WORK, LWORK, HB, INFO )
  IMPLICIT NONE
  !
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER, intent(in) :: UPLO
  INTEGER,   intent(out):: INFO
  INTEGER,   intent(in) :: LWORK, N, HB
  !     ..
  !     .. Array Arguments ..
  DOUBLE PRECISION, intent(out) :: D( * ), E( * )
  COMPLEX*16,     intent(inout) :: ABP( * ), WORK( * )
  COMPLEX*16,     intent(out)   :: TAU( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  ZHHTR4 reduces a complex Hermitian matrix A to real symmetric
  !  tridiagonal form T by an unitary similarity transformation:
  !  Q**T * A * Q = T.
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER*1
  !          = 'U':  Upper triangle of A is stored;
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  ABP     (input/output) COMPLEX*16 array, 
  !          dimension (n*(n+1)/2+(n/hb+1)*hb*(hb-1)/2)
  !          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !          N-by-N upper triangular part of A contains the upper
  !          triangular part of the matrix A, and the strictly lower
  !          triangular part of A is not referenced.  
  !          On exit, if UPLO = 'U', the diagonal and first superdiagonal
  !          of A are overwritten by the corresponding elements of the
  !          tridiagonal matrix T, and the elements above the first
  !          superdiagonal, with the array TAU, represent the orthogonal
  !          matrix Q as a product of elementary reflectors.
  !          See Further Details.
  !
  !  D       (output) DOUBLE PRECISION array, dimension (N)
  !          The diagonal elements of the tridiagonal matrix T:
  !          D(i) = A(i,i).
  !
  !  E       (output) DOUBLE PRECISION array, dimension (N-1)
  !          The off-diagonal elements of the tridiagonal matrix T:
  !          E(i) = A(i,i+1) if UPLO = 'U'.
  !
  !  TAU     (output) COMPLEX*16 array, dimension (N-1)
  !          The scalar factors of the elementary reflectors (see Further
  !          Details).
  !
  !  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
  !
  !  LWORK   (input) INTEGER
  !          LWORK = N*HB, where HB is the blocksize.
  !
  !  HB      (input) INTEGER
  !          The blocksize used in packed storage scheme.
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !
  !  Further Details
  !  ===============
  !
  !  If UPLO = 'U', the matrix Q is represented as a product of elementary
  !  reflectors
  !
  !     Q = H(n-1) . . . H(2) H(1).
  !
  !  Each H(i) has the form
  !
  !     H(i) = I - tau * v * v'
  !
  !  where tau is a real scalar, and v is a real vector with
  !  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
  !  A(1:i-1,i+1), and tau in TAU(i).
  !
  !  The contents of A on exit are illustrated by the following examples
  !  with n = 5:
  !
  !  if UPLO = 'U':               
  !
  !    (  d   e   v2  v3  v4 )     
  !    (      d   e   v3  v4 )    
  !    (          d   e   v4 )   
  !    (              d   e  )  
  !    (                  d  ) 
  !
  !  where d and e denote diagonal and off-diagonal elements of T, and vi
  !  denotes an element of the vector defining H(i).
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  COMPLEX*16   ONE
  PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            I, IB, J, LDWORK, NB
  INTEGER            I1I, I1J, IJI, IJJ, IJ1J, LDI, LDJ
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           ZLATD4, ZHER2K, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            ILAENV
  real               dtime2
  integer            IJ2K
  EXTERNAL           LSAME, ILAENV, IJ2K
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LWORK.LT.(N*HB) ) THEN
     INFO = -8
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZHHTR4', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) THEN
     WORK( 1 ) = 1
     RETURN
  END IF
  !
  !     Determine the block size.
  !
  !      lasttime=0
  !      write(*,*) 'ilaenv(DSYTRD)=', NB
  !
  NB=MIN(HB,N)
  LDWORK = N
  !
  IF( UPPER ) THEN
     !
     !
     !         T1=0
     !	 T2=0
     DO I= NB*((N-1)/NB)+1, 1, -NB    ! 20
        IB = MIN(NB, N-I+1)
        LDI = I+IB-1
        !           write(*,*) 'I=',I,' IB=',IB,' N=',N,' NB=',NB,
        !     $                ' LDWORK',LDWORK
        !
        !           Reduce columns i:i+ib-1 to tridiagonal form and form the
        !           matrix W which is needed to update the unreduced part of
        !           the matrix
        !
        CALL ZLATD4( UPLO, LDI, IB, ABP, E, TAU, WORK, LDWORK, NB )
        !STOP
        !	    T1 = T1 + dtime2(lasttime)
        !
        !           Update the unreduced submatrix A(1:i-1,1:i-1), using an
        !           update of the form:  A := A - V*W' - W*V'
        !
        I1I = IJ2K(UPLO, 1, I, N, NB)
        DO J = 1, I-1, NB   ! 5
           LDJ = J+NB-1
           !              write(*,*) 'ZHER2K I=',I,' J=',J
           I1J = IJ2K(UPLO, 1, J, N, NB)
           IJI = IJ2K(UPLO, J, I, N, NB)
           IJJ = IJ2K(UPLO, J, J, N, NB)
           IF( J.NE.1 ) THEN
              CALL ZGEMM( 'NoTranspose', 'Conjugate', J-1, NB, IB,-ONE, ABP(I1I), LDI, WORK(J), LDWORK,ONE, ABP(I1J), LDJ)
              CALL ZGEMM( 'NoTranspose', 'Conjugate', J-1, NB, IB,-ONE, WORK, LDWORK, ABP(IJI), LDI, ONE, ABP(I1J), LDJ)
           END IF
           CALL ZHER2K( UPLO, 'NoTranspose', NB, IB,-ONE, ABP(IJI), LDI, WORK(J), LDWORK, ONE, ABP(IJJ), LDJ)
           
        ENDDO  !5          CONTINUE
        !	    T2 = T2 + dtime2(lasttime)
        !
        !           Copy superdiagonal elements back into A, and diagonal
        !           elements into D
        !
        DO J = I, I + IB - 1    ! 10
           !	       write(*,*) 'Copy: I=',I,' J=',J
           IJJ = IJ2K(UPLO, J, J, N, NB)
           IF( J.GT.1 ) THEN
              IJ1J = IJ2K(UPLO, J-1, J, N, NB)
              ABP(IJ1J) = E( J-1 )
           END IF
           D( J ) = ABP( IJJ )
        ENDDO !10         CONTINUE
     ENDDO !20    CONTINUE
  END IF
  !
  !      write(*,*) '           DLATD4: ', T1
  !      write(*,*) '           DSYR2K: ', T2
  RETURN
  !
  !     End of ZHHTR4
  !
END SUBROUTINE ZHHTR4


SUBROUTINE ZHHTRD(UPLO,N,ABP,D,E,TAU,WORK,LWORK,HB,INFO)
  IMPLICIT NONE
  !
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER, intent(in) :: UPLO
  INTEGER,   intent(out):: INFO
  INTEGER,   intent(in) :: LWORK, N, HB
  !     ..
  !     .. Array Arguments ..
  DOUBLE PRECISION, intent(out) :: D( * ), E( * )
  COMPLEX*16,       intent(out) :: TAU( * )
  COMPLEX*16,     intent(inout) :: ABP( * ), WORK( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  ZHHTRD reduces a complex Hermitian matrix A to real symmetric
  !  tridiagonal form T by an unitary similarity transformation:
  !  Q**T * A * Q = T.
  !
  !  The matrix A is stored in packed mode and will be restored in packed 
  !  mode on exit. 
  !  Computation is done using block packed storage mode.
  !
  !  Note: this routine has three parameters which do not exist in the 
  !        call to ZHPTRD:
  !        WORK, LWORK, HB
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER*1
  !          = 'U':  Upper triangle of A is stored;
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  ABP     (input/output) COMPLEX*16 array, 
  !          dimension (n*(n+1)/2+(n/hb+1)*hb*(hb-1)/2)
  !          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !          N-by-N upper triangular part of A contains the upper
  !          triangular part of the matrix A, and the strictly lower
  !          triangular part of A is not referenced.  
  !          On exit, if UPLO = 'U', the diagonal and first superdiagonal
  !          of A are overwritten by the corresponding elements of the
  !          tridiagonal matrix T, and the elements above the first
  !          superdiagonal, with the array TAU, represent the orthogonal
  !          matrix Q as a product of elementary reflectors.
  !          See Further Details.
  !
  !  D       (output) DOUBLE PRECISION array, dimension (N)
  !          The diagonal elements of the tridiagonal matrix T:
  !          D(i) = A(i,i).
  !
  !  E       (output) DOUBLE PRECISION array, dimension (N-1)
  !          The off-diagonal elements of the tridiagonal matrix T:
  !          E(i) = A(i,i+1) if UPLO = 'U'.
  !
  !  TAU     (output) COMPLEX*16 array, dimension (N-1)
  !          The scalar factors of the elementary reflectors (see Further
  !          Details).
  !
  !  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
  !
  !  LWORK   (input) INTEGER
  !          LWORK >= N*HB.
  !
  !  HB      (input) INTEGER
  !          The blocksize used in packed storage scheme.
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !
  !  Further Details
  !  ===============
  !
  !  If UPLO = 'U', the matrix Q is represented as a product of elementary
  !  reflectors
  !
  !     Q = H(n-1) . . . H(2) H(1).
  !
  !  Each H(i) has the form
  !
  !     H(i) = I - tau * v * v'
  !
  !  where tau is a real scalar, and v is a real vector with
  !  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
  !  A(1:i-1,i+1), and tau in TAU(i).
  !
  !  The contents of A on exit are illustrated by the following examples
  !  with n = 5:
  !
  !  if UPLO = 'U':               
  !
  !    (  d   e   v2  v3  v4 )     
  !    (      d   e   v3  v4 )    
  !    (          d   e   v4 )   
  !    (              d   e  )  
  !    (                  d  ) 
  !
  !  where d and e denote diagonal and off-diagonal elements of T, and vi
  !  denotes an element of the vector defining H(i).
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  DOUBLE PRECISION   ONE
  PARAMETER          ( ONE = 1.0D0 )
  COMPLEX*16         CONE
  PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
  INTEGER            MINHB
  PARAMETER          ( MINHB = 1 )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            NB
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           XERBLA, ZHER2K, ZHHTR4, ZLATD4
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            ILAENV
  real               dtime2
  integer            IJ2K
  EXTERNAL           LSAME, ILAENV, IJ2K
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LWORK.LT.(N*HB) ) THEN
     INFO = -8
  ELSE IF( HB.LT.1 ) THEN
     INFO = -9
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZHHTRD', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) THEN
     WORK( 1 ) = 1
     RETURN
  END IF
  !
  !     Determine the block size.
  !
  NB=MIN(HB,N)
  IF( NB.LE.MINHB ) THEN
     CALL ZHPTRD( UPLO, N, ABP, D, E, TAU, INFO )
  ELSE IF( UPPER ) THEN
     !
     !      convert from LAPACK upper triangular packed storage format
     !                     into a upper triangular block packed matrix
     !
     !	call CPUTIM(CP(1))
     call zp2bp('U',N,ABP,NB)
     !	call CPUTIM(CP(2))
     !      
     !      compute the tridiagonalization of A
     !
     call ZHHTR4( UPLO, N, ABP, D, E, TAU, WORK, LWORK, NB, INFO )
     !STOP
     !      call CPUTIM(CP(3))
     !
     !      convert back to upper triangular packed storage format
     !
     call zbp2p('U',N,ABP,NB)
     !      call CPUTIM(CP(4))
     
     !
     !      timing output
     !
     !         DO 30 j = 1, 3
     !           CP(j) = CP(j+1)-CP(j)
     ! 30      CONTINUE
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Packing' , CP(1) , 1E-6*Flop/CP(1)
     !         Flop = (N*DBLE(N+1)*(2*N+1))/3
     !         WRITE (*,1001) 'Tridiagonalization' , CP(2) , 1E-6*Flop/CP(2)
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Unpacking' , CP(3) , 1E-6*Flop/CP(3)
     ! 1001    FORMAT (1X,'DPHTRF(',A,') :',t30,f7.3:t40,f8.3,' Mflops')
     !	 write (*,*) 'INFO = ', INFO
     !
  END IF
  RETURN
  !	
  !	Result stored using LAPACK packed format
  !     End of ZHHTRD
  !
END SUBROUTINE ZHHTRD

SUBROUTINE ZLATD4( UPLO, N, NB, ABP, E, TAU, W, LDW, HB )
  IMPLICIT NONE
  !
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER, intent(in) ::  UPLO
  INTEGER,   intent(in) ::  LDW, N, NB, HB
  !     ..
  !     .. Array Arguments ..
  DOUBLE PRECISION, intent(out) :: E( * )
  COMPLEX*16,      intent(inout):: ABP( * )
  COMPLEX*16,      intent(out)  :: TAU( * ), W( LDW, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  ZLATD4 reduces NB rows and columns of a real symmetric matrix A to
  !  symmetric tridiagonal form by an orthogonal similarity
  !  transformation Q' * A * Q, and returns the matrices V and W which are
  !  needed to apply the transformation to the unreduced part of A.
  !
  !  If UPLO = 'U', ZLATD4 reduces the last NB rows and columns of a
  !  matrix, of which the upper triangle is supplied;
  !
  !  This is an auxiliary routine called by DSYTRD.
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER
  !          Specifies whether the upper or lower triangular part of the
  !          symmetric matrix A is stored:
  !          = 'U': Upper triangular
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.
  !
  !  NB      (input) INTEGER
  !          The number of rows and columns to be reduced.
  !
  !  ABP     (input/output) COMPLEX*16 array, 
  !          dimension (1:n*(n+1)/2+(n/hb+1)*hb*(hb-1)/2)
  !          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !          n-by-n upper triangular part of A contains the upper
  !          triangular part of the matrix A, and the strictly lower
  !          triangular part of A is not referenced.  
  !          On exit:
  !          if UPLO = 'U', the last NB columns have been reduced to
  !            tridiagonal form, with the diagonal elements overwriting
  !            the diagonal elements of A; the elements above the diagonal
  !            with the array TAU, represent the orthogonal matrix Q as a
  !            product of elementary reflectors.
  !          See Further Details.
  !
  !  E       (output) DOUBLE PRECISION array, dimension (N-1)
  !          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
  !          elements of the last NB columns of the reduced matrix;
  !
  !  TAU     (output) COMPLEX*16 array, dimension (N-1)
  !          The scalar factors of the elementary reflectors, stored in
  !          TAU(n-nb:n-1) if UPLO = 'U'.
  !          See Further Details.
  !
  !  W       (output) COMPLEX*16 array, dimension (LDW,NB)
  !          The n-by-nb matrix W required to update the unreduced part
  !          of A.
  !
  !  HB      (input) INTEGER
  !          The blocksize used in packed storage scheme.
  !
  !  LDW     (input) INTEGER
  !          The leading dimension of the array W. LDW >= max(1,N).
  !
  !  Further Details
  !  ===============
  !
  !  If UPLO = 'U', the matrix Q is represented as a product of elementary
  !  reflectors
  !
  !     Q = H(n) H(n-1) . . . H(n-nb+1).
  !
  !  Each H(i) has the form
  !
  !     H(i) = I - tau * v * v'
  !
  !  where tau is a real scalar, and v is a real vector with
  !  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
  !  and tau in TAU(i-1).
  !
  !  The elements of the vectors v together form the n-by-nb matrix V
  !  which is needed, with W, to apply the transformation to the unreduced
  !  part of the matrix, using a symmetric rank-2k update of the form:
  !  A := A - V*W' - W*V'.
  !
  !  The contents of A on exit are illustrated by the following example
  !  with n = 5 and nb = 2:
  !
  !    (  a   a   a   v4  v5 )
  !    (      a   a   v4  v5 )  
  !    (          a   1   v5 ) 
  !    (              d   1  )
  !    (                  d  )
  !
  !  where d denotes a diagonal element of the reduced matrix, a denotes
  !  an element of the original matrix that is unchanged, and vi denotes
  !  an element of the vector defining H(i).
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  COMPLEX*16   ZERO, ONE, HALF
  PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ), HALF = ( 0.5D+0, 0.0D+0 ) )
  !     ..
  !     .. Local Scalars ..
  INTEGER            I, IW, LDI, LDJ
  COMPLEX*16         ALPHA
  INTEGER            I1I, I1I1, I1J, III, IJI, IJJ, II1I, III1, J, JB 
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           ZAXPY, ZGEMV, ZLARFG, ZSCAL, ZHEMV
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  !COMPLEX*16         ZDOTC
  EXTERNAL           LSAME!, ZDOTC
  INTEGER            IJ2K
  EXTERNAL		 IJ2K
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MIN
  
  !     ..
  !     .. Executable Statements ..
  !
  !     Quick return if possible
  !
  IF( N.LE.0 ) RETURN
  !
  IF( LSAME( UPLO, 'U' ) ) THEN
     !
     !        Reduce last NB columns of upper triangle
     !
     !	 write(*,*) 'ZLATD4: N=',N,' NB=',NB
     LDI = N
     DO I = N, N - NB + 1, -1   ! 10
        IW = I - N + NB
        I1I = IJ2K( UPLO, 1, I, N, HB )
        II1I = IJ2K( UPLO, I-1, I, N, HB )
        III = IJ2K(UPLO, I, I, N, HB )
        IF( I.LT.N ) THEN
           I1I1 = IJ2K( UPLO, 1, I+1, N, HB )
           III1 = IJ2K( UPLO, I, I+1, N, HB )
           !
           !              Update A(1:i,i)
           !
           ABP( III ) = DBLE( ABP( III ) )
           CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
           CALL ZGEMV( 'No transpose', I, N-I, -ONE, ABP(I1I1), LDI, W( I, IW+1 ), LDW, ONE, ABP(I1I), 1 )
           CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
           CALL ZLACGV( N-I, ABP( III1 ), LDI )
           CALL ZGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ), LDW, ABP(III1), LDI, ONE, ABP(I1I), 1 )
           CALL ZLACGV( N-I, ABP( III1 ), LDI )
           ABP( III ) = DBLE( ABP( III ) )
        END IF
        IF( I.GT.1 ) THEN
           !
           !              Generate elementary reflector H(i) to annihilate
           !              A(1:i-2,i)
           !
           ALPHA = ABP( II1I )
           !	       print *, 'Alpha=',alpha
           CALL ZLARFG( I-1, ALPHA, ABP(I1I), 1, TAU( I-1 ) )
           E( I-1 ) = ALPHA
           ABP(II1I) = ONE
           !
           !              Compute W(1:i-1,i)
           !
           
           DO J = 1, I-1
              W(J, IW) = ZERO
           ENDDO
           DO J = 1, I-1, HB
              JB = MIN( HB, I-J)
              IF( JB.EQ.HB ) THEN
                 LDJ = J+JB-1
              ELSE
                 LDJ = LDI
              END IF
              I1J = IJ2K( UPLO, 1, J, N, HB )
              IJI = IJ2K( UPLO, J, I, N, HB ) 
              IJJ = IJ2K( UPLO, J, J, N, HB ) 
              IF( J.NE.1) THEN 
                 CALL ZGEMV( 'NoTranspose', J-1, JB, ONE, ABP(I1J), LDJ, ABP(IJI), 1, ONE, W(1, IW), 1) 
                 CALL ZGEMV( 'Conjugate Transpose', J-1, JB, ONE, ABP(I1J), LDJ, ABP(I1I), 1, ONE, W(J, IW), 1)
              END IF
              !	         write(*,*) 'HB=',HB,' J=',J,' JB=',JB,' I=',I
              CALL ZHEMV( 'U', JB, ONE, ABP(IJJ), LDJ, ABP(IJI), 1, ONE, W( J, IW), 1)
           ENDDO
           
           
           IF( I.LT.N ) THEN
              CALL ZGEMV( 'Conjugate Transpose', I-1, N-I, ONE, W( 1, IW+1 ), LDW, ABP(I1I), 1, ZERO, W( I+1, IW ), 1 )
              CALL ZGEMV( 'No transpose', I-1, N-I, -ONE, ABP(I1I1), LDI, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 )
              CALL ZGEMV( 'Conjugate Transpose', I-1, N-I, ONE, ABP(I1I1), LDI, ABP(I1I), 1, ZERO, W( I+1, IW ), 1 )
              CALL ZGEMV( 'No transpose', I-1, N-I, -ONE, W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, W( 1, IW ), 1 )
           END IF
           
           CALL ZSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
           ALPHA  = -HALF*TAU(i-1)*dot_product(W(1:i-1,iw), ABP(i1i:i1i+i-2))
           CALL ZAXPY( I-1, ALPHA, ABP(I1I), 1, W( 1, IW ), 1 )
        END IF
        !
     ENDDO !10      CONTINUE
  END IF
  !
  RETURN
  !
  !     End of ZLATD4
  !
END SUBROUTINE ZLATD4

SUBROUTINE ZP2BP( UPLO, N, AP, HB )
  IMPLICIT NONE
  !
  !  ZP2BP transforms a packed stored matrix AP into a block packed matrix.
  !
  !
  !     .. Scalar Arguments ..
  CHARACTER*1        UPLO
  INTEGER            N, HB
  COMPLEX*16         AP(*)
  !     
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix AP is stored as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of AP
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of AP
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A. 
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - COMPLEX*16
  !           On entry, AP stores the matrix in packed storage format.
  !           In order to store the (not referenced elements in the 
  !           block diagonal), AP must be at least of size 
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !           On exit, AP stores the matrix in block packed storage.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            IJ2K
  EXTERNAL           LSAME, IJ2K
  !     .. External Subroutines ..
  EXTERNAL           XERBLA
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            INFO, I, I1, I2, I3, J
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  UPPER = LSAME( UPLO, 'U' )
  !
  INFO = 0
  IF(( .NOT.UPPER).AND.( .NOT.LSAME( UPLO , 'L' ) ) )THEN
     INFO = 1
  ELSE IF( N  .LT.0 )THEN
     INFO = 2
  ELSE IF( HB .LT.0 ) THEN
     INFO = 4
  END IF
  !
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'ZP2BP ', INFO )
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF ( UPPER ) THEN
     I1 = N*(N+1)/2
     I2 = IJ2K('U',N,N,N,HB)
     I3 = 0
     DO J = N,1,-1      ! 20
        DO I = J,1,-1   ! 10
           AP(I2) = AP(I1)
           I1 = I1 - 1
           I2 = I2 - 1
        ENDDO !10         CONTINUE
        IF ((J-1) .EQ. ((N/HB)*HB)) I3 = HB-1
        I3 = MOD (I3+1, HB)
        DO I = 0, I3-1  ! 15
           AP(I2-I) = 0.0
        ENDDO ! 15       CONTINUE
        I2 = I2 - I3
     ENDDO !  20    CONTINUE
  ELSE
     I1 = N*(N+1)/2
     I2 = IJ2K('L',N,N,N,HB)
     I3 = MOD(N-1,HB) 
     DO J = N,1,-1      ! 40
        DO I = N,J,-1   ! 30
           AP(I2) = AP(I1)
           I1 = I1 -1
           I2 = I2 -1
        ENDDO !30         CONTINUE
        DO I=0,I3-1  ! 35
           AP(I2-I)=0.0
        ENDDO !35            CONTINUE
        I2 = I2 - I3
        I3 = MOD(I3 - 1+HB, HB)
     ENDDO !   40    CONTINUE
  END IF
  !
  !     End of ZP2BP
  !
END SUBROUTINE ZP2BP

SUBROUTINE ZPHTR4( N, ABP, HB, INFO )
  IMPLICIT NONE
  !
  !
  !
  !     .. Scalar Arguments ..
  INTEGER            INFO, N, HB
  !     ..
  !     .. Array Arguments ..
  COMPLEX*16         ABP( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  ZPHTR4 computes an U**T*U factorization of a real symmetric
  !  positive definite matrix A.
  !
  !  This is the hierarchically blocked right-looking (k-variant) of
  !  the algorithm.
  !
  !  Arguments
  !  =========
  !
  !  N       (input) INTEGER
  !          The number of rows and columns of the matrix A.  N >= 0.
  !
  !  ABP     (input/output) COMPLEX*16 array, dimension (N*N/2 + HB*???)
  !          Only the block diagonal and the upper triangular part of A
  !          are stored columnwise (block packed storage). All blocks
  !          of the block diagonal are of order HB (the order of the last
  !          one may be less).
  !          On entry, the hermitian matrix A.  Only the upper triangular
  !          part of A is referenced.
  !          On exit, the factor U from the Cholesky factorization.
  !          The matrix U is stored using block packed storage.
  !
  !  HB      (input) INTEGER
  !          The blocksize for the hierarchically blocked algorithm.
  !          HB determines the order of the diagonal blocks. HB > 0.
  !          = 1  : use packed storage format
  !          >= N : use conventional storage format
  !          else : use block packed storage
  !
  !  INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order K is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  COMPLEX*16         ONE
  PARAMETER          ( ONE = (1.0D+0, 0) )
  !     ..
  !     .. Local Scalars ..
  INTEGER            K, KK, KKK, KB, KKB, KKKB
  INTEGER            KP, KPEND
  INTEGER            KTMP1, KTMP2, KTMP3
  !     ..
  !     .. External Subroutines and Functions ..
  INTEGER            IJ2K
  EXTERNAL           ZTRSM, IJ2K
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MIN
  !     ..
  !     .. Executable Statements ..
  !
  INFO = 0
  IF( HB.EQ.1 ) THEN
     CALL ZPPTRF( 'U', N, ABP, INFO )
  ELSEIF ( HB.GE.N ) THEN
     CALL ZPOTRF( 'U', N, ABP, N, INFO )
  ELSE
     DO K = 1, N, HB   ! 10
        KB = MIN( HB, N-K+1 )
        KTMP1 = IJ2K( 'U', K, K, N, HB )
        CALL ZPOTRF( 'U', KB, ABP( KTMP1 ), K+KB-1, INFO )
        IF( K+KB.LE.N ) THEN
           KPEND = ((N-K-KB+1)+(HB-1))/HB -1
           DO KP = 0, KPEND   ! 7
              KK = K+KB + KP*HB
              KKB = MIN( HB, N-KK+1 )
              KTMP2 = IJ2K( 'U', K, KK, N, HB )
              CALL ZTRSM( 'Left', 'U', 'ConjugateTranspose','Non-unit',KB, KKB, ONE, ABP(KTMP1), K+KB-1,ABP( KTMP2 ), KK+KKB-1 )
           ENDDO !    7          CONTINUE
           ! DOACROSS Local (KP,KK,KKB,KKK,KKKB,KTMP1,KTMP2,KTMP3)
           DO KP = 0, KPEND    ! 5
              KK = K+KB + KP*HB
              KKB = MIN( HB, N-KK+1 )
              KTMP1 = IJ2K( 'U', K, KK, N, HB )
              KTMP2 = IJ2K( 'U', KK, KK, N, HB )
              CALL ZHERK( 'U', 'ConjugateTranspose', KKB,KB, -ONE, ABP( KTMP1 ), KK+KKB-1, ONE,ABP( KTMP2 ), KK+KKB-1  )
              IF( KK+KKB.LE.N ) THEN
                 DO KKK = KK+KKB, N, HB    ! 3
                    KKKB = MIN( HB, N-KKK+1 )
                    KTMP1 = IJ2K( 'U', K, KKK, N, HB )
                    KTMP2 = IJ2K( 'U', K, KK, N, HB )
                    KTMP3 = IJ2K( 'U', KK, KKK, N, HB )
                    CALL ZGEMM( 'ConjugateTranspose','No transpose',KKB, KKKB, KB, -ONE, ABP( KTMP2 ),KK+KKB-1,ABP ( KTMP1 ), KKK+KKKB-1, ONE, ABP( KTMP3 ),KKK+KKKB-1 )
                 ENDDO !3                   CONTINUE
              END IF
           ENDDO !    5          CONTINUE
        END IF
        IF( INFO.NE.0 ) GO TO 20
     ENDDO ! 10    CONTINUE
  END IF
  RETURN
20 CONTINUE
  write(6,*) 'K in zphtr4',k,info
  INFO = INFO + K - 1
  RETURN
  !
  !     End of ZPHTR4
  !
END SUBROUTINE ZPHTR4

SUBROUTINE ZPHTRF( UPLO, N, AP, HB, INFO ) 
  !
  !
  CHARACTER      UPLO
  INTEGER        INFO, N, HB
  COMPLEX*16     AP( * )
  !       Hierarchically blocked version of ZPPTRF
  !
  !
  !	Parameters
  !
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower triangular
  !           part of the matrix AP is stored as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of AP
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of AP
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the symmetric matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - COMPLEX*16
  !           On entry, AP stores the matrix in packed storage format.
  !           In order to store the (not referenced elements in the
  !           block diagonal), AP must be at least of size
  !           N*(N+1)/2 + (N/HB+1)*HB(HB-1)/2.
  !           On exit, AP stores the matrix in block packed storage.
  !
  !  HB     - INTEGER.
  !           On entry, HB specifies the block size of the diagonal blocks.
  !           HB must be at least zero.
  !           Unchanged on exit.
  !
  !  INFO  - (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order K is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !  =====================================================================
  !           
  Integer            j
  REAL*8             CP(4), FLOP
  logical UPPER, LSAME
  integer ilaenv
  external ilaenv, LSAME
  !
  !     Test the input parameters.
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( HB.LT.0 ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZPHTRF', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) RETURN
  
  IF(HB.le.1) THEN
     call ZPPTRF( UPLO, N, AP, INFO ) 
  ELSE
     call CPUTIM(CP(1))
     !
     !      call the subroutines
     !			
     !
     !      convert from LAPACK upper triangular packed storage format
     !                     into a upper triangular block packed matrix
     !
     call zp2bp('U',N,AP,HB)
     call CPUTIM(CP(2))
     !
     !      compute an U**T*U factorization
     !      using the hierarchically blocked k-version of the UTU Algorithm
     !
     call ZPHTR4(N,AP,HB,INFO)
     call CPUTIM(CP(3))
     IF (INFO .NE. 0) THEN
        write(6,*) 'Info in zphtrf',info
        CALL OUTERR('zphtrf','ZPHTR4 aborted unsuccessfully.')
        GOTO 999
     ENDIF
     !
     !      convert back to upper triangular packed storage format
     !
     call zbp2p('U',N,AP,HB)
     call CPUTIM(CP(4))
     !
     !      timing output
     !
     !         DO 30 j = 1, 3
     !           CP(j) = CP(j+1)-CP(j)
     ! 30      CONTINUE
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Packing' , CP(1) , 1E-6*Flop/CP(1)
     !         Flop = (N*DBLE(N+1)*(2*N+1))/6
     !         WRITE (*,1001) 'Cholesky' , CP(2) , 1E-6*Flop/CP(2)
     !         Flop = (N*DBLE(N+1))/2
     !         WRITE (*,1001) 'Unpacking' , CP(3) , 1E-6*Flop/CP(3)
     ! 1001    FORMAT (1X,'ZPPTRFB(',A,') :',t30,f7.3:t40,f8.3,' Mflops')
     !	 write (*,*) 'INFO = ', INFO
  END IF
  RETURN
  !	
  !	Result stored using LAPACK packed format
  !
  !
999 STOP 'LAPW1 - Error'
  !
END SUBROUTINE ZPHTRF
