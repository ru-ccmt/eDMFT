subroutine sort_to_diagonal(A, ind, n)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(n,n)
  INTEGER, intent(out)    :: ind(n)
  INTEGER, intent(in)     :: n
  ! locals
  COMPLEX*16 :: Anew(n,n)
  REAL*8     :: column(n), m
  INTEGER    :: wind(n,n), cind(n), unique(n), towork(n), degenerate(n), present_(n), missing(n), tocorrect(n)
  INTEGER    :: i, j, current, t_size, d_size, m_size, ti, tii, ic, c
  LOGICAL    :: ind_in_deg
  wind=0
  do i=1,N
     cind = (/(j,j=1,n)/)
     column = abs(A(:,i))
     call SSORT(column, cind, n)
     wind(i,:) = cind
  enddo
  ind = wind(:,1)
  !print *, 'ind=', ind-1

  unique=(/(1,j=1,n)/)
  
  do i=1,n
     current = ind(i)
     do j=i+1,n
        if (ind(j) .eq. current) then
           unique(i)=0
           unique(j)=0
        endif
     enddo
  enddo
  !print *, 'unique=', unique

  t_size=0
  do i=1,n
     if (unique(i).eq.0) then
        t_size = t_size+1
        towork(t_size) = i
     endif
  enddo

  d_size=0
  do i=1,n
     if (unique(i).eq.0) then
        ind_in_deg=.false.
        do j=1,d_size
           if (ind(i).eq.degenerate(j)) ind_in_deg=.true.
        enddo
        if (.not. ind_in_deg) then
           d_size = d_size + 1
           degenerate(d_size) = ind(i)
        endif
     endif
  enddo
  !print *, 'degenerate:', degenerate(:d_size)
  
  present_=0
  do i=1,n
     do j=1,n
        if (ind(j).eq.i) present_(i)=1
     enddo
  enddo
  !print *, 'present=', present_
  
  m_size=0
  do i=1,n
     if (present_(i).eq.0) then
        m_size = m_size+1
        missing(m_size) = i
     endif
  enddo
  !print *, 'missing are:', missing(:m_size)
  
  tocorrect(:m_size) = missing(:m_size)
  tocorrect(m_size+1:m_size+d_size) = degenerate(:d_size)
  
  !print *, 'correct=', tocorrect(:m_size+d_size)
  !print *, 'towork=', towork(:t_size)
  
  do ic=1,m_size+d_size
     c = tocorrect(ic)
     ti = towork(1)
     tii=1
     m = abs(A(c,ti))
     do i=2,t_size
        if (abs(A(c,towork(i))) .gt. m) then
           ti = towork(i)
           tii=i
           m = abs(A(c,ti))
        endif
     enddo
     ind(ti)=c

     ! towork.remove(ti)
     do i=tii,t_size-1
        towork(i) = towork(i+1)
     enddo
     t_size = t_size-1
  enddo
  !print *, 'final ind=', ind-1
  
  !do i=1,n
  !   Anew(:,ind(i)) = A(:,i)
  !enddo
  !A(:,:) = Anew(:,:)
  
end subroutine sort_to_diagonal



SUBROUTINE SSORT(X, IY, N)
  IMPLICIT NONE
  !    Example of a Selection Sort   Using a Fortran 90 Intrinsic Function
  !***BEGIN PROLOGUE  SSORT
  !***PURPOSE  Sort an array and make the same interchanges in
  !            an auxiliary array.  The array is sorted in
  !            decreasing order.
  !***TYPE      SINGLE PRECISION
  !***KEYWORDS  SORT, SORTING
  !
  !   Description of Parameters
  !      X - array of values to be sorted   (usually abscissas)
  !      IY - array to be carried with X (all swaps of X elements are
  !          matched in IY .  After the sort IY(J) contains the original
  !          postition of the value X(J) in the unsorted X array.
  !      N - number of values in array X to be sorted
  INTEGER, intent(in) :: N
  !     .. Array Arguments ..  -----NOTE the 2 new ways of declaring array size
  REAL*8, intent(inout)  :: X(N)
  INTEGER, intent(inout) :: IY(N)
  !     .. Local Scalars ..
  REAL TEMP
  INTEGER I, ISWAP(1), ITEMP, ISWAP1
  !     .. Intrinsic Functions ..
  INTRINSIC MAXLOC
  !    MAXLOC is a FORTRAN 90 function that returns the index value for the
  !    maximum element in the array
  !***FIRST EXECUTABLE STATEMENT  SSORT
  !
  DO I=1,N-1
     ISWAP=MAXLOC(X(I:N))
     ISWAP1=ISWAP(1)+I-1
     IF(ISWAP1.NE.I) THEN
        TEMP=X(I)
        X(I)=X(ISWAP1)
        X(ISWAP1)=TEMP
        ITEMP=IY(I)
        IY(I)=IY(ISWAP1)
        IY(ISWAP1)=ITEMP
     ENDIF
  END DO
  RETURN
END SUBROUTINE SSORT
