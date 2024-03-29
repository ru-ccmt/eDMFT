! @Copyright 2007 Kristjan Haule
! 

SUBROUTINE MakeTanMesh(x, N, tanc, tanw, b0, b1)
  IMPLICIT NONE
  REAL*8, intent(out) :: x(2*N)
  INTEGER, intent(in) :: N
  REAL*8, intent(in)  :: tanc, tanw, b0, b1
  ! locals
  INTEGER :: i
  REAL*8 :: du, b1n, m0
  if (.NOT.(b0<b1)) print *,  "Relation must hold: b0<b1!"
  if (.NOT.(b0<tanw .AND. tanw<b1)) print *, "Relation mesu hold: b0<tanw<b1!"
  if (.NOT.(b0>0))  print *, "b0 must be positive!"
  du = atan(((tanc-b0)/tanw))
  b1n = atan((b1-tanc)/tanw)+du
  !m0 = [tanc + tanw * tan(b1n*(i-1)/(N-2)-du) for i in range(1,N)]
  do i=1,N
     m0 = tanc + tanw * tan(b1n*(i-1)/(N-1)-du)
     x(N+i) = m0
     x(N-(i-1)) = -m0
  end do
END SUBROUTINE MakeTanMesh

SUBROUTINE CombineMesh(x3, ni, om0, om, x, nom, nx)
  IMPLICIT NONE
  REAL*8, intent(out):: x3(nom+nx)
  INTEGER, intent(out):: ni
  REAL*8, intent(in) :: om0, om(nom), x(nx)
  INTEGER, intent(in):: nom, nx
  ! locals
  INTEGER :: i
  x3=0.0
  ni=0
  do i=1,nom
     if (om(i) >= om0-x(nx) .AND. om(i) <= om0-x(1) ) then
        ni = ni+1
        x3(ni) = om(i)
     end if
  end do
  do i=1,nx
     if (x(i)+om0 >= om(1) .AND. x(i)+om0 <= om(nom)) then
        ni = ni+1
        x3(ni) = x(i)+om0
     endif
  enddo
  
  CALL SHELL(ni,x3(1:ni))
  
END SUBROUTINE CombineMesh

!*****************************************************
!* Sorts an array ARR of length N in ascending order *
!*            by the Shell-Mezgar method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table ARR                  *
!*          ARR	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    ARR   table sorted in ascending order    *
!*                                                   *
!* NOTE: The Shell method is a N^3/2 routine and can *
!*       be used for relatively large arrays.        *
!*****************************************************         
SUBROUTINE SHELL(N,ARR)
  IMPLICIT NONE
  real*8, intent(inout) :: ARR(N)
  INTEGER, intent(in)   :: N
  !
  real*8 :: t, ALN2I, TINY
  INTEGER :: LOGNB2, m, nn, k, i, j, l
  parameter(ALN2I=1./0.69314718,TINY=1.E-5)
  !
  LOGNB2=INT(ALOG(FLOAT(N))*ALN2I+TINY)
  m=n
  do nn=1,LOGNB2
     m=m/2; k=n-m
     do j=1,k
        i=j
10      continue
        l=i+m
        if(ARR(l).LT.ARR(i)) then
           t=ARR(i)
           ARR(i)=ARR(l)
           ARR(l)=t
           i=i-m
           if(i.GE.1) GOTO 10
        end if
     end do
  end do
  return
END SUBROUTINE SHELL
