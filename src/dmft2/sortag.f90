!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
! 
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
! 
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
! 
!
      SUBROUTINE SORTAG (A, N, TAG)
      implicit real*4 (a-h,o-z)
!///
!       SUBROUTINE SORTAG
!
!        DECK USED - CDC 6500, IBM360/67
!
!        SORTAG
!
!        PURPOSE
!
!           TO SORT A VECTOR INTO INCREASING ORDER FROM A(1) TO A(N),
!           A MAY BE TYPE REAL OR TYPE INTEGER.  VECTOR TAG IS PER-
!           MUTED THE SAME AS VECTOR A.
!
!        USAGE
!
!           CALL SORTAG(A,N,TAG)
!
!        DESCRIPTION OF PARAMETERS
!
!           A - THE NAME OF THE N-VECTOR TO BE SORTED.
!               IF A IS TYPE REAL THEN EACH OF ITS COMPONENTS MUST BE
!               NORMALIZED FORM.
!           N - THE NUMBER OF ELEMENTS IN THE VECTOR TO BE SORTED.
!           TAG - THE NAME OF THE N-VECTOR CONTAINING THE TAG FIELDS.
!
!        REMARKS
!
!           THE PROCEDURE REQUIRES TWO ADDITIONAL ARRAYS IU(K) AND
!           IL(K) WHICH PERMIT SORTING UP TO 2**(K+1)-1 ELEMENTS.
!
!        EXAMPLE
!
!           IF N RANDOM NUMBERS, UPON GENERATION, WERE ASSIGNED
!           CONSECUTIVE TAGS OF 1 THROUGH N, AND THE ARRAY OF RANDOM
!           NUMBERS SORTED, THE TAG ARRAY COULD THEN BE USED TO
!           DETERMINE THAT THE RANDOM NUMBER NOW IN POSITION A(I) WAS
!           ORIGINALLY IN POSITION A(J) BECAUSE TAG(I) = J.
!
!        METHOD
!
!           'AN EFFICIENT ALGORITHM FOR SORTING WITH MINIMAL STORAGE'
!           BY RICHARD C. SINGLETON.  PREPARED FOR INSTITUTE RESEARCH
!           AND DEVELOPMENT.  STANFORD RESEARCH INSTITUTE PROJECT
!           387531-132 SEPTEMBER 1968.
!
!     ..................................................................
!      
!      RECEIVED FROM BAUKE DIJKSTR - DEC 1985 - VAX VERSION.
!
!=========================================================================
!      VERSION 23 DECEMBER 1985
!=========================================================================
!\\\

      PARAMETER (IKEY=28)
      DIMENSION A(N),IU(IKEY),IL(IKEY),TAG(N)
      INTEGER TAG,TG
!      INTEGER A,T,TT,TAG,TG

      IF(N .GE. 2**(IKEY+1))THEN
        WRITE(*,*)' TOO MANY ITEMS TO BE SORTED'
        WRITE(*,*)' INCREASE PARAMETER IKEY IN CFBLIB-ROUTINE SORTAG'
        STOP
      ENDIF

      M=1
      I=1
      J=N
    5 IF(I .GE. J) GO TO 70
   10 K=I
      IJ=(J+I)/2
      T=A(IJ)
      IF(A(I) .LE. T) GO TO 20
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(I)
      TAG(I)=TG
   20 L=J
      IF(A(J) .GE. T) GO TO 40
      A(IJ)=A(J)
      A(J)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(J)
      TAG(J)=TG
      IF(A(I) .LE. T) GO TO 40
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(I)
      TAG(I)=TG
      GO TO 40
   30 A(L)=A(K)
      A(K)=TT
      TG=TAG(L)
      TAG(L)=TAG(K)
      TAG(K)=TG
   40 L=L-1
      IF(A(L) .GT. T) GO TO 40
      TT=A(L)
   50 K=K+1
      IF(A(K) .LT. T) GO TO 50
      IF(K .LE. L) GO TO 30
      IF(L-I .LE. J-K) GO TO 60
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 80
   60 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 80
   70 M=M-1
      IF(M.EQ.0) RETURN
      I=IL(M)
      J=IU(M)
80    IF(J-I.GE.1) GOTO 10
      IF(I.EQ.1) GOTO 5
      I=I-1
   90 I=I+1
      IF(I .EQ. J) GO TO 70
      T=A(I+1)
      IF(A(I) .LE. T) GO TO 90
      TG=TAG(I+1)
      K=I
  100 A(K+1)=A(K)
      TAG(K+1)=TAG(K)
      K=K-1
      IF(T .LT. A(K)) GO TO 100
      A(K+1)=T
      TAG(K+1)=TG
      GO TO 90

      END

