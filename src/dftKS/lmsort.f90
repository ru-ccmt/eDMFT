SUBROUTINE LMSORT(NATO,NGAUNT,JATOM,LQIND,LQNS,GFAC,L0PTR, L0FLD,M0PTR,M0FLD,LPPTR,LPFLD,MPPTR,MPFLD,LMXPTR,LMXFLD,GNTFLD)
  IMPLICIT NONE
  !        Arguments
  INTEGER,    intent(in) :: NATO, NGAUNT, JATOM
  INTEGER,    intent(in) :: LQIND(NATO)
  INTEGER,    intent(in) :: LQNS(6,NGAUNT,NATO)
  COMPLEX*16, intent(in) :: GFAC(NGAUNT,NATO)
  !
  INTEGER,    intent(out):: L0PTR(NGAUNT+2),  L0FLD(NGAUNT)
  INTEGER,    intent(out):: M0PTR(NGAUNT+2),  M0FLD(NGAUNT)
  INTEGER,    intent(out):: LPPTR(NGAUNT+2),  LPFLD(NGAUNT)
  INTEGER,    intent(out):: MPPTR(NGAUNT+2),  MPFLD(NGAUNT)
  INTEGER,    intent(out):: LMXPTR(NGAUNT+2), LMXFLD(NGAUNT)
  COMPLEX*16, intent(out):: GNTFLD(NGAUNT)
  !
  !     ..................................................................
  !
  ! 2.     PURPOSE
  !           Convert the LP LL L0 LMX MP M0 combinations (sorted in that
  !           order) as given by array LQNS for a given atom JATOM into
  !           data structures suitable for the following processing
  !           scheme:
  !
  !           DO L0's 
  !              DO M0's associated with L0
  !                 DO LP's associated with (L0,M0)
  !                    DO MP's associated with (L0,M0,LP)
  !                       DO LMX's associated with (L0,M0,LP,MP)
  !
  !                          ...
  !
  !                       ENDDO
  !                    ENDDO
  !                 ENDDO
  !              ENDDO
  !           ENDDO
  !
  ! 3.     USAGE
  !        ARGUMENT-DESCRIPTION
  !           NATO   - INTEGER value                               (input)
  !                    Maximum number of useable non-equivalent atoms.
  !                    Used for dimensioning of LQIND, LQNS and GFAC.
  !
  !           NGAUNT - INTEGER value                               (input)
  !                    Maximum number of usable gaunt factors.
  !                    Used for dimensioning purposes.
  !
  !           JATOM  - INTEGER value                               (input)
  !                    Index of the actual non-equivalent atom for which
  !                    the restructuring should take place.
  !                    Constraint:  (1 .LE. JATOM .LE. NATO)
  !
  !           LQIND  - INTEGER array, dimension (NATO)             (input)
  !                    LQIND(JATOM) gives the number of l,m,l',m',L,M
  !                    combinations for atom JATOM
  !                    that means the number of gaunt factors for atom i
  !
  !           LQNS   - INTEGER array, dimension (6,NGAUNT,NATO)    (input)
  !                    LQNS(:,j,JATOM) ... l,m,LM,l',m' combination
  !                    corresponding to the gaunt-factor GFAC(j,JATOM)
  !                    for atom JATOM
  !                       LQNS(1,1:LQIND(i),i) ... l'+1 for atom i
  !                       LQNS(2,1:LQIND(i),i) ... L
  !                       LQNS(3,1:LQIND(i),i) ... l+1 for atom i
  !                       LQNS(4,1:LQIND(i),i) ... index of LM for atom i
  !                       LQNS(5,1:LQIND(i),i) ... m' for atom i
  !                       LQNS(6,1:LQIND(i),i) ... m for atom i
  !                    (1 .LE. j .LE. LQIND(i))
  !
  !           GFAC   - COMPLEX*16 array, dimension (NGAUNT,NATO)   (input)
  !                    GFAC(:,JATOM) ... gaunt factors for atom JATOM
  !
  !           L0PTR  - INTEGER array, dimension (NGAUNT+2)        (output)
  !                    points to a range of l values stored in L0FLD
  !                    L0PTR(1) provides the index of the first applicable
  !                       l value stored in L0FLD (L0FLD(L0PTR(1))
  !                    L0PTR(2)-1 gives the index of the last applicable
  !                       l value stored in L0FLD (L0FLD(L0PTR(2)-1))
  !                    L0PTR(3) .EQ. NIL  (NIL = -1 ... null pointer)
  !
  !           L0FLD  - INTEGER array, dimension (NGAUNT)          (output)
  !                    stores the l values sorted in ascending order
  !                    L0FLD(L0PTR(1)) ..... first (lowest) applicable l
  !                    L0FLD(L0PTR(2)-1) ... last (highest) applicable l
  !
  !           M0PTR  - INTEGER array, dimension (NGAUNT+2)        (output)
  !                    points to ranges of m values stored in M0FLD, each
  !                    associated with a specific l value;
  !                    given l = L0FLD(I):
  !                    M0PTR(I) provides the index of the first applicable
  !                       m value associated with l
  !                    M0PTR(I+1)-1 gives the index of the last applicable
  !                       m value associated with l
  !
  !           M0FLD  - INTEGER array, dimension (NGAUNT)          (output)
  !                    stores the m values arranged in ranges, each
  !                    associated with a specific l and sorted in
  !                    ascending order;
  !                    given l = L0FLD(I):
  !                    M0FLD(M0PTR(I)) ....... first (lowest) applicable m
  !                    M0FLD(M0PTR(I+1)-1) ... last (highest) applicable m
  !
  !           LPPTR  - INTEGER array, dimension (NGAUNT+2)        (output)
  !                    points to ranges of l' values stored in LPFLD, each
  !                    associated with a specific l,m combination;
  !                    given m = M0FLD(I):
  !                    LPPTR(I) provides the index of the first applicable
  !                       l' value associated with l,m
  !                    LPPTR(I+1)-1 gives the index of the last applicable
  !                       l' value associated with l,m
  !
  !           LPFLD  - INTEGER array, dimension (NGAUNT)          (output)
  !                    stores the l' values arranged in ranges, each
  !                    associated with a specific l,m and sorted in
  !                    ascending order;
  !                    given m = M0FLD(I):
  !                    LPFLD(LPPTR(I)) ...... first (lowest) applicable l'
  !                    LPFLD(LPPTR(I+1)-1) .. last (highest) applicable l'
  !
  !           MPPTR  - INTEGER array, dimension (NGAUNT+2)        (output)
  !                    points to ranges of m' values stored in MPFLD, each
  !                    associated with a specific l,m,l' combination;
  !                    given l' = LPFLD(I):
  !                    MPPTR(I) provides the index of the first applicable
  !                       m' value associated with l,m,l'
  !                    MPPTR(I+1)-1 gives the index of the last applicable
  !                       m' value associated with l,m,l'
  !
  !           MPFLD  - INTEGER array, dimension (NGAUNT)          (output)
  !                    stores the m' values arranged in ranges, each
  !                    associated with a specific l,m,l' and sorted in
  !                    ascending order;
  !                    given l' = LPFLD(I):
  !                    MPFLD(MPPTR(I)) ...... first (lowest) applicable m'
  !                    MPFLD(MPPTR(I+1)-1) .. last (highest) applicable m'
  !
  !           LMXPTR - INTEGER array, dimension (NGAUNT+2)        (output)
  !                    points to ranges of lmx values stored in LMXFLD,
  !                    each associated with a specific l,m,l',m'
  !                    combination;
  !                    given m' = MPFLD(I):
  !                    LMXPTR(I) provides the index of the first
  !                       applicable lmx value associated with l,m,l',m'
  !                    LMXPTR(I+1)-1 gives the index of the last
  !                       applicable lmx value associated with l,m,l',m'
  !
  !           LMXFLD - INTEGER array, dimension (NGAUNT)          (output)
  !                    stores the lmx values arranged in ranges, each
  !                    associated with a specific l,m,l',m' and sorted in
  !                    ascending order;
  !                    given m' = MPFLD(I):
  !                    LMXFLD(LMXPTR(I))     first (lowest) applicable lmx
  !                    LMXFLD(LMXPTR(I+1)-1) last (highest) applicable lmx
  !
  !           GNTFLD - COMPLEX*16 array, dimension (NGAUNT)       (output)
  !                    stores the gaunt factors associated with a specific
  !                    l,m,l',m',L,M combination (represented by a
  !                    l,m,l',m',lmx combination);
  !                    given lmx = LMXFLD(I):
  !                    GNTFLD(I) ... gaunt factor for l,m,l',m',L,M
  !           
  !        USED SUBROUTINES (DIRECTLY CALLED)
  !           none
  !
  !        INDIRECTLY CALLED SUBROUTINES
  !           none
  !
  !        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
  !           none
  !
  !        INPUT/OUTPUT (READ/WRITE)
  !           none
  !
  !        MACHINENDEPENDENT PROGRAMPARTS
  !           COMPLEX*16 declarations for GNTFLD and GFAC are used
  !
  ! 4.     REMARKS
  !           The following gives a short overview of the generated data
  !           structures:
  !
  !           L0FLD(L0PTR(1)), ... ,L0FLD(L0PTR(2)-1)
  !               l (L0) values (sorted ascending for increasing indices)
  !           M0FLD(M0PTR(J)), ... ,M0FLD(M0PTR(J+1)-1)
  !               m (M0) values for a given l, where
  !               index J = L0PTR(1)
  !           LPFLD(LPPTR(J)), ... ,LPFLD(LPPTR(J+1)-1)
  !               l' (LP) values for a given l,m combination, where
  !               index J = M0PTR(L0PTR(1))
  !           MPFLD(MPPTR(J)), ... ,MPFLD(MPPTR(J+1)-1)
  !               m' (LP) values for a given l,m,l' combination, where
  !               index J = LPPTR(M0PTR(L0PTR(1)))
  !           LMXFLD(LMXPTR(J)), ... ,LMXFLD(LMXPTR(J+1)-1)
  !               indices to lm values (LMX) for a  l,m,l'm' combination,
  !               where index J = MPPTR(LPPTR(M0PTR(L0PTR(1))))
  !           GNTFLD(K)
  !               gaunt factor associated with a l,m,l',m',LM combination
  !               index K is the same as the index used to access LMX from
  !               LMXFLD(K)
  !
  ! 5.     METHOD
  !           1. Sort l,m,l',m',lmx using a primitive insert sort
  !              algorithm, putting the values in the respective *FLD
  !              arrays.
  !              Now L0FLD(I),M0FLD(I),LPFLD(I),MPFLD(I),LMXFLD(I),
  !              GNTFLD(I) for (1 .LE. I .LE. LQIND(JATOM)) gives a
  !              l,m,l',m',lmx combination and the associated gaunt value.
  !              But it would not be easily possible to loop over all
  !              say m' values associated with a l,m,l' combination.
  !              Therefore,
  !           2. For each *FLD determine the appropriate ranges of values
  !              storing these ranges in the respective *PTR array as
  !              indicated in the PARAMETERS section above. During the
  !              range determination process all duplicate values will
  !              be removed.
  !        Parameters
  INTEGER            NIL
  PARAMETER          (NIL = -1)
  !
  !        Local Scalars
  !
  INTEGER            I, IDX, L0, LMX, LP, M0, MP, PTR, PTRBEG
  !
  IF (LQIND(JATOM) .LT. 1) THEN
     L0PTR(1) = 1
     L0PTR(2) = NIL
     GOTO 999
  ENDIF
  M0PTR(1) = 1
  LPPTR(1) = 1
  MPPTR(1) = 1
  LMXPTR(1) = 1
  L0FLD(1) = LQNS(3,1,JATOM) - 1
  M0FLD(1) = LQNS(6,1,JATOM)
  LPFLD(1) = LQNS(1,1,JATOM) - 1
  MPFLD(1) = LQNS(5,1,JATOM)
  LMXFLD(1) = LQNS(4,1,JATOM)
  GNTFLD(1) = GFAC(1,JATOM)
  DO I = 2, LQIND(JATOM)  !        Sort next L0, M0, LP, MP, LM tupel into existing fields
     L0 = LQNS(3,I,JATOM) - 1
     M0 = LQNS(6,I,JATOM)
     LP = LQNS(1,I,JATOM) - 1
     MP = LQNS(5,I,JATOM)
     LMX = LQNS(4,I,JATOM)
     DO IDX = I - 1, 1, -1  ! 10
        IF ((L0 .LT. L0FLD(IDX)) .OR. (L0 .EQ. L0FLD(IDX)).AND.((M0 .LT. M0FLD(IDX)) .OR. (M0 .EQ. M0FLD(IDX)).AND.((LP .LT. LPFLD(IDX)) .OR. (LP .EQ. LPFLD(IDX)).AND.((MP .LT. MPFLD(IDX)) .OR. (MP .EQ. MPFLD(IDX)).AND.(LMX .LT. LMXFLD(IDX)))))) THEN
           L0FLD(IDX+1) = L0FLD(IDX)
           M0FLD(IDX+1) = M0FLD(IDX)
           LPFLD(IDX+1) = LPFLD(IDX)
           MPFLD(IDX+1) = MPFLD(IDX)
           LMXFLD(IDX+1) = LMXFLD(IDX)
           GNTFLD(IDX+1) = GNTFLD(IDX)
        ELSE
           EXIT
        ENDIF
     ENDDO                ! 10
     L0FLD(IDX+1) = L0
     M0FLD(IDX+1) = M0
     LPFLD(IDX+1) = LP
     MPFLD(IDX+1) = MP
     LMXFLD(IDX+1) = LMX
     GNTFLD(IDX+1) = GFAC(I,JATOM)
  ENDDO
  !        determine which L0's are used and remove duplicate L0's
  !        (generate L0PTR)
  PTR = 1
  L0 = L0FLD(1)
  L0PTR(1) = 1
  M0PTR(1) = 1
  DO I = 2, LQIND(JATOM)
     IF (L0FLD(I) .GT. L0) THEN
        PTR = PTR + 1
        L0 = L0FLD(I)
        L0FLD(PTR) = L0
        M0PTR(PTR) = I
     ENDIF
  ENDDO
  L0PTR(2) = PTR + 1
  L0PTR(3) = NIL
  M0PTR(PTR+1) = LQIND(JATOM) + 1
  M0PTR(PTR+2) = NIL
  !    determine which M0's belong to a given L0 and remove
  !    duplicate M0's (generate M0PTR)
  PTR = 1
  DO I = 1, LQIND(JATOM)
     IF (M0PTR(I+1) .EQ. NIL) EXIT
     M0 = M0FLD(M0PTR(I))
     M0FLD(PTR) = M0
     LPPTR(PTR) = M0PTR(I)
     PTRBEG = PTR
     DO IDX = M0PTR(I) + 1, M0PTR(I+1) - 1
        IF (M0FLD(IDX) .GT. M0) THEN
           PTR = PTR + 1
           M0 = M0FLD(IDX)
           M0FLD(PTR) = M0
           LPPTR(PTR) = IDX
        ENDIF
     ENDDO
     PTR = PTR + 1
     M0PTR(I) = PTRBEG
  ENDDO
  M0PTR(I) = PTR
  LPPTR(PTR) = LQIND(JATOM) + 1
  LPPTR(PTR+1) = NIL
  !        determine which LP's belong to a given L0,M0 combination
  !        and remove duplicate LP's (generate LPPTR)
  PTR = 1
  DO I = 1, LQIND(JATOM)
     IF (LPPTR(I+1) .EQ. NIL) EXIT
     LP = LPFLD(LPPTR(I))
     LPFLD(PTR) = LP
     MPPTR(PTR) = LPPTR(I)
     PTRBEG = PTR
     DO IDX = LPPTR(I) + 1, LPPTR(I+1) - 1
        IF (LPFLD(IDX) .GT. LP) THEN
           PTR = PTR + 1
           LP = LPFLD(IDX)
           LPFLD(PTR) = LP
           MPPTR(PTR) = IDX
        ENDIF
     ENDDO
     PTR = PTR + 1
     LPPTR(I) = PTRBEG
  ENDDO
  LPPTR(I) = PTR
  MPPTR(PTR) = LQIND(JATOM) + 1
  MPPTR(PTR+1) = NIL
  !        determine which MP's belong to a given L0,M0,LP combination
  !        and remove duplicate MP's (generate MPPTR)
  PTR = 1
  DO I = 1, LQIND(JATOM)
     IF (MPPTR(I+1) .EQ. NIL) EXIT
     MP = MPFLD(MPPTR(I))
     MPFLD(PTR) = MP
     LMXPTR(PTR) = MPPTR(I)
     PTRBEG = PTR
     DO IDX = MPPTR(I) + 1, MPPTR(I+1) - 1
        IF (MPFLD(IDX) .GT. MP) THEN
           PTR = PTR + 1
           MP = MPFLD(IDX)
           MPFLD(PTR) = MP
           LMXPTR(PTR) = IDX
        ENDIF
     ENDDO
     PTR = PTR + 1
     MPPTR(I) = PTRBEG
  ENDDO
  MPPTR(I) = PTR
  LMXPTR(PTR) = LQIND(JATOM) + 1
  LMXPTR(PTR+1) = NIL
  !        determine which LM indices belong to a given L0,M0,LP,MP
  !        combination; also associate a gaunt factor with a given
  !        L0,M0,LP,MP,LMX combination
  !        (generate LMXPTR and GNTFLD)
  !        Notice: There are no duplicate LMX or gaunt values
  !                because LMX is the last field in the sorted
  !                L0,M0,LP,MP,LMX combinations.
  PTR = 1
  DO I = 1, LQIND(JATOM)
     IF (LMXPTR(I+1) .EQ. NIL) EXIT
     LMX = LMXFLD(LMXPTR(I))
     LMXFLD(PTR) = LMX
     GNTFLD(PTR) = GNTFLD(LMXPTR(I))
     PTRBEG = PTR
     DO IDX = LMXPTR(I) + 1, LMXPTR(I+1) - 1
        IF (LMXFLD(IDX) .GT. LMX) THEN
           PTR = PTR + 1
           LMX = LMXFLD(IDX)
           LMXFLD(PTR) = LMX
           GNTFLD(PTR) = GNTFLD(IDX)
        ENDIF
     ENDDO
     PTR = PTR + 1
     LMXPTR(I) = PTRBEG
  ENDDO
  LMXPTR(I) = PTR
  !
999 RETURN
END SUBROUTINE LMSORT
