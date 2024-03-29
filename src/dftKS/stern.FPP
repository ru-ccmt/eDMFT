SUBROUTINE STERN(NST,IORD,IMAT,G,TAU,STG,TAUP)
  IMPLICIT NONE
  !INCLUDE 'param.inc'
  !        Arguments
  INTEGER            IORD, NST
  INTEGER            G(3), STG(3,IORD), IMAT(3,3,IORD)
  DOUBLE PRECISION   TAU(3,IORD)
#ifdef _REAL  
  DOUBLE PRECISION TAUP(IORD)
#else  
  DOUBLE COMPLEX TAUP(IORD)
#endif
  !..................................................................
  !   Stern generates the star of rec lattice vector G(I).
  !   The star vectors are stored in STG, the star-size in NST,
  !   IMAT contains the symmetry-matrices.
  !..................................................................
  !   Parameters
  COMPLEX*16         IMAG
  PARAMETER          (IMAG = (0.0D+0,1.0D+0))
  !        Locals
  INTEGER            I, M
  DOUBLE PRECISION  ::  Pi, TK, Greal(3)
  INTEGER :: ind(iord)
  LOGICAL :: FoundStarMember
  !        Intrinsic Functions
  INTRINSIC          ATAN, DCMPLX, EXP
  !
  pi=acos(-1.d0)
  Greal(:) = G(:)
  NST = 0
  !         start loop over all symmetry operations
  DO i=1,iord
     TK = dot_product(Greal(:),TAU(:,I))*2.0d0*pi
     STG(:,I) = matmul(imat(:,:,I),G(:))
     FoundStarMember=.False.
     DO M = 1, NST !        proof, if the vector STG(J,I) is a new starmember or not
        IF( STG(1,M).NE.STG(1,I) .OR. STG(2,M).NE.STG(2,I) .OR. STG(3,M).NE.STG(3,I) ) CYCLE
        ! if we come here, previous stg(:,m) is exactly equal to current stg(:,i). Hence stg(:,i) is part of the same star.
        ind(M) = ind(M) + 1
#ifdef _REAL           
        TAUP(M) = TAUP(M) + COS(TK)
#else           
        TAUP(M) = TAUP(M) + DCMPLX(COS(TK),SIN(TK))
#endif
        FoundStarMember=.True.
     ENDDO

     if (.not.FoundStarMember) then
        ! We come here only when we start a new star member
        NST = NST + 1
        STG(:,NST) = STG(:,I)
        ind(NST) = 1
#ifdef _REAL        
        TAUP(NST) = COS(TK)
#else        
        TAUP(NST) = DCMPLX(COS(TK),SIN(TK))
#endif
     endif
  ENDDO
  TAUP(1:NST) = TAUP(1:NST)/ind(1:NST)
  RETURN
END SUBROUTINE
