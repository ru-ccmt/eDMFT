SUBROUTINE STERN(G,NST,STG,TAUP, iz, tau, iord, Qcomplex)
  ! Generates the start of reciprocal vector G
  ! Input:
  !     G(3)           -- reciprocal vector
  !     iz(3,3,iord)   -- 3x3 matrices for group operations (rotations, reflections...)
  !     tau(3,iord)    -- vector shifts for non-simmorphic symmetries
  !     Qcomplex       -- calculation is real when inversion is present and there is no spin-orbit. Otherwise we work with complex phase factors.
  ! Output:
  !     nst            -- number of star members
  !     stg            -- the first star member, i.e., iz*G
  !     taup           -- sum of phase factors for non-simmorphic symmetries : sum_{star_members=j} e^{2*pi*i*G*tau_j}, otherwise just degeneracy of star member
  !
  IMPLICIT NONE
  !
  INTEGER, intent(out)      :: NST,STG(3,iord)  
  COMPLEX*16, intent(out)   :: TAUP(iord)
  INTEGER, intent(in)       :: G(3)             
  INTEGER, intent(in)       :: iord
  INTEGER, intent(in)       :: iz(3,3,iord)
  REAL*8, intent(in)        :: tau(3,iord)
  LOGICAL, intent(in)       :: Qcomplex
  ! locals
  INTEGER :: i, M
  REAL*8  :: TK, pi, Greal(3)
  INTEGER :: ind(iord)
  LOGICAL :: FoundStarMember
  !---------------------------------------------------------------------  
  !
  pi=acos(-1.d0)
  NST=0
  Greal(:) = G(:)
  do i=1,iord
     TK = dot_product(Greal(:),TAU(:,I))*2.0d0*pi
     STG(:,I) = matmul(IZ(:,:,I),G(:))
     FoundStarMember=.False.
     DO M=1,NST
        IF( STG(1,M).NE.STG(1,I) .OR. STG(2,M).NE.STG(2,I) .OR. STG(3,M).NE.STG(3,I) ) CYCLE
        ! if we come here, previous stg(:,m) is exactly equal to current stg(:,i). Hence stg(:,i) is part of the same star.
        ind(M)=ind(M)+1
        if (Qcomplex) then
           TAUP(M) = TAUP(M) + dcmplx(cos(TK),sin(TK))
        else
           TAUP(M) = TAUP(M) + cos(TK)
        endif
        FoundStarMember=.True.
     ENDDO

     if (.not.FoundStarMember) then
        ! We come here only when we start a new star member
        NST=NST+1            ! next start member started
        STG(:,NST)=STG(:,I)  ! Symmetry related G
        ind(NST)=1           ! How many members of the star we have
        if (Qcomplex) then
           TAUP(NST) = dcmplx(cos(TK),sin(TK))
        else
           TAUP(NST) = cos(TK)
        endif
     endif
  enddo

  TAUP(1:NST)=TAUP(1:NST)/IND(1:NST)                                            
  RETURN                                                            
END SUBROUTINE STERN
