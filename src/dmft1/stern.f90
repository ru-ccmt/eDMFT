SUBROUTINE STERN(G,NST,STG,gind, iz, tau, iord)
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
  INTEGER, intent(out)      :: NST
  REAL*8, intent(out)       :: STG(3,iord)  
  INTEGER, intent(out)      :: gind(iord)
  REAL*8,  intent(in)       :: G(3)             
  INTEGER, intent(in)       :: iord
  INTEGER, intent(in)       :: iz(3,3,iord)
  REAL*8, intent(in)        :: tau(3,iord)
  ! locals
  INTEGER :: i, M
  REAL*8  :: tauc(3,iord), pi, Greal(3)
  INTEGER :: ind(iord)
  LOGICAL :: FoundStarMember
  !---------------------------------------------------------------------  
  !
  pi=acos(-1.d0)
  NST=0
  Greal(:) = G(:)
  gind(:)=0
  ind(:)=0
  do i=1,iord
     tauc(:,I) = TAU(:,I)
     STG(:,I) = matmul(IZ(:,:,I),G(:))
     FoundStarMember=.False.
     DO M=1,NST
        IF( abs(STG(1,M)-STG(1,I)).gt.1d-10 .OR. abs(STG(2,M)-STG(2,I)).gt.1d-10 .OR. abs(STG(3,M)-STG(3,I)).gt.1.d-10 .or. sum(abs(tauc(:,m)-tauc(:,i))).gt.1d-6) CYCLE
        !IF( abs(STG(1,M)-STG(1,I)).gt.1d-10 .OR. abs(STG(2,M)-STG(2,I)).gt.1d-10 .OR. abs(STG(3,M)-STG(3,I)).gt.1.d-10) CYCLE
        ! if we come here, previous stg(:,m) is exactly equal to current stg(:,i). Hence stg(:,i) is part of the same star.
        ind(M)=ind(M)+1
        !gind(M) = i
        FoundStarMember=.True.
     ENDDO

     if (.not.FoundStarMember) then
        ! We come here only when we start a new star member
        NST=NST+1            ! next start member started
        STG(:,NST)=STG(:,I)  ! Symmetry related G
        tauc(:,NST)=tauc(:,i)
        ind(NST)=1           ! How many members of the star we have
        gind(NST) = i
     endif
  enddo
  
  RETURN                                                            
END SUBROUTINE STERN
