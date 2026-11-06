INTEGER Function FindRotatedAtom(iz,tau,jatom,latom)
  !USE sym2,   ONLY: tau, iz
  USE structure, ONLY: pos, mult, lattic
  IMPLICIT NONE
  INTEGER, intent(in) :: jatom, latom
  INTEGER, intent(in) :: iz(3,3)
  REAL*8, intent(in)  :: tau(3)
  !
  REAL*8  :: Rn(3), diff(3), diff2(3)
  INTEGER :: index, m
  REAL*8, PARAMETER :: toler = 1.D-7
  REAL*8 :: toler2
  toler2=1.5d0*toler            

  Rn = matmul( transpose(iz(:,:)),pos(:,latom))
  Rn(:) = Rn(:) + tau(:)
  
  Rn(1) = MOD(Rn(1) + 1.d0, 1.d0)
  Rn(2) = MOD(Rn(2) + 1.d0, 1.d0)
  Rn(3) = MOD(Rn(3) + 1.d0, 1.d0)

  
  index = sum(mult(1:jatom-1))
  DO m=1,mult(jatom)
     index = index + 1
     diff(:) = MOD( abs( Rn(:) - pos(:,index))+toler, 1.d0) - toler
     !print *, 'm=', m, 'index=', index, 'diff=', diff
     if ( maxval(abs(diff)).lt.toler2 ) then
        FindRotatedAtom = index
        return
     endif
     !....check positions for centered lattices
     if(lattic(1:1).eq.'B') then
        diff(:) = mod(diff(:)+0.5d0+toler, 1.d0) - toler
        if ( maxval(abs(diff)).lt.toler2 ) then
           FindRotatedAtom = index
           return
        endif
     endif
     if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
        diff2(:) = diff(:)
        diff2(1) = mod(diff(1)+0.5d0+toler,1.d0)-toler
        diff2(2) = mod(diff(2)+0.5d0+toler,1.d0)-toler
        !print *, 'm=', m, 'index=', index, 'diff2=', diff2
        if ( maxval(abs(diff2)).lt.toler2 ) then
           FindRotatedAtom = index
           return
        endif
     endif
     if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
        diff2(:) = diff(:)
        diff2(1) = mod( diff(1)+0.5d0+toler,1.d0)-toler
        diff2(3) = mod( diff(3)+0.5d0+toler,1.d0)-toler
        !print *, 'm=', m, 'index=', index, 'diff3=', diff2
        if (maxval(abs(diff2)).lt.toler2 ) then
           FindRotatedAtom = index
           return
        endif
     endif
     if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
        diff2(:) = diff(:)
        diff2(2) = mod(diff(2)+0.5d0+toler,1.d0)-toler
        diff2(3) = mod(diff(3)+0.5d0+toler,1.d0)-toler
        !print *, 'm=', m, 'index=', index, 'diff4=', diff2
        if ( maxval(abs(diff2)).lt.toler2 ) then
           FindRotatedAtom = index
           return
        endif
     end if
  ENDDO
  WRITE(6,*) 'ERROR : Could not find equivalent atom in SymmetrizeDM. jatom=', jatom, 'latom=', latom, ' Boiling out'
  FindRotatedAtom = latom
end Function FindRotatedAtom

INTEGER Function FindRotatedCase(latom_rotated)
  use dmfts, ONLY : natom, iatom
  IMPLICIT NONE
  INTEGER, intent(in) :: latom_rotated
  !
  INTEGER :: icase_rotated
  do icase_rotated=1,natom
     if ( iatom(icase_rotated).eq.latom_rotated ) then
        FindRotatedCase = icase_rotated
        return
     endif
  enddo
  WRITE(6,*) 'ERROR : Could not find equivalent icase in SymmetrizeDM. latom_rotated=', latom_rotated
  FindRotatedCase = 1
end Function FindRotatedCase


SUBROUTINE CmpSymmetryTransformation(ig, RT, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2, Qprint)
  USE structure, ONLY: BR1, rotij
  USE sym2,  ONLY: tau, iz
  USE dmfts, ONLY: iso, natom, nl, ll, cix, iatom, isort, lmaxp, maxdim, ncix, crotloc
  USE defs,  ONLY: IMAG
  IMPLICIT NONE
  INTEGER, intent(in)      :: ig ! The index of the group operation
  COMPLEX*16, intent(out)  :: RT(maxdim,maxdim,ncix)  ! Symmetrization matrix for this group operation
  ! other needed inputs
  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmaxp+1), nindo(norbitals)
  INTEGER, intent(in)      :: norbitals, maxdim2
  LOGICAL, intent(in)      :: Qprint
  ! external
  REAL*8 :: detx
  INTEGER :: FindRotatedAtom
  INTEGER :: FindRotatedCase
  !
  REAL*8     :: BR1inv(3,3), Id(3,3), TR(3,3), Det
  REAL*8     :: local_axis_latom(3,3), local_axis_rotated(3,3), local_axis_rotated_inv(3,3), Tac(3,3)
  INTEGER    :: icase, lcase, icix, nind, iorb, ni, iorb_rotated, iorb1, iorb2, i1, i2, jatom, latom, nind1, nind2, ind1, ind2
  INTEGER    :: l, cdim, i, j, latom_rotated, icase_rotated
  COMPLEX*16 :: Ds(2,2), tmp(maxdim2,maxdim2)
  COMPLEX*16, allocatable :: Dr(:,:), Drs(:,:)
  COMPLEX*16 :: Drx(maxdim2,maxdim2), CfXC(maxdim2,maxdim2)
  !
  call inv_3x3(BR1,BR1inv)
  Id(:,:)=0.d0
  Id(1,1)=1.d0
  Id(2,2)=1.d0
  Id(3,3)=1.d0

  ! TR <=  BR1 * iz * BR1inv, is the rotation due to the group operation written in global coordinate system
  TR(:,:) = iz(:,:,ig)
  TR = matmul(BR1, matmul(TR, BR1inv) )
  
  if (Qprint) then
     WRITE(6,*) 'ig=', ig
     WRITE(6,*) 'iz='
     DO i=1,3
        WRITE(6,'(3I3,1x)') iz(i,:,ig)
     enddo
  endif
  
  RT(:,:,:) = 0.d0
  DO icase=1,natom
     latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
     jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )
     
     latom_rotated = FindRotatedAtom(iz(:,:,ig),tau(:,ig),jatom,latom)
     icase_rotated = FindRotatedCase(latom_rotated)
     !print *, 'icase=', icase, 'icase_rotated=', icase_rotated
     
     local_axis_latom   = matmul(crotloc(:,:,icase),         matmul(BR1, matmul(rotij(:,:,latom),         BR1inv)))
     local_axis_rotated = matmul(crotloc(:,:,icase_rotated), matmul(BR1, matmul(rotij(:,:,latom_rotated), BR1inv)))
     call inv_3x3(local_axis_rotated, local_axis_rotated_inv)
     
     Tac = matmul(local_axis_latom, matmul(TR, local_axis_rotated_inv))
     Det = detx(Tac)
     if (Det < 0.d0) then
        Tac = -Tac
     endif

     !if (Qprint) then
     !  WRITE(6,*) 'icase=', icase, 'Tac='
     !  DO i=1,3
     !     WRITE(6,'(3F10.3,1x)') Tac(i,:)
     !  enddo
     !endif
     
     do lcase=1,nl(icase)
        icix = cix(icase,lcase)
        if ( icix.EQ.0 ) CYCLE
        
        l = ll(icase,lcase)
        nind = (2*l+1)*iso
        iorb = iorbital(icase,lcase)
        iorb_rotated = iorbital(icase_rotated,lcase)
        
        allocate( Dr(nind,nind) )
        Dr = 0.d0
        if (sum(abs(Tac-Id)).LT.1e-10) then
           do i=1,nind   ! Identity requires no rotation
              Dr(i,i) = 1.d0
           enddo
        else
           ni = 2*l+1
           call GetWignerOrbitalMatrix(Dr(:ni,:ni), l, Tac)
           if (iso.eq.2) then ! Spin orbit requires simultaneous spin-rotations
              allocate( Drs(ni,ni) )
              Drs(:,:) = Dr(:ni,:ni)
              call GetWignerSpinMatrix(Ds, 1.d0/2.d0, Tac)
              Dr(:ni,    :ni) = Ds(1,1) * Drs(:,:)
              Dr(:ni,  ni+1:) = Ds(1,2) * Drs(:,:)
              Dr(ni+1:,  :ni) = Ds(2,1) * Drs(:,:)
              Dr(ni+1:,ni+1:) = Ds(2,2) * Drs(:,:)
              deallocate( Drs )
           endif
        endif
        
        if (Det.lt.0 .and. mod(l,2).eq.1) then  ! We have inversion, hence Y_{lm} -> (-1)^l Y_{lm}
           Dr(:,:) = Dr(:,:) * (-1)**l
        endif
        
        do iorb1=1,norbitals
           if (cix_orb(iorb) .NE. cix_orb(iorb1) ) CYCLE
           nind1 = nindo(iorb1)
           !!!!gfortran bug discovered in 2018. It can not handle temporary variables
           ! tmp(:nind1,:nind) = matmul( conjg(cfX(:nind1,:nind,iorb1,iorb)), Dr(:nind,:nind) )
           cfXC(:nind1,:nind) = conjg(cfX(:nind1,:nind,iorb1,iorb))
           tmp(:nind1,:nind) = matmul( cfXC(:nind1,:nind), Dr(:nind,:nind) )
           !!!!gfortran bug discovered in 2018. It can not handle temporary variables
           do iorb2=1,norbitals
              if (cix_orb(iorb_rotated) .NE. cix_orb(iorb2) ) CYCLE
              nind2 = nindo(iorb2)
              Drx(:nind1,:nind2) = matmul(tmp(:nind1,:nind), transpose(cfX(:nind2,:nind,iorb2,iorb_rotated)))                 
              do ind1=1,nind1
                 do ind2=1,nind2
                    i1 = iSx(ind1,iorb1)
                    i2 = iSx(ind2,iorb2)
                    RT(i1,i2,icix) = RT(i1,i2,icix) + Drx(ind1,ind2)
                 enddo
              enddo
           enddo
        enddo
        deallocate( Dr )
     enddo
  ENDDO
     
  if (Qprint) then
     do icix=1,ncix
        cdim = cixdim(icix)
        WRITE(6,*) 'icix=', icix
        do i=1,cdim
           do j=1,cdim
              WRITE(6,'(2F14.7,2x)',advance='no') RT(i,j,icix)
           enddo
           WRITE(6,*)
        enddo
     end do
  endif
END SUBROUTINE CmpSymmetryTransformation


!===============================================================================
! SYMMETRIZATION OF THE K-GLOBAL PROJECTOR OVERLAP
!===============================================================================
!
! The irreducible-k weight used by dmft2 already contains the multiplicity of
! the k-star.  That scalar multiplicity is sufficient for a scalar quantity,
! but not for
!
!                     S = sum_k w_k P_k^+ P_k,
!
! because S is a matrix in the correlated-orbital basis.  If g is a crystal
! symmetry operation and R_g is its representation in that basis, the band
! gauge cancels in P^+P and a missing member of the star is
!
!                     S(g k) = R_g^+ S(k) R_g.
!
! Consequently the full-star result can be reconstructed after the
! irreducible-k sum by applying the group projector below.  This is the same
! operation used by dmft1/SymmetrizeOverlap.  CmpSymmetryTransformation is
! deliberately reused here because it contains all of the nontrivial pieces:
! equivalent-atom mappings, local correlated axes, cfX transformations,
! spin rotations in the SO case, and orbital parity under inversion.
!
SUBROUTINE SymmetrizeOverlap(Olapm, cfX, cix_orb, cixdim, iSx, &
     iorbital, nindo, norbitals, maxdim2)
  USE sym2,   ONLY: iord
  USE dmfts,  ONLY: maxdim, ncix, natom, lmaxp
  USE com_mpi,ONLY: myrank, master
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix)
  COMPLEX*16, intent(in)    :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  INTEGER, intent(in)       :: cix_orb(norbitals), cixdim(ncix)
  INTEGER, intent(in)       :: iSx(maxdim2,norbitals)
  INTEGER, intent(in)       :: iorbital(natom,lmaxp+1), nindo(norbitals)
  INTEGER, intent(in)       :: norbitals, maxdim2
  ! locals
  COMPLEX*16, allocatable :: Olap0(:,:,:), Olaps(:,:,:), RT(:,:,:)
  REAL*8 :: change, herm_before, herm_after
  INTEGER :: ig, icix, cdim, i, j

  allocate(Olap0(maxdim,maxdim,ncix), Olaps(maxdim,maxdim,ncix), &
       RT(maxdim,maxdim,ncix))
  Olap0=Olapm
  Olaps=0.d0

  do ig=1,iord
     call CmpSymmetryTransformation(ig, RT, cfX, cix_orb, cixdim, &
          iSx, iorbital, nindo, norbitals, maxdim2, .False.)
     do icix=1,ncix
        cdim=cixdim(icix)
        Olaps(:cdim,:cdim,icix)=Olaps(:cdim,:cdim,icix)+ &
             matmul(conjg(transpose(RT(:cdim,:cdim,icix))), &
             matmul(Olap0(:cdim,:cdim,icix), &
             RT(:cdim,:cdim,icix)))/dble(iord)
     enddo
  enddo

  ! Keep the treatment of separate, symmetry-equivalent correlated blocks
  ! identical to dmft1.  For the common one-cix-per-sort case this is a no-op.
  call AverageEquivalentCixMatrices(Olaps,cixdim,maxdim,ncix)

  ! S and every derivative of S are Hermitian.  Enforce this only after the
  ! group average, so the assignment removes roundoff rather than hiding a
  ! problem in an individual k contribution.
  herm_before=0.d0
  do icix=1,ncix
     cdim=cixdim(icix)
     do j=1,cdim
        do i=1,cdim
           herm_before=max(herm_before, &
                abs(Olaps(i,j,icix)-conjg(Olaps(j,i,icix))))
        enddo
     enddo
     Olaps(:cdim,:cdim,icix)=0.5d0*(Olaps(:cdim,:cdim,icix)+ &
          conjg(transpose(Olaps(:cdim,:cdim,icix))))
  enddo

  change=maxval(abs(Olaps-Olap0))
  herm_after=0.d0
  do icix=1,ncix
     cdim=cixdim(icix)
     do j=1,cdim
        do i=1,cdim
           herm_after=max(herm_after, &
                abs(Olaps(i,j,icix)-conjg(Olaps(j,i,icix))))
        enddo
     enddo
  enddo
  Olapm=Olaps

  if (myrank.eq.master) then
     write(6,'(A,3(1x,ES14.6))') ':OVERLAP_SYM_SCALAR:',change, &
          herm_before,herm_after
  endif
  deallocate(RT,Olaps,Olap0)
END SUBROUTINE SymmetrizeOverlap


! Symmetrize the normalization functional gradient over the same missing
! k-star members that define the global projector overlap.  The accumulated
! object is momentBar=conjg(N), where N=sum_k U_k^+ C_k.  N is a general
! complex matrix (not a density matrix), so it must not be Hermitianized.
!
! This operation is deliberately diagnostic/default-off for now.  The
! production result is unchanged unless
!
!   EDMFT_STRESS_NORM_SYMMETRY=APPLY
!
! is set.  DIAG evaluates and reports the group projection without applying
! it.  This lets the finite-difference regression establish whether the
! missing adjoint k-star projection is the hydrostatic stress discrepancy
! before it becomes a production correction.
SUBROUTINE SymmetrizeNormalizationMoment(momentBar,cfX,cix_orb,cixdim, &
     iSx,iorbital,nindo,norbitals,maxdim2)
  USE sym2,   ONLY: iord
  USE dmfts,  ONLY: maxdim,ncix,natom,lmaxp
  USE com_mpi,ONLY: myrank,master
  IMPLICIT NONE
  COMPLEX*16,intent(inout) :: momentBar(maxdim,maxdim,ncix)
  COMPLEX*16,intent(in) :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  INTEGER,intent(in) :: cix_orb(norbitals),cixdim(ncix)
  INTEGER,intent(in) :: iSx(maxdim2,norbitals)
  INTEGER,intent(in) :: iorbital(natom,lmaxp+1),nindo(norbitals)
  INTEGER,intent(in) :: norbitals,maxdim2
  COMPLEX*16,allocatable :: Nraw(:,:,:),Nsym(:,:,:),RT(:,:,:)
  CHARACTER(len=32) :: mode
  INTEGER :: status,i,letter,ig,icix,cdim
  REAL*8 :: rawNorm,differenceNorm,relativeDifference,maxChange

  mode='OFF'
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_NORM_SYMMETRY',mode, &
       STATUS=status)
  if (status.ne.0 .or. len_trim(mode).eq.0) mode='OFF'
  mode=adjustl(mode)
  do i=1,len_trim(mode)
     letter=iachar(mode(i:i))
     if (letter.ge.iachar('a') .and. letter.le.iachar('z')) &
          mode(i:i)=achar(letter-iachar('a')+iachar('A'))
  enddo
  if (trim(mode).eq.'OFF') return
  if (trim(mode).ne.'DIAG' .and. trim(mode).ne.'APPLY') then
     WRITE(6,'(A,A)') 'ERROR: unknown EDMFT_STRESS_NORM_SYMMETRY=', &
          trim(mode)
     STOP 'invalid normalization-moment symmetry mode'
  endif

  allocate(Nraw(maxdim,maxdim,ncix),Nsym(maxdim,maxdim,ncix), &
       RT(maxdim,maxdim,ncix))
  ! momentBar is stored conjugated to match the projector-gradient convention.
  ! Undo that convention before applying N(gk)=R_g^+ N(k) R_g; conjugate the
  ! group-projected result only when returning to l2main.
  Nraw=conjg(momentBar)
  Nsym=0.d0
  do ig=1,iord
     CALL CmpSymmetryTransformation(ig,RT,cfX,cix_orb,cixdim,iSx, &
          iorbital,nindo,norbitals,maxdim2,.False.)
     do icix=1,ncix
        cdim=cixdim(icix)
        Nsym(:cdim,:cdim,icix)=Nsym(:cdim,:cdim,icix)+ &
             matmul(conjg(transpose(RT(:cdim,:cdim,icix))), &
             matmul(Nraw(:cdim,:cdim,icix), &
             RT(:cdim,:cdim,icix)))/dble(iord)
     enddo
  enddo
  CALL AverageEquivalentCixMatrices(Nsym,cixdim,maxdim,ncix)

  rawNorm=sqrt(sum(abs(Nraw)**2))
  differenceNorm=sqrt(sum(abs(Nsym-Nraw)**2))
  relativeDifference=differenceNorm/max(rawNorm,tiny(1.d0))
  maxChange=maxval(abs(Nsym-Nraw))
  if (myrank.eq.master) WRITE(6,'(A,1X,A,3(1X,ES20.11))') &
       ':DMFT_NORMALIZATION_MOMENT_SYMMETRY:',trim(mode), &
       rawNorm,relativeDifference,maxChange
  if (trim(mode).eq.'APPLY') momentBar=conjg(Nsym)

  deallocate(RT,Nsym,Nraw)
END SUBROUTINE SymmetrizeNormalizationMoment


! Symmetrize a Cartesian derivative of the overlap.  nq=3 denotes a polar
! vector (the force/displacement response), while nq=6 denotes the symmetric
! strain tensor in the ordering (xx,yy,zz,yz,xz,xy).  The orbital and
! Cartesian transformations must be applied as one direct-product group
! operation; symmetrizing the six local matrices independently is incorrect.
!
! For strain,
!
!   S^sym_{,ab} = 1/|G| sum_g R_g^+ [sum_cd Q_g(a,c)Q_g(b,d)S_{,cd}] R_g,
!
! where Q_g=BR1*iz_g*BR1^{-1} is the operation in global Cartesian
! coordinates.  Off-diagonal entries are tensor components, not engineering
! shear components, so no factor of two is introduced.
!
SUBROUTINE SymmetrizeOverlapDerivative(dS, nq, cfX, cix_orb, cixdim, &
     iSx, iorbital, nindo, norbitals, maxdim2)
  USE structure,ONLY: BR1
  USE sym2,    ONLY: iord, iz
  USE dmfts,   ONLY: maxdim, ncix, natom, lmaxp
  USE com_mpi, ONLY: myrank, master
  IMPLICIT NONE
  INTEGER, intent(in)       :: nq, norbitals, maxdim2
  COMPLEX*16, intent(inout) :: dS(maxdim,maxdim,ncix,nq)
  COMPLEX*16, intent(in)    :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  INTEGER, intent(in)       :: cix_orb(norbitals), cixdim(ncix)
  INTEGER, intent(in)       :: iSx(maxdim2,norbitals)
  INTEGER, intent(in)       :: iorbital(natom,lmaxp+1), nindo(norbitals)
  ! locals
  INTEGER, parameter :: voigt_a(6)=(/1,2,3,2,1,1/)
  INTEGER, parameter :: voigt_b(6)=(/1,2,3,3,3,2/)
  COMPLEX*16, allocatable :: dS0(:,:,:,:), dSs(:,:,:,:), rotated(:,:,:,:)
  COMPLEX*16, allocatable :: RT(:,:,:)
  REAL*8 :: BR1inv(3,3), Qcart(3,3), coeff
  REAL*8 :: change, herm_before, herm_after
  INTEGER :: ig, icix, cdim, iq, jq, a, b, c, d, i, j

  if (nq.ne.3 .and. nq.ne.6) then
     write(6,*) 'ERROR SymmetrizeOverlapDerivative: nq must be 3 or 6, nq=',nq
     stop 'Invalid Cartesian rank in SymmetrizeOverlapDerivative'
  endif

  allocate(dS0(maxdim,maxdim,ncix,nq), dSs(maxdim,maxdim,ncix,nq), &
       rotated(maxdim,maxdim,ncix,nq), RT(maxdim,maxdim,ncix))
  dS0=dS
  dSs=0.d0
  call inv_3x3(BR1,BR1inv)

  do ig=1,iord
     call CmpSymmetryTransformation(ig, RT, cfX, cix_orb, cixdim, &
          iSx, iorbital, nindo, norbitals, maxdim2, .False.)

     ! This is the same Cartesian operation used to construct the Wigner and
     ! local-axis transformation inside CmpSymmetryTransformation.
     Qcart=matmul(BR1,matmul(dble(iz(:,:,ig)),BR1inv))
     rotated=0.d0

     if (nq.eq.3) then
        ! Polar-vector response: dS'_a=sum_c Q_ac dS_c.
        do iq=1,3
           do jq=1,3
              rotated(:,:,:,iq)=rotated(:,:,:,iq)+ &
                   Qcart(iq,jq)*dS0(:,:,:,jq)
           enddo
        enddo
     else
        ! Expand the six stored symmetric components implicitly.  For c/=d,
        ! both dS_cd and dS_dc equal the one stored tensor component.
        do iq=1,6
           a=voigt_a(iq)
           b=voigt_b(iq)
           do jq=1,6
              c=voigt_a(jq)
              d=voigt_b(jq)
              if (c.eq.d) then
                 coeff=Qcart(a,c)*Qcart(b,d)
              else
                 coeff=Qcart(a,c)*Qcart(b,d)+ &
                      Qcart(a,d)*Qcart(b,c)
              endif
              rotated(:,:,:,iq)=rotated(:,:,:,iq)+ &
                   coeff*dS0(:,:,:,jq)
           enddo
        enddo
     endif

     do iq=1,nq
        do icix=1,ncix
           cdim=cixdim(icix)
           dSs(:cdim,:cdim,icix,iq)=dSs(:cdim,:cdim,icix,iq)+ &
                matmul(conjg(transpose(RT(:cdim,:cdim,icix))), &
                matmul(rotated(:cdim,:cdim,icix,iq), &
                RT(:cdim,:cdim,icix)))/dble(iord)
        enddo
     enddo
  enddo

  do iq=1,nq
     call AverageEquivalentCixMatrices(dSs(:,:,:,iq),cixdim,maxdim,ncix)
  enddo

  herm_before=0.d0
  do iq=1,nq
     do icix=1,ncix
        cdim=cixdim(icix)
        do j=1,cdim
           do i=1,cdim
              herm_before=max(herm_before,abs(dSs(i,j,icix,iq)- &
                   conjg(dSs(j,i,icix,iq))))
           enddo
        enddo
        dSs(:cdim,:cdim,icix,iq)=0.5d0*(dSs(:cdim,:cdim,icix,iq)+ &
             conjg(transpose(dSs(:cdim,:cdim,icix,iq))))
     enddo
  enddo

  change=maxval(abs(dSs-dS0))
  herm_after=0.d0
  do iq=1,nq
     do icix=1,ncix
        cdim=cixdim(icix)
        do j=1,cdim
           do i=1,cdim
              herm_after=max(herm_after,abs(dSs(i,j,icix,iq)- &
                   conjg(dSs(j,i,icix,iq))))
           enddo
        enddo
     enddo
  enddo
  dS=dSs

  if (myrank.eq.master) then
     if (nq.eq.3) then
        write(6,'(A,3(1x,ES14.6))') ':OVERLAP_SYM_VECTOR:',change, &
             herm_before,herm_after
     else
        write(6,'(A,3(1x,ES14.6))') ':OVERLAP_SYM_STRAIN:',change, &
             herm_before,herm_after
     endif
  endif
  deallocate(RT,rotated,dSs,dS0)
END SUBROUTINE SymmetrizeOverlapDerivative


! Average separate cix blocks that represent symmetry-equivalent atoms.  The
! orbital rotation has already put every block into the convention expected by
! dmft1, so this last average is meaningful and reproduces its historical
! treatment of equivalent correlated blocks.
SUBROUTINE AverageEquivalentCixMatrices(A,cixdim,maxdim,ncix)
  IMPLICIT NONE
  INTEGER, intent(in)       :: maxdim, ncix, cixdim(ncix)
  COMPLEX*16, intent(inout) :: A(maxdim,maxdim,ncix)
  INTEGER :: ncix_unique, size_cix_unique(ncix), cix_unique(ncix,ncix)
  INTEGER :: igroup, j, icix, cdim
  COMPLEX*16, allocatable :: tmp(:,:)

  call FindUniqueCix(ncix_unique,size_cix_unique,cix_unique)
  do igroup=1,ncix_unique
     if (size_cix_unique(igroup).le.1) cycle
     icix=cix_unique(igroup,1)
     cdim=cixdim(icix)
     allocate(tmp(cdim,cdim))
     tmp=0.d0
     do j=1,size_cix_unique(igroup)
        icix=cix_unique(igroup,j)
        if (cixdim(icix).ne.cdim) then
           write(6,*) 'ERROR: equivalent cix blocks have different dimensions'
           stop 'Inconsistent equivalent cix dimensions'
        endif
        tmp=tmp+A(:cdim,:cdim,icix)
     enddo
     tmp=tmp/dble(size_cix_unique(igroup))
     do j=1,size_cix_unique(igroup)
        icix=cix_unique(igroup,j)
        A(:cdim,:cdim,icix)=tmp
     enddo
     deallocate(tmp)
  enddo
END SUBROUTINE AverageEquivalentCixMatrices


SUBROUTINE FindUniqueCix(ncix_unique,size_cix_unique,cix_unique)
  USE dmfts, ONLY: natom, ncix, iatom, isort, cix
  IMPLICIT NONE
  INTEGER, intent(out) :: ncix_unique, size_cix_unique(ncix)
  INTEGER, intent(out) :: cix_unique(ncix,ncix)
  INTEGER :: cix_2_jatom(ncix), group_sort(ncix)
  INTEGER :: icix, icase, latom, jatom, igroup
  LOGICAL :: found

  ! Defaulting to a unique negative label avoids undefined grouping if a cix
  ! is intentionally not attached to the first correlated shell of a case.
  do icix=1,ncix
     cix_2_jatom(icix)=-icix
  enddo
  do icix=1,ncix
     do icase=1,natom
        latom=iatom(icase)
        jatom=isort(latom)
        if (icix.eq.cix(icase,1)) then
           cix_2_jatom(icix)=jatom
           exit
        endif
     enddo
  enddo

  ncix_unique=0
  size_cix_unique=0
  cix_unique=0
  group_sort=0
  do icix=1,ncix
     found=.False.
     do igroup=1,ncix_unique
        if (group_sort(igroup).eq.cix_2_jatom(icix)) then
           size_cix_unique(igroup)=size_cix_unique(igroup)+1
           cix_unique(igroup,size_cix_unique(igroup))=icix
           found=.True.
           exit
        endif
     enddo
     if (.not.found) then
        ncix_unique=ncix_unique+1
        group_sort(ncix_unique)=cix_2_jatom(icix)
        size_cix_unique(ncix_unique)=1
        cix_unique(ncix_unique,1)=icix
     endif
  enddo
END SUBROUTINE FindUniqueCix


SUBROUTINE SymmetrizeDensityMatrix(DM, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2)
  USE sym2,  ONLY: iord
  USE dmfts, ONLY: maxdim, ncix, natom, lmaxp
  USE defs,  ONLY: IMAG
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: DM(maxdim,maxdim,ncix)
  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmaxp+1), nindo(norbitals)
  INTEGER, intent(in)      :: norbitals, maxdim2
  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  ! locals
  COMPLEX*16 :: RT(maxdim,maxdim,ncix)
  INTEGER    :: icix, ig
  INTEGER    :: cdim
  COMPLEX*16, allocatable :: aDM(:,:,:)
  
  allocate( aDM(maxdim,maxdim,ncix) )
  aDM(:,:,:) = 0.d0
  
  do ig=1,iord
     Call CmpSymmetryTransformation(ig, RT, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2, .False.)
     do icix=1,ncix
        cdim = cixdim(icix)
        aDM(:,:,icix) = aDM(:,:,icix) + matmul(conjg(transpose(RT(:cdim,:cdim,icix))), matmul(DM(:cdim,:cdim,icix), RT(:cdim,:cdim,icix)))/iord
     end do
  enddo
  DM(:,:,:) = aDM(:,:,:)
END SUBROUTINE SymmetrizeDensityMatrix

SUBROUTINE SymmetrizeLocalQuantities(Gdloc0, Edimp0, Gdloc, Edimp, cfX, cix_orb, cixdim, iSx, iorbital, nindo, lg_deg, norbitals, maxdim2, nomega)
  USE sym2,  ONLY: iord
  USE dmfts, ONLY: maxdim, ncix, natom, lmaxp, ntcix, Sigind_orig
  USE defs,  ONLY: IMAG
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: Gdloc0(maxdim,maxdim,ncix,nomega), Edimp0(maxdim,maxdim,ncix)
  COMPLEX*16, intent(out)  :: Gdloc(ntcix,nomega)
  REAL*8, intent(out)      :: Edimp(ntcix)
  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmaxp+1), nindo(norbitals), lg_deg(ntcix)
  INTEGER, intent(in)      :: norbitals, maxdim2, nomega
  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  ! locals
  INTEGER    :: icix, ig, ip, iq, iom, it2
  INTEGER    :: cdim
  COMPLEX*16 :: gc(maxdim,maxdim)
  COMPLEX*16, allocatable :: RT(:,:,:,:)
  
  allocate( RT(maxdim,maxdim,ncix,iord) )
  
  do ig=1,iord
     Call CmpSymmetryTransformation(ig, RT(:,:,:,ig), cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2, .False.)
  enddo

  Gdloc(:,:) = 0.d0
  Edimp=0.d0
  do icix=1,ncix
     cdim = cixdim(icix)
     do iom=1,nomega
        ! Symmetrize the local Green's function
        gc(:,:)=0.d0
        do ig=1,iord
           gc(:cdim,:cdim) = gc(:cdim,:cdim) + matmul(conjg(transpose(RT(:cdim,:cdim,icix,ig))), matmul(Gdloc0(:cdim,:cdim,icix,iom), RT(:cdim,:cdim,icix,ig)))/iord
        end do
        Gdloc0(:cdim,:cdim,icix,iom) = gc(:cdim,:cdim)
     
        do ip=1,cdim
           do iq=1,cdim
              it2 = Sigind_orig( ip, iq, icix )
              if (it2.gt.0) then
                 Gdloc(it2,iom) = Gdloc(it2,iom) + gc(ip,iq)/lg_deg(it2)
              endif
           enddo
        enddo
     enddo
     
     ! Symmetrize the impurity levels
     gc(:,:)=0.d0
     do ig=1,iord
        gc(:cdim,:cdim) = gc(:cdim,:cdim) + matmul(conjg(transpose(RT(:cdim,:cdim,icix,ig))), matmul(Edimp0(:cdim,:cdim,icix), RT(:cdim,:cdim,icix,ig)))/iord
     end do
     Edimp0(:cdim,:cdim,icix) = gc(:cdim,:cdim)
     do ip=1,cdim
        do iq=1,cdim
           it2 = Sigind_orig( ip, iq, icix )
           if (it2.gt.0) then
              Edimp(it2) = Edimp(it2) + dble(gc(ip,iq))/lg_deg(it2)
           endif
        enddo
     enddo
  enddo
  
  deallocate( RT )
END SUBROUTINE SymmetrizeLocalQuantities

SUBROUTINE GetNds(Nds, DM, cixdim, lg_deg)
  use dmfts, ONLY: ncix, ntcix, maxdim, Sigind_orig
  IMPLICIT NONE
  REAL*8, intent(out)    :: Nds(ntcix)
  COMPLEX*16, intent(in) :: DM(maxdim,maxdim,ncix)
  INTEGER, intent(in)    :: cixdim(ncix), lg_deg(ntcix)
  ! locals
  INTEGER :: cdim, i, j, icix, it2
  !
  Nds(:)=0.d0
  do icix=1,ncix
     cdim = cixdim(icix)
     do i=1,cdim
        do j=1,cdim
           it2 = Sigind_orig( i, j, icix )
           if (it2.gt.0) then
              Nds(it2) = Nds(it2) + dble(DM(i,j, icix)) !/lg_deg(it2)
              !Nds(it2,2:4) = Nds(it2,2:4) + (DM(i,j,icix,2:4,1)-DM(i,j,icix,1,2:4))*IMAG
           endif
        enddo
     enddo
  enddo
  !do it2=1,ntcix
  !   Nds(it2) = Nds(it2)/deg(it2)
  !enddo
END SUBROUTINE GetNds

                       


!SUBROUTINE SymmetrizeDensityMatrix_old(DM, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2)
!  USE structure, ONLY: BR1, rotij
!  USE sym2,  ONLY: tau, iz, iord
!  USE dmfts, ONLY: iso, natom, nl, ll, cix, iatom, isort, lmaxp, maxdim, ncix, crotloc
!  USE defs,  ONLY: IMAG
!  IMPLICIT NONE
!  COMPLEX*16, intent(inout):: DM(maxdim,maxdim,ncix)
!  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmaxp+1), nindo(norbitals)
!  INTEGER, intent(in)      :: norbitals, maxdim2
!  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
!  ! external
!  REAL*8 :: detx
!  INTEGER :: FindRotatedAtom
!  INTEGER :: FindRotatedCase
!  !
!  COMPLEX*16 :: RT(maxdim,maxdim,ncix)
!  REAL*8     :: BR1inv(3,3), Id(3,3), TR(3,3), Det
!  REAL*8     :: local_axis_latom(3,3), local_axis_rotated(3,3), local_axis_rotated_inv(3,3), Tac(3,3)
!  INTEGER    :: icase, lcase, icix, nind, iorb, ni, iorb_rotated, iorb1, iorb2, i1, i2, jatom, latom, nind1, nind2, ig, ind1, ind2
!  INTEGER    :: l, cdim, i, j, latom_rotated, icase_rotated
!  COMPLEX*16 :: Ds(2,2), tmp(maxdim2,maxdim2)
!  COMPLEX*16, allocatable :: Dr(:,:), Drs(:,:)
!  COMPLEX*16 :: Drx(maxdim2,maxdim2)
!  !
!  COMPLEX*16, allocatable :: aDM(:,:,:)
!  
!  call inv_3x3(BR1,BR1inv)
!  Id(:,:)=0.d0
!  Id(1,1)=1.d0
!  Id(2,2)=1.d0
!  Id(3,3)=1.d0
!
!  allocate( aDM(maxdim,maxdim,ncix) )
!  aDM(:,:,:) = 0.d0
!  
!  do ig=1,iord
!     ! TR <=  BR1 * iz * BR1inv, is the rotation due to the group operation written in global coordinate system
!     TR(:,:) = iz(:,:,ig)
!     TR = matmul(BR1, matmul(TR, BR1inv) )
!     
!     WRITE(6,*) 'ig=', ig
!     WRITE(6,*) 'iz='
!     DO i=1,3
!        WRITE(6,'(3I3,1x)') iz(i,:,ig)
!     enddo
!     
!     RT(:,:,:) = 0.d0
!     DO icase=1,natom
!        latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
!        jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )
!        
!        latom_rotated = FindRotatedAtom(iz(:,:,ig),tau(:,ig),jatom,latom)
!        icase_rotated = FindRotatedCase(latom_rotated)
!        !print *, 'icase=', icase, 'icase_rotated=', icase_rotated
!        
!        local_axis_latom   = matmul(crotloc(:,:,icase),         matmul(BR1, matmul(rotij(:,:,latom),         BR1inv)))
!        local_axis_rotated = matmul(crotloc(:,:,icase_rotated), matmul(BR1, matmul(rotij(:,:,latom_rotated), BR1inv)))
!        call inv_3x3(local_axis_rotated, local_axis_rotated_inv)
!        
!        Tac = matmul(local_axis_latom, matmul(TR, local_axis_rotated_inv))
!        Det = detx(Tac)
!        if (Det < 0.d0) then
!           Tac = -Tac
!        endif
!        
!        !WRITE(6,*) 'icase=', icase, 'Tac='
!        !DO i=1,3
!        !   WRITE(6,'(3F10.3,1x)') Tac(i,:)
!        !enddo
!        
!        do lcase=1,nl(icase)
!           icix = cix(icase,lcase)
!           if ( icix.EQ.0 ) CYCLE
!
!           l = ll(icase,lcase)
!           nind = (2*l+1)*iso
!           iorb = iorbital(icase,lcase)
!           iorb_rotated = iorbital(icase_rotated,lcase)
!
!           allocate( Dr(nind,nind) )
!           Dr = 0.d0
!           if (sum(abs(Tac-Id)).LT.1e-10) then
!              do i=1,nind   ! Identity requires no rotation
!                 Dr(i,i) = 1.d0
!              enddo
!           else
!              ni = 2*l+1
!              call GetWignerOrbitalMatrix(Dr(:ni,:ni), l, Tac)
!              if (iso.eq.2) then ! Spin orbit requires simultaneous spin-rotations
!                 allocate( Drs(ni,ni) )
!                 Drs(:,:) = Dr(:ni,:ni)
!                 call GetWignerSpinMatrix(Ds, 1.d0/2.d0, Tac)
!                 Dr(:ni,    :ni) = Ds(1,1) * Drs(:,:)
!                 Dr(:ni,  ni+1:) = Ds(1,2) * Drs(:,:)
!                 Dr(ni+1:,  :ni) = Ds(2,1) * Drs(:,:)
!                 Dr(ni+1:,ni+1:) = Ds(2,2) * Drs(:,:)
!                 deallocate( Drs )
!              endif
!           endif
!        
!           if (Det.lt.0 .and. mod(l,2).eq.1) then  ! We have inversion, hence Y_{lm} -> (-1)^l Y_{lm}
!              Dr(:,:) = Dr(:,:) * (-1)**l
!           endif
!
!           do iorb1=1,norbitals
!              if (cix_orb(iorb) .NE. cix_orb(iorb1) ) CYCLE
!              nind1 = nindo(iorb1)
!              tmp(:nind1,:nind) = matmul( conjg(cfX(:nind1,:nind,iorb1,iorb)), Dr(:nind,:nind) )
!              do iorb2=1,norbitals
!                 if (cix_orb(iorb_rotated) .NE. cix_orb(iorb2) ) CYCLE
!                 nind2 = nindo(iorb2)
!                 Drx(:nind1,:nind2) = matmul(tmp(:nind1,:nind), transpose(cfX(:nind2,:nind,iorb2,iorb_rotated)))                 
!                 do ind1=1,nind1
!                    do ind2=1,nind2
!                       i1 = iSx(ind1,iorb1)
!                       i2 = iSx(ind2,iorb2)
!                       RT(i1,i2,icix) = RT(i1,i2,icix) + Drx(ind1,ind2)
!                    enddo
!                 enddo
!              enddo
!           enddo
!           deallocate( Dr )
!        enddo
!     ENDDO
!
!     do icix=1,ncix
!        cdim = cixdim(icix)
!        WRITE(6,*) 'icix=', icix
!        do i=1,cdim
!           do j=1,cdim
!              WRITE(6,'(2F14.7,2x)',advance='no') RT(i,j,icix)
!           enddo
!           WRITE(6,*)
!        enddo
!        aDM(:,:,icix) = aDM(:,:,icix) + matmul(conjg(transpose(RT(:cdim,:cdim,icix))), matmul(DM(:cdim,:cdim,icix), RT(:cdim,:cdim,icix)))/iord
!     end do
!  enddo
!  DM(:,:,:) = aDM(:,:,:)
!END SUBROUTINE SymmetrizeDensityMatrix_Old
!
!
!SUBROUTINE SymmetrizeTest(cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2)
!  USE structure, ONLY: BR1, rotij
!  USE sym2,  ONLY: tau, iz, iord
!  USE dmfts, ONLY: iso, natom, nl, ll, cix, iatom, isort, lmaxp, maxdim, ncix, crotloc
!  USE defs,  ONLY: IMAG
!  IMPLICIT NONE
!  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmaxp+1), nindo(norbitals)
!  INTEGER, intent(in)      :: norbitals, maxdim2
!  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
!  ! external
!  REAL*8 :: detx
!  INTEGER :: FindRotatedAtom
!  INTEGER :: FindRotatedCase
!  !
!  COMPLEX*16 :: RT(maxdim,maxdim,ncix)
!  REAL*8     :: BR1inv(3,3), Id(3,3), tmp3(3,3), TR(3,3), Det
!  REAL*8     :: local_axis_latom(3,3), local_axis_rotated(3,3), local_axis_rotated_inv(3,3), local_axis_latom_inv(3,3)
!  INTEGER    :: icase, lcase, icix, nind, iorb, ni, iorb_rotated, iorb1, iorb2, i1, i2, jatom, latom, nind1, nind2, ig, ind1, ind2
!  INTEGER    :: l, cdim, i, j, latom_rotated, icase_rotated, ii
!  COMPLEX*16 :: Ds(2,2)
!  COMPLEX*16, allocatable :: Dr(:,:), Dr2(:,:)
!  COMPLEX*16 :: Drx(maxdim2,maxdim2)
!  !
!  COMPLEX*16 :: T2C(5,5), Dcubic(5,5), Dx(5,5)
!  REAL*8 :: s2
!  !
!  s2 = 1.d0/sqrt(2.)
!  T2C(:,:) = 0.d0
!  T2C(1,2) = s2      ; T2C(1,4) = -s2
!  T2C(2,2) = s2*IMAG ; T2C(2,4) = s2*IMAG;
!  T2C(3,1) = s2*IMAG ; T2C(3,5) = -s2*IMAG;
!  T2C(4,3) = 1.0d0
!  T2C(5,1) = s2      ; T2C(5,5) = s2;
!
!  !call INVERSSYMDEF(BR1,BR1inv)
!  call inv_3x3(BR1,BR1inv)
!  Id(:,:)=0.d0
!  Id(1,1)=1.d0
!  Id(2,2)=1.d0
!  Id(3,3)=1.d0
!
!
!  l=2
!  allocate( Dr(2*l+1,2*l+1) )
!  DO icase=1,natom
!     latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
!     jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )
!     iorb = iorbital(icase,1)
!
!           
!     local_axis_latom   = matmul(crotloc(:,:,icase),         matmul(BR1, matmul(rotij(:,:,latom),         BR1inv)))
!     call inv_3x3(local_axis_latom, local_axis_latom_inv)
!
!     Det = detx(local_axis_latom)  ! possible inversion
!     WRITE(6,*) 'local_axis_latom with Det=', Det
!     DO i=1,3
!        WRITE(6,'(3F10.3,1x)') local_axis_latom(i,:)
!     enddo
!     
!     call GetWignerOrbitalMatrix(Dr, l, local_axis_latom*Det)
!     
!     WRITE(6,*) 'The corresponding Wigner rotation is'
!     DO i=1,5
!        WRITE(6,'(10F10.3,1x)') Dr(i,:)
!     enddo
!     
!     WRITE(6,*) 'The corresponding Wigner rotation in cubic harmonics is'
!     Dcubic = matmul(matmul(T2C, Dr), transpose(conjg(T2C)))
!
!     DO i=1,5
!        WRITE(6,'(10F10.3,1x)') Dcubic(i,:)
!     ENDDO
!     
!     do ii=1,1
!        !Dx = matmul( matmul(  cfX(:,:,ii,iorb), conjg(transpose(T2C)) ), Dcubic )
!        !Dx = matmul(matmul( conjg(cfX(:,:,ii,iorb)), Dr ), transpose(T2C))
!        Dx = matmul(matmul( cfX(:,:,ii,iorb), conjg(Dr) ), transpose(conjg(T2C)))
!        WRITE(6,*) 'cfX * Dr in cubic harmonics for iorb=', iorb
!        DO i=1,5
!           WRITE(6,'(10F10.3,1x)') Dx(i,:)*2.d0
!        ENDDO
!        WRITE(6,*)
!     enddo
!  ENDDO
!  deallocate( Dr )
!  return
!
!
!  
!  do ig=1,iord
!     ! RT^+ N RT
!     TR(:,:) = iz(:,:,ig)
!     Det = detx(TR)  ! possible inversion
!     TR = TR*Det     ! pure rotation
!     tmp3 = matmul(BR1, TR)
!     TR = matmul(tmp3,BR1inv)
!
!     WRITE(6,*) 'ig=', ig
!     WRITE(6,*) 'iz='
!     DO i=1,3
!        WRITE(6,'(3I3,1x)') iz(i,:,ig)
!     enddo
!     
!     if (sum(abs(TR-Id)).LT.1e-10) TR(:,:)=Id(:,:)
!     RT(:,:,:) = 0.d0
!     DO icase=1,natom
!        latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
!        jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )
!
!        latom_rotated = FindRotatedAtom(iz(:,:,ig),tau(:,ig),jatom,latom)
!        icase_rotated = FindRotatedCase(latom_rotated)
!        print *, 'icase=', icase, 'icase_rotated=', icase_rotated
!        
!        local_axis_latom   = matmul(crotloc(:,:,icase),         matmul(BR1, matmul(rotij(:,:,latom),         BR1inv)))
!        local_axis_rotated = matmul(crotloc(:,:,icase_rotated), matmul(BR1, matmul(rotij(:,:,latom_rotated), BR1inv)))
!        call inv_3x3(local_axis_rotated, local_axis_rotated_inv)
!        call inv_3x3(local_axis_latom, local_axis_latom_inv)
!        
!        TR = matmul(local_axis_latom, matmul(TR, local_axis_rotated_inv))
!        !TR = matmul(local_axis_latom_inv, matmul(TR, local_axis_rotated))
!        WRITE(6,*) 'icase=', icase, 'TR='
!        DO i=1,3
!           WRITE(6,'(3F10.3,1x)') TR(i,:)
!        enddo
!        
!        do lcase=1,nl(icase)
!           icix = cix(icase,lcase)
!           if ( icix.EQ.0 ) CYCLE
!
!           l = ll(icase,lcase)
!           nind = (2*l+1)*iso
!           iorb = iorbital(icase,lcase)
!           allocate( Dr(nind,nind) )
!           Dr=0.d0
!           ni = 2*l+1
!           call GetWignerOrbitalMatrix(Dr(:ni,:ni), l, TR)
!           if (iso.eq.2) then
!              allocate( Dr2(ni,ni) )
!              Dr2(:,:) = Dr(:ni,:ni)
!              call GetWignerSpinMatrix(Ds, 1.d0/2.d0, TR)
!              Dr(:ni,:ni) = Ds(1,1) * Dr2(:,:)
!              Dr(ni:,ni:) = Ds(2,2) * Dr2(:,:)
!              Dr(:ni,ni:) = Ds(1,2) * Dr2(:,:)
!              Dr(ni:,:ni) = Ds(2,1) * Dr2(:,:)
!              deallocate( Dr2 )
!           endif
!           iorb_rotated = iorbital(icase_rotated,lcase)
!
!           if (Det.lt.0 .and. mod(l,2).eq.1) then  ! We have inversion, hence Y_{lm} -> (-1)^l Y_{lm}
!              Dr(:,:) = Dr(:,:) * (-1)**l
!           endif
!           
!           do iorb1=1,norbitals
!              if (cix_orb(iorb) .NE. cix_orb(iorb1) ) CYCLE
!              nind1 = nindo(iorb1)
!              do iorb2=1,norbitals
!                 if (cix_orb(iorb_rotated) .NE. cix_orb(iorb2) ) CYCLE
!                 nind2 = nindo(iorb2)
!                 Drx(:nind1,:nind2) = matmul( conjg(cfX(:nind1,:nind,iorb1,iorb)) , matmul( Dr , transpose(cfX(:nind2,:nind,iorb2,iorb_rotated)) ) )
!                 do ind1=1,nind1
!                    do ind2=1,nind2
!                       i1 = iSx(ind1,iorb1)
!                       i2 = iSx(ind2,iorb2)
!                       RT(i1,i2,icix) = RT(i1,i2,icix) + Drx(ind1,ind2)
!                    enddo
!                 enddo
!              enddo
!           enddo
!           
!           Det = detx(local_axis_latom)  ! possible inversion
!           WRITE(6,*) 'local_axis_latom with Det=', Det
!           DO i=1,3
!              WRITE(6,'(3F10.3,1x)') local_axis_latom(i,:)
!           enddo
!           call GetWignerOrbitalMatrix(Dr(:ni,:ni), l, local_axis_latom*Det)
!           WRITE(6,*) 'The corresponding Wigner rotation is'
!           DO i=1,5
!              WRITE(6,'(10F10.3,1x)') Dr(i,:)
!           enddo
!           WRITE(6,*) 'The corresponding Wigner rotation in cubic harmonics is'
!           Dcubic = matmul(matmul(T2C, Dr), transpose(conjg(T2C)))
!           DO i=1,5
!              WRITE(6,'(10F10.3,1x)') Dcubic(i,:)
!           ENDDO
!           do iorb1=1,4
!              Dcubic = matmul( matmul(  cfX(:nind,:nind1,iorb,iorb1) , Dr), conjg(transpose(T2C)) )
!              WRITE(6,*) 'cfX * Dr in cubic harmonics for iorb=', iorb1
!              DO i=1,5
!                 WRITE(6,'(10F10.3,1x)') Dcubic(i,:)*2.d0
!              ENDDO
!              WRITE(6,*)
!           enddo
!           
!           deallocate( Dr )
!        enddo
!     ENDDO
!
!     do icix=1,ncix
!        cdim = cixdim(icix)
!        WRITE(6,*) 'icix=', icix
!        do i=1,cdim
!           do j=1,cdim
!              WRITE(6,'(2F14.7,2x)',advance='no') RT(i,j,icix)
!           enddo
!           WRITE(6,*)
!        enddo
!     end do
!  enddo
!END SUBROUTINE SymmetrizeTest
!
!
!SUBROUTINE SymmetrizeTest2(cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2)
!  USE structure, ONLY: BR1, rotij
!  USE sym2,  ONLY: tau, iz, iord
!  USE dmfts, ONLY: iso, natom, nl, ll, cix, iatom, isort, lmaxp, maxdim, ncix, crotloc
!  USE defs,  ONLY: IMAG
!  IMPLICIT NONE
!  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmaxp+1), nindo(norbitals)
!  INTEGER, intent(in)      :: norbitals, maxdim2
!  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
!  ! external
!  REAL*8 :: detx
!  INTEGER :: FindRotatedAtom
!  INTEGER :: FindRotatedCase
!  !
!  COMPLEX*16 :: RT(maxdim,maxdim,ncix)
!  REAL*8     :: BR1inv(3,3), Id(3,3), TR(3,3), Det
!  REAL*8     :: local_axis_latom(3,3), local_axis_rotated(3,3), local_axis_rotated_inv(3,3), Tac(3,3) !Ta1(3,3), Ta2(3,3)!, local_axis_latom_inv(3,3)
!  INTEGER    :: icase, lcase, icix, nind, iorb, ni, iorb_rotated, iorb1, iorb2, i1, i2, jatom, latom, nind1, nind2, ig, ind1, ind2
!  INTEGER    :: l, cdim, i, j, latom_rotated, icase_rotated
!  COMPLEX*16 :: Ds(2,2), tmp(maxdim2,maxdim2)
!  COMPLEX*16, allocatable :: Dr(:,:), Drs(:,:)
!  COMPLEX*16 :: Drx(maxdim2,maxdim2)
!  !
!  !INTEGER :: ii
!  !REAL*8  :: Det1, Det2
!  !COMPLEX*16, allocatable :: Dr1(:,:), Dr2(:,:)!, DD(:,:), Drc(:,:) !, C1t(:,:), C2t(:,:), C1(:,:,:), C2(:,:,:)
!  !
!  COMPLEX*16 :: T2C(5,5)!, Df(5,5)
!  REAL*8 :: s2
!
!  s2 = 1.d0/sqrt(2.)
!  
!  T2C(:,:) = 0.d0
!  T2C(1,2) = s2      ; T2C(1,4) = -s2
!  T2C(2,2) = s2*IMAG ; T2C(2,4) = s2*IMAG;
!  T2C(3,1) = s2*IMAG ; T2C(3,5) = -s2*IMAG;
!  T2C(4,3) = 1.0d0
!  T2C(5,1) = s2      ; T2C(5,5) = s2;
!  
!  call inv_3x3(BR1,BR1inv)
!  Id(:,:)=0.d0
!  Id(1,1)=1.d0
!  Id(2,2)=1.d0
!  Id(3,3)=1.d0
!  
!  do ig=1,iord
!     ! TR <=  BR1 * iz * BR1inv, is the rotation due to the group operation written in global coordinate system
!     TR(:,:) = iz(:,:,ig)
!     TR = matmul(BR1, matmul(TR, BR1inv) )
!     
!     WRITE(6,*) 'ig=', ig
!     WRITE(6,*) 'iz='
!     DO i=1,3
!        WRITE(6,'(3I3,1x)') iz(i,:,ig)
!     enddo
!     
!     RT(:,:,:) = 0.d0
!     DO icase=1,natom
!        latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
!        jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )
!        
!        latom_rotated = FindRotatedAtom(iz(:,:,ig),tau(:,ig),jatom,latom)
!        icase_rotated = FindRotatedCase(latom_rotated)
!        print *, 'icase=', icase, 'icase_rotated=', icase_rotated
!        
!        local_axis_latom   = matmul(crotloc(:,:,icase),         matmul(BR1, matmul(rotij(:,:,latom),         BR1inv)))
!        local_axis_rotated = matmul(crotloc(:,:,icase_rotated), matmul(BR1, matmul(rotij(:,:,latom_rotated), BR1inv)))
!        call inv_3x3(local_axis_rotated, local_axis_rotated_inv)
!        
!        Tac = matmul(local_axis_latom, matmul(TR, local_axis_rotated_inv))
!        Det = detx(Tac)
!        if (Det < 0.d0) then
!           Tac = -Tac
!        endif
!        
!        WRITE(6,*) 'icase=', icase, 'Tac='
!        DO i=1,3
!           WRITE(6,'(3F10.3,1x)') Tac(i,:)
!        enddo
!        
!        do lcase=1,nl(icase)
!           icix = cix(icase,lcase)
!           if ( icix.EQ.0 ) CYCLE
!
!           l = ll(icase,lcase)
!           nind = (2*l+1)*iso
!           iorb = iorbital(icase,lcase)
!           iorb_rotated = iorbital(icase_rotated,lcase)
!
!           allocate( Dr(nind,nind) )
!           Dr = 0.d0
!           if (sum(abs(Tac-Id)).LT.1e-10) then
!              do i=1,nind
!                 Dr(i,i) = 1.d0
!              enddo
!           else
!              ni = 2*l+1
!              call GetWignerOrbitalMatrix(Dr(:ni,:ni), l, Tac)
!              if (iso.eq.2) then
!                 allocate( Drs(ni,ni) )
!                 Drs(:,:) = Dr(:ni,:ni)
!                 call GetWignerSpinMatrix(Ds, 1.d0/2.d0, Tac)
!                 Dr(:ni,:ni) = Ds(1,1) * Drs(:,:)
!                 Dr(ni:,ni:) = Ds(2,2) * Drs(:,:)
!                 Dr(:ni,ni:) = Ds(1,2) * Drs(:,:)
!                 Dr(ni:,:ni) = Ds(2,1) * Drs(:,:)
!                 deallocate( Drs )
!              endif
!           endif
!        
!           if (Det.lt.0 .and. mod(l,2).eq.1) then  ! We have inversion, hence Y_{lm} -> (-1)^l Y_{lm}
!              Dr(:,:) = Dr(:,:) * (-1)**l
!           endif
!
!           do iorb1=1,norbitals
!              if (cix_orb(iorb) .NE. cix_orb(iorb1) ) CYCLE
!              nind1 = nindo(iorb1)
!              tmp(:nind1,:nind) = matmul( conjg(cfX(:nind1,:nind,iorb1,iorb)), Dr(:nind,:nind) )
!              do iorb2=1,norbitals
!                 if (cix_orb(iorb_rotated) .NE. cix_orb(iorb2) ) CYCLE
!                 nind2 = nindo(iorb2)
!                 Drx(:nind1,:nind2) = matmul(tmp(:nind1,:nind), transpose(cfX(:nind2,:nind,iorb2,iorb_rotated)))                 
!                 do ind1=1,nind1
!                    do ind2=1,nind2
!                       i1 = iSx(ind1,iorb1)
!                       i2 = iSx(ind2,iorb2)
!                       RT(i1,i2,icix) = RT(i1,i2,icix) + Drx(ind1,ind2)
!                    enddo
!                 enddo
!              enddo
!           enddo
!           deallocate( Dr )
!        enddo
!     ENDDO
!
!     do icix=1,ncix
!        cdim = cixdim(icix)
!        WRITE(6,*) 'icix=', icix
!        do i=1,cdim
!           do j=1,cdim
!              WRITE(6,'(2F14.7,2x)',advance='no') RT(i,j,icix)
!           enddo
!           WRITE(6,*)
!        enddo
!     end do
!  enddo
!END SUBROUTINE SymmetrizeTest2

