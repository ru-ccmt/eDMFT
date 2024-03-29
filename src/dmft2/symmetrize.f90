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

