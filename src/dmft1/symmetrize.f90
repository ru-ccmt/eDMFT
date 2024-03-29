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
  use case, ONLY : natom, iatom
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
  USE structure, ONLY: BR1, rotij, tau, iz
  USE case, ONLY: natom, nl, ll, cix, iatom, isort, maxdim, ncix, crotloc
  use com, ONLY: iso
  USE defs,  ONLY: IMAG
  USE param, ONLY: lmax2
  IMPLICIT NONE
  INTEGER, intent(in)      :: ig ! The index of the group operation
  COMPLEX*16, intent(out)  :: RT(maxdim,maxdim,ncix)  ! Symmetrization matrix for this group operation
  ! other needed inputs
  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmax2+1), nindo(norbitals)
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
  COMPLEX*16 :: Drx(maxdim2,maxdim2)
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


     do lcase=1,nl(icase)
     
        local_axis_latom   = matmul(crotloc(:,:,lcase,icase),         matmul(BR1, matmul(rotij(:,:,latom),         BR1inv)))
        local_axis_rotated = matmul(crotloc(:,:,lcase,icase_rotated), matmul(BR1, matmul(rotij(:,:,latom_rotated), BR1inv)))
        call inv_3x3(local_axis_rotated, local_axis_rotated_inv)
     
        Tac = matmul(local_axis_latom, matmul(TR, local_axis_rotated_inv))
        Det = detx(Tac)
        if (Det < 0.d0) then
           Tac = -Tac
        endif

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
           tmp(:nind1,:nind) = matmul( conjg(cfX(:nind1,:nind,iorb1,iorb)), Dr(:nind,:nind) )
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

SUBROUTINE SymmetrizeLocalQuantities(Gmloc, Eimpm, Olapm, s_oo, sigma, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, noccur, maxdim2, nomega)
  USE structure,  ONLY: iord
  USE case, ONLY: maxdim, ncix, natom, maxsize, csize, Sigind
  USE defs,  ONLY: IMAG
  USE param, ONLY: lmax2
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: Gmloc(maxdim,maxdim,ncix,nomega), Eimpm(maxdim,maxdim,ncix), Olapm(maxdim,maxdim,ncix)
  COMPLEX*16, intent(inout):: s_oo(maxsize,ncix), sigma(maxsize,ncix,nomega)
  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmax2+1), nindo(norbitals), noccur(maxsize,ncix)
  INTEGER, intent(in)      :: norbitals, maxdim2, nomega
  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  ! locals
  INTEGER    :: icix, ig, ip, iq, iom, it, i, j
  INTEGER    :: cdim
  COMPLEX*16 :: gc(maxdim,maxdim)
  !
  INTEGER    :: ncix_unique
  INTEGER    :: size_cix_unique(ncix)
  INTEGER    :: cix_unique(ncix,ncix)
  COMPLEX*16, allocatable :: tmp(:,:)
  !
  COMPLEX*16, allocatable :: RT(:,:,:,:), sig(:,:)
  
  allocate( RT(maxdim,maxdim,ncix,iord) )
  allocate( sig(maxdim,maxdim) )
  
  do ig=1,iord
     Call CmpSymmetryTransformation(ig, RT(:,:,:,ig), cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2, .False.)
  enddo

  do icix=1,ncix
     cdim = cixdim(icix)
     do iom=1,nomega
        ! Symmetrize the local Green's function
        gc(:,:)=0.d0
        do ig=1,iord
           gc(:cdim,:cdim) = gc(:cdim,:cdim) + matmul(conjg(transpose(RT(:cdim,:cdim,icix,ig))), matmul(Gmloc(:cdim,:cdim,icix,iom), RT(:cdim,:cdim,icix,ig)))/iord
        end do
        Gmloc(:cdim,:cdim,icix,iom) = gc(:cdim,:cdim)

        ! Symmetrize the self-energy
        sig=0.d0
        DO ip=1,cdim
           do iq=1,cdim
              it = abs(Sigind(ip,iq,icix))
              if (it.ne.0) then
                 sig(ip,iq) = sigma(it,icix,iom)
              endif
           enddo
        ENDDO
        gc(:,:)=0.d0
        do ig=1,iord
           gc(:cdim,:cdim) = gc(:cdim,:cdim) + matmul(conjg(transpose(RT(:cdim,:cdim,icix,ig))), matmul(sig(:cdim,:cdim), RT(:cdim,:cdim,icix,ig)))/iord
        end do
        sigma(:,icix,iom)=0.d0
        DO ip=1,cdim
           do iq=1,cdim
              it = abs(Sigind(ip,iq,icix))
              if (it.ne.0) then
                 sigma(it,icix,iom) = sigma(it,icix,iom) + gc(ip,iq)
              endif
           enddo
        ENDDO
        DO it=1,csize(icix)
           sigma(it, icix, iom) = sigma(it,icix,iom)/noccur(it,icix)
        ENDDO
        
     enddo
     
     ! Symmetrize the impurity levels
     gc(:,:)=0.d0
     do ig=1,iord
        gc(:cdim,:cdim) = gc(:cdim,:cdim) + matmul(conjg(transpose(RT(:cdim,:cdim,icix,ig))), matmul(Eimpm(:cdim,:cdim,icix), RT(:cdim,:cdim,icix,ig)))/iord
     end do
     Eimpm(:cdim,:cdim,icix) = gc(:cdim,:cdim)
     ! Symmetrize the overlap
     gc(:,:)=0.d0
     do ig=1,iord
        gc(:cdim,:cdim) = gc(:cdim,:cdim) + matmul(conjg(transpose(RT(:cdim,:cdim,icix,ig))), matmul(Olapm(:cdim,:cdim,icix), RT(:cdim,:cdim,icix,ig)))/iord
     end do
     Olapm(:cdim,:cdim,icix) = gc(:cdim,:cdim)

     ! Symmetrize the sigma infinity s_oo
     sig=0.d0
     DO ip=1,cdim
        do iq=1,cdim
           it = abs(Sigind(ip,iq,icix))
           if (it.ne.0) then
              sig(ip,iq) = s_oo(it,icix)
           endif
        enddo
     ENDDO
     gc(:,:)=0.d0
     do ig=1,iord
        gc(:cdim,:cdim) = gc(:cdim,:cdim) + matmul(conjg(transpose(RT(:cdim,:cdim,icix,ig))), matmul(sig(:cdim,:cdim), RT(:cdim,:cdim,icix,ig)))/iord
     end do
     s_oo(:,icix)=0.d0
     DO ip=1,cdim
        do iq=1,cdim
           it = abs(Sigind(ip,iq,icix))
           if (it.ne.0) then
              s_oo(it,icix) = s_oo(it,icix) + gc(ip,iq)
           endif
        enddo
     ENDDO
     DO it=1,csize(icix)
        s_oo(it, icix) = s_oo(it,icix)/noccur(it,icix)
     ENDDO
     
  enddo
  
  deallocate( RT )

  ! Now symmetrize over equivalent cix
  call FindUniqueCix(ncix_unique, size_cix_unique, cix_unique)
  do i=1,ncix_unique
     if (size_cix_unique(i).eq.1) cycle ! no symmetrization necessary
     icix = cix_unique(i,1)
     cdim = cixdim(icix)
     allocate( tmp(cdim,cdim) )

     !Gmloc(:cdim,:cdim,icix,iom)
     do iom=1,nomega
        tmp=0.d0
        do j=1,size_cix_unique(i)
           icix = cix_unique(i,j)
           tmp(:,:) = tmp(:,:) + Gmloc(:cdim,:cdim,icix,iom)
        enddo
        tmp(:,:) = tmp(:,:)/size_cix_unique(i)
        do j=1,size_cix_unique(i)
           icix = cix_unique(i,j)
           Gmloc(:cdim,:cdim,icix,iom) = tmp(:,:)
        enddo
     enddo
     !Eimpm(:cdim,:cdim,icix)
     tmp=0.d0
     do j=1,size_cix_unique(i)
        icix = cix_unique(i,j)
        tmp(:,:) = tmp(:,:) + Eimpm(:cdim,:cdim,icix)
     enddo
     tmp(:,:) = tmp(:,:)/size_cix_unique(i)
     do j=1,size_cix_unique(i)
        icix = cix_unique(i,j)
        Eimpm(:cdim,:cdim,icix) = tmp(:,:)
     enddo
     !Olapm(:cdim,:cdim,icix)
     tmp=0.d0
     do j=1,size_cix_unique(i)
        icix = cix_unique(i,j)
        tmp(:,:) = tmp(:,:) + Olapm(:cdim,:cdim,icix)
     enddo
     tmp(:,:) = tmp(:,:)/size_cix_unique(i)
     do j=1,size_cix_unique(i)
        icix = cix_unique(i,j)
        Olapm(:cdim,:cdim,icix) = tmp(:,:)
     enddo
     
     deallocate( tmp )
  enddo

  
END SUBROUTINE SymmetrizeLocalQuantities

SUBROUTINE SymmetrizeLocalQuantities2(g_inf, g_ferm, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2, nomega)
  USE structure,  ONLY: iord
  USE case, ONLY: maxdim, ncix, natom
  USE defs,  ONLY: IMAG
  USE param, ONLY: lmax2
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: g_inf(maxdim,maxdim,ncix,nomega), g_ferm(maxdim,maxdim,ncix)
  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmax2+1), nindo(norbitals)
  INTEGER, intent(in)      :: norbitals, maxdim2, nomega
  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  ! locals
  INTEGER    :: icix, ig, iom, i, j
  INTEGER    :: cdim
  INTEGER    :: ncix_unique
  INTEGER    :: size_cix_unique(ncix)
  INTEGER    :: cix_unique(ncix,ncix)
  COMPLEX*16 :: gc(maxdim,maxdim)
  COMPLEX*16, allocatable :: RT(:,:,:,:), tmp(:,:)
  
  allocate( RT(maxdim,maxdim,ncix,iord) )
  
  do ig=1,iord
     Call CmpSymmetryTransformation(ig, RT(:,:,:,ig), cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2, .False.)
  enddo

  do icix=1,ncix
     cdim = cixdim(icix)
     do iom=1,nomega
        ! Symmetrize the local Green's function
        gc(:,:)=0.d0
        do ig=1,iord
           gc(:cdim,:cdim) = gc(:cdim,:cdim) + matmul(conjg(transpose(RT(:cdim,:cdim,icix,ig))), matmul(g_inf(:cdim,:cdim,icix,iom), RT(:cdim,:cdim,icix,ig)))/iord
        end do
        g_inf(:cdim,:cdim,icix,iom) = gc(:cdim,:cdim)
     enddo
     gc(:,:)=0.d0
     do ig=1,iord
        gc(:cdim,:cdim) = gc(:cdim,:cdim) + matmul(conjg(transpose(RT(:cdim,:cdim,icix,ig))), matmul(g_ferm(:cdim,:cdim,icix), RT(:cdim,:cdim,icix,ig)))/iord
     end do
     g_ferm(:cdim,:cdim,icix) = gc(:cdim,:cdim)
  enddo
  
  deallocate( RT )

  ! Now symmetrize over equivalent cix
  call FindUniqueCix(ncix_unique, size_cix_unique, cix_unique)
  do i=1,ncix_unique
     if (size_cix_unique(i).eq.1) cycle ! no symmetrization necessary
     icix = cix_unique(i,1)
     cdim = cixdim(icix)
     allocate( tmp(cdim,cdim) )

     !g_inf(:cdim,:cdim,icix,iom)
     do iom=1,nomega
        tmp=0.d0
        do j=1,size_cix_unique(i)
           icix = cix_unique(i,j)
           tmp(:,:) = tmp(:,:) + g_inf(:cdim,:cdim,icix,iom)
        enddo
        tmp(:,:) = tmp(:,:)/size_cix_unique(i)
        do j=1,size_cix_unique(i)
           icix = cix_unique(i,j)
           g_inf(:cdim,:cdim,icix,iom) = tmp(:,:)
        enddo
     enddo
     !g_ferm(:cdim,:cdim,icix)
     tmp=0.d0
     do j=1,size_cix_unique(i)
        icix = cix_unique(i,j)
        tmp(:,:) = tmp(:,:) + g_ferm(:cdim,:cdim,icix)
     enddo
     tmp(:,:) = tmp(:,:)/size_cix_unique(i)
     do j=1,size_cix_unique(i)
        icix = cix_unique(i,j)
        g_ferm(:cdim,:cdim,icix) = tmp(:,:)
     enddo
     
     deallocate( tmp )
  enddo
  
END SUBROUTINE SymmetrizeLocalQuantities2


SUBROUTINE SymmetrizeOverlap(Olapm, cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2)
  USE structure,  ONLY: iord
  USE case, ONLY: maxdim, ncix, natom
  USE defs,  ONLY: IMAG
  USE param, ONLY: lmax2
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: Olapm(maxdim,maxdim,ncix)
  INTEGER, intent(in)      :: cix_orb(norbitals), cixdim(ncix), iSx(maxdim2, norbitals), iorbital(natom,lmax2+1), nindo(norbitals)
  INTEGER, intent(in)      :: norbitals, maxdim2
  COMPLEX*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  ! locals
  INTEGER    :: icix, ig, i, j
  INTEGER    :: cdim
  INTEGER    :: ncix_unique
  INTEGER    :: size_cix_unique(ncix)
  INTEGER    :: cix_unique(ncix,ncix)
  !
  COMPLEX*16, allocatable :: tmp(:,:)
  COMPLEX*16, allocatable :: Olaps(:,:,:), RT(:,:,:)
  
  allocate( Olaps(maxdim,maxdim,ncix), RT(maxdim,maxdim,ncix) )

  Olaps(:,:,:)=0.d0
  do ig=1,iord
     Call CmpSymmetryTransformation(ig, RT(:,:,:), cfX, cix_orb, cixdim, iSx, iorbital, nindo, norbitals, maxdim2, .False.)
     do icix=1,ncix
        cdim = cixdim(icix)
        Olaps(:cdim,:cdim,icix) = Olaps(:cdim,:cdim,icix) + matmul(conjg(transpose(RT(:cdim,:cdim,icix))), matmul(Olapm(:cdim,:cdim,icix), RT(:cdim,:cdim,icix)))/iord
     enddo
  enddo

  ! Now symmetrize over equivalent cix
  call FindUniqueCix(ncix_unique, size_cix_unique, cix_unique)
  do i=1,ncix_unique
     if (size_cix_unique(i).eq.1) cycle ! no symmetrization necessary
     icix = cix_unique(i,1)
     cdim = cixdim(icix)
     allocate( tmp(cdim,cdim) )
     tmp=0.d0
     do j=1,size_cix_unique(i)
        icix = cix_unique(i,j)
        tmp(:,:) = tmp(:,:) + Olaps(:cdim,:cdim,icix)
     enddo
     tmp(:,:) = tmp(:,:)/size_cix_unique(i)
     do j=1,size_cix_unique(i)
        icix = cix_unique(i,j)
        Olaps(:cdim,:cdim,icix) = tmp(:,:)
     enddo
     deallocate( tmp )
  enddo
  
  Olapm(:,:,:) = Olaps(:,:,:)
  deallocate( RT, Olaps )
END SUBROUTINE SymmetrizeOverlap

SUBROUTINE FindUniqueCix(ncix_unique, size_cix_unique, cix_unique)
  use case, ONLY: natom, ncix, iatom, isort, cix
  IMPLICIT NONE
  INTEGER, intent(out) :: ncix_unique
  INTEGER, intent(out) :: size_cix_unique(ncix)
  INTEGER, intent(out) :: cix_unique(ncix,ncix)
  !
  INTEGER :: cix_2_jatom(ncix)
  INTEGER :: icix, jatom_unique
  INTEGER :: icase, jatom, latom

  do icix=1,ncix
     DO icase=1,natom
        latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
        jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )
        if (icix .eq. cix(icase,1)) then
           cix_2_jatom(icix) = jatom
           exit
        endif
     ENDDO
  enddo
  
  ncix_unique=1
  size_cix_unique(:)=0
  cix_unique(:,:)=0
  icix=1
  jatom_unique = cix_2_jatom(icix)
  do icix=1,ncix
     if ( cix_2_jatom(icix) .eq. jatom_unique ) then
        size_cix_unique(ncix_unique) = size_cix_unique(ncix_unique) + 1
        cix_unique( ncix_unique, size_cix_unique(ncix_unique) ) = icix
        cycle
     endif
     jatom_unique = cix_2_jatom(icix)
     ncix_unique = ncix_unique + 1
     size_cix_unique(ncix_unique) = size_cix_unique(ncix_unique) + 1
     cix_unique( ncix_unique, size_cix_unique(ncix_unique) ) = icix
  enddo

  !print *, 'icix=', icix-1, 'ncix_unique=', ncix_unique
  !print *, 'size_cix_unique: '
  !do i=1,ncix_unique
  !   print *, size_cix_unique(i)
  !enddo
  !print *, 'cix_unique: '
  !do i=1,ncix_unique
  !   do j=1,size_cix_unique(i)
  !      print *, i, j, cix_unique(i,j)
  !   enddo
  !enddo
  
END SUBROUTINE FindUniqueCix
