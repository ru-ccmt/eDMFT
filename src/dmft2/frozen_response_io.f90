! Optional interfaces for the frozen-self-energy linear-response diagnostic.
!
! These routines do not alter a normal dmft2 calculation.  The Python response
! driver differentiates the central lattice resolvent, including the moving
! finite band window and the global inverse-square-root normalization.  Its
! Hermitian, SIGNED density matrix is then reconstructed using the same LAPW
! radial/angular/FFT machinery as an ordinary density.  No physical occupation
! cutoff may be applied: positive and negative response eigenvalues cancel.
!
! The input includes the central eigenvectors, not just eigenvalues.  Equal
! eigenvalues do not identify a gauge inside a degenerate subspace.  Refusing
! a different vector file prevents a silently gauge-dependent charge response.
SUBROUTINE ReadFrozenResponseDensity(prefix,ikp,n0,iso,nmat,nume, &
     first,last,E,As,rho)
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  CHARACTER(*),intent(in)::prefix
  INTEGER,intent(in)::ikp,n0,iso,nmat,nume,first,last
  REAL*8,intent(in)::E(nume)
  COMPLEX*16,intent(in)::As(nmat,nume,iso)
  COMPLEX*16,intent(out)::rho(last-first+1,last-first+1)
  CHARACTER(8)::magic
  CHARACTER(2048)::path
  INTEGER::u,ios,h(7),n
  REAL*8,allocatable::er(:)
  COMPLEX*16,allocatable::ar(:,:,:)
  n=last-first+1
  write(path,'(A,A,I8.8,A)') trim(prefix),'.k',ikp,'.bin'
  open(newunit=u,file=trim(path),status='old',access='stream', &
       form='unformatted',iostat=ios)
  if(ios.ne.0) error stop 'response density: cannot open input'
  read(u,iostat=ios) magic,h
  if(ios.ne.0) error stop 'response density: truncated header'
  if(magic.ne.'EDMRHO01') error stop 'response density: invalid magic'
  if(any(h.ne.(/1,ikp,n0,iso,n,first,last/))) &
       error stop 'response density: incompatible basis dimensions'
  allocate(er(n),ar(n0,n,iso))
  read(u,iostat=ios) er,ar,rho
  close(u)
  if(ios.ne.0) error stop 'response density: truncated matrices'
  if(maxval(abs(er-E(first:last))).gt.1.d-10) &
       error stop 'response density: eigenvalues differ from reference'
  if(maxval(abs(ar-As(1:n0,first:last,:))).gt.1.d-10) &
       error stop 'response density: eigenvector gauge differs from reference'
  if(.not.all(ieee_is_finite(dble(rho)))) error stop 'nonfinite response density'
  if(.not.all(ieee_is_finite(aimag(rho)))) error stop 'nonfinite response density'
  if(maxval(abs(rho-conjg(transpose(rho)))).gt.1.d-10) &
       error stop 'non-Hermitian response density'
  rho=0.5d0*(rho+conjg(transpose(rho)))
  deallocate(er,ar)
END SUBROUTINE ReadFrozenResponseDensity

! Export the actual quadrature, self-energy, and normalization symmetry map.
! A Python reimplementation must not guess the orbital symmetry from Sigind:
! overlap symmetry is a group average, followed by equivalent-site averaging.
! Exporting its action on real and imaginary matrix units also captures the
! final Hermitian projection, which is real-linear rather than complex-linear.
SUBROUTINE ExportFrozenResponseMetadata(prefix,omega,weights,npomega, &
     sigma,sinf,nomega,maxsize,mu,temp,SOlapm,cfX,cix_orb,cixdim, &
     iSx,iorbital,nindo,norbitals,maxdim2)
  USE dmfts,ONLY:maxdim,ncix,natom,lmaxp,Sigind
  USE com_mpi,ONLY:myrank,master
  IMPLICIT NONE
  CHARACTER(*),intent(in)::prefix
  INTEGER,intent(in)::npomega,nomega,maxsize,norbitals,maxdim2
  REAL*8,intent(in)::omega(nomega),weights(npomega),mu,temp
  COMPLEX*16,intent(in)::sigma(maxsize,ncix,nomega),sinf(maxsize,ncix)
  COMPLEX*16,intent(in)::SOlapm(maxdim,maxdim,ncix)
  COMPLEX*16,intent(in)::cfX(maxdim2,maxdim2,norbitals,norbitals)
  INTEGER,intent(in)::cix_orb(norbitals),cixdim(ncix),iSx(maxdim2,norbitals)
  INTEGER,intent(in)::iorbital(natom,lmaxp+1),nindo(norbitals)
  COMPLEX*16,allocatable::local(:,:,:),unit(:,:,:)
  INTEGER::u,iw,ic,i,j,part
  if(myrank.ne.master.or.len_trim(prefix).eq.0)return
  allocate(local(maxdim,maxdim,ncix),unit(maxdim,maxdim,ncix))
  open(newunit=u,file=trim(prefix)//'.meta',status='replace', &
       access='stream',form='unformatted')
  write(u)'EDMRSP01',2,maxdim,ncix,maxdim2,norbitals,npomega
  write(u)mu,temp,omega(:npomega),weights
  write(u)cixdim,cix_orb,nindo,iSx,Sigind,SOlapm
  call BuildLocalSigmaMatrix(local,sinf,Sigind,cixdim,maxdim,maxsize,ncix)
  write(u)local
  do iw=1,npomega
     call BuildLocalSigmaMatrix(local,sigma(:,:,iw),Sigind,cixdim, &
          maxdim,maxsize,ncix)
     write(u)local
  enddo
  do part=1,2
     do ic=1,ncix
        do j=1,maxdim
           do i=1,maxdim
              unit=0.d0
              if(i.le.cixdim(ic).and.j.le.cixdim(ic))then
                 if(part.eq.1)unit(i,j,ic)=1.d0
                 if(part.eq.2)unit(i,j,ic)=(0.d0,1.d0)
              endif
              call SymmetrizeOverlap(unit,cfX,cix_orb,cixdim,iSx, &
                   iorbital,nindo,norbitals,maxdim2)
              write(u)unit
           enddo
        enddo
     enddo
  enddo
  close(u)
  deallocate(local,unit)
  write(6,'(A,1X,A)')':FROZEN_RESPONSE_METADATA:',trim(prefix)//'.meta'
END SUBROUTINE ExportFrozenResponseMetadata

! Central k-point data for a standalone response calculation.  The saved U
! contains all stored rows of the raw overlap AFTER normalization, before
! truncation to W.  Python recovers raw P through the saved invertible X.
! Only W is embedded in the central Dyson equation.  Outside rows are needed
! for dPi_W and must never be embedded as an enlarged impurity window.
SUBROUTINE ExportFrozenResponseK(prefix,ikp,n0,iso,nmat,nume,nstored, &
     first,last,wfirst,wlast,maxdim2,norbitals,E,As,k3,k3lo,nnlo,U,wkp)
  IMPLICIT NONE
  CHARACTER(*),intent(in)::prefix
  INTEGER,intent(in)::ikp,n0,iso,nmat,nume,nstored,first,last,wfirst,wlast
  INTEGER,intent(in)::maxdim2,norbitals,nnlo,k3(3,nmat),k3lo(3,nnlo)
  REAL*8,intent(in)::E(nume),wkp
  COMPLEX*16,intent(in)::As(nmat,nume,iso),U(nstored,maxdim2,norbitals)
  INTEGER::fh
  CHARACTER(2048)::path
  if(len_trim(prefix).eq.0)return
  write(path,'(A,A,I8.8,A)')trim(prefix),'.k',ikp,'.bin'
  open(newunit=fh,file=trim(path),status='replace',access='stream',form='unformatted')
  write(fh)'EDMRPK01',1,ikp,n0,iso,nstored,first,last,wfirst,wlast, &
       maxdim2,norbitals,nnlo
  write(fh)wkp,E(first:last),As(1:n0,first:last,:),U
  write(fh)k3(:,1:n0-nnlo),k3lo
  close(fh)
END SUBROUTINE ExportFrozenResponseK
