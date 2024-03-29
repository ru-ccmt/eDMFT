! @Copyright 2007 Kristjan Haule
! 

SUBROUTINE Create_Atom_Arrays(csort, maxucase, isort, iatom, nat, natm, natom)
  USE com_mpi, ONLY: Qprint
  IMPLICIT NONE
  INTEGER, intent(out) :: csort(nat), maxucase
  INTEGER, intent(in)  :: isort(natm), iatom(natom), nat, natm, natom
  !local variables
  INTEGER :: icase, iucase
  csort=0
  maxucase=1
  do icase=1,natom
     if (csort(isort(iatom(icase))).EQ.0) then
        csort(isort(iatom(icase))) = maxucase
        maxucase = maxucase + 1
     endif
  enddo
  maxucase = maxucase-1
  
  if (Qprint) then
     WRITE(6,*)
     WRITE(6,*) '********** Information about the atomic index arrays ************'
     WRITE(6,*) 'iatom=', iatom
     WRITE(6,*) 'isort=', isort
     WRITE(6,*) 'csort=', csort
     WRITE(6,*) 'maxucase=', maxucase
  
     DO icase=1,natom
        iucase = csort(isort(iatom(icase)))
        WRITE(6, '(A,6I6,1x)') 'icase, iucase, jatom, latom=', icase, iucase, isort(iatom(icase)), iatom(icase)
     ENDDO
  endif
  
END SUBROUTINE Create_Atom_Arrays

SUBROUTINE Create_Orbital_Arrays(iorbital, norbitals, maxdim2, nl, ll, cix, natom, lmax2, iso)
  USE com_mpi, ONLY: Qprint
  IMPLICIT NONE
  INTEGER, intent(out) :: iorbital(natom,lmax2+1), maxdim2, norbitals
  INTEGER, intent(in)  :: nl(natom), ll(natom,4), cix(natom,4), natom, lmax2, iso
  !----- locals
  INTEGER :: icase, lcase, l1, nind, icix
  !----------- Find index to atom/L named iorbital(icase,lcase) -------------------------------!
  iorbital=0
  maxdim2=0   !-- maximum size of the matrix over all atoms and l's requested in the input ----!
              !---- do not confuse with maxdim, which is maximum size for correlated orbitals -!
  norbitals=0 !-- number of atom/l cases, called orbitals -------------------------------------!
  do icase=1,natom
     do lcase=1,nl(icase)
        norbitals = norbitals+1
        iorbital(icase,lcase)=norbitals  !-- index to the orbital number ----!
        icix = cix(icase,lcase)
        
        l1=ll(icase,lcase)
        nind=(2*l1+1)*iso
        maxdim2 = max(maxdim2, nind )

     enddo
  enddo  

  if (Qprint) then
     WRITE(6,*)
     WRITE(6,*) '*****  Arrays for correlated blocks ******'
     WRITE(6,*) 'norbitals=', norbitals
     WRITE(6,*) 'maxdim2=', maxdim2
     DO icase=1,natom
        DO lcase=1,nl(icase)
           WRITE(6,'(A,I2,A,I2,A,I4)') 'iorbital(',icase,',',lcase,')=', iorbital(icase,lcase)
        ENDDO
     ENDDO
     DO icase=1,natom
        DO lcase=1,nl(icase)
           WRITE(6,'(A,I2,A,I2,A,I4)') 'cix(',icase,',',lcase,')=', cix(icase,lcase)
        ENDDO
     ENDDO
  endif
  !----------- Find index to atom/L named iorbital(icase,lcase) --------
END SUBROUTINE Create_Orbital_Arrays
  

SUBROUTINE Create_Other_Arrays(cixdim, iSx, noccur, nindo, cix_orb, cfX, CF, nl, ll, cix, iorbital, csize, csizes, Sigind, iso, natom, maxdim, lmax2, ncix, maxsize, norbitals, maxdim2)
  USE com_mpi, ONLY: Qprint
  IMPLICIT NONE
  INTEGER, intent(out)    :: cixdim(ncix), iSx(maxdim2, norbitals), noccur(maxsize,ncix), nindo(norbitals), cix_orb(norbitals), csizes(ncix)
  COMPLEX*16, intent(out) :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  COMPLEX*16, intent(in)  :: CF(maxdim,maxdim,ncix)
  INTEGER, intent(in)     :: nl(natom), ll(natom,4), cix(natom,4), iorbital(natom,lmax2+1), csize(ncix)
  INTEGER, intent(in)     :: Sigind(maxdim,maxdim,ncix)
  INTEGER, intent(in)     :: natom, maxdim, lmax2, ncix, iso, maxsize, norbitals, maxdim2
  !----- locals
  INTEGER :: icase, lcase, iorb, icix, ip, iq, it, l1, nind1, nind2, ip1, iorb1, iorb2, ind1, ind2, wcsize

  cixdim=0
  iSx=0
  DO icase=1,natom            
     do lcase=1,nl(icase)     
        icix = cix(icase,lcase)
        l1 = ll(icase,lcase)
        iorb = iorbital(icase,lcase)
        nind1 = (2*l1+1)*iso
        nindo(iorb) = nind1

        cix_orb(iorb) = icix
        
        if ( icix.EQ.0 ) CYCLE
        do ip1=1,nind1
           cixdim(icix) = cixdim(icix) + 1
           iSx(ip1,iorb) = cixdim(icix)
           if (Qprint) WRITE(6,'(A,7I5)') 'icase,lcase,icix,iorb,nind1,ip1,iSx=', icase, lcase, icix, iorb, nind1,ip1,cixdim(icix)
        enddo
     enddo
  ENDDO

  !print *, 'DEBUGGING=', cixdim, Sigind

  noccur=0
  DO icix=1,ncix
     wcsize=0
     DO ip=1,cixdim(icix)
        DO iq=1,cixdim(icix)
           it = abs(Sigind(ip,iq,icix))
           if (Sigind(ip,iq,icix).gt.0) then
              if (noccur(it,icix).eq.0) wcsize=wcsize+1
           endif
           if (it.gt.0)  noccur(it,icix) = noccur(it,icix) + 1
        ENDDO
     ENDDO
     csizes(icix)=wcsize
  ENDDO

  DO iorb1=1,norbitals
     nind1 = nindo(iorb1)
     if ( cix_orb(iorb1).EQ.0 ) CYCLE
     DO iorb2=1,norbitals
        nind2 = nindo(iorb2)
        if ( cix_orb(iorb1).NE.cix_orb(iorb2) ) CYCLE
        icix = cix_orb(iorb1)
        do ind1=1,nind1
           do ind2=1,nind2
              !WRITE(*,'(I2,I2,I3,I3,I3,I3,2x,I3,2f14.5)') iorb1, iorb2, ind1, ind2, iSx(ind1,iorb1), iSx(ind2,iorb2), icix, CF( iSx(ind1,iorb1),iSx(ind2,iorb2), icix )
              cfX(ind1,ind2,iorb1,iorb2) = CF( iSx(ind1,iorb1),iSx(ind2,iorb2), icix )
           enddo
        enddo
     ENDDO
  ENDDO
  
  
  if (Qprint) then
     do icix=1,ncix
        WRITE(6,'(A,I2,A,I2)') 'cixdim(', icix, ')=', cixdim(icix)
        DO it=1,csize(icix)
           WRITE(6,'(A,I3,A,I2,A,I4)') 'noccur(',it,',',icix,')=', noccur(it,icix)
        ENDDO
     enddo
     do iorb=1,norbitals
        do ip1=1,nindo(iorb)
           WRITE(6,'(A,I2,A,I2,A,I3)') 'iSx(', ip1, ',', iorb, ')=', iSx(ip1,iorb)
        enddo
     enddo
     do iorb=1,norbitals
        WRITE(6,'(A,I2,A,I3)') 'nindo(', iorb, ')=', nindo(iorb)
     enddo
     do iorb1=1,norbitals
        do iorb2=1,norbitals
           if ( cix_orb(iorb1).EQ.0 ) CYCLE
           if ( cix_orb(iorb1).NE.cix_orb(iorb2) ) CYCLE
           nind1 = nindo(iorb1)
           nind2 = nindo(iorb2)
           WRITE(6,'(A,I2,A,I2)') 'CF for iorb1=', iorb1, ' iorb2=', iorb2
           do ind1=1,nind1
              do ind2=1,nind2
                 WRITE(6,'(f14.8,1x,f14.8,3x)', advance='no') real(cfX(ind1,ind2,iorb1,iorb2)), aimag(cfX(ind1,ind2,iorb1,iorb2))
              enddo
              WRITE(6,*)
           enddo
        enddo
     enddo
  endif
END SUBROUTINE Create_Other_Arrays


SUBROUTINE PrintSomeArrays(filename, nat, iso, norbitals, ncix, natom, nkpt, nmat, nume, Qcomplex, lmax2, maxdim2, maxdim, maxsize, nindo, cixdim, nl, ll, cix, iorbital, csize, iSx, Sigind, EF, VOL)
  use structure, ONLY : POS, rotij, tauij
  use case, ONLY : iatom, ifirst, shft
  IMPLICIT NONE
  CHARACTER*100, intent(in) :: filename
  INTEGER, intent(in) :: nat, iso, norbitals, ncix, natom, nkpt, nmat, nume, lmax2, maxdim2, maxdim, maxsize
  LOGICAL, intent(in) :: Qcomplex
  INTEGER, intent(in) :: nindo(norbitals), cixdim(ncix), nl(natom), ll(natom,4), cix(natom,4), iorbital(natom,lmax2+1)
  INTEGER, intent(in) :: csize(ncix), iSx(maxdim2,norbitals), Sigind(maxdim,maxdim,ncix)
  REAL*8, intent(in)  :: EF, VOL
  ! locals
  INTEGER :: fh, icase, lcase, iorb, ip1, icix, ip, iq, latom, lfirst
  REAL*8 :: POST(3)
  fh=996
  open(fh,file=TRIM(filename),status='unknown', form='formatted')
  WRITE(fh, *) 'nat, iso, norbitals, ncix, natom'
  WRITE(fh, *) nat, iso, norbitals, ncix, natom
  WRITE(fh, *) 'nkpt, nmat, nume'
  WRITE(fh, *) nkpt, nmat, nume
  WRITE(fh, *) 'Qcomplex'
  WRITE(fh, *) Qcomplex
  WRITE(fh, *) 'lmax2, maxdim2, maxdim, maxsize'
  WRITE(fh, *) lmax2, maxdim2, maxdim, maxsize
  WRITE(fh, *) 'nindo'
  WRITE(fh, *) (nindo(iorb),iorb=1,norbitals)
  WRITE(fh, *) 'cidim'
  WRITE(fh, *) (cixdim(icix),icix=1,ncix)
  WRITE(fh, *) 'nl'
  WRITE(fh, *) (nl(icase), icase=1,natom)
  WRITE(fh, *) 'll, iorbital, cix, nind'
  do icase=1,natom
     do lcase=1,nl(icase)
        WRITE(fh, *) ll(icase,lcase)
        WRITE(fh, *) iorbital(icase,lcase)
        WRITE(fh, *) cix(icase,lcase)
        WRITE(fh, *) (2*ll(icase,lcase)+1)*iso          ! nind
     enddo
  enddo
  
  WRITE(fh, *) 'csize'
  WRITE(fh, *) (csize(icix), icix=1,ncix)
  
  WRITE(fh, *) 'iSx'
  do iorb=1,norbitals
     do ip1=1,nindo(iorb)
        WRITE(fh, *) iSx(ip1,iorb)
     enddo
  enddo
  WRITE(fh, *) 'Sigind'
  DO icix=1,ncix
     DO ip=1,cixdim(icix)
        DO iq=1,cixdim(icix)
           WRITE(fh, *) abs(Sigind(ip,iq,icix))
        ENDDO
     ENDDO
  ENDDO
  WRITE(fh, *) 'EF, VOL'
  WRITE(fh, *) EF, VOL


  DO icase=1,natom  !--------------- over all atoms requested in the input ------------------------!
     latom = iatom(icase) 
     lfirst = ifirst(latom)
     do lcase=1,nl(icase)
        icix = cix(icase,lcase)
        if (icix>0) then
           POST = matmul(POS(:,lfirst), ROTIJ(:,:,latom))
           POST = POST + TAUIJ(:,latom) + shft(:,latom)
           WRITE(fh,*) (POST(ip),ip=1,3)
        endif
     enddo
  ENDDO
  
END SUBROUTINE PrintSomeArrays

SUBROUTINE CompressSigmaTransformation2(STrans, DMFTU, Sigind, iSx, cix, iorbital, ll, nl, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: STrans(maxsize,ncix,nbands,nbands)
  COMPLEX*16, intent(in)  :: DMFTU(nbands,maxdim2,norbitals)
  INTEGER, intent(in)     :: iSx(maxdim2, norbitals)
  INTEGER, intent(in)     :: Sigind(maxdim,maxdim,ncix), cix(natom,4)
  INTEGER, intent(in)     :: iorbital(natom,lmax2+1), ll(natom,4), nl(natom)
  INTEGER, intent(in)     :: natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize
  !----- locals
  INTEGER    :: i, j, icase, jcase, l1case, l2case, iorb1, iorb2, it, ind1, nind1, ind2, nind2, icix, l1, l2
  STrans=0
  DO icase=1,natom      
     do l1case=1,nl(icase) 
        icix = cix(icase,l1case)
        if ( icix.EQ.0 ) CYCLE
        l1 = ll(icase,l1case)
        nind1 = (2*l1+1)*iso
        iorb1 = iorbital(icase,l1case)
        DO jcase=1,natom
           do l2case=1,nl(jcase)
              if ( cix(jcase,l2case).NE.icix ) CYCLE
              l2 = ll(jcase,l2case)
              nind2 = (2*l2+1)*iso
              iorb2 = iorbital(jcase,l2case)
              do ind1=1,nind1
                 do ind2=1,nind2
                    it = abs(Sigind( iSx(ind1,iorb1),iSx(ind2,iorb2), icix ))
                    
                    if (it.gt.0) then
                       !WRITE(*,'(A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2)') 'icase=',icase, 'jcase=',jcase, 'l1case=',l1case, 'l2case=',l2case, 'iorb1=',iorb1, 'iorb2=',iorb2, 'ind1=',ind1, 'ind2=',ind2, 'iSx1=', iSx(ind1,iorb1), 'iSx2=', iSx(ind2,iorb2), 'it=',it
                       
                       DO i=1,nbands                     ! over bands-1
                          do j=1,nbands                  ! over bands-2
                             STrans(it,icix,i,j) = STrans(it,icix,i,j) + conjg(DMFTU(j,ind2,iorb2))*DMFTU(i,ind1,iorb1)
                          enddo
                       ENDDO
                       
                    endif
                 enddo
              enddo
           enddo
        ENDDO
     enddo
  ENDDO
END SUBROUTINE CompressSigmaTransformation2


SUBROUTINE AddSigma_optimized2(gij, sigma, STrans, csize, sign, nbands, ncix, maxsize)
  !-- Add's self-energy to inverse of the Green's function in band representation
  !-- 50% of all time is spend in this subroutine, the rest 50 % in CmpGk
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: gij(nbands,nbands)
  COMPLEX*16, intent(in)   :: sigma(maxsize,ncix)
  COMPLEX*16, intent(in)  :: STrans(maxsize,ncix,nbands,nbands)
  INTEGER, intent(in)      :: csize(ncix), sign, nbands, ncix, maxsize
  !----- locals
  COMPLEX*16 :: Sigmaij
  INTEGER    :: i, j, icix, it
  DO i=1,nbands                     ! over bands-1
     do j=1,nbands                  ! over bands-2
        Sigmaij=0
        DO icix=1,ncix
           do it=1,csize(icix)
              Sigmaij = Sigmaij + STrans(it,icix,i,j)*sigma(it,icix)
           enddo
        ENDDO
        gij(i,j) = gij(i,j) + sign*Sigmaij
     enddo
  ENDDO
END SUBROUTINE AddSigma_optimized2


SUBROUTINE CompressGcTransformation4(GTrans, DMFTrans, nbands, nindo, norbitals, maxdim2)
  IMPLICIT NONE
  COMPLEX*16, intent(out)  :: GTrans(nbands,nbands,norbitals)
  COMPLEX*16, intent(in)   :: DMFTrans(nbands,nbands,maxdim2,maxdim2,norbitals)
  INTEGER,    intent(in)   :: nbands, nindo(norbitals), norbitals, maxdim2
  !----- locals
  COMPLEX*16 :: csum
  INTEGER    :: i, j, iorb, ind1, nind
  DO i=1,nbands                     ! over bands-1
     do j=1,nbands                  ! over bands-2
        DO iorb=1,norbitals
           nind = nindo(iorb)

           csum=0
           do ind1=1,nind
              csum = csum + DMFTrans(i,j,ind1,ind1,iorb)
           enddo
           GTrans(i,j,iorb) = csum
           
        ENDDO
     enddo
  ENDDO
END SUBROUTINE CompressGcTransformation4

SUBROUTINE CmpGkc2(gmk, gij, DMFTU, iSx, iorbital, ll, nl, cix, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: gmk(maxdim,maxdim,ncix)
  COMPLEX*16, intent(in) :: gij(nbands,nbands)
  COMPLEX*16, intent(in) :: DMFTU(nbands,maxdim2,norbitals)
  INTEGER, intent(in)     :: iSx(maxdim2, norbitals)
  INTEGER, intent(in)     :: iorbital(natom,lmax2+1), ll(natom,4), nl(natom), cix(natom,4)
  INTEGER, intent(in)     :: natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands
  ! locals
  COMPLEX*16, allocatable :: tgk(:,:), tmp(:,:)
  INTEGER    :: icase, jcase, l1case, l2case, icix, l1, l2, nind1, nind2, iorb1, iorb2, ind1, ind2
  gmk=0
  allocate( tmp(maxdim2,nbands), tgk(maxdim2,maxdim2) )
  DO icase=1,natom      
     do l1case=1,nl(icase) 
        icix = cix(icase,l1case)
        if ( icix.EQ.0 ) CYCLE
        l1 = ll(icase,l1case)
        nind1 = (2*l1+1)*iso
        iorb1 = iorbital(icase,l1case)
        DO jcase=1,natom
           do l2case=1,nl(jcase)
              if ( cix(jcase,l2case).NE.icix ) CYCLE
              l2 = ll(jcase,l2case)
              nind2 = (2*l2+1)*iso
              iorb2 = iorbital(jcase,l2case)
              !print *, 'iorb1,iorb2=', iorb1, iorb2, maxdim2, nind1, nind2, nbands, maxdim2
              
              call zgemm('C','N', nind1, nbands, nbands, (1.d0,0.d0), DMFTU(:,:,iorb1), nbands, gij(:,:),nbands, (0.d0,0.d0), tmp(:,:),maxdim2)
              call zgemm('N','N', nind1, nind2, nbands, (1.d0,0.d0), tmp,maxdim2, DMFTU(:,:,iorb2),nbands, (0.d0,0.d0), tgk,maxdim2)
              
              !print *, 'tgk1=', tgk(:nind1,:nind2)
              !tgk(:nind1,:nind2) = matmul(matmul(transpose(conjg(DMFTU(:,:nind1,iorb1))), gij),DMFTU(:,:nind2,iorb2))
              !print *, 'tgk2=', tgk(:nind1,:nind2)

              do ind1=1,nind1
                 do ind2=1,nind2
                    gmk( iSx(ind1,iorb1), iSx(ind2,iorb2), icix ) = tgk(ind1,ind2)
                 enddo
              enddo
           enddo
        ENDDO
     enddo
  ENDDO
  deallocate( tmp, tgk )
END SUBROUTINE CmpGkc2


SUBROUTINE CmpGknc(gk, gij, GTrans, nindo, norbitals, ncix, nbands)
  ! Creates local green's function
  ! 50% of all time is spend in this subroutine, the rest 50 % in AddSigma
  IMPLICIT NONE
  COMPLEX*16, intent(out):: gk(norbitals)
  COMPLEX*16, intent(in) :: gij(nbands,nbands)
  COMPLEX*16, intent(in) :: GTrans(nbands,nbands,norbitals)
  INTEGER, intent(in)    :: nindo(norbitals), ncix, nbands, norbitals
  !------ locals
  INTEGER    :: i, j, nind, iorb
  COMPLEX*16 :: gc

  gk=0
  DO iorb=1,norbitals
     nind = nindo(iorb)
     gc=0
     if (ncix.gt.0) then    ! off-diagonal g in band representation
        DO i=1,nbands                      ! over bands-1
           do j=1,nbands                   ! over bands-2
              gc = gc + gij(i,j)*GTrans(i,j,iorb)
           enddo
        ENDDO
     else ! for non-correlated case, green's function is diagonal in band representation
        DO i=1,nbands      
           gc = gc + gij(i,i)*GTrans(i,i,iorb)
        ENDDO
     endif
     gk(iorb)=gc !/nind
  ENDDO
  
END SUBROUTINE CmpGknc



SUBROUTINE GetImpurityLevels2(Eimpmk, Olapmk, DMFTU, STrans, E, EF, s_oo, nemin, nume, iSx, iorbital, ll, nl, csize, cix, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Eimpmk(maxdim,maxdim,ncix)
  COMPLEX*16, intent(out) :: Olapmk(maxdim,maxdim,ncix)
  COMPLEX*16, intent(out) :: DMFTU(nbands,maxdim2,norbitals), STrans(maxsize,ncix,nbands,nbands)
  REAL*8,     intent(in)  :: E(nume), EF                         ! eigenvalues, chemical potential
  COMPLEX*16, intent(in)  :: s_oo(maxsize,ncix)
  INTEGER, intent(in)     :: nemin, nume, iSx(maxdim2, norbitals)
  INTEGER, intent(in)     :: iorbital(natom,lmax2+1), ll(natom,4), nl(natom), csize(ncix), cix(natom,4)
  INTEGER, intent(in)     :: natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize
  !--------- locals
  COMPLEX*16, allocatable :: tgk(:,:), tmp(:,:), gij(:,:), olp(:,:)
  INTEGER    :: i, icase, jcase, l1case, l2case, icix, l1, l2, nind1, nind2, iorb1, iorb2, ind1, ind2

  
  allocate( tmp(maxdim2,nbands), tgk(maxdim2, maxdim2), gij(nbands,nbands), olp(maxdim2,maxdim2) )

  gij=0
  DO i=1,nbands
     gij(i,i) = (E(i+nemin-1)-EF)
  ENDDO  
  
  CALL AddSigma_optimized2(gij, s_oo, STrans, csize, 1, nbands, ncix, maxsize)
  
  Eimpmk=0
  Olapmk=0
  DO icase=1,natom      
     do l1case=1,nl(icase) 
        icix = cix(icase,l1case)
        if ( icix.EQ.0 ) CYCLE
        l1 = ll(icase,l1case)
        nind1 = (2*l1+1)*iso
        iorb1 = iorbital(icase,l1case)
        DO jcase=1,natom
           do l2case=1,nl(jcase)
              if ( cix(jcase,l2case).NE.icix ) CYCLE
              l2 = ll(jcase,l2case)
              nind2 = (2*l2+1)*iso
              iorb2 = iorbital(jcase,l2case)

              call zgemm('C','N', nind1, nind2, nbands, (1.d0,0.d0), DMFTU(:,:,iorb1),nbands, DMFTU(:,:,iorb2),nbands, (0.d0,0.d0), olp,maxdim2)
              call zgemm('C','N', nind1, nbands, nbands, (1.d0,0.d0), DMFTU(:,:,iorb1), nbands, gij(:,:),nbands, (0.d0,0.d0), tmp(:,:),maxdim2)
              call zgemm('N','N', nind1, nind2,  nbands, (1.d0,0.d0), tmp,maxdim2, DMFTU(:,:,iorb2), nbands, (0.d0,0.d0), tgk,maxdim2)
              
              do ind1=1,nind1
                 do ind2=1,nind2
                    Eimpmk( iSx(ind1,iorb1), iSx(ind2,iorb2), icix ) = tgk(ind1,ind2)
                    Olapmk( iSx(ind1,iorb1), iSx(ind2,iorb2), icix ) = olp(ind1,ind2)
                 enddo
              enddo
              
           enddo
        ENDDO
     enddo
  ENDDO
  deallocate( gij, tmp, tgk, olp )
END SUBROUTINE GetImpurityLevels2

SUBROUTINE CountFrequencyShape(n0_om, nom_max, beta, omega, nomega)
  IMPLICIT NONE
  INTEGER, intent(out):: n0_om, nom_max
  REAL*8,  intent(out):: beta
  REAL*8,  intent(in) :: omega(nomega)
  INTEGER, intent(in) :: nomega
  !
  INTEGER :: iom
  REAL*8  :: PI
  PI=ACOS(-1.0D0)
  beta = PI/omega(1)
  do iom=1,nomega
     if (abs((omega(iom)*beta/PI+1)/2. - iom).GT.1e-1) EXIT
  enddo
  n0_om = iom-1
  if (n0_om.EQ.0) n0_om=1
  nom_max = nint((omega(nomega)*beta/PI+1)/2.)

END SUBROUTINE CountFrequencyShape


SUBROUTINE GetDelta2(Deltac, Glc, Eimpc, Olapc, logGloc, matsubara, omega, sigma, s_oo, gamma, gmloc, Eimpm, Olapm, Sigind, csize, cixdim, noccur, ncix, maxsize, maxdim, nomega, projector,ComputeLogGloc)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: Deltac(maxsize, ncix,nomega), Glc(maxsize,ncix,nomega)
  REAL*8, intent(out)    :: Eimpc(maxsize, ncix), Olapc(maxsize,ncix)
  REAL*8, intent(out)    :: logGloc(ncix)
  LOGICAL, intent(in)    :: matsubara
  REAL*8,  intent(in)    :: omega(nomega)
  COMPLEX*16, intent(in) :: sigma(maxsize,ncix,nomega), s_oo(maxsize,ncix)                                                                                                   
  REAL*8, intent(in)     :: gamma
  COMPLEX*16, intent(in) :: gmloc(maxdim,maxdim,ncix,nomega)
  COMPLEX*16, intent(inout) :: Eimpm(maxdim,maxdim,ncix), Olapm(maxdim,maxdim,ncix)
  INTEGER, intent(in)    :: Sigind(maxdim,maxdim,ncix), csize(ncix), cixdim(ncix)
  INTEGER, intent(in)    :: noccur(maxsize,ncix)
  INTEGER, intent(in)    :: ncix, maxsize, maxdim, nomega, projector
  LOGICAL, intent(in)    :: ComputeLogGloc
  !----- externals
  REAL*8 ::  FreeE0
  LOGICAL :: Off_diagonalc
  !----- locals
  real*8, PARAMETER       :: Ry2eV = 13.60569193
  COMPLEX*16, allocatable:: gcx(:,:), Eimpx(:,:), deltx(:,:), sigcx(:,:), sig_oo(:,:), gxc(:), Eim0x(:,:), olapx(:,:), tmp(:,:)
  COMPLEX*16 :: xomega, zsum
  INTEGER    :: icix, ip, iq, it, iom, cixdm, cixdms, n0_om, nom_max, n0, left
  complex*16 :: IMAG
  REAL*8 :: PI, beta, omw, ywr, ywi, ypwr, ypwi, yppwr, yppwi
  REAL*8, allocatable :: x(:), yr(:,:), yi(:,:), yppr(:,:), yppi(:,:)
  COMPLEX*16, ALLOCATABLE :: evals(:) 
  REAL*8, ALLOCATABLE     :: Einf(:)
  INTEGER, ALLOCATABLE    :: cind(:), cini(:), Sigini(:,:)
  INTEGER :: fh_Eimpx, iomega0
  LOGICAL :: Off_diagonal, write_Eimpx
  !---------------------------
  ! Example for index arrays:
  !---------------------------
  !  Sigind=
  !  [[1 0 0 0 0 0 0 0 0 0]
  !   [0 2 0 0 0 0 0 0 0 0]
  !   [0 0 0 0 0 0 0 0 0 0]
  !   [0 0 0 3 0 0 0 0 0 0]
  !   [0 0 0 0 4 0 0 0 0 0]
  !   [0 0 0 0 0 5 0 0 0 0]
  !   [0 0 0 0 0 0 6 0 0 0]
  !   [0 0 0 0 0 0 0 0 0 0]
  !   [0 0 0 0 0 0 0 0 0 0]
  !   [0 0 0 0 0 0 0 0 0 0]
  !  cixdm = 10
  !  cixdms = 6
  !  cind=[1,2,0,3,4,5,6,0,0,0]
  !  cini=[1,2,4,5,6,7]
  !  Sigini=[[1,0,0,0,0,0],[0,2,0,0,0,0],[0,0,4,0,0,0],[0,0,0,5,0,0],[0,0,0,0,6,0],[0,0,0,0,0,7]]


  !-----------------------
  DATA IMAG/(0.0D0,1.0D0)/
  !-----------------------
  write_Eimpx = .True.
  fh_Eimpx = 6
  !filename = 'Eimpx.dat'
  !open(fh_Eimpx,file=TRIM(filename),status='unknown')
  if (write_Eimpx) WRITE(fh_Eimpx,*) 'Full matrix of impurity levels follows'

  Deltac=0
  Olapc=0
  Glc=0
  Eimpc=0
  
  IF (matsubara) THEN
     PI=ACOS(-1.0D0)
     CALL CountFrequencyShape(n0_om, nom_max, beta, omega, nomega)
  END IF
  iomega0=1
  do iom=2,nomega
     if (abs(omega(iom)).LT.abs(omega(iomega0))) iomega0=iom
  enddo
  
  DO icix=1,ncix

     cixdm = cixdim(icix)
     allocate( cind(cixdm) )
     
     ! If Sigind(i,i)=0, we eliminate the i-th column and rown, because such orbital should be treated as non-correlated
     cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated
     cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block
     DO ip=1,cixdm
        it = Sigind(ip,ip,icix)
        if (it.gt.0) then
           cixdms = cixdms + 1
           cind(ip) = cixdms
        endif
     ENDDO

     !print *, 'icix=', icix, 'cixdms=', cixdms
     if (cixdms.EQ.0) then
        deallocate( cind )
        ! Even though all states are projected out, we still want to have GF printed.
        do iom=1,nomega !------- for all frequencies --!
           !---- packing Gloc to vector form

           DO ip=1,cixdm
              do iq=1,cixdm
                 it = abs(Sigind(ip,iq,icix))
                 if (it.ne.0) then
                    Glc( it, icix, iom ) =  Glc( it, icix, iom ) + gmloc(ip,iq,icix,iom)
                 endif
              enddo
           ENDDO
           DO it=1,csize(icix)
              Glc(it, icix, iom) = Glc(it, icix, iom)/noccur(it,icix)
           ENDDO
        enddo

        CYCLE
     endif

     allocate( cini(cixdms), Sigini(cixdms,cixdms) )
     do ip=1,cixdm
        if (cind(ip).gt.0) cini(cind(ip))=ip
     enddo
     DO ip=1,cixdms
        do iq=1,cixdms
           Sigini(ip,iq) = abs(Sigind(cini(ip),cini(iq),icix)) ! This abs is probably not needed, since we do not expect negative Sigind on the off-diagonal components
        enddo
     ENDDO
     !print *, 'cind=', cind
     !print *, 'cini=', cini
     !print *, 'Sigini=', Sigini

     allocate( sigcx(cixdms,cixdms), gcx(cixdms,cixdms), Eimpx(cixdms,cixdms), Eim0x(cixdms,cixdms), deltx(cixdms,cixdms), sig_oo(cixdms,cixdms), Einf(cixdms), olapx(cixdms,cixdms), tmp(cixdms,cixdms) )
     
     
     !---- Olap to vector form, and s_oo to matrix form
     DO ip=1,cixdm
        do iq=1,cixdm
           it = abs(Sigind(ip,iq,icix))
           if (it.gt.0) then
              Olapc(it, icix) =  Olapc(it, icix) + real(Olapm(ip,iq,icix))
           endif
        enddo
     ENDDO
     DO it=1,csize(icix)
        Olapc(it, icix) = Olapc(it, icix)/noccur(it,icix)
     ENDDO
     
     ! Set s_oo
     sig_oo=0
     DO ip=1,cixdms
        do iq=1,cixdms
           it = Sigini(ip,iq)
           if (it.gt.0) then
              sig_oo(ip,iq) = s_oo(it,icix)
           endif
        enddo
     ENDDO
     ! Set olapx, Eim0x
     DO ip=1,cixdms
        do iq=1,cixdms
           olapx(ip,iq) = Olapm(cini(ip),cini(iq),icix)
           Eim0x(ip,iq) = Eimpm(cini(ip),cini(iq),icix)
           if (projector.LT.0 .and. ip.NE.iq) then
              ! Set off-diagonal terms to zero. DMFT-SCC is taken to be purely diagonal
              olapx(ip,iq)=0.0
              Eim0x(ip,iq)=0.0
           end if
        enddo
     enddo

     ! DMFT SCC in non-orthogonal base:
     ! 1/(O*omega - E_imp - sigma - Delta) = \sum_k G_k 
     ! Olapm = sum_k U^+_k U_k
     ! Eimpm = sum_k U^+_k (eps_k-mu) U_k
     ! O^{-1}/omega + O^{-1}E_{imp}O^{-1}/omega^2 + ... = Olapm/omega + Eimpm/omega^2 + ...
     ! O = Olapm^{-1}
     ! E_{imp} = O * Eimpm * O
     ! 
     CALL zinv(olapx, cixdms)
     tmp = matmul(olapx,Eim0x)
     Eimpx = matmul( tmp, olapx) !-- normalize by overlap
     CALL GetEinf(Einf, Eimpx, Sigini, noccur(:,icix), cixdms, csize(icix), maxsize)
     Eimpx = Eimpx - sig_oo !--- subtracting S_oo from Eimp such that the real impurity level is Eimpc-Edc
     
     !---- Eimpc to vector form
     DO ip=1,cixdms
        do iq=1,cixdms
           it = Sigini(ip,iq)
           if (it.gt.0) then
              Eimpc(it, icix) =  Eimpc(it, icix) + real(Eimpx(ip,iq))
           endif
        enddo
     ENDDO
     DO it=1,csize(icix)
        Eimpc(it, icix) = Eimpc(it, icix)/noccur(it,icix)
     ENDDO
     
     
     if (write_Eimpx) then
        WRITE(fh_Eimpx,*) 'icix=', icix
        DO ip=1,cixdms
           DO iq=1,cixdms
              WRITE(fh_Eimpx,'(f14.8,1x,f14.8,3x)',advance='no') real(Eimpx(ip,iq))*Ry2eV, aimag(Eimpx(ip,iq))*Ry2eV
           ENDDO
           WRITE(fh_Eimpx,*)
        ENDDO
        WRITE(fh_Eimpx,*)
     end if

     do iom=1,nomega !------- for all frequencies --!
        
        IF (matsubara) THEN
           xomega = omega(iom)*IMAG
        ELSE
           xomega = omega(iom)
        ENDIF
        !---- unpacking sigma to matrix form
        !gcx = gmloc(:cixdm,:cixdm,icix,iom)
        sigcx=0
        gcx=0
        DO ip=1,cixdms
           do iq=1,cixdms
              !!! Sigma_loc
              it = Sigini(ip,iq)
              if (it.gt.0) then
                 sigcx(ip,iq) = sigma( it, icix, iom )
              endif
              !!! G_loc
              gcx(ip,iq) = gmloc(cini(ip),cini(iq),icix,iom)
           enddo
        ENDDO

        if (projector.LT.0) then
           ! Set off-diagonal terms to zero. DMFT-SCC is taken to be purely diagonal
           DO ip=1,cixdms
              do iq=1,cixdms
                 if (ip.NE.iq) then
                    gcx(ip,iq)=0.0
                 end if
              enddo
           enddo
        endif

        CALL zinv(gcx,cixdms)
        
        ! DMFT SCC in non-orthogonal base:
        ! 1/(O*omega - E_imp - sigma - Delta) = \sum_k G_k==G_loc
        ! Delta = O*omega - E_{imp} - sigma - G_loc^{-1}
        !--- DMFT self-conistency condition
        !deltx = Olapm(:cixdm,:cixdm,icix)*xomega - Eimpx - gcx - sigcx 
        deltx = olapx*xomega - Eimpx - gcx - sigcx 
        do ip=1,cixdms
           deltx(ip,ip) = deltx(ip,ip) + (0.d0, 1.d0)*gamma
        enddo


        if (iom==iomega0 .and. write_Eimpx) then
           tmp(:,:) = Eimpx(:,:) + 0.5*(deltx(:,:)+transpose(conjg(deltx(:,:))))
           WRITE(fh_Eimpx,*) 'icix=', icix, 'at omega=0'
           DO ip=1,cixdms
              DO iq=1,cixdms
                 WRITE(fh_Eimpx,'(f14.8,1x,f14.8,3x)',advance='no') real(tmp(ip,iq))*Ry2eV, aimag(tmp(ip,iq))*Ry2eV
              ENDDO
              WRITE(fh_Eimpx,*)
           ENDDO
           WRITE(fh_Eimpx,*)
        endif
        
        !---- packing delta to vector form
        DO ip=1,cixdms
           do iq=1,cixdms
              it = Sigini(ip,iq)
              if (it.gt.0) then
                 Deltac( it, icix, iom ) =  Deltac( it, icix, iom ) + deltx(ip,iq)
              endif
           enddo
        ENDDO
        !---- packing Gloc to vector form
        DO ip=1,cixdm
           do iq=1,cixdm
              it = abs(Sigind(ip,iq,icix))
              if (it.ne.0) then
                 Glc( it, icix, iom ) =  Glc( it, icix, iom ) + gmloc(ip,iq,icix,iom)
              endif
           enddo
        ENDDO
        DO it=1,csize(icix)
           Deltac(it, icix, iom) = Deltac(it, icix, iom)/noccur(it,icix)
           Glc(it, icix, iom) = Glc(it, icix, iom)/noccur(it,icix)
        ENDDO
     enddo

     logGloc(:)=0.0
     IF (matsubara .and. ComputeLogGloc) then  !!! This is new stuff on Tr(log(Gloc)) for free energy
        !!! At the moment this works only on imaginary axis. It should be extended to the real axis too.
        
        Off_diagonal = Off_diagonalc(Sigind, cixdm, maxdim)
        n0 = nomega-n0_om+1
        allocate( x(n0), yr(n0,csize(icix)), yi(n0,csize(icix)), yppr(n0,csize(icix)), yppi(n0,csize(icix)) )
        allocate( gxc(csize(icix)) )
        x(:) = omega(n0_om:)
        DO it=1,csize(icix)
           yr(:,it) = dble(Glc(it,icix,n0_om:))
           call spline_cubic_set( n0, x, yr(:,it), 0, 0.0D+00, 0, 0.0D+00, yppr(:,it) )
           yi(:,it) = dimag(Glc(it,icix,n0_om:))
           call spline_cubic_set( n0, x, yi(:,it), 0, 0.0D+00, 0, 0.0D+00, yppi(:,it) )
        ENDDO

        ALLOCATE( evals(cixdms) )

        logGloc(icix)=0.0
        DO iom=1,n0_om-1
           CALL GlocEigvals(evals, Glc(:,icix,iom), Sigini, Off_diagonal, cixdms, csize(icix))
           zsum=0.0
           do it=1,cixdms
              zsum = zsum + log(-evals(it)) + log(Einf(it)-omega(iom)*IMAG)
           enddo
           logGloc(icix) = logGloc(icix) + dble(zsum)*2/beta
        ENDDO
        left=0
        DO iom=n0_om,nom_max
           omw = (2*iom-1)*PI/beta
           DO it=1,csize(icix)
              call spline_cubic_val2(n0, x, yr(:,it), yppr(:,it), left, omw, ywr, ypwr, yppwr )
              call spline_cubic_val2(n0, x, yi(:,it), yppi(:,it), left, omw, ywi, ypwi, yppwi )
              gxc(it) = ywr+ywi*IMAG
           ENDDO
           CALL GlocEigvals(evals, gxc, Sigini, Off_diagonal, cixdms, csize(icix))
           zsum=0.0
           do it=1,cixdms
              zsum = zsum + log(-evals(it)) + log(Einf(it)-omw*IMAG)
           enddo
           logGloc(icix) = logGloc(icix) + dble(zsum)*2/beta
        ENDDO
        do it=1,cixdms
           logGloc(icix) = logGloc(icix) + FreeE0(Einf(it),1./beta)
        enddo
        WRITE(*,*) 'icix=', icix, 'logGloc=', logGloc(icix)
        
        DEALLOCATE( evals )
        deallocate( gxc )
        deallocate( x, yr, yi, yppr, yppi )
     END IF
     
     deallocate( sigcx, gcx, Eimpx, Eim0x, deltx, sig_oo, Einf, olapx, tmp )
     deallocate( cini, Sigini, cind )
  ENDDO
  !close(fh_Eimpx)
END SUBROUTINE GetDelta2

REAL*8 FUNCTION FreeE0(Energy, Temperature)
  IMPLICIT NONE
  REAL*8, intent(in) :: Energy, Temperature
  if (Energy/Temperature.LT.-200) then
     FreeE0 = Energy
     return
  endif
  if (Energy/Temperature.GT.200) then
     FreeE0 = 0.0
     return
  endif
  FreeE0 = -Temperature*log(1.+exp(-Energy/Temperature))
END FUNCTION FreeE0


SUBROUTINE GetEinf(Einf, Eimpx, Sigini, noccur, cixdms, csize, maxsize)
  IMPLICIT NONE
  INTEGER, intent(in)    :: cixdms
  REAL*8, intent(out)    :: Einf(cixdms)
  COMPLEX*16, intent(in) :: Eimpx(cixdms,cixdms)
  INTEGER, intent(in)    :: Sigini(cixdms,cixdms), noccur(maxsize), csize, maxsize
  ! locals
  REAL*8, allocatable :: Einf0(:)
  INTEGER :: ip, it
  
  ALLOCATE( Einf0(csize) )
  DO ip=1,cixdms
     it = Sigini(ip,ip)
     if (it.gt.0) then
        Einf0(it) = Einf0(it) + dble(Eimpx(ip,ip))
     endif
  ENDDO
  DO it=1,csize
     Einf0(it) = Einf0(it)/noccur(it)
  ENDDO
  
  Einf=0
  DO ip=1,cixdms
     it = Sigini(ip,ip)
     if (it.gt.0) then
        Einf(ip) = Einf0(it)
     endif
  ENDDO
  DEALLOCATE( Einf0 )
END SUBROUTINE GetEinf

LOGICAL FUNCTION Off_diagonalc(Sigind, cixdm, maxdim)
  IMPLICIT NONE
  !
  INTEGER, intent(in) :: Sigind(maxdim,maxdim), cixdm, maxdim
  LOGICAL :: Off_diagonal
  INTEGER :: ip, iq
  Off_diagonal=.False.
  DO ip=1,cixdm
     do iq=1,cixdm
        if (ip.NE.iq .AND. Sigind(ip,iq).GT.0) then
           Off_diagonal=.true.
           exit
        endif
     enddo
  ENDDO
  Off_diagonalc = Off_diagonal
  WRITE(*,*) 'Off_diagonal=', Off_diagonal
END FUNCTION Off_diagonalc

SUBROUTINE GlocEigvals(evals, gxc, Sigini, Off_diagonal, cixdms, csize)
  IMPLICIT NONE
  INTEGER, intent(in)     :: cixdms, csize
  COMPLEX*16, intent(out) :: evals(cixdms)
  COMPLEX*16, intent(in)  :: gxc(csize)
  INTEGER, intent(in) :: Sigini(cixdms,cixdms)
  LOGICAL, intent(in) :: Off_diagonal
  !! locals
  CHARACTER*1, PARAMETER :: jobvl = "N"  ! Jop parameter for zgeev
  CHARACTER*1, PARAMETER :: jobvr = "N"  ! Job parameter for zgeev
  COMPLEX*16, ALLOCATABLE :: cworkvec(:), evx(:,:), gcx(:,:)
  REAL*8, ALLOCATABLE     :: rworkvec(:)
  INTEGER :: ierr, ip, iq, it
  
  if (Off_diagonal) then
     ALLOCATE( cworkvec(4*cixdms), evx(cixdms,cixdms), gcx(cixdms,cixdms), rworkvec(2*cixdms) )
     gcx(:,:)=0.0
     DO ip=1,cixdms
        do iq=1,cixdms
           it = Sigini(ip,iq)
           if (it.gt.0) gcx(ip,iq) = gxc(it)
        enddo
     ENDDO
     CALL zgeev(jobvl,jobvr,cixdms,gcx,cixdms,evals,evx,cixdms,evx,cixdms,cworkvec,4*cixdms,rworkvec,ierr)
     IF(ierr.NE.0) THEN
        WRITE(6,*) 'Error code of zgeev ', ierr
        WRITE(0,*) 'Error in dia_gho! Stopp the code!'
     ENDIF
     DEALLOCATE( cworkvec, evx, gcx, rworkvec )
  else
     evals(:)=0
     DO ip=1,cixdms
        it = Sigini(ip,ip)
        if (it.gt.0) evals(ip) = gxc(it)
     ENDDO
  endif
END SUBROUTINE GlocEigvals
  
SUBROUTINE RenormalizeTransK(DMFTU, cix_orb, cixdim, nindo, iSx, Sigind, projector, nbands, maxdim2, norbitals, maxdim, ncix)
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: DMFTU(nbands,maxdim2,norbitals)
  INTEGER, intent(in)       :: cix_orb(norbitals), nindo(norbitals), iSx(maxdim2,norbitals), cixdim(ncix), Sigind(maxdim,maxdim,ncix)
  INTEGER, intent(in)       :: nbands, maxdim2, norbitals, maxdim, ncix, projector
  ! locals
  INTEGER :: iorb1, icix, nind1, ind1, cixdm
  REAL*8  :: olocef
  !
  INTEGER :: cixdms, cixdms_m, info, ip, lwork, lrwork, i
  INTEGER,    allocatable :: cind(:), cini(:), iwork(:)
  COMPLEX*16, allocatable :: Ucix(:,:), work(:), Uw(:,:), Vw(:,:)
  REAL*8,     allocatable :: rwork(:), ws(:)
  !
  DO icix=1,ncix
     cixdm = cixdim(icix)
     
     ! If we project out some orbitals, for example eg-orbitals, we might have
     ! Sigind= [[0,0,0,0,0], [0,0,0,0,0], [0,0,1,0,0], [0,0,0,2,0], [0,0,0,0,3]]
     ! then
     !     cind[1:5] = [0,0,1,2,3] and
     !     cini[1:3] = [3,4,5]
     allocate( cind(cixdm) )
     cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated                                                                                    
     cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block                                                                                         
     DO ip=1,cixdm
        if (Sigind(ip,ip,icix) .ne. 0) then
           cixdms = cixdms + 1
           cind(ip) = cixdms
        endif
     ENDDO
     
     allocate( cini(cixdms))
     do ip=1,cixdm
        if (cind(ip).gt.0) cini(cind(ip))=ip
     enddo
     ! Now cini[1:3] = [3,4,5] contains the small index of non-zero components
     
     ! If we have cluster-DMFT calculations, we need several orbitals combined into cix block
     allocate( Ucix(nbands,cixdms) )
     Ucix(:,:) = 0.d0
     DO iorb1=1,norbitals
        if ( cix_orb(iorb1).NE.icix ) CYCLE
        nind1 = nindo(iorb1)
        do ind1=1,nind1
           ip = iSx(ind1,iorb1)
           if (cind(ip).gt.0) Ucix(:,cind(ip)) = DMFTU(:,ind1,iorb1)
        enddo
     ENDDO

     cixdms_m = min(cixdms,nbands) ! should normally be cixdms, as the number of bands should be larger
     allocate( ws(cixdms_m), Uw(nbands,cixdms_m), Vw(cixdms_m,cixdms) )
     
     lwork = 2*cixdms*(cixdms+nbands)
     lrwork = 7*cixdms*(cixdms + 1)
     allocate( work(lwork), rwork(lrwork), iwork(8*cixdms) )

     !do i=1,nbands
     !   WRITE(6,'(A)', advance='no') "["
     !   do ip=1,cixdms
     !      WRITE(6, '(f14.10,"+",f8.5,3x,"*1j, ")',advance='no') real(Ucix(i,ip)), aimag(Ucix(i,ip))
     !   enddo
     !   WRITE(6,*) "],"
     !enddo

     
     call ZGESDD('S', nbands, cixdms, Ucix, nbands, ws, Uw, nbands, Vw, cixdms_m, work, lwork, rwork, iwork, info )
     if (info .ne. 0) then
        print *, 'SVD decomposition of the projector failed. Info-zgesdd=', info
        if (info.lt.0) print *, 'The ', abs(info),' th argument had an illegal value.'
        if (info.gt.0) print *, 'The updating process of DBDSDC did not converge.'
     endif
     call zgemm('N', 'N', nbands, cixdms, cixdms_m, (1.d0,0.d0), Uw, nbands, Vw, cixdms_m, (0.d0,0.d0), Ucix, nbands)


     !WRITE(6,*)
     !do i=1,nbands
     !   do ip=1,cixdms
     !      WRITE(6, '(f14.10,1x,f8.5,3x)',advance='no') real(Ucix(i,ip)), aimag(Ucix(i,ip))
     !   enddo
     !   WRITE(6,*)
     !enddo
     print *, 'Singular values=', ws

     
     deallocate( work, rwork, iwork )
     deallocate( ws, Uw, Vw )
     
     DO iorb1=1,norbitals
        if ( cix_orb(iorb1).NE.icix ) CYCLE
        nind1 = nindo(iorb1)
        do ind1=1,nind1
           ip = iSx(ind1,iorb1)
           if (cind(ip).gt.0) DMFTU(:,ind1,iorb1) = Ucix(:,cind(ip))
        enddo
     ENDDO
     
     deallocate( cind, cini )
     deallocate( Ucix )
  ENDDO
END SUBROUTINE RenormalizeTransK

SUBROUTINE RenormalizeTrans(DMFTU, Olapm0, SOlapm, cix_orb, cixdim, nindo, iSx, projector, nbands, maxdim2, norbitals, maxdim, ncix, SIMPLE)
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: DMFTU(nbands,maxdim2,norbitals)
  COMPLEX*16, intent(in)    :: Olapm0(maxdim,maxdim,ncix), SOlapm(maxdim,maxdim,ncix)
  INTEGER, intent(in)       :: cix_orb(norbitals), nindo(norbitals), iSx(maxdim2,norbitals), cixdim(ncix)
  INTEGER, intent(in)       :: nbands, maxdim2, norbitals, maxdim, ncix, projector
  LOGICAL, intent(in)       :: SIMPLE
  ! locals
  INTEGER :: iorb1, icix, nind1, ind1, cixdm
  REAL*8  :: olocef
  COMPLEX*16, allocatable :: tmp3(:,:), Ucix(:,:), Ucix2(:,:)
  
  if (SIMPLE) then
     DO iorb1=1,norbitals
        icix = cix_orb(iorb1)
        if ( icix.EQ.0 ) CYCLE
        nind1 = nindo(iorb1)

        do ind1=1,nind1
           olocef = 1/sqrt(real(Olapm0( iSx(ind1,iorb1), iSx(ind1,iorb1), icix )))
           DMFTU(:,ind1,iorb1) = DMFTU(:,ind1,iorb1) * olocef
        enddo
     ENDDO
  else if (.True.) then
     DO icix=1,ncix
        cixdm = cixdim(icix)
        allocate( Ucix(nbands,cixdm), Ucix2(nbands,cixdm) )
        Ucix(:,:) = 0.d0
        DO iorb1=1,norbitals
           if ( cix_orb(iorb1).NE.icix ) CYCLE
           nind1 = nindo(iorb1)
           do ind1=1,nind1
              Ucix(:,iSx(ind1,iorb1)) = DMFTU(:,ind1,iorb1)
           enddo
        ENDDO
        call zgemm('N','N', nbands, cixdm, cixdm, (1.d0,0.d0), Ucix, nbands, SOlapm(:,:,icix), maxdim, (0.d0,0.d0), Ucix2, nbands)
        DO iorb1=1,norbitals
           if ( cix_orb(iorb1).NE.icix ) CYCLE
           nind1 = nindo(iorb1)
           do ind1=1,nind1
              DMFTU(:,ind1,iorb1) = Ucix2(:,iSx(ind1,iorb1))
           enddo
        ENDDO
        deallocate( Ucix, Ucix2 )
     ENDDO
  else
     allocate( tmp3(nbands,maxdim2) )
     DO iorb1=1,norbitals
        icix = cix_orb(iorb1)
        if ( icix.EQ.0 ) CYCLE
        nind1 = nindo(iorb1)

        call zgemm('N','N', nbands, nind1, nind1, (1.d0,0.d0), DMFTU(:,:,iorb1),nbands, SOlapm(:,:,icix),maxdim, (0.d0,0.d0), tmp3, nbands)
        DMFTU(:,:nind1,iorb1) = tmp3(:,:nind1)
     ENDDO
     deallocate( tmp3 )
  endif
     
  
END SUBROUTINE RenormalizeTrans

SUBROUTINE read_overlap_from_file(info, SOlapm, cixdim, maxdim, ncix)
  USE com_mpi,ONLY: myrank, master
  IMPLICIT NONE
  INTEGER,    intent(out) :: info
  COMPLEX*16, intent(out) :: SOlapm(maxdim,maxdim,ncix)
  INTEGER,    intent(in)  :: cixdim(ncix), maxdim, ncix
  ! locals
  LOGICAL :: there
  INTEGER :: cixdm, ind1, ind2, icix

  INQUIRE( FILE='SOlapm.dat', EXIST=there)
  if (.not.there) then
     info=1
     return
  endif
  
  open(996, FILE='SOlapm.dat', status='old')
  DO icix=1,ncix
     cixdm = cixdim(icix)
     READ(996,*,ERR=991) !icix
     do ind1=1,cixdm
        do ind2=1,cixdm
           READ(996,'(2F20.16)',advance='no',ERR=991)  SOlapm(ind1,ind2,icix)
        enddo
        READ(996,*)
     enddo
  ENDDO
  close(996)
  info=0
  
  if (myrank.eq.master) then
     WRITE(6,*) 'SOlapm succesfully read from file'
     WRITE(6,*) 'Z will be used:'
     DO icix=1,ncix
        cixdm = cixdim(icix)
        do ind1=1,cixdm
           WRITE(6,'(F16.10)',advance='no')  1./dble(SOlapm(ind1,ind1,icix))
        enddo
        WRITE(6,*)
     ENDDO
     WRITE(6,*)
  endif

  return
991 CONTINUE
  info=2
  close(996)
END SUBROUTINE read_overlap_from_file
