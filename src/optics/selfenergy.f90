! @Copyright 2007 Kristjan Haule

FUNCTION CountSelfenergy(fh_sig, ncorr)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh_sig, ncorr
  INTEGER :: CountSelfEnergy
  ! locals
  REAL*8 :: sigma_vec(2*ncorr), ome
  INTEGER:: nom, ios, i
  nom=0
  ios=0
  DO
     ! read single complex number for each correlation index
     READ(fh_sig,*,IOSTAT=ios) ome,(sigma_vec(i),i=1,2*ncorr) 
     IF (ios.NE.0) EXIT ! If the read statement was not successful, exit
     nom=nom+1
  ENDDO
  IF (ios.GT.0) print *, 'ERROR Readings self energy file', fh_sig, ios, nom
  CountSelfenergy = nom
  REWIND(fh_sig)
  RETURN
END FUNCTION CountSelfenergy

SUBROUTINE ReadSelfenergy(fh_sig, sigma, omega, gammac, ncorr, nom, maxsize)
  IMPLICIT NONE
  INTEGER, intent(in)     :: fh_sig, ncorr, nom, maxsize
  REAL*8, intent(out)     :: omega(nom)
  REAL*8, intent(in)      :: gammac
  COMPLEX*16, intent(out) :: sigma(nom,maxsize)
  ! locals
  COMPLEX*16 :: s_oo(maxsize)
  REAL*8 :: sigma_vec(2*ncorr), ome
  INTEGER:: ios, iom, i
  sigma=0
  ios=0
  DO iom=1,nom
     ! read single complex number for each correlation index
     READ(fh_sig,*,IOSTAT=ios) ome,(sigma_vec(i),i=1,2*ncorr) 
     IF(ios.NE.0) then
        print *, 'ios=', ios, 'in reading self-energy' !EXIT ! If the read statement was not successful, exit
        exit
     endif
     omega(iom) = ome
     DO i=1,ncorr
        IF (ABS(sigma_vec(2*i)).LT.gammac) sigma_vec(2*i) = -gammac ! Set minimum broadening for correlated orbitals
        sigma(iom,i) = dcmplx(sigma_vec(2*i-1),sigma_vec(2*i))
     ENDDO
  ENDDO
  

  READ(fh_sig,*,IOSTAT=ios) ome,(sigma_vec(i),i=1,2*ncorr) 

  IF(ios.NE.0) print *, 'ERROR Readings s_oo', fh_sig, ios
  DO i=1,ncorr
     s_oo(i) = dcmplx(sigma_vec(2*i-1),sigma_vec(2*i))
  ENDDO

  IF (ios.GT.0) print *, 'ERROR Readings self energy file', fh_sig, ios
END SUBROUTINE ReadSelfenergy


SUBROUTINE CompressSigmaTransformation2(STrans, DMFTU, Sigind, iSx, cix, csize, iorbital, ll, nl, natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: STrans(maxsize,ncix,nbands,nbands)
  COMPLEX*16, intent(in)  :: DMFTU(nbands,maxdim2,norbitals)
  INTEGER, intent(in)     :: iSx(maxdim2, norbitals)
  INTEGER, intent(in)     :: Sigind(maxdim,maxdim,ncix), cix(natom,4), csize(ncix)
  INTEGER, intent(in)     :: iorbital(natom,lmax2+1), ll(natom,4), nl(natom)
  INTEGER, intent(in)     :: natom, iso, ncix, maxdim, maxdim2, lmax2, norbitals, nbands, maxsize
  !----- locals
  INTEGER    :: i, j, icase, jcase, lcase, l1case, l2case, iorb, iorb1, iorb2, it, ind1, nind1, ind2, nind2, icix, l1, l2, ip1, ip2
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
                    it = Sigind( iSx(ind1,iorb1),iSx(ind2,iorb2), icix )
                    
                    if (it.gt.0) then
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
  COMPLEX*16, intent(in)   :: STrans(maxsize,ncix,nbands,nbands)
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

