module sym2
integer,parameter :: iord=2
end module
module dmfts
integer,parameter :: maxdim=3,ncix=2,natom=2,lmaxp=0
end module
module com_mpi
integer :: myrank=1,master=0
end module
module test_rotations
complex*16 :: matrices(3,3,2,2)
end module
subroutine set_test_rotations(r)
use test_rotations
complex*16,intent(in) :: r(3,3,2,2)
matrices=r
end subroutine
subroutine CmpSymmetryTransformation(ig,RT,cfX,cix_orb,cixdim,iSx, &
 iorbital,nindo,norbitals,maxdim2,verbose)
use test_rotations
integer :: ig,norbitals,maxdim2,cix_orb(norbitals),cixdim(2), &
 iSx(maxdim2,norbitals),iorbital(2,1),nindo(norbitals)
complex*16 :: RT(3,3,2),cfX(maxdim2,maxdim2,norbitals,norbitals)
logical :: verbose
RT=matrices(:,:,:,ig)
end subroutine
subroutine FindUniqueCix(ng,sizes,groups)
integer,intent(out) :: ng,sizes(2),groups(2,2)
ng=1; sizes=0; groups=0
sizes(1)=2; groups(1,1)=1; groups(1,2)=2
end subroutine
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


SUBROUTINE SymmetrizeNormalizationAdjoint(Qbar,cfX,cix_orb,cixdim, &
     iSx,iorbital,nindo,norbitals,maxdim2)
  USE sym2,   ONLY: iord
  USE dmfts,  ONLY: maxdim,ncix,natom,lmaxp
  USE com_mpi,ONLY: myrank,master
  IMPLICIT NONE
  COMPLEX*16,intent(inout) :: Qbar(maxdim,maxdim,ncix)
  COMPLEX*16,intent(in) :: cfX(maxdim2,maxdim2,norbitals,norbitals)
  INTEGER,intent(in) :: cix_orb(norbitals),cixdim(ncix)
  INTEGER,intent(in) :: iSx(maxdim2,norbitals)
  INTEGER,intent(in) :: iorbital(natom,lmaxp+1),nindo(norbitals)
  INTEGER,intent(in) :: norbitals,maxdim2
  COMPLEX*16,allocatable :: Qraw(:,:,:),Qsym(:,:,:),RT(:,:,:)
  CHARACTER(len=32) :: mode
  INTEGER :: status,i,letter,ig,icix,cdim
  REAL*8 :: rawNorm,differenceNorm,relativeDifference,maxChange

  mode='APPLY'
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_NORM_SYMMETRY',mode, &
       STATUS=status)
  if (status.ne.0 .or. len_trim(mode).eq.0) mode='APPLY'
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
     STOP 'invalid normalization-adjoint symmetry mode'
  endif

  allocate(Qraw(maxdim,maxdim,ncix),Qsym(maxdim,maxdim,ncix), &
       RT(maxdim,maxdim,ncix))
  ! BuildNormalizationAdjoint returns conjg(Q), with Q Hermitian and
  ! dGamma_norm = 2 Re Tr(Q dS). Undo this storage convention first.
  ! The FORWARD overlap map is Herm o EquivalentAverage o OrbitalAverage:
  !     S = Herm E [sum_g R_g^dagger Sraw R_g / |G|].
  ! Reverse these operations for its exact ADJOINT. In particular the
  ! equivalent-cix average precedes the orbital rotations here, and the
  ! products are R_g Q R_g^dagger, not R_g^dagger Q R_g. No commutation of
  ! X, Q or equivalent-block rotations is assumed. This remains correct for
  ! a saved normalization X which is only approximately symmetry invariant.
  Qraw=conjg(Qbar)
  do icix=1,ncix
     cdim=cixdim(icix)
     Qraw(:cdim,:cdim,icix)=0.5d0*(Qraw(:cdim,:cdim,icix)+ &
          conjg(transpose(Qraw(:cdim,:cdim,icix))))
  enddo
  CALL AverageEquivalentCixMatrices(Qraw,cixdim,maxdim,ncix)
  Qsym=0.d0
  do ig=1,iord
     CALL CmpSymmetryTransformation(ig,RT,cfX,cix_orb,cixdim,iSx, &
          iorbital,nindo,norbitals,maxdim2,.False.)
     do icix=1,ncix
        cdim=cixdim(icix)
        Qsym(:cdim,:cdim,icix)=Qsym(:cdim,:cdim,icix)+ &
             matmul(RT(:cdim,:cdim,icix), &
             matmul(Qraw(:cdim,:cdim,icix), &
             conjg(transpose(RT(:cdim,:cdim,icix)))))/dble(iord)
     enddo
  enddo

  rawNorm=sqrt(sum(abs(Qbar)**2))
  differenceNorm=sqrt(sum(abs(Qsym-conjg(Qbar))**2))
  relativeDifference=differenceNorm/max(rawNorm,tiny(1.d0))
  maxChange=maxval(abs(Qsym-conjg(Qbar)))
  if (myrank.eq.master) WRITE(6,'(A,1X,A,3(1X,ES20.11))') &
       ':DMFT_NORMALIZATION_ADJOINT_SYMMETRY:',trim(mode), &
       rawNorm,relativeDifference,maxChange
  if (trim(mode).eq.'APPLY') Qbar=conjg(Qsym)

  deallocate(RT,Qsym,Qraw)
END SUBROUTINE SymmetrizeNormalizationAdjoint


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




SUBROUTINE BuildNormalizationAdjoint(Qbar,momentBar,X,Sigind,cixdim, &
     maxdim,ncix)
  IMPLICIT NONE
  INTEGER,intent(in) :: maxdim,ncix,cixdim(ncix), &
       Sigind(maxdim,maxdim,ncix)
  COMPLEX*16,intent(in) :: momentBar(maxdim,maxdim,ncix), &
       X(maxdim,maxdim,ncix)
  ! Qbar=conjg(Q) is returned because all projector-gradient contractions in
  ! l2main use the conjugated-gradient convention.
  COMPLEX*16,intent(out) :: Qbar(maxdim,maxdim,ncix)
  INTEGER :: icix,ip,iq,na,info,lwork,lrwork,liwork
  INTEGER,allocatable :: active(:),iwork(:)
  REAL*8,allocatable :: xeig(:),rwork(:)
  REAL*8 :: lij,directTest,adjointTest,adjointScale
  LOGICAL :: projectorDiagnostic
  CHARACTER(LEN=64) :: diagnosticEnvironment
  INTEGER :: diagnosticStatus
  COMPLEX*16,allocatable :: Xvec(:,:),Nmat(:,:),Np(:,:),Q0p(:,:), &
       Qp(:,:),Qmat(:,:),work(:),dStest(:,:),dStestPrime(:,:), &
       BtestPrime(:,:),Btest(:,:),Mtest(:,:)

  Qbar=0.d0
  diagnosticEnvironment=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_PROJECTOR_DIAG', &
       diagnosticEnvironment,STATUS=diagnosticStatus)
  projectorDiagnostic=diagnosticStatus.eq.0 .and. &
       len_trim(diagnosticEnvironment).gt.0 .and. &
       trim(diagnosticEnvironment).ne.'0'
  do icix=1,ncix
     allocate(active(cixdim(icix)))
     na=0
     do ip=1,cixdim(icix)
        if (any(Sigind(ip,1:cixdim(icix),icix).gt.0) .or. &
             any(Sigind(1:cixdim(icix),ip,icix).gt.0)) then
           na=na+1
           active(na)=ip
        endif
     enddo
     if (na.eq.0) then
        deallocate(active)
        cycle
     endif

     allocate(Xvec(na,na),Nmat(na,na),Np(na,na),Q0p(na,na), &
          Qp(na,na),Qmat(na,na),xeig(na))
     do ip=1,na
        do iq=1,na
           Xvec(ip,iq)=X(active(ip),active(iq),icix)
           Nmat(ip,iq)=conjg(momentBar(active(ip),active(iq),icix))
        enddo
     enddo

     lwork=4*na+na*na
     lrwork=10*na+2*na*na+1
     liwork=3+10*na
     allocate(work(lwork),rwork(lrwork),iwork(liwork))
     CALL ZHEEVD('V','U',na,Xvec,na,xeig,work,lwork,rwork,lrwork, &
          iwork,liwork,info)
     if (info.ne.0 .or. minval(xeig).le.1.d-14) then
        WRITE(6,*) 'ERROR: normalization adjoint diagonalization failed', &
             icix,info,minval(xeig)
        STOP 'BuildNormalizationAdjoint failed'
     endif

     Np=matmul(conjg(transpose(Xvec)),matmul(Nmat,Xvec))
     do ip=1,na
        do iq=1,na
           lij=-xeig(iq)/(1.d0/xeig(ip)+1.d0/xeig(iq))
           Q0p(ip,iq)=lij*Np(ip,iq)
        enddo
     enddo
     ! delta S is Hermitian.  Only the Hermitian part of the unconstrained
     ! adjoint can contribute to its real Frobenius contraction.
     Qp=0.5d0*(Q0p+conjg(transpose(Q0p)))
     Qmat=matmul(Xvec,matmul(Qp,conjg(transpose(Xvec))))
     do ip=1,na
        do iq=1,na
           Qbar(active(ip),active(iq),icix)=conjg(Qmat(ip,iq))
        enddo
     enddo

     if (projectorDiagnostic) then
        ! Verify the adjoint without reference to any numerical stress.  For
        ! a deterministic Hermitian trial dS, form B=X^{-1}dX with the same
        ! Sylvester map used by ComputeNormalizationDerivative.  The direct
        ! projector contraction and the adjoint contraction must agree:
        !
        !   2 Re sum(M*B) = 2 Re sum(Qbar*dS).
        !
        ! This catches a misplaced transpose or conjugation in the global
        ! normalization response even when X happens to be nearly diagonal.
        allocate(dStest(na,na),dStestPrime(na,na), &
             BtestPrime(na,na),Btest(na,na),Mtest(na,na))
        do ip=1,na
           do iq=1,na
              dStest(ip,iq)=cmplx(0.1d0*dble(ip+iq), &
                   0.03d0*dble(ip-iq),kind=8)
              Mtest(ip,iq)=momentBar(active(ip),active(iq),icix)
           enddo
        enddo
        dStestPrime=matmul(conjg(transpose(Xvec)), &
             matmul(dStest,Xvec))
        do ip=1,na
           do iq=1,na
              lij=-xeig(iq)/(1.d0/xeig(ip)+1.d0/xeig(iq))
              BtestPrime(ip,iq)=lij*dStestPrime(ip,iq)
           enddo
        enddo
        Btest=matmul(Xvec,matmul(BtestPrime, &
             conjg(transpose(Xvec))))
        directTest=2.d0*dble(sum(Mtest*Btest))
        adjointTest=2.d0*dble(sum(conjg(Qmat)*dStest))
        adjointScale=max(abs(directTest),abs(adjointTest),tiny(1.d0))
        WRITE(6,'(A,I4,3ES22.12)') ':DMFT_NORMALIZATION_ADJOINT: ', &
             icix,directTest,adjointTest, &
             abs(directTest-adjointTest)/adjointScale
        deallocate(dStest,dStestPrime,BtestPrime,Btest,Mtest)
     endif

     deallocate(active,Xvec,Nmat,Np,Q0p,Qp,Qmat,xeig,work,rwork,iwork)
  enddo
END SUBROUTINE BuildNormalizationAdjoint


  SUBROUTINE ProjectorGlobalAxes(axes,xrot,br,rotfrac)
    IMPLICIT NONE
    REAL*8,intent(out) :: axes(3,3)
    REAL*8,intent(in) :: xrot(3,3),br(3,3),rotfrac(3,3)
    REAL*8 :: brinv(3,3)
    call inv_3x3(br,brinv)
    axes=matmul(xrot,matmul(br,matmul(rotfrac,brinv)))
  END SUBROUTINE ProjectorGlobalAxes

SUBROUTINE inv_3x3(A,Ainv)
  IMPLICIT NONE
  REAL*8, intent(in)  :: A(3,3)
  REAL*8, intent(out) :: Ainv(3,3)
  ! locals
  REAL*8 :: det
  det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) &
      - A(3,1)*A(2,2)*A(1,3) - A(1,1)*A(3,2)*A(2,3) - A(2,1)*A(1,2)*A(3,3)
  Ainv(1,1) =(   A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / det
  Ainv(1,2) =( - A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / det
  Ainv(1,3) =(   A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / det
  Ainv(2,1) =( - A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / det
  Ainv(2,2) =(   A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / det
  Ainv(2,3) =( - A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / det
  Ainv(3,1) =(   A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / det
  Ainv(3,2) =( - A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / det
  Ainv(3,3) =(   A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / det
  RETURN
END SUBROUTINE inv_3x3

