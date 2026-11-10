module g2_angular_test_module
use iso_fortran_env,only: dp=>real64
implicit none
integer,parameter :: axis_a(6)=[1,2,3,2,1,1],axis_b(6)=[1,2,3,3,3,2]
complex(dp),parameter :: one=(1._dp,0._dp),zero=(0._dp,0._dp)
contains
  subroutine gauss_legendre_stress(n,x,w)
    integer,intent(in) :: n
    real(dp),intent(out) :: x(n),w(n)
    real(dp) :: z,old,p0,p1,p2,deriv,pi
    integer :: i,j,it
    pi=acos(-1._dp)
    do i=1,(n+1)/2
       z=cos(pi*(real(i,dp)-0.25_dp)/(real(n,dp)+0.5_dp))
       do it=1,64
          p0=1._dp; p1=z
          do j=2,n
             p2=((2*j-1)*z*p1-(j-1)*p0)/j
             p0=p1; p1=p2
          enddo
          deriv=n*(z*p1-p0)/(z*z-1._dp)
          old=z; z=old-p1/deriv
          if(abs(z-old)<4*epsilon(z)) exit
       enddo
       x(i)=-z; x(n+1-i)=z
       w(i)=2._dp/((1._dp-z*z)*deriv*deriv); w(n+1-i)=w(i)
    enddo
  end subroutine
  subroutine g2_angular_tables(lmax,lpot,mpot,phase_plus,phase_minus,ag,ar,aa)
    ! <Y_i|Vharm|Y_j>, <Y_i|n_a n_b Vharm|Y_j>, and
    ! <A_ab Y_i|Vharm|Y_j>.  The derivative is on the BRA in the last
    ! expression.  Keeping it as a scalar-product integral eliminates the
    ! otherwise easy-to-miss complex conjugation in three-Gaunt recoupling.
    integer,intent(in) :: lmax,lpot,mpot
    complex(dp),intent(in) :: phase_plus,phase_minus
    complex(dp),intent(out) :: ag((lmax+1)**2,(lmax+1)**2), &
         ar((lmax+1)**2,(lmax+1)**2,6),aa((lmax+1)**2,(lmax+1)**2,6)
    integer :: nt,np,nq,nh,lymax,it,ip,q,l,m,idx,jdx,v,t,iq,a,b
    real(dp) :: phi,pi,nvec(3)
    real(dp),allocatable :: x(:),w(:),xyz(:,:),qw(:)
    complex(dp),allocatable :: y(:,:),dy(:,:,:),weighted(:,:),yl(:),vp(:)
    complex(dp) :: ca(3,-1:1,-1:1)
    pi=acos(-1._dp); nh=(lmax+1)**2
    ! Products have degree at most 2*lmax+lpot+2.  The extra margin
    ! also covers top-l angular derivatives; no lmax+1 intermediate is cut.
    nt=2*lmax+lpot+4; np=2*(2*lmax+lpot+2)+1; nq=nt*np
    lymax=max(lmax+1,lpot)
    allocate(x(nt),w(nt),xyz(nq,3),qw(nq),y(nq,nh),dy(nq,nh,6), &
         weighted(nq,nh),yl((lymax+1)**2),vp(nq))
    call gauss_legendre_stress(nt,x,w)
    dy=zero; q=0
    do it=1,nt
       do ip=1,np
          q=q+1; phi=2*pi*(ip-1)/np
          nvec=[sqrt(max(0._dp,1-x(it)**2))*cos(phi), &
                sqrt(max(0._dp,1-x(it)**2))*sin(phi),x(it)]
          xyz(q,:)=nvec; qw(q)=w(it)*2*pi/np
          call YLM(nvec,lymax,yl,(lymax+1)**2)
          y(q,:)=yl(:nh)
          vp(q)=phase_plus*yl(lpot*(lpot+1)+mpot+1)
          if(mpot/=0) vp(q)=vp(q)+phase_minus*yl(lpot*(lpot+1)-mpot+1)
          do l=0,lmax
             do m=-l,l
                idx=l*(l+1)+m+1
                call c_factor(l,m,ca)
                do v=-1,1,2
                   if(l+v<0) cycle
                   do t=-1,1
                      if(abs(m+t)>l+v) cycle
                      jdx=(l+v)*(l+v+1)+m+t+1
                      do iq=1,6
                         a=axis_a(iq); b=axis_b(iq)
                         dy(q,idx,iq)=dy(q,idx,iq)+ &
                              0.5_dp*(nvec(a)*ca(b,v,t)+nvec(b)*ca(a,v,t))*yl(jdx)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    do idx=1,nh
       weighted(:,idx)=y(:,idx)*vp*qw
    enddo
    call zgemm('C','N',nh,nh,nq,one,y,nq,weighted,nq,zero,ag,nh)
    do iq=1,6
       call zgemm('C','N',nh,nh,nq,one,dy(:,:,iq),nq,weighted,nq,zero,aa(:,:,iq),nh)
    enddo
    do iq=1,6
       do idx=1,nh
          weighted(:,idx)=y(:,idx)*vp*qw*xyz(:,axis_a(iq))*xyz(:,axis_b(iq))
       enddo
       call zgemm('C','N',nh,nh,nq,one,y,nq,weighted,nq,zero,ar(:,:,iq),nh)
    enddo
    ! Numerical integration of exactly vanishing coefficients produces tiny
    ! noise; remove only roundoff, not physically small radial amplitudes.
    where(abs(ag)<2.e-13_dp) ag=zero
    where(abs(ar)<2.e-13_dp) ar=zero
    where(abs(aa)<2.e-13_dp) aa=zero
    deallocate(x,w,xyz,qw,y,dy,weighted,yl,vp)
  end subroutine
end module
module fact
real*8 :: fct(100)
end module
module com_mpi
logical :: Qprint=.false.
integer :: myrank=1,master=0
end module
module dmft_stress
real*8 :: stress_norm_response_sign=1d0
end module
module param
integer,parameter :: lmax2n=8
end module
module defs
complex*16,parameter :: zeroc=(0d0,0d0)
end module
module w_atpar
type cache_t
 logical :: ready=.false.
 integer :: nterm=0
 integer,allocatable :: ih(:),ibra(:),iket(:)
 complex*16,allocatable :: radial(:,:,:),angular(:,:,:)
end type
type(cache_t) :: nsp_geometry(1)
complex*16 :: VPts(2,2,1,1),WPts(2,2,1,1)
! Match production one-based bounds, including the first-record test.
integer :: lpv(1,1)=0,lv(1,1)=0,llv(1,1)=0,mmv(1,1)=0,Nhmx(1)
end module
SUBROUTINE BuildWindowResponseMatrices(rhoResponse,eResponse, &
     projectorGradient,DMFTUfull,Aweight,AEweight,E,nemin,nemax, &
     DM_nemin,DM_nemax,nbands_full,nbands,maxdim2,norbitals,maxdim, &
     ncix,cix_orb,cixdim,nindo,iSx,Sigind,gradientDirect,useDirectC)
  USE com_mpi, ONLY: Qprint
  IMPLICIT NONE
  INTEGER, intent(in) :: nemin, nemax, DM_nemin, DM_nemax
  INTEGER, intent(in) :: nbands_full, nbands, maxdim2, norbitals, maxdim, ncix
  INTEGER, intent(in) :: cix_orb(norbitals), cixdim(ncix), nindo(norbitals)
  INTEGER, intent(in) :: iSx(maxdim2,norbitals), Sigind(maxdim,maxdim,ncix)
  REAL*8, intent(in) :: E(*)
  COMPLEX*16, intent(in) :: DMFTUfull(nbands_full,maxdim2,norbitals)
  COMPLEX*16, intent(in) :: Aweight(nbands,nbands), AEweight(nbands,nbands)
  COMPLEX*16, intent(out) :: rhoResponse(nbands_full,nbands_full)
  COMPLEX*16, intent(out) :: eResponse(nbands_full,nbands_full)
  COMPLEX*16, intent(out) :: projectorGradient(nbands,maxdim2,norbitals)
  ! Direct C retains components in ker(U) which cannot be reconstructed
  ! from B=U*C.  The pseudoinverse below is retained as a diagnostic/fallback
  ! for full-rank test callers, not as a definition of the gradient.
  LOGICAL,intent(in) :: useDirectC
  COMPLEX*16,intent(in) :: gradientDirect(nbands,maxdim2,norbitals)
  INTEGER :: icix, ip, iq, iorb, ind, ia, iw, in, iw0, iw1, ja
  INTEGER :: nactive, icolumn, info, lwork, lrwork, liwork, nsmall, ndegenerate
  INTEGER, allocatable :: activeColumn(:,:), iwork(:)
  REAL*8 :: eigmax, cutoff, denominator, epsn, responseNorm
  REAL*8 :: sourceNorm, reconstructionResidual, smallestKept
  REAL*8, allocatable :: seig(:), rwork(:)
  LOGICAL :: projectorDiagnostic
  LOGICAL :: includeWithinWindow
  CHARACTER(LEN=64) :: diagnosticEnvironment
  INTEGER :: diagnosticStatus
  COMPLEX*16 :: response
  COMPLEX*16, allocatable :: Uall(:,:), Uwindow(:,:), Svec(:,:), Sinv(:,:)
  COMPLEX*16, allocatable :: B(:,:), reconstructedB(:,:), localPart(:,:), &
       Kfull(:,:), work(:)
  COMPLEX*16, parameter :: cone=(1.d0,0.d0), czero=(0.d0,0.d0)

  rhoResponse=czero
  eResponse=czero
  projectorGradient=czero
  diagnosticEnvironment=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_PROJECTOR_DIAG', &
       diagnosticEnvironment,STATUS=diagnosticStatus)
  projectorDiagnostic=diagnosticStatus.eq.0 .and. &
       len_trim(diagnosticEnvironment).gt.0 .and. &
       trim(diagnosticEnvironment).ne.'0'
  diagnosticEnvironment=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_MOVING_BAND', &
       diagnosticEnvironment,STATUS=diagnosticStatus)
  includeWithinWindow=diagnosticStatus.eq.0 .and. &
       len_trim(diagnosticEnvironment).gt.0 .and. &
       trim(diagnosticEnvironment).ne.'0'
  iw0=DM_nemin-nemin+1
  iw1=DM_nemax-nemin+1
  if (iw0.lt.1 .or. iw1.gt.nbands_full .or. iw1-iw0+1.ne.nbands) then
     WRITE(6,*) 'ERROR: inconsistent fixed DMFT window in BuildWindowResponseMatrices', &
          nemin,nemax,DM_nemin,DM_nemax,nbands_full,nbands
     STOP 'BuildWindowResponseMatrices: inconsistent band bounds'
  endif

  ! Keep every local state which occurs in at least one self-energy matrix
  ! element.  Usually testing Sigind(ip,ip) is sufficient, but the row/column
  ! test also covers a purely off-diagonal block.
  allocate(activeColumn(maxdim,ncix))
  activeColumn=0
  nactive=0
  do icix=1,ncix
     do ip=1,cixdim(icix)
        if (any(Sigind(ip,1:cixdim(icix),icix).gt.0) .or. &
            any(Sigind(1:cixdim(icix),ip,icix).gt.0)) then
           nactive=nactive+1
           activeColumn(ip,icix)=nactive
        endif
     enddo
  enddo
  if (nactive.eq.0) then
     deallocate(activeColumn)
     return
  endif

  allocate(Uall(nbands_full,nactive),Uwindow(nbands,nactive))
  Uall=czero
  do iorb=1,norbitals
     icix=cix_orb(iorb)
     if (icix.eq.0) cycle
     do ind=1,nindo(iorb)
        ip=iSx(ind,iorb)
        icolumn=activeColumn(ip,icix)
        if (icolumn.gt.0) Uall(:,icolumn)=DMFTUfull(:,ind,iorb)
     enddo
  enddo
  Uwindow=Uall(iw0:iw1,:)

  allocate(Svec(nactive,nactive),Sinv(nactive,nactive),seig(nactive))
  call zgemm('C','N',nactive,nactive,nbands,cone,Uwindow,nbands,Uwindow,nbands, &
       czero,Svec,nactive)
  lwork=2*nactive+nactive*nactive
  lrwork=5*nactive+2*nactive*nactive+1
  liwork=3+5*nactive
  allocate(work(lwork),rwork(lrwork),iwork(liwork))
  CALL ZHEEVD('V','U',nactive,Svec,nactive,seig,work,lwork,rwork,lrwork,iwork,liwork,info)
  if (info.ne.0) then
     WRITE(6,*) 'ERROR: projector metric diagonalization failed in window response',info
     STOP 'BuildWindowResponseMatrices: ZHEEVD failed'
  endif

  eigmax=maxval(seig)
  cutoff=max(1.d-13,1.d-11*eigmax)
  smallestKept=huge(1.d0)
  Sinv=Svec
  nsmall=0
  do ip=1,nactive
     if (seig(ip).gt.cutoff) then
        Sinv(:,ip)=Sinv(:,ip)/seig(ip)
        smallestKept=min(smallestKept,seig(ip))
     else
        Sinv(:,ip)=czero
        nsmall=nsmall+1
     endif
  enddo
  Sinv=matmul(Sinv,conjg(transpose(Svec)))

  allocate(B(nbands,nbands),localPart(nactive,nbands),Kfull(nbands_full,nbands))
  do iw=1,nbands
     do in=1,nbands
        ! eps multiplies the Green's function from the left.
        B(iw,in)=AEweight(iw,in)-E(DM_nemin+iw-1)*Aweight(iw,in)
     enddo
  enddo
  call zgemm('C','N',nactive,nbands,nbands,cone,Uwindow,nbands,B,nbands, &
       czero,localPart,nactive)
  localPart=matmul(Sinv,localPart)
  if (projectorDiagnostic .and. Qprint) then
     ! B=rhoE-epsilon*rho must lie in the column space of the normalized
     ! projector U_W if it is to be represented as U_W C.  This check is
     ! independent of the subsequent Sternheimer denominator and therefore
     ! cleanly distinguishes an algebra/convention error in C from an
     ! incomplete complementary-band sum.
     allocate(reconstructedB(nbands,nbands))
     reconstructedB=matmul(Uwindow,localPart)
     sourceNorm=sqrt(sum(abs(B)**2))
     reconstructionResidual=sqrt(sum(abs(B-reconstructedB)**2))/ &
          max(sourceNorm,tiny(1.d0))
     WRITE(6,'(A,4ES18.8,I6)') ':DMFT_PROJECTOR_RECONSTRUCTION: ', &
          sourceNorm,reconstructionResidual,eigmax,smallestKept,nsmall
     deallocate(reconstructedB)
  endif
  if (useDirectC) then
     if (projectorDiagnostic .and. Qprint) then
        sourceNorm=0.d0
        reconstructionResidual=0.d0
        do iorb=1,norbitals
           icix=cix_orb(iorb)
           if (icix.le.0) cycle
           do ind=1,nindo(iorb)
              icolumn=activeColumn(iSx(ind,iorb),icix)
              if (icolumn.le.0) cycle
              sourceNorm=sourceNorm+sum(abs(gradientDirect(:,ind,iorb))**2)
              reconstructionResidual=reconstructionResidual+sum(abs( &
                   gradientDirect(:,ind,iorb)-localPart(icolumn,:))**2)
           enddo
        enddo
        WRITE(6,'(A,2ES22.12,I6)') ':DMFT_C_DIRECT_VS_PINV: ', &
             sqrt(sourceNorm),sqrt(reconstructionResidual/max(sourceNorm,tiny(1.d0))),nsmall
     endif
     do iorb=1,norbitals
        icix=cix_orb(iorb)
        if (icix.le.0) cycle
        do ind=1,nindo(iorb)
           icolumn=activeColumn(iSx(ind,iorb),icix)
           if (icolumn.gt.0) localPart(icolumn,:)=gradientDirect(:,ind,iorb)
        enddo
     enddo
  else if (nsmall.gt.0) then
     WRITE(6,*) 'ERROR: rank-deficient projector requires the direct functional gradient C.'
     STOP 'A pseudoinverse cannot recover the null-space part of C'
  endif

  ! localPart(m,n)=sum_omega [Sigma U^+ G]_{m n}, including the
  ! irreducible-k weight.  It is the one-particle gradient conjugate to the
  ! explicit projector derivative: dGamma_proj=2 Re sum_mn C_mn dU_nm.
  ! This is an energy derivative. l2main converts it to raw STR_VAL
  ! with a minus sign; mixer reports +dGamma/d(strain)/Omega.
  ! (The corresponding atomic force has the opposite sign.)
  do iorb=1,norbitals
     icix=cix_orb(iorb)
     if (icix.eq.0) cycle
     do ind=1,nindo(iorb)
        ip=iSx(ind,iorb)
        icolumn=activeColumn(ip,icix)
        if (icolumn.gt.0) projectorGradient(:,ind,iorb)= &
             localPart(icolumn,:)
     enddo
  enddo
  call zgemm('N','N',nbands_full,nbands,nactive,cone,Uall,nbands_full, &
       localPart,nactive,czero,Kfull,nbands_full)

  ndegenerate=0
  do ia=1,nbands_full
     if (ia.ge.iw0 .and. ia.le.iw1) cycle
     do iw=1,nbands
        in=iw0+iw-1
        epsn=E(DM_nemin+iw-1)
        denominator=epsn-E(nemin+ia-1)
        if (abs(denominator).le.1.d-10) then
           ndegenerate=ndegenerate+1
           cycle
        endif
        response=Kfull(ia,iw)/denominator
        ! U_nl=<psi_n|phi_l>; if |psi_n>'=sum_a |psi_a> R_an,
        ! then U'_nl=sum_a conjg(R_an)*U_al.  Thus K_an multiplies
        ! conjg(H'_an-e_n O'_an), and belongs in the (a,n) block.
        ! Reversing these two entries violates arbitrary KS-phase invariance.
        rhoResponse(ia,in)=response
        rhoResponse(in,ia)=conjg(response)
        eResponse(ia,in)=epsn*response
        eResponse(in,ia)=conjg(eResponse(ia,in))
     enddo
  enddo

  if (includeWithinWindow) then
     ! Complete the band connection inside the fixed-index window without
     ! dividing two individually large terms by a small band splitting.
     ! Since B=rhoE-epsilon*rho and rhoE is Hermitian, every pair obeys
     !
     !   B_an-B_na^*=(epsilon_n-epsilon_a) rho_an,
     !
     ! while the energy-weighted ordered-pair sum is exactly rhoE_an.  The
     ! former implementation reconstructed these identities from two
     ! Kfull/denominator terms.  That is algebraically correct away from a
     ! degeneracy but catastrophically ill conditioned in FeSe: each term
     ! reached hundreds of GPa before cancellation.  Aweight and AEweight
     ! are the same already frequency-summed Hermitian matrices, so insert
     ! their off-diagonal blocks directly.  This is gauge covariant inside a
     ! degenerate block and makes the MOVING_BAND representation an exact
     ! diagnostic identity rather than a near-degenerate perturbation solve.
     do iw=1,nbands
        in=iw0+iw-1
        do ja=1,nbands
           if (ja.eq.iw) cycle
           ia=iw0+ja-1
           rhoResponse(ia,in)=rhoResponse(ia,in)+Aweight(ja,iw)
           eResponse(ia,in)=eResponse(ia,in)+AEweight(ja,iw)
        enddo
        ! Generalized parallel transport has
        ! A_n^+ O A_{n,q}=-O_{q,nn}/2, so the direct Klein derivative is
        ! -Re(K_nn) O_{q,nn}.  CMP_DMFT_STRESS, however, accumulates WIEN's
        ! intermediate STR_VAL, whose sign is the negative of that direct
        ! derivative because mixer finally forms -STR_VAL.  Insert
        ! +Herm(K_nn), not -Herm(K_nn), in the energy-response matrix.
        ! K_nn should be real at stationarity; taking its real part prevents
        ! insignificant Matsubara roundoff from entering a real stress.
        eResponse(in,in)=eResponse(in,in)+dble(Kfull(in,iw))
     enddo
  endif

  responseNorm=sqrt(sum(abs(rhoResponse)**2))
  if (Qprint .and. (nsmall.gt.0 .or. ndegenerate.gt.0)) then
     WRITE(6,'(A,2I6,A,ES12.4)') 'Window-response pseudoinverse: small/degenerate=', &
          nsmall,ndegenerate,'  ||rho_resp||=',responseNorm
  endif

  deallocate(activeColumn,Uall,Uwindow,Svec,Sinv,seig,work,rwork,iwork)
  deallocate(B,localPart,Kfull)
END SUBROUTINE BuildWindowResponseMatrices



SUBROUTINE ComputeNormalizationDerivative(dXrel,dS,X,Sigind,cixdim, &
     maxdim,ncix,SIMPLE,nq)
  IMPLICIT NONE
  INTEGER, intent(in) :: maxdim,ncix,nq,cixdim(ncix), &
       Sigind(maxdim,maxdim,ncix)
  LOGICAL, intent(in) :: SIMPLE
  COMPLEX*16, intent(in)  :: dS(maxdim,maxdim,ncix,nq), &
       X(maxdim,maxdim,ncix)
  COMPLEX*16, intent(out) :: dXrel(maxdim,maxdim,ncix,nq)
  INTEGER :: icix, cixdm, cixdms, ip, iq, idir, info
  INTEGER :: lwork, lrwork, liwork
  INTEGER, allocatable :: cini(:), iwork(:)
  REAL*8, allocatable :: xeig(:), rwork(:)
  COMPLEX*16, allocatable :: Xvec(:,:), dSt(:,:), Bt(:,:), Bout(:,:), work(:)

  dXrel=0.d0
  if (SIMPLE) then
     do icix=1,ncix
        do ip=1,cixdim(icix)
           if (Sigind(ip,ip,icix).eq.0) cycle
           do idir=1,nq
              ! X=S^{-1/2}: X^{-1} X_A = -1/2 X^2 S_A.
              dXrel(ip,ip,icix,idir) = -0.5d0*X(ip,ip,icix)**2*dS(ip,ip,icix,idir)
           enddo
        enddo
     enddo
     return
  endif

  do icix=1,ncix
     cixdm=cixdim(icix)
     allocate(cini(cixdm))
     cixdms=0
     do ip=1,cixdm
        if (Sigind(ip,ip,icix).ne.0) then
           cixdms=cixdms+1
           cini(cixdms)=ip
        endif
     enddo
     if (cixdms.eq.0) then
        deallocate(cini)
        cycle
     endif

     lwork=4*cixdms+cixdms*cixdms
     lrwork=10*cixdms+2*cixdms*cixdms+1
     liwork=3+10*cixdms
     allocate(Xvec(cixdms,cixdms), dSt(cixdms,cixdms), Bt(cixdms,cixdms), Bout(cixdms,cixdms))
     allocate(xeig(cixdms), work(lwork), rwork(lrwork), iwork(liwork))
     do ip=1,cixdms
        do iq=1,cixdms
           Xvec(ip,iq)=X(cini(ip),cini(iq),icix)
        enddo
     enddo
     CALL ZHEEVD('V','U',cixdms,Xvec,cixdms,xeig,work,lwork,rwork,lrwork,iwork,liwork,info)
     if (info.ne.0 .or. minval(xeig).le.1.d-14) then
        write(6,*) 'ERROR: diagonalization of S^{-1/2} failed in normalization derivative',icix,info
        STOP 'ComputeNormalizationDerivative failed'
     endif

     do idir=1,nq
        do ip=1,cixdms
           do iq=1,cixdms
              Bout(ip,iq)=dS(cini(ip),cini(iq),icix,idir)
           enddo
        enddo
        dSt=matmul(conjg(transpose(Xvec)),matmul(Bout,Xvec))
        do ip=1,cixdms
           do iq=1,cixdms
              ! Solve Y X_A+X_A Y=-X S_A X with Y=X^{-1}, then
              ! return the relative derivative B_A=X^{-1}X_A.  The actual
              ! normalized-projector derivative is Ubar_A=...+Ubar B_A.
              Bt(ip,iq)=-xeig(iq)*dSt(ip,iq)/(1.d0/xeig(ip)+1.d0/xeig(iq))
           enddo
        enddo
        Bout=matmul(Xvec,matmul(Bt,conjg(transpose(Xvec))))
        do ip=1,cixdms
           do iq=1,cixdms
              dXrel(cini(ip),cini(iq),icix,idir)=Bout(ip,iq)
           enddo
        enddo
     enddo
     deallocate(Xvec,dSt,Bt,Bout,xeig,work,rwork,iwork,cini)
  enddo
END SUBROUTINE ComputeNormalizationDerivative





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





SUBROUTINE BuildNormalizationWindowResponseMatrices(rhoResponse,eResponse, &
     normalizationGradient,Pfull,Qbar,wg,E,nemin,nemax,DM_nemin,DM_nemax,nbands_full,nbands, &
     maxdim2,norbitals,maxdim,ncix,cix_orb,cixdim,nindo,iSx)
  USE com_mpi,ONLY: Qprint
  USE dmft_stress,ONLY: stress_norm_response_sign
  IMPLICIT NONE
  INTEGER,intent(in) :: nemin,nemax,DM_nemin,DM_nemax,nbands_full,nbands
  INTEGER,intent(in) :: maxdim2,norbitals,maxdim,ncix
  INTEGER,intent(in) :: cix_orb(norbitals),cixdim(ncix),nindo(norbitals), &
       iSx(maxdim2,norbitals)
  REAL*8,intent(in) :: wg,E(*)
  COMPLEX*16,intent(in) :: Pfull(nbands_full,maxdim2,norbitals), &
       Qbar(maxdim,maxdim,ncix)
  COMPLEX*16,intent(out) :: rhoResponse(nbands_full,nbands_full), &
       eResponse(nbands_full,nbands_full), &
       normalizationGradient(nbands,maxdim2,norbitals)
  INTEGER :: iw0,iw1,icix,cixdm,iorb,ind,ia,iw,in,ndegenerate, &
       metricStatus
  REAL*8 :: denominator,epsn,responseNorm
  COMPLEX*16 :: response
  COMPLEX*16,allocatable :: Pall(:,:),Pwin(:,:),Gbar(:,:),Kfull(:,:)
  CHARACTER(len=64) :: metricEnvironment
  LOGICAL :: includeWindowMetric,includeWithinWindow

  rhoResponse=0.d0; eResponse=0.d0; normalizationGradient=0.d0
  iw0=DM_nemin-nemin+1
  iw1=DM_nemax-nemin+1
  if (iw0.lt.1 .or. iw1.gt.nbands_full .or. iw1-iw0+1.ne.nbands) &
       STOP 'BuildNormalizationWindowResponseMatrices: band bounds'
  allocate(Kfull(nbands_full,nbands))
  Kfull=0.d0

  ! Metric contribution omitted by the original
  ! complementary-band construction.  The generalized LAPW eigenvectors
  ! obey A_W^+ O A_W=1.  In parallel-transport gauge their within-window
  ! response consequently contains
  !
  !       A_W^+ O A_{W,q} = -1/2 O_{q,WW}.
  !
  ! For P=A_W^+ Q this produces the gauge-invariant contribution
  !
  !       S_q(metric) = - sum_k w_k P_k^+ O_{q,WW} P_k
  !
  ! to S=sum_k w_k P_k^+P_k and hence to X=S^(-1/2).  The historical
  ! routine below explicitly skips every row inside W, so it cannot create
  ! this term.  Kfull is the adjoint band-space kernel for the normalization
  ! response.  Notice the sign: the generalized
  ! eigenvectors satisfy
  !
  !       delta P_W = -1/2 O_{q,WW} P_W
  !
  ! in parallel-transport gauge.  Hence
  !
  !       2 Re sum(Gbar delta P_W)
  !         = -Tr[Herm(P_W Gbar^T) O_{q,WW}].
  !
  ! CMP_DMFT_STRESS writes the WIEN intermediate STR_VAL, not the final
  ! physical stress.  Its energy-response argument contributes
  ! +Tr[eResponse O_q] to STR_VAL, and mixer subsequently negates STR_VAL.
  ! Therefore the direct Klein contribution
  ! -Tr[Herm(K_WW) O_q] requires eResponse=+Herm(K_WW) here.  Using the
  ! opposite sign confuses the direct-functional convention with WIEN's
  ! stored valence convention and applies the mixer sign twice.
  !
  ! Include this connection in the fixed-coefficient derivative by default.
  ! Independent generalized-eigenproblem and basis-change tests verify its
  ! sign. It would be double-counted only if P' already included the moving
  ! orthonormal-frame metric response; cmp_overlap does not include it.
  metricEnvironment=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_NORM_WINDOW_METRIC', &
       metricEnvironment,STATUS=metricStatus)
  ! P' is constructed at fixed LAPW coefficients.  Its nonorthogonal
  ! window metric must be differentiated as well: -P^dagger O'_WW P.
  ! Set the diagnostic explicitly to 0 only to reproduce the old omission.
  includeWindowMetric=.true.
  if (metricStatus.eq.0 .and. trim(metricEnvironment).eq.'0') &
       includeWindowMetric=.false.
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_MOVING_BAND', &
       metricEnvironment,STATUS=metricStatus)
  ! A unitary rotation inside W cancels identically from
  ! S=sum_k P_k^+P_k.  The normalization derivative therefore needs the
  ! generalized-metric connection, but no separate within-W Sternheimer
  ! rotation.  (The fixed-X projector branch handles its own band
  ! connection.)  Keeping these notions separate prevents counting the
  ! same within-window rotation twice in MOVING_BAND diagnostics.
  includeWithinWindow=.false.
  if (metricStatus.eq.0 .and. len_trim(metricEnvironment).gt.0 .and. &
       trim(metricEnvironment).ne.'0') includeWindowMetric=.true.

  do icix=1,ncix
     cixdm=cixdim(icix)
     allocate(Pall(nbands_full,cixdm),Pwin(nbands,cixdm), &
          Gbar(nbands,cixdm))
     Pall=0.d0
     do iorb=1,norbitals
        if (cix_orb(iorb).ne.icix) cycle
        do ind=1,nindo(iorb)
           Pall(:,iSx(ind,iorb))=Pfull(:,ind,iorb)
        enddo
     enddo
     Pwin=Pall(iw0:iw1,:)
     ! Qbar=conjg(Q), hence conjg(Pwin)*Qbar=conjg(Pwin*Q).
     ! The factor two follows from differentiating both factors in P^+P;
     ! wg is exactly the normalization weight used in cmp_overlap.
     Gbar=2.d0*wg*matmul(conjg(Pwin),Qbar(1:cixdm,1:cixdm,icix))
     ! Store the conjugated gradient in the same sparse orbital layout used
     ! by the main projector gradient.  Its direct contraction is
     !
     !   2 Re sum(normalizationGradient * delta P),
     !
     ! where delta P is the raw-projector response.  This provides a central
     ! Sternheimer oracle which does not pass the normalization response
     ! through Blaha's rearranged DFT kernels.
     do iorb=1,norbitals
        if (cix_orb(iorb).ne.icix) cycle
        do ind=1,nindo(iorb)
           normalizationGradient(:,ind,iorb)=Gbar(:,iSx(ind,iorb))
        enddo
     enddo
     Kfull=Kfull+matmul(Pall,transpose(Gbar))
     deallocate(Pall,Pwin,Gbar)
  enddo

  if (includeWindowMetric) then
     eResponse(iw0:iw1,iw0:iw1)=eResponse(iw0:iw1,iw0:iw1)+ &
          0.5d0*(Kfull(iw0:iw1,:)+ &
          conjg(transpose(Kfull(iw0:iw1,:))))
     if (Qprint) WRITE(6,'(A,ES22.12)') &
          ':DMFT_NORM_WINDOW_METRIC_KERNEL: ', &
          sqrt(sum(abs(eResponse(iw0:iw1,iw0:iw1))**2))
  endif

  ndegenerate=0
  do ia=1,nbands_full
     if (ia.ge.iw0 .and. ia.le.iw1) cycle
     do iw=1,nbands
        in=iw0+iw-1
        epsn=E(DM_nemin+iw-1)
        denominator=epsn-E(nemin+ia-1)
        if (abs(denominator).le.1.d-10) then
           ndegenerate=ndegenerate+1
           cycle
        endif
        response=Kfull(ia,iw)/denominator
        ! Same covariant (a,n) orientation as the direct-projector response;
        ! the derivative of the projector bra conjugates R_an.
        rhoResponse(ia,in)=response
        rhoResponse(in,ia)=conjg(response)
        eResponse(ia,in)=epsn*response
        eResponse(in,ia)=conjg(eResponse(ia,in))
     enddo
  enddo
  if (includeWithinWindow) then
     do iw=1,nbands
        in=iw0+iw-1
        epsn=E(DM_nemin+iw-1)
        do ia=iw0,iw1
           if (ia.eq.in) cycle
           denominator=epsn-E(nemin+ia-1)
           if (abs(denominator).le.1.d-10) then
              ndegenerate=ndegenerate+1
              cycle
           endif
           response=Kfull(ia,iw)/denominator
           rhoResponse(ia,in)=rhoResponse(ia,in)+response
           rhoResponse(in,ia)=rhoResponse(in,ia)+conjg(response)
           eResponse(ia,in)=eResponse(ia,in)+epsn*response
           eResponse(in,ia)=eResponse(in,ia)+conjg(epsn*response)
        enddo
     enddo
  endif

  ! These are physical H'/O' response matrices. CMP_DMFT_STRESS already
  ! supplies the WIEN intermediate minus sign and mixer reverses it once.
  ! Consequently the physical response needs no further sign reversal.
  ! The multiplier is retained only for explicit legacy diagnostics.
  rhoResponse=stress_norm_response_sign*rhoResponse
  eResponse=stress_norm_response_sign*eResponse
  responseNorm=sqrt(sum(abs(rhoResponse)**2))
  if (Qprint .and. ndegenerate.gt.0) WRITE(6,'(A,I6,A,ES12.4)') &
       'Normalization-window response: degenerate=',ndegenerate, &
       '  ||rho_resp||=',responseNorm
  deallocate(Kfull)
END SUBROUTINE BuildNormalizationWindowResponseMatrices


SUBROUTINE AccumulateDirectProjectorPair(CPair,G,Uactive,localSigma, &
     cixdim,nbands,maxdim,ncix,UdagG,UdagGminus,product,productMinus)
  IMPLICIT NONE
  INTEGER,intent(in) :: nbands,maxdim,ncix,cixdim(ncix)
  COMPLEX*16,intent(in) :: G(nbands,nbands), &
       Uactive(nbands,maxdim,ncix),localSigma(maxdim,maxdim,ncix)
  COMPLEX*16,intent(inout) :: CPair(maxdim,nbands,ncix)
  COMPLEX*16,intent(inout) :: UdagG(maxdim,nbands), &
       UdagGminus(maxdim,nbands),product(maxdim,nbands), &
       productMinus(maxdim,nbands)
  INTEGER :: icix,nloc
  COMPLEX*16,parameter :: cone=(1.d0,0.d0),czero=(0.d0,0.d0)

  do icix=1,ncix
     nloc=cixdim(icix)
     if (nloc.le.0) cycle
     CALL zgemm('C','N',nloc,nbands,nbands,cone, &
          Uactive(1,1,icix),nbands,G,nbands,czero,UdagG,maxdim)
     CALL zgemm('N','N',nloc,nbands,nloc,cone, &
          localSigma(1,1,icix),maxdim,UdagG,maxdim,czero,product,maxdim)

     CALL zgemm('C','C',nloc,nbands,nbands,cone, &
          Uactive(1,1,icix),nbands,G,nbands,czero,UdagGminus,maxdim)
     CALL zgemm('C','N',nloc,nbands,nloc,cone, &
          localSigma(1,1,icix),maxdim,UdagGminus,maxdim,czero, &
          productMinus,maxdim)
     CPair(1:nloc,:,icix)=CPair(1:nloc,:,icix)+ &
          product(1:nloc,:)+productMinus(1:nloc,:)
  enddo
END SUBROUTINE AccumulateDirectProjectorPair

SUBROUTINE AccumulateDirectProjectorSpectral(CPair,G,Uactive,localSigma, &
     cixdim,nbands,maxdim,ncix,UdagG,UdagGminus,product,productMinus)
  IMPLICIT NONE
  INTEGER,intent(in) :: nbands,maxdim,ncix,cixdim(ncix)
  COMPLEX*16,intent(in) :: G(nbands,nbands), &
       Uactive(nbands,maxdim,ncix),localSigma(maxdim,maxdim,ncix)
  COMPLEX*16,intent(inout) :: CPair(maxdim,nbands,ncix)
  COMPLEX*16,intent(inout) :: UdagG(maxdim,nbands), &
       UdagGminus(maxdim,nbands),product(maxdim,nbands), &
       productMinus(maxdim,nbands)
  INTEGER :: icix,nloc
  COMPLEX*16,parameter :: cone=(1.d0,0.d0),czero=(0.d0,0.d0)

  do icix=1,ncix
     nloc=cixdim(icix)
     if (nloc.le.0) cycle
     CALL zgemm('C','N',nloc,nbands,nbands,cone, &
          Uactive(1,1,icix),nbands,G,nbands,czero,UdagG,maxdim)
     CALL zgemm('N','N',nloc,nbands,nloc,cone, &
          localSigma(1,1,icix),maxdim,UdagG,maxdim,czero,product,maxdim)

     CALL zgemm('C','C',nloc,nbands,nbands,cone, &
          Uactive(1,1,icix),nbands,G,nbands,czero,UdagGminus,maxdim)
     CALL zgemm('C','N',nloc,nbands,nloc,cone, &
          localSigma(1,1,icix),maxdim,UdagGminus,maxdim,czero, &
          productMinus,maxdim)
     CPair(1:nloc,:,icix)=CPair(1:nloc,:,icix)+ &
          -(product(1:nloc,:)-productMinus(1:nloc,:))/cmplx(0.d0,2.d0*acos(-1.d0),8)
  enddo
END SUBROUTINE AccumulateDirectProjectorSpectral



  subroutine c_alpha_m(c_matrix)
     implicit none
     complex*16         :: c_matrix(3,-1:1)
     complex*16, parameter :: imag  = (0.D0,1.d0)
     complex*16, parameter :: zeroc = (0.d0,0.d0)
     real*8, parameter     :: pi    = acos(-1.d0)
     real*8, parameter     :: two   = 2.D0
     real*8, parameter     :: three = 3.d0
     real*8, parameter     :: twopi = 2.D0*pi
    
       c_matrix(1,-1) = sqrt(twopi/three)
       c_matrix(1,0)  = zeroc
       c_matrix(1,1)  = -sqrt(twopi/three)

       c_matrix(2,-1) = imag*sqrt(twopi/three)    !our work doesn't have minus, kuppel does
       c_matrix(2,0)  = zeroc
       c_matrix(2,1)  = imag*sqrt(twopi/three)    !our work doesn't have minus, kuppel does

       c_matrix(3,-1) = zeroc
       c_matrix(3,0)  = sqrt(two*twopi/three)
       c_matrix(3,1)  = zeroc
  
  RETURN
  end subroutine c_alpha_m
  
      
  
      
  subroutine C_factor(ll,mm, c_ast_lm)
     !=====subroutine to calculate c_ast(l,m) which being used to calculate derivative of the spherical hermonics==========     
     implicit none
     integer,    intent(in)   :: ll,mm     
     complex*16, intent(out)  :: c_ast_lm(3,-1:1,-1:1)
     complex*16, parameter    :: zeroc = (0.d0,0.d0)
     complex*16, parameter    :: imag  = (0.D0,1.D0)
     integer,    parameter    :: one = 1.D0, two = 2.D0, three = 3.D0, f4r = 4.D0      
     integer                  :: index1, ss, tt, index, alpha, beta
     real*8                   :: delta(3,3), C1, C2, FF(6)

                    
     C1  = (ll-mm+one)*(ll-mm+two)
     C2  = f4r*(two*ll+one)*(two*ll+three)
     FF(1) = sqrt(C1/C2)

     C1 = (ll+mm-one)*(ll+mm)
     C2 = f4r*(two*ll-one)*(two*ll+one)
     FF(2) = sqrt(C1/C2)

     C1 = (ll+mm+one)*(ll+mm+two)
     C2 = f4r*(two*ll+one)*(two*ll+three)
     FF(3) = sqrt(C1/C2)

     C1 = (ll-mm-one)*(ll-mm)
     C2 = f4r*(two*ll-one)*(two*ll+one)
     FF(4) = sqrt(C1/C2)

     C1 = (ll+one)*(ll+one) - mm*mm
     C2 = (two*ll+one)*(two*ll+three)
     FF(5) = sqrt(C1/C2)

     C1 = ll*ll - mm*mm
     C2 = (two*ll-one)*(two*ll+one)
     FF(6) = sqrt(C1/C2)
     
                                                  
     delta  = 0
     DO alpha = 1, 3
        DO beta = 1, 3
           if(alpha .eq. beta) then
               delta(alpha,beta) = 1.d0
           else
               delta(alpha,beta) = 0.d0
           endif
        ENDDO
     ENDDO
     
     c_ast_lm(:,0,:) = zeroc
         
     DO alpha = 1, 3
           
        c_ast_lm(alpha,-1,-1) = -(ll+1)*(delta(alpha,1)+imag*delta(alpha,2))*FF(2)
        c_ast_lm(alpha,-1, 0) =  (ll+1)*delta(alpha,3)*FF(6)
        c_ast_lm(alpha,-1, 1) =  (ll+1)*(delta(alpha,1)-imag*delta(alpha,2))*FF(4)
              
        c_ast_lm(alpha, 1,-1) = -ll*(delta(alpha,1)+imag*delta(alpha,2))*FF(1)
        c_ast_lm(alpha, 1, 0) = -ll*delta(alpha,3)*FF(5)
        c_ast_lm(alpha, 1, 1) =  ll*(delta(alpha,1) - imag*delta(alpha,2))*FF(3)
                           
     ENDDO     
  return
  end subroutine C_factor





SUBROUTINE build_nsp_geometry(jatom)
  ! Precompute the k-independent three-Gaunt tensors in the fifth and sixth
  ! lines of Eq. (G2).  The cache is sparse in (radial channel,m_bra,m_ket)
  ! and is owned by w_atpar together with the case.nsh matrices whose labels
  ! define it.
  use param, only: lmax2n
  use defs, only: zeroc
  use com_mpi, only: myrank,master
  use w_atpar, only: VPts,WPts,lpv,lv,llv,mmv,Nhmx,nsp_geometry
  implicit none
  REAL*8, parameter :: hf=0.5d0
  INTEGER, intent(in) :: jatom
  INTEGER :: ih,l,lp,mbra,mket,index1,index2,maxterm,iterm
  INTEGER :: alpha,beta,sp,tp,t,u,v,v1,lint,nu
  INTEGER, external :: NOTRI
  REAL*8, external :: GAUNT
  REAL*8 :: gaunt1,gaunt_pot,gaunt_couple
  REAL*8 :: channel_diff,channel_scale
  LOGICAL :: duplicate_radial_channel
  COMPLEX*16 :: c_am(3,-1:1),ca_st_lm(3,-1:1,-1:1)
  COMPLEX*16 :: radial_geom(3,3),angular_geom(3,3)
  COMPLEX*16 :: radial_factor,angular_factor

  IF (nsp_geometry(jatom)%ready) RETURN
  call c_alpha_m(c_am)

  ! First determine a strict upper bound.  The same radial channel is
  ! repeated consecutively in case.nsh for the original Hamiltonian-allowed
  ! magnetic pairs.  Its phase-weighted matrices and (l,L,M,l') labels are
  ! identical, so comparing with the preceding record removes the replica.
  maxterm=0
  DO ih=1,Nhmx(jatom)
     l=lpv(ih,jatom)
     lp=lv(ih,jatom)
     duplicate_radial_channel=.FALSE.
     IF (ih.GT.1) THEN
     IF (lpv(ih-1,jatom).EQ.l .AND. &
          lv(ih-1,jatom).EQ.lp .AND. &
          llv(ih-1,jatom).EQ.llv(ih,jatom) .AND. &
          mmv(ih-1,jatom).EQ.mmv(ih,jatom)) THEN
        channel_diff=SUM(ABS(VPts(:,:,ih-1,jatom)-VPts(:,:,ih,jatom)))+ &
             SUM(ABS(WPts(:,:,ih-1,jatom)-WPts(:,:,ih,jatom)))
        channel_scale=SUM(ABS(VPts(:,:,ih,jatom)))+ &
             SUM(ABS(WPts(:,:,ih,jatom)))
        duplicate_radial_channel= &
             channel_diff.LE.1.d-12*(1.d0+channel_scale)
     ENDIF
     ENDIF
     IF (.NOT.duplicate_radial_channel) &
          maxterm=maxterm+(2*l+1)*(2*lp+1)
  ENDDO

  allocate(nsp_geometry(jatom)%ih(maxterm))
  allocate(nsp_geometry(jatom)%ibra(maxterm))
  allocate(nsp_geometry(jatom)%iket(maxterm))
  allocate(nsp_geometry(jatom)%radial(3,3,maxterm))
  allocate(nsp_geometry(jatom)%angular(3,3,maxterm))
  nsp_geometry(jatom)%radial=zeroc
  nsp_geometry(jatom)%angular=zeroc

  iterm=0
  DO ih=1,Nhmx(jatom)
     l=lpv(ih,jatom)
     lp=lv(ih,jatom)
     duplicate_radial_channel=.FALSE.
     IF (ih.GT.1) THEN
     IF (lpv(ih-1,jatom).EQ.l .AND. &
          lv(ih-1,jatom).EQ.lp .AND. &
          llv(ih-1,jatom).EQ.llv(ih,jatom) .AND. &
          mmv(ih-1,jatom).EQ.mmv(ih,jatom)) THEN
        channel_diff=SUM(ABS(VPts(:,:,ih-1,jatom)-VPts(:,:,ih,jatom)))+ &
             SUM(ABS(WPts(:,:,ih-1,jatom)-WPts(:,:,ih,jatom)))
        channel_scale=SUM(ABS(VPts(:,:,ih,jatom)))+ &
             SUM(ABS(WPts(:,:,ih,jatom)))
        duplicate_radial_channel= &
             channel_diff.LE.1.d-12*(1.d0+channel_scale)
     ENDIF
     ENDIF
     IF (duplicate_radial_channel) CYCLE

     DO mbra=-l,l
        index1=l*(l+1)+mbra+1
        CALL c_factor(l,mbra,ca_st_lm)
        DO mket=-lp,lp
           index2=lp*(lp+1)+mket+1
           radial_geom=zeroc
           DO sp=0,2
              IF (NOTRI(sp,1,1).LT.0) CYCLE
              DO tp=-sp,sp
                 DO t=-1,1
                    DO u=-1,1
                       IF (-tp+t+u.NE.0) CYCLE
                       gaunt1=GAUNT(sp,1,1,tp,t,u)
                       IF (ABS(gaunt1).LT.1.d-30) CYCLE
                       DO v1=0,llv(ih,jatom)+lp
                          IF (NOTRI(v1,llv(ih,jatom),lp).LT.0) CYCLE
                          IF (NOTRI(l,sp,v1).LT.0) CYCLE
                          DO v=-v1,v1
                             IF (-v+mmv(ih,jatom)+mket.NE.0) CYCLE
                             gaunt_pot=GAUNT(v1,llv(ih,jatom),lp,v, &
                                  mmv(ih,jatom),mket)
                             IF (ABS(gaunt_pot).LT.1.d-30) CYCLE
                             IF (-mbra+tp+v.NE.0) CYCLE
                             gaunt_couple=GAUNT(l,sp,v1,mbra,tp,v)
                             IF (ABS(gaunt_couple).LT.1.d-30) CYCLE
                             DO beta=1,3
                                DO alpha=1,3
                                   radial_factor=c_am(alpha,t)* &
                                        c_am(beta,u)*gaunt1* &
                                        gaunt_pot*gaunt_couple
                                   radial_geom(alpha,beta)= &
                                        radial_geom(alpha,beta)+radial_factor
                                ENDDO
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO

           angular_geom=zeroc
           DO v=-1,1,2
              ! l+1 is an intermediate harmonic even at the radial cutoff.
              IF (l+v.LT.0) CYCLE
              DO t=-1,1
                 IF (mbra+t.LT.-(l+v) .OR. mbra+t.GT.l+v) CYCLE
                 DO tp=-1,1
                    DO lint=0,llv(ih,jatom)+lp
                       IF (NOTRI(lint,llv(ih,jatom),lp).LT.0) CYCLE
                       IF (NOTRI(lint,1,l+v).LT.0) CYCLE
                       DO nu=-lint,lint
                          IF (-nu+mmv(ih,jatom)+mket.NE.0) CYCLE
                          gaunt_pot=GAUNT(lint,llv(ih,jatom),lp,nu, &
                               mmv(ih,jatom),mket)
                          IF (ABS(gaunt_pot).LT.1.d-30) CYCLE
                          IF (-nu+tp+mbra+t.NE.0) CYCLE
                          gaunt_couple=GAUNT(lint,1,l+v,nu,tp,mbra+t)
                          IF (ABS(gaunt_couple).LT.1.d-30) CYCLE
                          DO beta=1,3
                             DO alpha=1,3
                                angular_factor=hf*( &
                                     conjg(c_am(alpha,tp))* &
                                     conjg(ca_st_lm(beta,v,t))+ &
                                     conjg(c_am(beta,tp))* &
                                     conjg(ca_st_lm(alpha,v,t)))
                                angular_geom(alpha,beta)= &
                                     angular_geom(alpha,beta)+ &
                                     angular_factor*gaunt_pot*gaunt_couple
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO

           IF (MAX(MAXVAL(ABS(radial_geom)), &
                MAXVAL(ABS(angular_geom))).LT.1.d-30) CYCLE
           iterm=iterm+1
           nsp_geometry(jatom)%ih(iterm)=ih
           nsp_geometry(jatom)%ibra(iterm)=index1
           nsp_geometry(jatom)%iket(iterm)=index2
           nsp_geometry(jatom)%radial(:,:,iterm)=radial_geom
           nsp_geometry(jatom)%angular(:,:,iterm)=angular_geom
        ENDDO
     ENDDO
  ENDDO

  nsp_geometry(jatom)%nterm=iterm
  nsp_geometry(jatom)%ready=.TRUE.
  IF (myrank.EQ.master) WRITE(6,'(A,3I8)') &
       ':DMFT_G2_GEOMETRY_CACHE ',jatom,iterm,maxterm
END SUBROUTINE build_nsp_geometry


subroutine init_audit_fact()
use fact
implicit none
integer :: i
fct=0d0
fct(1)=1d0
do i=3,99,2
 fct(i)=fct(i-2)*dble((i-1)/2)
enddo
end subroutine
subroutine angular_probe(lb,mb,lk,mk,lp,mp,rad,ang)
use w_atpar
implicit none
integer,intent(in) :: lb,mb,lk,mk,lp,mp
complex*16,intent(out) :: rad(3,3),ang(3,3)
integer :: it
if (allocated(nsp_geometry(1)%ih)) deallocate(nsp_geometry(1)%ih, &
 nsp_geometry(1)%ibra,nsp_geometry(1)%iket,nsp_geometry(1)%radial, &
 nsp_geometry(1)%angular)
nsp_geometry(1)%ready=.false.
lpv(1,1)=lb;lv(1,1)=lk;llv(1,1)=lp;mmv(1,1)=mp;Nhmx(1)=1
VPts=1d0;WPts=1d0
call build_nsp_geometry(1)
rad=0d0;ang=0d0
do it=1,nsp_geometry(1)%nterm
 if(nsp_geometry(1)%ibra(it)/=lb*(lb+1)+mb+1) cycle
 if(nsp_geometry(1)%iket(it)/=lk*(lk+1)+mk+1) cycle
 rad=nsp_geometry(1)%radial(:,:,it)
 ang=nsp_geometry(1)%angular(:,:,it)
enddo
end subroutine

subroutine quadrature_probe(lb,mb,lk,mk,lp,mp,rad,ang)
use g2_angular_test_module
implicit none
integer,intent(in) :: lb,mb,lk,mk,lp,mp
complex(dp),intent(out) :: rad(3,3),ang(3,3)
complex(dp),allocatable :: ag(:,:),ar(:,:,:),aa(:,:,:)
integer :: lmax,nh,i,j,q,a,b
lmax=max(lb,lk);nh=(lmax+1)**2
allocate(ag(nh,nh),ar(nh,nh,6),aa(nh,nh,6))
call g2_angular_tables(lmax,lp,mp,one,zero,ag,ar,aa)
i=lb*(lb+1)+mb+1;j=lk*(lk+1)+mk+1
rad=zero;ang=zero
do q=1,6
 a=axis_a(q);b=axis_b(q)
 rad(a,b)=ar(i,j,q);rad(b,a)=ar(i,j,q)
 ang(a,b)=aa(i,j,q);ang(b,a)=aa(i,j,q)
enddo
deallocate(ag,ar,aa)
end subroutine
