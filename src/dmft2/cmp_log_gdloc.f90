! @Copyright 2007 Kristjan Haule
! 
SUBROUTINE cmpLogGdloc(logG, eimp_nd, eimp_nd2, DeltaG, forb, TrGSigVdc, Gdloc, Edimp, s_oo, DM, Nds, Temperature, Sigind, Sigind_orig, cixdim, ncix, maxdim, maxsize, ntcix, sigma, nomega, csize, fh_sig, nipc, SOlapm)
  USE defs,  ONLY: pi, IMAG
  USE splines, ONLY: zspline3, zspline3der
  USE muzero,ONLY: nomq, jomq, iomq, womq, n0_om
  USE dmfts, ONLY: iso, natom, iatom, isort, cix, nl
  USE structure, ONLY: natm
  IMPLICIT NONE
  REAL*8, intent(out)    :: logG, eimp_nd, eimp_nd2, DeltaG, forb(3,natm), TrGSigVdc
  COMPLEX*16, intent(in) :: Gdloc(ntcix,nomega,nipc), s_oo(maxsize,ncix)
  REAL*8, intent(in)     :: Edimp(ntcix), Nds(nipc,ntcix)
  COMPLEX*16, intent(in) :: DM(maxdim,maxdim,ncix)
  INTEGER, intent(in)    :: Sigind(maxdim,maxdim,ncix), Sigind_orig(maxdim,maxdim,ncix), cixdim(ncix), csize(ncix)
  INTEGER, intent(in)    :: maxsize, ncix, maxdim, ntcix
  COMPLEX*16, intent(in) :: sigma(maxsize,ncix,nomega)
  INTEGER, intent(in)    :: nomega, nipc, fh_sig
  COMPLEX*16, intent(in) :: SOlapm(maxdim,maxdim,ncix)
  REAL*8 :: Temperature
  ! Functions
  Interface
     real*8 FUNCTION FreeE0(Energy, Temp)
       real*8, intent(in) :: Energy, Temp
     end FUNCTION FreeE0
     real*8 Function ReturnHighFrequency(A,C,eps,Temperature)
       real*8, intent(in) :: A, C, eps, Temperature
     end Function ReturnHighFrequency
     REAL*8 Function ferm(x)
       REAL*8, intent(in) :: x
     end Function ferm
     Function NumericSecondDeriv(func,omega,ii,nomega) result(df2)
       INTEGER, intent(in)    :: nomega, ii
       COMPLEX*16, intent(in) :: func(nomega)
       REAL*8, intent(in)    :: omega(nomega)
       COMPLEX*16 :: df2
     End Function NumericSecondDeriv
  end interface
  ! locals
  real*8, PARAMETER       :: Ry2eV = 13.60569253d0
  COMPLEX*16, allocatable  :: gc(:,:), geig(:), Delta(:,:), Delta_all(:), dDelta_all(:), CC(:), GDC(:), Gdloc_(:), sigma_(:), Gdloc__(:,:), sec_(:)
  REAL*8, allocatable      :: epsinf(:), ndtot(:), omega_all(:)
  INTEGER, allocatable     :: Sigini(:,:), Sigini_orig(:,:)
  INTEGER, allocatable     :: icx_ind(:), it_ind(:)
  COMPLEX*16 :: omn, sec, df2, df3
  REAL*8     :: tlogGD, tlogG0, tlogGX, A, C, A2, C2, correct, GDd, C0, second_order_correction, renormalize
  INTEGER    :: deg(ntcix), id
  INTEGER    :: icix, ip, iq, ndim, cixdm, it, it2, iom, cixdms, im, nom_all, iatm, icase, lcase, jatom, latom
  LOGICAL    :: DIAGONAL
  REAL*8     :: omega(nomega), Bi(ntcix)
  CHARACTER*3 :: sitx
  REAL*8      :: GS_dynamic(0:3,ntcix), GS_static(0:3,ntcix)
  COMPLEX*16  :: sigma0(maxsize,ncix,nomega)
  CHARACTER*1 :: label(0:3)=['c', 'x','y','z']
  LOGICAL,parameter :: Olap_Renormalize=.True.
  !LOGICAL,parameter :: Olap_Renormalize=.False.
  
  TrGSigVdc = 0.d0
  
  C0=1.
  DIAGONAL=.True.
  do icix=1,ncix
     ndim = cixdim(icix)
     do ip=1,ndim
        do iq=1,ndim
           it2 = Sigind_orig(ip,iq,icix)
           if (it2.gt.0 .and. ip.ne.iq) then
              DIAGONAL=.False.
           endif
        enddo
     enddo
  enddo

  do iom=1,nomega
     im=jomq(iom)
     omega(iom) = (2*im-1)*pi*Temperature
  enddo
  nom_all = jomq(nomega)
  allocate(omega_all(nom_all))
  do iom=1,nom_all
     omega_all(iom) = (2*iom-1)*pi*Temperature 
  enddo

  ! Here we reread self-energy and we do not use any extra broadening this time.
  ! We also subtract s_oo
  !if (nipc.eq.4) then
  do icix=1,ncix
     rewind(fh_sig+icix)
     CALL ReadSelfenergy(fh_sig+icix, sigma0(:,icix,:), omega, s_oo(:,icix), 0.d0, csize(icix), nomega, maxsize)
     do it=1,maxsize
        sigma0(it,icix,:) = sigma0(it,icix,:)-s_oo(it,icix)
     enddo
  enddo
  !endif
  
  forb(:,:)=0.d0
  logG=0
  DeltaG=0
  eimp_nd=0
  eimp_nd2=0
  if (DIAGONAL) then

     allocate( Gdloc_(nom_all) )
     Gdloc_=0

     ! Get impurity energy for approximation of G(infinity)
     allocate(epsinf(ntcix), ndtot(ntcix))
     allocate( icx_ind(ntcix), it_ind(ntcix) )
     epsinf=0
     deg=0
     ndtot=0
     do icix=1,ncix
        do ip=1,cixdim(icix)
           it = Sigind(ip,ip,icix)
           it2 = Sigind_orig(ip,ip,icix)
           if (it2.gt.0 .and. it.gt.0) then
              epsinf(it2) = Edimp(it2)+s_oo(it,icix)
              deg(it2) = deg(it2)+1
              ndtot(it2) = ndtot(it2)+DM(ip,ip,icix)
              Bi(it2) = -aimag(sigma(it,icix,nomega))*omega(nomega)
              icx_ind(it2)=icix
              it_ind(it2)=it
           endif
        enddo
     enddo

     allocate(Delta(nomega,ntcix))
     allocate(Delta_all(nom_all), dDelta_all(nom_all), CC(nom_all), GDC(nom_all))

     GS_dynamic(:,:) = 0.d0
     GS_static(:,:) = 0.d0
     
     do it2=1,ntcix
        if (deg(it2)==0) CYCLE
        icix = icx_ind(it2)
        it = it_ind(it2)

        ! Interpolating Gdloc on all Matsubara points ( stored in Gdloc_ )
        Gdloc_(:n0_om)=Gdloc(it2,:n0_om,1) ! The first few points should not be interpolated
        if (n0_om.lt.nomega) then
           df2 = NumericSecondDeriv(Gdloc(it2,:,1),omega,n0_om,nomega)
           df3 = -2./(IMAG*omega(nomega)-epsinf(it2))**3
           Gdloc_(n0_om:) = zspline3(omega(n0_om:),Gdloc(it2,n0_om:,1),omega_all(n0_om:),df2,df3) ! The tail needs to be interpolated
        endif
        if (.false.) then
           write(sitx,'(I3)') it2
           open(999,file=TRIM('gd_debug.')//trim(ADJUSTL(sitx)),status='unknown')
           do iom=1,nom_all
              write(999,'(F18.12,1x,2F18.12,2x,2F18.12,1x,2F18.12)') omega_all(iom)*Ry2eV, Gdloc_(iom)/Ry2eV, 1/(omega_all(iom)*IMAG-epsinf(it2))/Ry2eV, dble(log(-Gdloc_(iom))), dble(log(-1/(omega_all(iom)*IMAG-epsinf(it2))))
           enddo
           close(999)

           open(999,file=TRIM('gd2_debug.')//trim(ADJUSTL(sitx)),status='unknown')
           write(999,*) '#  nipc=', nipc
           do iom=1,nomega
              write(999,'(F18.12,3x)', advance='no') omega(iom)*Ry2eV
              do ip=1,nipc
                 write(999,'(2F18.12,3x)', advance='no') Gdloc(it2,iom,ip)/Ry2eV
              enddo
              write(999,'(2F18.12,3x)', advance='no') sigma(it,icix,iom)*Ry2eV
              write(999,*)
           enddo
           close(999)
        endif

        
        tlogGD=0
        tlogG0=0
        sec=0
        do iom=1,nom_all
           omn = IMAG*(2*iom-1)*pi*Temperature 
           tlogGD = tlogGD + dble(log(-Gdloc_(iom)))
           tlogG0 = tlogG0 + dble(log(-1/(omn-epsinf(it2))))
           sec = sec + 1./((omn-C0)*(omn-epsinf(it2)))
        enddo
        !print *, 'epsinf(it2)=', epsinf(it2)
        tlogGX = FreeE0(epsinf(it2),Temperature)
        second_order_correction = Bi(it2)* ( (ferm(C0/Temperature)-ferm(epsinf(it2)/Temperature))/(C0-epsinf(it2)) - 2*Temperature*dble(sec) )
        logG = logG + (2*Temperature*(tlogGD-tlogG0) + tlogGX + second_order_correction )*deg(it2)
        eimp_nd = eimp_nd + Edimp(it2)*ndtot(it2)     ! Tr(eimp*n)
        eimp_nd2 = eimp_nd2 + epsinf(it2)*ndtot(it2)  ! Tr((eimp+s_oo)*n)
        !print *, 'sginf=', second_order_correction*deg(it2)*Ry2eV
        if (.false.) then
           WRITE(6,'(I3,1x,A,F13.8,1x,A,F13.8,1x,A,F13.8,1x,A,F13.8,1x,A,F13.8,1x,A,F13.8,1x,A,F16.8,1x,A,F16.8,1x,A,F13.8)') it2,'logG=', logG, 'ndtot=', ndtot(it2),'F0=', tlogGX, 'Eimp+soo=', epsinf(it2), 'Edimp=', Edimp(it2), 'sndocor=', second_order_correction, 'tlogGD=', tlogGD, 'tlogG0=', tlogG0, 'sec=', dble(sec)
        endif
        
        do iom=1,nomega
           im=jomq(iom)
           omn = IMAG*(2*im-1)*pi*Temperature 
           Delta(iom,it2) = omn-1/Gdloc(it2,iom,1)-Edimp(it2)-sigma(it,icix,iom)
        enddo
        call zspline3der(omega, Delta(:,it2), omega_all, Delta_all, dDelta_all)

        CC = Delta_all(:)-omega_all(:)*dDelta_all(:)
        CALL GetHighFrequency(A,C,CC,omega_all,nom_all)
        !print *, 'A=', A, 'C=', C, 'E=', epsinf(it2)
        
        GDC = Gdloc_(:)*CC - A/((omega_all(:)*IMAG-epsinf(it2))*(omega_all(:)*IMAG-C))
        correct = ReturnHighFrequency(A,C,epsinf(it2),Temperature)
        !print *, 'correct=', correct
        GDd = 2*Temperature*dble(sum(GDC)) + correct
        DeltaG = DeltaG + GDd*deg(it2) 

        ! NEW FOR FORCES
        ! Interpolating sigma on all Matsubara points
        allocate( sigma_(nom_all) )
        sigma_(:n0_om)= sigma0(it,icix,:n0_om) ! The first few points should not be interpolated
        if (n0_om.lt.nomega) then
           df2 = NumericSecondDeriv(sigma0(it,icix,:),omega,n0_om+1,nomega)
           df3 = NumericSecondDeriv(sigma0(it,icix,:),omega,nomega-1,nomega)
           sigma_(n0_om:) = zspline3(omega(n0_om:),sigma0(it,icix,n0_om:), omega_all(n0_om:), df2, df3) ! The tail needs to be interpolated
        endif
        CALL GetHighFrequency(A,C,sigma_,omega_all,nom_all)
        C=1.d0
        if (.false.) then
           write(sitx,'(I3)') it2
           open(999,file=TRIM('sig_debug.')//trim(ADJUSTL(sitx)),status='unknown')
           write(999,'(A,F12.8,1x,A,F12.8)') '# A=', A, 'C=', C
           do iom=1,nom_all
              write(999,'(F18.12,1x,2F18.12,2x)') omega_all(iom), sigma_(iom)
           enddo
           close(999)
        endif
        
        do ip=1,nipc
           if (ip.eq.1) then
              A2 = 1.d0
              C2 = epsinf(it2)
           else
              ! Interpolating vector form of Gdloc on all Matsubara points
              Gdloc_(:n0_om) = Gdloc(it2,:n0_om,ip) ! The first few points should not be interpolated
              if (n0_om.lt.nomega) then
                 df2 = NumericSecondDeriv(Gdloc(it2,:,ip),omega,n0_om+1,nomega)
                 df3 = NumericSecondDeriv(Gdloc(it2,:,ip),omega,nomega-1,nomega)
                 Gdloc_(n0_om:) = zspline3(omega(n0_om:),Gdloc(it2,n0_om:,ip),omega_all(n0_om:),df2,df3) ! The tail needs to be interpolated
              endif
              CALL GetHighFrequency(A2,C2,Gdloc_,omega_all,nom_all)
           endif
           
           GDC(:) = Gdloc_(:) * sigma_(:) - A*A2/((omega_all(:)*IMAG-C2)*(omega_all(:)*IMAG-C))
           correct = ReturnHighFrequency(A*A2, C, C2, Temperature)
           
           GS_dynamic(ip-1,it2) = GS_dynamic(ip-1,it2) + ( 2*Temperature*dble(sum(GDC)) + correct )*deg(it2)*2./iso
           GS_static (ip-1,it2) = GS_static(ip-1,it2)  + Nds(ip,it2)*dble(s_oo(it,icix))
           
           if (.false.) then
              write(sitx,'(I1,A,I1)') it2,'.',ip
              open(999,file=TRIM('gd3_debug.')//trim(ADJUSTL(sitx)),status='unknown')
              write(999,*) '#  nipc=', nipc, 'deg=', deg(it2), 'A=', A2, 'C=', C2, 'correct=', correct, 'dynamic=', GS_dynamic(ip-1,it2), 'static=', GS_static (ip-1,it2)
              do iom=1,nom_all
                 write(999,'(F18.12,1x,2F18.12,2x,F18.12)') omega_all(iom)*Ry2eV, Gdloc_(iom)/Ry2eV, dble(GDC(iom))
              enddo
              close(999)
           endif
        enddo
        deallocate( sigma_)
        if (.false.) then
           write(sitx,'(I3)') it2
           open(888+it2,file=TRIM('Deltax.')//trim(ADJUSTL(sitx)),status='unknown')
           do iom=1,nom_all
              write(888+it2,'(f10.6,1x)',advance='no') (2*iom-1)*pi*Temperature*Ry2eV
              write(888+it2,'(f10.6,1x,f10.6,2x,f10.6,1x,f10.6,2x)',advance='no') dble(Delta_all(iom))*Ry2eV, aimag(Delta_all(iom))*Ry2eV, dble(dDelta_all(iom))*Ry2eV, aimag(dDelta_all(iom))*Ry2eV
              write(888+it2,'(f10.6,1x,f10.6,2x)',advance='no') dble(GDC(iom)), aimag(GDC(iom))
              write(888+it2,*)
           enddo
           close(888+it2)
        endif
     enddo
     ! GS_dynamic(1:3,it2)
     ! GS_static (1:3,it2)
     TrGSigVdc = sum(GS_static(0,:)+GS_dynamic(0,:))

     WRITE(6,*) 'Nds for each orbital:'
     WRITE(6,'(A)',advance='no') '   c:'
     do it2=1,ntcix
        WRITE(6,'(f14.7,1x)',advance='no') Nds(1,it2)
     enddo
     WRITE(6,*)
     if (nipc.eq.4) then
        WRITE(6,'(A)',advance='no') '   x:'
        do it2=1,ntcix
           WRITE(6,'(f14.7,1x)',advance='no') Nds(2,it2)
        enddo
        WRITE(6,*)
        WRITE(6,'(A)',advance='no') '   y:'
        do it2=1,ntcix
           WRITE(6,'(f14.7,1x)',advance='no') Nds(3,it2)
        enddo
        WRITE(6,*)
        WRITE(6,'(A)',advance='no') '   z:'
        do it2=1,ntcix
           WRITE(6,'(f14.7,1x)',advance='no') Nds(4,it2)
        enddo
        WRITE(6,*)
     endif
     WRITE(6,*) '-Tr(G*(s_oo-Vdc)) for each orbital in mRy/Br:'
     do ip=1,nipc
        WRITE(6,'(A,A,A)',advance='no') '   ',label(ip-1),':'
        do it2=1,ntcix
           WRITE(6,'(f14.7,1x)',advance='no') -GS_static(ip-1,it2)*1000
        enddo
        WRITE(6,*)
     enddo
     WRITE(6,'(A)') '-Tr(G*(Sigma-s_oo)) for each orbital in mRy/Br:'
     do ip=1,nipc
        WRITE(6,'(A,A,A)',advance='no') '   ',label(ip-1),':'
        do it2=1,ntcix
           WRITE(6,'(f14.7,1x)',advance='no') -GS_dynamic(ip-1,it2)*1000
        enddo
        WRITE(6,*)
     enddo
     WRITE(6,'(A)') '-Tr(G*(Sigma-Vdc)) for each orbital in mRy/Br:'
     do ip=1,nipc
        WRITE(6,'(A,A,A)',advance='no') '   ',label(ip-1),':'
        do it2=1,ntcix
           WRITE(6,'(f14.7,1x)',advance='no') -(GS_static(ip-1,it2)+GS_dynamic(ip-1,it2))*1000
        enddo
        WRITE(6,*)
     enddo
     if (Olap_Renormalize) then
     WRITE(6,'(A)') '-Tr(G*(Sigma-Vdc))*Olap for each orbital in mRy/Br:'
     do ip=1,nipc
        WRITE(6,'(A,A,A)',advance='no') '   ',label(ip-1),':'
        do it2=1,ntcix
           icix = icx_ind(it2)
           it = it_ind(it2)
           renormalize = 1/(SOlapm(it,it,icix)*SOlapm(it,it,icix))
           WRITE(6,'(f14.7,1x)',advance='no') -(GS_static(ip-1,it2)+GS_dynamic(ip-1,it2))*renormalize*1000
        enddo
        WRITE(6,*)
     enddo
     endif
     
     if (nipc.eq.4) then
        DO icase=1,natom ! correlated atoms only
           latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
           jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )
           do lcase=1,nl(icase) 
              icix = cix(icase,lcase)
              if ( icix.EQ.0 ) CYCLE
              do ip=1,cixdim(icix)
                 do iq=1,cixdim(icix)
                    renormalize = 1.d0
                    if (Olap_Renormalize) renormalize = 1/(SOlapm(ip,ip,icix)*SOlapm(iq,iq,icix))
                    it2 = Sigind_orig(ip,iq,icix)
                    if (it2.gt.0) then
                       forb(1:3,latom) = forb(1:3,latom) - (GS_static(1:3,it2)+GS_dynamic(1:3,it2))*renormalize/deg(it2)
                    endif
                 enddo
              enddo
           enddo
        ENDDO
        
        WRITE(6,*) 'forb[mRy/Br]='
        do ip=1,3
           WRITE(6,'(A,A,A)',advance='no') '    ',label(ip),':'
           DO icase=1,natom
              latom = iatom(icase)
              WRITE(6,'(f16.10,1x)',advance='no') forb(ip,latom)*1000
           enddo
           WRITE(6,*)
        enddo
     endif
     
     deallocate(Delta_all, dDelta_all, CC, GDC)
     
     if (.False.) then
        open(888,file='Delta.xx',status='unknown')
        do iom=1,nomega
           im=jomq(iom)
           write(888,'(f10.6,1x)',advance='no') (2*im-1)*pi*Temperature*Ry2eV
           do it2=1,ntcix
              write(888,'(f10.6,1x,f10.6,2x)',advance='no') dble(Delta(iom,it2))*Ry2eV, aimag(Delta(iom,it2))*Ry2eV
           enddo
           write(888,*)
        enddo
        close(888)
     endif
     deallocate(Delta)
     !print *, 'eimp_nd=', eimp_nd
     deallocate( icx_ind, it_ind )
     deallocate( epsinf, ndtot )
     deallocate( Gdloc_ )
  else
  
     allocate( Gdloc__(ntcix,nom_all) )
     
     do icix=1,ncix
        cixdm = cixdim(icix)
         
        allocate( Sigini(cixdm,cixdm), Sigini_orig(cixdm,cixdm) )
        CALL Get_Sigini(cixdms, Sigini, Sigini_orig, Sigind(:,:,icix), Sigind_orig(:,:,icix), cixdm, maxdim)
        
        allocate( gc(cixdms,cixdms), geig(cixdms), epsinf(cixdms) )
        allocate( sec_(cixdms) )
        epsinf=0
        do ip=1,cixdms
           it = Sigini(ip,ip)
           it2 = Sigini_orig(ip,ip)
           if (it2.gt.0 .and. it.gt.0) then
              epsinf(ip) = Edimp(it2)+s_oo(it,icix)
              Bi(ip) = -aimag(sigma(it,icix,nomega))*omega(nomega)
           endif
        enddo
        
        !print *, 'cmp_log_gdloc'
        !do ip=1,cixdms
        !   do iq=1,cixdms
        !      it2 = Sigini_orig(ip,iq)
        !      if (it2.gt.0) then
        !         print *, ip, iq, it2, Gdloc(it2, 1, 1)
        !      endif
        !   enddo
        !enddo
        
        ! Interpolating Gdloc on all Matsubara points ( stored in Gdloc_ )
        Gdloc__(:,:n0_om)=Gdloc(:,:n0_om,1) ! The first few points should not be interpolated
        
        do ip=1,cixdms
           do iq=1,cixdms
              it2 = Sigini_orig(ip,iq)
              id=0
              if (ip.eq.iq) id=1
              if (it2.gt.0 .and. n0_om.lt.nomega) then
                 df2 = NumericSecondDeriv(Gdloc(it2,:,1),omega,n0_om,nomega)
                 df3 = -2.d0*id/(IMAG*omega(nomega)-epsinf(it2))**3
                 Gdloc__(it2,n0_om:) = zspline3(omega(n0_om:),Gdloc(it2,n0_om:,1),omega_all(n0_om:),df2,df3) ! The tail needs to be interpolated
              endif
           enddo
        enddo

        tlogGD=0.d0
        tlogG0=0.d0
        sec_(:)=0.d0
        do iom=1,nom_all
           omn = IMAG*(2*iom-1)*pi*Temperature 
           gc(:,:)=0.d0
           do ip=1,cixdms
              do iq=1,cixdms
                 it2 = Sigini_orig(ip,iq)
                 if (it2.gt.0) then
                    gc(ip,iq) = Gdloc__(it2,iom)
                 endif
              enddo
           enddo
           CALL eigvals0(gc, geig, cixdms) ! eigenvalues and eigenvectors
           
           do ip=1,cixdms
              tlogGD = tlogGD + dble(log(-geig(ip)))
              tlogG0 = tlogG0 + dble(log(-1/(omn-epsinf(ip))))
              sec_(ip) = sec_(ip) + 1./((omn-C0)*(omn-epsinf(ip)))
           enddo
        enddo

        tlogGX=0.d0
        second_order_correction=0.d0
        do ip=1,cixdms
           tlogGX = tlogGX + FreeE0(epsinf(ip),Temperature)
           second_order_correction = Bi(ip)* ( (ferm(C0/Temperature)-ferm(epsinf(ip)/Temperature))/(C0-epsinf(ip)) - 2*Temperature*dble(sec_(ip)) )
        enddo
        logG = logG + 2*Temperature*(tlogGD-tlogG0) + tlogGX + second_order_correction
        !print *, 'tlogGD=', tlogGD, 'tlogG0=', tlogG0, 'tlogGX=', tlogGX, 'second_order=', second_order_correction
        deallocate( sec_ )
        deallocate( Sigini, Sigini_orig )
        deallocate( gc, geig, epsinf )
     enddo
     deallocate( Gdloc__ )
  endif
  deallocate(omega_all)
  !print *, 'logG_final[Ry]=', logG
  !print *, 'logG_final[eV]=', logG*Ry2eV
END SUBROUTINE cmpLogGdloc

SUBROUTINE Get_Sigini(cixdms, Sigini, Sigini_orig, Sigind, Sigind_orig, cixdm, maxdim)
  IMPLICIT NONE
  INTEGER, intent(in)  :: Sigind(maxdim,maxdim), Sigind_orig(maxdim,maxdim), cixdm, maxdim
  INTEGER, intent(out) :: cixdms, Sigini(cixdm,cixdm), Sigini_orig(cixdm,cixdm)
  ! locals
  INTEGER :: cind(cixdm), it, ip, iq, cini(cixdm)
  
  cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated
  cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block
  DO ip=1,cixdm
     it = Sigind_orig(ip,ip)
     if (it.gt.0) then
        cixdms = cixdms + 1
        cind(ip) = cixdms
     endif
  ENDDO
     
  do ip=1,cixdm
     if (cind(ip).gt.0) cini(cind(ip))=ip
  enddo
  Sigini=0
  Sigini_orig=0
  DO ip=1,cixdms
     do iq=1,cixdms
        Sigini(ip,iq) = abs(Sigind(cini(ip),cini(iq)))           
        Sigini_orig(ip,iq) = abs(Sigind_orig(cini(ip),cini(iq))) 
     enddo
  ENDDO
END SUBROUTINE Get_Sigini

SUBROUTINE GetHighFrequency(A,C,CC,omega,N)
  ! Approximates CC ~  A/(i*om-C)
  IMPLICIT NONE
  REAL*8, intent(out)    :: A, C
  COMPLEX*16, intent(in) :: CC(N)
  REAL*8, intent(in)     :: omega(N)
  INTEGER, intent(in)    :: N
  !print *, 'CC(N)=', CC(N)
  if (abs(CC(N)).lt.1e-6) then
     A=0.d0
     C=1.d0
  else
     A = 1./aimag(1/(CC(N)*omega(N)))
     C = -A*dble(1./CC(N))
  endif
END SUBROUTINE GetHighFrequency

REAL*8 Function ReturnHighFrequency(A,C,eps,Temperature)
  IMPLICIT NONE
  REAL*8, intent(in) :: A, C, eps, Temperature
  Interface
     REAL*8 Function ferm(x)
       REAL*8, intent(in) :: x
     end Function ferm
  end Interface
  if (abs(A).gt.1e-10) then
     ReturnHighFrequency = A*(ferm(eps/Temperature)-ferm(C/Temperature))/(eps-C)
  else
     ReturnHighFrequency = 0.d0
  endif
  return
END Function ReturnHighFrequency

Function NumericSecondDeriv(funcy,omega,ii,nomega) result(df2)
  ! Numerically determines the derivative of function func(omega) at point ii
  IMPLICIT NONE
  COMPLEX*16 :: df2
  INTEGER, intent(in)   :: nomega, ii
  COMPLEX*16, intent(in):: funcy(nomega)
  REAL*8, intent(in)    :: omega(nomega)
  ! locals
  COMPLEX*16 :: f0, f1, f2
  REAL*8     :: dx0, dx1, dx2
  if (ii.le.1 .or. ii.ge.nomega) print *, 'Can not determine derivative outside the range of a function'
  f0 = funcy(ii-1)
  f1 = funcy(ii)
  f2 = funcy(ii+1)
  dx0 = omega(ii)-omega(ii-1)
  dx1 = omega(ii+1)-omega(ii)
  df2 = ((f2-f1)/dx1 - (f1-f0)/dx0)/(0.5*(dx0+dx1))
End Function NumericSecondDeriv
