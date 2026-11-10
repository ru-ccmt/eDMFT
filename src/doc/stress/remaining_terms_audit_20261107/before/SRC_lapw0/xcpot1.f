      SUBROUTINE XCPOT1 (luse2,LM,JATOM,LMMAX,LR2V,JSPIN,sigmaj,&
                         imfield,bextene,iprint,switch1,sigma_str)
      use struct, only: nat, iatnr, jri, mult
      use densit, only: lmmaxx,clmsp
      use work, only: exc,vxc,vxcold,vxcnonloc,r,vtau,vlapl,v,slater_tau_lm
      use vresp1, only: vclmsp
      use parallel, only :nstart,nstop,myid,abort_parallel
      USE xcparam, only: GGA_SWITCH, MGGA_TAU_SWITCH, MGGA_LAPL_SWITCH, VRESP_SWITCH, VX_SWITCH, &
      EX_SWITCH, READ_TAU_SWITCH, VDISP_SWITCH, VHALF_SWITCH, ID_FUNC_X, xc1, xc2, xc3, xc4, xcconst
      use TestCases,only : Lebedev, LEECE
!      use gga_str_kb, only: sig_str_u,sig_str_d, sig_str_ud  !,sigma_str    !for the stress tensor
      use str_switch, only : SW_STRESS, STRESS
!                                                                       
!     XCPOT1 CALCULATES THE EXCHANGE-CORRELATION POTENTIAL INSIDE       
!     THE SPHERES.                                                      
!     AFTER SUMATION OF  CHARGE DENSITY  BY CLM(R) FOR   THETA       
!     VALUES PER MESH-POINT, THE ( CUBIC ) HARMONIC-EXPANSION  OF       
!     EXCHANGE-POTENTIAL IS GENERATED AND STORED IN ARRAY VEX.
!                                                                       
!                                                                       
!     IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
!
      INCLUDE 'param.inc'

!     arguments
      integer, intent(in) :: luse2, iprint, jatom, imfield, jspin
      integer, intent(in) :: LM(2,NCOM+3,NAT),LMMAX(NAT)
      real*8, intent(in)  :: bextene
      logical, intent(in) :: lr2v
      character*5, intent(in):: switch1
      real*8, intent(out) :: sigmaj(2)
      real*8, intent(out) :: sigma_str(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:3)

!     local
      LOGICAL        :: robfit, spin, nospin
      real*8         :: YY(NCOM,nspd),YYB(nspd,NCOM),ROTA1(3)
      REAL*8            bbu(NCOM),bbbu(NCOM),bbd(NCOM), bbbd(NCOM), &
           vtau_bu(ncom),vtau_bd(ncom),vlapl_bu(ncom),vlapl_bd(ncom)
      COMPLEX*16, allocatable :: YL(:), dyt(:),dyt2(:), dyp(:),dyp2(:), dytp(:)
      !COMPLEX*16     :: YL((lmax2X+1)*(lmax2X+1)),&
      !     dyt((lmax2X+1)*(lmax2X+1)),dyt2((lmax2X+1)*(lmax2X+1)), &
      !     dyp((lmax2X+1)*(lmax2X+1)),dyp2((lmax2X+1)*(lmax2X+1)), &
      !     dytp((lmax2X+1)*(lmax2X+1))
      !real*8 :: spt(3,nspd),weight(nspd), &
      !     vxcsu(nspd),vxcsd(nspd), &
      !     vxu(nspd),vxd(nspd),exu(nspd),exd(nspd), &
      !     vtau_rad_up(nspd),vtau_rad_dn(nspd), &
      !     vlapl_rad_up(nspd),vlapl_rad_dn(nspd),&
      !     dcsp(ncom,nrad,2),ddcsp(ncom,nrad,2), &
      !     dydt(ncom,nspd),dydt2(ncom,nspd), &
      !     dydp(ncom,nspd),dydp2(ncom,nspd),dydtp(ncom,nspd), &
      !     qrhou(ncom),qrhod(ncom), &
      !     qsin(nspd),qcos(nspd)
      real*8, allocatable :: spt(:,:),weight(:),vxcsu(:),vxcsd(:),vxu(:),vxd(:),exu(:),exd(:)
      real*8, allocatable :: vtau_rad_up(:),vtau_rad_dn(:),vlapl_rad_up(:),vlapl_rad_dn(:)
      real*8, allocatable :: dcsp(:,:,:),ddcsp(:,:,:),dydt(:,:),dydt2(:,:),dydp(:,:),dydp2(:,:)
      real*8, allocatable :: dydtp(:,:),qrhou(:),qrhod(:),qsin(:),qcos(:)
!     derivatives of kinetic energy density
      real*8   ::  dtau_dr(ncom,nrad,2),d2tau_dr2(ncom,nrad,2), &
           dtau_dt_u, dtau_dt_d, dtau_dt_t, d2tau_dt2_u , d2tau_dt2_d, &
           d2tau_dt2_t, dtau_dp_u, dtau_dp_d, dtau_dp_t, d2tau_dp2_u, &
           d2tau_dp2_d, d2tau_dp2_t, d2tau_dtdp_u, d2tau_dtdp_d, &
           d2tau_dtdp_t, dtau_dr_u, dtau_dr_d, dtau_dr_t, d2tau_drdt_u, &
           d2tau_drdt_d, d2tau_drdt_t, d2tau_drdp_u, d2tau_drdp_d, &
           d2tau_drdp_t, d2tau_d2r_t, dtau_dx, dtau_dy, dtau_dz, &
           dtau_mag, d2tau, d2tau_x, d2tau_y, d2tau_z, dtau_d2tau, d2tau_d2r_u, &
           dtau_dx_u, dtau_dy_u, dtau_dz_u, dtau_mag_u, d2tau_u, &
           d2tau_x_u, d2tau_y_u, d2tau_z_u, dtau_d2tau_u, d2tau_d2r_d, &
           dtau_dx_d, dtau_dy_d, dtau_dz_d, dtau_mag_d, d2tau_d, d2tau_x_d, &
           d2tau_y_d, d2tau_z_d, dtau_d2tau_d, dtau_drho, dtau_u_drho_u, &
           dtau_u_drho_d, dtau_d_drho_u, dtau_d_drho_d, d2tau_dr2_u, &
           d2tau_dr2_d, d2tau_dr2_t, sig_str_u,sig_str_d,sig_str_ud

      integer  ::  quadtest,iatj,l,ll,iatest,ncomx,j,index,llmm,luse,&
           npt1,lm1,ind,init,ir,k,irsigmaup,irsigmadn
      real*8   :: dmaxup,dmaxdn,sigmamaxup,sigmamaxdn, &
           sum0,rp,rpp,ft,fu,fd,tautot,tauup,taudn,vxuold,vxdold,vxtold, &
           vxunonloc,vxdnonloc,gx,gxd,gxu,gy,gyd,gyu,gmag,gmagd,gmagu, &
           g2,g2d,g2u,ggx,ggxd,ggxu,ggy,ggyd,ggyu,ggz,ggzd,ggzu, &
           gdggt,gdggu,gdggd,vxcd,vxcu,vhup,vhdn,vhtot,gz,gzu,gzd, &
           d2ndru,d2ndrd,d2ndrt,dndtu,dndtd,dndtt,d2ndtu,d2ndtd,d2ndtt, &
           dndpu,dndpd,dndpt,d2ndpu,d2ndpd,d2ndpt,dndtpu,dndtpd,dndtpt, &
           dndrtu,dndrtd,dndrtt,dndrpu,dndrpd,dndrpt,dndru,dndrd,dndrt, &
           excd,excu,vxcua,vxcda,excua,excda,test,sigma,y,t, &
           vxtnonloc, vtauc_up, vtauc_dn, vlaplc_up, vlaplc_dn, &
           vtau_ua,vtau_da,vlapl_ua,vlapl_da
!          ,tmp,c1,c2,c3,c4, tmp1,tmp2,tmp3,tmp4,
      real*8   :: pi,rtpi,wnorm,sqfpi
      
!      parameters
       real*8    :: cutmuffin, thrd
       logical   :: debug
       parameter (thrd=1.D0/3.D0)
       parameter (debug=.false. )
       parameter (CutMuffin=1D-8)  !-12)

!      Function for Kahan summation
       real*8 ddotk

!                                                                       

!      ***** EECE Patching *****
       real*8 U(1),VT(1),SWORK(12*nspd), Sing(nspd)
       real*8 du_h(nspd), dd_h(nspd), tmpup(nspd), tmpdn(nspd)
       logical UseLM(ncom), UseInterp, UseLDA, there, EECE_SWITCH
       character*6 EECE_Show
       real*4,allocatable :: rwork(:)
       real*8,allocatable :: tmp(:,:)
       integer,allocatable:: itmp(:)
       integer nn,izer,i,NCOND,LWORK,INFO
       real*8 RUp,RDn,RCond
       
!******* for the stress tensor
       REAL*8 :: gxyzu1(nspd),gxyzd1(nspd),gxyzud1(nspd)
!
!       STRESS = .false.
!       IF(SW_STRESS .eq. 'YES') STRESS = .true.

       allocate(YL((lmax2X+1)*(lmax2X+1)),dyt((lmax2X+1)*(lmax2X+1)),dyt2((lmax2X+1)*(lmax2X+1)))
       allocate(dyp((lmax2X+1)*(lmax2X+1)),dyp2((lmax2X+1)*(lmax2X+1)), dytp((lmax2X+1)*(lmax2X+1)))
       allocate(spt(3,nspd),weight(nspd),vxcsu(nspd),vxcsd(nspd),vxu(nspd),vxd(nspd),exu(nspd),exd(nspd))
       allocate(vtau_rad_up(nspd),vtau_rad_dn(nspd),vlapl_rad_up(nspd),vlapl_rad_dn(nspd),dcsp(ncom,nrad,2))
       allocate(ddcsp(ncom,nrad,2),dydt(ncom,nspd),dydt2(ncom,nspd),dydp(ncom,nspd),dydp2(ncom,nspd))
       allocate(dydtp(ncom,nspd),qrhou(ncom),qrhod(ncom),qsin(nspd),qcos(nspd))
       
       UseInterp= .false.
       UseLDA   = .false.
       EECE_SWITCH=.false.
       if(switch1.eq.' EECE') EECE_SWITCH=.true.
       if(switch1.eq.' HYBR') EECE_SWITCH=.true.

       if(EECE_SWITCH)then

!       Experimental Options
        inquire(file='EECE_LDA',exist=there)
        if(there)UseLDA=.true.
        inquire(file='EECE_INT',exist=there)
        if(there)UseInterp=.true.

       endif

       UseLM=.true.

       if(EECE_SWITCH)then
        iatj = IATNR( JATOM )   
!       Set odd values to zero
        do L=1,LMMAX( JATOM )
                LL=IABS(LM(1,L,JATOM)) 
                if( 2*(LL/2) .ne. LL) UseLM(L)=.false.
!       Cubic -3,2
                if(iatj.GT.0)then
                        if(lm(1,L,jatom).EQ.-3.AND.lm(2,L,jatom).EQ.2)UseLM(L)=.true.
                endif
!               write(6,*)L,LL,LM(1:2,L,JATOM),UseLM(L)
        enddo
       endif
!
!       **** EECE End ****

       dmaxup=0
       dmaxdn=0
       sigmamaxup=0.d0
       sigmamaxdn=0.d0
       irsigmaup = jri(jatom)
       irsigmadn = jri(jatom)

!      Spin simplifications
        if(jspin .eq. 2)then
                spin=.true.
        else
                spin=.false.
        endif
!       Debug
        if(debug)spin=.true.
        nospin=.not. spin
!---------------------------------------------------------------------  
!.....DEFINE CONSTANTS, INITIALIZE VARIABLES 
      IF(.NOT.ALLOCATED(vclmsp)) allocate(vclmsp(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2))
         pi=acos(-1.d0)
         rtpi=sqrt(4.D0*pi)
         wnorm=4.D0*pi
       1980 FORMAT(3X)                                                        
        iatest=-1

 2000 FORMAT(16X,I2//)                                                  
 2010 FORMAT(16X,I2,5X,I2/)                                             
 2020 FORMAT(3X,5E13.7)                                                 
 2021 FORMAT(3X,4E19.12)                                                 
 2030 FORMAT(///)                                                       
 2031 FORMAT(/)                                                         
!                                                                       
!                                                                       
!fit      aa2=0.d0
      NCOMX=NCOM                                                        
      SQFPI=1.D0/rtpi                                                
!
      INDEX=1
      DO J=1,JATOM-1
        INDEX=INDEX+MULT(J)
      ENDDO
      LLMM=LMMAX( JATOM )
      iatj = IATNR( JATOM )                                                 
!************************************************
!
!     generates Gauss-Legendre points on sphere
!
      quadtest = 0
      LUSE=abs(LUSE2)
      if(LUSE .gt. 100)then
        quadtest=LUSE/100
        LUSE = LUSE - quadtest*100
      endif
      do L=1,LLMM
        LUSE=MAX(LUSE,abs(lm(1,L,JATOM)))
      enddo               
      luse=min(luse,lmax2X) 

!     Conditioning tests for DyScO3 indicates this is far better
      if(EECE_SWITCH .and. LEECE )Lebedev   = .true.
!        Lebedev   = .false.  !added and comment Kamal 2021 August.

      call gpoint(spt,weight,npt1,LUSE)
      if(myid .eq. 0 .and. Lebedev .and. jatom.eq.1 )then
          write(21,*)'Lebedev grid of ',NPT1
          write(6,*)'Lebedev grid of ',NPT1
      endif
!
!-----------
!     Since EECE is symmetric, only need half the points
      if(EECE_SWITCH  .and.  Lebedev) then
!       Some care here
!
!       Phase 1, sort to then just use z >= 0
        allocate (rwork(npt1), itmp(npt1) )        
        do J=1,npt1
                rwork(j) = spt(3,j)
                itmp(J)    = j
        enddo
        call sortag(rwork,npt1,itmp)

        deallocate(rwork) ; allocate(tmp(4,npt1) )
        tmp(1:3,1:npt1) = spt(1:3,1:npt1) 
        tmp(4,1:npt1)   = weight(1:npt1)
        nn=0
        izer=-9
        do j=npt1,1,-1
                i=itmp(j)
                t=tmp(3,i)
                if(t .gt. 1D-12)then
!                       Store, double weight
                        nn=nn+1
                        spt(1:3,nn) = tmp(1:3,i)
                        weight(nn)  = tmp(4,i)*2.D0
!                       write(6,'(a,I4,5ES12.4)')'Stored ',nn,spt(1:3,nn),weight(nn)
                else if(t .gt. -1D-12)then
!                       Store, weight for the moment (should reduce z=0 case)
                        if(izer .lt. 0)izer = nn
                        nn=nn+1
                        spt(1:3,nn) = tmp(1:3,i)
                        weight(nn)  = tmp(4,i)*2.D0
!                       write(6,'(a,I4,5ES12.4)')'Store0 ',nn,spt(1:3,nn),weight(nn)
                endif
        enddo
        npt1=nn
!
!       Now remove duplicates from izer to npt1 with z=0
        tmp(1:3,izer:npt1) = spt(1:3,izer:npt1)
        tmp(4  ,izer:npt1) = weight (izer:npt1)
        tmp(4  ,npt1)=-100.D0

        do J=izer,npt1-1
                if(tmp(4,j) .lt. -1.D0)cycle
                do k=j+1,npt1
                        if(     (abs(spt(1,j)+spt(1,k)) .lt. 1D-8 ) .and. &
                                (abs(spt(2,j)+spt(2,k)) .lt. 1D-8 )) &
                                tmp(4,k) = -100.D0
                enddo
        enddo
        nn=izer-1
        do j=izer,npt1
                if(tmp(4,j) .gt. -0.5D0)then
                        nn=nn+1
                        spt(:,nn)=tmp(1:3,j)
                        weight(nn)=tmp(4,j)
!                       write(6,'(a,I4,5ES12.4)')'Keep ',j,spt(1:3,j),weight(j)
                endif
        enddo                                 
!       if(myid .eq. 0)write(21,2121)npt1,nn
!2121   format('Lebedev Quadrature of ',i4,' Reduced to ',i4)
        npt1=nn

        deallocate (tmp, itmp)

      else if( EECE_SWITCH)then
        npt1=npt1/2
        weight(1:npt1)=weight(1:npt1)*2.0D0
      endif
!-----------
!     Note, we don't need LUSE later to be big
      
      LUSE=0
      do LM1=1,LLMM
        luse = max( iabs(LM(1,LM1,JATOM)), luse)
      enddo
      LUse=LUse+2
      Sum0 = 0.D0
      do 1 ind=1,npt1
         rota1(1)=spt(1,ind)
         rota1(2)=spt(2,ind)
         rota1(3)=spt(3,ind)
         CALL YLM (ROTA1,luse,YL)                                          
!
         if(GGA_SWITCH.eqv..true.) then
            call dylm(rota1,luse,yl,dyt,dyt2,dyp,dyp2,dytp,qsin(ind),qcos(ind))
         endif
!
         DO LM1=1,LMMAX(JATOM)                                             
            CALL SUML(LM1,YY(LM1,IND),YL,LM(1,1,JATOM))                      
         end DO
!
         sum0 = sum0 + YY(1,ind)*weight(ind)

         if(GGA_SWITCH.eqv..true.) then
            do lm1=1,LMMAX(JATOM)
               call suml(lm1,dydt(lm1,ind),dyt,lm(1,1,jatom))
               call suml(lm1,dydt2(lm1,ind),dyt2,lm(1,1,jatom))
               call suml(lm1,dydp(lm1,ind),dyp,lm(1,1,jatom))
               call suml(lm1,dydp2(lm1,ind),dyp2,lm(1,1,jatom))
               call suml(lm1,dydtp(lm1,ind),dytp,lm(1,1,jatom))
            end do
         endif
!
         IF(iatj.GT.0.AND.LLMM.GT.1) THEN                  
            CALL SUMFAC(llmm,YY,IND,LM(1,1,jatom),LMMAX(JATOM))
!
!         if(iex.gt.10) then
            if(GGA_SWITCH.eqv..true.) then
               call sumfac(llmm,dydt,ind,LM(1,1,jatom),LMMAX(JATOM))
               call sumfac(llmm,dydt2,ind,LM(1,1,jatom),LMMAX(JATOM))
               call sumfac(llmm,dydp,ind,LM(1,1,jatom),LMMAX(JATOM))
               call sumfac(llmm,dydp2,ind,LM(1,1,jatom),LMMAX(JATOM))
               call sumfac(llmm,dydtp,ind,LM(1,1,jatom),LMMAX(JATOM))
            endif
!        
         ENDIF
!                                                                       
  1   CONTINUE
      
!                                                                       
!========================================================
!         derivatives dCLM/dR and ddCLM/ddR
!========================================================
      if(GGA_SWITCH.eqv..true.) then
         call dfundr(jatom,jri(jatom),lmmax(jatom),clmsp,dcsp,ddcsp,spin,LM)
      end if
!========================================================
!
!                                                                       
!========================================================
!         derivatives d(SLATER_TAU_LM)/dR
!                 and dd(SLATER_TAU_LM)/ddR
!========================================================
      if(MGGA_TAU_SWITCH) then
         call dfundr(jatom,jri(jatom),lmmax(jatom),slater_tau_lm,dtau_dr,d2tau_dr2,spin,LM)
      end if
!========================================================
!
!     Establish Weight*YY, only values used here
      NCond = 0
      DO LM1=1,LLMM
        if(UseLM(LM1))then
                NCOND=NCOND+1
                YYB(1:NPT1,NCOND) = Weight(1:NPT1)*YY(LM1,1:NPT1)
        else
!               Zero ones not used
                dydt(lm1,:)     =0.D0
                dydt2(lm1,:)    =0.D0
                dydp(lm1,:)     =0.D0
                dydp2(lm1,:)    =0.D0
                dydtp(lm1,:)    =0.D0
                yy(lm1,:)       =0.D0
                dcsp(lm1,:,:)   =0.D0
                ddcsp(lm1,:,:)  =0.D0
        endif
      enddo
!
      LWORK=12*NSPD
!     Calculate Condition Number from singular values
      call DGESVD( 'N', 'N', NPT1, NCOND, YYB, nspd, Sing, U, 1, VT, 1, &
                     SWORK, LWORK, INFO )
      RUp = Maxval(Sing(1:NCOND))
      RDn = MinVal(Sing(1:NCOND))
      RCond = RUp/RDn
      write(6,'(A,ES12.5)')' Condition Number of XCPOT1 ',RCond


!     Store YY*weight product, we don't need to keep doing it.
      DO L=1,LMMAX(JATOM)
         YYB(1:npt1,L)=YY(L,1:npt1)*weight(1:npt1)
      enddo
      
!
!     LOOP OVER ALL RADIAL MESH POINTS                                  
!
      init=1
!
!$omp parallel do &
!$omp& private(bbbu,bbbd,rp,rpp,k, &
!$omp& qrhou,qrhod,ft,fu,fd,tautot,tauup,taudn,vxuold,vxdold,vxtold, &
!$omp& vxunonloc,vxdnonloc,gx,gxd,gxu,gy,gyd,gyu,gmag,gmagd,gmagu, &
!$omp& g2,g2d,g2u,ggx,ggxd,ggxu,ggy,ggyd,ggyu,ggz,ggzd,ggzu, &
!$omp& gdggt,gdggu,gdggd,vxcd,vxcu,vhup,vhdn,vhtot,gz,gzu,gzd, &
!$omp& d2ndru,d2ndrd,d2ndrt,dndtu,dndtd,dndtt,d2ndtu,d2ndtd,d2ndtt, &
!$omp& dndpu,dndpd,dndpt,d2ndpu,d2ndpd,d2ndpt,dndtpu,dndtpd,dndtpt, &
!$omp& dndrtu,dndrtd,dndrtt,dndrpu,dndrpd,dndrpt,dndru,dndrd,dndrt, &
!$omp& dtau_dt_u, dtau_dt_d, dtau_dt_t, d2tau_dt2_u , d2tau_dt2_d, &
!$omp& d2tau_dt2_t, dtau_dp_u, dtau_dp_d, dtau_dp_t, d2tau_dp2_u, &
!$omp& d2tau_dp2_d, d2tau_dp2_t, d2tau_dtdp_u, d2tau_dtdp_d, &
!$omp& d2tau_dtdp_t, dtau_dr_u, dtau_dr_d, dtau_dr_t, d2tau_drdt_u, &
!$omp& d2tau_drdt_d, d2tau_drdt_t, d2tau_drdp_u, d2tau_drdp_d, &
!$omp& d2tau_drdp_t, d2tau_d2r_t, dtau_dx, dtau_dy, dtau_dz, &
!$omp& dtau_mag, d2tau, d2tau_x, d2tau_y, d2tau_z, dtau_d2tau, d2tau_d2r_u, &
!$omp& dtau_dx_u, dtau_dy_u, dtau_dz_u, dtau_mag_u, d2tau_u, &
!$omp& d2tau_x_u, d2tau_y_u, d2tau_z_u, dtau_d2tau_u, d2tau_d2r_d, &
!$omp& dtau_dx_d, dtau_dy_d, dtau_dz_d, dtau_mag_d, d2tau_d, d2tau_x_d, &
!$omp& d2tau_y_d, d2tau_z_d, dtau_d2tau_d, dtau_drho, dtau_u_drho_u, &
!$omp& dtau_u_drho_d, dtau_d_drho_u, dtau_d_drho_d, d2tau_dr2_u, &
!$omp& d2tau_dr2_d, d2tau_dr2_t, &
!$omp& excd,excu,l,lm1,vxcua,vxcda,excua,excda,test,sigma,exu,exd,y,t, &
!$omp& vtauc_up,vtauc_dn,vlaplc_up,vlaplc_dn, &
!$omp& vtau_rad_up,vtau_rad_dn,vlapl_rad_up,vlapl_rad_dn, &
!$omp& vtau_ua,vtau_da,vlapl_ua,vlapl_da, gxyzu1, gxyzd1, gxyzud1, &
!$omp& sig_str_u,sig_str_d,sig_str_ud) &
!$omp& lastprivate(vxcsu,vxcsd,bbu,bbd, &
!$omp& vtau_bu,vtau_bd,vlapl_bu,vlapl_bd) &
!$omp& copyin(xc1,xc2,xc3,xc4,xcconst) &
!$omp& reduction(max:dmaxup,dmaxdn) schedule(guided)
! !!!$omp& tmp,c1,c2,c3,c4, tmp1,tmp2,tmp3,tmp4
      DO 34 ir=1,JRI(JATOM)               
!
!if(ir.lt.20) print*,r(ir),dcsp(1,ir,1),ddcsp(1,ir,1)
!if(ir.gt.JRI(JATOM)-20) print*,r(ir),dcsp(1,ir,1),ddcsp(1,ir,1)
         bbu=0.0D0
         bbbu=0.0D0
         vxcsu=0.0D0
         exu=0.0D0
         vtau_bu=0.0D0
         vlapl_bu=0.0D0
         if(spin)then
            bbd=0.0D0
            bbbd=0.0D0
            vxcsd=0.0D0
            exd=0.0D0
            vtau_bd=0.0D0
            vlapl_bd=0.0D0
         endif
!
         rp=r(ir)
         rpp=1.D0/(r(ir)*r(ir))

         DO 77 k=1,npt1
            
            qrhou(1:LMMAX(JATOM))=clmsp(ir,1:LMMAX(JATOM),jatom,1)
            if(spin)qrhod(1:LMMAX(JATOM))=clmsp(ir,1:LMMAX(JATOM),jatom,2)
            call srolyl(qrhou,qrhod,yy(1,k),sqfpi,iatj,llmm, &
                 fu,fd,ft,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
            ft=ft*RPP
            fu=fu*RPP
            fd=fd*RPP
!
            if(VRESP_SWITCH) then
               
               qrhou(1:LMMAX(JATOM))=vclmsp(ir,1:LMMAX(JATOM),jatom,1)
               if(spin)qrhod(1:LMMAX(JATOM))=vclmsp(ir,1:LMMAX(JATOM),jatom,2)
               
               call srolyl(qrhou,qrhod,yy(1,k),sqfpi,iatj,llmm, &
                    tauup,taudn,tautot,LM(1,1,jatom),LMMAX(JATOM),nospin)
               tautot=-tautot*RPP
               tauup=-tauup*RPP
               taudn=-taudn*RPP
               
               qrhou(1:LMMAX(JATOM))=v(ir,1:LMMAX(JATOM),jatom)
               if(spin)qrhod(1:LMMAX(JATOM))=v(ir,1:LMMAX(JATOM),jatom)
               qrhou(1)=qrhou(1)*RTPI
               if(spin)qrhod(1)=qrhod(1)*RTPI
               call srolyl(qrhou,qrhod,yy(1,k),sqfpi,iatj,llmm, &
                    vhup,vhdn,vhtot,LM(1,1,jatom),LMMAX(JATOM),nospin)
!.......v was not multiplied with r**2
               if ((VX_SWITCH(1:7) .eq. 'VX_GLLB') .or. (EX_SWITCH(1:12) .eq. 'EX_DELTAGLLB')) then
                  tautot=-tautot
                  tauup=-tauup
                  taudn=-taudn
               else ! subtract v_hartree*rho from sum_i of epsilon_i*rho_i
                  tautot=tautot-vhtot*ft
                  tauup=tauup -vhup*fu
                  taudn=taudn -vhdn*fd
               endif
!           Slater form of kinetic energy density is read from file
            else if(read_tau_switch) then
               qrhou(1:lmmax(jatom)) = slater_tau_lm(ir,1:lmmax(jatom),jatom,1)
               if(spin)qrhod(1:LMMAX(JATOM)) = slater_tau_lm(ir,1:lmmax(jatom),jatom,2)

               call srolyl(qrhou,qrhod,yy(1,k),sqfpi,iatj,llmm, &
                    tauup,taudn,tautot,lm(1,1,jatom),lmmax(jatom),nospin)
               
               tautot=tautot*RPP
               tauup=tauup*RPP
               taudn=taudn*RPP
            endif
!
!--------------------------------
!        if (iex .eq. 28 .or. iex .eq. 30) then
        if((VX_SWITCH.eq.'VX_BR89') .or. (VX_SWITCH.eq.'VX_SLATER') .or. &
           (VX_SWITCH.eq.'VX_MBJ') .or. (VX_SWITCH.eq.'VX_SMBJ') .or. &
           (VX_SWITCH.eq.'VX_GBJ') .or. (VX_SWITCH.eq.'VX_KLI') .or. &
           (VX_SWITCH.eq.'VX_R2V') .or. (VDISP_SWITCH.eq.'VDISP_NLVDW') .or. &
           (VHALF_SWITCH.ne.'VHALF_NONE') .or. (VX_SWITCH.eq.'VX_LMBJ') .or. &
           (id_func_x .eq. 207) .or. (id_func_x .eq. 208) .or. (id_func_x .eq. 209)) then

          qrhou(1:LMMAX(JATOM)) = vxcold(ir,1:LMMAX(JATOM),jatom,1)
          if (spin) qrhod(1:LMMAX(JATOM)) = vxcold(ir,1:LMMAX(JATOM),jatom,2)

          call srolyl(qrhou,qrhod,yy(1,k),1d0,iatj,llmm, &
                      vxuold,vxdold,vxtold,LM(1,1,jatom),LMMAX(JATOM),nospin)

          vxuold = vxuold*rpp
          vxdold = vxdold*rpp

        endif

        if((VX_SWITCH.eq.'VX_SLATER') .or. (VX_SWITCH.eq.'VX_SMBJ') .or. &
           (VX_SWITCH.eq.'VX_KLI') .or. (EX_SWITCH(1:12).eq.'EX_DELTAGLLB') .or. &
           (VX_SWITCH.eq.'VX_LMBJ')) then

          qrhou(1:LMMAX(JATOM)) = vxcnonloc(ir,1:LMMAX(JATOM),jatom,1)
          if (spin) qrhod(1:LMMAX(JATOM)) = vxcnonloc(ir,1:LMMAX(JATOM),jatom,2)

          call srolyl(qrhou,qrhod,yy(1,k),1d0,iatj,llmm, &
                      vxunonloc,vxdnonloc,vxtnonloc,LM(1,1,jatom),LMMAX(JATOM),nospin)

          vxunonloc = vxunonloc*rpp
          vxdnonloc = vxdnonloc*rpp

        endif

!        if(iex.gt.10) then
        if(GGA_SWITCH.eqv..true.) then
!
           qrhou(1:LMMAX(JATOM))=ddcsp(1:LMMAX(JATOM),ir,1)
           if(spin)qrhod(1:LMMAX(JATOM))=ddcsp(1:LMMAX(JATOM),ir,2)
           call srolyl(qrhou,qrhod,yy(1,k),sqfpi,iatj,llmm, &
                d2ndru,d2ndrd,d2ndrt,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           qrhou(1:LMMAX(JATOM))=clmsp(ir,1:LMMAX(JATOM),jatom,1)*rpp
           if(spin)qrhod(1:LMMAX(JATOM))=clmsp(ir,1:LMMAX(JATOM),jatom,2)*rpp
           
           call srolyl(qrhou,qrhod,dydt(1,k),sqfpi,iatj,llmm, &
                dndtu,dndtd,dndtt,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           call srolyl(qrhou,qrhod,dydt2(1,k),sqfpi,iatj,llmm, &
                d2ndtu,d2ndtd,d2ndtt,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           call srolyl(qrhou,qrhod,dydp(1,k),sqfpi,iatj,llmm, &
                dndpu,dndpd,dndpt,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           call srolyl(qrhou,qrhod,dydp2(1,k),sqfpi,iatj,llmm, &
                d2ndpu,d2ndpd,d2ndpt,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           call srolyl(qrhou,qrhod,dydtp(1,k),sqfpi,iatj,llmm, &
                dndtpu,dndtpd,dndtpt,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           if(spin)qrhod(1:LMMAX(JATOM))=dcsp(1:LMMAX(JATOM),ir,2)
           call srolyl(dcsp(1,ir,1),qrhod,dydt(1,k),sqfpi,iatj,llmm, &
                dndrtu,dndrtd,dndrtt,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           call srolyl(dcsp(1,ir,1),qrhod,dydp(1,k),sqfpi,iatj,llmm, &
                dndrpu,dndrpd,dndrpt,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           call srolyl(dcsp(1,ir,1),qrhod,yy(1,k),sqfpi,iatj,llmm, &
                dndru,dndrd,dndrt,LM(1,1,jatom),LMMAX(JATOM),nospin)

!
!--------------------------------         
        endif

        if(GGA_SWITCH.eqv..true.) then

           call grans(rp,qsin(k),qcos(k),dndrt,dndtt,dndpt, &
                d2ndrt,d2ndtt,d2ndpt,dndrtt,dndrpt,dndtpt, &
                gx,gy,gz,gmag,g2,ggx,ggy,ggz,gdggt,ir)
           call grans(rp,qsin(k),qcos(k),dndru,dndtu,dndpu, &
                d2ndru,d2ndtu,d2ndpu,dndrtu,dndrpu,dndtpu, &
                gxu,gyu,gzu,gmagu,g2u,ggxu,ggyu,ggzu,gdggu,ir)
           if(spin)then
              call grans(rp,qsin(k),qcos(k),dndrd,dndtd,dndpd, &
                   d2ndrd,d2ndtd,d2ndpd,dndrtd,dndrpd,dndtpd, &
                   gxd,gyd,gzd,gmagd,g2d,ggxd,ggyd,ggzd,gdggd,ir)
           else
              gxd=gxu
              gyd=gyu
              gzd=gzu
              gmagd=gmagu
              g2d=g2u
              ggxd=ggxu
              ggyd=ggyu
              ggzd=ggzu
              gdggd=gdggu
           endif
        endif

        if(MGGA_TAU_SWITCH) then
!          Compute derivatives of tau
           
!          tau file is required
           if (.not. read_tau_switch) then
              call outerr('LAPW0', 'Tau file required for tau-dependent meta-gga`s')
              call abort_parallel
              stop "tau file not provided"
           end if
           
!           second derivatives not needed
           qrhou(1:LMMAX(JATOM))=d2tau_dr2(1:LMMAX(JATOM),ir,1)
           if(spin)qrhod(1:LMMAX(JATOM))=d2tau_dr2(1:LMMAX(JATOM),ir,2)
           call srolyl(qrhou,qrhod,yy(1,k),sqfpi,iatj,llmm, &
                       d2tau_dr2_u,d2tau_dr2_d,d2tau_dr2_t,LM(1,1,jatom),LMMAX(JATOM),nospin)

           qrhou(1:LMMAX(JATOM))=slater_tau_lm(ir,1:LMMAX(JATOM),jatom,1)*rpp
           if(spin)qrhod(1:LMMAX(JATOM))=slater_tau_lm(ir,1:LMMAX(JATOM),jatom,2)*rpp

!          derivatives w.r.t. theta
           call srolyl(qrhou,qrhod,dydt(1,k),sqfpi,iatj,llmm, &
                dtau_dt_u,dtau_dt_d,dtau_dt_t,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           call srolyl(qrhou,qrhod,dydt2(1,k),sqfpi,iatj,llmm, &
                d2tau_dt2_u,d2tau_dt2_d,d2tau_dt2_t,LM(1,1,jatom),LMMAX(JATOM),nospin)
           
!          derivatives w.r.t. phi
           call srolyl(qrhou,qrhod,dydp(1,k),sqfpi,iatj,llmm, &
                dtau_dp_u,dtau_dp_d,dtau_dp_t,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
           call srolyl(qrhou,qrhod,dydp2(1,k),sqfpi,iatj,llmm, &
                d2tau_dp2_u,d2tau_dp2_d,d2tau_dp2_t,LM(1,1,jatom),LMMAX(JATOM),nospin)
!
!          mixed second derivatives w.r.t. theta, phi
           call srolyl(qrhou,qrhod,dydtp(1,k),sqfpi,iatj,llmm, &
                d2tau_dtdp_u,d2tau_dtdp_d,d2tau_dtdp_t,LM(1,1,jatom),LMMAX(JATOM),nospin)
           
!          derivatives w.r.t. radial coordinate
           qrhou(1:LMMAX(JATOM))=dtau_dr(1:LMMAX(JATOM),ir,1)
           if(spin)qrhod(1:LMMAX(JATOM))=dtau_dr(1:LMMAX(JATOM),ir,2)
!          derivative w.r.t r
           call srolyl(qrhou,qrhod,yy(1,k),sqfpi,iatj,llmm, &
                dtau_dr_u,dtau_dr_d,dtau_dr_t,LM(1,1,jatom),LMMAX(JATOM),nospin)
!          mixed second derivative r and theta
           call srolyl(qrhou,qrhod,dydt(1,k),sqfpi,iatj,llmm, &
                d2tau_drdt_u,d2tau_drdt_d,d2tau_drdt_t,LM(1,1,jatom),LMMAX(JATOM),nospin)
!          mixed second derivative r and phi
           call srolyl(qrhou,qrhod,dydp(1,k),sqfpi,iatj,llmm, &
                d2tau_drdp_u,d2tau_drdp_d,d2tau_drdp_t,LM(1,1,jatom),LMMAX(JATOM),nospin)

!--------------------------------         
!         
!          laplacian and gradient 
           call grans(rp,qsin(k),qcos(k),dtau_dr_t,dtau_dt_t,dtau_dp_t, &
                d2tau_d2r_t,d2tau_dt2_t,d2tau_dp2_t,d2tau_drdt_t,d2tau_drdp_t, d2tau_dtdp_t, &
                dtau_dx, dtau_dy, dtau_dz, dtau_mag, &
                d2tau, d2tau_x, d2tau_y, d2tau_z, dtau_d2tau, ir)
           call grans(rp,qsin(k),qcos(k),dtau_dr_u,dtau_dt_u,dtau_dp_u, &
                d2tau_d2r_u,d2tau_dt2_u,d2tau_dp2_u,d2tau_drdt_u,d2tau_drdp_u, d2tau_dtdp_u, &
                dtau_dx_u, dtau_dy_u, dtau_dz_u, dtau_mag_u, &
                d2tau_u, d2tau_x_u, d2tau_y_u, d2tau_z_u, dtau_d2tau_u, ir)
           if(spin)then
              call grans(rp,qsin(k),qcos(k),dtau_dr_d,dtau_dt_d,dtau_dp_d, &
                   d2tau_d2r_d,d2tau_dt2_d,d2tau_dp2_d,d2tau_drdt_d,d2tau_drdp_d, d2tau_dtdp_d, &
                   dtau_dx_d, dtau_dy_d, dtau_dz_d, dtau_mag_d, &
                   d2tau_d, d2tau_x_d, d2tau_y_d, d2tau_z_d, dtau_d2tau_d, ir)
           else
              dtau_dx_d=dtau_dx_u
              dtau_dy_d=dtau_dy_u
              dtau_dz_d=dtau_dz_u
              dtau_mag_d=dtau_mag_u
              d2tau_d=d2tau_u
              d2tau_x_d=d2tau_x_u
              d2tau_y_d=d2tau_y_u
              d2tau_z_d=d2tau_z_u
           endif
 
!          calculate dot product grad(tau).grad(rho)
           dtau_drho = dtau_dx*gx + dtau_dy*gy + dtau_dz*gz
           dtau_u_drho_u = dtau_dx_u*gxu + dtau_dy_u*gyu + dtau_dz_u*gzu
           dtau_u_drho_d = dtau_dx_u*gxd + dtau_dy_u*gyd + dtau_dz_u*gzd
           dtau_d_drho_u = dtau_dx_d*gxu + dtau_dy_d*gyu + dtau_dz_d*gzu
           dtau_d_drho_d = dtau_dx_d*gxd + dtau_dy_d*gyd + dtau_dz_d*gzd
        end if

        if(MGGA_LAPL_SWITCH) then
!          Compute derivatives of laplacian of density
           call outerr('LAPW0', 'Laplacian-dependent meta-gga`s&
                & are not yet implemented')
           call abort_parallel
           stop "not yet implemented"
        end if
     
        CALL VXCLM2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU,&
             FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,&
             VXCU,VXCD,excu,excd,&
             vtauc_up, vtauc_dn, vlaplc_up, vlaplc_dn,&
             dtau_drho,dtau_u_drho_u,dtau_u_drho_d,dtau_d_drho_u,dtau_d_drho_d,&
             VXUOLD,VXDOLD,VXUNONLOC,VXDNONLOC,&
             GGX,GGY,GGZ,GDGGT,gdggu,gdggd,tauup,taudn,&
             ir,CutMuffin,JSPIN,JATOM,sig_str_u,sig_str_d,sig_str_ud)
!     if(ir.ge.780.and. k.lt.3) then
!          write(6,*)'xcpot1',sig_str_u,sig_str_d,sig_str_ud       
!          write(6,*)FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU
!          write(6,*)FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD
!          write(6,*)   VXCU,VXCD,excu,excd
!          write(6,*) GGX,GGY,GGZ,GDGGT,gdggu,gdggd
!  endif
        !vtauc_up = 0.0d0
        !vtauc_dn = 0.0d0
!        
!                 
         IF(STRESS) THEN
           gxyzu1(k) = sig_str_u
           IF(spin) then
              gxyzd1(k)  = sig_str_d  
              gxyzud1(k) = sig_str_ud
           ELSE 
              gxyzd1(k)  = gxyzu1(k)              
              gxyzud1(k) = 0.d0
           ENDIF             
         ENDIF  
!
        if(imfield.gt.0) then
           vxcu=vxcu - Bextene
           vxcd=vxcd + Bextene
        endif
!
!       if(ir.eq.jri(jatom)) then
         vxcsu(k)=vxcu
         vxcsd(k)=vxcd
         exu  (k)=excu
         exd  (k)=excd
         
         if (MGGA_TAU_SWITCH) then
            vtau_rad_up(k) = vtauc_up
            vtau_rad_dn(k) = vtauc_dn
!            print*, vtauc_up
         end if
         
         if (MGGA_LAPL_SWITCH) then
            vlapl_rad_up(k) = vlaplc_up
            vlapl_rad_dn(k) = vlaplc_dn
         end if
         
         if(ir.eq.jri(jatom).and.iprint.eq.2) then
!this write will kill threaded performance if iprint=2, and possibly
!prevent really optimized scheduling even if not
                 write(789,'(29f10.5)') spt(1,k),spt(2,k),spt(3,k),fu,GXU,GYU,GZU,&
                                 GMAGU,G2U,GGXU,GGYU,GGZU,fd,GXD,GYD,GZD,GMAGD,&
                                 G2D,GGXD,GGYD,GGZD,excu,excd,vxcu,vxcd
         endif
! *** LDM Change
77      continue
!                                                                       
! Kahan summation dot product to evaluate the mean
        VXCUA = ddotk(npt1,YYB,1,VXCSU,1)                
        if(spin)VXCDA = ddotk(npt1,YYB,1,VXCSD,1)
        EXCUA = ddotk(npt1,YYB,1,EXU  ,1)
        if(spin)EXCDA = ddotk(npt1,YYB,1,EXD  ,1)

! Subtract the mean (largest term)
! Incorrect mean was used before, does not matter much but this is more rigorous
        VXCSU = VXCSU - VXCUA/Sum0
        if(spin)VXCSD = VXCSD - VXCDA/Sum0
        EXU   = EXU   - EXCUA /Sum0
        if(spin)EXD   = EXD   - EXCDA /Sum0

! Kahan summation dot product after subtraction of mean (large component)
        DO L=1,LMMAX(JATOM)                                                   
                if(UseLM(L))then
                        BBU(L)  = ddotk(npt1,YYB(1,L),1,VXCSU,1)
                        if(spin)BBD(L)  = ddotk(npt1,YYB(1,L),1,VXCSD,1)
                        BBBU(L) = ddotk(npt1,YYB(1,L),1,EXU  ,1)
                        if(spin)BBBD(L) = ddotk(npt1,YYB(1,L),1,EXD  ,1)
                endif
        enddo

!      Add back the base term (Note Sum0 does not belong here)
       bbu(1)=bbu(1)+vxcuA 
       if(spin)bbd(1)=bbd(1)+vxcdA
       bbbu(1)=bbbu(1)+excuA
       if(spin)bbbd(1)=bbbd(1)+excdA

! Add back the mean for diagnostics, slightly sloppy
        VXCSU = VXCSU + VXCUA/Sum0
        if(spin)VXCSD = VXCSD + VXCDA/Sum0

       if(debug)then
          if(nospin)then
            DO LM1=1,LLMM
               bbd (lm1)=bbu (lm1)
               bbbd(lm1)=bbbu(lm1)
            enddo
          endif

          do LM1=1,LLMM
                write(6,1551)IR,LM1,bbd(LM1),bbu(LM1),bbbu(LM1),bbbd(LM1)
          enddo
       endif

1551   format(2I5,4D26.18)



!
!    Integration over the angular grid for mgga parts
!
       if (MGGA_TAU_SWITCH) then
          vtau_ua = ddotk(npt1,yyb,1,vtau_rad_up,1)
          if(spin) vtau_da = ddotk(npt1,yyb,1,vtau_rad_dn,1)
          vtau_rad_up = vtau_rad_up - vtau_ua/sum0
          vtau_rad_dn = vtau_rad_dn - vtau_da/sum0

! Kahan summation dot product after subtraction of mean (large component)
          do L=1,LMMAX(JATOM)                                                   
             if(UseLM(L)) then
                vtau_bu(L)  = ddotk(npt1,YYB(1,L),1,vtau_rad_up,1)
                if(spin) vtau_bd(L)  = ddotk(npt1,YYB(1,L),1,vtau_rad_dn,1)
             endif
          enddo

!      Add back the base term (Note Sum0 does not belong here)
          vtau_bu(1)=vtau_bu(1)+vtau_ua 
          if(spin)vtau_bd(1)=vtau_bd(1)+vtau_da

! Add back the mean for diagnostics, slightly sloppy
          vtau_rad_up = vtau_rad_up + vtau_ua/sum0
          vtau_rad_dn = vtau_rad_dn + vtau_da/sum0

          if(debug) then
             if(nospin) vtau_bd=vtau_bu
             do LM1=1,LLMM
                write(6,1552)IR,LM1,vtau_bd(LM1),vtau_bu(LM1)
             enddo
          endif

1552      format(2I5,2D26.18)
          
     endif

     if (MGGA_LAPL_SWITCH) then
         vlapl_ua = ddotk(npt1,yyb,1,vlapl_rad_up,1)
          if(spin) vlapl_da = ddotk(npt1,yyb,1,vlapl_rad_dn,1)
          vlapl_rad_up = vlapl_rad_up - vlapl_ua/sum0
          vlapl_rad_dn = vlapl_rad_dn - vlapl_da/sum0

! Kahan summation dot product after subtraction of mean (large component)
          do L=1,LMMAX(JATOM)                                                   
             if(UseLM(L)) then
                vlapl_bu(L)  = ddotk(npt1,YYB(1,L),1,vlapl_rad_up,1)
                if(spin) vlapl_bd(L)  = ddotk(npt1,YYB(1,L),1,vlapl_rad_dn,1)
             endif
          enddo

!      Add back the base term (Note Sum0 does not belong here)
          vlapl_bu(1)=vlapl_bu(1)+vlapl_ua 
          if(spin)vlapl_bd(1)=vlapl_bd(1)+vlapl_da

! Add back the mean for diagnostics, slightly sloppy
          vlapl_rad_up = vlapl_rad_up + vlapl_ua/sum0
          vlapl_rad_dn = vlapl_rad_dn + vlapl_da/sum0

          if(debug) then
             if(nospin) vlapl_bd=vlapl_bu
             do LM1=1,LLMM
                write(6,1552)IR,LM1,vlapl_bd(LM1),vlapl_bu(LM1)
             enddo
          endif
     endif
     
!===========for the stress tensor=====
!      IF(STRESS) THEN
!           sumU = ddotk(npt1,yyb,1,gxyzu1,1)
!           if(spin) sumD = ddotk(npt1,yyb,1,gxyzd1,1)
!           if(spin) sumUD = ddotk(npt1,yyb,1,gxyzud1,1)
!           gxyzu1(1:npt1) = gxyzu1(1:npt1) - sumU   /sum0
!           gxyzd1(1:npt1) = gxyzd1(1:npt1) - sumD   /sum0
!           gxyzud1(1:npt1)= gxyzud1(1:npt1) - sumUD /sum0
!           do L =1, LMMAX(JATOM)
!              if(UseLM(L)) then
!                sigma_up(L) = ddotk(npt1,YYB(1,L),1,gxyzu1,1)
!                if(spin) sigma_dn(L) = ddotk(npt1,YYB(1,L),1,gxyzd1,1)
!                if(spin) sigma_updn(L) = ddotk(npt1,YYB(1,L),1,gxyzud1,1)
!              endif
!           enddo
!           sigma_up(1) = sigma_up(1) + sumU
!           if(spin) sigma_dn(1) = sigma_dn(1)+sumD
!           if(spin) sigma_updn(1) = sigma_updn(1)+sumUD
!           if(debug) then
!              if(nospin) sigma_dn   = sigma_up
!              if(nospin) sigma_updn = 0.d0
!           endif
!      ENDIF
     
!
!================= VXC, UP ==================
      SIGMA=0.D0
      DO k=1,npt1
        TEST=0.D0
!$omp simd reduction(+:TEST)
        DO LM1=1,LLMM
          TEST=TEST+bbu(LM1)*YY(LM1,k)
        ENDDO
        if(dmaxup.lt.abs(vxcsu(k)-TEST)) dmaxup= abs(vxcsu(k)-TEST)
        SIGMA=SIGMA+weight(k)*(vxcsu(k)-TEST)**2
      ENDDO
      SIGMA=SQRT(SIGMA/wnorm)
      IF(iatj.LT.0) THEN
          DO LM1=1,LLMM
            VXC(ir,LM1,JATOM,1)=bbu(LM1)
          ENDDO
      ELSE
          CALL cub_xc_back(vxc(1,1,jatom,1),ir,llmm, &
                        lm(1,1,jatom),bbu,lmmax(jatom))
      ENDIF                                                         
!$omp critical (sigmamaxup_critical)
      if(sigmamaxup.lt.sigma) then
        sigmamaxup=sigma
        irsigmaup=ir
      endif
!$omp end critical (sigmamaxup_critical)
      if(ir.eq.jri(jatom)) sigmaj(1)=sigma
!
!================= VXC, DN ==================
!
      if(spin)then
        SIGMA=0.D0
        DO k=1,npt1
          TEST=0.D0
!$omp simd reduction(+:TEST)
          DO LM1=1,LLMM
            TEST=TEST+bbd(LM1)*YY(LM1,k)
          ENDDO
          if(dmaxdn.lt.abs(vxcsd(k)-TEST)) dmaxdn= abs(vxcsd(k)-TEST)
          SIGMA=SIGMA+weight(k)*(vxcsd(k)-TEST)**2
        ENDDO
        SIGMA=SQRT(SIGMA/wnorm)
        IF(iatj.LT.0) THEN
          DO LM1=1,LLMM
            VXC(ir,LM1,JATOM,2)=bbd(LM1)
          ENDDO
        ELSE
          CALL cub_xc_back(vxc(1,1,jatom,2),ir,llmm, &
                           lm(1,1,jatom),bbd,lmmax(jatom))
        ENDIF

!$omp critical (sigmamaxdn_critical)
        if(sigmamaxdn.lt.sigma) then
          sigmamaxdn=sigma
          irsigmadn=ir
        endif
!$omp end critical (sigmamaxdn_critical)
        if(ir.eq.jri(jatom)) sigmaj(2)=sigma
      endif
!                                                                       
!================= EXC, UP ==================
!
      IF(iatj.LT.0) THEN                                        
       DO 19 LM1=1,LLMM                                                  
  19   EXC(ir,LM1,JATOM,1)=bbbu(LM1)                                   
      ELSE                                                              
         CALL cub_xc_back(exc(1,1,jatom,1),ir,llmm, &
                        lm(1,1,jatom),bbbu,lmmax(jatom))
      ENDIF                                                         
!                                                                       
!================= EXC, DN ==================
!
      if(spin)then
        IF(iatj.LT.0) THEN
                DO 5019 LM1=1,LLMM                                            
 5019           EXC(ir,LM1,JATOM,2)=bbbd(LM1)                                
        ELSE                                                              
                CALL cub_xc_back(exc(1,1,jatom,2),ir,llmm, &
                        lm(1,1,jatom),bbbd,lmmax(jatom))
        ENDIF                                                         
      endif

!                                                                       
!================= VTAU ==================
!
      IF (MGGA_TAU_SWITCH) THEN
         IF(iatj.LT.0) THEN                                        
            DO  LM1=1,LLMM
               vtau(ir,LM1,JATOM,1)=vtau_bu(LM1)
            end DO
         ELSE                                                              
            CALL cub_xc_back(vtau(1,1,jatom,1),ir,llmm, &
                 lm(1,1,jatom),vtau_bu,lmmax(jatom))
         ENDIF

         if(spin)then
            IF(iatj.LT.0) THEN
               DO LM1=1,LLMM                                            
                  vtau(ir,LM1,JATOM,2)=vtau_bd(LM1)
               end DO
            ELSE                                                              
               CALL cub_xc_back(vtau(1,1,jatom,2),ir,llmm, &
                    lm(1,1,jatom),vtau_bd,lmmax(jatom))
            ENDIF
         endif
      END IF
      
      
!                                                                       
!================= VLAPL ==================
!
      IF (MGGA_LAPL_SWITCH) THEN
         IF(iatj.LT.0) THEN                                        
            DO  LM1=1,LLMM
               vlapl(ir,LM1,JATOM,1)=vlapl_bu(LM1)
            end DO
         ELSE                                                              
            CALL cub_xc_back(vlapl(1,1,jatom,1),ir,llmm, &
                 lm(1,1,jatom),vlapl_bu,lmmax(jatom))
         ENDIF

         if(spin)then
            IF(iatj.LT.0) THEN
               DO LM1=1,LLMM                                            
                  vlapl(ir,LM1,JATOM,2)=vlapl_bd(LM1)
               end DO
            ELSE                                                              
               CALL cub_xc_back(vlapl(1,1,jatom,2),ir,llmm, &
                    lm(1,1,jatom),vlapl_bd,lmmax(jatom))
            ENDIF
         endif
      END IF

!
!
!=============FOR GGA STRESS====
!
      IF(STRESS) THEN
        If(GGA_SWITCH.eqv..true.) THEN 
        
         call kahan_summ_gxyz1(gxyzu1,gxyzd1,gxyzud1,npt1,LLMM,weight,YY, &
              sum0,spin,iatj,nat,LM,LMMAX,jatom,ir,debug,sigma_str) 

!               print*, sum0,gxyzu1(npt1),weight(npt1)
!           IF(iatj.LT.0) THEN
!             DO lm1 = 1, llmm
!                sigma_str(ir,lm1,jatom,1) = sigma_up(lm1)         
!             ENDDO
!           ELSE
!             call cub_xc_back(sigma_str(1,1,jatom,1),ir,llmm,&
!                   lm(1,1,jatom),sigma_up,lmmax(jatom))      
!           ENDIF
!           
!           IF(spin) THEN
!             IF(iatj.LT.0) THEN
!               DO lm1 = 1, llmm
!                  sigma_str(ir,lm1,jatom,2) = sigma_dn(lm1)
!                  sigma_str(ir,lm1,jatom,3) = sigma_updn(lm1)
!               ENDDO
!             ELSE
!               call cub_xc_back(sigma_str(1,1,jatom,2),ir,llmm,&
!                   lm(1,1,jatom),sigma_dn,lmmax(jatom))        
!               call cub_xc_back(sigma_str(1,1,jatom,2),ir,llmm,&
!                   lm(1,1,jatom),sigma_updn,lmmax(jatom))            
!             ENDIF
!           ENDIF
!           
!           IF(nospin) THEN
!             sigma_str(ir,1:lmmax(jatom),jatom,2) = sigma_str(ir,1:lmmax(jatom),jatom,1)      
!             sigma_str(ir,1:lmmax(jatom),jatom,3) = 0.d0
!           ENDIF
          
        ENDIF        
      ENDIF
      
34    CONTINUE          
!********************************************************
!      LOOP (34) OVER ALL RADIAL MESH COMPLETED
!                                                                       
!
!================ OUTPUT ===================
!
! Fitting is not (yet) done for mgga parts of potentials!!

      IF(iatj.LT.0) WRITE(6,1540)
      DO k=1,npt1
        TEST=0.D0
!$omp simd reduction(+:TEST)
        DO LM1=1,LLMM
          TEST=TEST+bbu(LM1)*YY(LM1,k)
        ENDDO
        if((k.eq.1).or.((iprint.gt.0).and.(MOD(k,2).EQ.0))) then
          IF(iatj.LT.0) THEN
            WRITE(6,1550)k,vxcsu(k),TEST,vxcsu(k)-TEST
          ELSE
            WRITE(6,*) k,vxcsu(k),TEST
          ENDIF
        endif
      ENDDO
      WRITE(6,'(a,i4,a,10x,2e20.5)') 'ATOM',JATOM, &
              '  AT RMT: SIGMA OF V-XC FIT: ',SIGMAJ(1)
      WRITE(6,'(a,i4,a,i4,")=",f6.3,a,2e13.5)') 'ATOM',JATOM, &
              '  LARGEST SIGMA AT R(',irsigmaup,r(irsigmaup),' and MAX DIFF:'&
              ,sigmamaxup,dmaxup

      if(spin)then
        IF(iatj.LT.0) WRITE(6,1540)
        DO k=1,npt1
          TEST=0.D0
!$omp simd reduction(+:TEST)
          DO LM1=1,LLMM
            TEST=TEST+bbd(LM1)*YY(LM1,k)
          ENDDO
          if((k.eq.1).or.((iprint.gt.0).and.(MOD(k,2).EQ.0))) then
            IF(iatj.LT.0) THEN
              WRITE(6,1550)k,vxcsd(k),TEST,vxcsd(k)-TEST
            ELSE
              WRITE(6,*) k,vxcsd(k),TEST
            ENDIF
          endif
        ENDDO
        WRITE(6,'(a,i4,a,10x,2e20.5)') 'ATOM',JATOM, &
                '  AT RMT: SIGMA OF V-XC FIT: ',SIGMAJ(2)
        WRITE(6,'(a,i4,a,i3,")=",f6.3,a,2e13.5)') 'ATOM',JATOM, &
                '  LARGEST SIGMA AT R(',irsigmadn,r(irsigmadn),' and MAX DIFF:'&
                ,sigmamaxdn,dmaxdn
      endif
!     For unpolarized case, copy values
      if(nospin)then
        sigmaj(2)=sigmaj(1)
        ir=jri(jatom)
        exc(1:ir,1:lmmax(jatom),JATOM,2)=exc(1:ir,1:lmmax(jatom),JATOM,1)
        vxc(1:ir,1:lmmax(jatom),JATOM,2)=vxc(1:ir,1:lmmax(jatom),JATOM,1)
        if(MGGA_TAU_SWITCH) vtau(:,:,:,2) = vtau(:,:,:,1)
        if(MGGA_LAPL_SWITCH) vlapl(:,:,:,2) = vlapl(:,:,:,1)
      endif

!      DO  LM1=1,LLMM
!      if(lm(1,lm1,jatom).eq.1) then
!              vxc(1:ir,lm1,JATOM,1:2)=0.d0
!      endif
!      if(lm(1,lm1,jatom).eq.2) then
!              vxc(1:ir,lm1,JATOM,1:2)=0.d0
!      endif
!      if(lm(1,lm1,jatom).eq.3) then
!              vxc(1:ir,lm1,JATOM,1:2)=0.d0
!      endif
!      enddo

 1530 FORMAT(3X,'CONDITION PARAMETER ATOM',I3,2X,A10,2X,'RCOND= ',       &
             F8.5)                                                      
 1540 FORMAT(3X,'INDEX',6X,'V-XC',7X,'FIT',5X,'DIFFERENCE')             
 1550 FORMAT(I8,2F10.5,F15.5)                                           


      deallocate(spt,weight,vxcsu,vxcsd,vxu,vxd,exu,exd)
      deallocate(vtau_rad_up,vtau_rad_dn,vlapl_rad_up,vlapl_rad_dn,dcsp)
      deallocate(ddcsp,dydt,dydt2,dydp,dydp2)
      deallocate(dydtp,qrhou,qrhod,qsin,qcos)
      deallocate(YL,dyt,dyt2)
      deallocate(dyp,dyp2,dytp)
      
      RETURN                                                            
      END
                                                               
