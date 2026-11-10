SUBROUTINE PWXAD5(I1,I2,I3,tvec,am, &
     fft,tauup,taudn,ifft1,ifft2,ifft3,vxu,vxd,exu,exd,vxuold,vxdold,&
     vxunonloc,vxdnonloc,vtau_up,vtau_dn,vlapl_up,vlapl_dn,fft_ext,&
     meta_fft,meta_fft_ext,JSPIN,&
     gga_int_up_11, gga_int_up_12, gga_int_up_13 , &
     gga_int_up_21, gga_int_up_22, gga_int_up_23 , &
     gga_int_up_31, gga_int_up_32, gga_int_up_33,  &
     gga_int_dn_11, gga_int_dn_12, gga_int_dn_13 , &
     gga_int_dn_21, gga_int_dn_22, gga_int_dn_23 , &
     gga_int_dn_31, gga_int_dn_32, gga_int_dn_33     )
! -------------------------------------------------------------------------
      USE xcparam, only: gga_switch,mgga_tau_switch,mgga_lapl_switch
      USE PWTest, only: pwindices, indmax
      use parallel, only: abort_parallel
      use str_switch, only : SW_STRESS
!      use gga_str_kb, only: sig_str_u,sig_str_d

!     IMPLICIT REAL*8(A-H,O-Z)
      implicit none

!     arguments
      integer, intent(in) :: i1,i2,i3,ifft1,ifft2,ifft3,jspin
      real*8, intent(in)  :: vxuold,vxdold, &
           vxunonloc,vxdnonloc
      real*8, intent(out) :: vxu,vxd,exu,exd, &
           vtau_up,vtau_dn,vlapl_up,vlapl_dn
      real*8, intent(inout) :: tauup,taudn
      real*8, intent(out) :: gga_int_up_11, gga_int_up_12, gga_int_up_13 , &
                             gga_int_up_21, gga_int_up_22, gga_int_up_23 , &
                             gga_int_up_31, gga_int_up_32, gga_int_up_33,  &
                             gga_int_dn_11, gga_int_dn_12, gga_int_dn_13 , &
                             gga_int_dn_21, gga_int_dn_22, gga_int_dn_23 , &
                             gga_int_dn_31, gga_int_dn_32, gga_int_dn_33
!     reciprocal vectors
      ! TVEC(1,I)=1.D0
      ! TVEC(2,I)=RVEC(1,I)
      ! TVEC(3,I)=RVEC(2,I)
      ! TVEC(4,I)=RVEC(3,I)
      ! TVEC(5,I)=RVEC(1,I)**2
      ! TVEC(6,I)=RVEC(2,I)**2
      ! TVEC(7,I)=RVEC(3,I)**2
      ! TVEC(8,I)=RVEC(2,I)*RVEC(3,I)
      ! TVEC(9,I)=RVEC(1,I)*RVEC(3,I)
      ! TVEC(10,I)=RVEC(1,I)*RVEC(2,I)

      real*8, intent(in)  :: tvec(10,50)
      
!     products of tvec summed over quadrature points
      real*8, intent(in)  :: am(125)

!     fft_arrays as input (for smoothing)
      complex*16, intent(in) ::  fft(ifft1,ifft2,ifft3,2)
      complex*16, intent(in) ::  fft_ext(ifft1,ifft2,2,2)

      complex*16, intent(in) ::  meta_fft(ifft1,ifft2,ifft3,2)
      complex*16, intent(in) ::  meta_fft_ext(ifft1,ifft2,2,2)

!     local variables and parameters
      real*8 :: BVT,BV1,BV2
      dimension  BVT(10),BV1(10),BV2(10)
      real*8 :: twobyrt3, tol, chargecut
      parameter (twobyrt3=1.1547005383792517D0, TOL=1.D-28)
      parameter (ChargeCut=1D-9)
      logical, parameter        :: test=.false.
!
!     Gerschberg-Stephanopolos (positivity) constraint
      logical, parameter        :: GSPW = .true.

!     derivatives of charge density
      real*8   :: du0,dd0,dup,ddn,gxu,gyu,gzu,gmagu,g2u,ggxu,ggyu,ggzu,fact, &
           gx,gy,gz,gmagt,ggx,ggy,ggz,gdggt,gxd,gyd,gzd,gmagd,g2d, &
           ggxd,ggyd,ggzd,gdggu,gdggd
      
!     derivatives of kinetic energy density
      real*8   ::  dtau_dx, dtau_dy, dtau_dz, &
           dtau_mag, d2tau, d2tau_x, d2tau_y, d2tau_z, &
           dtau_dx_u, dtau_dy_u, dtau_dz_u, dtau_mag_u, d2tau_u, &
           d2tau_dx2_u, d2tau_dy2_u, d2tau_dz2_u, &
           d2tau_dx2, d2tau_dy2, d2tau_dz2, &
           dtau_dx_d, dtau_dy_d, dtau_dz_d, dtau_mag_d, d2tau_d, d2tau_dx2_d, &
           d2tau_dy2_d, d2tau_dz2_d, dtau_drho, dtau_u_drho_u, &
           dtau_u_drho_d, dtau_d_drho_u, dtau_d_drho_d, &
           sig_str_u,sig_str_d,sig_str_ud
           
      integer  :: lx1,lx2,ind,ia,ja,ka
!
!   EXCHANGE-CORRELATION POTENTIAL SECTION
!
!  CALL CPUTIM(TIME0)
      DU0=fft(i1,i2,i3,1)
      DD0=fft(i1,i2,i3,2)
      
      if(GSPW)then
        DU0 = max(DU0,0.D0)
        DD0 = max(DD0,0.D0)
      endif

      IF(GGA_SWITCH .EQV. .true.) THEN
         
         DO LX1=1,10
            BV1(LX1)=0.D0
            BV2(LX1)=0.D0
         end DO
!
         DO IND=1,INDMAX
            ia=I1+PWIndices(1,IND)
            if(ia.lt.1    ) ia=ia+ifft1
            if(ia.gt.ifft1) ia=ia-ifft1
            
            ja=I2+PWIndices(2,IND)
            if(ja.lt.1    ) ja=ja+ifft2
            if(ja.gt.ifft2) ja=ja-ifft2
            
            ka=I3+PWIndices(3,IND)
            
            if(ka.lt.1) then
               DUP=fft_ext(ia,ja,1,1)
               DDN=fft_ext(ia,ja,1,2)
            else if(ka.gt.ifft3) then
               DUP=fft_ext(ia,ja,2,1)
               DDN=fft_ext(ia,ja,2,2)
            else
               DUP=fft(ia,ja,ka,1)
               DDN=fft(ia,ja,ka,2)
            endif
            
            if(GSPW)then
               DUP = max(DUP,0.D0)
               DDN = max(DDN,0.D0)
            endif

!       A Kahan summation could be used here, although there are probably
!       too few points for it to matter.
!
            DUP = DUP - DU0
            DDN = DDN - DD0
            
            DO LX2=1,10
               BV1(LX2)=BV1(LX2)+DUP*TVEC(LX2,IND)
               BV2(LX2)=BV2(LX2)+DDN*TVEC(LX2,IND)
            end DO
1819           format(3i3,11f10.7)

         ENDDO
!
         CALL CHSLV(AM,10,BV1)
         CALL CHSLV(AM,10,BV2) !BV = MATPROD(AM,BV)
         DO LX2=1,10
            BVT(LX2)=BV1(LX2)+BV2(LX2)
         end DO
         
         GXU=BV1(2)
         GYU=BV1(3)
         GZU=BV1(4)
         GMAGU=(GXU*GXU+GYU*GYU+GZU*GZU)
         G2U=2.D0*(BV1(5)+BV1(6)+BV1(7))
         if(abs(gmagu) .le. tol) then
            gmagu=0
            ggxu =twobyrt3*bv1(5)
            ggyu =twobyrt3*bv1(6)
            ggzu =twobyrt3*bv1(7)
            gdggu=gxu*ggxu+gyu*ggyu+gzu*ggzu
         else
            gmagu=sqrt(gmagu)
            FACT=1.D0/GMAGU
            GGXU=FACT*(2.D0*BV1(5)*GXU+BV1(10)*GYU+BV1(9)*GZU)
            GGYU=FACT*(BV1(10)*GXU+2.D0*BV1(6)*GYU+BV1(8)*GZU)
            GGZU=FACT*(BV1(9)*GXU+BV1(8)*GYU+2.D0*BV1(7)*GZU)
            gdggu=gxu*ggxu+gyu*ggyu+gzu*ggzu
         end if
         GX=BVT(2)
         GY=BVT(3)
         GZ=BVT(4)
         gmagt=(GX*GX+GY*GY+GZ*GZ)
         if(abs(gmagt) .le. tol) then
            gmagt=0.D0
            ggx =twobyrt3*bvt(5)
            ggy =twobyrt3*bvt(6)
            ggz =twobyrt3*bvt(7)
            gdggt=GX*GGX+GY*GGY+GZ*GGZ
         else
            gmagt=sqrt(gmagt)
            FACT=1.D0/gmagt
            GGX=FACT*(2.D0*BVT(5)*GX+BVT(10)*GY+BVT(9)*GZ)
            GGY=FACT*(BVT(10)*GX+2.D0*BVT(6)*GY+BVT(8)*GZ)
            GGZ=FACT*(BVT(9)*GX+BVT(8)*GY+2.D0*BVT(7)*GZ)
            GDGGT=GX*GGX+GY*GGY+GZ*GGZ
         end if
         GXD=BV2(2)
         GYD=BV2(3)
         GZD=BV2(4)
         GMAGD=(GXD*GXD+GYD*GYD+GZD*GZD)
         G2D=2.D0*(BV2(5)+BV2(6)+BV2(7))
         if(abs(gmagd) .le. tol) then
            gmagd=0.0D0
            ggxd =twobyrt3*bv2(5)
            ggyd =twobyrt3*bv2(6)
            ggzd =twobyrt3*bv2(7)
            gdggd=gxd*ggxd+gyd*ggyd+gzd*ggzd
         else
            gmagd=sqrt(gmagd)
            FACT=1.D0/GMAGD
            GGXD=FACT*(2.D0*BV2(5)*GXD+BV2(10)*GYD+BV2(9)*GZD)
            GGYD=FACT*(BV2(10)*GXD+2.D0*BV2(6)*GYD+BV2(8)*GZD)
            GGZD=FACT*(BV2(9)*GXD+BV2(8)*GYD+2.D0*BV2(7)*GZD)
            gdggd=gxd*ggxd+gyd*ggyd+gzd*ggzd
         endif

         if(mgga_tau_switch) then

            ! tauup=dble(meta_fft(i1,i2,i3,1))
            ! taudn=dble(meta_fft(i1,i2,i3,2))
            
            DO LX1=1,10
               BV1(LX1)=0.D0
               BV2(LX1)=0.D0
            end DO
!
            DO IND=1,INDMAX
               ia=I1+PWIndices(1,IND)
               if(ia.lt.1    ) ia=ia+ifft1
               if(ia.gt.ifft1) ia=ia-ifft1
            
               ja=I2+PWIndices(2,IND)
               if(ja.lt.1    ) ja=ja+ifft2
               if(ja.gt.ifft2) ja=ja-ifft2
               
               ka=I3+PWIndices(3,IND)
               
               if(ka.lt.1) then
                  DUP=meta_fft_ext(ia,ja,1,1)
                  DDN=meta_fft_ext(ia,ja,1,2)
               else if(ka.gt.ifft3) then
                  DUP=meta_fft_ext(ia,ja,2,1)
                  DDN=meta_fft_ext(ia,ja,2,2)
               else
                  DUP=meta_fft(ia,ja,ka,1)
                  DDN=meta_fft(ia,ja,ka,2)
               endif
            
               ! if(GSPW)then
               !    DUP = max(DUP,0.D0)
               !    DDN = max(DDN,0.D0)
               ! endif

               DUP = DUP - tauup
               DDN = DDN - taudn
            
               DO LX2=1,10
                  BV1(LX2)=BV1(LX2)+DUP*TVEC(LX2,IND)
                  BV2(LX2)=BV2(LX2)+DDN*TVEC(LX2,IND)
               end DO
               
            ENDDO
!
            CALL CHSLV(AM,10,BV1)
            CALL CHSLV(AM,10,BV2) !BV = MATPROD(AM,BV)
            DO LX2=1,10
               BVT(LX2)=BV1(LX2)+BV2(LX2)
            end DO
         
            dtau_dx_u=BV1(2)
            dtau_dy_u=BV1(3)
            dtau_dz_u=BV1(4)
            dtau_mag_u=(dtau_dx_u*dtau_dx_u+dtau_dy_u*dtau_dy_u+dtau_dz_u*dtau_dz_u)
            d2tau_u=2.D0*(BV1(5)+BV1(6)+BV1(7))
            if(abs(dtau_mag_u) .le. tol) then
               dtau_mag_u=0
               d2tau_dx2_u=twobyrt3*bv1(5)
               d2tau_dy2_u=twobyrt3*bv1(6)
               d2tau_dz2_u=twobyrt3*bv1(7)
            else
               dtau_mag_u=sqrt(dtau_mag_u)
               FACT=1.D0/dtau_mag_u
               d2tau_dx2_u=FACT*(2.D0*BV1(5)*dtau_dx_u+BV1(10)*dtau_dy_u+BV1(9)*dtau_dz_u)
               d2tau_dy2_u=FACT*(BV1(10)*dtau_dx_u+2.D0*BV1(6)*dtau_dy_u+BV1(8)*dtau_dz_u)
               d2tau_dz2_u=FACT*(BV1(9)*dtau_dx_u+BV1(8)*dtau_dy_u+2.D0*BV1(7)*dtau_dz_u)
            end if

            dtau_u_drho_u = dtau_dx_u*gxu + dtau_dy_u*gyu + dtau_dz_u*gzu
            dtau_u_drho_d = dtau_dx_u*gxd + dtau_dy_u*gyd + dtau_dz_u*gzd

            dtau_dx=BVT(2)
            dtau_dy=BVT(3)
            dtau_dz=BVT(4)
            dtau_mag=(dtau_dx*dtau_dx+dtau_dy*dtau_dy+dtau_dz*dtau_dz)
            if(abs(dtau_mag) .le. tol) then
               dtau_mag=0.D0
               d2tau_dx2=twobyrt3*bvt(5)
               d2tau_dy2=twobyrt3*bvt(6)
               d2tau_dz2=twobyrt3*bvt(7)
            else
               dtau_mag=sqrt(dtau_mag)
               FACT=1.D0/dtau_mag
               d2tau_dx2=FACT*(2.D0*BVT(5)*dtau_dx+BVT(10)*dtau_dy+BVT(9)*dtau_dz)
               d2tau_dy2=FACT*(BVT(10)*dtau_dx+2.D0*BVT(6)*dtau_dy+BVT(8)*dtau_dz)
               d2tau_dz2=FACT*(BVT(9)*dtau_dx+BVT(8)*dtau_dy+2.D0*BVT(7)*dtau_dz)
            end if
            dtau_drho = dtau_dx*gx + dtau_dy*gy + dtau_dz*gz

            dtau_dx_d=BV2(2)
            dtau_dy_d=BV2(3)
            dtau_dz_d=BV2(4)
            dtau_mag_d=(dtau_dx_d*dtau_dx_d+dtau_dy_d*dtau_dy_d+dtau_dz_d*dtau_dz_d)
            d2tau_d=2.D0*(BV2(5)+BV2(6)+BV2(7))
            if(abs(dtau_mag_d) .le. tol) then
               dtau_mag_d=0
               d2tau_dx2_d=twobyrt3*bv2(5)
               d2tau_dy2_d=twobyrt3*bv2(6)
               d2tau_dz2_d=twobyrt3*bv2(7)
            else
               dtau_mag_d=sqrt(dtau_mag_d)
               FACT=1.D0/dtau_mag_d
               d2tau_dx2_d=FACT*(2.D0*BV2(5)*dtau_dx_d+BV2(10)*dtau_dy_d+BV2(9)*dtau_dz_d)
               d2tau_dy2_d=FACT*(BV2(10)*dtau_dx_d+2.D0*BV2(6)*dtau_dy_d+BV2(8)*dtau_dz_d)
               d2tau_dz2_d=FACT*(BV2(9)*dtau_dx_d+BV2(8)*dtau_dy_d+2.D0*BV2(7)*dtau_dz_d)
            end if

            dtau_d_drho_u = dtau_dx_d*gxu + dtau_dy_d*gyu + dtau_dz_d*gzu
            dtau_d_drho_d = dtau_dx_d*gxd + dtau_dy_d*gyd + dtau_dz_d*gzd
            
         else if(mgga_lapl_switch) then
!          Compute derivatives of laplacian of density
            call outerr('LAPW0', 'Laplacian-dependent meta-gga`s&
                 & are not yet implemented')
            call abort_parallel
            stop "not yet implemented"
         end if
         
      ENDIF
!
8886  format (1(1x,i3),4(3x,f10.6,1x,f10.6))
      CALL VXCLM2(du0,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
           dd0,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD, &
           vxu,vxd,exu,exd,vtau_up,vtau_dn,vlapl_up,vlapl_dn, &
           dtau_drho,dtau_u_drho_u,dtau_u_drho_d,dtau_d_drho_u,dtau_d_drho_d,&
           vxuold,vxdold,vxunonloc,vxdnonloc,GGX,GGY,GGZ, &
           GDGGT,gdggu,gdggd,tauup,taudn, &
           i3+i2*1000+i1*1000000,chargeCut,JSPIN,0,sig_str_u,sig_str_d,sig_str_ud)
!if(i1.eq.5.and.i2.eq.5) write(6,*) 'pwxa ',i3,i3+i2*1000+i1*1000000,sig_str_u,vxu

!if(i2.eq.1.and.i3.eq.1) write(77,'(i5,f12.7,14f10.5)') i1,du0,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU,tauup,vtau_up,vxu,exu
!if(i2.eq.3.and.i3.eq.3) write(78,'(i5,f12.7,14f10.5)') i1,du0,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU,tauup,vtau_up,vxu,exu
!if(i2.eq.30.and.i3.eq.30) write(79,'(i5,f12.7,14f10.5)') i1,du0,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU,tauup,vtau_up,vxu,exu
!if(i2.eq.100.and.i3.eq.100) write(80,'(i5,f12.7,14f10.5)') i1,du0,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU,tauup,vtau_up,vxu,exu
!=================For GGA INT KAMAL===========
      IF(SW_STRESS .eq. 'YES') CALL GGA_INT_STRESS(GX,GY,GZ,GXU,GYU,GZU,GXD,GYD,GZD,&
                        JSPIN,du0,dd0,sig_str_u,sig_str_d,sig_str_ud, &
                        gga_int_up_11, gga_int_up_12, gga_int_up_13 , &
                        gga_int_up_21, gga_int_up_22, gga_int_up_23 , &
                        gga_int_up_31, gga_int_up_32, gga_int_up_33,  &
                        gga_int_dn_11, gga_int_dn_12, gga_int_dn_13 , &
                        gga_int_dn_21, gga_int_dn_22, gga_int_dn_23 , &
                        gga_int_dn_31, gga_int_dn_32, gga_int_dn_33 )
!=================end====================                
!if(i1.eq.5.and.i2.eq.5) write(6,*) 'pwxa5',i3,gga_int_up_11,sig_str_u,vxu
      RETURN
    END SUBROUTINE PWXAD5
    
SUBROUTINE GGA_INT_STRESS(GX,GY,GZ,GXU,GYU,GZU,GXD,GYD,GZD,JSPIN,du0,dd0,&
                        sig_str_u,sig_str_d,sig_str_ud, &
                        gga_int_up_11, gga_int_up_12, gga_int_up_13 , &
                        gga_int_up_21, gga_int_up_22, gga_int_up_23 , &
                        gga_int_up_31, gga_int_up_32, gga_int_up_33,  &
                        gga_int_dn_11, gga_int_dn_12, gga_int_dn_13 , &
                        gga_int_dn_21, gga_int_dn_22, gga_int_dn_23 , &
                        gga_int_dn_31, gga_int_dn_32, gga_int_dn_33 )
!This subroutine contains the interstitial part of the gga correction  --- KAMAL JUN 2021
!    use gga_str_kb, only: sig_str_u,sig_str_d
    IMPLICIT none
    real*8,intent(in)    :: GXU,GYU,GZU,GXD,GYD,GZD,GX,GY,GZ,du0,dd0,sig_str_u,sig_str_d,sig_str_ud
    integer,intent(in)   :: JSPIN
!     real*8, intent(out)  :: gga_inter(1:9,1:2)
    real*8, intent(out) :: gga_int_up_11, gga_int_up_12, gga_int_up_13 , &
                            gga_int_up_21, gga_int_up_22, gga_int_up_23 , &
                            gga_int_up_31, gga_int_up_32, gga_int_up_33,  &
                            gga_int_dn_11, gga_int_dn_12, gga_int_dn_13 , &
                            gga_int_dn_21, gga_int_dn_22, gga_int_dn_23 , &
                            gga_int_dn_31, gga_int_dn_32, gga_int_dn_33
    real*8               :: grad_rho(3,2)
    integer              :: alpha,beta,index1
    
    grad_rho(1,1) = GXU; grad_rho(2,1) = GYU; grad_rho(3,1) = GZU
    grad_rho(1,2) = GXD; grad_rho(2,2) = GYD; grad_rho(3,2) = GZD
    
!     index1 = 0
!     DO alpha = 1, 3
!        DO beta = 1, 3
!           index1 = index1 + 1
!             GGA_INTER(index1,1) = (grad_rho(alpha,1)+grad_rho(alpha,2))*&
!                                   (grad_rho(beta,1)+ grad_rho(beta,2))
!        ENDDO
!     ENDDO
!     
!     GGA_INTER(:,1) =  GGA_INTER(:,1)*sig_str_u/(du0+dd0)
!     IF(JSPIN .eq. 2) THEN
!       GGA_INTER(:,2) =  GGA_INTER(:,2)*sig_str_d
!     ELSE
!       GGA_INTER(:,2) = GGA_INTER(:,1) 
!     ENDIF      
    
    
    gga_int_up_11 = (grad_rho(1,1)+grad_rho(1,2))*&
                    (grad_rho(1,1)+ grad_rho(1,2))*sig_str_u/(du0+dd0)
    gga_int_up_12 = (grad_rho(1,1)+grad_rho(1,2))*&
                    (grad_rho(2,1)+ grad_rho(2,2))*sig_str_u/(du0+dd0)
    gga_int_up_13 = (grad_rho(1,1)+grad_rho(1,2))*&
                    (grad_rho(3,1)+ grad_rho(3,2))*sig_str_u/(du0+dd0)
    
    gga_int_up_21=  (grad_rho(2,1)+grad_rho(2,2))*&
                    (grad_rho(1,1)+ grad_rho(1,2))*sig_str_u/(du0+dd0)
    gga_int_up_22=  (grad_rho(2,1)+grad_rho(2,2))*&
                    (grad_rho(2,1)+ grad_rho(2,2))*sig_str_u/(du0+dd0)
    gga_int_up_23=  (grad_rho(2,1)+grad_rho(2,2))*&
                    (grad_rho(3,1)+ grad_rho(3,2))*sig_str_u/(du0+dd0) 
    
    gga_int_up_31 = (grad_rho(3,1)+grad_rho(3,2))*&
                    (grad_rho(1,1)+ grad_rho(1,2))*sig_str_u/(du0+dd0)
    gga_int_up_32 = (grad_rho(3,1)+grad_rho(3,2))*&
                    (grad_rho(2,1)+ grad_rho(2,2))*sig_str_u/(du0+dd0)
    gga_int_up_33 = (grad_rho(3,1)+grad_rho(3,2))*&
                    (grad_rho(3,1)+ grad_rho(3,2))*sig_str_u/(du0+dd0)
    
    IF( JSPIN .eq. 2) THEN
       gga_int_dn_11 = (grad_rho(1,1)+grad_rho(1,2))*&
                    (grad_rho(1,1)+ grad_rho(1,2))*sig_str_d/(du0+dd0)
       gga_int_dn_12 = (grad_rho(1,1)+grad_rho(1,2))*&
                    (grad_rho(2,1)+ grad_rho(2,2))*sig_str_d/(du0+dd0)
       gga_int_dn_13 = (grad_rho(1,1)+grad_rho(1,2))*&
                    (grad_rho(3,1)+ grad_rho(3,2))*sig_str_d/(du0+dd0)
    
       gga_int_dn_21=  (grad_rho(2,1)+grad_rho(2,2))*&
                    (grad_rho(1,1)+ grad_rho(1,2))*sig_str_d/(du0+dd0)
       gga_int_dn_22=  (grad_rho(2,1)+grad_rho(2,2))*&
                    (grad_rho(2,1)+ grad_rho(2,2))*sig_str_d/(du0+dd0)
       gga_int_dn_23=  (grad_rho(2,1)+grad_rho(2,2))*&
                    (grad_rho(3,1)+ grad_rho(3,2))*sig_str_d/(du0+dd0) 
    
       gga_int_dn_31 = (grad_rho(3,1)+grad_rho(3,2))*&
                    (grad_rho(1,1)+ grad_rho(1,2))*sig_str_d/(du0+dd0)
       gga_int_dn_32 = (grad_rho(3,1)+grad_rho(3,2))*&
                    (grad_rho(2,1)+ grad_rho(2,2))*sig_str_d/(du0+dd0)
       gga_int_dn_33 = (grad_rho(3,1)+grad_rho(3,2))*&
                    (grad_rho(3,1)+ grad_rho(3,2))*sig_str_d/(du0+dd0)          
    ELSE
        gga_int_dn_11  = gga_int_up_11
        gga_int_dn_12  = gga_int_up_12
        gga_int_dn_13  = gga_int_up_13
    
        gga_int_dn_21  = gga_int_up_21
        gga_int_dn_22  = gga_int_up_22
        gga_int_dn_23  = gga_int_up_23
    
        gga_int_dn_31 = gga_int_up_31
        gga_int_dn_32 = gga_int_up_32
        gga_int_dn_33 = gga_int_up_33
    ENDIF    
    
    
RETURN
END SUBROUTINE GGA_INT_STRESS
