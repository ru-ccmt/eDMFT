subroutine val_rela_kinstrmt(nemin,nemax,jatom,alm,blm,clm, &
                             buf_corr,buf_sr,buf_old,buf_surf)
   ! Scalar-relativistic correction to the valence kinetic stress inside
   ! one muffin-tin sphere.
   !
   ! WHY THIS ROUTINE IS NEEDED
   ! --------------------------
   ! The historical stress implementation evaluates the nonrelativistic
   ! operator p_alpha p_beta.  With CALC=RELA the radial functions and their
   ! normalization are scalar relativistic, but p_alpha p_beta is still the
   ! operator used by val_kinstrmt.  Merely including the small component in
   ! RINT13 does not turn that operator into the Koelling-Harmon operator.
   !
   ! In Rydberg units the upper component obeys the variational equation
   !
   !   H_SR(E) = -div[ A(r,E) grad ] + V(r),
   !   A(r,E)  = 1/M(r,E),
   !   M(r,E)  = 1 + (E - V(r))/(4*c_atomic**2).
   !
   ! A homogeneous strain changes the metric in the gradient bilinear.
   ! Consequently the scalar-relativistic kinetic stress is built from
   !
   !   Q_ab^SR = 1/4 Integral { (d_a psi)^* q_b
   !                           + q_a^* d_b psi
   !                           + (a <-> b) } d^3r,
   !
   ! where q=A grad(psi).  This symmetrized form is Hermitian also between
   ! radial basis functions generated at different linearization energies.
   ! For an energy derivative udot one must differentiate A as well:
   !
   !   q[udot] = A grad(udot) + A_E grad(u),
   !   A_E     = -1/(4*c_atomic**2*M**2).
   !
   ! The old RELA path corresponds, after partial integration, to the plain
   ! gradient bilinear of the large component plus 1/c_atomic**2 times the
   ! same bilinear of the stored small component.  We return the local,
   ! explicitly variational correction
   !
   !   buf_corr = Q_SR - Q_old,
   !
   ! as well as both terms separately.  This volume replacement is not the
   ! whole correction.  Partial integration produces an APW boundary term.
   ! The interstitial operator is nonrelativistic, whereas the radial flux on
   ! the muffin-tin side is A*d_r(psi).  The historical surface term uses
   ! d_r(psi) on both sides.  We therefore also return
   !
   !   buf_surf = R_MT^2 int dOmega n_a n_b psi^*
   !                              [A*d_r(psi)-d_r(psi)] .
   !
   ! This is a flux correction, not the fixed-RMT shape derivative.  The two
   ! have different origins and must never be combined.  l2main supplies
   ! occupation factors, spin factor two, equivalent-atom rotations, and the
   ! final 1/volume.
   !
   ! ANGULAR ALGEBRA
   ! ---------------
   ! Appendix C of Belbase, Troester, and Blaha gives
   !
   ! d_a[f_l(r)Y_lm] = rhat_a f'_l Y_lm
   !                  + f_l/r Sum_{s=+-1,t} c^a_st(l,m)Y_{l+s,m+t}.
   !
   ! hcoef below expands rhat_a Y_lm in the same output harmonics, while
   ! dcoef stores c^a_st.  Once both gradients are expressed in an
   ! orthonormal Y_LM basis, the angular integral is just a dot product.

   use param
   use com, only: rel
   use structure, only: r0,dx,jri
   use atspdt, only: el
   use xa, only: band_energy=>e
   use dmft_stress, only: stress_rela_setup_wall,stress_rela_band_wall
   use w_atpar, only: w_vr,w_rrad1,w_rrad2,w_rade1,w_rade2, &
                       w_radp1,w_radp2,w_radpe1,w_radpe2, &
                       w_a1lo,w_b1lo,w_ap1lo,w_bp1lo
   use lo, only: ilo_src=>ilo,loor_src=>loor,rlo_src=>rlo,elo_src=>elo
   implicit none

   integer, intent(in) :: nemin,nemax,jatom
   complex*16, intent(in) :: alm((lmax2+1)*(lmax2+1),nume), &
                             blm((lmax2+1)*(lmax2+1),nume), &
                             clm((lmax2+1)*(lmax2+1),nume,nloat)
   complex*16, intent(out) :: buf_corr(3,3,nume),buf_sr(3,3,nume), &
                              buf_old(3,3,nume),buf_surf(3,3,nume)

   integer :: lmax,nin,nout,nout2,imax,l,m,ll,mm,lm,lmout,v,t,alpha,beta
   integer :: ir,num,jlo,ich,nchan,diag_status,imap,small_scale_status, &
              small_surf_status,ensmall_status,ennat_status,massv_status
   real*8 :: c_atomic,c2inv,me,mass,ainv,ae,aee,r,energy,g0
   real*8 :: gaunt_value,h_norm,d_norm,hd_cross,max_h_err,max_d_err, &
             max_flux_err,max_flux_boundary_err,flux_scale
   real*8 :: sum_re,sum_im,old_small_scale,small_surf_scale, &
             ensmall_scale,small_norm,ennat_scale,massv_scale
   real*8 :: setup_wall_begin,band_wall_begin,band_wall_end
   real*8, external :: gaunt
   character(len=32) :: diag_env,small_scale_env,small_surf_env,ensmall_env, &
                        ennat_env,massv_env
   complex*16 :: c_am(3,-1:1),cfac(3,-1:1,-1:1),zval
   integer, allocatable :: nmap(:),outmap(:,:),nmap2(:),outmap2(:,:)
   complex*16, allocatable :: hmap(:,:,:),dmap(:,:,:),hmap2(:,:,:)
   complex*16, allocatable :: fvec(:),dfvec(:),svec(:),dsvec(:)
   complex*16, allocatable :: qdvec(:),quvec(:)
   complex*16, allocatable :: grad(:,:),aux(:,:),sgrad(:,:), &
                              npsi(:,:),ndflux(:,:),nsval(:,:),naux(:,:,:)
   complex*16, allocatable :: integ_sr(:,:,:),integ_old(:,:,:), &
                              integ_nat(:,:,:)
   real*8, allocatable :: rmesh(:),funct_re(:),funct_im(:), &
                          vry(:),dvry(:),massv_tensor(:,:,:)
   real*8, allocatable :: lf(:,:,:),ldf(:,:,:),sf(:,:,:),sdf(:,:,:)
   real*8, allocatable :: qd(:,:,:),qu(:,:,:),tmp(:)
   ! dmft2 caches every atom because its outer loop is over k points,
   ! whereas LAPW2 keeps only the current atom in module storage.  Copy the
   ! requested atom into local arrays so that the tested LAPW2 kernel below
   ! remains algebraically identical.
   real*8 :: vr(nrad)
   real*8 :: rrad1(nrad,0:lmax2),rrad2(nrad,0:lmax2)
   real*8 :: rade1(nrad,0:lmax2),rade2(nrad,0:lmax2)
   real*8 :: radp1(nrad,0:lmax2),radp2(nrad,0:lmax2)
   real*8 :: radpe1(nrad,0:lmax2),radpe2(nrad,0:lmax2)
   real*8 :: a1lo(nrad,nloat,0:lomax),b1lo(nrad,nloat,0:lomax)
   real*8 :: ap1lo(nrad,nloat,0:lomax),bp1lo(nrad,nloat,0:lomax)
   integer :: ilo(0:lmax2)
   logical :: loor(nloat,0:lmax2),rlo(nloat,0:lmax2)
   real*8 :: elo(0:lmax2,nloat)
   logical :: sdlo(nloat,0:lomax)

   buf_corr = (0.d0,0.d0)
   buf_sr   = (0.d0,0.d0)
   buf_old  = (0.d0,0.d0)
   buf_surf = (0.d0,0.d0)
   if(.not.rel) return

   call walltim(setup_wall_begin)

   vr=w_vr(:,jatom)
   rrad1=w_rrad1(:,:,jatom); rrad2=w_rrad2(:,:,jatom)
   rade1=w_rade1(:,:,jatom); rade2=w_rade2(:,:,jatom)
   radp1=w_radp1(:,:,jatom); radp2=w_radp2(:,:,jatom)
   radpe1=w_radpe1(:,:,jatom); radpe2=w_radpe2(:,:,jatom)
   a1lo=w_a1lo(:,:,:,jatom); b1lo=w_b1lo(:,:,:,jatom)
   ap1lo=w_ap1lo(:,:,:,jatom); bp1lo=w_bp1lo(:,:,:,jatom)
   ! dmft2 stores LO metadata only through LOMAX, while the WIEN2k kernel
   ! loops over all radial l channels through LMAX2.  Pad the absent higher-l
   ! channels explicitly; otherwise l>LOMAX reads beyond LOOR/ELO even when
   ! there is no local orbital in that channel.
   ilo=0; loor=.false.; rlo=.false.; elo=0.d0
   ilo(0:lomax)=ilo_src(0:lomax)
   loor(:,0:lomax)=loor_src(:,0:lomax)
   rlo(:,0:lomax)=rlo_src(:,0:lomax)
   elo(0:lomax,:)=elo_src(0:lomax,:)
   ! dftKS/dmft2 does not generate WIEN2k's second-derivative SDLO basis.
   ! Every cached local orbital is therefore an ordinary LO in this port.
   sdlo=.false.

   lmax = lmax2n
   nin  = (lmax+1)*(lmax+1)
   ! A Cartesian derivative of Y_lm contains l+1, so the output angular
   ! space needs one shell more than the input LAPW expansion.
   nout = (lmax+2)*(lmax+2)
   nout2 = (lmax+3)*(lmax+3)
   imax = jri(jatom)
   nchan = nloat+2             ! 1=u, 2=udot, 2+jlo=local orbital

   ! Each input Y_lm reaches at most six output harmonics: l+-1 and
   ! m+t, t=-1,0,1.  Store this sparse map explicitly.  A dense
   ! (nout x nin) contraction made the first diagnostic prototype roughly
   ! twenty-five times more expensive without changing any algebra.
   allocate(nmap(nin),outmap(6,nin))
   allocate(hmap(3,6,nin),dmap(3,6,nin))
   allocate(nmap2(nout),outmap2(6,nout),hmap2(3,6,nout))
   allocate(rmesh(imax))
   allocate(vry(imax),dvry(imax))
   allocate(lf(imax,0:lmax,nchan),ldf(imax,0:lmax,nchan))
   allocate(sf(imax,0:lmax,nchan),sdf(imax,0:lmax,nchan))
   allocate(qd(imax,0:lmax,nchan),qu(imax,0:lmax,nchan))

   nmap = 0
   outmap = 0
   hmap = (0.d0,0.d0)
   dmap = (0.d0,0.d0)
   nmap2 = 0
   outmap2 = 0
   hmap2 = (0.d0,0.d0)
   call c_alpha_m(c_am)

   ! Precompute the two angular maps in Eq. (C9).  GAUNT(L,l,1,M,m,t)
   ! is Integral Y_LM^* Y_lm Y_1t dOmega in WIEN2k's convention.
   do l=0,lmax
      do m=-l,l
         lm = l*(l+1)+m+1
         call c_factor(l,m,cfac)
         do v=-1,1,2
            ll=l+v
            if(ll.lt.0) cycle
            do t=-1,1
               mm=m+t
               if(abs(mm).gt.ll) cycle
               lmout=ll*(ll+1)+mm+1
               nmap(lm)=nmap(lm)+1
               imap=nmap(lm)
               outmap(imap,lm)=lmout
               gaunt_value=gaunt(ll,l,1,mm,m,t)
               do alpha=1,3
                  dmap(alpha,imap,lm)=cfac(alpha,v,t)
                  hmap(alpha,imap,lm)=c_am(alpha,t)*gaunt_value
               enddo
            enddo
         enddo
      enddo
   enddo

   ! A second multiplication map is needed for the full tensor moment
   ! int n_alpha*n_beta*rho_KH.  AUX already spans angular momenta through
   ! lmax+1, so multiplying it by one more direction cosine reaches
   ! lmax+2.  Keeping this map separate prevents accidental truncation of
   ! the relativistic density moment at the top angular channel.
   do l=0,lmax+1
      do m=-l,l
         lm=l*(l+1)+m+1
         do v=-1,1,2
            ll=l+v
            if(ll.lt.0) cycle
            do t=-1,1
               mm=m+t
               if(abs(mm).gt.ll) cycle
               lmout=ll*(ll+1)+mm+1
               nmap2(lm)=nmap2(lm)+1
               imap=nmap2(lm)
               outmap2(imap,lm)=lmout
               gaunt_value=gaunt(ll,l,1,mm,m,t)
               do alpha=1,3
                  hmap2(alpha,imap,lm)=c_am(alpha,t)*gaunt_value
               enddo
            enddo
         enddo
      enddo
   enddo

   ! Optional exact angular identities.  These are valuable guards against
   ! a phase or index-convention error in the new relativistic term:
   ! Sum_a ||rhat_a Y_lm||^2=1, Sum_a ||r d_aY_lm||^2=l(l+1), and their
   ! summed cross term vanishes.
   max_h_err=0.d0
   max_d_err=0.d0
   do lm=1,nin
      l=int(sqrt(dble(lm-1)))
      h_norm=0.d0
      d_norm=0.d0
      hd_cross=0.d0
      do alpha=1,3
         do imap=1,nmap(lm)
            h_norm=h_norm+abs(hmap(alpha,imap,lm))**2
            d_norm=d_norm+abs(dmap(alpha,imap,lm))**2
            hd_cross=hd_cross+real(conjg(hmap(alpha,imap,lm))* &
                                   dmap(alpha,imap,lm))
         enddo
      enddo
      max_h_err=max(max_h_err,abs(h_norm-1.d0))
      max_d_err=max(max_d_err,abs(d_norm-dble(l*(l+1))))
      max_d_err=max(max_d_err,abs(hd_cross))
   enddo
   diag_env=''
   call get_environment_variable('WIEN_STRESS_RELA_DIAG',diag_env, &
                                 status=diag_status)
   if(diag_status.eq.0 .and. len_trim(diag_env).gt.0 .and. &
      trim(diag_env).ne.'0') then
      write(6,'(A,2(1X,ES23.15))') ':DBG_RELA_ANGULAR:', &
                                    max_h_err,max_d_err
   endif

   do ir=1,imax
      rmesh(ir)=r0(jatom)*exp(dx(jatom)*dble(ir-1))
   enddo

   lf=0.d0; ldf=0.d0; sf=0.d0; sdf=0.d0; qd=0.d0; qu=0.d0
   c_atomic=137.0359895d0
   c2inv=1.d0/(c_atomic*c_atomic)
   me=0.25d0*c2inv

   ! Diagnostic scale of the historical small-component p_a p_b term.
   ! Production always uses one.  Setting the environment variable to zero
   ! is useful for determining whether a discrepancy belongs to the
   ! Koelling-Harmon large-component metric or to WIEN2k's separately stored
   ! radial small component.  It is deliberately a real scale (rather than a
   ! Boolean) so that linearity can be checked.  No result obtained with a
   ! non-unit value is suitable for production.
   old_small_scale=1.d0
   small_scale_env=''
   call get_environment_variable('WIEN_STRESS_RELA_OLD_SMALL_SCALE', &
                                 small_scale_env,status=small_scale_status)
   if(small_scale_status.eq.0 .and. len_trim(small_scale_env).gt.0) &
      read(small_scale_env,*,err=997) old_small_scale

   ! Diagnostic/validation switch for the surface left by the historical
   ! small-component second-derivative operator.  VAL_KINSTRMT applies
   ! -d_a d_b to both components through RINT13, whereas
   ! STRESS_APW_SURFACE contains only the large component that is matched to
   ! the interstitial plane wave.  Converting the old small-component term
   ! to the gradient form subtracted in Q_OLD therefore requires the exact
   ! integration-by-parts surface
   !
   !   S_ab^small = R_MT^2/(2 c^2) int dOmega
   !       [ n_a s^* d_b s + n_b s^* d_a s ].
   !
   ! Keep it opt-in until its sign and magnitude have been checked against
   ! the independent finite-strain pressure.  A value of one applies the
   ! derived term; zero reproduces the first RELA prototype.
   small_surf_scale=0.d0
   small_surf_env=''
   call get_environment_variable('WIEN_STRESS_RELA_SMALL_SURF_SCALE', &
                                 small_surf_env,status=small_surf_status)
   if(small_surf_status.eq.0 .and. len_trim(small_surf_env).gt.0) &
      read(small_surf_env,*,err=996) small_surf_scale

   ! Diagnostic for the leading generalized-overlap term of the
   ! Koelling-Harmon implementation.  WIEN2k normalizes each radial pair as
   ! int(P^2+Q^2/c^2)dr, so differentiating the implemented generalized
   ! eigenproblem produces an energy times small-component overlap term that
   ! is absent from the upper-component metric replacement above.  For an
   ! isotropic strain its tensor share is -epsilon*N_small/3 before the
   ! common factor two applied in L2MAIN.  The exact coefficient is one;
   ! SCALE is exposed only to measure 0 versus 1 while the complete tensor
   ! derivation is checked.
   ensmall_scale=0.d0
   ensmall_env=''
   call get_environment_variable('WIEN_STRESS_RELA_ENSMALL_SCALE', &
                                 ensmall_env,status=ensmall_status)
   if(ensmall_status.eq.0 .and. len_trim(ensmall_env).gt.0) &
      read(ensmall_env,*,err=995) ensmall_scale

   ! Diagnostic energy-dependent overlap of the eliminated Koelling-Harmon
   ! equation.  For
   !
   !     H(E) = -div[A(E,V) grad] + V,
   !
   ! the stationary norm is 1-H_E.  Its relativistic tensor is therefore
   !
   !     N^E_ab = 1/(4 c^2) Re int (A d_a psi)^* (A d_b psi) d^3r.
   !
   ! This is *not* the radial B^2/c^2 norm tested above: it contains the
   ! complete angular gradient and can consequently be much larger for
   ! semicore p and d channels.  It would enter if LAPW2 differentiated the
   ! eliminated, energy-dependent upper-component equation directly.
   ! WIEN2k does not implement that generalized eigenproblem: LAPW1 builds
   ! H*c=epsilon*S*c with the stored large-plus-small RINT13 overlap, and the
   ! existing Pulay terms already differentiate that S.  Adding N^E_ab here
   ! therefore differentiates a second overlap which is not present in the
   ! implemented matrix problem.  The FeSe RELA-minus-NREL fixed-G test
   ! confirms the double count: unit scale adds about 0.825 GPa to the
   ! relativistic pressure error.  Production scale is zero.  The switch is
   ! retained only for reproducing and auditing the eliminated-equation
   ! diagnostic; a nonzero value is not a production setting.
   ennat_scale=0.d0
   ennat_env=''
   call get_environment_variable('WIEN_STRESS_RELA_ENNAT_SCALE', &
                                 ennat_env,status=ennat_status)
   if(ennat_status.eq.0 .and. len_trim(ennat_env).gt.0) &
      read(ennat_env,*,err=994) ennat_scale

   ! Trace diagnostic for the potential part generated by the same
   ! energy-dependent mass.  The eigenvalue derivative of H(E,V) couples
   ! dV to the natural Koelling-Harmon density
   !
   !   rho_KH = |psi|^2 + (4 c^2)^-1 |A grad psi|^2,
   !
   ! whereas WIEN2k's electrostatic cancellation uses the implemented
   ! two-component density |psi|^2+|B|^2/c^2.  Their norm difference is
   ! small, but under an affine dilation it is multiplied by r*dV/dr,
   ! which is large in the nuclear region.  The residual contribution to
   ! the total-energy derivative must be represented by the opposite sign
   ! in STR_VAL because mixer forms STR_DENS-STR_VAL-STR_CORE.
   ! The one-half below is not an adjustable coefficient.  L2MAIN applies
   ! a common factor two to buffers representing the two metric derivatives
   ! of p_alpha*p_beta.  Delta-rho*dV has only one strain derivative, so its
   ! buffer must compensate that common factor by 1/2.
   ! The FeSe fixed-G finite-strain test reduces the missing RELA pressure
   ! from about 26 GPa to below 0.4 GPa when this derived term and N^E_ab are
   ! combined.  Keep unit scale as the production value; the environment
   ! variable is a diagnostic multiplier, not a fitted coefficient.
   massv_scale=1.d0
   massv_env=''
   call get_environment_variable('WIEN_STRESS_RELA_MASSV_SCALE', &
                                 massv_env,status=massv_status)
   if(massv_status.eq.0 .and. len_trim(massv_env).gt.0) &
      read(massv_env,*,err=993) massv_scale

   ! VRY is the spherical potential in Rydberg.  A centered derivative on
   ! the logarithmic mesh is more reliable here than imposing an even/odd
   ! origin continuation on the Coulomb singularity.
   do ir=1,imax
      vry(ir)=2.d0*vr(ir)/rmesh(ir)
   enddo
   dvry(1)=(vry(2)-vry(1))/(rmesh(2)-rmesh(1))
   do ir=2,imax-1
      dvry(ir)=(vry(ir+1)-vry(ir-1))/(rmesh(ir+1)-rmesh(ir-1))
   enddo
   dvry(imax)=(vry(imax)-vry(imax-1))/(rmesh(imax)-rmesh(imax-1))

   ! Construct values and ordinary radial derivatives of every radial basis
   ! function.  LAPW stores r*u(r) and r*du(r)/dr, hence the division by
   ! RMESH.  It is essential to use ATPAR's RADP/RADPE/AP1LO arrays here.
   ! Re-differentiating r*u/r with DERGL looks innocuous but loses accuracy
   ! precisely at the small radii where the relativistic mass differs most
   ! from one.  More importantly, it violates OUTWIN's exact first-order
   ! Koelling-Harmon flux identity and produced a large, spurious residual in
   ! the FeSe finite-difference test.
   do l=0,lmax
      lf(:,l,1)=rrad1(1:imax,l)/rmesh
      lf(:,l,2)=rade1(1:imax,l)/rmesh
      sf(:,l,1)=rrad2(1:imax,l)/rmesh
      sf(:,l,2)=rade2(1:imax,l)/rmesh
      ldf(:,l,1)=radp1(1:imax,l)/rmesh
      ldf(:,l,2)=radpe1(1:imax,l)/rmesh
      sdf(:,l,1)=radp2(1:imax,l)/rmesh
      sdf(:,l,2)=radpe2(1:imax,l)/rmesh
      do jlo=1,ilo(l)
         if(.not.loor(jlo,l)) cycle
         ich=2+jlo
         lf(:,l,ich)=a1lo(1:imax,jlo,l)/rmesh
         sf(:,l,ich)=b1lo(1:imax,jlo,l)/rmesh
         ldf(:,l,ich)=ap1lo(1:imax,jlo,l)/rmesh
         sdf(:,l,ich)=bp1lo(1:imax,jlo,l)/rmesh
      enddo
   enddo

   ! q=A grad(u) for u, and q=d[A grad(u)]/dE for udot.  The potential
   ! stored in VR is r*V_Hartree; 2*VR/r is therefore V in Rydberg units,
   ! consistent with EL and ELO.
   !
   ! Do NOT reconstruct the radial component as A*u' from RADP below.
   ! OUTWIN integrates the coupled Koelling-Harmon equations and stores the
   ! second radial component B explicitly.  Its first equation is the exact
   ! mesh identity
   !
   !              A(r,E) u'(r,E) = 2 B(r,E)/r .
   !
   ! The same identity can be differentiated with respect to the Rydberg
   ! energy once or twice.  Thus 2*SF is also the complete radial flux for
   ! udot and for an SDLO; it already contains every A_E and A_EE term.
   ! Reconstructing it from a separately stored/numerically differentiated
   ! large-component slope produced 5--7 percent pointwise errors near the
   ! nucleus in RELA mode, even though the two forms coincided accidentally
   ! at R_MT.  The angular component has no analogous stored flux and is
   ! therefore still evaluated explicitly with A, A_E, and A_EE.
   do l=0,lmax
      do ir=1,imax
         r=rmesh(ir)
         mass=1.d0+(el(l)-2.d0*vr(ir)/r)*me
         ainv=1.d0/mass
         ae=-me/(mass*mass)
         qd(ir,l,1)=2.d0*sf(ir,l,1)
         qu(ir,l,1)=ainv*lf(ir,l,1)/r
         qd(ir,l,2)=2.d0*sf(ir,l,2)
         qu(ir,l,2)=ainv*lf(ir,l,2)/r+ae*lf(ir,l,1)/r
      enddo

      do jlo=1,ilo(l)
         if(.not.loor(jlo,l)) cycle
         ich=2+jlo
         energy=elo(l,jlo)
         do ir=1,imax
            r=rmesh(ir)
            mass=1.d0+(energy-2.d0*vr(ir)/r)*me
            ainv=1.d0/mass
            ae=-me/(mass*mass)
            aee=2.d0*me*me/(mass*mass*mass)

            if(sdlo(jlo,l) .and. abs(energy-el(l)).lt.1.d-6) then
               ! a1lo is d^2u/dE_Ry^2 for an SDLO.  Differentiate A twice;
               ! omitting these two terms would lose the relativistic
               ! response of the energy-linearized local orbital.
               qd(ir,l,ich)=2.d0*sf(ir,l,ich)
               qu(ir,l,ich)=ainv*lf(ir,l,ich)/r+ &
                    2.d0*ae*lf(ir,l,2)/r+aee*lf(ir,l,1)/r
            else
               ! Ordinary local orbital (or the conservative fallback for
               ! an uncommon SDLO whose parent energy differs from EL).
               qd(ir,l,ich)=2.d0*sf(ir,l,ich)
               qu(ir,l,ich)=ainv*lf(ir,l,ich)/r
            endif
         enddo
         if(rlo(jlo,l) .and. diag_status.eq.0 .and. &
            len_trim(diag_env).gt.0 .and. trim(diag_env).ne.'0') then
            write(6,'(A,3I5)') ':DBG_RELA_RLO_FALLBACK:',jatom,l,jlo
         endif
      enddo
   enddo

   ! OUTWIN's first Koelling-Harmon equation gives the independent identity
   ! A*u'(r)=2*B_stored(r)/r for a primitive radial solution (WIEN2k stores
   ! c times the physical small component).  Check it before using A in a
   ! stress expression; a Hartree/Rydberg factor error would show up here
   ! immediately and would contaminate the full correction.
   max_flux_err=0.d0
   max_flux_boundary_err=0.d0
   do l=0,lmax
      do ir=1,imax
         ! SF=B_stored/r, so the right-hand side is simply 2*SF.
         flux_scale=max(1.d-14,abs(qd(ir,l,1)),abs(2.d0*sf(ir,l,1)))
         max_flux_err=max(max_flux_err, &
              abs(qd(ir,l,1)-2.d0*sf(ir,l,1))/flux_scale)
      enddo
      flux_scale=max(1.d-14,abs(qd(imax,l,1)), &
                     abs(2.d0*sf(imax,l,1)))
      max_flux_boundary_err=max(max_flux_boundary_err, &
           abs(qd(imax,l,1)-2.d0*sf(imax,l,1))/flux_scale)
   enddo
   if(diag_status.eq.0 .and. len_trim(diag_env).gt.0 .and. &
      trim(diag_env).ne.'0') then
      write(6,'(A,I5,2(1X,ES23.15))') ':DBG_RELA_FLUX:',jatom, &
                                      max_flux_err,max_flux_boundary_err
   endif

   ! Band loop.  For each radial point, form the large-component gradient,
   ! its scalar-relativistic flux q, and the old small-component gradient in
   ! the common Y_LM output basis.  Orthogonality then performs the complete
   ! angular integral without an angular mesh.
   ! Every band writes a distinct NUM slice of BUF_SR/BUF_OLD/BUF_CORR and
   ! BUF_SURF.  It is therefore an ideal OpenMP dimension, unlike the angular
   ! channels of the nonrelativistic G1 contraction, whose work is strongly
   ! imbalanced.  All arrays below are scratch storage for one band and must
   ! be private.  Allocating them once per thread (outside OMP DO) avoids both
   ! races and the much larger cost of allocating them at every radial point.
   ! The shared radial functions/maps above are read-only in this region.
   call walltim(band_wall_begin)
   stress_rela_setup_wall=stress_rela_setup_wall+ &
        band_wall_begin-setup_wall_begin
!$omp parallel default(shared) &
!$omp& private(num,ir,l,m,lm,lmout,imap,jlo,ich,alpha,beta,v,zval, &
!$omp&         sum_re,sum_im,small_norm,fvec,dfvec,svec,dsvec,qdvec,quvec, &
!$omp&         grad,aux,sgrad,npsi,ndflux,nsval,naux,integ_sr,integ_old, &
!$omp&         integ_nat,funct_re,funct_im,massv_tensor,tmp)
   allocate(fvec(nin),dfvec(nin),svec(nin),dsvec(nin))
   allocate(qdvec(nin),quvec(nin))
   allocate(grad(3,nout),aux(3,nout),sgrad(3,nout))
   allocate(npsi(3,nout),ndflux(3,nout))
   allocate(nsval(3,nout),naux(3,3,nout2))
   allocate(integ_sr(3,3,imax),integ_old(3,3,imax),integ_nat(3,3,imax))
   allocate(funct_re(imax),funct_im(imax),massv_tensor(3,3,imax),tmp(imax))
!$omp do schedule(dynamic,1)
   do num=nemin,nemax
      integ_sr=(0.d0,0.d0)
      integ_old=(0.d0,0.d0)
      if(ennat_scale.ne.0.d0) integ_nat=(0.d0,0.d0)
      if(ensmall_scale.ne.0.d0) tmp=0.d0
      if(massv_scale.ne.0.d0) massv_tensor=0.d0
      do ir=1,imax
         fvec=(0.d0,0.d0); dfvec=(0.d0,0.d0)
         svec=(0.d0,0.d0); dsvec=(0.d0,0.d0)
         qdvec=(0.d0,0.d0); quvec=(0.d0,0.d0)
         do l=0,lmax
            do m=-l,l
               lm=l*(l+1)+m+1
               fvec(lm)=alm(lm,num)*lf(ir,l,1)+ &
                        blm(lm,num)*lf(ir,l,2)
               dfvec(lm)=alm(lm,num)*ldf(ir,l,1)+ &
                         blm(lm,num)*ldf(ir,l,2)
               svec(lm)=alm(lm,num)*sf(ir,l,1)+ &
                        blm(lm,num)*sf(ir,l,2)
               dsvec(lm)=alm(lm,num)*sdf(ir,l,1)+ &
                         blm(lm,num)*sdf(ir,l,2)
               qdvec(lm)=alm(lm,num)*qd(ir,l,1)+ &
                         blm(lm,num)*qd(ir,l,2)
               quvec(lm)=alm(lm,num)*qu(ir,l,1)+ &
                         blm(lm,num)*qu(ir,l,2)
               do jlo=1,ilo(l)
                  if(.not.loor(jlo,l)) cycle
                  ich=2+jlo
                  fvec(lm)=fvec(lm)+clm(lm,num,jlo)*lf(ir,l,ich)
                  dfvec(lm)=dfvec(lm)+clm(lm,num,jlo)*ldf(ir,l,ich)
                  svec(lm)=svec(lm)+clm(lm,num,jlo)*sf(ir,l,ich)
                  dsvec(lm)=dsvec(lm)+clm(lm,num,jlo)*sdf(ir,l,ich)
                  qdvec(lm)=qdvec(lm)+clm(lm,num,jlo)*qd(ir,l,ich)
                  quvec(lm)=quvec(lm)+clm(lm,num,jlo)*qu(ir,l,ich)
               enddo
            enddo
         enddo

         grad=(0.d0,0.d0); aux=(0.d0,0.d0); sgrad=(0.d0,0.d0)
         do alpha=1,3
            do lm=1,nin
               do imap=1,nmap(lm)
                  lmout=outmap(imap,lm)
                  grad(alpha,lmout)=grad(alpha,lmout)+ &
                       hmap(alpha,imap,lm)*dfvec(lm)+ &
                       dmap(alpha,imap,lm)*fvec(lm)/rmesh(ir)
                  aux(alpha,lmout)=aux(alpha,lmout)+ &
                       hmap(alpha,imap,lm)*qdvec(lm)+ &
                       dmap(alpha,imap,lm)*quvec(lm)
                  sgrad(alpha,lmout)=sgrad(alpha,lmout)+ &
                       hmap(alpha,imap,lm)*dsvec(lm)+ &
                       dmap(alpha,imap,lm)*svec(lm)/rmesh(ir)
               enddo
            enddo
         enddo

         do alpha=1,3
            do beta=1,3
               zval=(0.d0,0.d0)
               do lmout=1,nout
                  ! Q_SR: Hermitian and strain-symmetrized large/flux
                  ! bilinear.  Q_old: the plain large gradient plus the
                  ! 1/c^2-weighted stored small-component gradient.
                  zval=zval+0.25d0*( &
                       conjg(grad(alpha,lmout))*aux(beta,lmout)+ &
                       conjg(aux(alpha,lmout))*grad(beta,lmout)+ &
                       conjg(grad(beta,lmout))*aux(alpha,lmout)+ &
                       conjg(aux(beta,lmout))*grad(alpha,lmout))
               enddo
               integ_sr(alpha,beta,ir)=rmesh(ir)*rmesh(ir)*zval

               zval=(0.d0,0.d0)
               do lmout=1,nout
                  zval=zval+0.5d0*( &
                       conjg(grad(alpha,lmout))*grad(beta,lmout)+ &
                       conjg(grad(beta,lmout))*grad(alpha,lmout)+ &
                       old_small_scale*c2inv* &
                       conjg(sgrad(alpha,lmout))*sgrad(beta,lmout)+ &
                       old_small_scale*c2inv* &
                       conjg(sgrad(beta,lmout))*sgrad(alpha,lmout))
               enddo
               integ_old(alpha,beta,ir)=rmesh(ir)*rmesh(ir)*zval

               ! The eliminated-equation overlap is a diagnostic which is
               ! disabled in production.  Avoid its full angular contraction
               ! unless the corresponding audit scale was explicitly set.
               if(ennat_scale.ne.0.d0) then
                  zval=(0.d0,0.d0)
                  do lmout=1,nout
                     zval=zval+0.5d0*me*( &
                          conjg(aux(alpha,lmout))*aux(beta,lmout)+ &
                          conjg(aux(beta,lmout))*aux(alpha,lmout))
                  enddo
                  integ_nat(alpha,beta,ir)=rmesh(ir)*rmesh(ir)*zval
               endif
            enddo
         enddo

         ! Stored scalar-relativistic small-component norm for this band.
         ! SVEC is B(r)/r in WIEN2k's convention and RINT13 weights B^2 by
         ! 1/c^2, hence r^2/c^2 times the angular norm is the radial density.
         if(ensmall_scale.ne.0.d0) then
            do lm=1,nin
               tmp(ir)=tmp(ir)+rmesh(ir)*rmesh(ir)*c2inv*abs(svec(lm))**2
            enddo
         endif

         ! MASSV is a production term, but its diagnostic scale may be set to
         ! zero.  In that audit mode skip the expensive second angular map
         ! entirely instead of computing a tensor which will be discarded.
         if(massv_scale.ne.0.d0) then
            nsval=(0.d0,0.d0)
            naux=(0.d0,0.d0)
            do alpha=1,3
               do lm=1,nin
                  do imap=1,nmap(lm)
                     lmout=outmap(imap,lm)
                     nsval(alpha,lmout)=nsval(alpha,lmout)+ &
                          hmap(alpha,imap,lm)*svec(lm)
                  enddo
               enddo
               do beta=1,3
                  do lm=1,nout
                     do imap=1,nmap2(lm)
                        lmout=outmap2(imap,lm)
                        naux(alpha,beta,lmout)=naux(alpha,beta,lmout)+ &
                             hmap2(alpha,imap,lm)*aux(beta,lm)
                     enddo
                  enddo
               enddo
            enddo
            do alpha=1,3
               do beta=1,3
                  zval=(0.d0,0.d0)
                  do v=1,3
                     do lmout=1,nout2
                        zval=zval+0.5d0*me*( &
                             conjg(naux(alpha,v,lmout))*naux(beta,v,lmout)+ &
                             conjg(naux(beta,v,lmout))*naux(alpha,v,lmout))
                     enddo
                  enddo
                  do lmout=1,nout
                     zval=zval-0.5d0*c2inv*( &
                          conjg(nsval(alpha,lmout))*nsval(beta,lmout)+ &
                          conjg(nsval(beta,lmout))*nsval(alpha,lmout))
                  enddo
                  massv_tensor(alpha,beta,ir)=rmesh(ir)*rmesh(ir)* &
                       rmesh(ir)*dvry(ir)*real(zval)
               enddo
            enddo
         endif
      enddo

      do alpha=1,3
         do beta=1,3
               funct_re=real(integ_sr(alpha,beta,:))
               funct_im=aimag(integ_sr(alpha,beta,:))
               call charge(rmesh,1.d0,1.d0,funct_re,dx(jatom),imax,sum_re)
               call charge(rmesh,1.d0,1.d0,funct_im,dx(jatom),imax,sum_im)
               buf_sr(alpha,beta,num)=dcmplx(sum_re,sum_im)

               funct_re=real(integ_old(alpha,beta,:))
               funct_im=aimag(integ_old(alpha,beta,:))
               call charge(rmesh,1.d0,1.d0,funct_re,dx(jatom),imax,sum_re)
               call charge(rmesh,1.d0,1.d0,funct_im,dx(jatom),imax,sum_im)
               buf_old(alpha,beta,num)=dcmplx(sum_re,sum_im)
               buf_corr(alpha,beta,num)=buf_sr(alpha,beta,num)- &
                                         buf_old(alpha,beta,num)

               if(ennat_scale.ne.0.d0) then
                  funct_re=real(integ_nat(alpha,beta,:))
                  funct_im=aimag(integ_nat(alpha,beta,:))
                  call charge(rmesh,1.d0,1.d0,funct_re,dx(jatom),imax,sum_re)
                  call charge(rmesh,1.d0,1.d0,funct_im,dx(jatom),imax,sum_im)
                  buf_corr(alpha,beta,num)=buf_corr(alpha,beta,num)+ &
                       ennat_scale*band_energy(num)*dcmplx(sum_re,sum_im)
               endif
         enddo
      enddo

      ! Correct the APW surface flux.  At this point FVEC, DFVEC and QDVEC
      ! still contain the values at the outermost radial point.  HMAP is the
      ! exact spherical-harmonic representation of n_alpha*Y_lm.  Thus NPSI
      ! represents n_alpha*psi and NDFLUX represents
      ! n_alpha*(A*d_r(psi)-d_r(psi)).  The factor R_MT^2 is the surface
      ! measure.  Symmetrizing alpha,beta makes the tensor Hermitian before
      ! the real occupation sum is taken in l2main.
      npsi=(0.d0,0.d0)
      ndflux=(0.d0,0.d0)
      do alpha=1,3
         do lm=1,nin
            do imap=1,nmap(lm)
               lmout=outmap(imap,lm)
               npsi(alpha,lmout)=npsi(alpha,lmout)+ &
                    hmap(alpha,imap,lm)*fvec(lm)
               ndflux(alpha,lmout)=ndflux(alpha,lmout)+ &
                    hmap(alpha,imap,lm)*(qdvec(lm)-dfvec(lm))
            enddo
         enddo
      enddo
      do alpha=1,3
         do beta=1,3
            zval=(0.d0,0.d0)
            do lmout=1,nout
               zval=zval+0.5d0*( &
                    conjg(npsi(alpha,lmout))*ndflux(beta,lmout)+ &
                    conjg(npsi(beta,lmout))*ndflux(alpha,lmout))
            enddo
            buf_surf(alpha,beta,num)=rmesh(imax)*rmesh(imax)*zval

         enddo
      enddo

      if(small_surf_scale.ne.0.d0) then
         ! SMALL COMPONENT: unlike the large component it has no APW
         ! continuation outside the sphere.  Rebuild NPSI from SVEC so it
         ! represents n_alpha*s and contract it with the complete small-
         ! component gradient already held in SGRAD.
         npsi=(0.d0,0.d0)
         do alpha=1,3
            do lm=1,nin
               do imap=1,nmap(lm)
                  lmout=outmap(imap,lm)
                  npsi(alpha,lmout)=npsi(alpha,lmout)+ &
                       hmap(alpha,imap,lm)*svec(lm)
               enddo
            enddo
         enddo
         do alpha=1,3
            do beta=1,3
               zval=(0.d0,0.d0)
               do lmout=1,nout
                  zval=zval+0.5d0*c2inv*( &
                       conjg(npsi(alpha,lmout))*sgrad(beta,lmout)+ &
                       conjg(npsi(beta,lmout))*sgrad(alpha,lmout))
               enddo
               buf_surf(alpha,beta,num)=buf_surf(alpha,beta,num)+ &
                    small_surf_scale*rmesh(imax)*rmesh(imax)*zval
            enddo
         enddo
      endif


      if(ensmall_scale.ne.0.d0) then
         call charge(rmesh,1.d0,1.d0,tmp,dx(jatom),imax,small_norm)
         do alpha=1,3
            buf_surf(alpha,alpha,num)=buf_surf(alpha,alpha,num)- &
                 ensmall_scale*band_energy(num)*small_norm/3.d0
         enddo
      endif

      if(massv_scale.ne.0.d0) then
         do alpha=1,3
            do beta=1,3
               funct_re=massv_tensor(alpha,beta,:)
               call charge(rmesh,1.d0,1.d0,funct_re,dx(jatom),imax,sum_re)
               buf_corr(alpha,beta,num)=buf_corr(alpha,beta,num)- &
                    0.5d0*massv_scale*sum_re
            enddo
         enddo
      endif
   enddo
!$omp end do
   deallocate(fvec,dfvec,svec,dsvec,qdvec,quvec)
   deallocate(grad,aux,sgrad,npsi,ndflux,nsval,naux)
   deallocate(integ_sr,integ_old,integ_nat,funct_re,funct_im, &
              massv_tensor,tmp)
!$omp end parallel
   call walltim(band_wall_end)
   stress_rela_band_wall=stress_rela_band_wall+ &
        band_wall_end-band_wall_begin

   if(diag_status.eq.0 .and. len_trim(diag_env).gt.0 .and. &
      trim(diag_env).ne.'0') then
      write(6,'(A,I5,4(1X,ES23.15))') ':DBG_RELA_KERNEL:',jatom, &
           sum(abs(buf_sr(:,:,nemin:nemax))), &
           sum(abs(buf_old(:,:,nemin:nemax))), &
           sum(abs(buf_corr(:,:,nemin:nemax))), &
           sum(abs(buf_surf(:,:,nemin:nemax)))
   endif

   deallocate(nmap,outmap,hmap,dmap,nmap2,outmap2,hmap2, &
              rmesh,vry,dvry)
   deallocate(lf,ldf,sf,sdf,qd,qu)
   return
993 write(6,*) 'Invalid WIEN_STRESS_RELA_MASSV_SCALE: ',trim(massv_env)
    stop 'WIEN_STRESS_RELA_MASSV_SCALE must be a real number'
994 write(6,*) 'Invalid WIEN_STRESS_RELA_ENNAT_SCALE: ',trim(ennat_env)
    stop 'WIEN_STRESS_RELA_ENNAT_SCALE must be a real number'
995 write(6,*) 'Invalid WIEN_STRESS_RELA_ENSMALL_SCALE: ', &
               trim(ensmall_env)
    stop 'WIEN_STRESS_RELA_ENSMALL_SCALE must be a real number'
996 write(6,*) 'Invalid WIEN_STRESS_RELA_SMALL_SURF_SCALE: ', &
               trim(small_surf_env)
    stop 'WIEN_STRESS_RELA_SMALL_SURF_SCALE must be a real number'
997 write(6,*) 'Invalid WIEN_STRESS_RELA_OLD_SMALL_SCALE: ', &
               trim(small_scale_env)
    stop 'WIEN_STRESS_RELA_OLD_SMALL_SCALE must be a real number'
end subroutine val_rela_kinstrmt



