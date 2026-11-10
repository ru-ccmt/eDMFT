      subroutine pbe_vsigma_builtin(rhoup,rhodn,gradup,graddn, &
           gradtot,mu,kappa,beta,corr_scale,jspin,cutup,cutdn, &
           vsigup,vsigdn,vsigud)
      implicit none
      real*8, intent(in) :: rhoup,rhodn,gradup,graddn,gradtot
      real*8, intent(in) :: mu,kappa,beta,corr_scale
      integer, intent(in) :: jspin
      logical, intent(in) :: cutup,cutdn
      real*8, intent(out) :: vsigup,vsigdn,vsigud

      real*8, parameter :: pi=3.1415926535897932384626433832795d0
      real*8, parameter :: pi32=29.608813203268075856503472999628d0
      real*8, parameter :: thrd=1.d0/3.d0
      real*8, parameter :: thrd2=2.d0/3.d0
      real*8, parameter :: thrd4=4.d0/3.d0
      real*8, parameter :: ax=-0.738558766382022405884230032680836d0
      real*8, parameter :: gamma=0.03109069086965489503494086371273d0
      real*8, parameter :: gamspin=0.5198420997897463295344212145565d0
      real*8, parameter :: fzz=8.d0/(9.d0*gamspin)

      real*8 :: rho2,fk,s,p0,exunif,vsxup,vsxdn
      real*8 :: rho,zeta,phi,phi3,ks,twoksphi,tcor,t2,t4
      real*8 :: rs,eu,eurs,ep,eprs,alfm,alfrsm,fspin,z4,ec
      real*8 :: delta,b,b2,q4,q5,q8,q9,vsc,corr_factor

      vsigup=0.d0
      vsigdn=0.d0
      vsigud=0.d0
      vsxup=0.d0
      vsxdn=0.d0

! Exact spin scaling of exchange:
! E_x[n_up,n_dn] = (E_x[2*n_up]+E_x[2*n_dn])/2.
! For F_x=1+kappa-kappa/(1+(mu/kappa)*s^2), differentiation
! with respect to |grad n_spin|^2 gives the expressions below.
      if((.not.cutup).and.rhoup.gt.1.d-18.and. &
           abs(kappa).gt.1.d-14) then
        rho2=2.d0*rhoup
        fk=(pi32*rho2)**thrd
        s=2.d0*gradup/(2.d0*fk*rho2)
        p0=1.d0+(mu/kappa)*s*s
        exunif=ax*rho2**thrd
        vsxup=exunif*mu/(4.d0*fk*fk*rhoup*p0*p0)
      endif
      if((.not.cutdn).and.rhodn.gt.1.d-18.and. &
           abs(kappa).gt.1.d-14) then
        rho2=2.d0*rhodn
        fk=(pi32*rho2)**thrd
        s=2.d0*graddn/(2.d0*fk*rho2)
        p0=1.d0+(mu/kappa)*s*s
        exunif=ax*rho2**thrd
        vsxdn=exunif*mu/(4.d0*fk*fk*rhodn*p0*p0)
      endif

! PBE correlation depends on the total-density invariant
! sigma=|grad(n_up+n_dn)|^2.  Use precisely the PW91-LDA correlation
! interpolation and constants used by CORPBE, then differentiate its
! logarithmic gradient correction H analytically.
      vsc=0.d0
      rho=rhoup+rhodn
      if(rho.gt.1.d-18.and.abs(beta).gt.1.d-14.and. &
           abs(corr_scale).gt.1.d-14) then
        zeta=(rhoup-rhodn)/rho
        phi=((1.d0+zeta)**thrd2+(1.d0-zeta)**thrd2)/2.d0
        phi3=phi**3
        fk=(pi32*rho)**thrd
        rs=(0.75d0/(pi*rho))**thrd
        call gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0, &
             1.6382d0,0.49294d0,1.d0,rs,eu,eurs)
        call gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0, &
             3.3662d0,0.62517d0,1.d0,rs,ep,eprs)
        call gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0, &
             0.88026d0,0.49671d0,1.d0,rs,alfm,alfrsm)
        fspin=((1.d0+zeta)**thrd4+(1.d0-zeta)**thrd4-2.d0)/gamspin
        z4=zeta**4
        ec=eu*(1.d0-fspin*z4)+ep*fspin*z4- &
             alfm*fspin*(1.d0-z4)/fzz

        ks=sqrt(4.d0*fk/pi)
        twoksphi=2.d0*ks*phi
        tcor=gradtot/(twoksphi*rho)
        t2=tcor*tcor
        t4=t2*t2
        delta=beta/gamma
        b=delta/(exp(-ec/(phi3*gamma))-1.d0)
        b2=b*b
        q4=1.d0+b*t2
        q5=1.d0+b*t2+b2*t4
        q8=q5*q5+delta*q4*q5*t2
        q9=1.d0+2.d0*b*t2

! PBEC suppresses H completely if both spin channels are cut, and uses
! H/2 if only one is cut.  Apply the same rule to dH/dsigma.
        corr_factor=1.d0
        if(cutup.and.cutdn) corr_factor=0.d0
        if(cutup.neqv.cutdn) corr_factor=0.5d0
        vsc=corr_scale*corr_factor*beta*phi3*q9/ &
             (q8*twoksphi*twoksphi*rho)
      endif

      if(jspin.eq.1) then
! In the non-spin-polarized interface the independent invariant is
! |grad rho|^2.  vsxup was evaluated using rho_up=rho/2 and is twice
! the derivative with respect to this total-density invariant.
        vsigup=0.5d0*vsxup+vsc
        vsigdn=vsigup
        vsigud=0.d0
      else
        vsigup=vsxup+vsc
        vsigdn=vsxdn+vsc
! sigma_total=sigma_uu+2*sigma_ud+sigma_dd.
        vsigud=2.d0*vsc
      endif
      return
      end
