! Complete local bilinear form of Blaha G2 for a general Hermitian density.
!
! No on-shell KS equation is used on the band coefficients.  A pseudo-band
! may be any complex linear combination, and rhoE is independent of rho.
! The cache includes u, udot and every genuine LO primitive.  Angular
! quadrature is polynomial-exact for the retained harmonics; it is performed
! once per atom, never inside a k/band/Matsubara loop.  This also generates
! strain-allowed potential channels which are absent from case.nsh because
! the ordinary Hamiltonian Gaunt coefficient is zero.
module stress_g2_matrices
  use iso_fortran_env, only: dp=>real64
  implicit none
  private
  public :: prepare_g2_matrices, release_g2_matrices, evaluate_g2
  public :: gauss_legendre_stress, g2_angular_tables
  integer,parameter :: axis_a(6)=[1,2,3,2,1,1], axis_b(6)=[1,2,3,3,3,2]
  complex(dp),parameter :: one=(1._dp,0._dp),zero=(0._dp,0._dp)
  type local_g2
     integer :: n=0
     integer,allocatable :: ly(:),primitive(:),l(:),u_partner(:)
     real(dp),allocatable :: energy(:)
     complex(dp),allocatable :: s0(:,:),sr(:,:,:),sa(:,:,:)
     complex(dp),allocatable :: v0(:,:),vr(:,:,:),va(:,:,:)
  end type
  type(local_g2),allocatable,save :: atoms(:)
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

  subroutine prepare_g2_matrices()
    use param,only: nrad,nloat,lmax2,lomax,ncom,filename_V_vns
    use com,only: nat,rel
    use com_mpi,only: myrank,master
    use structure,only: jri,r0,dx,iatnr
    use lo,only: ilo,loor,elo
    use atspdt,only: el
    use w_atpar,only: w_rrad1,w_rrad2,w_rade1,w_rade2,w_a1lo,w_b1lo, &
         w_radp1,w_radp2,w_radpe1,w_radpe2,w_ap1lo,w_bp1lo,w_lm,w_lmmax, &
         ri_mat,retrieve_w_atpar,Vts,lpv,lv,mpv,mv,Nhmx
    use readPotential,only: init_V_vns,read_V_vns,close_V_vns
    integer :: jatom,nr,nh,n,l,m,p,i,j,iq,ilm,lpot,mpot,lfirst,lmmax, &
         ri,nradial,ir,ih,idx,jdx,lns
    integer,allocatable :: rad_id(:,:),rl(:),rp(:),rid(:)
    real(dp),allocatable :: f(:,:),g(:,:),rfp(:,:),rgp(:,:),r(:), &
         vraw(:,:),potential(:,:),srad(:,:),drad(:,:),vmat(:,:),wmat(:,:),tmp1(:),tmp2(:)
    complex(dp),allocatable :: ag(:,:),ar(:,:,:),aa(:,:,:),vreference(:,:)
    complex(dp) :: pp,pm
    real(dp) :: denom,angular_trace
    real(dp),external :: rint13n
    if(allocated(atoms)) return
    allocate(atoms(nat)); nh=(lmax2+1)**2
    allocate(ag(nh,nh),ar(nh,nh,6),aa(nh,nh,6))
    call init_V_vns(filename_V_vns,9929)
    do jatom=1,nat
       call retrieve_w_atpar(jatom,lfirst,lmmax)
       nr=jri(jatom)
       ! The represented nonspherical Hamiltonian has its own angular
       ! cutoff (LNSMAX in lapw1). Preserve it; completing strain-allowed
       ! triples must not silently add Hamiltonian channels above that cutoff.
       lns=lmax2
       if(Nhmx(jatom)>0) lns=min(lmax2,max(maxval(lpv(:Nhmx(jatom),jatom)), &
            maxval(lv(:Nhmx(jatom),jatom))))
       allocate(rad_id(2+nloat,0:lmax2));rad_id=0
       nradial=0; n=0
       do l=0,lmax2
          do p=1,2+nloat
             if(p>2) then
                if(l>lomax) cycle
                if(p-2>ilo(l)) cycle
                if(.not.loor(p-2,l)) cycle
             endif
             nradial=nradial+1;rad_id(p,l)=nradial;n=n+(2*l+1)
          enddo
       enddo
       atoms(jatom)%n=n
       allocate(atoms(jatom)%ly(n),atoms(jatom)%primitive(n),atoms(jatom)%l(n), &
            atoms(jatom)%u_partner(n),atoms(jatom)%energy(n),rid(n),rl(nradial),rp(nradial))
       allocate(f(nr,nradial),g(nr,nradial),rfp(nr,nradial),rgp(nr,nradial), &
            r(nr),vraw(nrad,0:ncom-1),potential(nrad,ncom),srad(nradial,nradial), &
            drad(nradial,nradial),vmat(nradial,nradial),wmat(nradial,nradial),tmp1(nr),tmp2(nr))
       do ir=1,nr
          r(ir)=r0(jatom)*exp((ir-1)*dx(jatom))
       enddo
       n=0
       do l=0,lmax2
          do p=1,2+nloat
             ri=rad_id(p,l)
             if(ri==0) cycle
             rl(ri)=l;rp(ri)=p
             select case(p)
             case(1)
                f(:,ri)=w_rrad1(:nr,l,jatom);g(:,ri)=w_rrad2(:nr,l,jatom)
                rfp(:,ri)=r*w_radp1(:nr,l,jatom);rgp(:,ri)=r*w_radp2(:nr,l,jatom)
             case(2)
                f(:,ri)=w_rade1(:nr,l,jatom);g(:,ri)=w_rade2(:nr,l,jatom)
                rfp(:,ri)=r*w_radpe1(:nr,l,jatom);rgp(:,ri)=r*w_radpe2(:nr,l,jatom)
             case default
                f(:,ri)=w_a1lo(:nr,p-2,l,jatom);g(:,ri)=w_b1lo(:nr,p-2,l,jatom)
                rfp(:,ri)=r*w_ap1lo(:nr,p-2,l,jatom);rgp(:,ri)=r*w_bp1lo(:nr,p-2,l,jatom)
             end select
             do m=-l,l
                n=n+1;rid(n)=ri
                atoms(jatom)%ly(n)=l*(l+1)+m+1
                atoms(jatom)%primitive(n)=p;atoms(jatom)%l(n)=l
                atoms(jatom)%energy(n)=el(l)
                if(p>2) atoms(jatom)%energy(n)=elo(l,p-2)
             enddo
          enddo
       enddo
       atoms(jatom)%u_partner=0
       do i=1,n
          if(atoms(jatom)%primitive(i)/=2) cycle
          do j=1,n
             if(atoms(jatom)%primitive(j)/=1) cycle
             if(atoms(jatom)%ly(i)==atoms(jatom)%ly(j)) atoms(jatom)%u_partner(i)=j
          enddo
       enddo
       do j=1,nradial
          do i=1,nradial
             srad(i,j)=rint13n(nr,f(:,i),g(:,i),f(:,j),g(:,j),r0(jatom),dx(jatom),rel)
             ! rfp=r^2*u'_bra; the derivative is never inferred by transposing.
             drad(i,j)=rint13n(nr,rfp(:,i),rgp(:,i),f(:,j),g(:,j),r0(jatom),dx(jatom),rel)
          enddo
       enddo
       allocate(atoms(jatom)%s0(n,n),atoms(jatom)%sr(n,n,6),atoms(jatom)%sa(n,n,6), &
            atoms(jatom)%v0(n,n),atoms(jatom)%vr(n,n,6),atoms(jatom)%va(n,n,6))
       atoms(jatom)%s0=zero;atoms(jatom)%v0=zero;atoms(jatom)%vr=zero;atoms(jatom)%va=zero
       ! Vharm=1 gives both spherical strain angular tables.
       call g2_angular_tables(lmax2,0,0,cmplx(sqrt(4*acos(-1._dp)),0._dp,dp),zero,ag,ar,aa)
       do j=1,n
          do i=1,n
             idx=atoms(jatom)%ly(i);jdx=atoms(jatom)%ly(j)
             if(idx==jdx) atoms(jatom)%s0(i,j)= &
                  ri_mat(atoms(jatom)%primitive(i),atoms(jatom)%primitive(j),atoms(jatom)%l(i),jatom)
             atoms(jatom)%sr(i,j,:)=drad(rid(i),rid(j))*ar(idx,jdx,:)
             atoms(jatom)%sa(i,j,:)=srad(rid(i),rid(j))*aa(idx,jdx,:)
          enddo
       enddo
       call read_V_vns(vraw,jatom,nrad,nr,w_lmmax(jatom))
       potential=0._dp
       potential(:,2:w_lmmax(jatom))=vraw(:,1:w_lmmax(jatom)-1)
       if(iatnr(jatom)>0) call cubic_potential_components(potential,w_lm(:,:,jatom),nr,w_lmmax(jatom))
       do ilm=2,w_lmmax(jatom)
          if(maxval(abs(potential(:nr,ilm)))==0._dp) cycle
          lpot=abs(w_lm(1,ilm,jatom));mpot=w_lm(2,ilm,jatom)
          call potential_phases(w_lm(1,ilm,jatom),mpot,iatnr(jatom)>0,pp,pm)
          call g2_angular_tables(lmax2,lpot,mpot,pp,pm,ag,ar,aa)
          do i=1,nradial
             tmp1=f(:,i)*potential(:nr,ilm);tmp2=g(:,i)*potential(:nr,ilm)
             do j=1,nradial
                vmat(i,j)=rint13n(nr,tmp1,tmp2,f(:,j),g(:,j),r0(jatom),dx(jatom),rel)
             enddo
             tmp1=rfp(:,i)*potential(:nr,ilm);tmp2=rgp(:,i)*potential(:nr,ilm)
             do j=1,nradial
                wmat(i,j)=rint13n(nr,tmp1,tmp2,f(:,j),g(:,j),r0(jatom),dx(jatom),rel)
             enddo
          enddo
          do j=1,n
             if(atoms(jatom)%l(j)>lns) cycle
             do i=1,n
                if(atoms(jatom)%l(i)>lns) cycle
                idx=atoms(jatom)%ly(i);jdx=atoms(jatom)%ly(j)
                atoms(jatom)%v0(i,j)=atoms(jatom)%v0(i,j)+vmat(rid(i),rid(j))*ag(idx,jdx)
                atoms(jatom)%vr(i,j,:)=atoms(jatom)%vr(i,j,:)+wmat(rid(i),rid(j))*ar(idx,jdx,:)
                atoms(jatom)%va(i,j,:)=atoms(jatom)%va(i,j,:)+vmat(rid(i),rid(j))*aa(idx,jdx,:)
             enddo
          enddo
       enddo
       ! Independently reconstruct the ordinary Hamiltonian matrix from
       ! its .nsh records.  This checks real/cubic harmonic phases, radial
       ! normalization and compatibility with the vectors that supplied E_l.
       allocate(vreference(n,n));vreference=zero
       do ih=1,Nhmx(jatom)
          do j=1,n
             if(atoms(jatom)%ly(j)/=lv(ih,jatom)*(lv(ih,jatom)+1)+mv(ih,jatom)+1) cycle
             do i=1,n
                if(atoms(jatom)%ly(i)/=lpv(ih,jatom)*(lpv(ih,jatom)+1)+mpv(ih,jatom)+1) cycle
                vreference(i,j)=vreference(i,j)+ &
                     Vts(atoms(jatom)%primitive(i),atoms(jatom)%primitive(j),ih,jatom)
             enddo
          enddo
       enddo
       denom=sqrt(sum(abs(vreference)**2))
       if(denom>1.e-20_dp.and.myrank==master) write(6,'(A,I6,2ES18.8)') &
            ':DMFT_G2_NSH_CHECK: ',jatom,denom, &
            sqrt(sum(abs(atoms(jatom)%v0-vreference)**2))/denom
       deallocate(vreference)
       angular_trace=maxval(abs(sum(atoms(jatom)%sa(:,:,1:3),dim=3)))+ &
            maxval(abs(sum(atoms(jatom)%va(:,:,1:3),dim=3)))
       if(angular_trace>1.e-9_dp) error stop 'G2 angular trace identity failed'
       if(myrank==master) write(6,'(A,3I6,2ES18.8)') ':DMFT_G2_COMPLETE_CACHE: ',jatom,n,lns, &
            angular_trace,maxval(abs(atoms(jatom)%v0-conjg(transpose(atoms(jatom)%v0))))
       deallocate(rad_id,rl,rp,rid,f,g,rfp,rgp,r,vraw,potential,srad,drad,vmat,wmat,tmp1,tmp2)
    enddo
    call close_V_vns()
    deallocate(ag,ar,aa)
  end subroutine

  subroutine potential_phases(signed_l,m,cubic,pp,pm)
    integer,intent(in) :: signed_l,m
    logical,intent(in) :: cubic
    complex(dp),intent(out) :: pp,pm
    integer :: parity
    parity=1
    if(mod(m,2)/=0) parity=-1
    pp=one;pm=one
    if(cubic) then
       if(any(abs(signed_l)==[3,7,9])) then
          pp=cmplx(0._dp,-1._dp,dp);pm=-pp
       endif
    else if(m/=0) then
       if(signed_l<0) then
          pp=cmplx(0._dp,-real(parity,dp),dp)/sqrt(2._dp)
          pm=cmplx(0._dp,1._dp,dp)/sqrt(2._dp)
       else
          pp=real(parity,dp)/sqrt(2._dp);pm=1._dp/sqrt(2._dp)
       endif
    endif
  end subroutine

  subroutine cubic_potential_components(v,lm,nr,nlm)
    ! WIEN cubic harmonics are constrained real-harmonic combinations.
    ! Transform their radial coefficients once to the complex-Y convention.
    ! Coefficients match dftKS CBCOMB; the .nsh file itself remains unchanged.
    real(dp),intent(inout) :: v(:,:)
    integer,intent(in) :: lm(:,:),nr,nlm
    real(dp) :: c(0:10,0:10),a(3),factor,total(nr)
    integer :: i,j,count,l,m
    c=0._dp
    c(0,0)=1; c(3,2)=1
    c(4,0)=sqrt(7._dp/3)/2;c(4,4)=sqrt(5._dp/3)/2
    c(6,0)=sqrt(.5_dp)/2;c(6,2)=sqrt(11._dp)/4
    c(6,4)=-sqrt(7._dp/2)/2;c(6,6)=-sqrt(5._dp)/4
    c(7,2)=sqrt(13._dp/6)/2;c(7,6)=sqrt(11._dp/6)/2
    c(8,0)=sqrt(33._dp)/8;c(8,4)=sqrt(7._dp/3)/4;c(8,8)=sqrt(65._dp/3)/8
    c(9,2)=sqrt(3._dp)/4;c(9,4)=sqrt(17._dp/6)/2
    c(9,6)=-sqrt(13._dp)/4;c(9,8)=-sqrt(7._dp/6)/2
    c(10,0)=sqrt(65._dp/6)/8;c(10,2)=sqrt(247._dp/6)/8
    c(10,4)=-sqrt(11._dp/2)/4;c(10,6)=sqrt(19._dp/3)/16
    c(10,8)=-sqrt(187._dp/6)/8;c(10,10)=-sqrt(85._dp)/16
    i=1
    do while(i<=nlm)
       l=abs(lm(1,i));m=lm(2,i)
       select case(l)
       case(0)
          i=i+1;cycle
       case(3)
          v(:nr,i)=v(:nr,i)/sqrt(2._dp);i=i+1;cycle
       case(4,6,7,9)
          count=2
       case(8,10)
          count=3
       case default
          error stop 'Unsupported cubic harmonic in complete G2'
       end select
       if(i+count-1>nlm) error stop 'Incomplete cubic harmonic group'
       total=0._dp
       do j=1,count
          a(j)=c(l,m+4*(j-1));total=total+a(j)*v(:nr,i+j-1)
       enddo
       do j=1,count
          factor=sqrt(2._dp)
          if(j==1.and.m==0) factor=1._dp
          v(:nr,i+j-1)=total*a(j)/factor
       enddo
       i=i+count
    enddo
  end subroutine

  subroutine evaluate_g2(jatom,first,last,alm,blm,clm,salm,sblm,sclm,band_energy, &
       spherical,nsp_match,nsp_radial,nsp_angular,sph_match,sph_radial,sph_angular,with_nsp)
    use structure,only: rotloc
    integer,intent(in) :: jatom,first,last
    complex(dp),intent(in) :: alm(:,:),blm(:,:),clm(:,:,:),salm(:,:,:),sblm(:,:,:),sclm(:,:,:,:)
    real(dp),intent(in) :: band_energy(:)
    complex(dp),intent(out) :: spherical(:,:,:),nsp_match(:,:,:),nsp_radial(:,:,:),nsp_angular(:,:,:), &
         sph_match(:,:,:),sph_radial(:,:,:),sph_angular(:,:,:)
    logical,intent(in) :: with_nsp
    complex(dp),allocatable :: c(:,:),d(:,:,:),hc(:,:),tmp(:,:)
    complex(dp) :: tensor(3,3),local_tensor(3,3),value
    integer :: n,nb,i,j,ly,p,iq,a,b,k,band
    if(.not.allocated(atoms)) error stop 'G2 cache not prepared'
    n=atoms(jatom)%n;nb=last-first+1
    allocate(c(n,nb),d(n,nb,6),hc(n,nb),tmp(n,nb))
    do i=1,n
       ly=atoms(jatom)%ly(i);p=atoms(jatom)%primitive(i)
       select case(p)
       case(1)
          c(i,:)=alm(ly,first:last);d(i,:,:)=salm(ly,first:last,:)
       case(2)
          c(i,:)=blm(ly,first:last);d(i,:,:)=sblm(ly,first:last,:)
       case default
          c(i,:)=clm(ly,first:last,p-2);d(i,:,:)=sclm(ly,first:last,p-2,:)
       end select
       ! The shared builder returns physical s.  Blaha D=s+gamma*A/2.
       do iq=1,3
          d(i,:,iq)=d(i,:,iq)+0.5_dp*c(i,:)
       enddo
       ! Matching derivatives arrive in the reference Cartesian axes;
       ! all G2 integral tables use the local harmonic axes.  Convert here
       ! before symmetry_stress_rotij later converts the complete tensor back.
       do k=1,nb
          do iq=1,6
             a=axis_a(iq);b=axis_b(iq);tensor(a,b)=d(i,k,iq);tensor(b,a)=d(i,k,iq)
          enddo
          local_tensor=matmul(rotloc(:,:,jatom),matmul(tensor,transpose(rotloc(:,:,jatom))))
          do iq=1,6
             d(i,k,iq)=local_tensor(axis_a(iq),axis_b(iq))
          enddo
       enddo
       hc(i,:)=atoms(jatom)%energy(i)*c(i,:)
    enddo
    do i=1,n
       j=atoms(jatom)%u_partner(i)
       if(j>0) hc(j,:)=hc(j,:)+c(i,:)
    enddo
    do k=1,nb
       hc(:,k)=hc(:,k)-band_energy(first+k-1)*c(:,k)
    enddo
    spherical=zero;sph_match=zero;sph_radial=zero;sph_angular=zero
    nsp_match=zero;nsp_radial=zero;nsp_angular=zero
    call zgemm('N','N',n,nb,n,one,atoms(jatom)%s0,n,hc,n,zero,tmp,n)
    do iq=1,6
       a=axis_a(iq);b=axis_b(iq)
       do k=1,nb
          band=first+k-1;value=sum(conjg(d(:,k,iq))*tmp(:,k))
          sph_match(a,b,band)=value;sph_match(b,a,band)=value
       enddo
    enddo
    do iq=1,6
       a=axis_a(iq);b=axis_b(iq)
       call zgemm('N','N',n,nb,n,one,atoms(jatom)%sr(:,:,iq),n,hc,n,zero,tmp,n)
       do k=1,nb
          band=first+k-1;value=sum(conjg(c(:,k))*tmp(:,k))
          sph_radial(a,b,band)=value;sph_radial(b,a,band)=value
       enddo
       call zgemm('N','N',n,nb,n,one,atoms(jatom)%sa(:,:,iq),n,hc,n,zero,tmp,n)
       do k=1,nb
          band=first+k-1;value=sum(conjg(c(:,k))*tmp(:,k))
          sph_angular(a,b,band)=value;sph_angular(b,a,band)=value
       enddo
    enddo
    spherical=sph_match+sph_radial+sph_angular
    if(with_nsp) then
       call zgemm('N','N',n,nb,n,one,atoms(jatom)%v0,n,c,n,zero,tmp,n)
       do iq=1,6
          a=axis_a(iq);b=axis_b(iq)
          do k=1,nb
             band=first+k-1;value=sum(conjg(d(:,k,iq))*tmp(:,k))
             nsp_match(a,b,band)=value;nsp_match(b,a,band)=value
          enddo
       enddo
       do iq=1,6
          a=axis_a(iq);b=axis_b(iq)
          call zgemm('N','N',n,nb,n,one,atoms(jatom)%vr(:,:,iq),n,c,n,zero,tmp,n)
          do k=1,nb
             band=first+k-1;value=sum(conjg(c(:,k))*tmp(:,k))
             nsp_radial(a,b,band)=value;nsp_radial(b,a,band)=value
          enddo
          call zgemm('N','N',n,nb,n,one,atoms(jatom)%va(:,:,iq),n,c,n,zero,tmp,n)
          do k=1,nb
             band=first+k-1;value=sum(conjg(c(:,k))*tmp(:,k))
             nsp_angular(a,b,band)=value;nsp_angular(b,a,band)=value
          enddo
       enddo
    endif
    deallocate(c,d,hc,tmp)
  end subroutine

  subroutine release_g2_matrices()
    if(allocated(atoms)) deallocate(atoms)
  end subroutine
end module stress_g2_matrices
