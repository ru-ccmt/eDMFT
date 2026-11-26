! Muffin-tin scalar-relativistic stress in the primitive (l,m,radial-channel)
! basis. Algebra follows VAL_RELA_KINSTRMT in patched LAPW2. In dmft2 the
! outer loop is over k, so each rank retains one sparse operator per atom.
module stress_rela_operator
  implicit none
  private
  public :: rela_atom_operator, build_rela_operator, apply_rela_operator

  integer, parameter :: sr=1, old=2, surface=3, mass_v=4, energy_norm=5

  type rela_atom_operator
     integer, allocatable :: lm(:), channel(:), bra(:), ket(:)
     complex*16, allocatable :: kernel(:,:,:,:) ! (a,b,term,pair)
     real*8, allocatable :: small_overlap(:)
     real*8 :: small_energy_scale, natural_energy_scale, mass_scale
  end type

contains

  subroutine build_rela_operator(atom,r,dx,lf,ldf,sf,sdf,qd,qu,dv,ilo,loor, &
       nmap,outmap,hmap,dmap,nmap2,outmap2,hmap2,old_scale,surf_scale,en_scale,nat_scale,mass_scale)
    use param, only: nrad
    type(rela_atom_operator), intent(out) :: atom
    real*8, intent(in) :: r(:),dx,lf(:,0:,:),ldf(:,0:,:),sf(:,0:,:),sdf(:,0:,:), &
         qd(:,0:,:),qu(:,0:,:),dv(:),old_scale,surf_scale,en_scale,nat_scale,mass_scale
    integer, intent(in) :: ilo(0:),nmap(:),outmap(:,:),nmap2(:),outmap2(:,:)
    logical, intent(in) :: loor(:,0:)
    complex*16, intent(in) :: hmap(:,:,:),dmap(:,:,:),hmap2(:,:,:)
    real*8, parameter :: c2inv=1.d0/(137.0359895d0**2), me=0.25d0*c2inv
    real*8 :: mesh(nrad),work(nrad),edge_area
    real*8, allocatable :: grad(:,:,:),flux(:,:,:),small_grad(:,:,:),small(:,:),large_edge(:)
    real*8, allocatable :: integral(:,:,:,:,:),mass_small(:,:),small_overlap(:,:),pair_overlap(:)
    complex*16, allocatable :: angular(:,:,:,:),moment(:,:,:,:,:),kernel(:,:,:,:)
    complex*16 :: angular_pair(3,3,2,2),moment_pair(3,3,2,2),op(3,3,5)
    integer, allocatable :: channel_id(:,:),basis_radial(:),first(:),last(:),bra(:),ket(:)
    integer :: nr,lmax,nch,nlm,nout,nout2,nradial,nbasis,npair,l,ch,lm,radial_id
    integer :: p,q,i,j,a,b,v,h,k,out,dest,ip,iq,lp,lq,nterm,capacity

    atom%small_energy_scale=en_scale
    atom%natural_energy_scale=nat_scale
    atom%mass_scale=mass_scale
    nr=size(r)
    mesh(:nr)=r
    edge_area=r(nr)**2
    lmax=ubound(lf,2)
    nch=size(lf,3)
    nlm=size(nmap)
    nout=size(nmap2)
    nout2=(lmax+3)**2

    ! A radial channel is shared by all m of the same l. Only u, udot and
    ! the local orbitals actually present on this atom occur in the basis.
    allocate(channel_id(0:lmax,nch))
    channel_id=0
    nradial=0
    nbasis=0
    do l=0,lmax
       do ch=1,nch
          if (ch>2) then
             if (ch-2>ilo(l)) cycle
             if (.not.loor(ch-2,l)) cycle
          endif
          nradial=nradial+1
          channel_id(l,ch)=nradial
          nbasis=nbasis+2*l+1
       enddo
    enddo
    allocate(grad(nr,2,nradial),flux(nr,2,nradial),small_grad(nr,2,nradial),small(nr,nradial),large_edge(nradial))
    allocate(atom%lm(nbasis),atom%channel(nbasis),basis_radial(nbasis),first(nlm),last(nlm))
    p=0
    do lm=1,nlm
       l=int(sqrt(dble(lm-1)))
       first(lm)=p+1
       do ch=1,nch
          radial_id=channel_id(l,ch)
          if (radial_id==0) cycle
          p=p+1
          atom%lm(p)=lm
          atom%channel(p)=ch
          basis_radial(p)=radial_id
          if (lm/=l*l+1) cycle
          ! Cartesian gradients have two radial factors: derivative and f/r.
          grad(:,1,radial_id)=ldf(:,l,ch)
          grad(:,2,radial_id)=lf(:,l,ch)/r
          flux(:,1,radial_id)=qd(:,l,ch)
          flux(:,2,radial_id)=qu(:,l,ch)
          small_grad(:,1,radial_id)=sdf(:,l,ch)
          small_grad(:,2,radial_id)=sf(:,l,ch)/r
          small(:,radial_id)=sf(:,l,ch)
          large_edge(radial_id)=lf(nr,l,ch)
       enddo
       last(lm)=p
    enddo

    ! Integrate once with precisely the dmft2 logarithmic-mesh quadrature,
    ! including its endpoint rule on an even mesh. The five products below
    ! are grad-grad, grad-flux, small-small, flux-flux and r V' flux-flux.
    allocate(integral(2,2,nradial,nradial,5),mass_small(nradial,nradial),small_overlap(nradial,nradial))
    integral=0.d0
    mass_small=0.d0
    small_overlap=0.d0
    do q=1,nradial
       do p=1,nradial
          do j=1,2
             do i=1,2
                work(:nr)=r*r*grad(:,i,p)*grad(:,j,q)
                call charge(mesh,1.d0,1.d0,work,dx,nr,integral(i,j,p,q,1))
                work(:nr)=r*r*grad(:,i,p)*flux(:,j,q)
                call charge(mesh,1.d0,1.d0,work,dx,nr,integral(i,j,p,q,2))
                work(:nr)=r*r*small_grad(:,i,p)*small_grad(:,j,q)
                call charge(mesh,1.d0,1.d0,work,dx,nr,integral(i,j,p,q,3))
                if (nat_scale/=0.d0) then
                   work(:nr)=r*r*flux(:,i,p)*flux(:,j,q)
                   call charge(mesh,1.d0,1.d0,work,dx,nr,integral(i,j,p,q,4))
                endif
                if (mass_scale/=0.d0) then
                   work(:nr)=r*r*r*dv*flux(:,i,p)*flux(:,j,q)
                   call charge(mesh,1.d0,1.d0,work,dx,nr,integral(i,j,p,q,5))
                endif
             enddo
          enddo
          if (mass_scale/=0.d0) then
             work(:nr)=r*r*r*dv*small(:,p)*small(:,q)
             call charge(mesh,1.d0,1.d0,work,dx,nr,mass_small(p,q))
          endif
          if (en_scale/=0.d0) then
             work(:nr)=r*r*c2inv*small(:,p)*small(:,q)
             call charge(mesh,1.d0,1.d0,work,dx,nr,small_overlap(p,q))
          endif
       enddo
    enddo

    ! Angular coefficients in the output Y_LM basis. A second direction
    ! cosine is required for the n_a n_b moment of the mass-density term.
    allocate(angular(nout,3,2,nlm))
    angular=0.d0
    do lm=1,nlm
       do k=1,nmap(lm)
          out=outmap(k,lm)
          angular(out,:,1,lm)=hmap(:,k,lm)
          angular(out,:,2,lm)=dmap(:,k,lm)
       enddo
    enddo
    if (mass_scale/=0.d0) then
       allocate(moment(nout2,3,3,2,nlm))
       moment=0.d0
       do lm=1,nlm
          do k=1,nmap(lm)
             out=outmap(k,lm)
             do h=1,nmap2(out)
                dest=outmap2(h,out)
                do i=1,2
                   do v=1,3
                      moment(dest,:,v,i,lm)=moment(dest,:,v,i,lm)+hmap2(:,h,out)*angular(out,v,i,lm)
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    ! Dense upper bound is temporary for a single atom. Keep only nonzero
    ! pairs on the rank, with both bra-ket orders for complex coefficients.
    nterm=mass_v
    if (nat_scale/=0.d0) nterm=energy_norm
    capacity=nbasis**2
    allocate(kernel(3,3,nterm,capacity),bra(capacity),ket(capacity),pair_overlap(capacity))
    npair=0
    do iq=1,nlm
       do ip=1,nlm
          moment_pair=0.d0
          do j=1,2
             do i=1,2
                do b=1,3
                   do a=1,3
                      angular_pair(a,b,i,j)=0.5d0*(dot_product(angular(:,a,i,ip),angular(:,b,j,iq)) &
                                                   +dot_product(angular(:,b,i,ip),angular(:,a,j,iq)))
                      if (mass_scale/=0.d0) then
                         do v=1,3
                            moment_pair(a,b,i,j)=moment_pair(a,b,i,j)+0.5d0*( &
                                 dot_product(moment(:,a,v,i,ip),moment(:,b,v,j,iq)) &
                                +dot_product(moment(:,b,v,i,ip),moment(:,a,v,j,iq)))
                         enddo
                      endif
                   enddo
                enddo
             enddo
          enddo
          if (ip/=iq .and. all(angular_pair==0.d0) .and. all(moment_pair==0.d0)) cycle
          do q=first(iq),last(iq)
             do p=first(ip),last(ip)
                lp=basis_radial(p)
                lq=basis_radial(q)
                op=0.d0
                do j=1,2
                   do i=1,2
                      op(:,:,sr)=op(:,:,sr)+0.5d0*angular_pair(:,:,i,j)*(integral(i,j,lp,lq,2)+integral(j,i,lq,lp,2))
                      op(:,:,old)=op(:,:,old)+angular_pair(:,:,i,j)*(integral(i,j,lp,lq,1)+old_scale*c2inv*integral(i,j,lp,lq,3))
                      if (nat_scale/=0.d0) op(:,:,energy_norm)=op(:,:,energy_norm)+me*angular_pair(:,:,i,j)*integral(i,j,lp,lq,4)
                      if (mass_scale/=0.d0) op(:,:,mass_v)=op(:,:,mass_v)+me*moment_pair(:,:,i,j)*integral(i,j,lp,lq,5)
                   enddo
                enddo
                op(:,:,mass_v)=op(:,:,mass_v)-c2inv*angular_pair(:,:,1,1)*mass_small(lp,lq)
                op(:,:,surface)=edge_area*angular_pair(:,:,1,1)*large_edge(lp)*(flux(nr,1,lq)-grad(nr,1,lq))
                if (surf_scale/=0.d0) then
                   do j=1,2
                      op(:,:,surface)=op(:,:,surface)+surf_scale*edge_area*c2inv*small(nr,lp)* &
                              angular_pair(:,:,1,j)*small_grad(nr,j,lq)
                   enddo
                endif
                if (all(op==0.d0) .and. (ip/=iq .or. small_overlap(lp,lq)==0.d0)) cycle
                npair=npair+1
                bra(npair)=p
                ket(npair)=q
                kernel(:,:,:,npair)=op(:,:,:nterm)
                pair_overlap(npair)=merge(small_overlap(lp,lq),0.d0,ip==iq)
             enddo
          enddo
       enddo
    enddo
    atom%bra=bra(:npair)
    atom%ket=ket(:npair)
    atom%kernel=kernel(:,:,:,:npair)
    if (en_scale/=0.d0) atom%small_overlap=pair_overlap(:npair)
  end subroutine

  subroutine apply_rela_operator(atom,nemin,nemax,alm,blm,clm,energy,buf_corr,buf_sr,buf_old,buf_surf)
    type(rela_atom_operator), intent(in) :: atom
    integer, intent(in) :: nemin,nemax
    complex*16, intent(in) :: alm(:,:),blm(:,:),clm(:,:,:)
    real*8, intent(in) :: energy(:)
    ! Caller zeroes every output, including unused band slices.
    complex*16, intent(inout) :: buf_corr(:,:,:),buf_sr(:,:,:),buf_old(:,:,:),buf_surf(:,:,:)
    complex*16, allocatable :: coeff(:)
    complex*16 :: sums(3,3,5),product,overlap
    integer :: n,p,ch,k,t,a
    !$omp parallel default(shared) private(n,p,ch,k,t,a,coeff,sums,product,overlap)
    allocate(coeff(size(atom%lm)))
    !$omp do schedule(static)
    do n=nemin,nemax
       do p=1,size(coeff)
          ch=atom%channel(p)
          select case (ch)
          case (1)
             coeff(p)=alm(atom%lm(p),n)
          case (2)
             coeff(p)=blm(atom%lm(p),n)
          case default
             coeff(p)=clm(atom%lm(p),n,ch-2)
          end select
       enddo
       sums=0.d0
       overlap=0.d0
       do k=1,size(atom%bra)
          product=conjg(coeff(atom%bra(k)))*coeff(atom%ket(k))
          do t=1,size(atom%kernel,3)
             sums(:,:,t)=sums(:,:,t)+product*atom%kernel(:,:,t,k)
          enddo
          if (allocated(atom%small_overlap)) overlap=overlap+product*atom%small_overlap(k)
       enddo
       buf_sr(:,:,n)=sums(:,:,sr)
       buf_old(:,:,n)=sums(:,:,old)
       buf_surf(:,:,n)=sums(:,:,surface)
       buf_corr(:,:,n)=sums(:,:,sr)-sums(:,:,old)-0.5d0*atom%mass_scale*real(sums(:,:,mass_v),kind=8)
       ! Band energy is attached at application time. The same operator works
       ! for density pseudo-bands and for the two signed window responses.
       if (atom%natural_energy_scale/=0.d0) then
          buf_corr(:,:,n)=buf_corr(:,:,n)+atom%natural_energy_scale*energy(n)*sums(:,:,energy_norm)
       endif
       if (allocated(atom%small_overlap)) then
          do a=1,3
             buf_surf(a,a,n)=buf_surf(a,a,n)-atom%small_energy_scale*energy(n)*real(overlap,kind=8)/3.d0
          enddo
       endif
    enddo
    !$omp end do
    deallocate(coeff)
    !$omp end parallel
  end subroutine
end module stress_rela_operator
