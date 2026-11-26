! Used by "make check-rela-operator". Compare raw complex tensors before
! occupation, rotations, symmetry averaging, or summing equivalent atoms.
program test_rela_operator
  use param, only: lmax2,lomax,nloat,nrad,nume
  use fact, only: fct
  use com, only: rel
  use atspdt, only: el
  use structure, only: jri,dx,r0
  use lo, only: ilo,loor,rlo,elo
  use xa, only: e
  use w_atpar, only: w_allocate,rela_operator,w_vr,w_rrad1,w_rrad2,w_rade1,w_rade2, &
       w_radp1,w_radp2,w_radpe1,w_radpe2,w_a1lo,w_b1lo,w_ap1lo,w_bp1lo
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: nb=5,nlm=(lmax2+1)**2
  complex*16 :: a(nlm,nb),b(nlm,nb),c(nlm,nb,nloat),out(3,3,nb,4),reference(3,3,nb,4)
  real*8 :: f,r,x,g,max_error
  integer :: atom,l,n,lm,i,j,visit,u
  character(10) :: mode
  character(512) :: path
  call get_command_argument(1,mode)
  call get_command_argument(2,path)
  if (mode=='write') then
     open(newunit=u,file=path,access='stream',form='unformatted',status='replace')
  else if (mode=='check') then
     open(newunit=u,file=path,access='stream',form='unformatted',status='old')
  else
     error stop 'expected write or check'
  endif
  nume=nb
  rel=.true.
  f=1.d0
  do i=1,49
     fct(2*i-1)=f
     f=f*i
  enddo
  call w_allocate(1,2)
  allocate(jri(2),r0(2),dx(2),e(nb),rela_operator(2))
  jri=[48,53] ! Exercise both endpoint rules in native CHARGE.
  dx=[0.045d0,0.05d0]
  r0=[0.15d0,0.1d0]
  do atom=1,2
     do i=1,nrad
        r=r0(atom)*exp(dx(atom)*(min(i,jri(atom))-1))
        w_vr(i,atom)=-3.d0*atom*exp(-r)
        do l=0,lmax2
           g=r**(l+1)*exp(-r*(0.4d0+atom/10.d0))
           w_rrad1(i,l,atom)=g
           w_rrad2(i,l,atom)=0.3d0*g
           w_rade1(i,l,atom)=r*g
           w_rade2(i,l,atom)=-0.2d0*g*r
           w_radp1(i,l,atom)=g*(l/r-0.4d0-atom/10.d0)
           w_radp2(i,l,atom)=0.3d0*w_radp1(i,l,atom)
           w_radpe1(i,l,atom)=r*w_radp1(i,l,atom)+g
           w_radpe2(i,l,atom)=-0.2d0*w_radpe1(i,l,atom)
           if (l>lomax) cycle
           do j=1,nloat
              w_a1lo(i,j,l,atom)=g/(r+j)
              w_b1lo(i,j,l,atom)=-0.4d0*w_a1lo(i,j,l,atom)
              w_ap1lo(i,j,l,atom)=(w_radp1(i,l,atom)-w_a1lo(i,j,l,atom))/(r+j)
              w_bp1lo(i,j,l,atom)=-0.4d0*w_ap1lo(i,j,l,atom)
           enddo
        enddo
     enddo
  enddo
  max_error=0.d0
  do visit=1,4
     atom=1+mod(visit-1,2)
     ilo=0
     ilo(0:lomax)=atom
     loor=.false.
     loor(:atom,:)=.true.
     rlo=.false.
     el=0.2d0*atom
     elo=-0.5d0*atom
     do n=1,nb
        e(n)=visit*0.2d0-n*0.1d0
        do lm=1,nlm
           x=0.317d0*(n*lm+visit)
           a(lm,n)=dcmplx(sin(x),cos(1.7d0*x))/nlm
           b(lm,n)=dcmplx(cos(x),sin(1.3d0*x))/nlm
           do j=1,nloat
              c(lm,n,j)=a(lm,n)+dcmplx(0.d0,0.2d0*j)*b(lm,n)
           enddo
        enddo
     enddo
     call val_rela_kinstrmt(2,nb-1,atom,a,b,c,out(:,:,:,1),out(:,:,:,2),out(:,:,:,3),out(:,:,:,4))
     if (mode=='write') then
        write(u) out
     else
        read(u) reference
        if (.not.all(ieee_is_finite(abs(out-reference)))) error stop 'nonfinite tensor'
        if (maxval(abs(reference))<1.d-8) error stop 'empty reference tensor'
        max_error=max(max_error,maxval(abs(out-reference)/(1.d0+abs(reference))))
        if (max_error>2.d-12) error stop 'cached RELA tensor differs from radial integration'
        if (.not.allocated(rela_operator(atom)%kernel)) error stop 'cache did not run'
     endif
  enddo
  close(u)
  if (mode=='check') print '(A,ES12.4)','Maximum scaled complex tensor error ',max_error
end program
