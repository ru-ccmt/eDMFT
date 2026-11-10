program test_complete_g2
 use param
 use structure
 use lo
 use atspdt
 use w_atpar
 use stress_g2_matrices
 implicit none
 integer,parameter :: nb=3,nh=(lmax2+1)**2
 integer :: i,j,l,p,ir,n,unit,iq
 real(8) :: r(nrad),f(nrad,4,0:lmax2),deriv(nrad,4,0:lmax2),z(nrad),alpha
 real(8) :: benergy(nb)
 real(8),external :: rint13n
 complex(8) :: a(nh,nb),b(nh,nb),c(nh,nb,nloat),sa(nh,nb,6),sb(nh,nb,6),sc(nh,nb,nloat,6)
 complex(8) :: sph(3,3,nb),ns0(3,3,nb),nsr(3,3,nb),nsa(3,3,nb)
 complex(8) :: sp0(3,3,nb),spr(3,3,nb),spa(3,3,nb)
 r0(1)=2.1d0*exp(-(nrad-1)*dx(1));z=0
 rotloc=0;do i=1,3;rotloc(i,i,1)=1;enddo
 ! A nontrivial local frame also tests the Cartesian-to-local matching map.
 rotloc(1,1,1)=cos(.37d0);rotloc(2,2,1)=cos(.37d0)
 rotloc(1,2,1)=-sin(.37d0);rotloc(2,1,1)=sin(.37d0)
 do ir=1,nrad;r(ir)=r0(1)*exp((ir-1)*dx(1));enddo
 w_lm(:,1,1)=[0,0];w_lm(:,2,1)=[2,0]
 w_rrad2=0;w_rade2=0;w_b1lo=0;w_radp2=0;w_radpe2=0;w_bp1lo=0
 do l=0,lmax2
 el(l)=.3d0+.11d0*l
 do p=1,4
 alpha=.7d0+.2d0*p
 f(:,p,l)=r**(l+1)*exp(-alpha*r)
 deriv(:,p,l)=(l-alpha*r)*f(:,p,l)/r
 enddo
 w_rrad1(:,l,1)=f(:,1,l);w_rade1(:,l,1)=f(:,2,l)
 w_radp1(:,l,1)=deriv(:,1,l);w_radpe1(:,l,1)=deriv(:,2,l)
 do p=1,nloat
 w_a1lo(:,p,l,1)=f(:,p+2,l);w_ap1lo(:,p,l,1)=deriv(:,p+2,l)
 elo(l,p)=1.3d0+.21d0*l+.6d0*p
 enddo
 do i=1,4
 do j=1,4
 ri_mat(i,j,l,1)=rint13n(nrad,f(:,i,l),z,f(:,j,l),z,r0(1),dx(1),.false.)
 enddo
 enddo
 enddo
 call prepare_g2_matrices()
 n=atoms(1)%n
 open(newunit=unit,file='g2_fixture.bin',access='stream',form='unformatted',status='replace')
 write(unit) n,nb,nrad,r0(1),dx(1),rotloc(:,:,1)
 write(unit) atoms(1)%ly,atoms(1)%l,atoms(1)%primitive,atoms(1)%energy,atoms(1)%u_partner
 write(unit) atoms(1)%s0,atoms(1)%sr,atoms(1)%sa,atoms(1)%v0,atoms(1)%vr,atoms(1)%va
 do i=1,nh
 do j=1,nb
 a(i,j)=cmplx(sin(.17d0*i*j),cos(.13d0*(i+j)),8)
 b(i,j)=cmplx(cos(.19d0*i*j),sin(.31d0*(i+j)),8)*.3d0
 do p=1,nloat
 c(i,j,p)=cmplx(cos(.11d0*(i+j+p)),sin(.23d0*i*j*p),8)*.2d0
 enddo
 do iq=1,6
 sa(i,j,iq)=a(i,j)*cmplx(.1d0*iq,.03d0*(i-j),8)
 sb(i,j,iq)=b(i,j)*cmplx(.07d0*iq,.02d0*(j-i),8)
 do p=1,nloat
 sc(i,j,p,iq)=c(i,j,p)*cmplx(.05d0*(p+iq),.01d0*(j-i),8)
 enddo
 enddo
 enddo
 enddo
 benergy=[-.4d0,.2d0,1.1d0]
 call evaluate_g2(1,1,nb,a,b,c,sa,sb,sc,benergy,sph,ns0,nsr,nsa,sp0,spr,spa,.true.)
 write(unit) a,b,c,sa,sb,sc,benergy
 write(unit) sp0,spr,spa,ns0,nsr,nsa
 close(unit)
 call release_g2_matrices()
end program
