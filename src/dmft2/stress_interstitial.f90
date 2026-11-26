! Interstitial kinetic-stress kernels adapted from WIEN2k 24 SRC_lapw2.
! This file contains the complex-vector version used by dmft2.

SUBROUTINE step_fun(X1,Y1,Z1,N1,thetagp)       

       use structure, only: mult,pos,vol,rmt,br1
       use com, only: nat
       use defs, only: pi       
       
       implicit none
       integer,intent(in)   ::  N1 
       integer,intent(in)   :: X1(N1),Y1(N1),Z1(N1)
            
       real*8               :: gnorm(N1)
       real*8               :: A(3), gxyz(3)
       real*8               :: arg, q, vat(nat)
       integer              :: ig,ig1, iat, imu, index, n, k
       
      complex*16, intent(out) :: thetagp(N1)

complex*16 :: phs
                    
 thetagp = (0d0,0d0)
                                                           
       do iat = 1, nat
          vat(iat) = (4d0*pi/3.0)*rmt(iat)**3
       enddo
                  
        do ig = 1, N1
           gxyz(1) = X1(ig)
           gxyz(2) = Y1(ig)
           gxyz(3) = Z1(ig)
           A(1) = X1(ig)*BR1(1,1) + Y1(ig)*BR1(1,2) + Z1(ig)*BR1(1,3) 
           A(2) = X1(ig)*BR1(2,1) + Y1(ig)*BR1(2,2) + Z1(ig)*BR1(2,3)
           A(3) = X1(ig)*BR1(3,1) + Y1(ig)*BR1(3,2) + Z1(ig)*BR1(3,3)
           gnorm(ig) = sqrt(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
           index = 0
                 if(gnorm(ig) .eq. 0d0) then
                      thetagp(ig) = 1d0 - sum(mult(1:nat)*vat(1:nat))/vol                       
                 else
                 do iat = 1, nat                     
          phs = (0d0,0d0)                    
                    do imu = 1, mult(iat)
                       index = index + 1
                       arg  = (gxyz(1)*pos(1,index) + &
                               gxyz(2)*pos(2,index) + &
                               gxyz(3)*pos(3,index))*2d0*pi                        
             phs = phs + dcmplx(cos(-arg),sin(-arg))

                    enddo
                    q = gnorm(ig)*rmt(iat)                    
                    thetagp(ig) = thetagp(ig) + 3d0*vat(iat)* &
                                (q*cos(q) - sin(q))/(q*q*q)*phs/vol                      
                 enddo 
                 endif                           
        enddo                  
       ! step function which has been calculating here is the same as the one in the mixer but tricky part is that
       ! in the mixer, step function is actually vol time the thetag.
END


! Analytic derivative of the interstitial step-function overlap for a fixed,
! spherical muffin-tin partition.  For a reciprocal difference g,
!
!   Theta_g = delta_g0 + sum_a 3 V_a/Omega f(q_a) exp(-i g.tau_a),
!   f(q)=(q cos(q)-sin(q))/q^3,   q_a=|g| R_a.
!
! Fractional atomic positions keep the phase invariant under homogeneous
! strain.  The reciprocal length and cell volume give
!
!   |g|_,alpha beta = -g_alpha g_beta/|g|,
!   (1/Omega)_,alpha beta = -delta_alpha beta/Omega.
!
! This term is absent from Blaha's on-shell IBS kernel because the
! interstitial sqrt(Omega)*plane wave is strain invariant.  It is nevertheless
! compulsory when the overlap derivative is isolated off shell by M_Sigma.
SUBROUTINE step_fun_stress_derivative(X1,Y1,Z1,N1,dthetagp)
  use structure, only: mult,pos,vol,rmt,br1
  use com, only: nat
  use defs, only: pi
  implicit none
  integer,intent(in) :: N1
  integer,intent(in) :: X1(N1),Y1(N1),Z1(N1)
  complex*16,intent(out) :: dthetagp(3,3,N1)
  integer :: ig,iat,imu,index,alpha,beta
  real*8 :: gfrac(3),gcart(3),gnorm,arg,q,vat,fq,dfq,qder
  complex*16 :: phase,central

  dthetagp=(0.d0,0.d0)
  do ig=1,N1
     gfrac=(/dble(X1(ig)),dble(Y1(ig)),dble(Z1(ig))/)
     gcart=matmul(BR1,gfrac)
     gnorm=sqrt(dot_product(gcart,gcart))
     if (gnorm.lt.1.d-13) then
        do iat=1,nat
           vat=4.d0*pi*rmt(iat)**3/3.d0
           do alpha=1,3
              dthetagp(alpha,alpha,ig)=dthetagp(alpha,alpha,ig)+ &
                   dble(mult(iat))*vat/vol
           enddo
        enddo
        cycle
     endif

     index=0
     do iat=1,nat
        vat=4.d0*pi*rmt(iat)**3/3.d0
        phase=(0.d0,0.d0)
        do imu=1,mult(iat)
           index=index+1
           arg=2.d0*pi*dot_product(gfrac,pos(:,index))
           phase=phase+dcmplx(cos(-arg),sin(-arg))
        enddo
        q=gnorm*rmt(iat)
        fq=(q*cos(q)-sin(q))/(q*q*q)
        dfq=((3.d0-q*q)*sin(q)-3.d0*q*cos(q))/(q**4)
        central=3.d0*vat*fq*phase/vol
        do alpha=1,3
           do beta=1,3
              qder=-rmt(iat)*gcart(alpha)*gcart(beta)/gnorm
              dthetagp(alpha,beta,ig)=dthetagp(alpha,beta,ig)+ &
                   3.d0*vat*dfq*qder*phase/vol
              if (alpha.eq.beta) dthetagp(alpha,beta,ig)= &
                   dthetagp(alpha,beta,ig)-central
           enddo
        enddo
     enddo
  enddo
END SUBROUTINE step_fun_stress_derivative


! Contract the complete interstitial overlap derivative with one set of
! signed pseudo-bands.  BAND_WEIGHT already contains the k/spin/Matsubara
! weight.  DTHETA is a full Hermitian derivative, so unlike a one-sided
! muffin-tin Pulay kernel there is no extra factor two here.
SUBROUTINE val_overlap_interstitial_stress(isize,overlap_stress,nemin, &
     nemax,band_weight)
  use xa3, only: a,k3
  use param, only: nume
  implicit none
  integer,intent(in) :: isize,nemin,nemax
  real*8,intent(in) :: band_weight(nume)
  real*8,intent(out) :: overlap_stress(3,3)
  integer :: npair,ig,igp,index,num,alpha,beta
  integer,allocatable :: gx(:),gy(:),gz(:)
  complex*16,allocatable :: dtheta(:,:,:)
  complex*16 :: value

  npair=isize*isize
  allocate(gx(npair),gy(npair),gz(npair),dtheta(3,3,npair))
  index=0
  do ig=1,isize
     do igp=1,isize
        index=index+1
        gx(index)=k3(1,ig)-k3(1,igp)
        gy(index)=k3(2,ig)-k3(2,igp)
        gz(index)=k3(3,ig)-k3(3,igp)
     enddo
  enddo
  call step_fun_stress_derivative(gx,gy,gz,npair,dtheta)

  overlap_stress=0.d0
  do num=nemin,nemax
     do alpha=1,3
        do beta=1,3
           value=(0.d0,0.d0)
           index=0
           do ig=1,isize
              do igp=1,isize
                 index=index+1
                 value=value+conjg(a(ig,num))*a(igp,num)* &
                      dtheta(alpha,beta,index)
              enddo
           enddo
           overlap_stress(alpha,beta)=overlap_stress(alpha,beta)+ &
                band_weight(num)*dble(value)
        enddo
     enddo
  enddo
  deallocate(gx,gy,gz,dtheta)
END SUBROUTINE val_overlap_interstitial_stress
       

! For each plane wave G contract the mask with all the current pseudo-bands:
!   M(G,n) = sum_Gp Theta(G-Gp) A(Gp,n).
! The band sum of Re[A(G,n)* M(G,n)] (conjugating A) then multiplies the
! same Cartesian outer product for every tensor component. BAND_WEIGHT is
! supplied by the caller: it can be signed for complementary-window and
! projector-normalization terms. The slab bounds the workspace at O(N_G).
subroutine val_kin_interstitial(isize, kin_int_fn, nemin, nemax, band_weight)
  use xa3, only: a, k3
  use stress_pw, only: F1p
  use param, only: nume
  implicit none
  integer, intent(in) :: isize, nemin, nemax
  real*8, intent(in) :: band_weight(nume)
  real*8, intent(inout) :: kin_int_fn(3,3)
  integer, parameter :: block_size=32
  integer :: first, rows, row, axis, g, band, nbands, alpha, beta
  integer, allocatable :: delta_g(:,:,:)
  complex*16, allocatable :: theta(:,:), masked(:,:)
  complex*16, parameter :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
  real*8 :: weighted_density

  nbands=nemax-nemin+1
  if (nbands<=0 .or. isize<=0) return
  allocate(delta_g(isize,block_size,3), theta(isize,block_size), masked(block_size,nbands))
  do first=1,isize,block_size
     rows=min(block_size,isize-first+1)
     do axis=1,3
        do row=1,rows
           delta_g(:,row,axis)=k3(axis,first+row-1)-k3(axis,1:isize)
        enddo
     enddo
     call step_fun(delta_g(1,1,1),delta_g(1,1,2),delta_g(1,1,3),isize*rows,theta)
     ! theta is stored as (Gp,G), with argument G-Gp. Use its transpose
     ! without conjugating the mask, including for noncentrosymmetric cells.
     call zgemm('T','N',rows,nbands,isize,one,theta,isize,a(1,nemin),size(a,1),zero,masked,block_size)
     do row=1,rows
        g=first+row-1
        weighted_density=0.d0
        do band=nemin,nemax
           weighted_density=weighted_density+band_weight(band)* &
                real(conjg(a(g,band))*masked(row,band-nemin+1),kind=8)
        enddo
        do beta=1,3
           do alpha=1,3
              kin_int_fn(alpha,beta)=kin_int_fn(alpha,beta)+ &
                   2.d0*weighted_density*F1p(g,alpha)*F1p(g,beta)
           enddo
        enddo
     enddo
  enddo
end subroutine val_kin_interstitial
