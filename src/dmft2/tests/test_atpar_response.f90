program test_atpar_response
  use iso_fortran_env, only: real64
  use outwin_response
  implicit none
  integer, parameter :: dp=real64,n=401,lmax=3,nlo=1
  type(radial_pair_context) :: ctx
  real(dp) :: r(n),potential(n),r0,dx,energy,charge,error,largest,values(n,4,2)
  real(dp) :: el(0:lmax),elo(0:lmax,nlo)
  real(dp) :: u(n,0:lmax),du(n,0:lmax),b(n,0:lmax),db(n,0:lmax)
  real(dp) :: up(n,0:lmax),dup(n,0:lmax),bp(n,0:lmax),dbp(n,0:lmax)
  real(dp) :: alo_r(n,nlo,0:lmax),blo_r(n,nlo,0:lmax),aplo(n,nlo,0:lmax),bplo(n,nlo,0:lmax)
  real(dp) :: boundary(nlo+2,0:lmax),slope(nlo+2,0:lmax),overlap(nlo+2,nlo+2,0:lmax)
  real(dp) :: alo(0:lmax,nlo),blo(0:lmax,nlo),clo(0:lmax,nlo)
  logical :: loor(nlo,0:lmax),rlo(nlo,0:lmax),lapw(0:lmax),rel
  integer :: ilo(0:lmax),i,mode,l
  r0=1e-4_dp;dx=log(2.1_dp/r0)/(n-1);energy=.17_dp;charge=26;largest=0
  do i=1,n
     r(i)=r0*exp(dx*(i-1))
     potential(i)=-charge*exp(-r(i))+.15_dp*r(i)**2
  enddo
  loor=.false.;rlo=.false.;lapw=.true.;ilo=0;el=2*energy;elo=2*energy
  do mode=0,1
     rel=mode==1
     call ATPARN(u,du,b,db,up,dup,bp,dbp, &
          alo_r,blo_r,aplo,bplo,boundary,slope,overlap,alo,blo,clo, &
          2*potential,n,el,elo,loor,rlo,lapw,0,0,0,ilo,rel,charge, &
          .false.,n,r0,dx,r(n),lmax,lmax,nlo,n)
     do l=0,lmax
        call prepare_radial_pair(ctx,rel,potential,r0,dx,energy,real(l,dp),charge)
        values=radial_pair_value(ctx)
        error=max(maxval(abs(values(:,1,1)-u(:,l))),maxval(abs(values(:,2,1)-b(:,l))), &
                  maxval(abs(values(:,3,1)-up(:,l))),maxval(abs(values(:,4,1)-bp(:,l))), &
                  maxval(abs(values(:,1,2)-du(:,l))),maxval(abs(values(:,2,2)-db(:,l))), &
                  maxval(abs(values(:,3,2)-dup(:,l))),maxval(abs(values(:,4,2)-dbp(:,l))))
        largest=max(largest,error)
        write(*,'(A,L1,A,I2,A,ES15.7)') 'RELA=',rel,' l=',l,' ATPARN max abs difference=',error
     enddo
  enddo
  if(largest>2e-10_dp) error stop 'ATPAR central convention mismatch'
  write(*,'(A,ES15.7)') 'PASS central ATPARN closure=',largest
end program
