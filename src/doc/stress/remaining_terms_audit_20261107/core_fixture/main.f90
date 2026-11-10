program check_core
 use struct
 use core_nsp
 use core_sp
 implicit none
 real(8) :: r(1022),rho(1022),dx
 integer :: i,mode,l,m,ia
 allocate(lm11_st(2,264,1),lmmax_st(1),vns_st(1022,264,1))
 dx=log(1.d5)/880.d0
 rotij=0;rotloc=0
 do i=1,3
  rotloc(i,i,1)=1;rotij(i,i,:)=1
 enddo
 do i=1,881
  r(i)=1.d-5*exp((i-1)*dx);rho(i)=r(i)**2
 enddo
 do mode=0,10
  lm11_st=0;lmmax_st=0;vns_st=0;dv1p=0;sum_tot=0
  if(mode==0)then
   dv1p(1:881)=r(1:881)**2/2
   l=0;m=0
  else
   select case(mode)
   case(1);l=2;m=0
   case(2);l=2;m=1
   case(3);l=2;m=2
   case(4);l=-2;m=1
   case(5);l=-2;m=2
   case(6);l=1;m=0
   case(7);l=3;m=0
   case(8);l=4;m=0
   case(9);l=4;m=4
   case(10);l=-4;m=2
   end select
   lmmax_st=1;lm11_st(:,1,1)=[l,m]
   vns_st(1:881,1,1)=r(1:881)**abs(l)
  endif
  ia=0
  call core_corr_stress(1,881,r,rho,1,ia)
  write(*,'(a,3i4,9es24.15)') 'FIXTURE',mode,l,m,sum_tot(:,1)
 enddo
end program
