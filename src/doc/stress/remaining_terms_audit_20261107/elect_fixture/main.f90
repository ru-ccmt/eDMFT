program check
 use struct
 use densit
 use work
 implicit none
 integer :: lm(2,264,1),lmmax(1),mode,side,i,l,m
 real(8) :: r(881),xc(9),pi
 complex(8) :: inter(9)
 pi=acos(-1d0);dx=log(1d5)/880
 do i=1,881
 r(i)=r0(1)*exp((i-1)*dx(1))
 enddo
 rotloc=0;rotij=0
 do i=1,3
 rotloc(i,i,1)=1;rotij(i,i,1)=1
 enddo
 do side=1,3
 do mode=0,5
  lm=0;clm=0;v=0;exc=0;vxc=0;xc=0;inter=0
  clm(:,1,1)=4*pi*r*r
  v(:,1,1)=sqrt(4*pi)
  lmmax=1
  if(mode>0)then
   select case(mode)
   case(1);l=2;m=0
   case(2);l=2;m=1
   case(3);l=2;m=2
   case(4);l=-2;m=1
   case(5);l=-2;m=2
   end select
   lm(:,2,1)=[l,m];lm(:,3,1)=[l,-m];lmmax=2
   if(side==1.or.side==3)then
    v=0;v(:,2,1)=r*r
   endif
   if(side==2.or.side==3)then
    clm=0;clm(:,2,1)=r**4
   endif
  endif
  write(21,'(a,2i4)') 'CASE',side,mode
  call finl_elect_str(lm,lmmax,inter,0d0,xc)
 enddo
 enddo
end program
