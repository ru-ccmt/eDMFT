!BOP
!
! !ROUTINE: momradintc
!
! !INTERFACE:
subroutine momradintc(rel,iat,nln,is,iso)
  !
  ! !USES:
  use comi, nat=>natti		!LO ncg
  use core    
  use lolog1,   only: loor1
  use moments 
  use potnlc	
  use radfun
  INCLUDE 'param.inc'		!LO
  ! !INPUT PARAMETERS:
  integer(4), intent(in) :: iat  ! Runs over inequivalent atoms      
  integer(4), intent(in) :: nln    ! azimuthal quantum number of core state
  logical  rel
  integer(4), intent(in) :: is   
  ! !LOCAL VARIABLES:
  integer(4) :: k,i
  integer(4) :: jri1  ! number of radial mesh points
  integer(4) :: m,m0,npo2
  real(8) :: fl      ! angular momentum quantum number      
  real(8) :: iuup    ! \int(u_l du_l'/dr r^2 dr)      
  real(8) :: iuu    ! \int(u_l u_l' r dr)      
  real(8), dimension(4) :: c
  !LO-----------------------------
  !LO--------------------------------------NEW RADIAL FUNCTION--------------------
  real(8), allocatable :: ucl(:)
  real(8), allocatable :: uscl(:)
  real(8), allocatable :: ul(:)
  real(8), allocatable :: usl(:)
  real(8), allocatable :: upl(:)
  real(8), allocatable :: uspl(:)
  real(8), allocatable :: udl(:)
  real(8), allocatable :: usdl(:)
  real(8), allocatable :: udor(:)
  real(8), allocatable :: usdor(:)
  real(8), allocatable :: ulol(:)
  real(8), allocatable :: uslol(:)
  real(8), allocatable :: ulopl(:)
  real(8), allocatable :: uslopl(:)
  real(8), allocatable :: ulr(:)
  real(8), allocatable :: uslr(:)
  real(8), allocatable :: udlr(:)
  real(8), allocatable :: usdlr(:)
  real(8), allocatable :: ulolr(:)
  real(8), allocatable :: uslolr(:)
  real(8), allocatable :: ulolp(:)
  real(8), allocatable :: uslolp(:)
  real(8), allocatable :: udlp(:)
  real(8), allocatable :: usdlp(:)
  real(8), allocatable :: ulolrp(:)
  real(8), allocatable :: uslolrp(:)
  real(8), allocatable :: ulp(:)
  real(8), allocatable :: uslp(:)
  real(8), allocatable :: ulrp(:)
  real(8), allocatable :: uslrp(:)
  real(8), allocatable :: ucor(:)
  real(8), allocatable :: udlrp(:)
  real(8), allocatable :: usdlrp(:)
  real(8), allocatable :: uscor(:)
  ! !DEFINED PARAMETERS:
  integer(4), parameter :: np = 4 
  real(8), parameter :: two = 2.0d+0      
  real(8), parameter :: one = 1.0d+0      
  real(8), parameter :: cf2 = 1.0d-22      
  ! !EXTERNAL ROUTINES: 
  external rint13
  ! !INTRINSIC ROUTINES:  
  intrinsic iabs
  intrinsic isign
  !     Allocate necesary arrays
  !     
  jri1=jri(iat)
  allocate(ucl(jri1))
  allocate(uscl(jri1))
  allocate(ucor(jri1))
  allocate(uscor(jri1))
  allocate(ul(jri1))
  allocate(usl(jri1))
  allocate(upl(jri1))
  allocate(uspl(jri1))
  allocate(udl(jri1))
  allocate(usdl(jri1))
  allocate(udor(jri1))
  allocate(usdor(jri1))
  allocate(ulol(jri1))
  allocate(uslol(jri1))
  allocate(ulopl(jri1))
  allocate(uslopl(jri1))
  allocate(ulr(jri1))
  allocate(uslr(jri1))
  allocate(udlr(jri1))
  allocate(usdlr(jri1))
  allocate(ulolr(jri1))
  allocate(uslolr(jri1))
  allocate(ulolp(jri1))
  allocate(uslolp(jri1))
  allocate(udlp(jri1))
  allocate(usdlp(jri1))
  allocate(ulolrp(jri1))
  allocate(uslolrp(jri1))
  allocate(ulp(jri1))
  allocate(uslp(jri1))
  allocate(ulrp(jri1))
  allocate(uslrp(jri1))
  allocate(udlrp(jri1))
  allocate(usdlrp(jri1))
  !LO
  !LO_____Inizialize core integrals________________________________
  !
  ncormax = ncg(iat)
  npo2=np/2
  fl=dble(nln)
  !=========================================================================
  !               begin the calculation of iul1ul like integrals           !
  !=========================================================================
  !
  !      store the needed shared arrays in the local ones
  !?????
  ucl(1:jri1)=ucore_1(1:jri1,is)
  uscl(1:jri1)=uscore_1(1:jri1,is)
  ucor(1:jri1)=ucore_1(1:jri1,is)/rr(1:jri1,iat)
  uscor(1:jri1)=uscore_1(1:jri1,is)/rr(1:jri1,iat)	
  !
  if(nln.gt.0)then
     ulr(1:jri1)=rrad01(1:jri1,nln-1)/rr(1:jri1,iat)     ! u_A/r
     uslr(1:jri1)=rrad02(1:jri1,nln-1)/rr(1:jri1,iat)	 ! u_B/r/137
     ul(1:jri1)=rrad01(1:jri1,nln-1)                     ! u_A
     usl(1:jri1)=rrad02(1:jri1,nln-1)                    ! u_B/137
     udlr(1:jri1)=rade01(1:jri1,nln-1)/rr(1:jri1,iat)    ! udot_A/r
     usdlr(1:jri1)=rade02(1:jri1,nln-1)/rr(1:jri1,iat)   ! udor_B/r /137
     udl(1:jri1)=rade01(1:jri1,nln-1)                    ! udot_A
     usdl(1:jri1)=rade02(1:jri1,nln-1)                   ! udot_B/137
     upl(1)=0.0d0
     do m=2,jri1
        if(m.le.npo2)then
           m0=1
        elseif(m.gt.jri1-npo2)then
           m0=jri1-np+1
        else
           m0=m-npo2
        endif
        upl(m)=polynom(1,np,rr(m0,iat),ul(m0),c,rr(m,iat))-ulr(m)  ! upl = d(u_A)/dr-u_A/r
     enddo !m
     uspl(1)=0.0d0
     do m=2,jri1
        if(m.le.npo2)then
           m0=1
        elseif(m.gt.jri1-npo2)then
           m0=jri1-np+1
        else
           m0=m-npo2
        endif
        uspl(m)=polynom(1,np,rr(m0,iat),usl(m0),c,rr(m,iat))-uslr(m) ! uspl = (d(u_B)/dr-u_B/r)/137.
     enddo !m
     udor(1)=0.0d0
     do m=2,jri1
        if(m.le.npo2)then
           m0=1
        elseif(m.gt.jri1-npo2)then
           m0=jri1-np+1
        else
           m0=m-npo2
        endif
        udor(m)=polynom(1,np,rr(m0,iat),udl(m0),c,rr(m,iat))-udlr(m) ! udor = d(udot_A)/dr-udot_A/r
     enddo !m
     usdor(1)=0.0d0
     do m=2,jri1
        if(m.le.npo2)then
           m0=1
        elseif(m.gt.jri1-npo2)then
           m0=jri1-np+1
        else
           m0=m-npo2
        endif
        usdor(m)=polynom(1,np,rr(m0,iat),usdl(m0),c,rr(m,iat))-usdlr(m) ! usdor = (d(udot_B)/dr-udot_B/r)/137.
     enddo !m
     !        
     call rint13(rel,ucl,uscl,upl,uspl,iuup,iat) ! iuup = <u_Acore*r| d(u_A)/dr-u_A/r > + <u_Bcore*r| d(u_B)/dr-u_B/r > /137^2
     call rint13(rel,ucor,uscor,ul,usl,iuu,iat)  ! iuu  = <u_Acore  | u_A> + <u_Bcore|u_B>/137^2    
     iucl1ul(iat,is)=iuup-(fl-one)*iuu           ! iucl1ul = iuup - (l-1)*iuu
     call rint13(rel,ucl,uscl,udor,usdor,iuup,iat)! iuup = <u_Acore*r| d(udot_A)/dr-udot_A/r> + <u_Bcore*r|(d(udot_B)/dr-udot_B/r)>/137^2
     call rint13(rel,ucor,uscor,udl,usdl,iuu,iat) ! iuu = <u_Acore|udot_A> +  <u_Bcore|udot_B>/137^2
     iucl1udl(iat,is)=iuup-(fl-one)*iuu           ! iucl1udl = <u_Acore*r|d(udot_A)/dr-udot_A/r> + <u_Bcore*r|(d(udot_B)/dr-udot_B/r)>/137^2 - (l-1)*(<u_Acore|udot_A> +  <u_Bcore|udot_B>/137^2)
     
     !------------------------------------------------------        
     !       And now for local orbitals, l< lomax
     !------------------------------------------------------        
     if(nln.le.lomax+1)then
        if(loor1(nln-1))then
           ulol(1:jri1)=(a01lo(1:jri1,nln-1))
           uslol(1:jri1)=b01lo(1:jri1,nln-1)
           ulolr(1:jri1)=(a01lo(1:jri1,nln-1))/rr(1:jri1,iat)
           uslolr(1:jri1)=(b01lo(1:jri1,nln-1))/rr(1:jri1,iat)
           ulopl(1)=0.0d0
           do m=2,jri1
              if(m.le.npo2)then
                 m0=1
              elseif(m.gt.jri1-npo2)then
                 m0=jri1-np+1
              else
                 m0=m-npo2
              endif
              ulopl(m)=polynom(1,np,rr(m0,iat),ulol(m0),c,rr(m,iat))-ulolr(m)
           enddo !m
           
           uslopl(1)=0.0d0
           do m=2,jri1
              if(m.le.npo2)then
                 m0=1
              elseif(m.gt.jri1-npo2)then
                 m0=jri1-np+1
              else
                 m0=m-npo2
              endif
              uslopl(m)=polynom(1,np,rr(m0,iat),uslol(m0),c,rr(m,iat))-uslolr(m)
           enddo !m
           !            
           call rint13(rel,ucl,uscl,ulopl,uslopl,iuup,iat)
           call rint13(rel,ucor,uscor,ulol,uslol,iuu,iat)         
           iucl1ulol(iat,is)=iuup-(fl-one)*iuu
        endif ! loor 
     endif ! l <= lomax  
  endif ! l > 0
  !=========================================================================
  !                 end the calculation of iul1ul like integrals           !
  !=========================================================================
  
  !=========================================================================
  !               begin the calculation of iulul1 like integrals           !
  !=========================================================================
  !
  !      store the needed shared arrays in the local ones
  !       
  ulrp(1:jri1)=rrad01(1:jri1,nln+1)/rr(1:jri1,iat)
  uslrp(1:jri1)=rrad02(1:jri1,nln+1)/rr(1:jri1,iat)	
  ulp(1:jri1)=rrad01(1:jri1,nln+1)
  uslp(1:jri1)=rrad02(1:jri1,nln+1)
  udlrp(1:jri1)=rade01(1:jri1,nln+1)/rr(1:jri1,iat)
  usdlrp(1:jri1)=rade02(1:jri1,nln+1)/rr(1:jri1,iat)
  ul(1:jri1)=rrad01(1:jri1,nln+1)
  usl(1:jri1)=rrad02(1:jri1,nln+1)
  udl(1:jri1)=rade01(1:jri1,nln+1)
  usdl(1:jri1)=rade02(1:jri1,nln+1)
  udlp(1:jri1)=rade01(1:jri1,nln+1)
  usdlp(1:jri1)=rade02(1:jri1,nln+1)
  
  upl(1)=0.0d0	  
  do m=2,jri1
     if(m.le.npo2)then
        m0=1
     elseif(m.gt.jri1-npo2)then
        m0=jri1-np+1
     else
        m0=m-npo2
     endif
     upl(m)=polynom(1,np,rr(m0,iat),ulp(m0),c,rr(m,iat))-ulrp(m)
  enddo !m
  
  uspl(1)=0.0d0
  do m=2,jri1
     if(m.le.npo2)then
        m0=1
     elseif(m.gt.jri1-npo2)then
        m0=jri1-np+1
     else
        m0=m-npo2
     endif
     uspl(m)=polynom(1,np,rr(m0,iat),uslp(m0),c,rr(m,iat))-uslrp(m)
  enddo !m
  
  udor(1)=0.0d0
  do m=2,jri1
     if(m.le.npo2)then
        m0=1
     elseif(m.gt.jri1-npo2)then
        m0=jri1-np+1
     else
        m0=m-npo2
     endif
     udor(m)=polynom(1,np,rr(m0,iat),udlp(m0),c,rr(m,iat))-udlrp(m)
  enddo !m!
  
  usdor(1)=0.0d0
  do m=2,jri1
     if(m.le.npo2)then
        m0=1
     elseif(m.gt.jri1-npo2)then
        m0=jri1-np+1
     else
        m0=m-npo2
     endif
     usdor(m)=polynom(1,np,rr(m0,iat),usdlp(m0),c,rr(m,iat))-usdlrp(m)
  enddo !m 
  !        
  call rint13(rel,ucl,uscl,upl,uspl,iuup,iat)
  call rint13(rel,ucor,uscor,ul,usl,iuu,iat)
  iuclul1(iat,is)=iuup+(fl+two)*iuu
  call rint13(rel,ucl,uscl,udor,usdor,iuup,iat)
  call rint13(rel,ucor,uscor,udl,usdl,iuu,iat)
  iucludl1(iat,is)=iuup+(fl+two)*iuu
  
  !------------------------------------------------------        
  !       And now for local orbitals, l< lomax
  !------------------------------------------------------        
  if(nln.le.lomax-1)then
     if(loor1(nln+1))then
        ulolp(1:jri1)=(a01lo(1:jri1,nln+1))
        uslolp(1:jri1)=b01lo(1:jri1,nln+1)
        ulolrp(1:jri1)=(a01lo(1:jri1,nln+1))/rr(1:jri1,iat)
        uslolrp(1:jri1)=(b01lo(1:jri1,nln+1))/rr(1:jri1,iat)
	ulopl(1)=0.0d0
        do m=2,jri1
           if(m.le.npo2)then
              m0=1
           elseif(m.gt.jri1-npo2)then
              m0=jri1-np+1
           else
              m0=m-npo2
           endif
           ulopl(m)=polynom(1,np,rr(m0,iat),ulolp(m0),c,rr(m,iat))-ulolrp(m)
        enddo !m
	uslopl(1)=0.0d0
        do m=2,jri1
           if(m.le.npo2)then
              m0=1
           elseif(m.gt.jri1-npo2)then
              m0=jri1-np+1
           else
              m0=m-npo2
           endif
           uslopl(m)=polynom(1,np,rr(m0,iat),uslolp(m0),c,rr(m,iat))-uslolrp(m)
        enddo !m
        call rint13(rel,ucl,uscl,ulopl,uslopl,iuup,iat)
        call rint13(rel,ucor,uscor,ulolp,uslolp,iuu,iat)
        iuclulol1(iat,is)=iuup+(fl+two)*iuu
     endif ! loor
  endif ! l < lomax  
  !     Deallocate local arrays
  !
  deallocate(ucl)
  deallocate(uscl)
  deallocate(ucor)
  deallocate(uscor)
  deallocate(ul)
  deallocate(usl)
  deallocate(upl)
  deallocate(uspl)
  deallocate(udl)
  deallocate(usdl)
  deallocate(udor)
  deallocate(usdor)
  deallocate(ulol)
  deallocate(uslol)
  deallocate(ulopl)
  deallocate(uslopl)
  deallocate(ulr)
  deallocate(uslr)
  deallocate(udlr)
  deallocate(usdlr)
  deallocate(ulolr)
  deallocate(uslolr)
  deallocate(ulolp)
  deallocate(uslolp)
  deallocate(udlp)
  deallocate(usdlp)
  deallocate(ulolrp)
  deallocate(uslolrp)
  deallocate(ulp)
  deallocate(uslp)
  deallocate(ulrp)
  deallocate(uslrp)
  deallocate(udlrp)
  deallocate(usdlrp)
  
  !LO________POLYNOM FUNCTION ADDED________________
  ! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
  ! This file is distributed under the terms of the GNU Lesser General Public
  ! License. See the file COPYING for license details.
  !BOP
  ! !ROUTINE: polynom
  ! !INTERFACE:
contains
  real(8) function polynom(m,np,xa,ya,c,x)
    ! !INPUT/OUTPUT PARAMETERS:
    !   m  : order of derivative (in,integer)
    !   np : number of points to fit (in,integer)
    !   xa : abscissa array (in,real(np))
    !   ya : ordinate array (in,real(np))
    !   c  : work array (out,real(np))
    !   x  : evaluation abscissa (in,real)
    ! !DESCRIPTION:
    !   Fits a polynomial of order $n_p-1$ to a set of $n_p$ points. If $m\ge 0$ the
    !   function returns the $m$th derviative of the polynomial at $x$, while for
    !   $m<0$ the integral of the polynomial from the first point in the array to
    !   $x$ is returned.
    !
    ! !REVISION HISTORY:
    !   Created October 2002 (JKD)
    !EOP
    !BOC
    implicit none
    ! argmuments
    integer, intent(in) :: m
    integer, intent(in) :: np
    real(8), intent(in) :: xa(np)
    real(8), intent(in) :: ya(np)
    real(8), intent(out) :: c(np)
    real(8), intent(in) :: x
    ! local variables
    integer i,j,k
    real(8) x0,x1,x2,x3,y1,y2,y3
    real(8) t1,t2,t3,t4,t5,t6,t7,sum
    ! fast evaluations for small np
    select case(np)
    case(1)
       select case(m)
       case(:-1)
          polynom=ya(1)*(x-xa(1))
       case(0)
          polynom=ya(1)
       case default
          polynom=0.d0
       end select
       return
    case(2)
       c(2)=(ya(2)-ya(1))/(xa(2)-xa(1))
       t1=x-xa(1)
       select case(m)
       case(:-1)
          polynom=0.5d0*c(2)*t1**2+ya(1)*t1
       case(0)
          polynom=c(2)*t1+ya(1)
       case(1)
          polynom=c(2)
       case default
          polynom=0.d0
       end select
       return
    case(3)
       x1=xa(2)-xa(1)
       x2=xa(3)-xa(1)
       y1=ya(2)-ya(1)
       y2=ya(3)-ya(1)
       t1=1.d0/(x1*x2*(x2-x1))
       t2=x1*y2
       t3=x2*y1
       c(2)=t1*(x2*t3-x1*t2)
       c(3)=t1*(t2-t3)
       t1=x-xa(1)
       select case(m)
       case(:-1)
          polynom=(1.d0/3.d0)*c(3)*t1**3+0.5d0*c(2)*t1**2+ya(1)*t1
       case(0)
          polynom=c(3)*t1**2+c(2)*t1+ya(1)
       case(1)
          polynom=2.d0*c(3)*t1+c(2)
       case(2)
          polynom=2.d0*c(3)
       case default
          polynom=0.d0
       end select
       return
    case(4)
       x1=xa(2)-xa(1)
       x2=xa(3)-xa(1)
       x3=xa(4)-xa(1)
       y1=ya(2)-ya(1)
       y2=ya(3)-ya(1)
       y3=ya(4)-ya(1)
       t1=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
       t2=x1*x2*y3
       t3=x2*x3*y1
       t4=x3*x1*y2
       t5=x1**2
       t6=x2**2
       t7=x3**2
       c(2)=t1*(t3*(x3*t6-x2*t7)+t4*(x1*t7-x3*t5)+t2*(x2*t5-x1*t6))
       c(3)=t1*(t3*(t7-t6)+t4*(t5-t7)+t2*(t6-t5))
       c(4)=t1*(t3*(x2-x3)+t4*(x3-x1)+t2*(x1-x2))
       t1=x-xa(1)
       t2=t1**2
       select case(m)
       case(:-1)
          polynom=0.25d0*c(4)*t2**2+(1.d0/3.d0)*c(3)*t1*t2+0.5d0*c(2)*t2+ya(1)*t1
       case(0)
          polynom=c(4)*t1*t2+c(3)*t2+c(2)*t1+ya(1)
       case(1)
          polynom=3.d0*c(4)*t2+2.d0*c(3)*t1+c(2)
       case(2)
          polynom=6.d0*c(4)*t1+2.d0*c(3)
       case(3)
          polynom=6.d0*c(4)
       case default
          polynom=0.d0
       end select
       return
    end select
    if (np.le.0) then
       write(*,*)
       write(*,'("Error(polynom): np <= 0 : ",I8)') np
       write(*,*)
       stop
    end if
    if (m.ge.np) then
       polynom=0.d0
       return
    end if
    ! find the polynomial coefficients in divided differences form
    c(:)=ya(:)
    do i=2,np
       do j=np,i,-1
          c(j)=(c(j)-c(j-1))/(xa(j)-xa(j+1-i))
       end do
    end do
    ! special case m=0
    if (m.eq.0) then
       sum=c(1)
       t1=1.d0
       do i=2,np
          t1=t1*(x-xa(i-1))
          sum=sum+c(i)*t1
       end do
       polynom=sum
       return
    end if
    x0=xa(1)
    ! convert to standard form
    do j=1,np-1
       do i=1,np-j
          k=np-i
          c(k)=c(k)-(xa(k-j+1)-x0)*c(k+1)
       end do
    end do
    if (m.gt.0) then
       ! take the m'th derivative
       do j=1,m
          do i=m+1,np
             c(i)=c(i)*dble(i-j)
          end do
       end do
       t1=c(np)
       t2=x-x0
       do i=np-1,m+1,-1
          t1=t1*t2+c(i)
       end do
       polynom=t1
    else
       ! find the integral
       t1=c(np)/dble(np)
       t2=x-x0
       do i=np-1,1,-1
          t1=t1*t2+c(i)/dble(i)
       end do
       polynom=t1*t2
    end if
    return
  end function polynom
  !LO
  !LO
end subroutine momradintc
!EOC      
