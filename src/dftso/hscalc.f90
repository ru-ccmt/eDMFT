real*8 function hscalc(f_l1,f_l2,vrr,eev,l,mesh,dx,rnot)
  USE param
  implicit none
  real*8 vrr(nrad),eev
  integer l 
  real*8 res
  integer mesh
  real*8 dx,rnot
  
  real*8 f_l1(nrad),f_l2(nrad)      
  real*8 mfield(nrad),rmesh(nrad),fint(nrad)
  real*8 dv(nrad)
  integer i
  
  do i = 1,mesh
     rmesh(i)  = rnot*exp((i-1)*dx)
     mfield(i) = 0.5d0+1.d0/(2.d0*(2.d0*clight)**2)*(eev-vrr(i))  
  enddo

  call dergl2(vrr(1),dv(1),rnot,dx,mesh)

  do i=1,mesh
     fint(i)=dv(i)*f_l1(i)*f_l2(i)/(4.d0*mfield(i)**2*(2.d0*clight)**2)*2.d0/rmesh(i) 
  enddo
  
  call cali(fint,res,rmesh(1),dx,mesh) 
  hscalc=res
  return
end function hscalc
      
