program probe_spherical
 use iso_fortran_env, only: real64
 use lapw_radial_response
 use lapw_radial_hamiltonian
 implicit none
 integer,parameter::dp=real64
 type(radial_basis_context)::ctx
 integer::n,nlo,mode,ell,lapw,iu,ou,np,m
 real(dp)::r0,dx,rmt,e,z,de
 real(dp),allocatable::v(:),dv(:),b(:,:),h(:,:),o(:,:),zero(:),el(:),del(:)
 character(1024)::input,output
 call get_command_argument(1,input);call get_command_argument(2,output)
 open(newunit=iu,file=trim(input),status='old')
 read(iu,*)n,nlo,mode,ell,lapw,np
 read(iu,*)r0,dx,rmt,e,z,de
 m=np+nlo
 allocate(v(n),dv(n),b(2,np),h(m,m),o(m,m),zero(n),el(nlo),del(nlo));zero=0
 if(nlo>0)then
  read(iu,*)el;read(iu,*)del
 endif
 read(iu,*)v;read(iu,*)dv;read(iu,*)b;close(iu)
 call prepare_radial_basis(ctx,mode==1,v,r0,dx,rmt,e,el,lapw==1,real(ell,dp),z)
 call radial_hamiltonian_action(ctx,b,e,el,zero,dv,de,del,-dv,h,o)
 open(newunit=ou,file=trim(output),access='stream',form='unformatted',status='replace')
 write(ou)h,o;close(ou)
end program
