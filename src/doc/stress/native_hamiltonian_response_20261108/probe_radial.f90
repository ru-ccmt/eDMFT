program probe_radial
 use iso_fortran_env, only: real64
 use lapw_radial_response
 use lapw_radial_hamiltonian
 implicit none
 integer,parameter::dp=real64
 type(radial_basis_context)::ctx
 integer::n,np,mode,ell,lapw,iu,ou
 real(dp)::r0,dx,rmt,e,z
 real(dp),allocatable::v(:),dv(:),b(:,:),h(:,:),o(:,:),zero(:),el(:),de(:)
 character(1024)::input,output
 call get_command_argument(1,input);call get_command_argument(2,output)
 open(newunit=iu,file=trim(input),status='old')
 read(iu,*)n,np,mode,ell,lapw
 read(iu,*)r0,dx,rmt,e,z
 allocate(v(n),dv(n),b(2,np),h(np,np),o(np,np),zero(n),el(0),de(0));zero=0
 read(iu,*)v;read(iu,*)dv;read(iu,*)b;close(iu)
 call prepare_radial_basis(ctx,mode==1,v,r0,dx,rmt,e,el,lapw==1,real(ell,dp),z)
 call radial_hamiltonian_action(ctx,b,e,el,zero,dv,0._dp,de,-dv,h,o)
 open(newunit=ou,file=trim(output),access='stream',form='unformatted',status='replace')
 write(ou)h,o;close(ou)
end program
