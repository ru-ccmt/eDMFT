program export_radial
 use iso_fortran_env, only: real64
 use lapw_radial_response
 implicit none
 integer,parameter::dp=real64
 type(radial_basis_context)::ctx
 integer::n,nlo,mode,ell,lapw,iu,ou
 real(dp)::r0,dx,rmt,e,z,de
 real(dp),allocatable::v(:),dv(:),b(:,:),u(:,:,:),dm(:,:),el(:),del(:)
 character(1024)::input,output
 call get_command_argument(1,input);call get_command_argument(2,output)
 open(newunit=iu,file=trim(input),status='old')
 read(iu,*)n,nlo,mode,ell,lapw
 read(iu,*)r0,dx,rmt,e,z,de
 allocate(v(n),dv(n),b(2,0),u(n,4,2+nlo),dm(2+nlo,nlo),el(nlo),del(nlo))
 if(nlo>0)then
  read(iu,*)el;read(iu,*)del
 endif
 read(iu,*)v;read(iu,*)dv;close(iu)
 call prepare_radial_basis(ctx,mode==1,v,r0,dx,rmt,e,el,lapw==1,real(ell,dp),z)
 call radial_components_action(ctx,b,dv,de,del,u,dm)
 open(newunit=ou,file=trim(output),access='stream',form='unformatted',status='replace')
 write(ou)ctx%primitive,u,ctx%overlap_metric,ctx%projector_metric,radial_component_coefficients(ctx,b),dm
 close(ou)
end program
