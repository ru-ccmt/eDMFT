program radial_pullback
 use iso_fortran_env, only: real64
 use lapw_radial_response
 implicit none
 integer,parameter::dp=real64
 type(radial_basis_context)::ctx
 integer::n,nlo,mode,ell,lapw,iu,ou
 real(dp)::r0,dx,rmt,e,z,de,eb
 real(dp),allocatable::v(:),dv(:),b(:,:),gp(:,:,:),gm(:,:),el(:),del(:),vb(:),elb(:)
 character(1024)::input,covector,output
 call get_command_argument(1,input);call get_command_argument(2,covector);call get_command_argument(3,output)
 open(newunit=iu,file=trim(input),status='old')
 read(iu,*)n,nlo,mode,ell,lapw
 read(iu,*)r0,dx,rmt,e,z,de
 allocate(v(n),dv(n),b(2,0),gp(n,4,2+nlo),gm(2+nlo,nlo),el(nlo),del(nlo),vb(n),elb(nlo))
 if(nlo>0)then
  read(iu,*)el;read(iu,*)del
 endif
 read(iu,*)v;read(iu,*)dv;close(iu)
 call prepare_radial_basis(ctx,mode==1,v,r0,dx,rmt,e,el,lapw==1,real(ell,dp),z)
 open(newunit=iu,file=trim(covector),access='stream',form='unformatted',status='old')
 read(iu)gp,gm;close(iu)
 call radial_components_pullback(ctx,b,gp,gm,vb,eb,elb)
 open(newunit=ou,file=trim(output),access='stream',form='unformatted',status='replace')
 write(ou)vb,eb,elb;close(ou)
end program
