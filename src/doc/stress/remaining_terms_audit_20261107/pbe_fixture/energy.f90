subroutine test_energy(nu,nd,gu,gd,energy,q)
 implicit none
 real(8)::nu,nd,gu(3),gd(3),energy,q(3,3),u,v,g,exu,exd,vxu,vxd,ec,vcu,vcd,su,sd,sud
 logical::cu,cd
 u=sqrt(sum(gu*gu));v=sqrt(sum(gd*gd));g=sqrt(sum((gu+gd)**2))
 call pbex(nu,u,0d0,0d0,nd,v,0d0,0d0,exu,exd,vxu,vxd,1d-18,cu,cd)
 call pbec(nu,u,0d0,nd,v,0d0,g,0d0,ec,vcu,vcd,cu,cd)
 energy=2d0*(nu*exu+nd*exd+(nu+nd)*ec)
 call pbe_vsigma_builtin(nu,nd,u,v,g,.2195149727645171d0,.804d0,.06672455060314922d0,1d0,2,cu,cd,su,sd,sud)
 call gga_metric_tensor(2,gu,gd,su,sd,sud,q)
end subroutine
