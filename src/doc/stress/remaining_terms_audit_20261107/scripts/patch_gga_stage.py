from pathlib import Path
root=Path('/private/tmp/stress_remaining_audit_20261107/SRC_lapw0')
def change(name,old,new,count=1):
 p=root/name;s=p.read_text();assert s.count(old)==count,(name,s.count(old),old[:70]);p.write_text(s.replace(old,new))
# A small pure algebra kernel shared by real-space angular integration and
# the interstitial path. sig_* are the existing VXCLM2 (2/JSPIN)*Hartree
# derivatives; tensor is the positive metric contraction in Rydbergs.
(root/'gga_metric_tensor.f90').write_text('''! Metric derivative of a spin-resolved GGA energy density.
! With f in Ry/Bohr^3, define a_s=df/d|grad n_s|^2 and
! a_ud=df/d(grad n_up.grad n_dn).  The explicit gradient metric term
! in physical stress is -Q_ab/Omega, where
! Q_ab=2 a_up g_up,a g_up,b+2 a_dn g_dn,a g_dn,b
!      +a_ud(g_up,a g_dn,b+g_dn,a g_up,b).
! VXCLM2 returns sig=(2/JSPIN)*df_Ha/dsigma, so a=JSPIN*sig.
! In the unpolarized interface sig_up is instead df_Ry/d|grad n|^2.
! No occupations, density response, or extra density factor enter Q.
subroutine gga_metric_tensor(jspin,grad_up,grad_dn,sig_up,sig_dn,sig_ud,q)
  implicit none
  integer,intent(in) :: jspin
  real(8),intent(in) :: grad_up(3),grad_dn(3),sig_up,sig_dn,sig_ud
  real(8),intent(out) :: q(3,3)
  real(8) :: grad_total(3)
  integer :: a,b
  grad_total=grad_up+grad_dn
  do b=1,3
    do a=1,3
      if(jspin==1)then
        q(a,b)=2.d0*sig_up*grad_total(a)*grad_total(b)
      else
        q(a,b)=2.d0*dble(jspin)*(sig_up*grad_up(a)*grad_up(b) &
                              +sig_dn*grad_dn(a)*grad_dn(b)) &
              +dble(jspin)*sig_ud*(grad_up(a)*grad_dn(b)+grad_dn(a)*grad_up(b))
      endif
    enddo
  enddo
end subroutine gga_metric_tensor
''')
change('Makefile','kahan_summ_gxyz.o stress_gga.o','kahan_summ_gxyz.o stress_gga.o gga_metric_tensor.o')
change('modules.F','  real*8                  :: gga_fin_str(1:9)','''  real*8                  :: gga_fin_str(1:9)
! Full angular-grid integral of the gradient metric tensor, per atom type.
! Allocated only for a stress calculation; XCPOT1 writes disjoint atom slots.
  real(8),allocatable     :: gga_mt_metric(:,:)
  logical                :: gga_metric_full=.true.''')
change('lapw0.F','use gga_str_kb,  only: gga_int','use gga_str_kb,  only: gga_int, gga_mt_metric, gga_metric_full')
change('lapw0.F','      gga_int        = zeroc','''      ! Keep the former harmonic-diagonal GGA contraction available only
      ! for regression diagnostics. FULL evaluates the derived metric term.
      call init_gga_metric(nat)
      gga_int        = zeroc''')
# Kernel initialization lives in stress_gga so the existing module build order
# remains sufficient and no new module dependency is needed.
p=root/'stress_gga.F';s=p.read_text();s=s.replace('use gga_str_kb, only: gga_fin_str','use gga_str_kb, only: gga_fin_str, gga_mt_metric, gga_metric_full')
old='    DO 29 jatom = nstart(myid), nstop(myid)'
assert s.count(old)==1
new='''    if(gga_metric_full.and.allocated(gga_mt_metric))then
       ! Q was integrated on exactly the angular grid where VXCLM2 samples
       ! the energy density. Products of different harmonics and tangential
       ! density gradients are thus retained, without an intermediate
       ! projection/truncation of the nonlinear derivative df/dsigma.
       do jatom=nstart(myid),nstop(myid)
          do mu=1,mult(jatom)
             index_a=index_a+1
             buf_str4=gga_mt_metric(:,jatom)
             call symmetry_stress_rotij(index_a,jatom,buf_str4)
             gga_fin_str=gga_fin_str+buf_str4
          enddo
       enddo
#ifdef Parallel
       gga_fin_str1=gga_fin_str
       call mpi_reduce(gga_fin_str1,gga_fin_str,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif
       if(myid.eq.0)then
          write(6,'(a,9es24.15)') ':GGA_METRIC_MT_RY:',gga_fin_str
          write(6,'(a,9es24.15)') ':GGA_METRIC_INT_RY:',real(gga_str_int,8)
          gga_fin_str=gga_fin_str+real(gga_str_int,8)
       endif
       return
    endif
    DO 29 jatom = nstart(myid), nstop(myid)'''
s=s.replace(old,new)
s+='''
subroutine init_gga_metric(natom)
  use gga_str_kb, only: gga_mt_metric,gga_metric_full
  implicit none
  integer,intent(in) :: natom
  character(32) :: mode
  integer :: status
  if(allocated(gga_mt_metric))deallocate(gga_mt_metric)
  allocate(gga_mt_metric(9,natom));gga_mt_metric=0.d0
  mode=''
  call get_environment_variable('WIEN_STRESS_GGA_MODE',mode,status=status)
  gga_metric_full=trim(mode)/='LEGACY'
end subroutine init_gga_metric
''';p.write_text(s)
# XCPOT1 already has all derivatives in the local orthonormal spherical
# frame (GRANS); rotate that gradient before forming the Cartesian tensor.
change('xcpot1.f','      use str_switch, only : SW_STRESS, STRESS','''      use str_switch, only : SW_STRESS, STRESS
      use gga_str_kb, only: gga_mt_metric,gga_metric_full''')
change('xcpot1.f','       REAL*8 :: gxyzu1(nspd),gxyzd1(nspd),gxyzud1(nspd)','''       REAL*8 :: gxyzu1(nspd),gxyzd1(nspd),gxyzud1(nspd)
       real(8),allocatable :: metric_radial(:,:)
       real(8) :: er_st(3),et_st(3),ep_st(3),gu_st(3),gd_st(3),q_st(3,3)
       real(8) :: cp_st,sp_st,xy_st,metric_integral
       integer :: aa_st,bb_st,ii_st''')
change('xcpot1.f','       UseInterp= .false.','''       if(STRESS.and.GGA_SWITCH.and.gga_metric_full)then
          allocate(metric_radial(nrad,9));metric_radial=0.d0
       endif
       UseInterp= .false.''')
change('xcpot1.f','!$omp& sig_str_u,sig_str_d,sig_str_ud) &','''!$omp& sig_str_u,sig_str_d,sig_str_ud,er_st,et_st,ep_st,gu_st,gd_st, &
!$omp& q_st,cp_st,sp_st,xy_st,aa_st,bb_st,ii_st) &''')
anchor='''         IF(STRESS) THEN
           gxyzu1(k) = sig_str_u'''
block='''         if(STRESS.and.GGA_SWITCH.and.gga_metric_full)then
           ! GRANS returns components along (e_r,e_theta,e_phi), not the
           ! fixed Cartesian axes. Rotate before taking outer products.
           er_st=spt(:,k)
           xy_st=sqrt(er_st(1)**2+er_st(2)**2)
           cp_st=1.d0;sp_st=0.d0
           if(xy_st>1.d-14)then
             cp_st=er_st(1)/xy_st;sp_st=er_st(2)/xy_st
           endif
           et_st=[er_st(3)*cp_st,er_st(3)*sp_st,-xy_st]
           ep_st=[-sp_st,cp_st,0.d0]
           gu_st=gxu*er_st+gyu*et_st+gzu*ep_st
           gd_st=gxd*er_st+gyd*et_st+gzd*ep_st
           call gga_metric_tensor(jspin,gu_st,gd_st,sig_str_u,sig_str_d,sig_str_ud,q_st)
           ii_st=0
           do aa_st=1,3
             do bb_st=1,3
               ii_st=ii_st+1
               metric_radial(ir,ii_st)=metric_radial(ir,ii_st)+weight(k)*rp*rp*q_st(aa_st,bb_st)
             enddo
           enddo
         endif
         IF(STRESS) THEN
           gxyzu1(k) = sig_str_u'''
change('xcpot1.f',anchor,block)
change('xcpot1.f','      deallocate(spt,weight,vxcsu,vxcsd,vxu,vxd,exu,exd)','''      if(allocated(metric_radial))then
         ! The radial integrand already includes r^2 and the angular
         ! quadrature weights. Charge performs integral dr, not r^2 dr.
         do ii_st=1,9
           call charge(r,log(r(2)/r(1)),metric_radial(:,ii_st),1,jri(jatom),metric_integral)
           gga_mt_metric(ii_st,jatom)=metric_integral
         enddo
         deallocate(metric_radial)
      endif
      deallocate(spt,weight,vxcsu,vxcsd,vxu,vxd,exu,exd)''')
# Correct the spin-polarized interstitial contraction. Preserve the old
# return convention 2*n_up*up + 2*n_dn*dn used by XCPOT3, with identical
# two slots Q/(2*n_total). This avoids changing its FFT machinery.
p=root/'pwxad5.f';s=p.read_text();anchor='    real*8               :: grad_rho(3,2)';assert s.count(anchor)==1
s=s.replace(anchor,anchor+'\n    real(8) :: metric_q(3,3),metric_flat(9),density_total\n    logical :: metric_full\n    character(32) :: metric_mode\n    integer :: metric_status')
anchor='    grad_rho(1,1) = GXU; grad_rho(2,1) = GYU; grad_rho(3,1) = GZU';assert s.count(anchor)==1
# use module flag: avoid getenv per FFT point
s=s.replace('    logical :: metric_full\n    character(32) :: metric_mode\n    integer :: metric_status','')
start=s.index('SUBROUTINE GGA_INT_STRESS(')
pos=s.index('    IMPLICIT none',start);s=s[:pos]+s[pos:].replace('    IMPLICIT none','    use gga_str_kb, only: gga_metric_full\n    IMPLICIT none',1)
anchor='    grad_rho(1,2) = GXD; grad_rho(2,2) = GYD; grad_rho(3,2) = GZD'
assign=[]
for channel in ['up','dn']:
 for a in range(1,4):
  for b in range(1,4):assign.append(f'      gga_int_{channel}_{a}{b}=metric_q({a},{b})')
block=anchor+'''
    if(gga_metric_full)then
      call gga_metric_tensor(jspin,grad_rho(:,1),grad_rho(:,2),sig_str_u,sig_str_d,sig_str_ud,metric_q)
      density_total=du0+dd0
      if(density_total>1.d-30)then
        metric_q=metric_q/(2.d0*density_total)
      else
        metric_q=0.d0
      endif
'''+ '\n'.join(assign)+'''
      return
    endif
'''
assert s.count(anchor)==1;s=s.replace(anchor,block);p.write_text(s)
print('Staged full muffin-tin/spin-resolved GGA metric contraction')
