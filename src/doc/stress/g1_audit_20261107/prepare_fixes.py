from pathlib import Path
import re,shutil,hashlib,json,difflib
root=Path(__file__).resolve().parent
paths=[Path('/Users/haule/wien2k_24/SRC_lapw2/val_kinstrmt.F'),Path('/Users/haule/wien2k_24/SRC_lapw2/kin_str_integrals.F'),Path('/Users/haule/projects/eDMFT_new/src/dmft2/stress_kernels.f90'),Path('/Users/haule/projects/eDMFT_new/src/dmft2/stress_integrals.f90')]
comment='''    ! G1 product rule (audit 2026-11-07): F2 differentiates n_alpha only.
    ! The derivative of Y in D_beta(n_alpha*Y) belongs in F3 as well as
    ! the separate n_beta*D_alpha(Y) row. Thus F3 is
    !       (2*R1-R0) * (n_alpha*D_beta+n_beta*D_alpha)/2,
    ! R1=int r*u_bra*u'_ket dr, R0=int u_bra*u_ket dr.
    ! Halving the R1 coefficient was an incorrect earlier correction.
    ! In particular psi=r*Y_10 is proportional to z and has ZERO Hessian;
    ! the halved coefficient gives a spurious, traceless tensor instead.
    ! LO--LO products have no extra multiplicity: the ordered pair sum
    ! already contains every primitive pair, just as for u and udot.
    ! PREAUDIT explicitly reproduces the two old angular defects for
    ! diagnostics only; the physical/default scales are both one.
    g1_deriv_scale = 1.d0
    g1_lo_scale = 1.d0
    g1_legacy_env = ''
    call get_environment_variable('WIEN_STRESS_G1_PREAUDIT',g1_legacy_env, &
                                  status=g1_legacy_status)
    if(g1_legacy_status.eq.0 .and. len_trim(g1_legacy_env).gt.0 .and. &
       trim(g1_legacy_env).ne.'0') then
       g1_deriv_scale = 0.5d0
       g1_lo_scale = 2.d0
    endif
'''
helper='''
subroutine g1_radial_derivative(n,r,f,df)
  ! Stress-only derivative on the physical positive-r mesh.
  ! The old G1 path passed an UNINITIALIZED IsEven to DERGL, which
  ! actually uses the flag to mirror f across r=0. Even a specified
  ! parity is not generally valid for Coulomb-cusp/RELA radial functions,
  ! or for both a radial function and its derivative. No extrapolated
  ! negative-radius samples are needed: SPLINE already supplies a
  ! one-sided endpoint condition from positive-mesh divided differences.
  ! Keep the existing cubic interpolation and twice-differentiated
  ! strategy; this repair does NOT add physical radial response to strain.
  ! Generic DERGL and non-stress density/force routines are untouched.
  implicit none
  integer,intent(in) :: n
  real*8,intent(in) :: r(n),f(n)
  real*8,intent(out) :: df(n)
  real*8 :: c(n),d(n)
  if(n.lt.4) stop 'G1 radial derivative requires at least four points'
  call spline(n,r,f,df,c,d)
end subroutine g1_radial_derivative
'''
manifest={};patch=[]
for p in paths:
    kind='lapw2' if 'wien2k_24' in str(p) else 'dmft2'
    before=p.read_text();after=before
    target=root/'prepared'/kind/p.name;target.parent.mkdir(parents=True,exist_ok=True)
    backup=root/'before'/kind/p.name;backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,backup)
    if 'kinstrmt' in p.name or 'kernels' in p.name:
        start=after.index('    ! Eq. (G1)') if kind=='lapw2' else after.index('    ! Equation G1 carries')
        end=after.index('    DO l = 0, lmax',start)
        after=after[:start]+comment+'\n'+after[end:]
        after=after.replace(':: g1_deriv_scale',':: g1_deriv_scale, g1_lo_scale',1)
        if kind=='lapw2':
            for n in [1,2]:
                old=f'radial_buf{n} = radial_buf{n} + 2.d0* &'
                assert after.count(old)==1
                after=after.replace(old,f'radial_buf{n} = radial_buf{n} + g1_lo_scale* &')
        else:
            assert after.count('ohalf*g1_deriv_scale*gaunt1*2.d0* &')==1
            after=after.replace('ohalf*g1_deriv_scale*gaunt1*2.d0* &','ohalf*g1_deriv_scale*gaunt1*g1_lo_scale* &')
            old='2.d0*conjg( clm(index1,num,jlo) )*clm(index2,num,jlop)*&'
            assert after.count(old)==1
            after=after.replace(old,'g1_lo_scale*conjg( clm(index1,num,jlo) )*clm(index2,num,jlop)*&')
    else:
        match=re.search(r'(?ims)^subroutine kin_str_integrals\(.*?^end subroutine kin_str_integrals',after)
        assert match
        body=match.group()
        body=re.sub(r'(?im)^.*IsEven.*\n',lambda m:'' if not re.search(r'call dergl',m.group(),re.I) else m.group(),body)
        body,n=re.subn(r'(?i)call dergl(?:_stress)?\(([^\n]*?),(?:IsEven|\.true\.),g0\)',r'call g1_radial_derivative(\1)',body)
        assert n==24,(p,n)
        assert 'IsEven' not in body
        after=after[:match.start()]+body+after[match.end():]+helper
    target.write_text(after)
    manifest[str(p)]={'before_sha256':hashlib.sha256(before.encode()).hexdigest(),'after_sha256':hashlib.sha256(after.encode()).hexdigest(),'prepared':str(target),'backup':str(backup)}
    patch.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile=str(p)+'.before_g1_audit',tofile=str(p)))
(root/'g1-fixes.patch').write_text(''.join(patch));(root/'source_changes.json').write_text(json.dumps(manifest,indent=2)+'\n')
