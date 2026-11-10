from pathlib import Path
import difflib

root = Path(__file__).resolve().parent
p = root / 'SRC_lapw2/val_corr_stress.F'
old = p.read_text()
new = old.replace(
    '(elo(lp,jlop)-e(num))*jr3upu_22(jlo,jlop,l,lp)',
    '''! G2 differentiates the BRA radial function.  The integral
                                                 ! table stores <i|r d_r|j>, so transpose BOTH
                                                 ! LO and angular indices, while retaining the
                                                 ! KET linearization energy elo(lp,jlop).
                                                 (elo(lp,jlop)-e(num))*jr3upu_22(jlop,jlo,lp,l)''')
new = new.replace(
    '''                                     radial_factor = &
                                          conjg(c_am(alpha,t))* &
                                          conjg(c_am(beta,u))*gaunt1* &''',
    '''                                     ! This Gaunt ordering expands n_alpha*n_beta
                                     ! on the KET side of <Y_bra|...|Y_ket>.
                                     ! Conjugating c_am here reverses xy/yz signs.
                                     ! The spherical block uses another ordering
                                     ! and must NOT be changed in the same way.
                                     radial_factor = &
                                          c_am(alpha,t)* &
                                          c_am(beta,u)*gaunt1* &''')
new = new.replace(
    '                IF (l+v.LT.0 .OR. l+v.GT.lmax2n) CYCLE',
    '''                ! grad_Omega Y_l contains l+1 even for l=lmax2n.
                ! This is an intermediate angular harmonic, not an extra
                ! radial basis channel.  Cutting it spoils the exact zero
                ! trace n dot grad_Omega=0, including in TRACE mode.
                IF (l+v.LT.0) CYCLE''')
start = '''       IF (ih.GT.1 .AND. lpv(ih-1).EQ.l .AND. lv(ih-1).EQ.lp .AND. &'''
assert new.count(start) == 2
new = new.replace(start, '''       ! Fortran .AND. need not short circuit: guard ih-1 separately.
       IF (ih.GT.1) THEN
       IF (lpv(ih-1).EQ.l .AND. lv(ih-1).EQ.lp .AND. &''')
end = '''               channel_diff.LE.1.d-12*(1.d0+channel_scale)
       ENDIF
       IF (duplicate_radial_channel) CYCLE'''
assert new.count(end) == 2
new = new.replace(end, '''               channel_diff.LE.1.d-12*(1.d0+channel_scale)
       ENDIF
       ENDIF
       IF (duplicate_radial_channel) CYCLE''')
assert old != new
p.write_text(new)
(root/'lapw2-g2-audit.patch').write_text(''.join(difflib.unified_diff(
    old.splitlines(True), new.splitlines(True),
    fromfile='a/SRC_lapw2/val_corr_stress.F',
    tofile='b/SRC_lapw2/val_corr_stress.F')))
print('Patched four identified G2 defects in staged source.')
