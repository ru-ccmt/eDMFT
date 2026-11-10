from pathlib import Path
root=Path('/private/tmp/stress_remaining_audit_20261107')

def change(p,old,new,count=1):
 s=p.read_text();assert s.count(old)==count,(p.name,s.count(old),old[:80]);p.write_text(s.replace(old,new))
core=root/'SRC_lcore/core_corr_stress.f'
change(core,'         ll_p = lmtot(1,lm1p)','''         ! The sign of a WIEN lattice-harmonic l labels the sine/cosine
         ! combination, not angular momentum.  multfc_core has already
         ! encoded that sign in the complex coefficients.  Gaunt triangle
         ! rules and C_FACTOR require the physical, nonnegative l.
         ll_p = abs(lmtot(1,lm1p))''',2)
# The same convention occurs in the electrostatic angular contraction.
p=root/'SRC_lapw0/finl_elect_str.F'
s=p.read_text();s=s.replace('l0 = lmtot(1,lm1)','l0 = abs(lmtot(1,lm1))').replace('l0p = lmtot(1,lm1p)','l0p = abs(lmtot(1,lm1p))').replace('l0   = lmtot(1,lm1)','l0   = abs(lmtot(1,lm1))');p.write_text(s)
