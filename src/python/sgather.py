#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys,re,os
import optparse
import utils
from indmffile import Indmfl, ParsIndmfi
#from scipy import *
from numpy import array,zeros,loadtxt,savetxt,shape,pi,isfinite,all,abs,searchsorted,round
from functools import reduce

# --- frequency-mesh handling for the mixing step -----------------------------
# See docs: sgather used to pair the new and old self-energy BY ROW NUMBER and
# truncate with min(len_new,len_old).  When the solver returns a DENSER mesh than
# sig.inp (bctqmc Nambu at beta=150: solver 2385 rows vs 585 log-mesh rows), that
# silently amputated omega in [24.5, 99.84] and then stamped the truncated axis on
# top -- destroying 18 iterations with exit code 0.  Everything below aligns by
# exact integer Matsubara index instead, and keeps the SOLVER's mesh (the mode
# iNORG has always run in).

OM_SENTINEL = 5000.0      # the high-frequency marker rows carry (s_oo-Edc, 0)
SENTINEL_MIN = 1000.0     # anything above this is a marker, not a physical point
# The sentinel's dynamic part is only zero AFTER (s_oo-Edc) has been subtracted;
# before that it equals s_oo-Edc itself (measured 1.0018 for bctqmc, 0.056 for
# iNORG), so this check must never run before the subtraction.
SENTINEL_ZERO_TOL = 1e-8
# Tolerance for checking a RAW marker row against its own (s_oo-Edc, 0), before any
# averaging.  Loose enough for the solvers' '%20.15f' output (measured residual ~3e-16
# relative), tight enough that a genuinely wrong marker cannot hide.
SENTINEL_RAW_TOL = 1e-9
# Matsubara-exactness tolerance, RELATIVE ON OMEGA.  Do not re-express this as an
# absolute tolerance on n: the deviation in n grows with n, while the relative
# error in omega does not.  bctqmc's Nambu writer uses '%20.15f' (15 digits after
# the point), giving max relative error 4.8e-12 on real production data -- which
# is |n-round(n)| up to 1.2e-9 at beta=150.  A tolerance of 5e-10 on n therefore
# ABORTS on healthy data; 1e-9 relative on omega leaves ~200x margin.
OM_RTOL = 1e-9
NHEAD_CONSEC = 8          # how many leading indices must be consecutive from n=0.
                          # Deliberately small: sgather never learns the ssplit -n that
                          # produced the mesh, and a legitimate nom=50 log mesh has only
                          # 51 consecutive head points, so a large fixed value aborts
                          # healthy files.  8 is ample for its actual purpose (beta
                          # aliasing shows up at the first step: 1,4,7,... vs 0,1,2,...).
BRACKET_RATIO_MAX = 1.6   # max omega_R/omega_L for a single interpolation bracket
                          # (measured max in the production log tail: 1.4187).
                          # A policy limit on how far we will reconstruct, NOT a
                          # bound on the resulting error.
N_ENDPOINT_TRIM = 2       # Matsubara indices two log meshes may differ by at the top
NROW_ABORT_FRAC = 0.80    # kept rows < this * len(new mesh) -> abort


def strip_sentinel(rSigma, what):
    """Drop trailing omega=5000 marker rows (there can be more than one: the solver
    writes one and ssplit appends its own).  Must be called AFTER the (s_oo-Edc)
    subtraction, so a healthy marker reads as ~0 here.

    Matches the marker VALUE, not merely 'above SENTINEL_MIN': a physical
    high-frequency row whose dynamic part happened to fall under the zero tolerance
    must not be silently discarded.  Anything else above SENTINEL_MIN is rejected
    by check_raw_impurity rather than guessed at."""
    while shape(rSigma)[1] > 1 and abs(rSigma[0, -1] - OM_SENTINEL) < 1e-6:
        dyn = abs(rSigma[1:, -1])
        if len(dyn) and dyn.max() > SENTINEL_ZERO_TOL:
            raise RuntimeError(
                'sgather: %s omega=%g row is not a zero-dynamic sentinel '
                '(max|dyn|=%.3e > %.1e). Refusing to drop a row that carries data.'
                % (what, rSigma[0, -1], dyn.max(), SENTINEL_ZERO_TOL))
        print('sgather: stripping %s sentinel row omega=%g (max|dyn|=%.2e)'
              % (what, rSigma[0, -1], dyn.max() if len(dyn) else 0.0))
        rSigma = rSigma[:, :-1]
    return rSigma


def matsubara_index(om, beta):
    "n such that om = (2n+1)*pi/beta"
    return (om*beta/pi - 1.0)/2.0


def validate_mesh(om, beta, what, nhead=None):
    """Abort on anything that would make index alignment meaningless.

    nhead: how many leading indices must be consecutive from n=0.  Keep this SMALL.
    sgather is not told the ssplit -n/-t that produced the mesh, and the dense-head
    length cannot be inferred from the file: a legitimate create_log_mesh(nom=50,
    ntail=30) output has only 51 consecutive head points (ctqmc's own default nom is
    50, src/impurity/ctqmc.py), so any larger fixed requirement aborts a healthy file.
    A short prefix is sufficient for what this check is actually for -- beta aliasing
    shows up at the very first step (1,4,7,... instead of 0,1,2,...), not deep in the
    mesh.
    """
    if nhead is None:
        nhead = NHEAD_CONSEC
    if len(om) < 2:
        raise RuntimeError('sgather: %s has only %d frequency rows' % (what, len(om)))
    if not all(isfinite(om)):
        raise RuntimeError('sgather: %s frequency mesh contains NaN/Inf' % what)
    if not all(om[1:] > om[:-1]):
        raise RuntimeError('sgather: %s frequency mesh is not strictly increasing '
                           '(duplicate or out-of-order rows)' % what)
    # every point must be an exact fermionic Matsubara point of this beta
    nn = matsubara_index(om, beta)
    ni = round(nn)
    om_exact = (2.0*ni + 1.0)*pi/beta
    rel = abs(om - om_exact)/om_exact
    if rel.max() > OM_RTOL:
        ibad = int(rel.argmax())
        raise RuntimeError('sgather: %s omega[%d]=%.12g is not a Matsubara point of '
                           'beta=%.6g (n=%.9f, relative error %.3e > %.1e)'
                           % (what, ibad, om[ibad], beta, nn[ibad], rel.max(), OM_RTOL))
    # The head must be n = 0,1,2,...  This is the ONLY check that catches beta
    # aliasing: every beta=50 Matsubara point is also a beta=150 point
    # ((2n+1)/50 = (2m+1)/150 => m = 3n+1), and this project runs iNORG at
    # beta=50 and bctqmc at beta=150.  A beta=50 mesh read as beta=150 has head
    # indices 1,4,7,... and fails here immediately.
    ni = ni.astype(int)
    nchk = min(nhead, len(ni))
    if ni[0] != 0:
        raise RuntimeError('sgather: %s does not start at the first Matsubara point '
                           '(n[0]=%d, expected 0) -- wrong beta?' % (what, ni[0]))
    for i in range(1, nchk):
        if ni[i] != ni[i-1] + 1:
            raise RuntimeError('sgather: %s head is not a consecutive Matsubara run: '
                               'n[%d]=%d follows n[%d]=%d -- wrong beta (aliasing) or '
                               'corrupted mesh?' % (what, i, ni[i], i-1, ni[i-1]))
    return ni


def check_raw_impurity(om, data, s_oo, Edc, ncols_expected, filename, om_ref, nhead=None):
    """Validate ONE raw imp.?/sig.out before its columns are folded into rSigma.

    This has to happen per-file and unconditionally (also at mix=1): the aggregate
    rSigma carries only imp.0's axis (see the `om = loadtxt(imp.0).transpose()[0]`
    upstream), so once the columns of imp.1 have been added by row position, a
    disagreement between the two axes is unrecoverable -- and equal-length,
    different-index files do not even raise.  Likewise a NaN or a non-Matsubara
    axis at mix=1 would otherwise reach sig.inp with exit 0.
    """
    beta = pi/om[0]
    if not all(isfinite(data)):
        raise RuntimeError('sgather: %s contains NaN/Inf in its dynamic columns'
                           % filename)
    # The trailing omega=5000 marker is deliberately NOT a Matsubara point, so the
    # mesh checks run on the physical rows only.  The marker itself is checked below
    # (it must be the known value, and its dynamic part is verified once the
    # (s_oo-Edc) subtraction has been applied downstream, never before).
    nphys = len(om)
    while nphys > 1 and abs(om[nphys-1] - OM_SENTINEL) < 1e-6:
        nphys -= 1
    validate_mesh(om[:nphys], beta, 'impurity self-energy %s' % filename, nhead=nhead)
    if len(s_oo) != len(Edc):
        raise RuntimeError('sgather: %s header lengths disagree: len(s_oo)=%d, '
                           'len(Edc)=%d' % (filename, len(s_oo), len(Edc)))
    if len(s_oo) < ncols_expected:
        raise RuntimeError('sgather: %s header has %d entries but %d self-energy '
                           'components are expected from the Sigind of this impurity'
                           % (filename, len(s_oo), ncols_expected))
    if shape(data)[0] < 2*ncols_expected + 1:
        raise RuntimeError('sgather: %s has %d columns, too few for the %d self-energy '
                           'components expected from its Sigind'
                           % (filename, shape(data)[0], ncols_expected))
    # The marker of each RAW file is verified HERE, on that file alone, before any
    # averaging.  The post-aggregation zero-check cannot do this job: two impurities
    # contributing to the same column with opposite marker residuals would cancel
    # below tolerance in the aggregate and pass while each raw row is bad.
    if om[-1] >= SENTINEL_MIN:
        if abs(om[-1] - OM_SENTINEL) > 1e-6:
            raise RuntimeError('sgather: %s ends at omega=%g, which is above %g but is not '
                               'the omega=%g high-frequency marker'
                               % (filename, om[-1], SENTINEL_MIN, OM_SENTINEL))
        # the raw marker row holds (s_oo[i]-Edc[i], 0) per component, pre-subtraction
        for i in range(ncols_expected):
            re_got, im_got = data[1+2*i, -1], data[2+2*i, -1]
            re_exp = s_oo[i] - Edc[i]
            if abs(re_got - re_exp) > SENTINEL_RAW_TOL*(1.0 + abs(re_exp)) \
               or abs(im_got) > SENTINEL_RAW_TOL:
                raise RuntimeError('sgather: %s omega=%g marker row disagrees with its own '
                                   'header at component %d: got (%.6g, %.6g), expected '
                                   '(s_oo-Edc, 0) = (%.6g, 0)'
                                   % (filename, om[-1], i+1, re_got, im_got, re_exp))
    if om_ref is not None:
        if len(om) != len(om_ref) or not all(abs(om - om_ref) <= OM_RTOL*abs(om_ref)):
            raise RuntimeError('sgather: %s is on a different frequency mesh than the '
                               'first impurity (%d vs %d rows). Impurity self-energies '
                               'are combined by row, so their meshes must be identical.'
                               % (filename, len(om), len(om_ref)))


def check_raw_real_axis(om, data, s_oo, Edc, ncols_expected, filename, om_ref):
    """Structural validation of ONE raw real-axis imp.?/sig.out.

    Same aggregation hazard as the Matsubara case (columns are folded in by row onto
    imp.0's axis), minus everything Matsubara-specific: a real mesh spans negative to
    positive frequencies, so there is no index to align on and no marker row.
    """
    if len(om) < 2:
        raise RuntimeError('sgather: %s has only %d frequency rows' % (filename, len(om)))
    if not all(isfinite(data)):
        raise RuntimeError('sgather: %s contains NaN/Inf' % filename)
    if not all(om[1:] > om[:-1]):
        raise RuntimeError('sgather: %s real-frequency mesh is not strictly increasing '
                           '(duplicate or out-of-order rows)' % filename)
    if len(s_oo) != len(Edc):
        raise RuntimeError('sgather: %s header lengths disagree: len(s_oo)=%d, len(Edc)=%d'
                           % (filename, len(s_oo), len(Edc)))
    if len(s_oo) < ncols_expected:
        raise RuntimeError('sgather: %s header has %d entries but %d self-energy '
                           'components are expected from the Sigind of this impurity'
                           % (filename, len(s_oo), ncols_expected))
    if shape(data)[0] < 2*ncols_expected + 1:
        raise RuntimeError('sgather: %s has %d columns, too few for the %d self-energy '
                           'components expected from its Sigind'
                           % (filename, shape(data)[0], ncols_expected))
    if om_ref is not None:
        if len(om) != len(om_ref) or not all(abs(om - om_ref) <= 1e-8*(1.0 + abs(om_ref))):
            raise RuntimeError('sgather: %s is on a different real-frequency mesh than the '
                               'first impurity (%d vs %d rows). Impurity self-energies are '
                               'combined by row, so their meshes must be identical.'
                               % (filename, len(om), len(om_ref)))


def mix_on_matsubara(rSigma, rSigma_old, mix):
    """Mix new and old self-energies at the SAME frequency, output on the NEW mesh.

    The frequency axis is metadata and never participates in the mixing.  Old
    values are taken by exact integer Matsubara index where present, else
    interpolated in x=1/omega between the bracketing old points -- a convex
    combination with REAL, COLUMN-IDENTICAL weights, which is what keeps the
    Nambu tau3 relations (Sigma_hh=-conj(Sigma_pp), Phi^dag=conj(Phi)) exact and
    Im(Sigma)<=0 intact.  Never use per-column weights here.
    """
    om_new = rSigma[0]
    om_old = rSigma_old[0]
    if shape(rSigma)[0] != shape(rSigma_old)[0]:
        raise RuntimeError('sgather: new self-energy has %d columns but the old one has '
                           '%d -- refusing to mix different formats'
                           % (shape(rSigma)[0], shape(rSigma_old)[0]))
    if not all(isfinite(rSigma)) or not all(isfinite(rSigma_old)):
        raise RuntimeError('sgather: self-energy contains NaN/Inf')

    # beta from the first Matsubara point; both files must agree, otherwise we
    # would be averaging two temperatures.
    beta = pi/om_new[0]
    if abs(om_old[0] - om_new[0]) > OM_RTOL*om_new[0]:
        raise RuntimeError('sgather: first Matsubara frequency differs between the new '
                           '(%.12g, beta=%.6g) and old (%.12g, beta=%.6g) self-energy '
                           '-- refusing to mix two temperatures'
                           % (om_new[0], beta, om_old[0], pi/om_old[0]))
    nhead = min(NHEAD_CONSEC, len(om_new), len(om_old))
    n_new = validate_mesh(om_new, beta, 'new self-energy (imp.?/sig.out)', nhead=nhead)
    n_old = validate_mesh(om_old, beta, 'old self-energy (%s)' % 'sig.inp', nhead=nhead)

    # fast path: the two meshes are already the same
    if len(n_new) == len(n_old) and all(n_new == n_old):
        out = mix*rSigma + (1-mix)*rSigma_old
        out[0] = om_new
        print('sgather: mixing on identical mesh, %d rows (mix=%g)' % (len(n_new), mix))
        return out

    # mismatch path: align by integer index, output on the new mesh
    pos = {}
    for i in range(len(n_old)):
        pos[int(n_old[i])] = i
    n_old_max = int(n_old[-1])
    x_old = 1.0/om_old
    n_head_old = 1
    while n_head_old < len(n_old) and n_old[n_head_old] == n_old[n_head_old-1] + 1:
        n_head_old += 1

    keep = []
    old_rows = []
    n_exact = n_interp = n_trim = 0
    max_ratio = 1.0
    max_err = 0.0
    for i in range(len(n_new)):
        n = int(n_new[i])
        if n in pos:
            old_rows.append(rSigma_old[:, pos[n]])
            keep.append(i)
            n_exact += 1
        elif n > n_old_max:
            # the endpoints of two log meshes can differ by a point or two
            # (99.84 vs 99.88); more than that means the old file is truncated.
            if n - n_old_max <= 2:
                n_trim += 1
                continue
            raise RuntimeError('sgather: the old self-energy stops at Matsubara index %d '
                               '(omega=%.6g) while the solver reaches %d (omega=%.6g). '
                               'That is a truncated/poisoned sig.inp, not a mesh '
                               'refinement -- aborting instead of propagating it.'
                               % (n_old_max, om_old[-1], n, om_new[i]))
        else:
            j = int(searchsorted(n_old, n))          # n_old[j-1] < n < n_old[j]
            # A missing index strictly BELOW the end of the old file's own consecutive
            # prefix means that prefix has a hole in it -- corrupted, not a log tail.
            # Note this cannot be decided from j alone: the first legitimate tail gap
            # (dense head n=0..500, tail resuming at 503) and a genuine hole right
            # after the prefix (n=0..99 then 101) BOTH land on j == n_head_old.  The
            # index value is what separates them.
            if n < int(n_old[n_head_old-1]):
                raise RuntimeError('sgather: the old self-energy is missing Matsubara index '
                                   '%d inside its own dense head (n=0..%d) -- corrupted '
                                   'mesh, aborting.' % (n, int(n_old[n_head_old-1])))
            ratio = om_old[j]/om_old[j-1]
            if ratio > BRACKET_RATIO_MAX:
                raise RuntimeError('sgather: interpolation bracket omega=[%.6g, %.6g] is too '
                                   'wide (ratio %.4f > %.2f) to reconstruct the old '
                                   'self-energy at omega=%.6g -- aborting rather than '
                                   'guessing.' % (om_old[j-1], om_old[j], ratio,
                                                  BRACKET_RATIO_MAX, om_new[i]))
            if ratio > max_ratio: max_ratio = ratio
            # linear in x=1/omega: Im(Sigma) ~ 1/omega is exactly linear there, so
            # the 1/(i*omega) moment is represented exactly and Re ~ 1/omega**2
            # carries only second-order error.
            xL = x_old[j-1]; xR = x_old[j]
            wR = (1.0/om_new[i] - xL)/(xR - xL)
            old_rows.append((1.0-wR)*rSigma_old[:, j-1] + wR*rSigma_old[:, j])
            keep.append(i)
            n_interp += 1
            # Audit metadata ONLY -- the spread of the bracket endpoints, scaled by
            # w(1-w).  This is NOT a bound on the interpolation error: with equal
            # endpoints and an interior feature it reads zero while the true error is
            # arbitrary.  Measured on real data it reads 2.6e-2 where the true max
            # error is 1.9e-3.  Logged to spot pathological brackets, never enforced.
            err = abs(rSigma_old[1:, j] - rSigma_old[1:, j-1]).max()*wR*(1.0-wR)
            if err > max_err: max_err = err

    if not keep:
        raise RuntimeError('sgather: no frequency of the new self-energy could be matched '
                           'against the old one')
    rSigma = rSigma[:, keep]
    old_aligned = array(old_rows).transpose()
    out = mix*rSigma + (1-mix)*old_aligned
    out[0] = rSigma[0]                       # the new axis, never mixed

    # Range guards.  These are INDEX-based, not fractional: a fractional endpoint
    # test like om_out_max >= 0.98*om_old[-1] passes while losing ~47 of the top
    # Matsubara points at beta=150, and repeating that allowed 2% per iteration
    # compounds to ~18% after ten -- a slow-motion version of the very bug being
    # fixed.  In Matsubara index the tolerance is unambiguous: the only endpoint
    # discrepancy we accept is the couple of indices two log meshes can differ by.
    om_out_max = out[0, -1]
    n_out_max = int(round(matsubara_index(om_out_max, beta)))
    n_old_top = int(n_old[-1])
    if n_out_max < n_old_top - N_ENDPOINT_TRIM:
        raise RuntimeError('sgather: the mixed self-energy would stop at Matsubara index '
                           '%d (omega=%.6g), down from %d (omega=%.6g) in the previous '
                           'iteration -- more than the %d indices two meshes may '
                           'legitimately differ by. That is the high-frequency '
                           'truncation this guard exists to stop.'
                           % (n_out_max, om_out_max, n_old_top, om_old[-1],
                              N_ENDPOINT_TRIM))
    # Coverage of the NEW mesh: this is what detects rows being dropped relative to
    # what the solver actually produced.  min(new,old) would go vacuous exactly when
    # the new mesh is the shorter one.
    if shape(out)[1] < len(n_new) - N_ENDPOINT_TRIM:
        nlost = len(n_new) - shape(out)[1]
        if shape(out)[1] < NROW_ABORT_FRAC*len(n_new):
            raise RuntimeError('sgather: the mixed self-energy kept only %d of the %d '
                               'frequencies the solver produced (%d dropped) -- too much '
                               'of the mesh was lost, aborting.'
                               % (shape(out)[1], len(n_new), nlost))
        print('WARNING sgather: kept %d of the %d solver frequencies (%d dropped)'
              % (shape(out)[1], len(n_new), nlost))

    print('sgather: mesh-align mix=%g beta=%.6g n_old=%d n_new=%d n_out=%d '
          'om_max old=%.6g new=%.6g out=%.6g n_max old=%d out=%d '
          'exact=%d interp=%d trimmed=%d max_bracket_ratio=%.4f '
          'interp_endpoint_spread=%.3e'
          % (mix, beta, shape(rSigma_old)[1], len(n_new), shape(out)[1],
             om_old[-1], om_new[-1], om_out_max, n_old_top, n_out_max,
             n_exact, n_interp, n_trim, max_ratio, max_err))
    return out

#def union(data):
#    " Takes a union of array or list"
#    c = []
#    for d in data:
#        if d not in c:
#            c.append(d)
#    return c
#
#def SimplifySiginds(siginds):
#    " Takes dictionary of Sigind's and creates list or non-zero columns"
#    colsp={}
#    colsm={}
#    for icix in list(siginds.keys()):
#        Sigind = siginds[icix]
#        cols_all = union(array(Sigind).flatten())
#        cols_all = sorted(cols_all,key=lambda x: abs(x))
#        colp=[x for x in cols_all if x>0]
#        colm=[x for x in cols_all if x<0]
#        #col = sort(union(array(Sigind).flatten())).tolist()
#        #if 0 in col: col.remove(0)
#        colsp[icix] = colp
#        colsm[icix] = colm
#    return (colsp,colsm)

def SimplifySiginds(siginds):
    """Takes dictionary of Sigind's and creates list of non-zero values
    that appear in Siginds, cols[icix]=[1,2,3,...]
    """
    colsp={}
    colsm={}
    for icix in siginds:
        unique_values = set(siginds[icix].ravel())
        colsp[icix] = sorted([x for x in unique_values if x>0])
        colsm[icix] = sorted([x for x in unique_values if x<0], reverse=True)
        # should remove all negative indices!!!
    return (colsp,colsm)

if __name__=='__main__':
    """ Takes the self-energy files from all impurity problems, and combines
    them into a single self-energy file.
    """
    usage = """usage: %prog [ options ]

    The script takes the self-energy files from all impurity problems,
    and combines them into a single self-energy file.

    To give filename, you can use the following expressions:
      - word 'case', which will be replaced by current case, determined by
                    the presence of struct file.
      - '?', which will be replaced by icix.
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-o", "--osig", dest="osig", default='sig.inp', help="filename of the output self-energy file. Default: 'sig.inp'")
    parser.add_option("-i", "--isig", dest="isig", default='imp.?/sig.out', help="filename of the input self-energy from all impurity problems. Default: 'imp.?/sig.out'")
    parser.add_option("-l", "--lext", dest="m_extn", default='', help="For magnetic calculation, it can be 'dn'.")
    parser.add_option("-m", "--mix", dest="mix", type=float, default=1.0, help="Mixing parameter for self-energy. Default=1 -- no mixing")

    # Next, parse the arguments
    (options, args) = parser.parse_args()

    # A negative mix silently EXTRAPOLATES away from both self-energies; >1 does the
    # same on the other side.  run_dmft.py only ever passes 0<mix<1, but this script
    # is also run by hand.
    if not (0.0 <= options.mix <= 1.0):
        print('ERROR sgather: -m/--mix must be in [0,1], got %g' % options.mix)
        sys.exit(1)

    env = utils.W2kEnvironment()
    case = env.case
    
    options.osig = re.sub(r'case', case, options.osig)
    options.isig = re.sub(r'case', case, options.isig)
    
    print('sgather.py: case=%s, isig=%s, osig=%s' %  (case, options.isig, options.osig))
    

    inl = Indmfl(case)
    inl.read()
    if options.m_extn:
        inldn = Indmfl(case, 'indmfl'+options.m_extn)
        inldn.read()

    # Impurity quantities
    iSiginds = ParsIndmfi(case)
    icols,icolsm = SimplifySiginds(iSiginds)
    print('sgather: icols=', icols)

    # lattice quantities
    colsp,colsm = SimplifySiginds(inl.siginds)
    print('sgather: colsp=', colsp, 'colsm=', colsm)

    if options.m_extn:
        colspdn, colsmdn = SimplifySiginds(inldn.siginds)
        print('sgather: colspdn=', colspdn, 'colsmdn=', colsmdn)


    # these are columns we can fill in whit our impurity problems
    allcols = sorted( reduce(lambda x,y: x+y, list(icols.values())) )
    print('sgather: allcols=', allcols)
    #print('len(allcols)=', len(allcols))
    noccur = zeros(max(allcols),dtype=int)
    #print('len(noccur)=', len(noccur))
    # these columns can not be obtained from impurity problems. Need another way.
    missing_columns={}
    for icix in list(colsm.keys()):
        # This needs particular attention. Note that there are two types of negative cix indices :
        # a) If both positive and negative indices appear for the same orbital (same icix) then we expect
        #    to just merely shift the spectra that corresponds to the negative columns, and in this case
        #    we add such terms in s_oo and Edc.
        # b) If all indices of certain orbitals are negative, such orbital is treated as "open core" orbital,
        #    and is not computed by solving an impurity problem. Rather it is split by a fixed self-energy,
        #    stored in the file sfx.[icix]. This case needs particular attention in ssplit.py part, but
        #    not in sgather.py step, as such "open core" orbitals do not appear in sig.inp, and hence can
        #    be safely ignored in this step.
        if len(colsp[icix]) > 0:  # only if some orbitals are positive in this icix, such case is not "open core" case, and we need to treat it.
            for i in colsm[icix]: # missing columns do not appear in self-energy
                missing_columns[abs(i)]=1
    print('sgather: missing_columns=', missing_columns)

    for c in allcols:
        noccur[c-1]+=1
    print('sgather: noccur=', noccur)
    
    filename = re.sub(r'\?', str(0), options.isig)
    om = loadtxt(filename).transpose()[0]
    
    # Array of all sigmas
    rSigma=zeros(( len(allcols)*2+1, len(om) ),dtype=float)
    rs_oo=zeros( len(allcols)+len(missing_columns) ,dtype=float)
    rEdc=zeros( len(allcols)+len(missing_columns) ,dtype=float)
    rSigma[0] = om
    print('sgather: shape(rSigma)', shape(rSigma))
    print('sgather: shape(rs_oo)=', shape(rs_oo))

    # Reading self-energies from all impurity problems
    for icix in list(icols.keys()):
        # Processing Delta
        #mincol = min(icols[icix])
        #maxcol = max(icols[icix])
        
        filename = re.sub(r'\?', str(icix), options.isig)
        print(('sgather: icix=%d reading from: %s' % (icix, filename)), 'cols=', icols[icix])
        
        # Searching for s_oo and Edc
        with open(filename, 'r') as fh_sig:
            line1 = next(fh_sig)[1:]
            line2 = next(fh_sig)[1:]
            print('sgather: line1=', line1,end='')
            print('sgather: line2=', line2,end='')
            exec(line1.strip(), globals())
            exec(line2.strip(), globals())
        # reading sigma
        data = loadtxt(filename).transpose()

        # Validate THIS raw file before folding its columns in.  Unconditional in the
        # Matsubara case -- the mixing branch below only ever sees the aggregate on
        # imp.0's axis, so a per-impurity disagreement is unrecoverable by then.
        _om_ref = None if icix == list(icols.keys())[0] else om
        if inl.matsubara:
            check_raw_impurity(data[0], data, s_oo, Edc, len(icols[icix]), filename,
                               _om_ref)
        else:
            # Real axis: no Matsubara index to check, but the row-aggregation hazard is
            # identical -- an equal-length imp.1 on a different mesh would be stamped
            # onto imp.0's axis.  So the structural checks still run, unconditionally.
            check_raw_real_axis(data[0], data, s_oo, Edc, len(icols[icix]), filename,
                                _om_ref)

        for iic,c in enumerate(icols[icix]):
            imiss = len([x for x in missing_columns if x<=c])
            ic = c-imiss
            #print 'writting into ', 2*ic-1, 'and', 2*ic, 'from data at ', 2*(c-mincol+1)-1, 'and', 2*(c-mincol+1)
            print('sgather: writting into ', 2*ic-1, 'and', 2*ic, 'from data at ', 2*iic+1, 'and', 2*iic+2)
            #rSigma[2*ic-1] += data[2*(c-mincol+1)-1]*(1./noccur[c-1])
            #rSigma[2*ic]   += data[2*(c-mincol+1)  ]*(1./noccur[c-1])
            rSigma[2*ic-1] += data[2*iic+1]*(1./noccur[c-1])
            rSigma[2*ic]   += data[2*iic+2]*(1./noccur[c-1])
        
        for iic,c in enumerate(icols[icix]):
            #print 'writting s_oo into ', c-1, 'from', c-mincol
            print('writting s_oo into ', c-1, 'from', iic)
            rs_oo[c-1] += s_oo[iic]*(1./noccur[c-1])
            rEdc[c-1]  += Edc[iic]*(1./noccur[c-1])
        
    print('sgather: rs_oo=', rs_oo, 'rEdc=', rEdc)
    # The impurity self-energy has a static contribution
    # equal to s_oo-Edc
    # This self-energy will however go to zero in infinity. No static contribution
    for c in allcols:
        imiss = len([x for x in missing_columns if x<=c])
        ic = c-imiss
        rSigma[2*ic-1] -= (rs_oo[c-1]-rEdc[c-1])

    # Drop the omega=5000 marker row the solver appended.  Unconditional (also at
    # mix=1): ssplit.py adds its own sinf=[5000.] row, so leaving this one in gives
    # sig.inp1 two of them, and create_log_mesh then treats the marker as a data
    # point.  This has to come after the subtraction above -- before it, the row
    # legitimately holds s_oo-Edc and the zero-check would fire on every run.
    if inl.matsubara:
        rSigma = strip_sentinel(rSigma, 'new self-energy')


    #cm = colsm.values()
    #if options.m_extn:
    #    cm += colsmdn.values()
    #cm=array(union(cm)).flatten()
    #print 'cm=', cm

    if len(missing_columns)>0: # Some columns contain only s_oo and Edc (no dynamic component)
        # old-s_oo and old-Edc
        with open(options.osig,'r') as fi: # checking old s_oo and Edc from sig.inp
            exec(next(fi)[1:].strip(), globals())
            exec(next(fi)[1:].strip(), globals())
            
        for i in list(missing_columns.keys()):
            rs_oo[i-1] = s_oo[i-1]
            rEdc[i-1] = Edc[i-1]

    if options.mix!=1.0 and os.path.isfile(options.osig) and os.path.getsize(options.osig)>0:
        print('Mixing self-energy with mix=', options.mix)
        with open(options.osig,'r') as fi: # checking old s_oo and Edc
            exec(next(fi)[1:].strip(), globals())
            exec(next(fi)[1:].strip(), globals())
            
        rSigma_old = loadtxt(options.osig).transpose()
        print('shape sigma old', shape(rSigma_old))
        print('shape sigma new', shape(rSigma))

        if inl.matsubara:
            # The old file may also still carry a marker row (an older sig.inp written
            # before this fix, or one that came back through ssplit).  Its dynamic part
            # was already reduced by the s_oo-Edc of ITS own iteration, so it reads as
            # ~0 here, exactly like the new one.
            rSigma_old = strip_sentinel(rSigma_old, 'old self-energy')

            # Mix at equal frequency, aligned by exact integer Matsubara index, output
            # on the solver's mesh.  Replaces the old pair-by-row-number + min() logic,
            # which silently truncated the high-frequency tail whenever the solver mesh
            # was denser than sig.inp.
            rSigma = mix_on_matsubara(rSigma, rSigma_old, options.mix)
        else:
            # Real-frequency axis: there is no Matsubara index to align on, so only the
            # same-mesh case is well defined.  (This is the pre-existing behaviour minus
            # the min()-truncation, which was never right here either.)
            if shape(rSigma) != shape(rSigma_old):
                raise RuntimeError('sgather: real-axis self-energies have different shapes '
                                   '(new %s vs old %s); refusing to mix by row position'
                                   % (shape(rSigma), shape(rSigma_old)))
            if not all(abs(rSigma[0] - rSigma_old[0]) <= 1e-8*(1.0 + abs(rSigma_old[0]))):
                raise RuntimeError('sgather: real-axis self-energies are on different '
                                   'frequency meshes; refusing to mix by row position')
            new_om = rSigma[0].copy()
            rSigma = options.mix*rSigma + (1-options.mix)*rSigma_old
            rSigma[0] = new_om                  # the axis is metadata, never mixed

        # s_oo and Edc take the same alpha as the dynamic columns.  The stored
        # quantity is (Sigma - (s_oo-Edc)), so mixing all three with one alpha is
        # exactly equivalent to mixing the full self-energy -- with a different
        # alpha here the reconstruction downstream would be silently wrong.
        rs_oo[:] = options.mix*rs_oo[:]  + (1-options.mix)*array(s_oo)
        rEdc[:]  = options.mix*rEdc[:]   + (1-options.mix)*array(Edc)


    print('osig=', options.osig)
    # Write through a temporary file and rename: jobs in this campaign are killed
    # routinely, and a half-written sig.inp is unrecoverable.
    tmp = options.osig + '.tmp'
    with open(tmp, 'w') as fout:
        print('# s_oo=', rs_oo.tolist(), file=fout)
        print('# Edc=', rEdc.tolist(), file=fout)
        savetxt(fout,rSigma.transpose())
    os.replace(tmp, options.osig)
