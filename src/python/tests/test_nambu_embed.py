#!/usr/bin/env python
"""Standalone checks of the Nambu embedding / back-projection used by dmftgk.py --nambu.

Pure numpy, no MPI, no import of dmftgk.py (which pulls in mpi4py and pylab).
The two candidate embedding conventions are re-implemented here in a few lines and
compared against an EXACT Gor'kov Green's function built in the orbital basis with no
projector anywhere -- the oracle is deliberately not written in terms of the formulas
under test.

Run:  <python> src/python/tests/test_nambu_embed.py
Exit 0 = all checks pass; nonzero + a message on stderr = failure.
"""
import sys
import numpy as np


# ----------------------------------------------------------------------------------
# Model.  A 1D tight-binding chain with n orbitals per cell.  The real-space hoppings
# are drawn ONCE, in the constructor, and are REAL, so that H(-k) == conj(H(k)).
# Drawing randoms inside H(k) would give k and -k different Hamiltonians and silently
# invalidate every k -> -k statement below.
# ----------------------------------------------------------------------------------
class Model:
    def __init__(self, n, nc, seed, rot=None):
        self.n, self.nc = n, nc
        rng = np.random.default_rng(seed)
        a0 = rng.normal(size=(n, n))
        t0 = 0.5 * (a0 + a0.T)              # real symmetric on-site block
        t1 = rng.normal(size=(n, n))        # real nearest-neighbour hopping
        if rot is not None:
            # rot is REAL and orthogonal, so the rotated hoppings stay real and
            # H(-k) == conj(H(k)) survives the basis change.
            t0 = (rot @ t0 @ rot.conj().T).real
            t1 = (rot @ t1 @ rot.conj().T).real
        self.t0, self.t1 = t0, t1
        # local Nambu self-energy blocks; the physical relations are imposed by hand.
        A = rng.normal(size=(nc, nc)) + 1j * rng.normal(size=(nc, nc))
        self.S_pp = 0.5 * (A + A.T)                     # symmetric complex
        B = rng.normal(size=(nc, nc)) + 1j * rng.normal(size=(nc, nc))
        self.Phi = 0.25 * (B + B.T)                     # symmetric complex
        self.S_hh = -np.conj(self.S_pp)
        self.Phi_dag = np.conj(self.Phi)

    def Hk(self, k):
        """H(k) = t0 + t1 e^{ik} + t1^T e^{-ik}.  Hermitian, H(-k) = conj(H(k))."""
        return self.t0 + self.t1 * np.exp(1j * k) + self.t1.T * np.exp(-1j * k)

    def pad(self, M):
        out = np.zeros((self.n, self.n), dtype=complex)
        out[:self.nc, :self.nc] = M
        return out

    # ---- ORACLE: exact Nambu G in the orbital basis, no projector anywhere ---------
    def exact_blocks(self, k, z):
        n, nc = self.n, self.nc
        Ginv = np.zeros((2 * n, 2 * n), dtype=complex)
        Ginv[:n, :n] = z * np.eye(n) - self.Hk(k) - self.pad(self.S_pp)
        Ginv[n:, n:] = z * np.eye(n) + self.Hk(-k).T - self.pad(self.S_hh)
        Ginv[:n, n:] = -self.pad(self.Phi)
        Ginv[n:, :n] = -self.pad(self.Phi_dag)
        G = np.linalg.inv(Ginv)
        return (G[:nc, :nc], G[:nc, n:n + nc], G[n:n + nc, :nc], G[n:n + nc, n:n + nc])

    def projector(self, k):
        e, V = np.linalg.eigh(self.Hk(k))
        return e, np.conj(V[:self.nc, :])

    # ---- PROJECTOR ROUTE ----------------------------------------------------------
    def projector_blocks(self, k, z, conv):
        e, U = self.projector(k)
        nb = len(e)
        c = CONV[conv]
        Spp = embed(self.S_pp, U, c['emb']['pp'])
        Sph = embed(self.Phi, U, c['emb']['ph'])
        Shp = embed(self.Phi_dag, U, c['emb']['hp'])
        Shh = embed(self.S_hh, U, c['emb']['hh'])
        Ginv = np.zeros((2 * nb, 2 * nb), dtype=complex)
        Ginv[:nb, :nb] = np.diag(z - e) - Spp
        Ginv[nb:, nb:] = np.diag(z + e) - Shh
        Ginv[:nb, nb:] = -Sph
        Ginv[nb:, :nb] = -Shp
        G = np.linalg.inv(Ginv)
        return (back(G[:nb, :nb], U, c['bak']['pp']),
                back(G[:nb, nb:], U, c['bak']['ph']),
                back(G[nb:, :nb], U, c['bak']['hp']),
                back(G[nb:, nb:], U, c['bak']['hh']))

    def ksum_exact(self, z, nk=24):
        """Sum the exact local blocks over a BZ mesh closed under k -> -k."""
        ks = 2 * np.pi * np.arange(nk) / nk
        acc = [np.zeros((self.nc, self.nc), dtype=complex) for _ in range(4)]
        for k in ks:
            for a, b in zip(acc, self.exact_blocks(k, z)):
                a += b / nk
        return acc


# 'n' means U (not conjugated), 'c' means conj(U).
def embed(M, U, flags):
    a = U if flags[0] == 'n' else np.conj(U)
    b = U if flags[1] == 'n' else np.conj(U)
    return np.einsum('ab,am,bn->mn', M, a, b)


def back(Gblk, U, flags):
    a = U if flags[0] == 'n' else np.conj(U)
    b = U if flags[1] == 'n' else np.conj(U)
    return np.einsum('am,bn,mn->ab', a, b, Gblk)


CONV = {
    # ALPHA is what dmftgk.py does: ONE tensor (STrans / GTrans) for all four blocks.
    'ALPHA': dict(emb=dict(pp='nc', ph='nc', hp='nc', hh='nc'),
                  bak=dict(pp='cn', ph='cn', hp='cn', hh='cn')),
    # BETA is the known-wrong alternative, kept only as a negative control.
    'BETA':  dict(emb=dict(pp='nc', ph='nn', hp='cc', hh='cn'),
                  bak=dict(pp='cn', ph='cc', hp='nn', hh='nc')),
}

# The convention T1/T5d assert to be exact.  Flip this to 'BETA' to prove the test
# can actually fail (see the acceptance run in the report).
CONV_UNDER_TEST = 'ALPHA'


def maxdiff(a, b):
    return max(np.max(np.abs(x - y)) for x, y in zip(a, b))


FAILURES = []


def check(name, ok, msg):
    print(('PASS ' if ok else 'FAIL ') + name + ': ' + msg)
    if not ok:
        FAILURES.append(name + ': ' + msg)


# generic k -- deliberately NOT Gamma, and nonzero frequency: at Gamma with omega=0
# the two conventions agree to 4e-16 and the negative control would be vacuous.
K_GEN = 0.837
Z_MATS = 0.7j
Z_REAL = 0.35 + 0.02j

M = Model(n=5, nc=3, seed=20260815)


def test_T0_model():
    d = np.max(np.abs(M.Hk(-K_GEN) - np.conj(M.Hk(K_GEN))))
    check('T0a H(-k)==conj(H(k))', d < 1e-14, f'max|H(-k)-H(k)*| = {d:.3e}')
    dh = np.max(np.abs(M.Hk(K_GEN) - M.Hk(K_GEN).conj().T))
    check('T0b H hermitian', dh < 1e-14, f'max|H-H^dag| = {dh:.3e}')
    imu = np.max(np.abs(M.projector(K_GEN)[1].imag))
    check('T0c projector genuinely complex', imu > 0.1, f'max|Im U| = {imu:.3e}')


def test_T1_alpha_exact():
    for tag, z in (('matsubara', Z_MATS), ('real-axis', Z_REAL)):
        d = maxdiff(M.projector_blocks(K_GEN, z, CONV_UNDER_TEST), M.exact_blocks(K_GEN, z))
        check(f'T1 {CONV_UNDER_TEST}==exact ({tag})', d < 1e-12,
              f'max|{CONV_UNDER_TEST}-exact| = {d:.3e} (need <1e-12)')


def test_T2_beta_negative_control():
    for tag, z in (('matsubara', Z_MATS), ('real-axis', Z_REAL)):
        d = maxdiff(M.projector_blocks(K_GEN, z, 'BETA'), M.exact_blocks(K_GEN, z))
        check(f'T2 BETA!=exact ({tag})', d > 1e-3, f'max|BETA-exact| = {d:.3e} (need >1e-3)')


def test_T3_matsubara_fabrication_ok():
    gpp, f, fdag, ghh = M.ksum_exact(Z_MATS)
    d1 = np.max(np.abs(ghh - (-np.conj(gpp))))
    d2 = np.max(np.abs(fdag - np.conj(f)))
    check('T3 G_hh==-conj(G_pp) on k-sum (matsubara)', d1 < 1e-12,
          f'max|G_hh+G_pp*| = {d1:.3e} (need <1e-12)')
    check('T3 F_dag==conj(F) on k-sum (matsubara)', d2 < 1e-12,
          f'max|F_dag-F*| = {d2:.3e} (need <1e-12)')


def test_T3b_perk_violates():
    gpp, f, fdag, ghh = M.exact_blocks(K_GEN, Z_MATS)
    d = max(np.max(np.abs(ghh + np.conj(gpp))), np.max(np.abs(fdag - np.conj(f))))
    check('T3b per-k violates the relation', d > 1e-3,
          f'max per-k violation = {d:.3e} (so the check MUST stay outside the k-loop)')


def test_T3c_no_transpose():
    """F_dag is the ELEMENTWISE conjugate of F, not F^dagger.

    dmftgk.py used to assemble the local Nambu matrix with conj(flc.T); that is only
    the same thing when F happens to be symmetric (which it is when sigind is
    diagonal, e.g. the production plaquette K basis, but not in general).
    """
    gpp, f, fdag, ghh = M.ksum_exact(Z_MATS)
    asym = np.max(np.abs(f - f.T))
    d_el = np.max(np.abs(fdag - np.conj(f)))
    d_tr = np.max(np.abs(fdag - np.conj(f.T)))
    check('T3c F is not symmetric in general', asym > 1e-3, f'max|F-F.T| = {asym:.3e}')
    check('T3c F_dag == elementwise conj(F)', d_el < 1e-12, f'max|F_dag-conj(F)| = {d_el:.3e}')
    check('T3c F_dag != conj(F.T)', d_tr > 1e-3, f'max|F_dag-conj(F.T)| = {d_tr:.3e}')


def test_T4_realaxis_fabrication_wrong():
    w, eta = 0.35, 0.02
    gpp, f, fdag, ghh = M.ksum_exact(w + 1j * eta)
    gpp_m, f_m, _, _ = M.ksum_exact(-w + 1j * eta)
    d1 = np.max(np.abs(ghh + np.conj(gpp)))
    d2 = np.max(np.abs(fdag - np.conj(f)))
    check('T4 fabrication WRONG on real axis', max(d1, d2) > 1e-3,
          f'max|G_hh+G_pp*(w)| = {d1:.3e}, max|F_dag-F*(w)| = {d2:.3e} (need >1e-3)')
    r1 = np.max(np.abs(ghh - (-np.conj(gpp_m))))
    r2 = np.max(np.abs(fdag - np.conj(f_m)))
    check('T4 reflected rule G_hh(w)==-conj(G_pp(-w))', r1 < 1e-10,
          f'max|G_hh(w)+G_pp*(-w)| = {r1:.3e} (need <1e-10)')
    check('T4 reflected rule F_dag(w)==conj(F(-w))', r2 < 1e-10,
          f'max|F_dag(w)-F*(-w)| = {r2:.3e} (need <1e-10)')


def test_T5_plaquette():
    """Production geometry: 2x2 plaquette R -> cluster momentum K basis."""
    R = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    Kc = np.array([[0, 0], [np.pi, 0], [0, np.pi], [np.pi, np.pi]], dtype=float)
    W = np.exp(-1j * (Kc @ R.T)) / 2.0
    imw = np.max(np.abs(W.imag))
    check('T5a W real', imw < 1e-12, f'max|Im W| = {imw:.3e}')
    du = np.max(np.abs(W @ W.conj().T - np.eye(4)))
    check('T5b W unitary', du < 1e-14, f'max|W W^dag - I| = {du:.3e}')

    n = 6                                   # 4 plaquette sites + 2 uncorrelated bands
    Wfull = np.zeros((n, n), dtype=complex)
    Wfull[:4, :4] = W
    Wfull[4:, 4:] = np.eye(n - 4)
    Mk = Model(n=n, nc=4, seed=31415, rot=Wfull.real)
    d = np.max(np.abs(Mk.Hk(-K_GEN) - np.conj(Mk.Hk(K_GEN))))
    check('T5c rotated model keeps H(-k)==conj(H(k))', d < 1e-13, f'max diff = {d:.3e}')
    da = maxdiff(Mk.projector_blocks(K_GEN, Z_MATS, CONV_UNDER_TEST), Mk.exact_blocks(K_GEN, Z_MATS))
    check(f'T5d {CONV_UNDER_TEST}==exact in K basis', da < 1e-12,
          f'max|{CONV_UNDER_TEST}-exact| = {da:.3e} (need <1e-12)')
    db = maxdiff(Mk.projector_blocks(K_GEN, Z_MATS, 'BETA'), Mk.exact_blocks(K_GEN, Z_MATS))
    check('T5e BETA!=exact in K basis', db > 1e-3, f'max|BETA-exact| = {db:.3e} (need >1e-3)')


def main():
    print(f'convention under test: {CONV_UNDER_TEST}')
    test_T0_model()
    test_T1_alpha_exact()
    test_T2_beta_negative_control()
    test_T3_matsubara_fabrication_ok()
    test_T3b_perk_violates()
    test_T3c_no_transpose()
    test_T4_realaxis_fabrication_wrong()
    test_T5_plaquette()
    if FAILURES:
        print('\n%d CHECK(S) FAILED:' % len(FAILURES), file=sys.stderr)
        for f in FAILURES:
            print('  ' + f, file=sys.stderr)
        return 1
    print('\nAll checks passed.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
