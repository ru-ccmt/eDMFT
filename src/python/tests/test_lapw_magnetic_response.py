"""Shared-mu constraint and spin averaging, with independent branch oracles."""
import unittest
import numpy as np
from lapw_magnetic_response import CollinearComplement
from lapw_complement_response import complement_electronic_operator
from test_lapw_complement_response import fixture
from frozen_sigma_tangent import real_inner


class MagneticTests(unittest.TestCase):
    def setUp(self):
        self.up=fixture()[2];self.down=fixture()[2]
        self.combo=CollinearComplement(self.up,self.down)

    def perturbations(self,inputs):
        rng=np.random.default_rng(662)
        return inputs.pack(*(rng.normal(size=b.size) for b in inputs.blocks))

    def number(self,inputs,outputs,x,y):
        counts=[]
        for branch,ci,co,xx,yy in zip(self.combo.branches,inputs.blocks,outputs.blocks,inputs.unpack(x),outputs.unpack(y)):
            _,o,_,_,_=ci.unpack(xx); rr,_=co.unpack(yy)
            counts.append(np.trace(rr,axis1=1,axis2=2).real.sum()+real_inner(branch.rho,o))
        return np.array(counts)

    def test_shared_number_and_exact_transpose(self):
        op,inputs,outputs=self.combo.electronic_operator(fixed_number=True)
        x=self.perturbations(inputs)
        n=self.number(inputs,outputs,x,op.matvec(x))
        self.assertLess(abs(n.sum()),3e-13)
        self.assertGreater(abs(n[0]),1e-3) # Spin populations must be allowed to move.
        self.assertLess(op.check_transpose(samples=2),1e-12)

    def test_fixed_mu_equals_two_independent_branches(self):
        op,inputs,outputs=self.combo.electronic_operator(fixed_number=False)
        x=self.perturbations(inputs)
        expected=[]
        for b,xx in zip(self.combo.branches,inputs.unpack(x)):
            child,_,_=complement_electronic_operator(b,fixed_number=False)
            expected.append(child.matvec(xx))
        np.testing.assert_allclose(op.matvec(x),outputs.pack(*expected),atol=1e-14)
        self.assertLess(op.check_transpose(samples=1),1e-12)

    def test_spin_antisymmetric_shift_survives_fixed_N(self):
        op,inputs,outputs=self.combo.electronic_operator(fixed_number=True)
        directions=[]
        for branch,c,sign in zip(self.combo.branches,inputs.blocks,[1,-1]):
            nk,nr=branch.model.eps.shape
            h=sign*np.broadcast_to(np.eye(nr),(nk,nr,nr))
            q=[np.zeros((s.nq,s.nr),complex) for s in branch.solvers]
            directions.append(c.pack(h,np.zeros_like(h),np.zeros_like(branch.model.raw),q,q))
        x=inputs.pack(*directions);n=self.number(inputs,outputs,x,op.matvec(x))
        np.testing.assert_allclose(n,[-self.up.model.compressibility,self.down.model.compressibility],rtol=2e-13)

    def test_half_weighted_scalar_and_common_mu_guard(self):
        g,c=self.combo.functional_gradients()
        from lapw_strain_functional import lattice_functional_gradients
        child,cc=lattice_functional_gradients(self.up)
        np.testing.assert_allclose(g['lattice'],.5*c.pack(child['lattice'],child['lattice']))
        self.down.model.meta.mu+=1e-3
        with self.assertRaisesRegex(ValueError,'share mu'): CollinearComplement(self.up,self.down)

    def test_unequal_spin_populations_against_nonlinear_common_mu(self):
        """A ferromagnetic fixture must conserve only the combined population.

        Shift the down-spin self-energy, leaving its KS bands and mu shared.
        Independently solve the nonlinear number equation at +/- perturbation;
        this tests the shared-mu tangent without imposing spin compensation.
        """
        from dataclasses import replace
        from scipy.optimize import brentq
        from frozen_sigma_tangent import CachedFrozenSigmaResponse
        from lapw_complement_response import EmptyComplementResponse

        old = self.down.model
        shift = .25 * np.eye(2)[:, :, None]
        meta = replace(old.meta, sigma_inf=old.meta.sigma_inf + shift,
                       sigma=old.meta.sigma + shift[:, :, :, None])
        model = CachedFrozenSigmaResponse(meta, old.eps, old.raw,
                                         old.apply_x(old.raw, old.X), old.kw, old.W)
        down = EmptyComplementResponse(model, self.down.solvers, self.down.raw_q)
        combo = CollinearComplement(self.up, down)
        op, inputs, outputs = combo.electronic_operator(fixed_number=True)
        models = [b.model for b in combo.branches]
        identities = [np.broadcast_to(np.eye(m.eps.shape[1]),
                       m.eps.shape + (m.eps.shape[1],)) for m in models]
        number = lambda r: np.trace(r, axis1=1, axis2=2).real.sum()
        centers = [m.central_density() for m in models]
        self.assertGreater(abs(number(centers[0]) - number(centers[1])), .01)
        target = sum(number(r) for r in centers)
        shifts = [.4, -.2]
        zero_q = [[np.zeros((s.nq, s.nr), complex) for s in b.solvers]
                  for b in combo.branches]
        x = inputs.pack(*(c.pack(a * eye, np.zeros_like(eye),
                                 np.zeros_like(m.raw), q, q)
                          for c, a, eye, m, q in zip(inputs.blocks, shifts,
                                                     identities, models, zero_q)))

        def endpoint(step):
            def densities(dmu):
                return [m.density((step * a - dmu) * eye, 1., m.uc)
                        for m, a, eye in zip(models, shifts, identities)]
            dmu = brentq(lambda v: sum(number(r) for r in densities(v)) - target,
                         -.01, .01, xtol=1e-14)
            return densities(dmu)

        step = 1e-5
        plus, minus = endpoint(step), endpoint(-step)
        expected = outputs.pack(*(c.pack((p - n) / (2 * step), q)
                                  for c, p, n, q in zip(outputs.blocks, plus,
                                                       minus, zero_q)))
        np.testing.assert_allclose(op.matvec(x), expected, atol=2e-7, rtol=2e-6)

if __name__=='__main__':unittest.main()
