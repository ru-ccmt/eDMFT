from pathlib import Path
root=Path(__file__).parent/'stage'
p=root/'frozen_sigma_response.py';s=p.read_text()
s=s.replace('self, metadata, energies, raw, normalized, kw, window, *, parallel=SERIAL','self, metadata, energies, raw, normalized, kw, window, *, parallel=SERIAL,\n        diagnostic_gap_cutoff=None')
a='''        if self.minimum_gap < 1e-8:
            raise ValueError(f"window cuts degeneracy: {self.minimum_gap} Ry")'''
b='''        # Strict production behavior is unchanged. This opt-in only measures
        # the sensitivity of the existing hard window to omitting small-gap
        # cross-window response blocks. It neither smooths the projector nor
        # defines a unique derivative of a degenerate spectral cutoff.
        if diagnostic_gap_cutoff is not None and (
            not np.isfinite(diagnostic_gap_cutoff) or diagnostic_gap_cutoff <= 0
        ):
            raise ValueError("diagnostic gap cutoff must be finite and positive")
        self.diagnostic_gap_cutoff = diagnostic_gap_cutoff
        self.diagnostic_only = diagnostic_gap_cutoff is not None
        if self.minimum_gap < 1e-8 and not self.diagnostic_only:
            raise ValueError(f"window cuts degeneracy: {self.minimum_gap} Ry")

    def window_denominators(self):
        """One reciprocal-gap policy shared by forward and transpose actions.

        Infinite denominators implement a zero inverse on the selected blocks.
        Central U, Matsubara weights, k weights and direct density terms remain
        unchanged. A diagnostic result must never be published as validated
        stress, even if the threshold happens not to omit any blocks.
        """
        result = []
        for e in self.eps:
            d = e[self.W][None, :] - e[self.Q][:, None]
            if self.diagnostic_only:
                d = np.where(abs(d) <= self.diagnostic_gap_cutoff, np.inf, d)
            result.append(d)
        return result'''
assert a in s;s=s.replace(a,b)
s=s.replace('for k, e in enumerate(self.eps):\n            a = dh[k][np.ix_(Q, W)] / (e[W][None, :] - e[Q][:, None])','for k, denominator in enumerate(self.window_denominators()):\n            a = dh[k][np.ix_(Q, W)] / denominator')
p.write_text(s)
p=root/'frozen_sigma_tangent.py';s=p.read_text();a='self._denominators = [e[self.W][None, :] - e[self.Q][:, None] for e in self.eps]';assert a in s;s=s.replace(a,'self._denominators = self.window_denominators()');p.write_text(s)
