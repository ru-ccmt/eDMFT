from pathlib import Path
p=Path(__file__).parent/'stage'
f=p/'lapw_complement_response.py';s=f.read_text().replace('''        slope_tolerance=1e-12
''','''        slope_tolerance=1e-12,
        active_columns=None,
''',1)
needle='''        self.minimum_gap = low - float(np.max(e))
''';s=s.replace(needle,'''        # Mixing two thermally empty states outside the DMFT window has no
        # density/projector readout. Its eigenvector derivative is undefined
        # at a degeneracy, but the physical derivative is regular. The caller
        # identifies columns irrelevant to the projector; independently check
        # their thermal occupation and slope before excluding any solve.
        self.active_columns = np.ones(len(e), bool) if active_columns is None else np.asarray(active_columns)
        if self.active_columns.shape != e.shape or self.active_columns.dtype != bool or not np.any(self.active_columns):
            raise ValueError("invalid active complement columns")
        inactive = ~self.active_columns
        f_inactive = expit((self.mu-e[inactive])/self.temperature)
        self.inactive_charge_bound = float(2*np.sum(f_inactive))
        if self.inactive_charge_bound > charge_tolerance or self.inactive_charge_bound/self.temperature > slope_tolerance:
            raise ValueError("inactive retained columns are not thermally empty")
        self.minimum_gap_all = low - float(np.max(e))
        self.minimum_gap = low - float(np.max(e[self.active_columns]))
''',1)
s=s.replace('''        self.diagonal = d[:, None] - e[None, :]
''','''        safe_energies = np.where(self.active_columns, e, np.min(e[self.active_columns]))
        self.diagonal = d[:, None] - safe_energies[None, :]
''',1).replace('''        x = -(self.rotation.conj().T @ rhs)
''','''        x = -(self.rotation.conj().T @ (rhs * self.active_columns))
''',1).replace('''            np.linalg.norm(self.h @ tangent - tangent * self.energies + rhs)
            / max(np.linalg.norm(rhs), np.finfo(float).tiny)
''','''            np.linalg.norm((self.h @ tangent - tangent * self.energies + rhs) * self.active_columns)
            / max(np.linalg.norm(rhs * self.active_columns), np.finfo(float).tiny)
''',1).replace('''            minimum_gap_Ry=self.minimum_gap,
''','''            minimum_gap_Ry=self.minimum_gap,
            retained_endpoint_gap_Ry=self.minimum_gap_all,
            inactive_retained_columns=int(np.sum(~self.active_columns)),
            inactive_retained_charge_bound=self.inactive_charge_bound,
''',1)
# A constructed complement must never deactivate correlated columns, including
# unoccupied ones whose self-energy/projector derivative is still required.
needle='''        self.rho = model.central_density()
''';s=s.replace(needle,'''        for solver in self.solvers:
            if not np.all(solver.active_columns[model.W]):
                raise ValueError("DMFT window columns cannot be omitted from complement response")
'''+needle,1);f.write_text(s)
f=p/'lapw_valence_state.py';s=f.read_text().replace('from scipy.special import spherical_jn','from scipy.special import spherical_jn, expit');needle='''            temperature=model.meta.temperature,
        )
        bases.append(base)''';assert needle in s;s=s.replace(needle,'''            temperature=model.meta.temperature,
            active_columns=(np.isin(np.arange(nret), model.W) |
                (2*nret*expit((model.meta.mu-ev)/model.meta.temperature)
                 > min(1e-12, model.meta.temperature*1e-12))),
        )
        bases.append(base)''',1);f.write_text(s)
f=p/'lapw_total_stress.py';s=f.read_text();needle='''        residual = state.metadata["central_fixed_point_residual_norm"]
''';s=s.replace(needle,'''        if any(getattr(v.model, "diagnostic_only", False) for v in getattr(state.valence, "branches", [state.valence])):
            raise NotImplementedError("window-omission diagnostic cannot publish total stress")
'''+needle,1);f.write_text(s)
