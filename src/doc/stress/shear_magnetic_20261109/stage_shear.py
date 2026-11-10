from pathlib import Path
r=Path(__file__).parent/'stage'
p=r/'lapw_total_stress.py';s=p.read_text()
s=s.replace('Complete native frozen-Sigma stress for symmetry-enforced diagonal tensors.', 'Complete native frozen-Sigma stress in the full crystal-invariant tensor space.')
a=s.index('Supported publication:');b=s.index('"""',a)
s=s[:a]+'''The six Cartesian symmetric components are reduced only by the native charge
symmetry. Triclinic cells retain all six; monoclinic cells retain their allowed
shear. Non-invariant perturbations would need expanded density/star layouts,
but their first derivatives vanish for a scalar with the imposed symmetry.
The report explicitly identifies projection onto that imposed symmetry when
the supplied self-energy is not invariant. RMT and local frames remain fixed.
SOC and spin-dependent KS potentials are outside this implementation.
'''+s[b:]
s=s.replace('def diagonal_invariant_strains(structure):','def invariant_strains(structure):',1)
a=s.index('    """A complete strain basis');b=s.index('    s = structure',a)
s=s[:a]+'''    """Deterministic complete basis of crystal-invariant symmetric strains.

    Project all six Cartesian directions with the point-group Reynolds map.
    Keep hydrostatic strain first, then independent diagonal and shear images.
    Never infer missing stress components from a hydrostatic measurement.
    Fixed-point sources remain in the native symmetry-adapted coordinates:
    every returned e obeys R e R^T=e. For P1/-1 this is all six components.
    """
'''+s[b:]
s=s.replace('if invariant.shape[1] == 0 or np.linalg.norm(invariant[3:]) > 1e-9:\n        raise NotImplementedError("symmetry allows shear; expanded sources required")','if invariant.shape[1] == 0:\n        raise ValueError("point symmetry has no invariant hydrostatic direction")')
s=s.replace('for i in [2, 1, 0]:','for i in [2, 1, 0, 3, 4, 5]:',1).replace('incomplete invariant diagonal strain basis','incomplete invariant symmetric strain basis')
needle='\n\ndef stress_from_contractions'
s=s.replace(needle,'''\n\ndef diagonal_invariant_strains(structure):
    """Compatibility check for callers explicitly requesting a diagonal tensor."""
    directions = invariant_strains(structure)
    if any(np.linalg.norm(e - np.diag(np.diag(e))) > 1e-9 for e in directions):
        raise NotImplementedError("symmetry allows shear; use invariant_strains")
    return directions
'''+needle,1)
s=s.replace('        progress=None,\n    ):','        progress=None,\n        diagnostic=False,\n    ):',1)
a=s.index('        if hasattr(state.valence, "branches"):');b=s.index('        if not np.isfinite(fixed_point_tolerance)',a)
s=s[:a]+'''        branches = getattr(state.valence, "branches", [state.valence])
        omitted = any(getattr(v.model, "diagnostic_only", False) for v in branches)
        if omitted and not diagnostic:
            raise NotImplementedError("window-omission diagnostic cannot publish total stress")
        self.diagnostic_only = bool(diagnostic or omitted)
        residual = state.metadata["central_fixed_point_residual_norm"]
'''+s[b:]
s=s.replace('self.strains = diagonal_invariant_strains(', 'self.strains = invariant_strains(')
a=s.index('        meta = state.valence.model.meta');b=s.index('        self.sigma_symmetry_error',a)
old=s[a:b]
new=old.replace('        meta = state.valence.model.meta\n        errors = []\n','        errors = []\n        for branch in branches:\n            meta = branch.model.meta\n')
# Indent the complete per-spin symmetry test by one level.
rest=new.index('        for sigma in ')
new=new[:rest]+''.join('    '+line if line.strip() else line for line in new[rest:].splitlines(True))
s=s[:a]+new+s[b:]
s=s.replace('            full_tensor_within_symmetry=True,','            full_tensor_within_symmetry=True,\n            diagnostic_only=self.diagnostic_only,\n            stress_ready=not self.diagnostic_only,')
s=s.replace('scope="nonmagnetic fixed-N primitive orthogonal; symmetry-enforced diagonal tensor",','scope="complete symmetric tensor within imposed native charge symmetry; shared unpolarized KS potential",')
p.write_text(s)
p=r/'lapw_strain_field.py';s=p.read_text();a=s.index('"""');b=s.index('"""',a+3)
s=s[:a]+'''"""General symmetry-preserving strain of the native potential map.

All six symmetric Cartesian components are supported, including shear. The
native star/equivalent-atom layout can be reused precisely when e commutes
with its Cartesian point group. Symmetry-forbidden first stress derivatives
are zero; symmetry-breaking response/elastic constants need expanded layouts.
The density-coordinate transport is the independent scalar convention tr(e)/3.
"""'''+s[b+3:]
s=s.replace('Analytic potential source for diagonal cell strain plus scalar transport.','Analytic potential source for symmetric cell strain plus scalar transport.').replace('general diagonal e','general symmetric e').replace('"diagonal strain"','"symmetric strain"')
s += '\n# Historical function name retained for existing importers.\nsymmetric_strain_potential = diagonal_strain_potential\n'
p.write_text(s)
p=r/'lapw_magnetic_response.py';s=p.read_text().replace('fingerprints=up.fingerprints,','fingerprints=up.fingerprints | down.fingerprints,').replace('fingerprint_paths=up.fingerprint_paths,','fingerprint_paths=up.fingerprint_paths | down.fingerprint_paths,');p.write_text(s)
p=r/'stress_response_publish.py';s=p.read_text().replace('report.get("converged") is not True\n        or report.get("full_tensor_within_symmetry") is not True','report.get("converged") is not True\n        or report.get("full_tensor_within_symmetry") is not True\n        or report.get("diagnostic_only", False)');p.write_text(s)
