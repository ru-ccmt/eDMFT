# Cluster Hubbard interaction (`Uc`)

This document describes the scalar cluster-Hubbard interaction shared by the
BCTQMC and iNORG impurity wrappers.

## Parameter contract

For a cluster impurity with exactly one active orbital on every site, specifying

```python
iparams0 = {
    # ...
    "cluster": [True, "# cluster impurity"],
    "U":  [9.0, "# double-counting/default Slater U"],
    "Uc": [9.0, "# cluster-Hubbard interaction"],
    "J":  [1.0, "# double-counting/default Slater J"],
    "CoulombF": ["'Ising'", "# density-density interaction"],
}
```

selects the physical interaction

\[
H_{\mathrm{int}}=U_c\sum_R n_{R\uparrow}n_{R\downarrow}.
\]

`Uc` is the final interaction in the active one-band model. It is not a Slater
integral and is never added to `U`. In this mode `J` does not enter the impurity
Hamiltonian; the existing `U`/`J` (and optional `Ud`/`Jd`) double-counting
behavior is unchanged.

`Uc` may be any finite real scalar, including zero or a negative value. The
scalar interface intentionally rejects multi-orbital clusters because a single
number does not specify their inter-orbital, exchange, and pair-hopping terms.

## Compatibility rules

- BCTQMC cluster input with `Uc` uses `Uc` to generate the cluster cix.
- BCTQMC cluster input without `Uc` falls back to legacy `Ud` with a warning.
- BCTQMC cluster input with neither `Uc` nor legacy `Ud` fails instead of
  silently reusing a stale cix.
- iNORG cluster input without `Uc` retains the existing Slater `U`/`J` tensor.
- If both `Uc` and `Ud` are present, `Uc` controls the solver Hamiltonian and
  `Ud` retains its existing double-counting override meaning.
- `Uc` is only accepted for a cluster impurity. Supplying it for a single-site
  impurity is an error rather than a silently ignored parameter.
- `Ising`, `IsingS`, and `IsingK` are accepted and serialized as the equivalent
  density-density one-band Hubbard interaction. `Full` is rejected in `Uc`
  mode; the legacy generic Full-tensor path is not changed.

## BCTQMC CIX updates and restart

`UpdateAtom` retains its legacy responsibility: it requests regeneration of the
cluster CIX so that dynamic impurity levels can be updated.  Regenerating a CIX
does not by itself remove `status.*`, `sampled_data`, or `PPSigma.0`; compatible
Bold-CTQMC warm-start data therefore survives ordinary DMFT iterations with
`UpdateAtom=1`.

Restart files are cleared only when the interaction/CIX cache is missing or its
sidecar fingerprint does not match the requested `Uc`, cluster shape, mode,
basis, or `Sigind`.  This separates a routine CIX update from a physical or
structural change of the impurity problem.

## Four-site plaquette

For `Uc=9`, the site-basis tensor has only

\[
U_{aaaa}=9,\qquad a=0,1,2,3,
\]

with every other entry zero. Consequently

```text
Vij_abba = diag(9, 9, 9, 9)
Vij_abab = 0
```

and the plaquette K-basis vertex has magnitude

\[
U_c/N_{\mathrm{site}}=9/4=2.25
\]

for every momentum-conserving term. BCTQMC and iNORG use different raw spin
ordering and Nambu conventions, so raw line counts or signs should only be
compared after canonicalizing those conventions. The physical site-basis
Hamiltonian is identical.

By contrast, legacy iNORG input `U=9,J=1` without `Uc` uses Slater integrals.
For an active d-x2-y2 orbital its intra-orbital matrix element is

\[
9+\frac{8}{7}=\frac{71}{7}=10.142857\ldots,
\]

which is intentionally preserved for backward compatibility.
