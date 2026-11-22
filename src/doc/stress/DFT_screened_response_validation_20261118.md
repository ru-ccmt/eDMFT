# Final screened DFT stress in LAPW: implementation and validation

18 November 2026. The values below are the **mean tensile stress**
`tr(sigma)/3` in GPa, matching the sign of the current WIEN `:PRESSURE:` row.
The thermodynamic compression pressure has the opposite sign.

## The missing basis variable

WIEN LAPW1 replaces an energy parameter of `0.30` in `case.in1[c]` by
`E_ref - 0.20` or `E_ref + 0.20` Ry. At convergence `E_ref` is the Fermi
energy. This means that a strained, self-consistent calculation changes the
radial basis even at fixed Kohn–Sham potential. The original screened stress
held `E_ref` fixed. In SrVO3 the endpoint values of `E_ref` change at roughly
-3.17 Ry per unit hydrostatic lattice strain. Omitting that path left the
screened stress farther from the energy finite difference than native LAPW2.

Let `q = delta n` be the response of the material-coordinate charge and
`e = delta E_ref`. At fixed particle number, define `J` as the previously
implemented charge Jacobian, `p` as the charge response to a unit change in
`E_ref`, `m` as the Fermi-energy response to charge, and `m_e` as the
Fermi-energy response to a unit change in `E_ref`. The explicit strain sources
are `s` for charge and `t` for the Fermi energy. The coupled response is

```text
(I - J) q - p e       = s,
       -m q + (1-m_e)e = t.
```

For the finite-basis Mermin functional, let `a` be its explicit strain
derivative at fixed `n,E_ref`, and let `c,g` be its derivatives with respect
to `n,E_ref`. The full derivative is `D = a + c·q + g e`. The code solves one
adjoint for all symmetry-allowed strain directions:

```text
[(I-J)^T  -m^T] [lambda  ]   [c]
[ -p^T   1-m_e] [lambda_e] = [g],
D = a + lambda·s + lambda_e t.
```

The `0.30` channels are identified from the actual IN1 records, including
APW+lo and repeated local-orbital entries. Searched radial energies retain
their independently differentiated root response. Their slots are excluded
from the `E_ref` channel. The fixed-N Fermi derivative uses the same Mermin
occupations and overlap convention as the charge response.

## Matched material checks

All three structures in each triplet used `TEMPS 0.005`, fixed RMT and
RKMAX 8, and the same stored band count and reciprocal-vector **sets at
every k point**. The endpoint charge densities were reconverged after
switching from TETRA to TEMPS. The numerical result is the central difference
of converged `:ENE` values divided by `3 Omega` for an isotropic `a,b,c`
change of ±0.001. Native stress is from the central LAPW2/MIXER cycle;
screened stress is the final-only response on the same pre-MIXER charge.

| System | Mode | k points | Bands | Native | Screened | Numerical | Screened minus numerical |
|---|---|---:|---:|---:|---:|---:|---:|
| SrVO3 | NREL, cubic | 75 | 43 | 8.59353 | 8.36898 | 8.36772 | +0.00126 |
| FeSe | RELA, tetragonal | 140 | 66 | 12.16906 | 12.05719 | 12.05661 | +0.00057 |

Thus the absolute native-to-numerical discrepancy drops from **0.226 to
0.00126 GPa** for SrVO3 and from **0.112 to 0.000574 GPa** for FeSe. The
native 0.1–0.2 GPa scale is consequential for structural optimization;
the corrected discrepancies are below 0.002 GPa in these matched tests.

The SrVO3 endpoint charge distances were at most 5.5e-6; FeSe's were at
most 6e-7. These comparisons validate the coupled derivative in these two
metallic test cases. They do not establish accuracy for SOC, a moving RMT,
branch-changing radial searches, or an insulating Fermi-level policy.

The response contraction took about 205 s for SrVO3 and 257 s for FeSe with
four MPI ranks and one math-library thread per rank. Those times exclude the
isolated LAPW0/LAPW1/LAPW2/LCORE snapshot preparation.

## Execution and outputs

Set `STR` in `case.in0`, `STRT` (stress) or `STRFT` (stress and forces) in
`case.in2[c]`, and `TEMPS` on its occupation line. With `source ~/activate`,
`run_dft.py` or WIEN `run_lapw` invokes the screened response **once after
charge convergence**. Ordinary LAPW2 evaluates its native STR/STRF terms in
the final cycle. The response operates in a separate
`case.dft_response_runs/<id>/` directory, then appends the final corrected
stress rows to `case.scfm` and `case.scf`. The native rows remain available
above them. `case.output2_response` gives the comparison and a link to the
full JSON report. `mpi_prefix.dat` selects the k-point MPI launch.

An early `run_lapw -e <step>` exit does not launch the response. The WIEN
`run_lapw -f case` path now passes the chosen case to MIXER's definition and
convergence checks. The DFT Python merger no longer imports a stale
`Eorb_imp.dat` from a copied DMFT directory.

Reproducible isolated calculations and reports are under
`/private/tmp/edmft_dft_response/`. The original material directories were
not modified for these finite-difference tests.
