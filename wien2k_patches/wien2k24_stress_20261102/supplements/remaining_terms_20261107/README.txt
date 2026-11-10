LAPW0 / LCORE REMAINING-TERMS AUDIT SUPPLEMENT, 2026-11-07

Contains the signed-harmonic, full GGA tensor, and electrostatic precision fixes.
Read the detailed report under src/doc/stress/remaining_terms_audit_20261107/.

Check applicability:
  python3 apply_remaining_stress.py /path/to/WIENROOT
Apply after a successful check:
  python3 apply_remaining_stress.py /path/to/WIENROOT --apply

The script recognizes exact pre-audit Fortran source hashes, verifies the patch
in a temporary tree, and creates timestamped backups. It preserves site compiler
and library settings in Makefile, adding only gga_metric_tensor.o. It refuses
unknown or partially modified sources. It performs no build or installation.
Build SRC_lapw0 and SRC_lcore with the target site's normal settings, then install
the resulting executables. Use a sequential make when changing serial/MPI mode;
the upstream Makefile has incomplete module dependencies for a clean parallel
make (this is separate from MPI execution, which passed the diagnostic runs).

This supplements the older cumulative bundle. Its top-level installer/target
manifest are not rewritten here and do not incorporate all intervening G1/G2
and fixed-band-count audits. Do not interpret that older manifest as a full
inventory of the current laptop sources. No remote executable was changed.

Normal GGA stress uses the full tensor. The explicit diagnostic environment
WIEN_STRESS_GGA_MODE=LEGACY restores the old GGA contraction for comparison.
Core/electrostatic harmonic and precision corrections remain active in that mode.
The ordinary potential and core density are unchanged by these corrections.
The principal SrVO3 eDMFT stress mismatch remains unresolved; see the report.
