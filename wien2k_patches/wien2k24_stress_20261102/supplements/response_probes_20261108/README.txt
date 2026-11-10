Optional response probes, 8 November 2026

Adds native LAPW1 H/O dumps, low-l controls, and all-angular restoration for
an alternate radial basis. Ordinary execution is unchanged with flags absent.
See src/doc/stress/radial_all_l_20261108/REPORT.txt in the eDMFT repository for
results, conventions, limitations, and the unresolved RELA radial comparison.

This is an INCREMENTAL supplement to the exact development sources recorded
in source_changes.json. It is not a replacement for the cumulative stress
patch. The baseline includes earlier fixed-band/debug changes; verify hashes
before applying, and do not force a patch onto an unknown source version.
No compiler settings or Makefiles are changed. Build real and complex LAPW1
with the destination site's normal toolchain after applying it. The native
matrix dump currently requires a serial WIEN LAPW1 build.

The combined patch was applied to a disposable exact-baseline tree and all
resulting hashes verified. No remote deployment was performed in this pass.
