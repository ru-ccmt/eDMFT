# Complete WIEN2k-24 stress and response patch bundle

Updated 9 November 2026. The supported installer now reproduces the **current
36 audited Fortran sources** in LAPW0, LAPW1, LAPW2, LCORE and MIXER, plus the
two necessary Makefile object-list additions. All prerequisites of the old
supplements are included. No full licensed WIEN distribution is included.

This bundle packages existing, tested local source changes. It does not alter
the physics or rebuild the currently installed executables by itself. The
Python eDMFT screened-response engine is distributed separately in
`src/python/`, with native support in `src/dmft2/atpar/`.

## Install or upgrade

Use a **complete private WIEN2k-24 source/build tree**. The system installation
previously inspected on rupc had incomplete LCORE/MIXER sources. Do not replace
site compiler/MPI/library settings with the laptop Makefiles.

```bash
source ~/wmodules                         # on rupc
./apply_wien2k24_stress_patches.sh --check /path/to/private/wien2k_24
./apply_wien2k24_stress_patches.sh --apply /path/to/private/wien2k_24
./verify_wien2k24_stress_patches.sh /path/to/private/wien2k_24
```

The default is `--check`: it verifies the package, recognizes all source hashes,
applies candidate changes in a temporary tree, and checks every resulting hash
without modifying WIENROOT. `--apply` performs the same complete preflight,
creates timestamped backups and a transaction manifest, then replaces the
reviewed files. It refuses unknown source contents and source symlinks.

Recognition is **per file**, so documented mixtures of pristine, legacy
required/full, pre-cache, cached, and supplement/audit states can be upgraded
without guessing the historical application order. An unknown edit aborts the
entire preflight before any production source changes. Concurrent edits detected
after preflight also abort. Failed replacements are rolled back from backups.
Do not force a raw patch onto a tree rejected by the installer.

Only `gga_metric_tensor.o` in LAPW0 and `val_rela_kinstrmt.o` in LAPW2 are
inserted into the actual `OBJS` assignments. Comments, compiler flags, MPI
configuration, library paths, and other Makefile bytes are preserved. Existing
objects are not duplicated. Both objects are verified even when all source
hashes already match.

Optional commands, after a source check:

```bash
./apply_wien2k24_stress_patches.sh --apply --build /path/to/private/wien2k_24
./apply_wien2k24_stress_patches.sh --install /path/to/private/wien2k_24
./apply_wien2k24_stress_patches.sh --install-helpers ~/dbin /path/to/private/wien2k_24
```

Builds use sequential `make -j1` because upstream module dependencies are
incomplete. LAPW1 and LAPW2 real/complex targets are built in order. This does
not restrict MPI execution of the programs. Cross-compiler and cluster builds
still need validation on the destination machine. Never copy macOS binaries
or `.so` libraries onto Linux.

`WIEN_STRESS_BACKUP_DIR` can select the parent of the unique backup directory.
`PYTHON` can select the Python 3 interpreter. The installer uses the standard
library and the system `patch` program, with zero context fuzz for recognized
per-file transitions. The default-off native debug controls remain included;
`--with-debug` is accepted for older invocations.

## Corrections now included

| Program | Current corrections and capabilities |
|---|---|
| LAPW0 | Built-in PBE/PBEsol stress derivatives; complete GGA metric tensor; signed/cubic harmonic and precision fixes; correct spin/atom indexing; FFTW interface hardening; earlier allocation/declaration prerequisites |
| LAPW1 | Fixed-G diagnostic; complete explicit `FIXED_NBANDS` parsing and serial/parallel selection; matching radial data; native H/O exports and radial-response probe controls |
| LAPW2 | Final G1 product-rule and LO multiplicity corrections; deterministic radial derivatives; G2 LO-index/angular conjugation/intermediate-channel fixes; safe bounds; RELA corrections/cache; complete local-orbital radial diagnostic export |
| LCORE | Correct integration boundary, normalization, harmonic handling and precision; fresh RELA stress output |
| MIXER | Correct total stress sign and assembly; existing debug output |

**The earlier G1 half-coefficient correction was superseded.** The current
`val_kinstrmt.F` uses the audited `2*R1-R0` product rule and no extra LO--LO
multiplicity. Historical documentation claiming the earlier half coefficient
was correct must not be used as the current implementation specification.
`WIEN_STRESS_G1_PREAUDIT` exists only to reproduce the earlier defects.

`FIXED_NBANDS` is opt-in and separate from historical `nemax_it` behavior. The
bundle contains all four necessary LAPW1 source changes; merely applying the
old response-probe supplement did not supply the declaration/parser/solver
prerequisites. Source hashes pin the complete implementation.

Radial diagnostics, H/O dumps and fixed-G controls remain dormant unless
requested. The installed eDMFT response implementation is needed for the full
screened stress; this WIEN bundle alone is not a claim of complete DFT/eDMFT
stress accuracy for every cell, occupation policy or strain component.

## Additional compatibility changes

Two portable FFTW C-interface changes in `SRC_3ddens/3ddens.F` and
`SRC_hf/hf.F` are included as an optional group:

```bash
./apply_wien2k24_stress_patches.sh --check --with-fftw-compat /path/to/WIENROOT
./apply_wien2k24_stress_patches.sh --apply --with-fftw-compat /path/to/WIENROOT
./verify_wien2k24_stress_patches.sh --with-fftw-compat /path/to/WIENROOT
```

With that option, `--build`/`--install` additionally handle 3ddens and real/complex
HF. The default stress installation does not require those source directories.
LAPW0's own FFTW changes are already part of the mandatory set.

Local `ssh`/`rsh`, tty and shell-wrapper customizations and compiler Makefiles
are site settings, not portable stress fixes. A blank-line-only change in
`SRC_mixer/AtomsDiag.f` is also excluded. These exclusions are recorded in
`COVERAGE.json`; no functional stress correction is excluded for those reasons.

## Bundle contents and history

- `sources.json`: current hashes and reviewed older per-file versions.
- `TARGET_SOURCES.sha256`: exact mandatory target source hashes (36 files).
- `OPTIONAL_FFTW_TARGET_SOURCES.sha256`: the two optional target hashes.
- `BASE_SOURCES.sha256`: pristine sources; the two new sources have no base file.
- `MANIFEST.sha256`: all distributable bundle files except itself and caches.
- `patches/current/`: reviewable source-only cumulative upgrades from five
  historical states. Makefile object additions are handled by the installer.
- `patches/files/`: the exact per-file transitions used by the installer.
- `patches/build/`: reviewable pristine Makefile object-list patch.
- `tests/test_bundle.py`: reproducibility, idempotence, preservation and refusal
  tests against a supplied pristine tree; no installed source is changed.
- `VALIDATION.txt`: results of the local packaging validation.
- `history/` and `patches/history/20261105/`: earlier installer/documentation and
  cumulative patches, retained for provenance rather than deployment.

The old `patches/required`, `patches/optional/010*` and `supplements/` are retained
as historical evidence. **Do not apply them after the current installer.** Their
changes, including missing prerequisites, are already incorporated. The old
root `*_legacy.sh` entry points now delegate to the current installer so they
cannot accidentally install the superseded G1 implementation.

Historical numerical reports remain under `src/doc/stress/`. They refer to
specific geometries, density snapshots, cutoff values and stages of the audit;
packaging the sources does not create new physical validation results.
