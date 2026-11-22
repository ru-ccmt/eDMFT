# Complete WIEN2k-24 stress and response patch bundle

Updated 21 November 2026 for WIEN2k 24.1. The bundle covers the 36 audited
Fortran sources in LAPW0, LAPW1, LAPW2, LCORE and MIXER, the two Makefile
object-list additions, the final screened DFT stress hook in `run_lapw`, and
the sources and build procedure for four WIEN response executables. Two FFTW
interface changes for 3ddens and HF are available separately.

## Install or upgrade

Start in this bundle directory. Set `WIENROOT` to a **complete, private and
configured WIEN2k 24.1 source/build tree**. It must contain the component
sources, Makefiles and `SRC_Globals`. On RUPC, first run
`source ~/wmodules.sh` to load the compilers and libraries, then set `WIENROOT`
explicitly because that module selects the system installation. Set
`WIEN_DMFT_ROOT` to the directory containing the installed eDMFT programs.
For example, after installing eDMFT in `~/dbin`:

```bash
export WIEN_DMFT_ROOT="$HOME/dbin"
export WIENROOT=/path/to/private/configured/wien2k_24
./apply_wien2k24_stress_patches.sh --check --with-run-lapw \
    --install-response-workers "$WIEN_DMFT_ROOT" "$WIENROOT"
./apply_wien2k24_stress_patches.sh --apply --with-run-lapw --install \
    --install-response-workers "$WIEN_DMFT_ROOT" "$WIENROOT"
./verify_wien2k24_stress_patches.sh --with-run-lapw \
    --install-response-workers "$WIEN_DMFT_ROOT" "$WIENROOT"
```

The middle command updates the sources and driver, rebuilds the ordinary
executables, then builds the four response workers in a clean temporary tree.
The first and third commands are read-only source checks. Verification checks
the proposed worker source hashes. The install step also checks that all four
new binaries exist and are executable before replacing any workers in the
destination. Backups of replaced sources, scripts and workers have their own
transaction records. If a compiler or linker fails, correct the build settings
and repeat the install command.

For ordinary stress alone, omit the driver and worker options. Running the
installer with only the source-tree argument performs a read-only preflight.

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

Builds use `make -j1` because upstream Fortran module dependencies are
incomplete. LAPW1 and LAPW2 real/complex targets are built in order. On sites
with MPI-configured LAPW0 Makefiles, add `--with-lapw0-mpi` to the install
command to build and install `lapw0_mpi` as well. The four response workers
are serial native exporters. The eDMFT Python driver distributes the response
over k points. Build platform-specific executables and libraries on each host.

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
| LAPW2 | Final G1 product-rule and LO multiplicity corrections; deterministic radial derivatives; G2 LO-index/angular conjugation/intermediate-channel fixes; safe bounds; RELA corrections/cache; complete local-orbital radial diagnostic export; STRT and STRFT input modes |
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
- `patches/build/`: reviewable Makefile object-list patch and the separate XC/core response-worker patch.
- `tests/test_bundle.py`: reproducibility, idempotence, preservation and refusal
  tests against a supplied pristine tree; no installed source is changed.
- `VALIDATION.txt`: results of the local packaging validation.
- `history/` and `patches/history/20261105/`: earlier installer/documentation and
  cumulative patches, retained for provenance rather than deployment.

The legacy stress patches and older supplements are retained for provenance.
The current installer applies all selected changes itself. The
`*_legacy.sh` entry points delegate to it.

## Screened response and its dependencies

`--with-run-lapw` installs the final-step hook for the distributed WIEN2k
24.1 script or the previously supported private script. It updates the active
`run_lapw` when present, or the `SRC/run_lapw` template in an unpacked tree.
A site-modified script with an unknown checksum requires manual review against
`supplements/run_lapw_dft_response.patch`. In split installations, deploy the
updated driver and ordinary executables from the build tree to the runtime
tree using the site's usual procedure.

LAPW2 accepts `STRT` and `STRFT` as `STR` and `STRF` during charge
iterations. After convergence, `run_lapw` calls `dft_stress_step.py` once.
It uses `case.clmsum_old` to match the last vectors and potential. The eDMFT
DFT driver `run_dft.py` also supports this final step. Use `TEMPS` smearing
for a comparison to the derivative of the stationary Mermin free energy.

`--install-response-workers DIR` builds the following executables with the
configured native compiler. The XC/core exporter changes are applied **only
in the temporary worker build**.

| Executable installed in DIR | Data used by screened response |
|---|---|
| `lapw0_response` | Native XC potential, muffin-tin XC quadrature, and interstitial FFT data |
| `lcore_response` | Core radial mesh, potential and large/small components |
| `lapw1_response`, `lapw1c_response` | Real/complex LAPW1 Hamiltonian and overlap exports |

The source patch is `patches/build/response-workers.patch`. These exports were
originally documented under `src/doc/stress/xc_core_response_20261108/`.
That development archive is no longer needed to build the WIEN workers.

eDMFT installs its response drivers, `liblapw_response.so`, and `dmft2`.
For screened eDMFT calculations, the driver uses the installed `dmft2` for
both the ordinary charge calculation and the response export, with or without
spin–orbit coupling. The installation commands above supply the four additional
WIEN workers. Users testing a different eDMFT executable can set `solver.dmft2`
in the case response JSON to its path.

The optional `--install-helpers` installs two standalone radial diagnostics.

In matched DFT finite-difference tests, the mean-stress discrepancy decreased
from 0.226 to 0.00126 GPa for cubic SrVO3 (NREL) and from 0.112 to
0.000574 GPa for tetragonal FeSe (RELA). The response includes the
Fermi-level-relative IN1 radial-energy path. These are tested cases, not a
general error bound. See `src/doc/stress/DFT_screened_response_validation_20261118.md`
in the local development tree for inputs and response reports.
