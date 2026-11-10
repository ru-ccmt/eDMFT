Second-pass reproducibility
===========================
Read REPORT.txt and the accompanying self-contained LaTeX/PDF note first.
No production source was changed by this pass. Existing baseline dirty files
were retained. Scripts extract and compile the current source bodies, so
compare SOURCE_AND_BINARY_SHA256.json before interpreting a future rerun.

New independent tests:
  combined_operator_test.py   actual G1+G2 versus finite full-H/O endpoints
  check_matching.py          actual APW/LAPW/LO derivative versus endpoints
  check_surface.py           actual APW jump matrix versus angular quadrature
Rerun earlier controls:
  check_g1.py                 Cartesian polynomial Hessian, all complex entries
  check_interstitial.py      Fourier step and complex-density contraction
  check_radial.py             positive-mesh radial numerical differentiation

These audit scripts use local absolute source/compiler paths and require
NumPy, SciPy with sph_harm_y, SymPy, gfortran, and macOS Accelerate where used.
They write beside the script. COPY the audit folder to a fresh writable
scratch folder before running; do not run in the documentation directory.
Example (set audit_dir to this archive's absolute directory):
  scratch=$(mktemp -d /tmp/g1g2_repeat.XXXXXX)
  cp "$audit_dir"/*.py "$scratch/"
  cd "$scratch"
  /Users/haule/miniconda3-clean/envs/work311/bin/python combined_operator_test.py
  /Users/haule/miniconda3-clean/envs/work311/bin/python check_matching.py
  /Users/haule/miniconda3-clean/envs/work311/bin/python check_surface.py
  /Users/haule/miniconda3-clean/envs/work311/bin/python check_g1.py
  /Users/haule/miniconda3-clean/envs/work311/bin/python check_interstitial.py
  /Users/haule/miniconda3-clean/envs/work311/bin/python check_radial.py
The intentional pre-audit G1 negative control fails its oracle; default
production must pass. It is not a newly failing production test.

The existing production tests live in src/dmft2/tests/ with README.txt. Their
JSON results from this second pass are under existing/ here. Use a separate
EDMFT_STRESS_TEST_DIR for those routines to avoid module-file collisions.

run_dft.py/run_dmft.py reuse the first-pass /private/tmp fixtures and execute
only the installed current binaries, in new isolated copies. They do not
reconverge endpoints. summarize.py compares to the first-pass material data;
verify_mixer.py checks DFT assembly. These are session-specific harnesses,
not general material setup tools. The fresh small SCF records and original
endpoint reference JSON are preserved even if scratch data are removed.
