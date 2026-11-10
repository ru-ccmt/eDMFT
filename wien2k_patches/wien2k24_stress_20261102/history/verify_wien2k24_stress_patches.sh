#!/usr/bin/env bash
set -eu

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 WIENROOT" >&2
  exit 2
fi

wienroot=${1%/}
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

echo 'Checking exact patched-source checksums...'
(cd "$wienroot" && shasum -a 256 -c "$script_dir/TARGET_SOURCES.sha256")

check() {
  description=$1
  shift
  if "$@"; then
    printf 'OK   %s\n' "$description"
  else
    printf 'FAIL %s\n' "$description" >&2
    exit 1
  fi
}

check 'LAPW0 built-in PBE stress derivative' \
  grep -Fq 'call pbe_vsigma_builtin' "$wienroot/SRC_lapw0/vxclm2.F"
check 'LAPW0 RELA density stress is freshly written' \
  grep -Fq 'allowed mixer to reuse stale :STR_DENS' "$wienroot/SRC_lapw0/finl_elect_str.F"
check 'LCORE RELA core stress is freshly written' \
  grep -Fq 'RELA mode caused mixer to retain stale NREL data' "$wienroot/SRC_lcore/hfsd.f"
check 'MIXER physical stress sign' \
  grep -Eq 'final_stress\(ialpha,ibeta\)[[:space:]]*=[[:space:]]*-[[:space:]]*Tstress' "$wienroot/SRC_mixer/mixer.F"
check 'LAPW1 fixed-G diagnostic' \
  grep -Fq 'LAPW1_FIXGVEC_V1' "$wienroot/SRC_lapw1/coors.f"
check 'LAPW1 nonspherical radial matrices' \
  grep -Fq 'WPUU' "$wienroot/SRC_lapw1/forfhs.f"
check 'LAPW2 scalar-relativistic source is compiled' \
  grep -Fq 'val_rela_kinstrmt.o' "$wienroot/SRC_lapw2/Makefile"
check 'LAPW2 scalar-relativistic production mode' \
  grep -Fq 'WIEN_STRESS_RELA_MODE' "$wienroot/SRC_lapw2/symmetrization_stress.f"
check 'LAPW2 safe G2 trace mode' \
  grep -Fq 'WIEN_STRESS_G2_MODE' "$wienroot/SRC_lapw2/symmetrization_stress.f"
check 'LAPW2 RELA band OpenMP loop' \
  grep -Fq '!$omp do schedule(dynamic,1)' "$wienroot/SRC_lapw2/val_rela_kinstrmt.F"
check 'LAPW2 atom-dependent RELA operator cache' \
  grep -Fq 'module rela_operator_cache' "$wienroot/SRC_lapw2/val_rela_kinstrmt.F"
check 'LAPW2 cache tracks the radial quadrature rule' \
  grep -Fq 'use TestCases, only: SimpInt,ConOpt' "$wienroot/SRC_lapw2/val_rela_kinstrmt.F"
check 'LAPW2 matching-coefficient BLAS contraction' \
  grep -Fq "CALL zgemm('N','N',maxindex" "$wienroot/SRC_lapw2/val_corr_stress.F"
check 'LAPW2 interstitial BLAS convolution' \
  grep -Fq "call zgemm('T','N',isize" "$wienroot/SRC_lapw2/val_kin_interstitial.F"
check 'LAPW2 radial-response diagnostic' \
  grep -Fq 'WIEN_STRESS_RADIAL_DUMP' "$wienroot/SRC_lapw2/atpar.F"

env PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/wien_stress_pycache \
  python3 -m py_compile \
  "$script_dir/helpers/rmt_stress_correction.py" \
  "$script_dir/helpers/radial_response_fd.py"
echo 'OK   bundled Python diagnostics compile'
echo 'Complete WIEN2k-24 stress patch verified.'

