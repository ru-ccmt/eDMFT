#!/usr/bin/env bash
set -eu

if [ "${#}" -ne 1 ]; then
  echo "Usage: verify_wien2k24_stress_patches.sh WIENROOT" >&2
  exit 2
fi

wienroot=${1%/}
fail=0

check_present() {
  label=${1}
  pattern=${2}
  file=${3}
  if grep -Eq "${pattern}" "${file}"; then
    echo "PASS  ${label}"
  else
    echo "FAIL  ${label}" >&2
    fail=1
  fi
}

check_absent() {
  label=${1}
  pattern=${2}
  file=${3}
  if grep -Eq "${pattern}" "${file}"; then
    echo "FAIL  ${label}" >&2
    fail=1
  else
    echo "PASS  ${label}"
  fi
}

lapw0=${wienroot}/SRC_lapw0/vxclm2.F
lapw0_elect=${wienroot}/SRC_lapw0/finl_elect_str.F
lapw0_gga=${wienroot}/SRC_lapw0/stress_gga.F
lapw0_xc=${wienroot}/SRC_lapw0/xcpot3.F
lcore=${wienroot}/SRC_lcore/core_corr_stress.f
mixer=${wienroot}/SRC_mixer/mixer.F
lapw2=${wienroot}/SRC_lapw2/val_kinstrmt.F

check_present "lapw0 calls built-in PBE VSIGMA" "call[[:space:]]+pbe_vsigma_builtin" "${lapw0}"
check_present "lapw0 contains built-in PBE VSIGMA routine" "subroutine[[:space:]]+pbe_vsigma_builtin" "${lapw0}"
check_present "lcore stops core-correction integral at R_MT" "counter[[:space:]]*=[[:space:]]*nrp[[:space:]]*$" "${lcore}"
check_absent "lcore no longer integrates nrp+141" "counter[[:space:]]*=[[:space:]]*nrp[[:space:]]*\+[[:space:]]*141" "${lcore}"
check_present "lcore accumulates one normalized total" "sum_tot\(index2,jatom\)[[:space:]]*=.*cor_corr_tot\(index2\)" "${lcore}"
check_present "mixer converts Tstress to physical stress" "final_stress\(ialpha,ibeta\)[[:space:]]*=[[:space:]]*-[[:space:]]*Tstress\(index1\)" "${mixer}"
if grep -Fq "WIEN2k stress Eq. (G1) F3 factor correction" "${lapw2}" || \
   grep -Eq "g1_deriv_scale[[:space:]]*=[[:space:]]*0\.5d0" "${lapw2}"; then
  echo "PASS  lapw2 marks the Eq. (G1) F3 factor correction"
else
  echo "FAIL  lapw2 does not mark the Eq. (G1) F3 factor correction" >&2
  fail=1
fi
check_present "lapw0 preserves cubic l=8/10 c3" "c3[[:space:]]*=[[:space:]]*c_kub\(IABS\(lm1\(1,i\)\),lm1\(2,i\)\+8[[:space:]]*\)" "${lapw0_elect}"
check_present "lcore preserves cubic l=8/10 c3" "c3[[:space:]]*=[[:space:]]*c_kub\(IABS\(lm1\(1,i\)\),lm1\(2,i\)\+8[[:space:]]*\)" "${lcore}"
check_present "lapw0 spin-down XC stress uses spin-down Vxc" "RHOKSP\(J,2\)=VKEXC\(J,2\)-VKVXC\(J,2\)" "${lapw0_xc}"
check_present "lapw0 electrostatic MPI rotation index is pre-advanced" "DO jatom[[:space:]]*=[[:space:]]*1,[[:space:]]*nstart\(myid\)-1" "${lapw0_elect}"
check_present "lapw0 GGA MPI rotation index is pre-advanced" "DO jatom[[:space:]]*=[[:space:]]*1,[[:space:]]*nstart\(myid\)-1" "${lapw0_gga}"

g1_factor_count=$(grep -Ec 'ohalf\*(ohalf|g1_deriv_scale)\*gaunt1' "${lapw2}" || true)
if [ "${g1_factor_count}" -ge 4 ]; then
  echo "PASS  lapw2 applies the F3 factor to APW and all LO combinations"
else
  echo "FAIL  lapw2 has only ${g1_factor_count} of 4 corrected F3 factors" >&2
  fail=1
fi

if grep -Eq "WIEN_STRESS_DEBUG" "${wienroot}/SRC_lapw0/finl_elect_str.F" && grep -Eq "WIEN_STRESS_DEBUG" "${lcore}" && grep -Eq "DBG_STR_CORE_RADIAL" "${lcore}" && grep -Eq "WIEN_STRESS_DEBUG" "${wienroot}/SRC_lapw2/symmetrization_stress.f" && grep -Eq "DBG_STR_VAL" "${wienroot}/SRC_lapw2/symmetrization_stress.f" && grep -Eq "WIEN_STRESS_DEBUG" "${wienroot}/SRC_mixer/scfana.f"; then
  echo "INFO  optional WIEN_STRESS_DEBUG diagnostics are installed"
else
  echo "INFO  optional diagnostics are not fully installed"
fi

if [ -f "${wienroot}/SRC_lapw1/coors.f" ] && grep -Eq "LAPW1_FIXGVEC_V1" "${wienroot}/SRC_lapw1/coors.f" && grep -Eq "lapw1\.fixgvec\.read" "${wienroot}/SRC_lapw1/coors.f"; then
  echo "INFO  optional serial fixed-G LAPW1 diagnostic is installed"
else
  echo "INFO  optional fixed-G LAPW1 diagnostic is not installed"
fi

exit "${fail}"
