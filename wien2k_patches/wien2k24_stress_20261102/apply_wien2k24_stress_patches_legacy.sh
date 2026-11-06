#!/usr/bin/env bash
set -eu

usage() {
  cat <<'EOF'
Usage:
  apply_wien2k24_stress_patches.sh [options] WIENROOT

Options:
  --check         Check patch applicability only (default).
  --apply         Apply the required source patches.
  --with-debug    Also apply optional diagnostics, including fixed-G LAPW1.
  --build         Build lapw0, lcore, mixer, lapw2/lapw2c, and optional LAPW1 diagnostics.
  --install       Copy rebuilt executables into WIENROOT (implies --build).
  -h, --help      Show this help.

The script never invokes sudo.  For a system installation such as
/opt/apps/... run --check as an ordinary user first, then arrange the
required write permission before using --apply or --install.
EOF
}

mode=check
with_debug=0
do_build=0
do_install=0

while [ "${#}" -gt 0 ]; do
  case "${1}" in
    --check) mode=check ;;
    --apply) mode=apply ;;
    --with-debug) with_debug=1 ;;
    --build) do_build=1 ;;
    --install) do_build=1; do_install=1 ;;
    -h|--help) usage; exit 0 ;;
    --*) echo "Unknown option: ${1}" >&2; usage >&2; exit 2 ;;
    *) break ;;
  esac
  shift
done

if [ "${#}" -ne 1 ]; then
  usage >&2
  exit 2
fi

wienroot=${1%/}
script_dir=$(CDPATH= cd -- "$(dirname -- "${0}")" && pwd)

required_patches="
${script_dir}/patches/required/0001-lapw0-pbe-vsigma.patch
${script_dir}/patches/required/0002-lcore-core-stress.patch
${script_dir}/patches/required/0003-mixer-stress-sign.patch
${script_dir}/patches/required/0004-lapw2-g1-f3-factor.patch
${script_dir}/patches/required/0005-cubic-l8-l10-coefficient.patch
${script_dir}/patches/required/0006-lapw0-spin-and-parallel-index.patch
"
debug_patches="
${script_dir}/patches/optional/0101-lapw0-stress-debug.patch
${script_dir}/patches/optional/0102-lcore-stress-debug.patch
${script_dir}/patches/optional/0103-mixer-stress-debug.patch
${script_dir}/patches/optional/0104-lapw1-fixed-g-debug.patch
${script_dir}/patches/optional/0105-lapw2-stress-debug.patch
"

for path in "${wienroot}/SRC_lapw0/vxclm2.F" \
            "${wienroot}/SRC_lapw0/finl_elect_str.F" \
            "${wienroot}/SRC_lapw0/stress_gga.F" \
            "${wienroot}/SRC_lapw0/xcpot3.F" \
            "${wienroot}/SRC_lcore/core_corr_stress.f" \
            "${wienroot}/SRC_mixer/mixer.F" \
            "${wienroot}/SRC_lapw2/val_kinstrmt.F"
do
  if [ ! -f "${path}" ]; then
    echo "Missing expected WIEN2k source: ${path}" >&2
    exit 1
  fi
done

if [ "${with_debug}" -eq 1 ]; then
  for path in "${wienroot}/SRC_lapw1/coors.f" "${wienroot}/SRC_lapw2/symmetrization_stress.f"
  do
    if [ ! -f "${path}" ]; then
      echo "Missing expected WIEN2k source: ${path}" >&2
      exit 1
    fi
  done
fi

patch_state() {
  patch_file=${1}
  if patch --dry-run --batch --forward -p1 -d "${wienroot}" < "${patch_file}" >/dev/null 2>&1; then
    echo applicable
  elif patch --dry-run --batch --reverse -p1 -d "${wienroot}" < "${patch_file}" >/dev/null 2>&1; then
    echo applied
  elif patch_signature_present "$(basename "${patch_file}")"; then
    # A later optional patch, or an unrelated local edit, can change context
    # lines and prevent GNU patch from recognizing an earlier patch in reverse.
    # In that case use narrow semantic signatures for the known fixes.  The
    # complete set of signatures is checked again by the verifier below.
    echo applied
  else
    echo incompatible
  fi
}

patch_signature_present() {
  patch_name=${1}
  case "${patch_name}" in
    0001-lapw0-pbe-vsigma.patch)
      grep -Eiq "call[[:space:]]+pbe_vsigma_builtin" "${wienroot}/SRC_lapw0/vxclm2.F" &&
      grep -Eiq "subroutine[[:space:]]+pbe_vsigma_builtin" "${wienroot}/SRC_lapw0/vxclm2.F"
      ;;
    0002-lcore-core-stress.patch)
      grep -Eq "counter[[:space:]]*=[[:space:]]*nrp[[:space:]]*$" "${wienroot}/SRC_lcore/core_corr_stress.f" &&
      ! grep -Eq "counter[[:space:]]*=[[:space:]]*nrp[[:space:]]*\+[[:space:]]*141" "${wienroot}/SRC_lcore/core_corr_stress.f" &&
      grep -Eq "sum_tot\(index2,jatom\)[[:space:]]*=.*cor_corr_tot\(index2\)" "${wienroot}/SRC_lcore/core_corr_stress.f"
      ;;
    0003-mixer-stress-sign.patch)
      grep -Eq "final_stress\(ialpha,ibeta\)[[:space:]]*=[[:space:]]*-[[:space:]]*Tstress\(index1\)" "${wienroot}/SRC_mixer/mixer.F"
      ;;
    0004-lapw2-g1-f3-factor.patch)
      { grep -Fq "WIEN2k stress Eq. (G1) F3 factor correction" "${wienroot}/SRC_lapw2/val_kinstrmt.F" &&
        [ "$(grep -Fc 'ohalf*ohalf*gaunt1' "${wienroot}/SRC_lapw2/val_kinstrmt.F")" -ge 4 ]; } ||
      { grep -Eq "g1_deriv_scale[[:space:]]*=[[:space:]]*0\.5d0" "${wienroot}/SRC_lapw2/val_kinstrmt.F" &&
        [ "$(grep -Fc 'ohalf*g1_deriv_scale*gaunt1' "${wienroot}/SRC_lapw2/val_kinstrmt.F")" -ge 4 ]; }
      ;;
    0005-cubic-l8-l10-coefficient.patch)
      grep -Fq "c3 = c_kub(IABS(lm1(1,i)),lm1(2,i)+8 )" "${wienroot}/SRC_lapw0/finl_elect_str.F" &&
      grep -Fq "c3 = c_kub(IABS(lm1(1,i)),lm1(2,i)+8 )" "${wienroot}/SRC_lcore/core_corr_stress.f"
      ;;
    0006-lapw0-spin-and-parallel-index.patch)
      grep -Fq "RHOKSP(J,2)=VKEXC(J,2)-VKVXC(J,2)" "${wienroot}/SRC_lapw0/xcpot3.F" &&
      grep -Fq "DO jatom = 1, nstart(myid)-1" "${wienroot}/SRC_lapw0/finl_elect_str.F" &&
      grep -Fq "DO jatom = 1, nstart(myid)-1" "${wienroot}/SRC_lapw0/stress_gga.F"
      ;;
    0101-lapw0-stress-debug.patch)
      grep -Eq "WIEN_STRESS_DEBUG" "${wienroot}/SRC_lapw0/finl_elect_str.F"
      ;;
    0102-lcore-stress-debug.patch)
      grep -Eq "WIEN_STRESS_DEBUG" "${wienroot}/SRC_lcore/core_corr_stress.f"
      ;;
    0103-mixer-stress-debug.patch)
      grep -Eq "WIEN_STRESS_DEBUG" "${wienroot}/SRC_mixer/scfana.f"
      ;;
    0104-lapw1-fixed-g-debug.patch)
      grep -Eq "LAPW1_FIXGVEC_V1" "${wienroot}/SRC_lapw1/coors.f" &&
      grep -Eq "lapw1\.fixgvec\.read" "${wienroot}/SRC_lapw1/coors.f"
      ;;
    0105-lapw2-stress-debug.patch)
      grep -Eq "DBG_STR_VAL" "${wienroot}/SRC_lapw2/symmetrization_stress.f" &&
      grep -Eq "WIEN_STRESS_DEBUG" "${wienroot}/SRC_lapw2/symmetrization_stress.f"
      ;;
    *)
      return 1
      ;;
  esac
}

bad=0
required_pending=0
for patch_file in ${required_patches}; do
  state=$(patch_state "${patch_file}")
  printf '%-12s %s\n' "${state}" "$(basename "${patch_file}")"
  if [ "${state}" = applicable ]; then
    required_pending=1
  fi
  if [ "${state}" = incompatible ]; then
    bad=1
  fi
done

if [ "${with_debug}" -eq 1 ] && [ "${mode}" = check ]; then
  for patch_file in ${debug_patches}; do
    state=$(patch_state "${patch_file}")
    if [ "${state}" = incompatible ] && [ "${required_pending}" -eq 1 ]; then
      state=deferred
    elif [ "${state}" = incompatible ]; then
      bad=1
    fi
    printf '%-12s %s\n' "${state}" "$(basename "${patch_file}")"
  done
fi

if [ "${bad}" -ne 0 ]; then
  echo "At least one patch is incompatible with this source tree." >&2
  echo "Do not force it; compare that file with the WIEN2k-24 base first." >&2
  exit 1
fi

if [ "${mode}" = apply ]; then
  stamp=$(date +%Y%m%d_%H%M%S)
  backup_dir=${WIEN_STRESS_BACKUP_DIR:-"${wienroot}/stress_patch_backup_${stamp}"}
  mkdir -p "${backup_dir}/SRC_lapw0" "${backup_dir}/SRC_lcore" \
           "${backup_dir}/SRC_mixer" "${backup_dir}/SRC_lapw2"
  if [ "${with_debug}" -eq 1 ]; then
    mkdir -p "${backup_dir}/SRC_lapw1"
  fi
  cp -p "${wienroot}/SRC_lapw0/vxclm2.F" "${backup_dir}/SRC_lapw0/"
  cp -p "${wienroot}/SRC_lapw0/finl_elect_str.F" "${backup_dir}/SRC_lapw0/"
  cp -p "${wienroot}/SRC_lapw0/stress_gga.F" "${backup_dir}/SRC_lapw0/"
  cp -p "${wienroot}/SRC_lapw0/xcpot3.F" "${backup_dir}/SRC_lapw0/"
  cp -p "${wienroot}/SRC_lcore/core_corr_stress.f" "${backup_dir}/SRC_lcore/"
  cp -p "${wienroot}/SRC_mixer/mixer.F" "${backup_dir}/SRC_mixer/"
  cp -p "${wienroot}/SRC_mixer/scfana.f" "${backup_dir}/SRC_mixer/"
  cp -p "${wienroot}/SRC_lapw2/val_kinstrmt.F" "${backup_dir}/SRC_lapw2/"
  if [ "${with_debug}" -eq 1 ]; then
    cp -p "${wienroot}/SRC_lapw1/coors.f" "${backup_dir}/SRC_lapw1/"
    cp -p "${wienroot}/SRC_lapw2/symmetrization_stress.f" "${backup_dir}/SRC_lapw2/"
  fi
  echo "Backed up sources to ${backup_dir}"

  for patch_file in ${required_patches}; do
    state=$(patch_state "${patch_file}")
    if [ "${state}" = applicable ]; then
      patch --batch --forward -p1 -d "${wienroot}" < "${patch_file}"
    else
      echo "Already applied: $(basename "${patch_file}")"
    fi
  done

  if [ "${with_debug}" -eq 1 ]; then
    for patch_file in ${debug_patches}; do
      state=$(patch_state "${patch_file}")
      printf '%-12s %s\n' "${state}" "$(basename "${patch_file}")"
      if [ "${state}" = applicable ]; then
        patch --batch --forward -p1 -d "${wienroot}" < "${patch_file}"
      elif [ "${state}" = applied ]; then
        echo "Already applied: $(basename "${patch_file}")"
      else
        echo "Optional patch is incompatible after required fixes: ${patch_file}" >&2
        exit 1
      fi
    done
  fi
fi

if [ "${do_build}" -eq 1 ]; then
  make -C "${wienroot}/SRC_lapw0"
  make -C "${wienroot}/SRC_lcore"
  make -C "${wienroot}/SRC_mixer"
  make -C "${wienroot}/SRC_lapw2" real
  make -C "${wienroot}/SRC_lapw2" complex
  if [ "${with_debug}" -eq 1 ]; then
    make -C "${wienroot}/SRC_lapw1" real
    make -C "${wienroot}/SRC_lapw1" complex
  fi
fi

if [ "${do_install}" -eq 1 ]; then
  cp "${wienroot}/SRC_lapw0/lapw0" "${wienroot}/lapw0"
  cp "${wienroot}/SRC_lcore/lcore" "${wienroot}/lcore"
  cp "${wienroot}/SRC_mixer/mixer" "${wienroot}/mixer"
  cp "${wienroot}/SRC_lapw2/lapw2" "${wienroot}/lapw2"
  cp "${wienroot}/SRC_lapw2/lapw2c" "${wienroot}/lapw2c"
  if [ "${with_debug}" -eq 1 ]; then
    cp "${wienroot}/SRC_lapw1/lapw1" "${wienroot}/lapw1"
    cp "${wienroot}/SRC_lapw1/lapw1c" "${wienroot}/lapw1c"
    echo "Installed lapw0, lcore, mixer, lapw2, lapw2c, lapw1, and lapw1c into ${wienroot}"
  else
    echo "Installed lapw0, lcore, mixer, lapw2, and lapw2c into ${wienroot}"
  fi
fi

if [ "${mode}" = apply ] || [ "${required_pending}" -eq 0 ]; then
  "${script_dir}/verify_wien2k24_stress_patches.sh" "${wienroot}"
else
  echo "Compatibility check passed; required patches have not been applied."
fi
