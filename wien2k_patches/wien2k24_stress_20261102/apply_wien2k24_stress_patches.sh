#!/usr/bin/env bash
set -eu

usage() {
  cat <<'EOF'
Usage:
  apply_wien2k24_stress_patches.sh [options] WIENROOT

Options:
  --check                 Check applicability only (default).
  --apply                 Apply or upgrade to the complete stress sources.
  --with-debug            Accepted for compatibility; diagnostics are included
                          in the cumulative patch and remain default-off.
  --build                 Build all affected WIEN2k programs.
  --install               Install rebuilt executables (implies --build).
  --install-helpers DIR   Install the two Python diagnostics into DIR.
  -h, --help              Show this help.

The installer recognizes pristine WIEN2k-24, the former required-only and
required+debug patch sets, and the complete target before operator caching.
It never invokes sudo and refuses unknown or partially modified sources.
Site-specific Makefile compiler/library settings are preserved; the patch
only adds val_rela_kinstrmt.o to the LAPW2 object list.
EOF
}

mode=check
do_build=0
do_install=0
helper_dir=

while [ "$#" -gt 0 ]; do
  case "$1" in
    --check) mode=check ;;
    --apply) mode=apply ;;
    --with-debug) : ;;
    --build) do_build=1 ;;
    --install) do_build=1; do_install=1 ;;
    --install-helpers)
      shift
      [ "$#" -gt 0 ] || { echo '--install-helpers requires DIR' >&2; exit 2; }
      helper_dir=$1
      ;;
    -h|--help) usage; exit 0 ;;
    --*) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    *) break ;;
  esac
  shift
done

[ "$#" -eq 1 ] || { usage >&2; exit 2; }
wienroot=${1%/}
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
target_manifest=${script_dir}/TARGET_SOURCES.sha256
previous_manifest=${script_dir}/PRE_OPERATOR_TARGET_SOURCES.sha256

complete_patch=${script_dir}/patches/current/0001-wien2k24-complete-stress.patch
upgrade_required=${script_dir}/patches/current/0002-upgrade-required-to-complete.patch
upgrade_full=${script_dir}/patches/current/0003-upgrade-full-legacy-to-complete.patch
upgrade_cache=${script_dir}/patches/current/0004-upgrade-complete-to-operator-cache.patch

for file in "$target_manifest" "$previous_manifest" "$complete_patch" "$upgrade_required" "$upgrade_full" "$upgrade_cache"; do
  [ -f "$file" ] || { echo "Missing bundle file: $file" >&2; exit 1; }
done

target_matches() {
  (cd "$wienroot" && shasum -a 256 -c "$target_manifest" >/dev/null 2>&1) &&
    grep -Fq 'val_rela_kinstrmt.o' "$wienroot/SRC_lapw2/Makefile"
}

patch_applies() {
  patch --dry-run --batch --forward -p1 -d "$wienroot" < "$1" >/dev/null 2>&1
}

state=incompatible
selected_patch=
if target_matches; then
  state=applied
elif (cd "$wienroot" && shasum -a 256 -c "$previous_manifest" >/dev/null 2>&1) &&
    grep -Fq 'val_rela_kinstrmt.o' "$wienroot/SRC_lapw2/Makefile"; then
  state=complete-before-operator-cache
  selected_patch=$upgrade_cache
elif patch_applies "$complete_patch"; then
  state=pristine
  selected_patch=$complete_patch
elif patch_applies "$upgrade_full"; then
  state=legacy-full
  selected_patch=$upgrade_full
elif patch_applies "$upgrade_required"; then
  state=legacy-required
  selected_patch=$upgrade_required
fi

printf '%-18s %s\n' "$state" "$wienroot"
if [ "$state" = incompatible ]; then
  echo 'Source tree is neither pristine, a recognized legacy/previous-complete state, nor the complete target.' >&2
  echo 'Do not force the patch; compare it with BASE_SOURCES.sha256 and TARGET_SOURCES.sha256.' >&2
  exit 1
fi

if [ "$mode" = apply ] && [ "$state" != applied ]; then
  stamp=$(date +%Y%m%d_%H%M%S)
  backup_dir=${WIEN_STRESS_BACKUP_DIR:-"${wienroot}/stress_patch_backup_${stamp}"}
  mkdir -p "$backup_dir"

  # Back up every existing target file while preserving its relative path.
  while read -r checksum relative; do
    [ -n "$relative" ] || continue
    if [ -f "$wienroot/$relative" ]; then
      mkdir -p "$backup_dir/$(dirname "$relative")"
      cp -p "$wienroot/$relative" "$backup_dir/$relative"
    fi
  done < "$target_manifest"
  cp -p "$wienroot/SRC_lapw2/Makefile" "$backup_dir/SRC_lapw2.Makefile"
  echo "Backed up affected sources to $backup_dir"

  patch --batch --forward -p1 -d "$wienroot" < "$selected_patch"
  if ! target_matches; then
    echo 'Patch completed, but exact target verification failed.' >&2
    echo "Original files are in $backup_dir" >&2
    exit 1
  fi
  state=applied
  echo 'Complete WIEN2k-24 stress patch verified.'
fi

if [ "$do_build" -eq 1 ]; then
  [ "$state" = applied ] || {
    echo 'Use --apply before --build on an unpatched tree.' >&2
    exit 1
  }
  make -C "$wienroot/SRC_lapw0"
  make -C "$wienroot/SRC_lapw1" real
  make -C "$wienroot/SRC_lapw1" complex
  make -C "$wienroot/SRC_lcore"
  make -C "$wienroot/SRC_mixer"
  make -C "$wienroot/SRC_lapw2" real
  make -C "$wienroot/SRC_lapw2" complex
fi

if [ "$do_install" -eq 1 ]; then
  cp "$wienroot/SRC_lapw0/lapw0" "$wienroot/lapw0"
  cp "$wienroot/SRC_lapw1/lapw1" "$wienroot/lapw1"
  cp "$wienroot/SRC_lapw1/lapw1c" "$wienroot/lapw1c"
  cp "$wienroot/SRC_lcore/lcore" "$wienroot/lcore"
  cp "$wienroot/SRC_mixer/mixer" "$wienroot/mixer"
  cp "$wienroot/SRC_lapw2/lapw2" "$wienroot/lapw2"
  cp "$wienroot/SRC_lapw2/lapw2c" "$wienroot/lapw2c"
  echo 'Installed lapw0, lapw1, lapw1c, lcore, mixer, lapw2, and lapw2c.'
fi

if [ -n "$helper_dir" ]; then
  mkdir -p "$helper_dir"
  cp "$script_dir/helpers/rmt_stress_correction.py" "$helper_dir/"
  cp "$script_dir/helpers/radial_response_fd.py" "$helper_dir/"
  chmod 755 "$helper_dir/rmt_stress_correction.py" "$helper_dir/radial_response_fd.py"
  echo "Installed diagnostic helpers in $helper_dir"
fi

