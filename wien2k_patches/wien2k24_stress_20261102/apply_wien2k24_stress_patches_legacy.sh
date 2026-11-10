#!/usr/bin/env bash
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
echo 'Legacy entry point: upgrading to the current audited sources.' >&2
exec "$script_dir/apply_wien2k24_stress_patches.sh" "$@"
