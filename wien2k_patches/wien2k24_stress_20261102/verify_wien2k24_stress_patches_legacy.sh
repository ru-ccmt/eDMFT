#!/usr/bin/env bash
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "$script_dir/verify_wien2k24_stress_patches.sh" "$@"
