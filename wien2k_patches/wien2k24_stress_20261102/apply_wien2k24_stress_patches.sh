#!/usr/bin/env bash
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "${PYTHON:-python3}" "$script_dir/apply_wien2k24_stress.py" "$@"
