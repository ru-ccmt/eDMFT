#!/bin/sh
# Diagnostic only: link unmodified native MIXER quadrature/mask/normalization.
set -eu
native_wien_root=${WIENROOT:-/Users/haule/wien2k_24}
native_fortran=${FC:-/opt/homebrew/bin/gfortran}
cd "$(dirname "$0")"
"$native_fortran" -shared -fPIC -O1 -ffree-form -ffree-line-length-none \
  reference_modules.f90 \
  "$native_wien_root/SRC_mixer/charge.f" \
  "$native_wien_root/SRC_mixer/setn.f" \
  "$native_wien_root/SRC_mixer/normal.f" \
  reference_wrapper.f90 -o libmixer_normalization.dylib
