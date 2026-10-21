#!/usr/bin/env python3
import ast
import sys
from pathlib import Path
import numpy as np

if len(sys.argv) != 2:
    sys.exit(f"Usage: {Path(sys.argv[0]).name} imp.0/Sbulla.out")

source = Path(sys.argv[1])
output = Path(source.stem + ".inp")

with source.with_name("sig.out").open() as f:
    next(f)  # Skip original s_oo
    Edc = ast.literal_eval(next(f).split("=", 1)[1].strip())

data = np.loadtxt(source, ndmin=2)
s_oo = data[-1, 1::2].copy()

if data.shape[1] != 1 + 2 * len(Edc):
    raise ValueError("Self-energy columns do not match Edc")

data[:, 1::2] -= s_oo

with output.open("w") as f:
    print("# s_oo=", s_oo.tolist(), file=f)
    print("# Edc=", Edc, file=f)
    np.savetxt(f, data, fmt="%.15e")
