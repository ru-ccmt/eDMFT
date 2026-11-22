#!/usr/bin/env python3
"""Compare full LAPW1 solvers on an isolated snapshot of a material's inputs.

Example: python check_material_solvers.py CASE_DIRECTORY --executable ./lapw1c
The original calculation is read only. A small k list and all outputs are
created in a temporary directory, whose location is printed on completion.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import time

import numpy as np


def record(stream):
    head = stream.read(4)
    if not head:
        return None
    count, = struct.unpack("<i", head)
    body = stream.read(count)
    assert struct.unpack("<i", stream.read(4))[0] == count
    return body


def vectors(path, nat):
    result = []
    with path.open("rb") as stream:
        for _ in range(2*nat):
            record(stream)
        while True:
            head = record(stream)
            if head is None:
                return result
            *_, n, m, weight = struct.unpack("<ddd10siid", head)
            labels = record(stream)
            energies, columns = [], []
            for _ in range(m):
                band, energy = struct.unpack("<id", record(stream))
                data = record(stream)
                columns.append(np.frombuffer(data, "<f8" if len(data)==8*n else "<c16"))
                energies.append(energy)
            result.append((np.array(energies), np.column_stack(columns)))


def matrix(path):
    with path.open("rb") as stream:
        magic=stream.read(8)
        assert magic in (b"EDMFTOV1", b"EDMFTHV1")
        assert struct.unpack("<i",stream.read(4))[0] == 0x01020304
        header=np.fromfile(stream,dtype="<i4",count=8)
        n=int(header[3])
        stream.read(24+12*n)
        a=np.empty((n,n),complex)
        for j in range(n):
            a[:j+1,j]=np.fromfile(stream,dtype="<c16",count=j+1)
            a[j,:j]=a[:j,j].conj()
        return a


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    parser.add_argument("--executable",type=Path,required=True)
    parser.add_argument("--case")
    parser.add_argument("--kpoints",type=int,default=2)
    parser.add_argument("--ranks",type=int,default=1)
    args=parser.parse_args()
    if args.kpoints<1 or args.ranks<1:
        parser.error("--kpoints and --ranks must be positive")
    source=args.directory.resolve()
    case=args.case or next(source.glob("*.struct")).stem
    exe=args.executable.resolve()
    iscomplex=exe.name.endswith("c")
    stage=Path(tempfile.mkdtemp(prefix="dftks_material_"))
    base=stage/"input";base.mkdir()
    names=[case+ext for ext in (".struct",".vsp",".vns",".in1c" if iscomplex else ".in1")]
    hashes={}
    for name in names:
        src=source/name
        before=src.read_bytes()
        shutil.copy2(src,base/name)
        after=src.read_bytes()
        assert before==after==(base/name).read_bytes(), f"Input changed during snapshot: {src}"
        hashes[name]=hashlib.sha256(before).hexdigest()
    lines=(source/(case+".klist")).read_text().splitlines()
    (base/(case+".klist")).write_text("\n".join(lines[:args.kpoints])+"\nEND\n")
    nat=int((base/(case+".struct")).read_text().splitlines()[1][27:30])
    definitions={4:(case+".klist","old","formatted"),
                 5:(case+(".in1c" if iscomplex else ".in1"),"old","formatted"),
                 6:(case+".output1","unknown","formatted"),
                 10:(case+".vector","unknown","unformatted"),
                 11:(case+".energy","unknown","formatted"),
                 18:(case+".vsp","old","formatted"),
                 19:(case+".vns","old","formatted"),
                 20:(case+".struct","old","formatted"),
                 21:(case+".scf1","unknown","formatted"),
                 55:(case+".vec","unknown","formatted"),
                 71:(case+".nsh","unknown","unformatted")}
    (base/"lapw1.def").write_text("".join(
        f"{i}, '{f}', '{s}', '{fmt}', 0\n" for i,(f,s,fmt) in definitions.items()))
    env=dict(os.environ,OMP_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1",
             OPENBLAS_NUM_THREADS="1",DFTKS_OVERLAP_DUMP="operators")
    for key in list(env):
        if key.startswith("DFTKS_EIG") or key.startswith("DFTKS_DAVIDSON"):
            env.pop(key)
    cache=stage/"cache"
    results=[];reference=None
    for name,mode in [("gvx","GVX"),("gvd","GVD"),("seed","DAVIDSON"),
                      ("warm","DAVIDSON"),("auto","AUTO")]:
        directory=stage/name;shutil.copytree(base,directory)
        start=time.perf_counter()
        command=[str(exe),"lapw1.def"]
        if args.ranks>1:
            command=["mpirun","-n",str(args.ranks)]+command
        proc=subprocess.run(command,cwd=directory,
             env=dict(env,DFTKS_EIGENSOLVER=mode,DFTKS_EIG_CACHE=str(cache)),
             capture_output=True,text=True,timeout=180)
        elapsed=time.perf_counter()-start
        (directory/"console.log").write_text(proc.stdout+proc.stderr)
        assert proc.returncode==0, (name,proc.stdout,proc.stderr)
        points=vectors(directory/(case+".vector"),nat)
        if reference is None:
            reference=points
        maxerr=maxres=maxmetric=0.
        for k,((e,z),(eref,zref)) in enumerate(zip(points,reference),1):
            assert len(e)==len(eref)
            h=matrix(directory/f"operators.ham.k{k:08d}.bin")
            o=matrix(directory/f"operators.k{k:08d}.bin")
            hz,oz=h@z,o@z
            residual=np.linalg.norm(hz-oz*e,axis=0)/np.maximum(1.,
                 np.linalg.norm(hz,axis=0)+abs(e)*np.linalg.norm(oz,axis=0))
            maxres=max(maxres,float(np.max(residual)))
            maxmetric=max(maxmetric,float(np.max(abs(z.conj().T@o@z-np.eye(len(e))))))
            maxerr=max(maxerr,float(np.max(abs(e-eref))))
        assert len(points)==len(reference)==args.kpoints
        assert maxerr<2e-8 and maxres<2e-9 and maxmetric<2e-8, (name,maxerr,maxres,maxmetric)
        log=(directory/(case+".output1")).read_text()
        results.append(dict(run=name,mode=mode,wall_seconds=elapsed,
             band_counts=[len(e) for e,z in points],eigenvalue_error=maxerr,
             residual=maxres,metric_error=maxmetric,
             solver_log=[line.strip() for line in log.splitlines() if ":DFTKS_" in line]))
        print(name,results[-1],flush=True)
    (stage/"results.json").write_text(json.dumps(dict(source=str(source),
         executable=str(exe),input_sha256=hashes,results=results),indent=2))
    print(f"PASS: {stage / 'results.json'}")


if __name__=="__main__":
    main()
