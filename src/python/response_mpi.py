"""Choose response MPI launchers without making each case carry a JSON file.

The caller supplies DmftEnvironment.MPI2, which already honors mpi_prefix.dat2
before mpi_prefix.dat. Native Fortran workers inherit that exact argument list.
Python's mpi4py can be linked to a DIFFERENT MPI implementation (in particular,
conda MPICH versus Homebrew Open MPI on macOS). Reusing the wrong launcher can
start independent singletons or hang, so a plain local launch is translated to
a verified launcher from the active Python environment, retaining its rank count.

Host/scheduler/binding arguments must never disappear in that translation.
Matching MPI implementations retain the complete prefix; scheduler launchers
retain the site's prefix. Cross-implementation translation of a complex prefix
requires an explicit solver.mpi override. An explicit [] means serial. These
choices affect execution only, not the stress formulas or physical inputs.
"""
from functools import lru_cache
from pathlib import Path
import json
import os
import re
import shlex
import shutil
import subprocess
import sys


def argument_list(value, label):
    """Response subprocesses use argv, never shell expansion or evaluation."""
    if not isinstance(value, list) or any(not isinstance(x, str) or not x for x in value):
        raise ValueError(label + ' must be an explicit argument list')
    if any(any(c in x for c in '$;|><\n\r') for x in value):
        raise ValueError(label + ' must not contain shell expansions or operators')
    return list(value)


def _family(text):
    name = text.lower()
    if 'intel' in name and 'mpi' in name:
        return 'intel'
    if 'mvapich' in name:
        return 'mvapich'
    if 'open mpi' in name or 'openrte' in name or 'open-mpi' in name:
        return 'openmpi'
    if 'mpich' in name or 'hydra build details' in name:
        return 'mpich'
    return None


@lru_cache(maxsize=4)
def python_mpi_family(executable):
    """Probe in a short child WITHOUT MPI_Init in the launcher's parent.

    Initializing MPI here would create progress threads and make subsequent
    native subprocess/fork execution unsafe. MPI library identification is
    permitted before initialization; mpi4py's automatic init is disabled.
    """
    code = (
        'import mpi4py; mpi4py.rc.initialize=False; mpi4py.rc.finalize=False; '
        'from mpi4py import MPI; import json; '
        'print(json.dumps(MPI.get_vendor()[0]))'
    )
    env = os.environ | {'PYTHONDONTWRITEBYTECODE': '1'}
    try:
        result = subprocess.run([executable, '-c', code], capture_output=True,
                                text=True, env=env, timeout=10, check=True)
        vendor = json.loads(result.stdout.strip())
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        raise ValueError('Cannot identify mpi4py in the active Python environment; '
                         'install its compatible MPI runtime or set solver.mpi explicitly') from error
    family = _family(vendor)
    if family is None:
        raise ValueError('Unrecognized Python MPI implementation '+repr(vendor)+
                         '; set solver.mpi explicitly')
    return family


@lru_cache(maxsize=32)
def launcher_family(executable):
    try:
        result = subprocess.run([executable, '--version'], capture_output=True,
                                text=True, timeout=5, check=True)
    except (OSError, subprocess.SubprocessError):
        return None
    return _family(result.stdout + '\n' + result.stderr)


def _executable(name):
    expanded = os.path.expanduser(name)
    if '/' in expanded:
        p = Path(expanded)
        return str(p.resolve()) if p.is_file() and os.access(p, os.X_OK) else None
    return shutil.which(expanded)


def plain_local_ranks(command):
    """Recognize only a local mpiexec/mpirun plus one explicit rank count.

    This deliberately excludes hostfiles, placement, environment assignments,
    scheduler wrappers and MPMD syntax. Silently discarding any of them can
    turn a cluster job into hundreds of ranks on its login node.
    """
    if not command or Path(command[0]).name not in ('mpiexec', 'mpirun', 'mpiexec.hydra'):
        return None
    args = command[1:]
    number = None
    if len(args) == 2 and args[0] in ('-n', '-np', '--np', '--ntasks'):
        number = args[1]
    elif len(args) == 1:
        match = re.fullmatch(r'(?:-np|-n|--np=|--ntasks=)([0-9]+)', args[0])
        number = match[1] if match else None
    return int(number) if number and number.isdigit() and int(number) > 0 else None


def _matching_launcher(family):
    # Prefer the runtime installed beside this Python. PATH is a fallback,
    # but every candidate is checked; its filename alone proves nothing.
    directories = [Path(sys.executable).parent, Path(sys.prefix)/'bin']
    candidates = [str(d/name) for d in directories for name in ('mpiexec', 'mpirun')]
    candidates += [shutil.which(name) for name in ('mpiexec', 'mpirun')]
    seen = set()
    for candidate in candidates:
        if not candidate:
            continue
        executable = _executable(candidate)
        if executable and executable not in seen:
            seen.add(executable)
            if launcher_family(executable) == family:
                return executable
    raise ValueError('No launcher matching Python MPI ('+family+
                     ') was found; activate a compatible environment or set solver.mpi')


def _site_allocation():
    """A bare rank count can still mean a multi-node job through environment."""
    return any(os.environ.get(key) for key in (
        'SLURM_JOB_ID', 'SLURM_JOBID', 'PBS_JOBID', 'LSB_JOBID',
        'COBALT_JOBID', 'LOADL_STEP_ID', 'HYDRA_HOST_FILE', 'MPICH_HOSTFILE',
        'I_MPI_HYDRA_HOST_FILE', 'OMPI_MCA_orte_default_hostfile',
        'PRTE_MCA_prte_default_hostfile',
    ))


def automatic_python_launch(native):
    if not native:
        return []
    family = python_mpi_family(sys.executable)
    name = Path(native[0]).name
    # Schedulers launch the site's configured PMI/PMIx environment. Preserve
    # their allocation, task-count, affinity and transport options verbatim.
    if name in ('srun', 'aprun', 'jsrun', 'ibrun'):
        return list(native)
    if name not in ('mpiexec', 'mpirun', 'mpiexec.hydra'):
        raise ValueError('Cannot automatically adapt MPI wrapper '+repr(native[0])+
                         '; set solver.mpi to its Python-compatible argument list')
    executable = _executable(native[0])
    actual = launcher_family(executable) if executable else None
    ranks = plain_local_ranks(native)
    if actual == family:
        launch = list(native)
    elif ranks is not None and not _site_allocation():
        launch = [_matching_launcher(family), '-n', str(ranks)]
    else:
        raise ValueError('Native/Python MPI implementations differ or cannot be verified; '
                         'cannot translate host/binding options or a site allocation safely. '
                         'Set solver.mpi explicitly')
    # Local conda MPICH/OFI on this Mac spent two CPU cores per rank polling
    # while idle. The TCP provider removes that overhead, and passed both MPI
    # response tests and a full SrVO3 run. Do NOT force TCP on Linux clusters,
    # multi-host launches, explicit JSON commands or a user-selected provider.
    if (sys.platform == 'darwin' and family == 'mpich' and ranks is not None
            and not _site_allocation()
            and not os.environ.get('FI_PROVIDER')
            and not any('FI_PROVIDER' in word for word in launch)):
        launch[1:1] = ['-genv', 'FI_PROVIDER', 'tcp']
    return launch


def resolve_launches(prefix, solver, *, analytic=True):
    """Return Python and native argv, honoring independent explicit overrides.

    No prefix means serial by default. JSON keys are optional, and presence
    (including []) wins over inference. The native backend retains its previous
    launch semantics; only the analytic backend needs Python MPI detection.
    """
    # If both commands are supplied, an arbitrary historical shell prefix for
    # ordinary dmft2 need not be interpreted by this separate response path.
    need_prefix = 'mpi' not in solver or 'native_mpi' not in solver
    native_default = argument_list(shlex.split(prefix), 'MPI prefix') if need_prefix else []
    native = argument_list(solver.get('native_mpi', native_default), 'solver.native_mpi')
    if 'mpi' in solver:
        python = argument_list(solver['mpi'], 'solver.mpi')
    else:
        python = automatic_python_launch(native_default) if analytic else native_default
    return python, native
