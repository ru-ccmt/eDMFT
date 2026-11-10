"""K-distributed response, with replicated real-space Krylov coordinates.

Only k-dependent matrices, wavefunctions and frequency caches are partitioned.
Weights retain their GLOBAL normalization. Scalar/core/radial field coordinates
are replicated; their inner product must not be summed over MPI ranks again.
Collectives occur at explicit local/global map boundaries, including dX and dmu.
No mpi4py dependency is needed for an ordinary serial import or calculation.
"""

from __future__ import annotations
import numpy as np


class ResponseParallel:
    def __init__(self, comm=None):
        self.comm = comm
        self.rank = 0 if comm is None else comm.Get_rank()
        self.size = 1 if comm is None else comm.Get_size()
        self.root = self.rank == 0

    @classmethod
    def world(cls):
        from mpi4py import MPI

        return cls(MPI.COMM_WORLD)

    def reduce(self, value, operation="sum"):
        if self.size == 1:
            return value
        from mpi4py import MPI

        op = {"sum": MPI.SUM, "min": MPI.MIN, "max": MPI.MAX}[operation]
        if isinstance(value, np.ndarray):
            source = np.ascontiguousarray(value)
            result = np.empty_like(source)
            self.comm.Allreduce(source, result, op=op)
            return result
        return self.comm.allreduce(value, op=op)

    def sum(self, value):
        return self.reduce(value)

    def root_call(self, operation):
        """Root-only I/O with success/error broadcast before dependent work."""
        if self.size == 1:
            return operation()
        payload = None
        if self.root:
            try:
                payload = (True, operation())
            except Exception as error:
                payload = (False, f"{type(error).__name__}: {error}")
        ok, result = self.comm.bcast(payload, root=0)
        if not ok:
            raise RuntimeError(result)
        return result

    def fingerprints(self, paths):
        """Hash disjoint files on each rank; exchange only small digest maps.

        Stream files to bound temporary memory even for a very large vector.
        Header/inventory discovery stays small and does not read dense matrices.
        """
        import hashlib
        from pathlib import Path

        partial = {}
        for i, (name, path) in enumerate(sorted(paths.items())):
            if i % self.size != self.rank:
                continue
            digest = hashlib.sha256()
            with Path(path).open("rb") as stream:
                for block in iter(lambda: stream.read(4 * 1024**2), b""):
                    digest.update(block)
            partial[str(name)] = digest.hexdigest()
        return {k: v for part in self.gather(partial) for k, v in part.items()}

    def gather(self, value):
        return [value] if self.size == 1 else self.comm.allgather(value)

    def partition(self, costs):
        """Deterministic longest-job-first ownership; global IDs stay intact.

        A nonempty partition is required by compact matrix coordinates. Extra
        ranks are rejected explicitly, never silently left waiting in a solve.
        Preparation scales roughly with the cube of the LAPW matrix dimension;
        use that estimate, not renormalized Brillouin-zone integration weights.
        """
        costs = np.asarray(costs, dtype=float)
        if len(costs) < self.size or np.any(costs <= 0):
            raise ValueError(
                "response MPI ranks must not exceed the number of k points"
            )
        loads = np.zeros(self.size)
        owners = np.empty(len(costs), dtype=int)
        for k in np.argsort(-costs, kind="stable"):
            rank = int(np.argmin(loads))
            owners[k] = rank
            loads[rank] += costs[k]
        return np.flatnonzero(owners == self.rank)

    def check_transpose(
        self,
        operator,
        *,
        input_local=False,
        output_local=False,
        seed=59,
        samples=3,
        tolerance=1e-11,
    ):
        """Distributed dot test; count each replicated field exactly once."""
        rx = np.random.default_rng(seed + (1009 * self.rank if input_local else 0))
        ry = np.random.default_rng(seed + 1 + (1013 * self.rank if output_local else 0))
        errors = []

        def inner(a, b, local):
            value = float(a @ b)
            return self.sum(value) if local else value

        for unused in range(samples):
            x = rx.normal(size=operator.input_size)
            y = ry.normal(size=operator.output_size)
            ax, aty = operator.matvec(x), operator.rmatvec(y)
            lhs = inner(y, ax, output_local)
            rhs = inner(x, aty, input_local)
            scale = max(
                np.sqrt(inner(y, y, output_local) * inner(ax, ax, output_local)),
                np.sqrt(inner(x, x, input_local) * inner(aty, aty, input_local)),
                np.finfo(float).tiny,
            )
            errors.append(abs(lhs - rhs) / scale)
        error = max(errors)
        if error > tolerance:
            raise ValueError(
                f"{operator.name}: distributed transpose failed: {error:g}"
            )
        return error

    def local_to_global(self, operator):
        """Sum local output fields; transpose feeds the same covector to each k pool."""
        from stress_response_adjoint import PairedOperator

        return PairedOperator(
            operator.input_size,
            operator.output_size,
            lambda x: self.sum(operator.matvec(x)),
            operator.rmatvec,
            "sum_k " + operator.name,
        )

    def global_to_local(self, operator):
        """Replicated input, local k output; transpose sums local input covectors."""
        from stress_response_adjoint import PairedOperator

        return PairedOperator(
            operator.input_size,
            operator.output_size,
            operator.matvec,
            lambda y: self.sum(operator.rmatvec(y)),
            operator.name + " replicated input",
        )

    def sum_operator(self, operator):
        """Sum k-dependent contributions with BOTH input/output fields replicated."""
        return self.global_to_local(self.local_to_global(operator))


SERIAL = ResponseParallel()


def parallel_of(model):
    """Serial-compatible for independently constructed diagnostic states."""
    return getattr(model, "parallel", SERIAL)


def limit_blas_threads():
    """Prevent nested BLAS workers from oversubscribing a k-point MPI launch.

    Linux BLAS/OpenMP use the launch environment. Recent macOS Accelerate
    additionally provides a per-calling-thread API; VECLIB_MAXIMUM_THREADS
    alone does not reliably constrain its newer BLAS/LAPACK interface.
    This is an execution policy, not a change to numerical operations.
    """
    import sys

    if sys.platform != "darwin":
        return "launch environment"
    import ctypes

    try:
        library = ctypes.CDLL(
            "/System/Library/Frameworks/Accelerate.framework/Accelerate"
        )
        setter = library.BLASSetThreading
    except (OSError, AttributeError):
        return "launch environment (legacy Accelerate)"
    setter.argtypes = [ctypes.c_uint]
    setter.restype = ctypes.c_int
    if setter(1) != 0:
        raise RuntimeError("Accelerate rejected single-thread BLAS setting")
    return "Accelerate BLAS_THREADING_SINGLE_THREADED"
