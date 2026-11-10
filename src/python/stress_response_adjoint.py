#!/usr/bin/env python3
"""Paired real linear operators and screened six-direction stress contractions.

This layer never infers a transpose by probing one density coefficient at a
time. Each physical adapter must supply BOTH its action and exact discrete
transpose. Native LAPW spatial-potential/radial/core adapters must be validated
before this solver can replace the finite-difference WIEN response backend.

Coordinate convention: real Euclidean vectors. For Hermitian band matrices,
HermitianCoordinates is an isometry for Re Tr(A^dagger B). Spatial density
coefficients require their own explicit coordinates and all integration
weights in their forward/pullback maps; no continuum metric is assumed here.
"""
from __future__ import annotations

from dataclasses import dataclass
import time
from typing import Callable

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres

RY_BOHR3_TO_GPA = 14710.507641078
STRAIN_COMPONENTS = ("xx", "yy", "zz", "yz", "xz", "xy")


def _real_vector(x, size, name):
    x = np.asarray(x)
    if np.iscomplexobj(x) or x.shape != (size,) or not np.all(np.isfinite(x)):
        raise ValueError(f"{name} must be a finite real vector of length {size}")
    return np.asarray(x, float)


@dataclass(frozen=True)
class PairedOperator:
    input_size: int
    output_size: int
    action: Callable
    transpose: Callable
    name: str = "operator"

    def __post_init__(self):
        if self.input_size < 1 or self.output_size < 1:
            raise ValueError("operator dimensions must be positive")

    def matvec(self, x):
        return _real_vector(
            self.action(_real_vector(x, self.input_size, self.name)),
            self.output_size,
            self.name + " result",
        )

    def rmatvec(self, y):
        return _real_vector(
            self.transpose(_real_vector(y, self.output_size, self.name)),
            self.input_size,
            self.name + " transpose result",
        )

    def __matmul__(self, other):
        if other.output_size != self.input_size:
            raise ValueError("operator composition dimensions disagree")
        return PairedOperator(
            other.input_size,
            self.output_size,
            lambda x: self.matvec(other.matvec(x)),
            lambda y: other.rmatvec(self.rmatvec(y)),
            self.name + " @ " + other.name,
        )

    def __add__(self, other):
        if (self.input_size, self.output_size) != (other.input_size, other.output_size):
            raise ValueError("operator addition dimensions disagree")
        return PairedOperator(
            self.input_size,
            self.output_size,
            lambda x: self.matvec(x) + other.matvec(x),
            lambda y: self.rmatvec(y) + other.rmatvec(y),
            self.name + " + " + other.name,
        )

    def check_transpose(self, seed=59, samples=3, tolerance=1e-11):
        rng = np.random.default_rng(seed)
        errors = []
        for _ in range(samples):
            x, y = rng.normal(size=self.input_size), rng.normal(size=self.output_size)
            ax, aty = self.matvec(x), self.rmatvec(y)
            error = abs(np.dot(y, ax) - np.dot(x, aty))
            scale = max(
                np.linalg.norm(y) * np.linalg.norm(ax),
                np.linalg.norm(x) * np.linalg.norm(aty),
                np.finfo(float).tiny,
            )
            errors.append(float(error / scale))
        if max(errors) > tolerance:
            raise ValueError(f"{self.name}: transpose test failed: {max(errors):.3g}")
        return max(errors)


class HermitianCoordinates:
    """Independent n^2 real coordinates per Hermitian n x n matrix.

    Diagonal entries are unscaled; real/imaginary upper-triangle entries
    carry sqrt(2). This avoids redundant variables and shear-like factors
    of two in transpose tests and Krylov inner products.
    """

    def __init__(self, shape):
        self.shape = tuple(shape)
        if len(shape) < 2 or shape[-2] != shape[-1] or min(shape) < 1:
            raise ValueError("invalid Hermitian matrix shape")
        self.n = shape[-1]
        self.count = int(np.prod(shape[:-2]))
        self.size = self.count * self.n * self.n
        self.upper = np.triu_indices(self.n, 1)

    def pack(self, matrix):
        matrix = np.asarray(matrix)
        if matrix.shape != self.shape or not np.all(np.isfinite(matrix)):
            raise ValueError("invalid Hermitian coordinate input")
        if np.max(abs(matrix - matrix.swapaxes(-1, -2).conj())) > 1e-10 * max(
            1.0, np.max(abs(matrix))
        ):
            raise ValueError("coordinate input is not Hermitian")
        a = matrix.reshape((self.count, self.n, self.n))
        upper = a[:, self.upper[0], self.upper[1]]
        return np.concatenate(
            (
                np.diagonal(a, axis1=1, axis2=2).real,
                np.sqrt(2) * upper.real,
                np.sqrt(2) * upper.imag,
            ),
            axis=1,
        ).ravel()

    def unpack(self, vector):
        v = _real_vector(vector, self.size, "Hermitian coordinates").reshape(
            self.count, -1
        )
        a = np.zeros((self.count, self.n, self.n), complex)
        diag = np.arange(self.n)
        a[:, diag, diag] = v[:, : self.n]
        nupper = len(self.upper[0])
        upper = (
            v[:, self.n : self.n + nupper] + 1j * v[:, self.n + nupper :]
        ) / np.sqrt(2)
        a[:, self.upper[0], self.upper[1]] = upper
        a[:, self.upper[1], self.upper[0]] = upper.conj()
        return a.reshape(self.shape)


def electronic_operator(model, *, fixed_number=True, embedding=True):
    coordinates = HermitianCoordinates((*model.eps.shape, model.eps.shape[1]))
    return (
        PairedOperator(
            coordinates.size,
            coordinates.size,
            lambda x: coordinates.pack(
                model.response(coordinates.unpack(x), fixed_number, embedding)[0]
            ),
            lambda y: coordinates.pack(
                model.adjoint(coordinates.unpack(y), fixed_number, embedding)
            ),
            "frozen-Sigma electronic response",
        ),
        coordinates,
    )


class ResponseConvergenceError(RuntimeError):
    def __init__(self, report):
        super().__init__(
            f"{report['direction']} screened response did not converge: residual {report['residual_norm']:.6g}"
        )
        self.report = report


class ScreenedResponse:
    """Solve (I-J)x=b or its transpose without assuming J is symmetric.

    GMRES is appropriate for the non-self-adjoint normalized-projector map.
    Always re-evaluate the final TRUE residual. Failure never yields a stress
    result. No nonlinear charge loop or impurity solve occurs in this layer.
    """

    def __init__(self, jacobian):
        if jacobian.input_size != jacobian.output_size:
            raise ValueError("charge Jacobian must be square")
        self.jacobian = jacobian

    def solve(
        self,
        source,
        *,
        adjoint=False,
        rtol=1e-9,
        atol=1e-11,
        restart=30,
        max_cycles=20,
        initial=None,
    ):
        n = self.jacobian.input_size
        source = _real_vector(source, n, "response source")
        if (
            rtol < 0
            or atol < 0
            or not np.isfinite(rtol + atol)
            or restart < 1
            or max_cycles < 1
        ):
            raise ValueError("invalid screened solver controls")
        if initial is not None:
            initial = _real_vector(initial, n, "initial response")
        count, norms = 0, []
        started = time.monotonic()
        action = self.jacobian.rmatvec if adjoint else self.jacobian.matvec

        def apply(x):
            nonlocal count
            count += 1
            return x - action(x)

        operator = LinearOperator((n, n), matvec=apply, dtype=float)
        answer, info = gmres(
            operator,
            source,
            x0=initial,
            rtol=rtol,
            atol=atol,
            restart=min(restart, n),
            maxiter=max_cycles,
            callback=lambda norm: norms.append(float(norm)),
            callback_type="pr_norm",
        )
        residual = source - apply(answer)
        norm = float(np.linalg.norm(residual))
        threshold = max(atol, rtol * float(np.linalg.norm(source)))
        report = dict(
            converged=bool(info == 0 and norm <= threshold),
            gmres_info=int(info),
            direction="adjoint" if adjoint else "forward",
            residual_norm=norm,
            threshold=threshold,
            actions=count,
            inner_residual_history=norms,
            seconds=time.monotonic() - started,
        )
        if not report["converged"]:
            raise ResponseConvergenceError(report)
        return answer, report

    def contractions(self, gradient, sources, explicit=None, **controls):
        """One adjoint for any number of sources of ONE scalar observable.

        sources[:,nu]=b_nu, gradient=c, result_nu=a_nu+lambda^T b_nu.
        Distinct observable gradients require distinct solves (or a block
        solve); six unrelated gradients cannot share one adjoint.
        """
        sources = np.asarray(sources)
        if (
            np.iscomplexobj(sources)
            or sources.ndim != 2
            or sources.shape[0] != self.jacobian.input_size
            or not np.all(np.isfinite(sources))
        ):
            raise ValueError("sources must have shape (charge coordinates, directions)")
        count = sources.shape[1]
        explicit = (
            np.zeros(count)
            if explicit is None
            else _real_vector(explicit, count, "explicit derivatives")
        )
        lam, report = self.solve(gradient, adjoint=True, **controls)
        return explicit + lam @ sources, lam, report


def engineering_strain_bases():
    """dF/deta for eta=(exx,eyy,ezz,2eyz,2exz,2exy), at F=I."""
    basis = np.zeros((6, 3, 3))
    for i in range(3):
        basis[i, i, i] = 1
    for nu, (i, j) in enumerate([(1, 2), (0, 2), (0, 1)], 3):
        basis[nu, i, j] = basis[nu, j, i] = 0.5
    return basis


def tensile_tensor(derivatives, volume):
    """Energy derivatives in Ry -> tensile-positive GPa; engineering shears."""
    d = _real_vector(derivatives, 6, "six energy derivatives")
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError("volume must be positive")
    xx, yy, zz, yz, xz, xy = d * RY_BOHR3_TO_GPA / volume
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
