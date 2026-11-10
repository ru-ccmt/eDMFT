"""Fixed-geometry angular and band projection of native radial H/O blocks.

The caller supplies the actual native angular factors and ordered coefficient
rows. This adapter does not infer an atom's phases, rotate local orbitals, or
generate a strain. The same map is used for H and O and has an exact transpose.
"""

import numpy as np

from stress_response_adjoint import HermitianCoordinates, PairedOperator


class NativeBandProjection:
    """Project ONE radial channel after the native upper-triangle copy.

    Input is a real, potentially nonsymmetric n x n radial block in native
    coefficient order. ``angular`` includes phases, angular sums, volume and
    multiplicities. ``coefficients`` are the corresponding n rows of central
    KS eigenvectors, with NO density, occupation, spin or k-point weighting.

    For F(A)=triu(A)+triu(A,1)^dagger (real diagonal), the forward map is
        B(R) = C^dagger F(angular * R) C.
    In particular, do not replace F by (A+A^dagger)/2: native H can retain
    a finite-stencil Wronskian asymmetry before copying its upper triangle.

    A sector may cover only a subset of the LAPW coefficient rows; those
    rows must be supplied in the SAME increasing order as the native matrix.
    Contributions from different sectors are added after band projection.
    """

    def __init__(self, angular, coefficients):
        self.angular = np.array(angular, dtype=complex, copy=True)
        self.coefficients = np.array(coefficients, dtype=complex, copy=True)
        a, c = self.angular, self.coefficients
        if a.ndim != 2 or a.shape[0] != a.shape[1] or not a.shape[0]:
            raise ValueError("angular block must be nonempty and square")
        if c.ndim != 2 or c.shape[0] != a.shape[0] or not c.shape[1]:
            raise ValueError("coefficient rows do not match the radial block")
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(c)):
            raise ValueError("nonfinite angular factors or coefficients")
        self.n = len(a)
        self.coordinates = HermitianCoordinates((c.shape[1], c.shape[1]))
        self.operator = PairedOperator(
            self.n**2,
            self.coordinates.size,
            lambda x: self.coordinates.pack(self.project(x.reshape(a.shape))),
            lambda y: self.pullback(self.coordinates.unpack(y)).ravel(),
            "native triangle, angular factors and central band projection",
        )

    def project(self, radial):
        radial = np.asarray(radial)
        if (
            np.iscomplexobj(radial)
            or radial.shape != self.angular.shape
            or not np.all(np.isfinite(radial))
        ):
            raise ValueError("invalid real radial block")
        upper = np.triu(radial * self.angular, 1)
        native = upper + upper.conj().T
        np.fill_diagonal(native, np.diag(radial * self.angular).real)
        c = self.coefficients
        return c.conj().T @ native @ c

    def pullback(self, band_covector):
        # Validate Hermiticity using the SAME isometric coordinate convention.
        self.coordinates.pack(band_covector)
        c = self.coefficients
        full = c @ band_covector @ c.conj().T
        # Both mirrored entries contribute to the ONE independent upper entry.
        upper = 2 * np.triu(full, 1)
        np.fill_diagonal(upper, np.diag(full).real)
        return (self.angular.conj() * upper).real
