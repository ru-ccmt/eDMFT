"""Compose the native moving-basis charge response and its exact transpose.

The native adapters are explicit required inputs. This module cannot infer a
Hartree/XC, core, or basis-density response from band matrices alone, and never
replaces an absent adapter with a frozen approximation. The existing finite-
difference WIEN driver remains the production reference until these physical
adapters and six strain sources have passed native comparisons.
"""

from dataclasses import dataclass
import numpy as np

from stress_response_adjoint import (
    PairedOperator,
    HermitianCoordinates,
    ScreenedResponse,
)


class GeneralizedCoordinates:
    """Isometric real coordinates for the joint (dH,dO,dphi) perturbation.

    Hermitian matrices use independent coordinates with sqrt(2) on their
    off-diagonal entries. The unnormalized complex orbital overlaps use
    separate real/imaginary arrays, with no sqrt(2), spin or k weighting.
    """

    def __init__(self, band_shape, raw_shape):
        self.band = HermitianCoordinates(band_shape)
        self.raw_shape = tuple(raw_shape)
        if not self.raw_shape or min(self.raw_shape) < 1:
            raise ValueError("invalid raw-projector shape")
        self.raw_size = int(np.prod(self.raw_shape))
        self.size = 2 * self.band.size + 2 * self.raw_size

    def pack(self, h, o, raw):
        raw = np.asarray(raw)
        if raw.shape != self.raw_shape or not np.all(np.isfinite(raw)):
            raise ValueError("invalid raw-projector coordinates")
        return np.concatenate(
            (self.band.pack(h), self.band.pack(o), raw.real.ravel(), raw.imag.ravel())
        )

    def unpack(self, vector):
        v = np.asarray(vector)
        if np.iscomplexobj(v) or v.shape != (self.size,) or not np.all(np.isfinite(v)):
            raise ValueError("invalid generalized coordinates")
        n, r = self.band.size, self.raw_size
        return (
            self.band.unpack(v[:n]),
            self.band.unpack(v[n : 2 * n]),
            (v[2 * n : 2 * n + r] + 1j * v[2 * n + r :]).reshape(self.raw_shape),
        )


def generalized_electronic_operator(model, fixed_number=True):
    coordinates = GeneralizedCoordinates(
        (*model.eps.shape, model.eps.shape[1]), model.raw.shape
    )
    band = coordinates.band
    operator = PairedOperator(
        coordinates.size,
        band.size,
        lambda v: band.pack(
            model.generalized_response(
                *coordinates.unpack(v), fixed_number=fixed_number
            )[0]
        ),
        lambda y: coordinates.pack(
            *model.generalized_adjoint(band.unpack(y), fixed_number=fixed_number)
        ),
        "moving-basis frozen-Sigma electronic response",
    )
    return operator, coordinates


@dataclass
class ChargeResponseChain:
    """One prepared charge-map Jacobian with its reverse action.

    K: charge -> potential (Hartree/XC);
    B: potential -> joint band (dH,dO,dphi);
    E: joint band -> contravariant density coefficients;
    R: density coefficients -> spatial valence charge at fixed basis;
    D: potential -> direct density change of the basis at fixed rho0;
    C: potential -> core density;
    N: derivative of the implemented output-charge normalization.

       J = N [R E B + D + C] K.

    D is essential: E converts the changing overlap metric back to density
    coefficients. Only R E B + D represents a physical density change.
    All operators use REAL Euclidean pairings in their stated coordinates.
    Required quadrature/star/multiplicity factors belong in the native maps.
    E already includes spin and k weights; R must not apply them again.
    """

    potential: PairedOperator
    basis: PairedOperator
    electronic: PairedOperator
    density: PairedOperator
    basis_density: PairedOperator
    core: PairedOperator
    normalization: PairedOperator

    def __post_init__(self):
        if self.normalization.input_size != self.normalization.output_size:
            raise ValueError("charge normalization must preserve coordinate space")
        if self.potential.input_size != self.normalization.output_size:
            raise ValueError("input/output charge coordinates differ")
        self.valence = self.density @ self.electronic @ self.basis
        self.charge_from_potential = self.normalization @ (
            self.valence + self.basis_density + self.core
        )
        self.jacobian = self.charge_from_potential @ self.potential
        self.screened = ScreenedResponse(self.jacobian)

    def explicit_source(
        self, band_source, basis_density_source, core_source, potential_source
    ):
        """Assemble ONE physical strain source; geometry is supplied by caller.

        band_source includes explicit geometry changes of H/O/phi at fixed
        potential. basis_density_source includes coordinate/volume transport
        and explicit basis changes, core_source explicit core geometry, and
        potential_source the explicit geometric Hartree/XC/ionic potential.
        No shape change or strain is manufactured by this composition layer.
        """
        direct = self.density.matvec(self.electronic.matvec(band_source))
        for name, value in (("basis", basis_density_source), ("core", core_source)):
            value = np.asarray(value)
            if (
                np.iscomplexobj(value)
                or value.shape != direct.shape
                or not np.all(np.isfinite(value))
            ):
                raise ValueError(f"invalid {name} strain-density source")
            direct = direct + value
        return self.normalization.matvec(direct) + self.charge_from_potential.matvec(
            potential_source
        )

    def check_transposes(self, **controls):
        """Check leaves as well as the complete chain; not a physics closure test."""
        names = (
            "potential",
            "basis",
            "electronic",
            "density",
            "basis_density",
            "core",
            "normalization",
            "jacobian",
        )
        return {name: getattr(self, name).check_transpose(**controls) for name in names}
