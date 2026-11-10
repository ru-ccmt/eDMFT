"""Independent Coulomb benchmarks, normalization, and exact discrete adjoints."""

import numpy as np
import unittest
import functools
from lapw_core_response import (
    CorePotentialExtension,
    MatchedCoreState,
    AtomCoreResponse,
)


def cases(*values):
    def decorate(function):
        function.cases = values
        return function

    return decorate


def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    for name, function in sorted(globals().items()):
        if name.startswith("test_") and callable(function):
            if hasattr(function, "cases"):
                for value in function.cases:
                    suite.addTest(
                        unittest.FunctionTestCase(
                            functools.update_wrapper(
                                functools.partial(function, value), function
                            ),
                            description=f"{name} ({value})",
                        )
                    )
            else:
                suite.addTest(unittest.FunctionTestCase(function))
    return suite


def state(rela=False, z=5.0, occupation=2.0):
    r = np.geomspace(1e-5, 20.0, 1201)
    gamma = np.sqrt(1 - (z / 137.0359895) ** 2) if rela else 1.0
    op = MatchedCoreState(
        r,
        -z / r,
        -z * z / (1 + gamma),
        z=z,
        kappa=-1,
        match=int(np.argmin(abs(r - 0.3))),
        end=len(r) - 1,
        occupation=occupation,
        relativistic=rela,
        output_points=len(r),
    )
    return op, gamma


@cases(False, True)
def test_exact_coulomb_density(rela):
    op, gamma = state(rela)
    exact = op.r ** (2 * gamma) * np.exp(-10 * op.r)
    exact *= 2 / (op.weights @ exact)
    np.testing.assert_allclose(op.density, exact, atol=2e-8, rtol=2e-8)


@cases(False, True)
def test_core_normalization_and_adjoint(rela):
    op, _ = state(rela)
    response = op.action(0.1 * np.sin(op.r))
    assert abs(op.weights @ response) < 1e-13
    assert op.operator.check_transpose(samples=2) < 1e-12


@cases(False, True)
def test_uniform_potential_shift_localized_core(rela):
    op, _ = state(rela)
    # The asymptotic boundary is at a finite zero potential; exponential
    # localization makes this a valid gauge check to the stated tolerance.
    assert abs(op.energy_gradient @ np.ones(len(op.r)) - 1) < 1e-8
    assert np.linalg.norm(op.action(np.ones(len(op.r)))) < 1e-8


def test_nonrelativistic_coulomb_scaling():
    op, _ = state(False)
    # Finite R0 omits the derivative of the Coulomb origin parameter. In the
    # NREL r^1 origin this leaves a vanishing-grid error, tested at R0=1e-5.
    np.testing.assert_allclose(
        op.action(-1 / op.r), op.density * (3 / 5 - 2 * op.r), atol=2e-7, rtol=2e-6
    )
    assert abs(op.energy_gradient @ (-1 / op.r) + 5) < 2e-8


def test_extension_coulomb_and_transpose():
    r = np.geomspace(1e-4, 2.0, 201)
    ext = CorePotentialExtension(r, np.full(len(r), -10.0), shift_ry=0)
    assert np.max(abs(ext.central[: len(r)] + 5 / r)) < 1e-10
    assert ext.operator.check_transpose(samples=2) < 1e-13
    direction = 0.2 * r * np.cos(r)
    h = 1e-4
    p = CorePotentialExtension(r, -10 + h * direction)
    m = CorePotentialExtension(r, -10 - h * direction)
    np.testing.assert_allclose(
        ext.action(direction), (p.central - m.central) / (2 * h), atol=5e-8, rtol=1e-6
    )


def test_core_rejects_bad_indices_and_supercritical_origin():
    r = np.geomspace(1e-5, 20.0, 301)
    args = dict(
        z=1,
        kappa=-1,
        match=200,
        end=300,
        occupation=2,
        relativistic=True,
        output_points=301,
    )
    for change in [
        dict(end=301),
        dict(match=200.5),
        dict(output_points=0),
        dict(z=140),
        dict(occupation=np.nan),
        dict(relativistic="RELA"),
    ]:
        with unittest.TestCase().assertRaises(ValueError):
            MatchedCoreState(r, -1 / r, -0.5, **(args | change))
