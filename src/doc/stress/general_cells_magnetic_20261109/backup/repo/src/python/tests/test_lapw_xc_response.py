"""Analytic limits, independent perturbations, and discrete transpose checks."""

import numpy as np
import unittest
import functools
from scipy.special import eval_legendre
from lapw_xc_response import (
    unpolarized_xc,
    local_jacobian,
    NativeSplineDerivative,
    NativeDensityRadialDerivatives,
    MuffinTinXCResponse,
    InterstitialXCResponse,
)
from lapw_hartree_periodic import NativeStarTable, NativeHartreeReanalysis


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


def test_pbe_uniform_density_limit():
    n = np.geomspace(1e-4, 1e4, 100)
    lda = unpolarized_xc(n, 0, 0, 0, functional="LDA", charge_cut=1e-8)
    pbe = unpolarized_xc(n, 0, 0, 0, functional="PBE", charge_cut=1e-8)
    np.testing.assert_allclose(pbe, lda, rtol=4e-16, atol=4e-14)


@cases("LDA", "PBE")
def test_local_jacobian_finite_difference(functional):
    x = np.array(
        [[0.2, 0.3, -0.1, 0.12], [0.001, 0.01, 0.003, -0.001], [1e-5, 0.03, 0.05, 0.01]]
    )

    def f(a):
        return unpolarized_xc(*a.T, functional=functional, charge_cut=1e-9)

    val, j = local_jacobian(f, x)
    d = x * np.array([0.12, -0.2, 0.06, 0.13])
    h = 1e-4
    np.testing.assert_allclose(
        np.sum(j * d, axis=1),
        (f(x + h * d) - f(x - h * d)) / (2 * h),
        rtol=5e-8,
        atol=1e-9,
    )


def test_zero_density_and_functional_guard():
    for f in ["LDA", "PBE"]:
        assert unpolarized_xc(0, 0, 0, 0, functional=f, charge_cut=1e-8) == 0
    with unittest.TestCase().assertRaises(NotImplementedError):
        unpolarized_xc(1, 0, 0, 0, functional="SCAN", charge_cut=1e-8)


def test_spline_cubic_and_transpose():
    x = np.geomspace(0.03, 3.0, 50)
    op = NativeSplineDerivative(x)
    np.testing.assert_allclose(op.action(x**3), 3 * x * x, rtol=5e-13, atol=1e-12)
    rng = np.random.default_rng(4)
    a = rng.normal(size=(len(x), 3))
    b = rng.normal(size=a.shape)
    np.testing.assert_allclose(
        np.sum(b * op.action(a)), np.sum(a * op.pullback(b)), rtol=1e-13
    )


def test_radial_parity_and_origin_transpose():
    r = np.geomspace(0.001, 2.0, 100)
    op = NativeDensityRadialDerivatives(r, [0, 1, 2])
    rng = np.random.default_rng(1)
    a = rng.normal(size=(len(r), 3))
    b = rng.normal(size=a.shape)
    c = rng.normal(size=a.shape)
    d1, d2 = op.action(a)
    np.testing.assert_allclose(
        np.sum(b * d1 + c * d2), np.sum(a * op.pullback(b, c)), rtol=1e-12
    )


def mt(functional, center=None):
    x, w = np.polynomial.legendre.leggauss(24)
    sin = np.sqrt(1 - x * x)
    r = np.geomspace(0.001, 2.0, 101)
    ells = [0, 2]
    y = np.column_stack(
        [
            np.ones(len(x)) / np.sqrt(4 * np.pi),
            np.sqrt(5 / (4 * np.pi)) * (3 * x * x - 1) / 2,
        ]
    )
    dt = np.zeros_like(y)
    dtt = np.zeros_like(y)
    dt[:, 1] = -np.sqrt(5 / (4 * np.pi)) * 3 * x * sin
    dtt[:, 1] = np.sqrt(5 / (4 * np.pi)) * 3 * (1 - 2 * x * x)
    if center is None:
        center = np.column_stack(
            [
                4 * np.pi * r * r * (0.6 * np.exp(-r * r) + 0.2),
                0.01 * r**4 * np.exp(-r * r),
            ]
        )
    return (
        MuffinTinXCResponse(
            r,
            ells,
            center,
            y,
            2 * np.pi * w,
            functional=functional,
            sine=sin,
            cosine=x,
            angular_derivatives=[dt, dtt, y * 0, y * 0, y * 0],
        ),
        center,
    )


@cases("LDA", "PBE")
def test_mt_response_and_transpose(functional):
    op, c = mt(functional)
    d = c * np.array([0.2, -0.3])
    h = 1e-3
    p, _ = mt(functional, c + h * d)
    m, _ = mt(functional, c - h * d)
    np.testing.assert_allclose(
        op.action(d), (p.central - m.central) / (2 * h), rtol=2e-6, atol=2e-8
    )
    assert op.operator.check_transpose(samples=2) < 1e-12


def ir(functional, center=None, rho_min=(2e-8, 1e-9)):
    reps = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    stars = NativeStarTable(
        reps,
        np.array([np.eye(3, dtype=int), -np.eye(3, dtype=int)]),
        np.zeros((2, 3)),
        inversion=True,
    )
    length = np.array([8.0, 9.0, 10.0])
    rean = NativeHartreeReanalysis(
        stars,
        (12, 12, 12),
        np.diag(2 * np.pi / length),
        np.array([[0.2, 0.3, 0.4]]),
        np.array([1.0]),
        np.prod(length),
    )
    if center is None:
        center = np.array([0.3, 0.03, 0.01, -0.02], complex)
    return (
        InterstitialXCResponse(
            rean, center, length, functional=functional, rho_min=rho_min
        ),
        center,
    )


@cases("LDA", "PBE")
def test_ir_response_and_transpose(functional):
    op, c = ir(functional)
    d = np.array([0.01, -0.02, 0.005, 0.002], complex)
    h = 1e-3
    p, _ = ir(functional, c + h * d)
    m, _ = ir(functional, c - h * d)
    np.testing.assert_allclose(
        op.action(d), (p.central - m.central) / (2 * h), rtol=2e-7, atol=1e-10
    )
    assert op.operator.check_transpose(samples=2) < 1e-12


def test_vacuum_average_reverse():
    op, c = ir("LDA", rho_min=(0.16, 0.14))
    assert op.average_count > 0
    rng = np.random.default_rng(2)
    x = rng.normal(size=op.shape)
    y = rng.normal(size=op.shape)
    np.testing.assert_allclose(
        np.sum(op.average(x) * y), np.sum(x * op.average(y, True)), rtol=1e-13
    )
    assert op.operator.check_transpose(samples=2) < 1e-12
