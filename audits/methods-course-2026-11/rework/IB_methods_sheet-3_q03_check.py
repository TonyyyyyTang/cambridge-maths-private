#!/usr/bin/env python3
"""Computational verification for IB/methods/sheet-3/q03 (2025-26 Sheet 3, Q3)."""

import math


def gauss_legendre_8(func, a, b, panels=40):
    """Composite 8-point Gauss-Legendre quadrature on [a, b]."""
    if abs(b - a) < 1e-15:
        return 0.0
    nodes = [
        -0.9602898564975363,
        -0.7966664774136267,
        -0.5255324099163290,
        -0.1834346424956498,
        0.1834346424956498,
        0.5255324099163290,
        0.7966664774136267,
        0.9602898564975363,
    ]
    weights = [
        0.1012285362903763,
        0.2223810344533745,
        0.3137066458778873,
        0.3626837833783620,
        0.3626837833783620,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763,
    ]
    h = (b - a) / panels
    total = 0.0
    for i in range(panels):
        x0 = a + i * h
        mid = x0 + 0.5 * h
        half = 0.5 * h
        s = 0.0
        for nd, wt in zip(nodes, weights):
            s += wt * func(mid + half * nd)
        total += half * s
    return total


def G(x, xi):
    """Green's function for y'' - 2y' + y = f(x) with y(0) = y(1) = 0."""
    if x <= xi:
        return -x * (1.0 - xi) * math.exp(x - xi)
    return -xi * (1.0 - x) * math.exp(x - xi)


def G_tilde(x, xi):
    """Self-adjoint Green's function for (e^{-2x} y')' + e^{-2x} y = delta(x - xi)."""
    return math.exp(2.0 * xi) * G(x, xi)


def solve_via_G(x, f):
    """Compute y(x) = int_0^1 G(x, xi) f(xi) dxi, splitting at xi = x."""
    left = gauss_legendre_8(lambda xi: G(x, xi) * f(xi), 0.0, x, panels=30)
    right = gauss_legendre_8(lambda xi: G(x, xi) * f(xi), x, 1.0, panels=30)
    return left + right


def solve_explicit_formula(x, f):
    """Compute y(x) = -e^x [(1-x) int_0^x xi e^{-xi} f(xi) dxi + x int_x^1 (1-xi) e^{-xi} f(xi) dxi]."""
    i1 = gauss_legendre_8(lambda xi: xi * math.exp(-xi) * f(xi), 0.0, x, panels=30)
    i2 = gauss_legendre_8(lambda xi: (1.0 - xi) * math.exp(-xi) * f(xi), x, 1.0, panels=30)
    return -math.exp(x) * ((1.0 - x) * i1 + x * i2)


def check_greens_function():
    """Verify BCs, homogeneous ODE, continuity, jump condition, and symmetry properties."""
    h = 1e-4
    for xi in [0.15, 0.37, 0.72, 0.91]:
        # Boundary conditions
        assert abs(G(0.0, xi)) < 1e-15
        assert abs(G(1.0, xi)) < 1e-15

        # Homogeneous ODE away from xi
        for x in [0.5 * xi, 0.5 * (xi + 1.0)]:
            gp = (G(x + h, xi) - G(x - h, xi)) / (2 * h)
            gpp = (G(x + h, xi) - 2 * G(x, xi) + G(x - h, xi)) / (h ** 2)
            assert abs(gpp - 2 * gp + G(x, xi)) < 1e-6

        # Continuity at x = xi
        g_left = G(xi - 1e-10, xi)
        g_right = G(xi + 1e-10, xi)
        assert abs(g_right - g_left) < 1e-9

        # Jump condition [dG/dx]_{xi^-}^{xi^+} = 1
        eps = 1e-6
        gx_right = (-3 * G(xi, xi) + 4 * G(xi + eps, xi) - G(xi + 2 * eps, xi)) / (2 * eps)
        gx_left = (3 * G(xi, xi) - 4 * G(xi - eps, xi) + G(xi - 2 * eps, xi)) / (2 * eps)
        jump = gx_right - gx_left
        assert abs(jump - 1.0) < 1e-9, f"Jump mismatch at xi={xi}: {jump}"

        # Self-adjoint G_tilde jump condition = e^{2 xi}
        gt_right = (-3 * G_tilde(xi, xi) + 4 * G_tilde(xi + eps, xi) - G_tilde(xi + 2 * eps, xi)) / (2 * eps)
        gt_left = (3 * G_tilde(xi, xi) - 4 * G_tilde(xi - eps, xi) + G_tilde(xi - 2 * eps, xi)) / (2 * eps)
        assert abs((gt_right - gt_left) - math.exp(2.0 * xi)) < 1e-8

    # Check G is NOT symmetric (in fact G(x,xi) = e^{2(x-xi)} G(xi,x)), while G_tilde IS symmetric
    assert abs(G(0.2, 0.7) - G(0.7, 0.2)) > 0.05
    for x in [0.1, 0.35, 0.8]:
        for xi in [0.2, 0.65, 0.9]:
            assert abs(G(x, xi) - math.exp(2.0 * (x - xi)) * G(xi, x)) < 1e-15
            assert abs(G_tilde(x, xi) - G_tilde(xi, x)) < 1e-15


def check_exact_solutions():
    """Verify Green's function solution against analytic solutions for multiple f(x)."""
    cases = [
        (
            "f(x) = 1",
            lambda x: 1.0,
            lambda x: 1.0 - math.exp(x) + (1.0 - math.exp(-1.0)) * x * math.exp(x),
        ),
        (
            "f(x) = e^x sin(pi x)",
            lambda x: math.exp(x) * math.sin(math.pi * x),
            lambda x: -math.exp(x) * math.sin(math.pi * x) / (math.pi ** 2),
        ),
        (
            "f(x) = e^x (6x - 2)",
            lambda x: math.exp(x) * (6.0 * x - 2.0),
            lambda x: x * x * (x - 1.0) * math.exp(x),
        ),
        (
            "f(x) = e^x (12x^2 - 6x)",
            lambda x: math.exp(x) * (12.0 * x * x - 6.0 * x),
            lambda x: (x ** 4 - x ** 3) * math.exp(x),
        ),
    ]

    max_err = 0.0
    for name, f, y_exact in cases:
        for x in [0.0, 0.13, 0.42, 0.79, 1.0]:
            y1 = solve_via_G(x, f)
            y2 = solve_explicit_formula(x, f)
            ye = y_exact(x)
            err = max(abs(y1 - ye), abs(y2 - ye))
            max_err = max(max_err, err)
            assert err < 1e-13, f"{name} at x={x}: err={err}"

            # Also verify self-adjoint integral int_0^1 G_tilde(x,xi) e^{-2xi} f(xi) dxi
            y_sa = (
                gauss_legendre_8(lambda xi: G_tilde(x, xi) * math.exp(-2.0 * xi) * f(xi), 0.0, x, panels=30)
                + gauss_legendre_8(lambda xi: G_tilde(x, xi) * math.exp(-2.0 * xi) * f(xi), x, 1.0, panels=30)
            )
            assert abs(y_sa - ye) < 1e-13

    # Verify non-elementary forcing f(x) = cos(2x) + x^2 directly via ODE FD
    f4 = lambda x: math.cos(2.0 * x) + x * x
    assert abs(solve_explicit_formula(0.0, f4)) < 1e-15
    assert abs(solve_explicit_formula(1.0, f4)) < 1e-15
    h = 1e-4
    max_fd_res = 0.0
    for x in [0.2, 0.55, 0.85]:
        ym = solve_explicit_formula(x - h, f4)
        y0 = solve_explicit_formula(x, f4)
        yp = solve_explicit_formula(x + h, f4)
        d1 = (yp - ym) / (2 * h)
        d2 = (yp - 2 * y0 + ym) / (h ** 2)
        res = abs(d2 - 2 * d1 + y0 - f4(x))
        max_fd_res = max(max_fd_res, res)
        assert res < 1e-6

    return max_err, max_fd_res


def main():
    check_greens_function()
    max_err, max_fd_res = check_exact_solutions()
    print("IB_methods_sheet-3_q03_check.py: ALL CHECKS PASSED")
    print(f"  max exact vs Green-integral error:  {max_err:.3e}")
    print(f"  max FD ODE residual (cos(2x)+x^2):  {max_fd_res:.3e}")


if __name__ == "__main__":
    main()
