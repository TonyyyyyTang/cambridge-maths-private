"""Computational verification for IB/2009/P4/q5B (Lagrange multipliers & max distance)."""

import math
from fractions import Fraction


def f(x, y):
    return x * x + y * y


def g(x, y):
    return x * x + y * y + x * y - 4


def grad_f(x, y):
    return (2 * x, 2 * y)


def grad_g(x, y):
    return (2 * x + y, x + 2 * y)


def det_cond(x, y):
    fx, fy = grad_f(x, y)
    gx, gy = grad_g(x, y)
    return fx * gy - gx * fy


def test_exact_stationary_points():
    # Case y = -x: exact rational points (2, -2) and (-2, 2), lambda = 2
    for x, y in ((Fraction(2), Fraction(-2)), (Fraction(-2), Fraction(2))):
        assert g(x, y) == 0
        assert det_cond(x, y) == 0
        fx, fy = grad_f(x, y)
        gx, gy = grad_g(x, y)
        lam = Fraction(2)
        assert fx == lam * gx and fy == lam * gy
        assert f(x, y) == Fraction(8)

    # Case y = x: points (+-2/sqrt(3), +-2/sqrt(3)), lambda = 2/3
    c = 2.0 / math.sqrt(3.0)
    for sgn in (1.0, -1.0):
        x, y = sgn * c, sgn * c
        assert abs(g(x, y)) < 1e-14
        assert abs(det_cond(x, y)) < 1e-14
        fx, fy = grad_f(x, y)
        gx, gy = grad_g(x, y)
        lam = 2.0 / 3.0
        assert abs(fx - lam * gx) < 1e-14 and abs(fy - lam * gy) < 1e-14
        assert abs(f(x, y) - 8.0 / 3.0) < 1e-14

    print("Exact stationary point checks passed: f_max = 8 (d_max = 2*sqrt(2)), f_min = 8/3 (d_min = 2*sqrt(6)/3)")


def test_dense_ellipse_scan():
    # Exact parametrization of x^2 + y^2 + xy = 4 via rotated axes:
    # (3/2) u^2 + (1/2) v^2 = 4 => u = sqrt(8/3) cos(t), v = sqrt(8) sin(t)
    # x = (u + v) / sqrt(2), y = (u - v) / sqrt(2)
    N = 200000
    d_max = -1.0
    d_min = 1e9
    max_pt = None
    min_pt = None
    max_g_err = 0.0
    det_zero_count = 0

    for k in range(N):
        t = 2.0 * math.pi * k / N
        u = math.sqrt(8.0 / 3.0) * math.cos(t)
        v = math.sqrt(8.0) * math.sin(t)
        x = (u + v) / math.sqrt(2.0)
        y = (u - v) / math.sqrt(2.0)
        max_g_err = max(max_g_err, abs(g(x, y)))
        d = math.hypot(x, y)
        if d > d_max:
            d_max = d
            max_pt = (x, y)
        if d < d_min:
            d_min = d
            min_pt = (x, y)
        # Check determinant sign changes (zeros on the curve)
        t_next = 2.0 * math.pi * (k + 1) / N
        u2 = math.sqrt(8.0 / 3.0) * math.cos(t_next)
        v2 = math.sqrt(8.0) * math.sin(t_next)
        x2 = (u2 + v2) / math.sqrt(2.0)
        y2 = (u2 - v2) / math.sqrt(2.0)
        d1 = det_cond(x, y)
        d2 = det_cond(x2, y2)
        if d1 == 0.0 or (d1 < 0.0 < d2) or (d1 > 0.0 > d2):
            det_zero_count += 1

    assert max_g_err < 1e-13
    expected_d_max = 2.0 * math.sqrt(2.0)
    expected_d_min = 2.0 * math.sqrt(6.0) / 3.0
    assert abs(d_max - expected_d_max) < 1e-12, f"d_max={d_max}, expected={expected_d_max}"
    assert abs(d_min - expected_d_min) < 1e-12, f"d_min={d_min}, expected={expected_d_min}"
    assert det_zero_count == 4, f"Expected 4 stationary points, found {det_zero_count}"

    print(
        f"Dense scan ({N} pts) passed: max_g_err={max_g_err:.2e}, "
        f"d_max={d_max:.12f} at ({max_pt[0]:.6f}, {max_pt[1]:.6f}), "
        f"d_min={d_min:.12f} at ({min_pt[0]:.6f}, {min_pt[1]:.6f}), "
        f"stationary points found={det_zero_count}"
    )


if __name__ == "__main__":
    test_exact_stationary_points()
    test_dense_ellipse_scan()
    print("ALL CHECKS PASSED for IB_2009_P4_q5B")
