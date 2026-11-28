#!/usr/bin/env python3
"""Verification script for IB/methods/sheet-1/qx27 (2007-08 Sheet 1, Q11 / 2004/4/II/15A).

Checks:
1. Inverse FT of f_tilde(k) = (e^k - e^{-k}) 1_{|k|<=1} equals
   f(x) = (2i / (pi (1 + x^2))) (cosh(1) sin(x) - x cos(x) sinh(1)).
2. Strip Laplace solution psi(x, y) matches:
   - Fourier inversion integral of (e^{k(1-y)} - e^{-k(1-y)}) 1_{|k|<=1},
   - Closed form (2i / (pi ((1-y)^2 + x^2))) ((1-y) cosh(1-y) sin(x) - x cos(x) sinh(1-y)),
   - Holomorphic imaginary-part form (2i/pi) Im(sinh(z)/z) with z = (1-y) + i*x.
3. PDE nabla^2 psi = 0 in 0 < y < 1 via 4th-order finite differences.
4. Boundary conditions psi(x, 0) = f(x) and psi(x, 1) = 0, including removable singularity at (0, 1).
5. Parseval energy check in x for several y in [0, 1).
"""

import cmath
import math


def simpson_complex(func, a, b, n=1000):
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    total = func(a) + func(b)
    for j in range(1, n):
        x = a + j * h
        total += (4.0 if j % 2 == 1 else 2.0) * func(x)
    return total * (h / 3.0)


def f_exact(x):
    return (2j / (math.pi * (1.0 + x * x))) * (
        math.cosh(1.0) * math.sin(x) - x * math.cos(x) * math.sinh(1.0)
    )


def psi_closed(x, y):
    alpha = 1.0 - y
    denom = alpha * alpha + x * x
    if denom == 0.0:
        return 0.0j
    return (2j / (math.pi * denom)) * (
        alpha * math.cosh(alpha) * math.sin(x) - x * math.cos(x) * math.sinh(alpha)
    )


def psi_holomorphic(x, y):
    z = complex(1.0 - y, x)
    if abs(z) < 1e-12:
        w = 1.0 + z * z / 6.0 + (z ** 4) / 120.0
    else:
        w = cmath.sinh(z) / z
    return (2j / math.pi) * w.imag


def psi_fourier_integral(x, y):
    alpha = 1.0 - y
    return simpson_complex(
        lambda k: (math.exp(alpha * k) - math.exp(-alpha * k)) * cmath.exp(1j * k * x),
        -1.0, 1.0, n=4000
    ) / (2.0 * math.pi)


def d2_central_4th(func, t, h=1e-3):
    """4th-order central finite difference for second derivative."""
    return (-func(t + 2 * h) + 16 * func(t + h) - 30 * func(t)
            + 16 * func(t - h) - func(t - 2 * h)) / (12.0 * h * h)


def main():
    # 1. Part 1: Inverse FT check
    max_f_err = 0.0
    for x in (-10.0, -3.2, -1.0, -0.1, 0.0, 0.25, 1.0, 2.7, 8.5):
        num = psi_fourier_integral(x, 0.0)
        ex = f_exact(x)
        err = abs(num - ex)
        max_f_err = max(max_f_err, err)
        assert err < 1e-12, f"Part 1 mismatch at x={x}: {num} vs {ex}"

    # 2. Part 2: Closed form vs Fourier integral vs holomorphic form
    max_psi_int_err = 0.0
    max_psi_hol_err = 0.0
    for y in (0.0, 0.15, 0.5, 0.85, 0.999, 1.0):
        for x in (-7.0, -1.5, -0.001, 0.0, 0.001, 0.8, 4.2):
            val_closed = psi_closed(x, y)
            val_int = psi_fourier_integral(x, y)
            val_hol = psi_holomorphic(x, y)
            e1 = abs(val_closed - val_int)
            e2 = abs(val_closed - val_hol)
            max_psi_int_err = max(max_psi_int_err, e1)
            max_psi_hol_err = max(max_psi_hol_err, e2)
            assert e1 < 1e-11, f"psi int mismatch at ({x}, {y}): {val_closed} vs {val_int}"
            assert e2 < 1e-11, f"psi hol mismatch at ({x}, {y}): {val_closed} vs {val_hol}"

    # 3. Boundary conditions check
    max_bc0_err = 0.0
    max_bc1_err = 0.0
    for x in (-12.0, -2.5, -0.3, 0.0, 0.7, 5.0):
        e0 = abs(psi_closed(x, 0.0) - f_exact(x))
        e1 = abs(psi_closed(x, 1.0))
        max_bc0_err = max(max_bc0_err, e0)
        max_bc1_err = max(max_bc1_err, e1)
        assert e0 < 1e-15 and e1 < 1e-15

    # Removable singularity check near (0, 1): psi(x, y) ~ (2i / (3 pi)) * (1-y) * x
    for eps in (1e-2, 1e-3, 1e-4):
        x, y = 2.0 * eps, 1.0 - 3.0 * eps
        val = psi_closed(x, y)
        asymp = (2j / (3.0 * math.pi)) * (1.0 - y) * x
        assert abs(val - asymp) < 10.0 * (eps ** 4)

    # 4. PDE check: psi_xx + psi_yy = 0 via 4th-order finite differences
    max_pde_err = 0.0
    for y in (0.1, 0.3, 0.5, 0.7, 0.9):
        for x in (-4.0, -1.2, 0.0, 0.6, 2.5, 6.0):
            psi_xx = d2_central_4th(lambda u: psi_closed(u, y), x, h=3e-3)
            psi_yy = d2_central_4th(lambda v: psi_closed(x, v), y, h=3e-3)
            res = abs(psi_xx + psi_yy)
            max_pde_err = max(max_pde_err, res)
            assert res < 1e-10, f"PDE residual too large at ({x}, {y}): {res}"

    # 5. Transformed ODE check on [-1, 1] and Parseval check
    # (1/(2pi)) int_{-1}^1 |2 sinh(k alpha)|^2 dk = (2/pi) * (sinh(2 alpha)/(2 alpha) - 1)
    for y in (0.0, 0.4, 0.8):
        alpha = 1.0 - y
        rhs = (2.0 / math.pi) * (math.sinh(2.0 * alpha) / (2.0 * alpha) - 1.0)
        # Integrate |psi(x, y)|^2 on [-X, X]; tail is O(1/X)
        X = 500.0
        lhs = simpson_complex(lambda x: abs(psi_closed(x, y)) ** 2, -X, X, n=40000).real
        assert abs(lhs - rhs) < 2e-3, f"Parseval mismatch at y={y}: {lhs} vs {rhs}"

    print("IB_methods_sheet-1_qx27_check.py: ALL CHECKS PASSED")
    print(f"  Part 1 max inverse FT quadrature error:      {max_f_err:.3e}")
    print(f"  Part 2 max closed-vs-integral error:         {max_psi_int_err:.3e}")
    print(f"  Part 2 max closed-vs-holomorphic error:      {max_psi_hol_err:.3e}")
    print(f"  Boundary condition errors (y=0, y=1):        {max_bc0_err:.3e}, {max_bc1_err:.3e}")
    print(f"  Max 4th-order FD Laplace residual |nabla^2|: {max_pde_err:.3e}")


if __name__ == "__main__":
    main()
