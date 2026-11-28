#!/usr/bin/env python3
"""Verification script for IB/methods/sheet-1/qx24 (2007-08 Sheet 1, Q6).

Checks Fourier transforms (convention: f_tilde(k) = int_{-inf}^{inf} f(x) e^{-ikx} dx)
of 1, e^{iax}, sin(ax), cos(ax) supported on |x| < c.
"""

import cmath
import math


def simpson_complex(func, a, b, n=6000):
    """Composite Simpson's rule for complex-valued functions on [a, b]."""
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    total = func(a) + func(b)
    for j in range(1, n):
        x = a + j * h
        total += (4.0 if j % 2 == 1 else 2.0) * func(x)
    return total * (h / 3.0)


def sinc_unnorm(u):
    """sin(u)/u with removable singularity at u = 0."""
    if abs(u) < 1e-12:
        return 1.0 - (u * u) / 6.0 + (u ** 4) / 120.0
    return math.sin(u) / u


def ft1_exact(k, c):
    """FT of 1 on |x| < c."""
    return 2.0 * c * sinc_unnorm(k * c)


def ft2_exact(k, a, c):
    """FT of e^{iax} on |x| < c."""
    return 2.0 * c * sinc_unnorm((k - a) * c)


def ft3_exact_split(k, a, c):
    """FT of sin(ax) on |x| < c (split sinc form)."""
    return -1j * c * (sinc_unnorm((k - a) * c) - sinc_unnorm((k + a) * c))


def ft3_exact_combined(k, a, c):
    """FT of sin(ax) on |x| < c (combined fraction form for k != +/- a)."""
    return (2j * (k * math.cos(k * c) * math.sin(a * c)
                  - a * math.sin(k * c) * math.cos(a * c))) / (k * k - a * a)


def ft4_exact_split(k, a, c):
    """FT of cos(ax) on |x| < c (split sinc form)."""
    return c * (sinc_unnorm((k - a) * c) + sinc_unnorm((k + a) * c))


def ft4_exact_combined(k, a, c):
    """FT of cos(ax) on |x| < c (combined fraction form for k != +/- a)."""
    return (2.0 * (k * math.sin(k * c) * math.cos(a * c)
                   - a * math.cos(k * c) * math.sin(a * c))) / (k * k - a * a)


def main():
    max_err_quad = 0.0
    max_err_comb = 0.0
    max_err_euler = 0.0

    param_sets = [
        (0.0, 1.5),
        (1.0, 2.0),
        (2.3, 1.7),
        (-1.8, 0.9),
        (3.5, 2.4),
    ]

    for a, c in param_sets:
        k_vals = [0.0, a, -a, a + 1e-7, -a - 1e-7, 0.5, -1.2, 4.1, -3.7]
        for k in k_vals:
            # Numerical FTs via composite Simpson
            num1 = simpson_complex(lambda x: cmath.exp(-1j * k * x), -c, c)
            num2 = simpson_complex(lambda x: cmath.exp(1j * (a - k) * x), -c, c)
            num3 = simpson_complex(lambda x: math.sin(a * x) * cmath.exp(-1j * k * x), -c, c)
            num4 = simpson_complex(lambda x: math.cos(a * x) * cmath.exp(-1j * k * x), -c, c)

            ex1 = ft1_exact(k, c)
            ex2 = ft2_exact(k, a, c)
            ex3 = ft3_exact_split(k, a, c)
            ex4 = ft4_exact_split(k, a, c)

            for num, ex in ((num1, ex1), (num2, ex2), (num3, ex3), (num4, ex4)):
                err = abs(num - ex)
                max_err_quad = max(max_err_quad, err)
                assert err < 1e-11, f"Quad mismatch at a={a}, c={c}, k={k}: {num} vs {ex}"

            # Euler linearity check: FT(cos) + i FT(sin) = FT(e^{iax})
            err_e = abs((ex4 + 1j * ex3) - ex2)
            max_err_euler = max(max_err_euler, err_e)
            assert err_e < 1e-14

            # Combined fraction check away from k = +/- a
            if abs(k * k - a * a) > 1e-3:
                err_c3 = abs(ft3_exact_combined(k, a, c) - ex3)
                err_c4 = abs(ft4_exact_combined(k, a, c) - ex4)
                max_err_comb = max(max_err_comb, err_c3, err_c4)
                assert err_c3 < 1e-13 and err_c4 < 1e-13

        # Special-value checks
        assert abs(ft1_exact(0.0, c) - 2.0 * c) < 1e-14
        assert abs(ft2_exact(a, a, c) - 2.0 * c) < 1e-14
        if a != 0.0:
            expected_f3_a = -1j * (c - math.sin(2.0 * a * c) / (2.0 * a))
            expected_f4_a = c + math.sin(2.0 * a * c) / (2.0 * a)
            assert abs(ft3_exact_split(a, a, c) - expected_f3_a) < 1e-14
            assert abs(ft3_exact_split(-a, a, c) + expected_f3_a) < 1e-14
            assert abs(ft4_exact_split(a, a, c) - expected_f4_a) < 1e-14
            assert abs(ft4_exact_split(-a, a, c) - expected_f4_a) < 1e-14

        # Symmetry checks
        for k in (0.7, 2.9):
            assert abs(ft1_exact(-k, c) - ft1_exact(k, c)) < 1e-14
            assert abs(ft3_exact_split(-k, a, c) + ft3_exact_split(k, a, c)) < 1e-14
            assert abs(ft3_exact_split(k, a, c).real) < 1e-14
            assert abs(ft4_exact_split(-k, a, c) - ft4_exact_split(k, a, c)) < 1e-14

    # Parseval check: (1/(2pi)) int_{-inf}^{inf} |ft(k)|^2 dk = int_{-c}^{c} |f(x)|^2 dx
    # Use tail bound O(1/K) with K = 800 to verify within ~5e-3 (and exact energy values).
    a, c = 1.7, 1.3
    K = 600.0
    N = 120000
    p1_num = simpson_complex(lambda k: abs(ft1_exact(k, c)) ** 2, -K, K, N).real / (2.0 * math.pi)
    p2_num = simpson_complex(lambda k: abs(ft2_exact(k, a, c)) ** 2, -K, K, N).real / (2.0 * math.pi)
    p3_num = simpson_complex(lambda k: abs(ft3_exact_split(k, a, c)) ** 2, -K, K, N).real / (2.0 * math.pi)
    p4_num = simpson_complex(lambda k: abs(ft4_exact_split(k, a, c)) ** 2, -K, K, N).real / (2.0 * math.pi)

    p1_ex = 2.0 * c
    p2_ex = 2.0 * c
    p3_ex = c - math.sin(2.0 * a * c) / (2.0 * a)
    p4_ex = c + math.sin(2.0 * a * c) / (2.0 * a)

    parseval_errs = [abs(p1_num - p1_ex), abs(p2_num - p2_ex),
                     abs(p3_num - p3_ex), abs(p4_num - p4_ex)]
    assert max(parseval_errs) < 2e-3, f"Parseval errs too large: {parseval_errs}"

    print("IB_methods_sheet-1_qx24_check.py: ALL CHECKS PASSED")
    print(f"  Max quadrature error (all 4 transforms, all params/special k): {max_err_quad:.3e}")
    print(f"  Max split-vs-combined error (sin, cos):                        {max_err_comb:.3e}")
    print(f"  Max Euler identity error (FT(cos)+i*FT(sin)-FT(e^{{iax}})):     {max_err_euler:.3e}")
    print(f"  Parseval truncation errors (K={K:.0f}): {[round(e, 6) for e in parseval_errs]}")


if __name__ == "__main__":
    main()
