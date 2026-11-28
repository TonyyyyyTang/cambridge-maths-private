#!/usr/bin/env python3
"""Computational verification for IB/2014/P2/q5D (Fourier series of x(1-x) and term-by-term differentiation)."""

import math
from fractions import Fraction

# 8-point Gauss-Legendre nodes and weights on [-1, 1]
GL_NODES = [
    -0.96028985649753623168,
    -0.79666647741362673959,
    -0.52553240991632898582,
    -0.18343464249564980494,
    0.18343464249564980494,
    0.52553240991632898582,
    0.79666647741362673959,
    0.96028985649753623168,
]
GL_WEIGHTS = [
    0.10122853629037625915,
    0.22238103445337447054,
    0.31370664587788728734,
    0.36268378337836198297,
    0.36268378337836198297,
    0.31370664587788728734,
    0.22238103445337447054,
    0.10122853629037625915,
]


def quad(func, a, b, panels=40):
    """Composite 8-point Gauss-Legendre quadrature on [a, b]."""
    h = (b - a) / panels
    total = 0.0
    for j in range(panels):
        x0 = a + j * h
        mid = x0 + 0.5 * h
        half = 0.5 * h
        for t, w in zip(GL_NODES, GL_WEIGHTS):
            total += w * func(mid + half * t)
    return 0.5 * h * total


def f_periodic(x):
    u = x - math.floor(x)
    return u * (1.0 - u)


def f_prime_periodic(x):
    u = x - math.floor(x)
    return 1.0 - 2.0 * u


def main():
    # 1. Check Fourier coefficients of f(x) = x(1-x) on [0, 1)
    a0_num = 2.0 * quad(lambda x: x * (1.0 - x), 0.0, 1.0)
    assert abs(a0_num - 1.0 / 3.0) < 1e-14, f"a0 mismatch: {a0_num}"

    max_err_a = 0.0
    max_err_b = 0.0
    for n in range(1, 15):
        an_num = 2.0 * quad(lambda x: x * (1.0 - x) * math.cos(2.0 * math.pi * n * x), 0.0, 1.0, panels=60)
        bn_num = 2.0 * quad(lambda x: x * (1.0 - x) * math.sin(2.0 * math.pi * n * x), 0.0, 1.0, panels=60)
        an_exact = -1.0 / (math.pi ** 2 * n ** 2)
        max_err_a = max(max_err_a, abs(an_num - an_exact))
        max_err_b = max(max_err_b, abs(bn_num))
    assert max_err_a < 1e-13, f"a_n mismatch: {max_err_a}"
    assert max_err_b < 1e-13, f"b_n mismatch: {max_err_b}"

    # 2. Check pointwise convergence of Fourier series of f on [-1, 2]
    N = 4000
    test_xs = [-0.7, -0.25, 0.0, 0.1, 0.25, 0.5, 0.83, 1.0, 1.4]
    max_fs_err = 0.0
    for x in test_xs:
        s = 1.0 / 6.0 - sum(math.cos(2.0 * math.pi * n * x) / (math.pi ** 2 * n ** 2) for n in range(1, N + 1))
        err = abs(s - f_periodic(x))
        max_fs_err = max(max_fs_err, err)
    assert max_fs_err < 3e-5, f"FS convergence error too large: {max_fs_err}"

    # 3. Exact Parseval identity check using Fractions:
    # int_0^1 x^2(1-x)^2 dx = 1/3 - 2/4 + 1/5 = 1/30
    # (a0/2)^2 + (1/2) sum_{n>=1} a_n^2 = (1/6)^2 + (1/(2 pi^4)) * (pi^4 / 90) = 1/36 + 1/180 = 1/30
    lhs_parseval = Fraction(1, 3) - Fraction(2, 4) + Fraction(1, 5)
    rhs_parseval = Fraction(1, 6) ** 2 + Fraction(1, 2) * Fraction(1, 90)
    assert lhs_parseval == Fraction(1, 30) and rhs_parseval == Fraction(1, 30)

    # 4. Check Fourier coefficients of f'(x) = 1 - 2x match term-by-term derivative
    alpha0_num = 2.0 * quad(lambda x: 1.0 - 2.0 * x, 0.0, 1.0)
    assert abs(alpha0_num) < 1e-14

    max_err_alpha = 0.0
    max_err_beta = 0.0
    for n in range(1, 15):
        alphan_num = 2.0 * quad(lambda x: (1.0 - 2.0 * x) * math.cos(2.0 * math.pi * n * x), 0.0, 1.0, panels=60)
        betan_num = 2.0 * quad(lambda x: (1.0 - 2.0 * x) * math.sin(2.0 * math.pi * n * x), 0.0, 1.0, panels=60)
        # From IBP: alpha_n = 2*pi*n*b_n = 0, beta_n = -2*pi*n*a_n = 2/(pi*n)
        betan_exact = 2.0 / (math.pi * n)
        max_err_alpha = max(max_err_alpha, abs(alphan_num))
        max_err_beta = max(max_err_beta, abs(betan_num - betan_exact))
    assert max_err_alpha < 1e-13, f"alpha_n mismatch: {max_err_alpha}"
    assert max_err_beta < 1e-13, f"beta_n mismatch: {max_err_beta}"

    # 5. Check G(0) = 0 = 0.5 * (f'(0+) + f'(1-)) and interior convergence of G(x)
    G0 = (2.0 / math.pi) * sum(math.sin(2.0 * math.pi * n * 0.0) / n for n in range(1, 5000))
    midpoint = 0.5 * ((1.0 - 2.0 * 0.0) + (1.0 - 2.0 * 1.0))
    assert G0 == 0.0 and midpoint == 0.0

    # Check interior convergence of G(x) to 1 - 2x using Cesaro / large-N sum
    N_g = 20000
    max_g_err = 0.0
    for x in [0.15, 0.25, 0.5, 0.75, 0.85]:
        gx = (2.0 / math.pi) * sum(math.sin(2.0 * math.pi * n * x) / n for n in range(1, N_g + 1))
        max_g_err = max(max_g_err, abs(gx - f_prime_periodic(x)))
    assert max_g_err < 5e-4, f"G(x) convergence error too large: {max_g_err}"

    # 6. Contrast check: discontinuous periodic extension h(x)=x on [0,1) fails alpha_0 = 0
    # h'(x) = 1 has alpha_0 = 2, whereas derivative of constant term of h gives 0.
    h_alpha0 = 2.0 * quad(lambda x: 1.0, 0.0, 1.0)
    assert abs(h_alpha0 - 2.0) < 1e-14

    print("IB_2014_P2_q5D verification PASSED:")
    print(f"  a0 = {a0_num:.16f} (exact 1/3)")
    print(f"  max |a_n - (-1/(pi^2 n^2))| (n=1..14) = {max_err_a:.3e}")
    print(f"  max |b_n| (n=1..14) = {max_err_b:.3e}")
    print(f"  max |S_{N}(x) - f_per(x)| = {max_fs_err:.3e}")
    print(f"  Parseval exact: {lhs_parseval} == {rhs_parseval}")
    print(f"  max |alpha_n| = {max_err_alpha:.3e}, max |beta_n - 2/(pi n)| = {max_err_beta:.3e}")
    print(f"  G(0) = {G0}, midpoint 0.5*(f'(0+)+f'(1-)) = {midpoint}")
    print(f"  max |G_{N_g}(x) - f'_per(x)| on interior points = {max_g_err:.3e}")


if __name__ == "__main__":
    main()
