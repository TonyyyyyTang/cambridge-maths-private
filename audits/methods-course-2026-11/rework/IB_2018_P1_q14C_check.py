#!/usr/bin/env python3
"""Independent verification script for IB/2018/P1/q14C (Methods)."""

import cmath
import math


def simpson_complex(func, a, b, n):
    """Composite Simpson's rule for complex-valued functions (n even)."""
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    total = func(a) + func(b)
    for j in range(1, n):
        x = a + j * h
        total += (4 if j % 2 == 1 else 2) * func(x)
    return total * (h / 3.0)


def simpson_real(func, a, b, n):
    """Composite Simpson's rule for real-valued functions (n even)."""
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    total = func(a) + func(b)
    for j in range(1, n):
        x = a + j * h
        total += (4 if j % 2 == 1 else 2) * func(x)
    return total * (h / 3.0)


def check_convolution_theorem():
    """Verify F{f * g}(k) = F{f}(k) * F{g}(k) on two Gaussian test functions."""
    # f(x) = exp(-x^2), g(x) = exp(-2 (x - 0.4)^2)
    # Exact convolution: (f * g)(x) = sqrt(pi/3) * exp(-(2/3)(x - 0.4)^2)
    def f(x):
        return math.exp(-x * x)

    def g(x):
        return math.exp(-2.0 * (x - 0.4) ** 2)

    def fg_exact(x):
        return math.sqrt(math.pi / 3.0) * math.exp(-(2.0 / 3.0) * (x - 0.4) ** 2)

    # Check convolution integral matches exact formula
    for x in [-1.2, 0.0, 0.7, 1.5]:
        conv_num = simpson_real(lambda y, x=x: f(y) * g(x - y), -8.0, 8.0, 1000)
        assert abs(conv_num - fg_exact(x)) < 1e-12

    # Check F{f*g}(k) = F{f}(k) F{g}(k)
    max_err = 0.0
    for k in [-2.3, -0.5, 0.0, 1.1, 3.4]:
        f_hat = simpson_complex(lambda x, k=k: f(x) * cmath.exp(-1j * k * x), -8.0, 8.0, 1000)
        g_hat = simpson_complex(lambda x, k=k: g(x) * cmath.exp(-1j * k * x), -8.0, 8.0, 1000)
        fg_hat = simpson_complex(
            lambda x, k=k: fg_exact(x) * cmath.exp(-1j * k * x), -10.0, 10.0, 1200
        )
        err = abs(fg_hat - f_hat * g_hat)
        max_err = max(max_err, err)
        assert err < 1e-11, f"Convolution theorem mismatch at k={k}: err={err}"

    print(f"[PASS] Convolution theorem verified (max error = {max_err:.3e})")


def h_hat_exact(k):
    """Exact Fourier transform of h(x) = sin(x)/x^2."""
    if k <= -1.0:
        return 1j * math.pi
    elif k >= 1.0:
        return -1j * math.pi
    else:
        return -1j * math.pi * k


def h_hat_compact(k):
    """Equivalent compact form -(i*pi/2)*(|k+1| - |k-1|)."""
    return -0.5j * math.pi * (abs(k + 1.0) - abs(k - 1.0))


def tail_cos_over_x2(omega, M):
    """Asymptotic tail int_M^infty cos(omega x)/x^2 dx accurate to O(M^-5)."""
    if abs(omega) < 1e-14:
        return 1.0 / M
    return (
        -math.sin(omega * M) / (omega * M ** 2)
        + 2.0 * math.cos(omega * M) / ((omega ** 2) * (M ** 3))
        + 6.0 * math.sin(omega * M) / ((omega ** 3) * (M ** 4))
    )


def check_sinx_over_x2_ft():
    """Verify F{sin(x)/x^2}(k) directly by quadrature + tail and via k-space convolution."""
    # 1. Check equivalence of piecewise and compact formulas
    test_ks = [-3.5, -1.7, -1.0, -0.6, -0.15, 0.0, 0.25, 0.8, 1.0, 1.6, 4.2]
    for k in test_ks:
        assert abs(h_hat_exact(k) - h_hat_compact(k)) < 1e-15

    # 2. Direct quadrature of -2i * int_0^infty sin(kx) sin(x) / x^2 dx
    # Since sin(kx) sin(x) / x^2 is smooth at x=0 (limit = k) and decays like 1/x^2,
    # integrate on [0, M] with Simpson and add analytic integration-by-parts tail on [M, inf).
    M = 200.0
    max_direct_err = 0.0
    for k in test_ks:
        def integrand(x, k=k):
            if x == 0.0:
                return k
            return math.sin(k * x) * math.sin(x) / (x * x)

        # Ensure enough points per oscillation (max freq is |k| + 1)
        n_panels = max(4000, int(40 * M * (abs(k) + 1.0)))
        if n_panels % 2 == 1:
            n_panels += 1
        main_part = simpson_real(integrand, 0.0, M, n_panels)
        tail = 0.5 * (tail_cos_over_x2(k - 1.0, M) - tail_cos_over_x2(k + 1.0, M))
        num_val = -2.0j * (main_part + tail)
        exact_val = h_hat_exact(k)
        err = abs(num_val - exact_val)
        max_direct_err = max(max_direct_err, err)
        assert err < 1e-7, f"Direct FT mismatch at k={k}: num={num_val}, exact={exact_val}"

    # 3. Direct check of dual convolution integral -(i*pi/2) int_{-1}^1 sgn(k - q) dq
    def sgn(v):
        if v > 0:
            return 1.0
        elif v < 0:
            return -1.0
        return 0.0

    max_conv_err = 0.0
    for k in test_ks:
        # Split integral at q = k if k in (-1, 1) to avoid discontinuity inside Simpson panel
        if -1.0 < k < 1.0:
            val_int = simpson_real(lambda q: sgn(k - q), -1.0, k - 1e-12, 100) + simpson_real(
                lambda q: sgn(k - q), k + 1e-12, 1.0, 100
            )
        else:
            val_int = simpson_real(lambda q: sgn(k - q), -1.0 + 1e-12, 1.0 - 1e-12, 200)
        conv_val = -0.5j * math.pi * val_int
        err = abs(conv_val - h_hat_exact(k))
        max_conv_err = max(max_conv_err, err)
        assert err < 1e-10, f"k-space convolution mismatch at k={k}: err={err}"

    # 4. Derivative check: d/dk F{sin(x)/x^2}(k) = -i * F{sin(x)/x}(k) = -i*pi * 1_{|k|<1}
    dk = 1e-5
    for k in [-2.0, -0.5, 0.3, 2.5]:
        deriv = (h_hat_exact(k + dk) - h_hat_exact(k - dk)) / (2.0 * dk)
        expected_deriv = -1j * math.pi if abs(k) < 1.0 else 0.0j
        assert abs(deriv - expected_deriv) < 1e-10

    print(
        f"[PASS] F{{sin(x)/x^2}}(k) verified: direct quadrature max err = {max_direct_err:.3e}, "
        f"k-space convolution max err = {max_conv_err:.3e}"
    )


if __name__ == "__main__":
    check_convolution_theorem()
    check_sinx_over_x2_ft()
    print("ALL CHECKS PASSED FOR IB_2018_P1_q14C")
