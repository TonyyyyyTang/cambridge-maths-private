#!/usr/bin/env python3
"""Independent verification for IB/2023/P3/q14A (Fourier transforms & d'Alembert's solution)."""

import cmath
import math


def simpson_c(func, a, b, n=1200):
    """Composite Simpson's rule for complex-valued functions."""
    assert n % 2 == 0
    h = (b - a) / n
    s = func(a) + func(b)
    for k in range(1, n):
        x = a + k * h
        s += (4.0 if k % 2 == 1 else 2.0) * func(x)
    return s * h / 3.0


def ft(func, k, L=12.0, n=1000):
    """Forward Fourier transform: int f(x) e^{-ikx} dx."""
    return simpson_c(lambda x: func(x) * cmath.exp(-1j * k * x), -L, L, n=n)


def ift(func_k, x, K=18.0, n=1200):
    """Inverse Fourier transform: (1/(2pi)) int F(k) e^{ikx} dk."""
    return simpson_c(lambda k: func_k(k) * cmath.exp(1j * k * x), -K, K, n=n) / (2.0 * math.pi)


def main():
    # 1. Verify part (a): Convolution theorem h(x) = int f(y) g(x - y) dy
    def f1(x):
        return math.exp(-x * x) * math.cos(0.7 * x)

    def g1(x):
        return math.exp(-1.5 * (x - 0.4) ** 2)

    # Precompute FT on a grid? Direct calls are fast enough for a handful of x values.
    max_conv_err = 0.0
    for x in [-1.5, -0.3, 0.0, 0.8, 2.1]:
        conv_direct = simpson_c(lambda y: f1(y) * g1(x - y), -10.0, 10.0, n=800).real
        # Analytic FTs for speed & high precision:
        # f1(x) = 0.5 e^{-x^2}(e^{0.7ix} + e^{-0.7ix})
        # FT[e^{-a x^2}](k) = sqrt(pi/a) e^{-k^2/(4a)}
        def f1_hat(k):
            return 0.5 * math.sqrt(math.pi) * (
                math.exp(-((k - 0.7) ** 2) / 4.0) + math.exp(-((k + 0.7) ** 2) / 4.0)
            )

        def g1_hat(k):
            return math.sqrt(math.pi / 1.5) * math.exp(-(k * k) / 6.0) * cmath.exp(-1j * k * 0.4)

        h_ift = ift(lambda k: f1_hat(k) * g1_hat(k), x, K=16.0, n=1000)
        max_conv_err = max(max_conv_err, abs(h_ift - conv_direct))
    # Also verify our analytic f1_hat, g1_hat match numerical ft()
    for k in [-2.3, 0.0, 1.7]:
        assert abs(ft(f1, k) - 0.5 * math.sqrt(math.pi) * (
            math.exp(-((k - 0.7) ** 2) / 4.0) + math.exp(-((k + 0.7) ** 2) / 4.0)
        )) < 1e-11
        assert abs(ft(g1, k) - math.sqrt(math.pi / 1.5) * math.exp(-(k * k) / 6.0) * cmath.exp(-1j * k * 0.4)) < 1e-11
    assert max_conv_err < 1e-11, f"Convolution theorem error: {max_conv_err}"
    print(f"[PASS] Part (a) convolution theorem max err = {max_conv_err:.3e}")

    # 2. Verify part (b): p'(x) = g(x) => p_hat(k) = g_hat(k) / (i k)
    def p_test(x):
        return (1.0 + 0.5 * x) * math.exp(-x * x)

    def g_test(x):
        return (0.5 - 2.0 * x - x * x) * math.exp(-x * x)

    max_deriv_err = 0.0
    for k in [-3.1, -0.7, 1.2, 2.8]:
        p_hat = ft(p_test, k)
        g_hat = ft(g_test, k)
        max_deriv_err = max(max_deriv_err, abs(p_hat - g_hat / (1j * k)))
    assert max_deriv_err < 1e-11, f"Part (b) error: {max_deriv_err}"
    print(f"[PASS] Part (b) p_hat(k) = g_hat(k)/(ik) max err = {max_deriv_err:.3e}")

    # 3. Verify part (c): distributional inverse FTs of cos(ka) and sin(ka)
    # Pairing: int F^{-1}[phi](x) psi(x) dx == (1/(2pi)) int phi(k) psi_hat(-k) dk
    a_val = 1.7
    def psi(x):
        return math.exp(-(x - 0.3) ** 2)

    def psi_hat(k):
        return math.sqrt(math.pi) * math.exp(-(k * k) / 4.0) * cmath.exp(-1j * k * 0.3)

    lhs_cos = simpson_c(lambda k: math.cos(k * a_val) * psi_hat(-k), -16.0, 16.0, n=1200) / (2.0 * math.pi)
    rhs_cos = 0.5 * (psi(-a_val) + psi(a_val))  # from (1/2)[delta(x+a) + delta(x-a)]
    lhs_sin = simpson_c(lambda k: math.sin(k * a_val) * psi_hat(-k), -16.0, 16.0, n=1200) / (2.0 * math.pi)
    rhs_sin = (psi(-a_val) - psi(a_val)) / (2j)  # from (1/(2i))[delta(x+a) - delta(x-a)]
    assert abs(lhs_cos - rhs_cos) < 1e-11, f"cos delta err: {abs(lhs_cos - rhs_cos)}"
    assert abs(lhs_sin - rhs_sin) < 1e-11, f"sin delta err: {abs(lhs_sin - rhs_sin)}"
    print(f"[PASS] Part (c) inverse FT of cos(ka), sin(ka): errs = {abs(lhs_cos-rhs_cos):.2e}, {abs(lhs_sin-rhs_sin):.2e}")

    # 4. Verify part (d): Fourier solution == d'Alembert's solution & satisfies PDE + ICs
    def f_ic(x):
        return math.exp(-x * x)

    def g_ic(x):
        return math.exp(-((x - 0.5) ** 2))

    def f_ic_hat(k):
        return math.sqrt(math.pi) * math.exp(-(k * k) / 4.0)

    def g_ic_hat(k):
        return math.sqrt(math.pi) * math.exp(-(k * k) / 4.0) * cmath.exp(-1j * k * 0.5)

    def u_dalembert(x, t):
        return 0.5 * (f_ic(x + t) + f_ic(x - t)) + 0.5 * simpson_c(g_ic, x - t, x + t, n=1200).real

    def u_fourier(x, t):
        def integrand(k):
            sinc_t = t if abs(k) < 1e-12 else math.sin(k * t) / k
            return f_ic_hat(k) * math.cos(k * t) + g_ic_hat(k) * sinc_t
        return ift(integrand, x, K=18.0, n=1600).real

    max_dal_err = 0.0
    for x in [-1.2, 0.0, 0.9, 2.3]:
        for t in [0.4, 1.5, 2.7]:
            max_dal_err = max(max_dal_err, abs(u_dalembert(x, t) - u_fourier(x, t)))
    assert max_dal_err < 1e-11, f"Fourier vs d'Alembert err: {max_dal_err}"

    # Verify PDE u_tt = u_xx and ICs on u_dalembert
    h = 1e-3
    max_pde_err = 0.0
    for x in [-1.0, 0.3, 1.8]:
        for t in [0.5, 1.9]:
            utt = (-u_dalembert(x, t + 2 * h) + 16 * u_dalembert(x, t + h) - 30 * u_dalembert(x, t) + 16 * u_dalembert(x, t - h) - u_dalembert(x, t - 2 * h)) / (12 * h * h)
            uxx = (-u_dalembert(x + 2 * h, t) + 16 * u_dalembert(x + h, t) - 30 * u_dalembert(x, t) + 16 * u_dalembert(x - h, t) - u_dalembert(x - 2 * h, t)) / (12 * h * h)
            max_pde_err = max(max_pde_err, abs(utt - uxx))
    max_ic_err = 0.0
    for x in [-1.5, 0.0, 1.1]:
        u0 = u_dalembert(x, 0.0)
        ut0 = (-u_dalembert(x, 2 * h) + 8 * u_dalembert(x, h) - 8 * u_dalembert(x, -h) + u_dalembert(x, -2 * h)) / (12 * h)
        max_ic_err = max(max_ic_err, abs(u0 - f_ic(x)), abs(ut0 - g_ic(x)))
    assert max_pde_err < 1e-6 and max_ic_err < 1e-10
    print(f"[PASS] Part (d) Fourier==d'Alembert ({max_dal_err:.2e}), PDE ({max_pde_err:.2e}), ICs ({max_ic_err:.2e})")

    # 5. Verify part (e): f(x)=0, g(x)=x for |x|<=1 else 0
    def g_e(xi):
        return xi if abs(xi) <= 1.0 else 0.0

    def u_e_exact(x, t):
        """Exact piecewise solution via intersection [x-t, x+t] cap [-1, 1]."""
        lo = max(x - t, -1.0)
        hi = min(x + t, 1.0)
        if lo >= hi:
            return 0.0
        return 0.25 * (hi * hi - lo * lo)

    def u_e_5region(x, t):
        assert t > 1.0
        if x <= -1.0 - t:
            return 0.0
        elif x < 1.0 - t:
            return 0.25 * ((x + t) ** 2 - 1.0)
        elif x <= t - 1.0:
            return 0.0
        elif x < t + 1.0:
            return 0.25 * (1.0 - (x - t) ** 2)
        else:
            return 0.0

    # Check exact agreement across all 5 regions for t > 1, and verify u == 0 in |x| <= t - 1 and != 0 just outside
    for t in [1.2, 2.0, 3.5, 5.0]:
        for i in range(-100, 101):
            x = (t + 2.0) * i / 100.0
            v1 = u_e_exact(x, t)
            v2 = u_e_5region(x, t)
            assert abs(v1 - v2) < 1e-15
            if abs(x) <= t - 1.0:
                assert abs(v1) < 1e-15
            elif -1.0 - t + 1e-9 < x < 1.0 - t - 1e-9:
                assert v1 < 0.0
            elif t - 1.0 + 1e-9 < x < t + 1.0 - 1e-9:
                assert v1 > 0.0

    # Independent check against direct Fourier inversion of g_e_hat(k) * sin(kt)/k
    # g_e_hat(k) = int_{-1}^1 x e^{-ikx} dx = -2i (sin k - k cos k) / k^2
    def g_e_hat(k):
        if abs(k) < 1e-7:
            return -2j * (k / 3.0 - (k ** 3) / 30.0)
        return -2j * (math.sin(k) - k * math.cos(k)) / (k * k)

    max_fe_err = 0.0
    for (x, t) in [(0.0, 2.0), (0.5, 2.0), (-0.8, 2.0), (2.0, 2.0), (-2.0, 2.0), (3.5, 2.0)]:
        def integrand(k):
            sinc_t = t if abs(k) < 1e-12 else math.sin(k * t) / k
            return g_e_hat(k) * sinc_t
        val_ft = ift(integrand, x, K=120.0, n=6000).real
        val_ex = u_e_exact(x, t)
        max_fe_err = max(max_fe_err, abs(val_ft - val_ex))
    assert max_fe_err < 1e-3, f"Part (e) FT inversion err: {max_fe_err}"
    print(f"[PASS] Part (e) wedge |x| <= t - 1 (t > 1) verified; FT cross-check max err = {max_fe_err:.3e}")
    print("ALL CHECKS PASSED for IB/2023/P3/q14A.")


if __name__ == "__main__":
    main()
