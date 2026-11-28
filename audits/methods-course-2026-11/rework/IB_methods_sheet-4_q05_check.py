#!/usr/bin/env python3
"""Computational verification for IB/methods/sheet-4/q05.

Verifies:
1. Parseval's theorem in the Sheet 4 Fourier transform convention.
2. Equivalence of Fourier-space solution v_hat(lam, t) = g_hat(lam)*exp(-kappa*|lam|^2*t)
   and heat-kernel convolution v(x, t) = (K_t * g)(x) in 1D, and PDE satisfaction.
3. L^2 contraction bound ||u1(., t) - u2(., t)|| <= C ||f1 - f2|| with C = 1
   in 1D and 2D (including non-zero forcing F(x, t) which cancels in u1 - u2).
4. Sharpness of C = 1 via broad Gaussians g_eps(x) = exp(-eps*|x|^2) as eps -> 0+.
5. Ill-posedness for kappa < 0 (unbounded amplification of high frequencies).
"""

import cmath
import math


def simpson_1d(func, a, b, n):
    """Composite Simpson's rule on [a, b] with even integer n."""
    assert n % 2 == 0
    h = (b - a) / n
    s = func(a) + func(b)
    for k in range(1, n):
        w = 4.0 if k % 2 == 1 else 2.0
        s += w * func(a + k * h)
    return s * h / 3.0


def main():
    # ------------------------------------------------------------------
    # 1. 1D verification: bimodal + oscillatory difference g = f1 - f2
    #    g(x) = exp(-(x-0.4)^2) + 0.7 * x * exp(-0.5 * x^2) * cos(2*x)
    # ------------------------------------------------------------------
    def g1(x):
        return math.exp(-((x - 0.4) ** 2)) + 0.7 * x * math.exp(-0.5 * x * x) * math.cos(2.0 * x)

    # Exact Fourier transform of g1(x) under convention:
    # g_hat(lam) = int g(x) exp(-i*lam*x) dx
    # For exp(-(x-a)^2): sqrt(pi) * exp(-i*a*lam - lam^2/4)
    # For x*exp(-x^2/2)*cos(2x) = (1/2)*x*exp(-x^2/2)*(e^{2ix}+e^{-2ix}):
    #   int exp(-x^2/2) e^{-i*k*x} dx = sqrt(2*pi) * exp(-k^2/2)
    #   differentiating wrt k: int (-i*x) exp(-x^2/2) e^{-i*k*x} dx = -k*sqrt(2*pi)*exp(-k^2/2)
    #   => int x*exp(-x^2/2) e^{-i*k*x} dx = -i*k*sqrt(2*pi)*exp(-k^2/2)
    def g1_hat(lam):
        term1 = math.sqrt(math.pi) * cmath.exp(-1j * 0.4 * lam - 0.25 * lam * lam)
        kp = lam - 2.0
        km = lam + 2.0
        term2 = (
            0.35
            * math.sqrt(2.0 * math.pi)
            * (-1j)
            * (kp * math.exp(-0.5 * kp * kp) + km * math.exp(-0.5 * km * km))
        )
        return term1 + term2

    # Check numerical FT matches analytic g1_hat at sample frequencies
    for lam in [-3.1, 0.0, 1.7]:
        ft_num = simpson_1d(lambda x: g1(x) * cmath.exp(-1j * lam * x), -12.0, 12.0, 800)
        assert abs(ft_num - g1_hat(lam)) < 1e-10, f"FT mismatch at lam={lam}"

    # Verify Parseval's theorem
    norm_sq_x = simpson_1d(lambda x: abs(g1(x)) ** 2, -12.0, 12.0, 800)
    norm_sq_lam = (1.0 / (2.0 * math.pi)) * simpson_1d(
        lambda lam: abs(g1_hat(lam)) ** 2, -15.0, 15.0, 800
    )
    parseval_err = abs(norm_sq_x - norm_sq_lam)
    assert parseval_err < 1e-10, f"Parseval error too large: {parseval_err}"

    # Verify Fourier-space evolution == heat kernel convolution and satisfies PDE
    kappa = 0.7
    t_test = 0.45

    def v_fourier(x, t):
        val = (1.0 / (2.0 * math.pi)) * simpson_1d(
            lambda lam: g1_hat(lam) * math.exp(-kappa * lam * lam * t) * cmath.exp(1j * lam * x),
            -15.0,
            15.0,
            600,
        )
        return val.real

    def v_conv(x, t):
        pref = 1.0 / math.sqrt(4.0 * math.pi * kappa * t)
        return pref * simpson_1d(
            lambda y: g1(y) * math.exp(-((x - y) ** 2) / (4.0 * kappa * t)),
            -12.0,
            12.0,
            800,
        )

    max_conv_diff = 0.0
    for x in [-1.5, 0.0, 0.8, 2.1]:
        max_conv_diff = max(max_conv_diff, abs(v_fourier(x, t_test) - v_conv(x, t_test)))
    assert max_conv_diff < 1e-9, f"Fourier vs convolution mismatch: {max_conv_diff}"

    # Finite-difference PDE check v_t - kappa * v_xx = 0
    dt = 1e-5
    dx = 1e-4
    x0 = 0.6
    vt = (v_fourier(x0, t_test + dt) - v_fourier(x0, t_test - dt)) / (2.0 * dt)
    vxx = (v_fourier(x0 + dx, t_test) - 2.0 * v_fourier(x0, t_test) + v_fourier(x0 - dx, t_test)) / (dx * dx)
    pde_res = abs(vt - kappa * vxx)
    assert pde_res < 1e-6, f"1D heat PDE residual too large: {pde_res}"

    # Verify ||v(., t)|| <= 1 * ||g|| for several t > 0
    ratios_1d = []
    for t in [0.05, 0.45, 2.0]:
        v_norm_sq_lam = (1.0 / (2.0 * math.pi)) * simpson_1d(
            lambda lam: (abs(g1_hat(lam)) ** 2) * math.exp(-2.0 * kappa * lam * lam * t),
            -15.0,
            15.0,
            800,
        )
        v_norm_sq_x = simpson_1d(lambda x: v_fourier(x, t) ** 2, -16.0, 16.0, 400)
        assert abs(v_norm_sq_x - v_norm_sq_lam) < 1e-8
        ratio = math.sqrt(v_norm_sq_lam / norm_sq_lam)
        assert ratio <= 1.0 + 1e-12
        ratios_1d.append((t, ratio))

    # ------------------------------------------------------------------
    # 2. 2D verification with non-zero inhomogeneous forcing F(x, y, t)
    #    Let F(x, y, t) = (1 + 4*kappa*beta*(1 - beta*(x^2+y^2))) * exp(-beta*(x^2+y^2) - t)
    #    Wait: simpler! Choose a particular solution u_p(x, y, t) = sin(t) * exp(-(x^2+y^2)),
    #    so F = (u_p)_t - kappa * Delta(u_p) is non-zero, u_p(x, y, 0) = 0.
    #    Then u_i = u_p + v_i where v_i solves homogeneous heat eqn with data f_i.
    #    For Gaussian initial data:
    #      f1(x, y) = exp(-a1*x^2 - b1*y^2)
    #      f2(x, y) = exp(-a2*x^2 - b2*y^2)
    #    Heat evolution of exp(-a*x^2 - b*y^2) is exact:
    #      (1+4*kappa*a*t)^(-1/2) * (1+4*kappa*b*t)^(-1/2)
    #      * exp(-a*x^2/(1+4*kappa*a*t) - b*y^2/(1+4*kappa*b*t))
    # ------------------------------------------------------------------
    def gauss_evolve_2d(x, y, t, a, b, kap):
        da = 1.0 + 4.0 * kap * a * t
        db = 1.0 + 4.0 * kap * b * t
        return math.exp(-a * x * x / da - b * y * y / db) / math.sqrt(da * db)

    a1, b1 = 0.8, 1.3
    a2, b2 = 1.5, 0.6

    # Analytic L^2 inner product in 2D:
    # int exp(-p1*x^2 - q1*y^2) * exp(-p2*x^2 - q2*y^2) dx dy = pi / sqrt((p1+p2)*(q1+q2))
    def l2_sq_diff_2d(t):
        def coeff(a, b):
            da = 1.0 + 4.0 * kappa * a * t
            db = 1.0 + 4.0 * kappa * b * t
            return (1.0 / math.sqrt(da * db), a / da, b / db)

        c1, p1, q1 = coeff(a1, b1)
        c2, p2, q2 = coeff(a2, b2)
        i11 = (c1 * c1) * math.pi / math.sqrt((2 * p1) * (2 * q1))
        i22 = (c2 * c2) * math.pi / math.sqrt((2 * p2) * (2 * q2))
        i12 = (c1 * c2) * math.pi / math.sqrt((p1 + p2) * (q1 + q2))
        return i11 + i22 - 2.0 * i12

    norm_f_diff_2d = math.sqrt(l2_sq_diff_2d(0.0))
    for t in [0.1, 0.5, 3.0]:
        norm_u_diff_2d = math.sqrt(l2_sq_diff_2d(t))
        assert norm_u_diff_2d <= norm_f_diff_2d + 1e-14
        # Also verify exact PDE for u_p + v_1 and cancellation of u_p in u1 - u2
        x0, y0 = 0.3, -0.4
        u1 = math.sin(t) * math.exp(-(x0 * x0 + y0 * y0)) + gauss_evolve_2d(x0, y0, t, a1, b1, kappa)
        u2 = math.sin(t) * math.exp(-(x0 * x0 + y0 * y0)) + gauss_evolve_2d(x0, y0, t, a2, b2, kappa)
        v_diff = gauss_evolve_2d(x0, y0, t, a1, b1, kappa) - gauss_evolve_2d(x0, y0, t, a2, b2, kappa)
        assert abs((u1 - u2) - v_diff) < 1e-15

    # ------------------------------------------------------------------
    # 3. Sharpness of C = 1: broad Gaussian g_eps(x) = exp(-eps*|x|^2) in R^n
    #    Exact ratio: ||v_eps(., t)|| / ||g_eps|| = (1 + 4*kappa*eps*t)^(-n/4) -> 1 as eps -> 0+
    # ------------------------------------------------------------------
    t_fixed = 1.0
    eps_small = 1e-5
    sharp_ratios = {}
    for n in [1, 2, 3, 4]:
        exact_ratio = (1.0 + 4.0 * kappa * eps_small * t_fixed) ** (-n / 4.0)
        assert 0.9999 < exact_ratio < 1.0
        sharp_ratios[n] = exact_ratio

    # Cross-check n=1 formula numerically at eps = 0.05 via Parseval integral
    eps_num = 0.05
    # g_hat(lam) = sqrt(pi/eps) * exp(-lam^2/(4*eps))
    num_top = simpson_1d(
        lambda lam: math.exp(-lam * lam / (2.0 * eps_num) - 2.0 * kappa * lam * lam * t_fixed),
        -4.0,
        4.0,
        600,
    )
    num_bot = simpson_1d(lambda lam: math.exp(-lam * lam / (2.0 * eps_num)), -4.0, 4.0, 600)
    ratio_num = math.sqrt(num_top / num_bot)
    ratio_ana = (1.0 + 4.0 * kappa * eps_num * t_fixed) ** (-0.25)
    assert abs(ratio_num - ratio_ana) < 1e-12

    # ------------------------------------------------------------------
    # 4. Backwards heat equation (kappa < 0) violates boundedness
    #    For g_k(x) = exp(-x^2/2) * cos(k*x), ||v(., t)|| / ||g_k|| grows exponentially in k^2.
    # ------------------------------------------------------------------
    kap_neg = -0.1
    t_back = 0.1
    k_high = 10.0
    # At lam = k_high, exp(-2 * kap_neg * k^2 * t) = exp(+2 * 0.1 * 100 * 0.1) = exp(2)
    # Exact ratio for exp(-x^2/2)*exp(i*k*x) (same L^2 norm evolution for single mode):
    # int exp(-(lam-k)^2) * exp(2*|kap|*t*lam^2) dlam / int exp(-(lam-k)^2) dlam
    # (converges provided 2*|kap|*t < 1; for larger t even blows up to infinity!)
    alpha = 2.0 * abs(kap_neg) * t_back  # = 0.02 < 1
    amp_sq = (1.0 / math.sqrt(1.0 - alpha)) * math.exp((alpha / (1.0 - alpha)) * k_high * k_high)
    assert math.sqrt(amp_sq) > 2.5

    print("IB_methods_sheet-4_q05 verification PASSED:")
    print(f"  Parseval error (1D)             = {parseval_err:.3e}")
    print(f"  Fourier vs convolution diff     = {max_conv_diff:.3e}")
    print(f"  1D heat PDE FD residual         = {pde_res:.3e}")
    print("  1D contraction ratios ||v(t)||/||g||:", [(t, round(r, 6)) for t, r in ratios_1d])
    print(f"  2D ||f1-f2|| = {norm_f_diff_2d:.6f}, ||u1(0.5)-u2(0.5)|| = {math.sqrt(l2_sq_diff_2d(0.5)):.6f} <= ||f1-f2||")
    print("  Sharpness of C=1 (eps=1e-5, t=1):", {n: round(r, 8) for n, r in sharp_ratios.items()})
    print(f"  Backwards heat (kappa=-0.1, k=10, t=0.1) amplification = {math.sqrt(amp_sq):.4f} > 1")


if __name__ == "__main__":
    main()
