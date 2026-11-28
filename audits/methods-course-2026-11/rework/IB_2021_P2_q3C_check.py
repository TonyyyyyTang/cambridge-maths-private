#!/usr/bin/env python3
"""Independent verification for IB/2021/P2/q3C (Sturm-Liouville eigenvalue problem)."""

import cmath
import math


def simpson(func, a, b, n=2000):
    """Composite Simpson's rule on [a, b] with even n."""
    assert n % 2 == 0
    h = (b - a) / n
    s = func(a) + func(b)
    for k in range(1, n):
        x = a + k * h
        s += (4.0 if k % 2 == 1 else 2.0) * func(x)
    return s * h / 3.0


def y_n(n, x):
    return math.exp(-x) * math.sin(n * math.pi * x)


def lam_n(n):
    return 1.0 + (n * math.pi) ** 2


def w(x):
    return math.exp(2.0 * x)


def f_target(x):
    return x - x * x


def A_n_exact(n):
    """Closed-form evaluation of 2 * int_0^1 (x - x^2) e^x sin(n pi x) dx for cross-check."""
    alpha = complex(1.0, n * math.pi)
    val = ((alpha - 2.0) * cmath.exp(alpha) + (alpha + 2.0)) / (alpha ** 3)
    return 2.0 * val.imag


def main():
    # 1. Verify Sturm-Liouville identity: d/dx(e^{2x} y') == e^{2x} (y'' + 2 y')
    h = 1e-3
    def test_u(x):
        return math.sin(2.3 * x) + 0.7 * x ** 3 - 0.4 * math.exp(-0.5 * x)

    max_sl_err = 0.0
    for i in range(1, 20):
        x = i / 20.0
        # 4th-order central differences
        u_p = (-test_u(x + 2 * h) + 8 * test_u(x + h) - 8 * test_u(x - h) + test_u(x - 2 * h)) / (12 * h)
        u_pp = (-test_u(x + 2 * h) + 16 * test_u(x + h) - 30 * test_u(x) + 16 * test_u(x - h) - test_u(x - 2 * h)) / (12 * h * h)
        rhs = w(x) * (u_pp + 2.0 * u_p)

        def p_u_p(t):
            tp = (-test_u(t + 2 * h) + 8 * test_u(t + h) - 8 * test_u(t - h) + test_u(t - 2 * h)) / (12 * h)
            return w(t) * tp

        lhs = (-p_u_p(x + 2 * h) + 8 * p_u_p(x + h) - 8 * p_u_p(x - h) + p_u_p(x - 2 * h)) / (12 * h)
        max_sl_err = max(max_sl_err, abs(lhs - rhs))
    assert max_sl_err < 1e-7, f"SL identity error: {max_sl_err}"
    print(f"[PASS] Sturm-Liouville form identity max err = {max_sl_err:.3e}")

    # 2. Verify boundary conditions and eigenvalue equation for n = 1..6
    max_ode_err = 0.0
    for n in range(1, 7):
        assert abs(y_n(n, 0.0)) < 1e-15
        assert abs(y_n(n, 1.0)) < 1e-14
        lam = lam_n(n)
        for i in range(1, 20):
            x = i / 20.0
            yp = (-y_n(n, x + 2 * h) + 8 * y_n(n, x + h) - 8 * y_n(n, x - h) + y_n(n, x - 2 * h)) / (12 * h)
            ypp = (-y_n(n, x + 2 * h) + 16 * y_n(n, x + h) - 30 * y_n(n, x) + 16 * y_n(n, x - h) - y_n(n, x - 2 * h)) / (12 * h * h)
            res = ypp + 2.0 * yp + lam * y_n(n, x)
            max_ode_err = max(max_ode_err, abs(res))
    assert max_ode_err < 1e-6, f"ODE residual too large: {max_ode_err}"
    print(f"[PASS] Eigenvalue ODE & BCs (n=1..6) max residual = {max_ode_err:.3e}")

    # 3. Verify strictly increasing eigenvalues and no eigenvalues for lambda <= 1
    lams = [lam_n(n) for n in range(1, 10)]
    assert all(lams[i] < lams[i + 1] for i in range(len(lams) - 1))
    # For lambda < 1, solution with y(0)=0, y'(0)=1 is e^{-x} sinh(mu x)/mu; at x=1 this is > 0
    for lam in [-10.0, -2.0, -0.5, 0.0, 0.5, 0.999]:
        mu = math.sqrt(1.0 - lam)
        y1 = math.exp(-1.0) * math.sinh(mu) / mu
        assert y1 > 0.0
    # For lambda == 1, solution with y(0)=0, y'(0)=1 is x e^{-x}; at x=1 this is e^{-1} > 0
    assert math.exp(-1.0) > 0.0
    print(f"[PASS] No eigenvalues for lambda <= 1; lambda_1 = {lam_n(1):.6f}, lambda_2 = {lam_n(2):.6f}")

    # 4. Verify orthogonality with weight w(x) = e^{2x}
    max_orth_err = 0.0
    for m in range(1, 7):
        for n in range(1, 7):
            ip = simpson(lambda x: w(x) * y_n(m, x) * y_n(n, x), 0.0, 1.0, n=1000)
            expected = 0.5 if m == n else 0.0
            max_orth_err = max(max_orth_err, abs(ip - expected))
    assert max_orth_err < 1e-11, f"Orthogonality error: {max_orth_err}"
    print(f"[PASS] Weighted orthogonality <y_m, y_n>_w = (1/2) delta_mn max err = {max_orth_err:.3e}")

    # 5. Verify A_n formula: ratio of integrals == 2 * int_0^1 (x-x^2) e^x sin(n pi x) dx == closed form
    max_an_err = 0.0
    for n in range(1, 11):
        num = simpson(lambda x: w(x) * f_target(x) * y_n(n, x), 0.0, 1.0, n=2000)
        den = simpson(lambda x: w(x) * (y_n(n, x) ** 2), 0.0, 1.0, n=2000)
        an_ratio = num / den
        an_simp = 2.0 * simpson(lambda x: f_target(x) * math.exp(x) * math.sin(n * math.pi * x), 0.0, 1.0, n=2000)
        an_ex = A_n_exact(n)
        max_an_err = max(max_an_err, abs(an_ratio - an_simp), abs(an_simp - an_ex))
    assert max_an_err < 1e-10, f"A_n mismatch: {max_an_err}"
    print(f"[PASS] A_n integral vs closed form (n=1..10) max err = {max_an_err:.3e}")

    # 6. Verify eigenfunction expansion converges to f(x) = x - x^2 on [0, 1]
    N_terms = 200
    A_coeffs = [A_n_exact(n) for n in range(1, N_terms + 1)]
    max_recon_err = 0.0
    for i in range(0, 101):
        x = i / 100.0
        approx = sum(A_coeffs[n - 1] * y_n(n, x) for n in range(1, N_terms + 1))
        max_recon_err = max(max_recon_err, abs(approx - f_target(x)))
    assert max_recon_err < 1e-4, f"Series reconstruction error too large: {max_recon_err}"
    print(f"[PASS] Eigenfunction expansion (N={N_terms}) max error on [0,1] = {max_recon_err:.3e}")
    print("ALL CHECKS PASSED for IB/2021/P2/q3C.")


if __name__ == "__main__":
    main()
