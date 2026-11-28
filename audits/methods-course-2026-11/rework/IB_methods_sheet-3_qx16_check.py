#!/usr/bin/env python3
"""Verification script for IB/methods/sheet-3/qx16 (Modes on a string)."""

from fractions import Fraction
import math


def composite_simpson(f, a, b, n):
    assert n % 2 == 0
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * f(a + i * h)
    return s * h / 3.0


def main():
    # Arbitrary non-trivial physical parameters
    mu = 0.35
    T = 14.0
    ell = 2.5
    V = 1.8
    c = math.sqrt(T / mu)

    def y_t0(x):
        return (4.0 * V / (ell * ell)) * x * (ell - x)

    # 1. Verify Fourier sine coefficients omega_n * B_n by direct quadrature
    max_coeff_err = 0.0
    for n in range(1, 13):
        num_val = (2.0 / ell) * composite_simpson(
            lambda x, n=n: y_t0(x) * math.sin(n * math.pi * x / ell),
            0.0,
            ell,
            2000,
        )
        exact_val = (32.0 * V / (n**3 * math.pi**3)) if (n % 2 == 1) else 0.0
        max_coeff_err = max(max_coeff_err, abs(num_val - exact_val))
    assert max_coeff_err < 1e-10, f"Fourier coeff err {max_coeff_err}"
    print(f"[PASS] Fourier sine coefficients omega_n*B_n verified for n=1..12 (max err={max_coeff_err:.2e})")

    # Build modal representation (odd n up to N_max)
    N_max = 199
    odd_modes = list(range(1, N_max + 1, 2))

    def omega(n):
        return n * math.pi * c / ell

    def k_wave(n):
        return n * math.pi / ell

    def B_amp(n):
        return 32.0 * V * ell / (n**4 * math.pi**4 * c)

    def y_val(x, t):
        return sum(B_amp(n) * math.sin(omega(n) * t) * math.sin(k_wave(n) * x) for n in odd_modes)

    def y_t(x, t):
        return sum(omega(n) * B_amp(n) * math.cos(omega(n) * t) * math.sin(k_wave(n) * x) for n in odd_modes)

    def y_x(x, t):
        return sum(k_wave(n) * B_amp(n) * math.sin(omega(n) * t) * math.cos(k_wave(n) * x) for n in odd_modes)

    # 2. Check BCs, ICs, and PDE (via finite differences on y)
    for t in [0.0, 0.13, 0.47, 1.25]:
        assert abs(y_val(0.0, t)) < 1e-15
        assert abs(y_val(ell, t)) < 1e-13

    max_ic0_err = 0.0
    max_ic1_err = 0.0
    for j in range(1, 50):
        x = ell * j / 50.0
        max_ic0_err = max(max_ic0_err, abs(y_val(x, 0.0)))
        max_ic1_err = max(max_ic1_err, abs(y_t(x, 0.0) - y_t0(x)))
    assert max_ic0_err < 1e-15
    assert max_ic1_err < 1e-5  # truncation tail sum_{n>199} n^-3 ~ 1e-5
    print(f"[PASS] BCs and ICs verified (max |y(x,0)|={max_ic0_err:.2e}, max |y_t(x,0)-f(x)|={max_ic1_err:.2e})")

    # Finite-difference PDE check mu * y_tt = T * y_xx (using first 25 odd modes so smooth)
    modes_pde = list(range(1, 50, 2))

    def y_smooth(x, t):
        return sum(B_amp(n) * math.sin(omega(n) * t) * math.sin(k_wave(n) * x) for n in modes_pde)

    h_x, h_t = 1e-5, 1e-5
    max_pde_res = 0.0
    for x in [0.3, 1.1, 2.0]:
        for t in [0.1, 0.37, 0.9]:
            ytt = (y_smooth(x, t + h_t) - 2.0 * y_smooth(x, t) + y_smooth(x, t - h_t)) / (h_t * h_t)
            yxx = (y_smooth(x + h_x, t) - 2.0 * y_smooth(x, t) + y_smooth(x - h_x, t)) / (h_x * h_x)
            max_pde_res = max(max_pde_res, abs(mu * ytt - T * yxx))
    assert max_pde_res < 1e-4, f"PDE residual {max_pde_res}"
    print(f"[PASS] Wave PDE mu*y_tt = T*y_xx verified by finite differences (max res={max_pde_res:.2e})")

    # 3. Verify kinetic energy K(t), potential energy U(t), and energy conservation E(t) = E(0)
    E0_exact = (4.0 / 15.0) * mu * ell * V * V
    E0_quad = 0.5 * mu * composite_simpson(lambda x: y_t0(x) ** 2, 0.0, ell, 2000)
    assert abs(E0_quad - E0_exact) < 1e-12

    pref = 256.0 * mu * ell * V * V / (math.pi**6)
    max_K_err = 0.0
    max_U_err = 0.0
    max_E_err = 0.0
    for t in [0.0, 0.07, 0.21, 0.53, 1.17]:
        K_modal = pref * sum((math.cos(omega(n) * t) ** 2) / (n**6) for n in odd_modes)
        U_modal = pref * sum((math.sin(omega(n) * t) ** 2) / (n**6) for n in odd_modes)
        K_quad = 0.5 * mu * composite_simpson(lambda x, t=t: y_t(x, t) ** 2, 0.0, ell, 1000)
        U_quad = 0.5 * T * composite_simpson(lambda x, t=t: y_x(x, t) ** 2, 0.0, ell, 1000)
        max_K_err = max(max_K_err, abs(K_modal - K_quad))
        max_U_err = max(max_U_err, abs(U_modal - U_quad))
        max_E_err = max(max_E_err, abs((K_modal + U_modal) - E0_exact))
    assert max_K_err < 1e-8, f"K err {max_K_err}"
    assert max_U_err < 1e-8, f"U err {max_U_err}"
    assert max_E_err < 1e-11, f"E conservation err {max_E_err}"
    print(
        f"[PASS] Modal K(t), U(t) match spatial quadrature (max K err={max_K_err:.2e}, "
        f"max U err={max_U_err:.2e}) and E(t)=E(0)=4/15*mu*ell*V^2 (err={max_E_err:.2e})"
    )

    # 4. Exact rational algebra & numerical series check for sum_{n odd} 1/n^6 = pi^6 / 960
    int_u2_1mu2 = Fraction(1, 3) - Fraction(2, 4) + Fraction(1, 5)
    assert int_u2_1mu2 == Fraction(1, 30)
    E0_coeff = Fraction(8, 1) * int_u2_1mu2
    assert E0_coeff == Fraction(4, 15)
    ratio_pi6 = E0_coeff / Fraction(256, 1)
    assert ratio_pi6 == Fraction(1, 960)

    series_sum = sum(1.0 / (n**6) for n in range(1, 20000, 2))
    target_sum = (math.pi**6) / 960.0
    assert abs(series_sum - target_sum) < 1e-14
    print(
        f"[PASS] Series identity verified: sum_(n odd) 1/n^6 = {series_sum:.15f}, "
        f"pi^6/960 = {target_sum:.15f} (diff={abs(series_sum - target_sum):.2e})"
    )
    print("ALL CHECKS PASSED for IB_methods_sheet-3_qx16.")


if __name__ == "__main__":
    main()
