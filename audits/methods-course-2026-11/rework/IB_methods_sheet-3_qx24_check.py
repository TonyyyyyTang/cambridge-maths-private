#!/usr/bin/env python3
"""Verification script for IB/methods/sheet-3/qx24 (Gaussian FT and delta sequence)."""

import cmath
import math


def composite_simpson_complex(f, a, b, n):
    assert n % 2 == 0
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * f(a + i * h)
    return s * h / 3.0


def main():
    # Convention: f_tilde(k) = int_{-inf}^{inf} f(x) exp(-i k x) dx
    # (confirmed by Q7 on the same Michaelmas 2016 sheet).

    # 1. Verify FT of f(x) = exp(-n^2 (x - mu)^2) for various n (pos and neg), mu, k
    max_ft_err = 0.0
    test_cases = [
        (0.7, 0.0),
        (1.3, 2.1),
        (2.5, -1.4),
        (-1.8, 0.9),  # tests |n| when n < 0
    ]
    for n, mu in test_cases:
        abs_n = abs(n)
        half_w = 10.0 / abs_n  # exp(-100) ~ 3.7e-44
        for k in [-3.5, -1.0, 0.0, 0.8, 2.7]:
            num_ft = composite_simpson_complex(
                lambda x, n=n, mu=mu, k=k: math.exp(-(n**2) * (x - mu) ** 2) * cmath.exp(-1j * k * x),
                mu - half_w,
                mu + half_w,
                2000,
            )
            exact_ft = (math.sqrt(math.pi) / abs_n) * cmath.exp(-1j * k * mu - (k**2) / (4.0 * n**2))
            max_ft_err = max(max_ft_err, abs(num_ft - exact_ft))
    assert max_ft_err < 1e-11, f"FT error {max_ft_err}"
    print(f"[PASS] Gaussian FT verified across positive/negative n, mu, k (max err={max_ft_err:.2e})")

    # 2. Verify Fourier inversion f(x) = (1/(2pi)) int f_tilde(k) exp(+i k x) dk
    n, mu = 1.4, -0.75
    K_cut = 20.0 * n  # exp(-100)
    max_inv_err = 0.0
    for x in [-2.0, -0.75, 0.0, 0.5, 1.8]:
        f_inv = (1.0 / (2.0 * math.pi)) * composite_simpson_complex(
            lambda k, x=x: (math.sqrt(math.pi) / n)
            * cmath.exp(-1j * k * mu - (k**2) / (4.0 * n**2))
            * cmath.exp(1j * k * x),
            -K_cut,
            K_cut,
            2000,
        )
        f_exact = math.exp(-(n**2) * (x - mu) ** 2)
        max_inv_err = max(max_inv_err, abs(f_inv - f_exact))
    assert max_inv_err < 1e-11, f"Inversion error {max_inv_err}"
    print(f"[PASS] Fourier inversion formula verified (max err={max_inv_err:.2e})")

    # 3. Verify differentiation ODE d(g_tilde)/dk = -(k / (2 n^2)) g_tilde(k)
    n_ode = 1.6
    half_w = 10.0 / n_ode

    def g_tilde_num(k):
        return composite_simpson_complex(
            lambda x: math.exp(-(n_ode**2) * x * x) * cmath.exp(-1j * k * x),
            -half_w,
            half_w,
            1000,
        ).real

    hk = 1e-5
    max_ode_err = 0.0
    for k in [-2.5, -0.5, 0.0, 1.2, 3.0]:
        dg_dk = (g_tilde_num(k + hk) - g_tilde_num(k - hk)) / (2.0 * hk)
        rhs = -(k / (2.0 * n_ode * n_ode)) * g_tilde_num(k)
        max_ode_err = max(max_ode_err, abs(dg_dk - rhs))
    assert max_ode_err < 1e-8, f"ODE error {max_ode_err}"
    print(f"[PASS] Differentiation ODE dg_tilde/dk = -(k/(2n^2)) g_tilde verified (max err={max_ode_err:.2e})")

    # 4. Verify integral of delta_n(x) = 1 and FT delta_n_tilde(k) = exp(-k^2 / (4 n^2))
    max_area_err = 0.0
    max_d_ft_err = 0.0
    for n_val in [0.1, 0.5, 2.0, 10.0, 50.0]:
        hw = 10.0 / n_val
        area = composite_simpson_complex(
            lambda x, nv=n_val: (nv / math.sqrt(math.pi)) * math.exp(-(nv**2) * x * x),
            -hw,
            hw,
            2000,
        ).real
        max_area_err = max(max_area_err, abs(area - 1.0))

        for k in [0.0, 0.5 * n_val, 2.0 * n_val]:
            d_ft_num = composite_simpson_complex(
                lambda x, nv=n_val, k=k: (nv / math.sqrt(math.pi))
                * math.exp(-(nv**2) * x * x)
                * cmath.exp(-1j * k * x),
                -hw,
                hw,
                2000,
            )
            d_ft_ex = math.exp(-(k**2) / (4.0 * n_val**2))
            max_d_ft_err = max(max_d_ft_err, abs(d_ft_num - d_ft_ex))
    assert max_area_err < 1e-12
    assert max_d_ft_err < 1e-11
    print(
        f"[PASS] Area int delta_n(x) dx = 1 (max err={max_area_err:.2e}) "
        f"and delta_n_tilde(k) = exp(-k^2/(4n^2)) (max err={max_d_ft_err:.2e})"
    )

    # 5. Verify sifting limit int delta_n(x) phi(x) dx -> phi(0) and delta_n_tilde(k) -> 1 as n -> inf
    test_phis = [
        ("cos(3x) + x^2", lambda x: math.cos(3.0 * x) + x * x),
        ("exp(-x)*sin(2x + 0.5)", lambda x: math.exp(-x) * math.sin(2.0 * x + 0.5)),
        ("1/(1 + x^2)", lambda x: 1.0 / (1.0 + x * x)),
    ]
    for label, phi in test_phis:
        prev_err = float("inf")
        for n_val in [5.0, 20.0, 100.0, 500.0]:
            hw = 10.0 / n_val
            val = composite_simpson_complex(
                lambda x, nv=n_val, phi=phi: (nv / math.sqrt(math.pi)) * math.exp(-(nv**2) * x * x) * phi(x),
                -hw,
                hw,
                2000,
            ).real
            err = abs(val - phi(0.0))
            assert err < prev_err
            prev_err = err
        assert prev_err < 1e-4, f"Sifting convergence too slow for {label}: {prev_err}"

    # Uniform convergence of delta_n_tilde(k) -> 1 on [-10, 10]
    n_large = 1000.0
    max_k_err = max(abs(math.exp(-(k**2) / (4.0 * n_large**2)) - 1.0) for k in [-10.0, -5.0, 0.0, 5.0, 10.0])
    assert max_k_err < 3e-5
    print(
        f"[PASS] Distributional limit delta_n -> delta verified on 3 test functions "
        f"and max |delta_1000_tilde(k) - 1| on [-10,10] = {max_k_err:.2e}"
    )
    print("ALL CHECKS PASSED for IB_methods_sheet-3_qx24.")


if __name__ == "__main__":
    main()
