#!/usr/bin/env python3
"""Independent verification for IB/2023/P2/q14A (Laplace's equation in polar coordinates & annulus)."""

import math


def simpson(func, a, b, n=2000):
    assert n % 2 == 0
    h = (b - a) / n
    s = func(a) + func(b)
    for k in range(1, n):
        x = a + k * h
        s += (4.0 if k % 2 == 1 else 2.0) * func(x)
    return s * h / 3.0


def polar_laplacian_fd(phi, r, th, hr=1e-3, hth=1e-3):
    """4th-order central finite-difference polar Laplacian: (1/r)(r phi_r)_r + (1/r^2) phi_{th,th}."""
    phi_r = (-phi(r + 2 * hr, th) + 8 * phi(r + hr, th) - 8 * phi(r - hr, th) + phi(r - 2 * hr, th)) / (12 * hr)
    phi_rr = (-phi(r + 2 * hr, th) + 16 * phi(r + hr, th) - 30 * phi(r, th) + 16 * phi(r - hr, th) - phi(r - 2 * hr, th)) / (12 * hr * hr)
    phi_tt = (-phi(r, th + 2 * hth) + 16 * phi(r, th + hth) - 30 * phi(r, th) + 16 * phi(r, th - hth) - phi(r, th - 2 * hth)) / (12 * hth * hth)
    return phi_rr + phi_r / r + phi_tt / (r * r)


def f_specific(th):
    th = th % (2.0 * math.pi)
    if th < math.pi:
        return 0.5 * math.pi - th
    return -1.5 * math.pi + th


def main():
    # 1. Verify general polar Laplace basis modes satisfy nabla^2 phi = 0
    basis = [
        lambda r, t: 1.0,
        lambda r, t: math.log(r),
    ]
    for n in range(1, 5):
        basis.extend([
            lambda r, t, n=n: (r ** n) * math.cos(n * t),
            lambda r, t, n=n: (r ** (-n)) * math.cos(n * t),
            lambda r, t, n=n: (r ** n) * math.sin(n * t),
            lambda r, t, n=n: (r ** (-n)) * math.sin(n * t),
        ])
    max_basis_err = 0.0
    for fn in basis:
        for r in [1.3, 2.7, 5.1]:
            for t in [0.2, 1.7, 4.3]:
                max_basis_err = max(max_basis_err, abs(polar_laplacian_fd(fn, r, t)))
    assert max_basis_err < 1e-6, f"Basis Laplacian residual: {max_basis_err}"
    print(f"[PASS] Part (a) basis functions satisfy nabla^2 phi = 0, max err = {max_basis_err:.3e}")

    # 2. Verify part (b)(i) for a general finite cosine series f(theta) = sum_{n=1}^5 A_n cos(n theta)
    A_test = {1: 0.7, 2: -1.3, 3: 0.45, 4: -0.2, 5: 0.11}

    def f_gen(t):
        return sum(A * math.cos(n * t) for n, A in A_test.items())

    def phi_gen_form1(r, t):
        s = math.log(r) - 1.0
        for n, A in A_test.items():
            s += A * (r ** n + math.exp(2.0 * n) * (r ** (-n))) / (math.exp(2.0 * n) + 1.0) * math.cos(n * t)
        return s

    def phi_gen_form2(r, t):
        s = math.log(r) - 1.0
        lr = math.log(r)
        for n, A in A_test.items():
            s += A * (math.cosh(n * (lr - 1.0)) / math.cosh(n)) * math.cos(n * t)
        return s

    e2 = math.exp(2.0)
    max_equiv_err = 0.0
    max_pde_gen = 0.0
    max_bc_gen = 0.0
    for r in [1.2, 2.5, 4.0, 6.5]:
        for t in [0.0, 0.6, 1.9, 3.5, 5.2]:
            max_equiv_err = max(max_equiv_err, abs(phi_gen_form1(r, t) - phi_gen_form2(r, t)))
            max_pde_gen = max(max_pde_gen, abs(polar_laplacian_fd(phi_gen_form1, r, t)))
    for k in range(32):
        t = 2.0 * math.pi * k / 32.0
        max_bc_gen = max(
            max_bc_gen,
            abs(phi_gen_form1(1.0, t) - (f_gen(t) - 1.0)),
            abs(phi_gen_form1(e2, t) - (f_gen(t) + 1.0)),
        )
    assert max_equiv_err < 1e-14
    assert max_pde_gen < 1e-6
    assert max_bc_gen < 1e-13
    print(f"[PASS] Part (b)(i) general solution: form equiv={max_equiv_err:.2e}, PDE={max_pde_gen:.2e}, BCs={max_bc_gen:.2e}")

    # 3. Verify part (b)(ii) Fourier coefficients of specific f(theta)
    a0 = (simpson(f_specific, 0.0, math.pi) + simpson(f_specific, math.pi, 2.0 * math.pi)) / (2.0 * math.pi)
    assert abs(a0) < 1e-14, f"a0 = {a0}"

    max_coeff_err = 0.0
    for n in range(1, 13):
        An_num = (
            simpson(lambda t: (0.5 * math.pi - t) * math.cos(n * t), 0.0, math.pi, n=4000)
            + simpson(lambda t: (-1.5 * math.pi + t) * math.cos(n * t), math.pi, 2.0 * math.pi, n=4000)
        ) / math.pi
        Bn_num = (
            simpson(lambda t: (0.5 * math.pi - t) * math.sin(n * t), 0.0, math.pi, n=4000)
            + simpson(lambda t: (-1.5 * math.pi + t) * math.sin(n * t), math.pi, 2.0 * math.pi, n=4000)
        ) / math.pi
        An_ex = (4.0 / (math.pi * n * n)) if (n % 2 == 1) else 0.0
        max_coeff_err = max(max_coeff_err, abs(An_num - An_ex), abs(Bn_num))
    assert max_coeff_err < 1e-12, f"Fourier coeff error: {max_coeff_err}"
    print(f"[PASS] Part (b)(ii) Fourier coeffs (n=1..12) max err = {max_coeff_err:.3e}")

    # Parseval check
    lhs_parseval = (
        simpson(lambda t: f_specific(t) ** 2, 0.0, math.pi)
        + simpson(lambda t: f_specific(t) ** 2, math.pi, 2.0 * math.pi)
    ) / math.pi
    rhs_parseval = math.pi ** 2 / 6.0
    assert abs(lhs_parseval - rhs_parseval) < 1e-12
    print(f"[PASS] Parseval identity: (1/pi) int f^2 = {lhs_parseval:.10f} == pi^2/6")

    # 4. Verify explicit series solution in part (b)(ii) satisfies PDE in interior and BCs on r=1, r=e^2
    def phi_spec(r, t, n_odd_terms=100):
        lr = math.log(r)
        s = lr - 1.0
        for k in range(n_odd_terms):
            n = 2 * k + 1
            # Use cosh form to avoid float overflow for large n
            # cosh(n*(lr-1))/cosh(n) = (e^{n*(lr-2)} + e^{-n*lr}) / (1 + e^{-2n})
            radial = (math.exp(n * (lr - 2.0)) + math.exp(-n * lr)) / (1.0 + math.exp(-2.0 * n))
            s += (4.0 / (math.pi * n * n)) * radial * math.cos(n * t)
        return s

    max_pde_spec = 0.0
    for r in [1.3, 2.2, 3.8, 5.9]:
        for t in [0.15, 0.9, 2.1, 3.8, 5.5]:
            max_pde_spec = max(max_pde_spec, abs(polar_laplacian_fd(lambda rr, tt: phi_spec(rr, tt, 120), r, t)))
    assert max_pde_spec < 1e-6, f"Interior PDE residual: {max_pde_spec}"

    max_bc_spec = 0.0
    for j in range(64):
        t = 2.0 * math.pi * j / 64.0
        err1 = abs(phi_spec(1.0, t, 2000) - (f_specific(t) - 1.0))
        err2 = abs(phi_spec(e2, t, 2000) - (f_specific(t) + 1.0))
        max_bc_spec = max(max_bc_spec, err1, err2)
    assert max_bc_spec < 5e-4, f"Boundary mismatch: {max_bc_spec}"
    print(f"[PASS] Part (b)(ii) explicit solution: interior PDE={max_pde_spec:.2e}, BC err (2000 odd terms)={max_bc_spec:.2e}")
    print("ALL CHECKS PASSED for IB/2023/P2/q14A.")


if __name__ == "__main__":
    main()
