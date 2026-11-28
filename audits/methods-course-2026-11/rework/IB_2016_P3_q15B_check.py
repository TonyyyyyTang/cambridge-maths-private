#!/usr/bin/env python3
"""Independent verification script for IB/2016/P3/q15B (Methods)."""

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


def check_part_a():
    """Verify Fourier transform of exp(-a^2 x^2) = (sqrt(pi)/a) exp(-k^2/(4 a^2))."""
    max_err = 0.0
    for a in [0.5, 1.0, 1.7, 2.5]:
        L = 10.0 / a  # exp(-100) ~ 3.7e-44
        for k in [0.0, 0.7, -1.3, 3.2, -5.0]:
            num = simpson_complex(
                lambda x, a=a, k=k: math.exp(-(a * x) ** 2) * cmath.exp(-1j * k * x),
                -L,
                L,
                4000,
            )
            exact = (math.sqrt(math.pi) / a) * math.exp(-(k ** 2) / (4.0 * a ** 2))
            err = abs(num - exact)
            max_err = max(max_err, err)
            assert err < 1e-11, f"Part (a) mismatch for a={a}, k={k}: err={err}"
    print(f"[PASS] Part (a): Gaussian FT verified (max error = {max_err:.3e})")


def check_part_b():
    """Verify diffusion equation heat-kernel solution and Property P."""
    D = 1.3
    alpha = 1.0

    # Test 1: Box initial data g(y) = 1 on [-1, 1], 0 outside
    def theta_box_quad(x, t):
        pref = 1.0 / math.sqrt(4.0 * math.pi * D * t)
        return pref * simpson_real(
            lambda y: math.exp(-((x - y) ** 2) / (4.0 * D * t)),
            -alpha,
            alpha,
            1200,
        )

    def theta_box_erf(x, t):
        s = math.sqrt(4.0 * D * t)
        return 0.5 * (math.erf((x + alpha) / s) - math.erf((x - alpha) / s))

    max_erf_err = 0.0
    for x in [-3.0, -0.5, 0.0, 0.8, 2.5]:
        for t in [0.05, 0.2, 1.0]:
            err = abs(theta_box_quad(x, t) - theta_box_erf(x, t))
            max_erf_err = max(max_erf_err, err)
            assert err < 1e-11, f"Box heat kernel mismatch at x={x}, t={t}: err={err}"

    # Verify PDE theta_t = D theta_xx via central differences on exact erf solution
    max_pde_err = 0.0
    dt = 1e-5
    dx = 1e-4
    for x in [-2.3, -0.4, 0.7, 1.9]:
        for t in [0.15, 0.6]:
            u_t = (theta_box_erf(x, t + dt) - theta_box_erf(x, t - dt)) / (2.0 * dt)
            u_xx = (
                theta_box_erf(x + dx, t)
                - 2.0 * theta_box_erf(x, t)
                + theta_box_erf(x - dx, t)
            ) / (dx ** 2)
            err = abs(u_t - D * u_xx)
            max_pde_err = max(max_pde_err, err)
            assert err < 1e-6, f"Diffusion PDE residual too large: {err}"

    # Verify initial condition recovery as t -> 0+
    t_small = 1e-6
    assert abs(theta_box_erf(0.3, t_small) - 1.0) < 1e-12
    assert abs(theta_box_erf(2.0, t_small) - 0.0) < 1e-12

    # Test 2: Smooth compactly supported positive bump g(y) = (1 - y^2)^2 on [-1, 1]
    def g_bump(y):
        return (1.0 - y * y) ** 2 if abs(y) < alpha else 0.0

    integral_g = simpson_real(g_bump, -alpha, alpha, 2000)
    assert abs(integral_g - 16.0 / 15.0) < 1e-12

    # Verify Property P and strict lower bound for large beta, small epsilon
    for beta, eps in [(3.0, 0.1), (5.0, 0.2), (-4.0, 0.15)]:
        pref = 1.0 / math.sqrt(4.0 * math.pi * D * eps)
        val = pref * simpson_real(
            lambda y: g_bump(y) * math.exp(-((beta - y) ** 2) / (4.0 * D * eps)),
            -alpha,
            alpha,
            600,
        )
        lower_bound = (
            pref
            * math.exp(-((abs(beta) + alpha) ** 2) / (4.0 * D * eps))
            * integral_g
        )
        assert val > 0.0, f"Property P failed at beta={beta}, eps={eps}"
        assert val >= lower_bound * (1.0 - 1e-12), (
            f"Lower bound violated: val={val}, lb={lower_bound}"
        )

    print(
        f"[PASS] Part (b): Heat kernel vs erf max err = {max_erf_err:.3e}, "
        f"PDE residual = {max_pde_err:.3e}, Property P + lower bound verified"
    )


def check_part_c():
    """Verify wave equation d'Alembert solution and failure of Property P."""
    alpha = 1.0

    # C^2 compactly supported initial data phi(x) = (1 - x^2)^3 on [-1, 1]
    def phi(x):
        return (1.0 - x * x) ** 3 if abs(x) < alpha else 0.0

    for c in [1.5, -0.8, 0.0]:
        def u(x, t, c=c):
            return 0.5 * (phi(x - c * t) + phi(x + c * t))

        # Check initial data u(x,0) = phi(x) and u_t(x,0) = 0
        dt = 1e-5
        for x in [-1.5, -0.6, 0.0, 0.4, 1.3]:
            assert abs(u(x, 0.0) - phi(x)) < 1e-15
            u_t0 = (u(x, dt) - u(x, -dt)) / (2.0 * dt)
            assert abs(u_t0) < 1e-12

        # Check wave PDE u_tt = c^2 u_xx away from support boundaries x +- c*t = +-1
        dx = 1e-4
        max_wave_err = 0.0
        for x, t in [(0.2, 0.3), (-0.5, 0.2), (1.1, 0.4)]:
            # Skip points too close to x +- c*t = +-1 where 3rd derivative jumps
            if min(
                abs(abs(x - c * t) - alpha), abs(abs(x + c * t) - alpha)
            ) < 0.05:
                continue
            u_tt = (u(x, t + dt) - 2.0 * u(x, t) + u(x, t - dt)) / (dt ** 2)
            u_xx = (u(x + dx, t) - 2.0 * u(x, t) + u(x - dx, t)) / (dx ** 2)
            err = abs(u_tt - (c ** 2) * u_xx)
            max_wave_err = max(max_wave_err, err)
            assert err < 1e-5, f"Wave PDE residual too large for c={c}: {err}"

        # Check failure of Property P: whenever |beta| > alpha + |c| * eps, u(beta, eps) == 0
        for eps in [1e-4, 0.05, 0.3, 1.0]:
            for margin in [1e-6, 0.1, 5.0]:
                for sign in [-1.0, 1.0]:
                    beta = sign * (alpha + abs(c) * eps + margin)
                    assert u(beta, eps) == 0.0, (
                        f"Expected u({beta}, {eps}) == 0 for c={c}, got {u(beta, eps)}"
                    )

        # Also check that inside the domain of influence u can be positive
        assert u(abs(c) * 0.5, 0.5) > 0.0

    print("[PASS] Part (c): d'Alembert PDE/ICs and failure of Property P verified")


if __name__ == "__main__":
    check_part_a()
    check_part_b()
    check_part_c()
    print("ALL CHECKS PASSED FOR IB_2016_P3_q15B")
