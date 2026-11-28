#!/usr/bin/env python3
"""Verification script for IB/methods/sheet-3/qx07 (Brachistochrone)."""

import math


def solve_theta_for_ratio(ratio):
    """Solve (theta - sin(theta)) / (1 - cos(theta)) = ratio for theta in (0, 2*pi)."""
    lo, hi = 1e-12, 2.0 * math.pi - 1e-12
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        val = (mid - math.sin(mid)) / (1.0 - math.cos(mid))
        if val < ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def S_calibration(x, y, g=9.81):
    """Hamilton-Jacobi calibration function S(x, y) for y > 0."""
    if x <= 0.0:
        return math.sqrt(2.0 * y / g)
    theta = solve_theta_for_ratio(x / y)
    a = y / (1.0 - math.cos(theta))
    return theta * math.sqrt(a / g)


def composite_simpson(f, a, b, n):
    assert n % 2 == 0
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * f(a + i * h)
    return s * h / 3.0


def main():
    g = 9.81
    l = 3.7  # arbitrary positive distance
    a = l / (2.0 * math.pi)
    T_exact = math.sqrt(2.0 * math.pi * l / g)

    # 1. Verify Beltrami first integral and full 2nd-order Euler-Lagrange ODE on cycloid
    # F(y, y') = sqrt((1 + y'^2)/(2 g y))
    # Full EL ODE: 2 y y'' + (1 + y'^2) = 0
    max_beltrami_err = 0.0
    max_el_err = 0.0
    for k in range(1, 100):
        th = 2.0 * math.pi * k / 100.0
        y = a * (1.0 - math.cos(th))
        dx_dth = a * (1.0 - math.cos(th))
        dy_dth = a * math.sin(th)
        yp = dy_dth / dx_dth  # cot(th/2)
        assert abs(yp - 1.0 / math.tan(0.5 * th)) < 1e-12

        # y'' = (d(yp)/dth) / (dx/dth) = (-1/(2 sin^2(th/2))) / (2a sin^2(th/2))
        ypp = (-0.5 / (math.sin(0.5 * th) ** 2)) / dx_dth

        beltrami_err = abs(y * (1.0 + yp**2) - 2.0 * a)
        el_err = abs(2.0 * y * ypp + (1.0 + yp**2))
        max_beltrami_err = max(max_beltrami_err, beltrami_err)
        max_el_err = max(max_el_err, el_err)
    assert max_beltrami_err < 1e-11
    assert max_el_err < 1e-10
    print(
        f"[PASS] Cycloid satisfies Beltrami (max err={max_beltrami_err:.2e}) "
        f"and full EL ODE 2yy''+1+y'^2=0 (max err={max_el_err:.2e})"
    )

    # 2. Verify vertical tangents at A and B and endpoints
    x_A, y_A = a * (0.0 - math.sin(0.0)), a * (1.0 - math.cos(0.0))
    x_B, y_B = a * (2.0 * math.pi - math.sin(2.0 * math.pi)), a * (1.0 - math.cos(2.0 * math.pi))
    assert abs(x_A) < 1e-15 and abs(y_A) < 1e-15
    assert abs(x_B - l) < 1e-14 and abs(y_B) < 1e-14
    eps = 1e-6
    assert 1.0 / math.tan(0.5 * eps) > 1e5
    assert 1.0 / math.tan(0.5 * (2.0 * math.pi - eps)) < -1e5
    print("[PASS] Endpoints (0,0), (l,0) and vertical tangents at A, B verified")

    # 3. Direct numerical quadrature of travel time along cycloid
    # (a) Parametric quadrature avoiding exact 0 and 2pi endpoints
    def integrand_th(th):
        one_minus_cos = 2.0 * (math.sin(0.5 * th) ** 2)
        dx = a * one_minus_cos
        dy = a * math.sin(th)
        y = a * one_minus_cos
        return math.sqrt((dx * dx + dy * dy) / (2.0 * g * y))

    T_num_th = composite_simpson(integrand_th, 1e-10, 2.0 * math.pi - 1e-10, 2000)
    assert abs(T_num_th - T_exact) < 1e-9

    # (b) Direct x-integral via substitution x = l * sin^2(u) to resolve integrable x^(-1/3) singularity
    # Note near x=0: y ~ (9/2 a)^(1/3) x^(2/3), sqrt(1+y'^2)/sqrt(y) ~ x^(-2/3), integrable!
    # Use x = l * (3u^2 - 2u^3) or x = l * u^3 / (u^3 + (1-u)^3) for smooth endpoint clustering.
    def theta_minus_sin(th):
        if th < 1e-3:
            t2 = th * th
            return (th * t2) * (1.0 / 6.0 - t2 / 120.0 + (t2 * t2) / 5040.0)
        return th - math.sin(th)

    def eval_cycloid_y_yp_at_x(x_val):
        # Solve a*(th - sin(th)) = x_val for x_val in (0, l/2]
        lo, hi = 0.0, math.pi
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if a * theta_minus_sin(mid) < x_val:
                lo = mid
            else:
                hi = mid
        th = 0.5 * (lo + hi)
        return 2.0 * a * (math.sin(0.5 * th) ** 2), 1.0 / math.tan(0.5 * th)

    def x_integrand_u(u):
        # By symmetry T = 2 * int_0^{l/2} F(y, y') dx; set x = (l/2) * u^3 so dx = (3l/2) u^2 du
        # Near x=0, F ~ x^(-2/3) ~ u^(-2), so F * dx/du is smooth and finite as u -> 0+!
        x_val = 0.5 * l * (u**3)
        dx_du = 1.5 * l * (u**2)
        y_val, yp_val = eval_cycloid_y_yp_at_x(x_val)
        return 2.0 * math.sqrt((1.0 + yp_val**2) / (2.0 * g * y_val)) * dx_du

    T_num_x = composite_simpson(x_integrand_u, 1e-12, 1.0, 400)
    assert abs(T_num_x - T_exact) < 1e-9, f"T_num_x={T_num_x}, T_exact={T_exact}"
    print(f"[PASS] Travel time quadrature: exact={T_exact:.10f}, param={T_num_th:.10f}, x-int={T_num_x:.10f}")

    # 4. Verify Hamilton-Jacobi calibration |grad S|^2 = 1/(2gy) and exact gradient formula
    h = 1e-6
    max_hj_err = 0.0
    for x_val in [0.1, 0.5, 1.5, 3.0]:
        for y_val in [0.1, 0.4, 1.2, 2.5]:
            Sx = (S_calibration(x_val + h, y_val, g) - S_calibration(x_val - h, y_val, g)) / (2.0 * h)
            Sy = (S_calibration(x_val, y_val + h, g) - S_calibration(x_val, y_val - h, g)) / (2.0 * h)
            th = solve_theta_for_ratio(x_val / y_val)
            a_loc = y_val / (1.0 - math.cos(th))
            Sx_ex = 1.0 / (2.0 * math.sqrt(g * a_loc))
            Sy_ex = (1.0 / math.tan(0.5 * th)) / (2.0 * math.sqrt(g * a_loc))
            assert abs(Sx - Sx_ex) < 1e-7
            assert abs(Sy - Sy_ex) < 1e-7
            hj_err = abs((Sx * Sx + Sy * Sy) - 1.0 / (2.0 * g * y_val))
            max_hj_err = max(max_hj_err, hj_err)
    assert max_hj_err < 1e-7
    print(f"[PASS] Hamilton-Jacobi calibration S(x,y)=theta*sqrt(a/g) verified (max |grad S|^2 err={max_hj_err:.2e})")

    # 5. Compare against competitor tunnels and perturbations of the cycloid
    # (a) Perturbations in theta-parametric form: y(th) = a(1-cos th) + eps*sin^2(th/2)*sin(k*th/2)
    for k in [1, 2, 3]:
        for eps_p in [-0.15, -0.05, 0.05, 0.15]:
            def pert_int(th):
                s_half = math.sin(0.5 * th)
                one_minus_cos = 2.0 * s_half * s_half
                dx = a * one_minus_cos
                y = a * one_minus_cos + eps_p * a * (s_half**2) * math.sin(0.5 * k * th)
                dy = a * math.sin(th) + eps_p * a * (
                    s_half * math.cos(0.5 * th) * math.sin(0.5 * k * th)
                    + 0.5 * k * (s_half**2) * math.cos(0.5 * k * th)
                )
                return math.sqrt((dx * dx + dy * dy) / (2.0 * g * y))

            T_pert = composite_simpson(pert_int, 1e-8, 2.0 * math.pi - 1e-8, 1000)
            assert T_pert > T_exact + 1e-6, f"Perturbation k={k}, eps={eps_p} beat cycloid!"

    # (b) Semicircular tunnel: x = (l/2)(1 - cos phi), y = (l/2) sin phi
    T_semi = composite_simpson(
        lambda ph: (0.5 * l) / math.sqrt(g * l * math.sin(ph)),
        1e-8,
        math.pi - 1e-8,
        2000,
    )
    # Better accuracy for endpoint singularity phi^(-1/2) via phi = pi * sin^2(u)
    T_semi = composite_simpson(
        lambda u: (0.5 * l) * (math.pi * 2.0 * math.sin(u) * math.cos(u)) / math.sqrt(g * l * math.sin(math.pi * math.sin(u) ** 2)),
        1e-8,
        0.5 * math.pi - 1e-8,
        1000,
    )
    # (c) Optimal V-tunnel (depth l/2)
    T_V = 2.0 * math.sqrt(2.0 * l / g)
    # (d) Parabolic tunnel y = 4H x(l-x)/l^2 at H = 2a = l/pi (same max depth as cycloid)
    H = l / math.pi
    T_parab = composite_simpson(
        lambda u: math.sqrt(1.0 + (4.0 * H * (1.0 - 2.0 * math.sin(u) ** 2) / l) ** 2)
        / math.sqrt(2.0 * g * 4.0 * H * (math.sin(u) ** 2) * (math.cos(u) ** 2))
        * (2.0 * l * math.sin(u) * math.cos(u)),
        1e-8,
        0.5 * math.pi - 1e-8,
        1000,
    )
    assert T_exact < T_semi < T_V
    assert T_exact < T_parab
    print(
        f"[PASS] Strictly minimal vs competitors (times in units of sqrt(l/g)): "
        f"cycloid={T_exact/math.sqrt(l/g):.5f}, semicircle={T_semi/math.sqrt(l/g):.5f}, "
        f"parabola={T_parab/math.sqrt(l/g):.5f}, V-tunnel={T_V/math.sqrt(l/g):.5f}"
    )
    print("ALL CHECKS PASSED for IB_methods_sheet-3_qx07.")


if __name__ == "__main__":
    main()
