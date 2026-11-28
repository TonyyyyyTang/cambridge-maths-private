#!/usr/bin/env python3
"""Independent verification for IB/2024/P3/q5B (Laplace's equation in an annulus)."""

import math


def polar_laplacian_fd(u, r, th, hr=1e-3, hth=1e-3):
    """4th-order central finite-difference polar Laplacian."""
    u_r = (-u(r + 2 * hr, th) + 8 * u(r + hr, th) - 8 * u(r - hr, th) + u(r - 2 * hr, th)) / (12 * hr)
    u_rr = (-u(r + 2 * hr, th) + 16 * u(r + hr, th) - 30 * u(r, th) + 16 * u(r - hr, th) - u(r - 2 * hr, th)) / (12 * hr * hr)
    u_tt = (-u(r, th + 2 * hth) + 16 * u(r, th + hth) - 30 * u(r, th) + 16 * u(r, th - hth) - u(r, th - 2 * hth)) / (12 * hth * hth)
    return u_rr + u_r / r + u_tt / (r * r)


def solve_radial_ode_fd(a, b, n_mode=2, N=400):
    """Independent 2nd-order finite-difference solver for r^2 R'' + r R' - n^2 R = 0, R(a)=0, R(b)=1."""
    h = (b - a) / N
    # Interior nodes r_i = a + i*h, i = 1 .. N-1
    m = N - 1
    lower = [0.0] * m
    diag = [0.0] * m
    upper = [0.0] * m
    rhs = [0.0] * m
    for idx in range(m):
        i = idx + 1
        r = a + i * h
        c_minus = r * r / (h * h) - r / (2.0 * h)
        c_zero = -2.0 * r * r / (h * h) - float(n_mode * n_mode)
        c_plus = r * r / (h * h) + r / (2.0 * h)
        diag[idx] = c_zero
        if idx > 0:
            lower[idx] = c_minus
        if idx < m - 1:
            upper[idx] = c_plus
        else:
            rhs[idx] -= c_plus * 1.0  # R(b) = 1
    # Thomas algorithm
    for idx in range(1, m):
        w = lower[idx] / diag[idx - 1]
        diag[idx] -= w * upper[idx - 1]
        rhs[idx] -= w * rhs[idx - 1]
    R = [0.0] * m
    R[-1] = rhs[-1] / diag[-1]
    for idx in range(m - 2, -1, -1):
        R[idx] = (rhs[idx] - upper[idx] * R[idx + 1]) / diag[idx]
    return [a + (idx + 1) * h for idx in range(m)], R


def u_explicit(a, b, r, th):
    return (b ** 2 / (b ** 4 - a ** 4)) * (r ** 2 - (a ** 4) / (r ** 2)) * math.cos(2.0 * th)


def main():
    # 1. Verify general separable modes satisfy Laplace's equation
    modes = [
        lambda r, t: 1.0,
        lambda r, t: math.log(r),
    ]
    for n in range(1, 5):
        modes.extend([
            lambda r, t, n=n: (r ** n) * math.cos(n * t),
            lambda r, t, n=n: (r ** (-n)) * math.cos(n * t),
            lambda r, t, n=n: (r ** n) * math.sin(n * t),
            lambda r, t, n=n: (r ** (-n)) * math.sin(n * t),
        ])
    max_mode_err = 0.0
    for fn in modes:
        for r in [1.1, 2.3]:
            for t in [0.4, 2.2, 4.9]:
                max_mode_err = max(max_mode_err, abs(polar_laplacian_fd(fn, r, t)))
    assert max_mode_err < 1e-6, f"Mode residual: {max_mode_err}"
    print(f"[PASS] General separation modes satisfy nabla^2 u = 0, max err = {max_mode_err:.3e}")

    # 2. Verify explicit solution across multiple annuli (a, b)
    for (a, b) in [(1.0, 2.0), (0.5, 3.0), (2.0, 5.0)]:
        max_pde = 0.0
        for i in range(1, 10):
            r = a + (b - a) * i / 10.0
            for j in range(12):
                th = 2.0 * math.pi * j / 12.0
                max_pde = max(max_pde, abs(polar_laplacian_fd(lambda rr, tt: u_explicit(a, b, rr, tt), r, th)))
        assert max_pde < 1e-6, f"PDE error for ({a},{b}): {max_pde}"

        max_bc_a = 0.0
        max_bc_b = 0.0
        for j in range(64):
            th = 2.0 * math.pi * j / 64.0
            max_bc_a = max(max_bc_a, abs(u_explicit(a, b, a, th)))
            max_bc_b = max(max_bc_b, abs(u_explicit(a, b, b, th) - math.cos(2.0 * th)))
        assert max_bc_a < 1e-14 and max_bc_b < 1e-14

        # Cross-check radial factor against independent tridiagonal FD ODE solver
        rs, R_fd = solve_radial_ode_fd(a, b, n_mode=2, N=500)
        max_ode_diff = 0.0
        for r, rv in zip(rs, R_fd):
            r_exact = (b ** 2 / (b ** 4 - a ** 4)) * (r ** 2 - (a ** 4) / (r ** 2))
            max_ode_diff = max(max_ode_diff, abs(rv - r_exact))
        assert max_ode_diff < 5e-5, f"FD radial ODE mismatch for ({a},{b}): {max_ode_diff}"
        print(f"[PASS] Annulus ({a}, {b}): PDE={max_pde:.2e}, BC_a={max_bc_a:.2e}, BC_b={max_bc_b:.2e}, FD_ODE_diff={max_ode_diff:.2e}")

    print("ALL CHECKS PASSED for IB/2024/P3/q5B.")


if __name__ == "__main__":
    main()
