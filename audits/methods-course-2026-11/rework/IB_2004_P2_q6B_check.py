#!/usr/bin/env python3
"""Independent verification for IB/2004/P2/q6B."""

import math


def phi_exact(r, theta):
    if r < 1.0:
        return -0.25 * (r ** 2) * math.cos(2.0 * theta) - 0.125 * (r ** 4) * math.cos(4.0 * theta)
    else:
        return -0.25 * (r ** -2) * math.cos(2.0 * theta) - 0.125 * (r ** -4) * math.cos(4.0 * theta)


def dphi_dr_exact(r, theta, side):
    if side == "-":
        return -0.5 * r * math.cos(2.0 * theta) - 0.5 * (r ** 3) * math.cos(4.0 * theta)
    else:
        return 0.5 * (r ** -3) * math.cos(2.0 * theta) + 0.5 * (r ** -5) * math.cos(4.0 * theta)


def check_laplace_polar_and_cartesian():
    h = 1e-4
    max_polar_err = 0.0
    max_cart_err = 0.0
    for r in (0.15, 0.4, 0.75, 0.92, 1.08, 1.5, 2.7, 5.0):
        for k in range(13):
            theta = 2.0 * math.pi * k / 13.0 + 0.123
            # Polar 4th-order central finite differences
            pr = (
                -phi_exact(r + 2 * h, theta)
                + 8 * phi_exact(r + h, theta)
                - 8 * phi_exact(r - h, theta)
                + phi_exact(r - 2 * h, theta)
            ) / (12 * h)
            prr = (
                -phi_exact(r + 2 * h, theta)
                + 16 * phi_exact(r + h, theta)
                - 30 * phi_exact(r, theta)
                + 16 * phi_exact(r - h, theta)
                - phi_exact(r - 2 * h, theta)
            ) / (12 * h * h)
            ptt = (
                -phi_exact(r, theta + 2 * h)
                + 16 * phi_exact(r, theta + h)
                - 30 * phi_exact(r, theta)
                + 16 * phi_exact(r, theta - h)
                - phi_exact(r, theta - 2 * h)
            ) / (12 * h * h)
            lap_polar = prr + pr / r + ptt / (r * r)
            max_polar_err = max(max_polar_err, abs(lap_polar))
            assert abs(lap_polar) < 1e-7

            # Cartesian 4th-order central finite differences
            x = r * math.cos(theta)
            y = r * math.sin(theta)

            def phi_xy(xx, yy):
                return phi_exact(math.hypot(xx, yy), math.atan2(yy, xx))

            pxx = (
                -phi_xy(x + 2 * h, y)
                + 16 * phi_xy(x + h, y)
                - 30 * phi_xy(x, y)
                + 16 * phi_xy(x - h, y)
                - phi_xy(x - 2 * h, y)
            ) / (12 * h * h)
            pyy = (
                -phi_xy(x, y + 2 * h)
                + 16 * phi_xy(x, y + h)
                - 30 * phi_xy(x, y)
                + 16 * phi_xy(x, y - h)
                - phi_xy(x, y - 2 * h)
            ) / (12 * h * h)
            lap_cart = pxx + pyy
            max_cart_err = max(max_cart_err, abs(lap_cart))
            assert abs(lap_cart) < 1e-7

    print(
        f"[PASS] Laplace's equation verified in polar (max err={max_polar_err:.2e}) "
        f"& Cartesian (max err={max_cart_err:.2e})"
    )


def check_boundary_and_matching_conditions():
    # Asymptotics as r -> 0 and r -> infinity
    for k in range(36):
        theta = 2.0 * math.pi * k / 36.0
        assert abs(phi_exact(1e-6, theta)) < 1e-12
        assert abs(phi_exact(1e6, theta)) < 1e-12

    # Continuity and derivative jump at r = 1
    h = 1e-6
    max_cont_err = 0.0
    max_jump_exact_err = 0.0
    max_jump_fd_err = 0.0
    for k in range(100):
        theta = 2.0 * math.pi * k / 100.0
        phi_minus = phi_exact(1.0 - 1e-15, theta)
        phi_plus = phi_exact(1.0 + 1e-15, theta)
        max_cont_err = max(max_cont_err, abs(phi_plus - phi_minus))
        assert abs(phi_plus - phi_minus) < 1e-12

        # Exact derivative jump
        target = math.cos(2.0 * theta) + math.cos(4.0 * theta)
        jump_exact = dphi_dr_exact(1.0, theta, "+") - dphi_dr_exact(1.0, theta, "-")
        max_jump_exact_err = max(max_jump_exact_err, abs(jump_exact - target))
        assert abs(jump_exact - target) < 1e-14

        # 2nd-order one-sided finite difference jump
        phi_1 = -0.25 * math.cos(2.0 * theta) - 0.125 * math.cos(4.0 * theta)
        dp_plus = (-3.0 * phi_1 + 4.0 * phi_exact(1.0 + h, theta) - phi_exact(1.0 + 2 * h, theta)) / (2.0 * h)
        dp_minus = (3.0 * phi_1 - 4.0 * phi_exact(1.0 - h, theta) + phi_exact(1.0 - 2 * h, theta)) / (2.0 * h)
        max_jump_fd_err = max(max_jump_fd_err, abs((dp_plus - dp_minus) - target))
        assert abs((dp_plus - dp_minus) - target) < 1e-8

    print(
        f"[PASS] Asymptotics, continuity (err={max_cont_err:.2e}), "
        f"& jump condition (exact err={max_jump_exact_err:.2e}, FD err={max_jump_fd_err:.2e}) verified"
    )


def check_independent_single_layer_potential():
    # Independent check via 2D Green's function (single-layer potential on r'=1):
    # G(x, x') = (1/(2*pi)) ln |x - x'| satisfies nabla^2 G = delta(x - x'),
    # so phi(r,theta) = (1/(2*pi)) int_0^{2pi} ln(sqrt(r^2 + 1 - 2 r cos(theta - t))) (cos 2t + cos 4t) dt
    N = 256
    dt = 2.0 * math.pi / N
    max_err = 0.0
    for r in (0.2, 0.5, 0.8, 1.25, 2.0, 4.0):
        for theta in (0.0, 0.37, 1.2, math.pi / 2, 2.5, 5.1):
            s = 0.0
            for k in range(N):
                t = k * dt
                dist2 = r * r + 1.0 - 2.0 * r * math.cos(theta - t)
                f = math.cos(2.0 * t) + math.cos(4.0 * t)
                s += 0.5 * math.log(dist2) * f
            phi_sl = s * dt / (2.0 * math.pi)
            err = abs(phi_sl - phi_exact(r, theta))
            max_err = max(max_err, err)
            assert err < 1e-13
    print(f"[PASS] Independent 2D single-layer Green's integral matches exact solution (max err={max_err:.2e})")


if __name__ == "__main__":
    check_laplace_polar_and_cartesian()
    check_boundary_and_matching_conditions()
    check_independent_single_layer_potential()
    print("ALL CHECKS PASSED FOR IB/2004/P2/q6B")
