#!/usr/bin/env python3
"""Independent verification script for IB/2017/P3/q15A (Methods)."""

from fractions import Fraction
import math


def simpson_real(func, a, b, n):
    """Composite Simpson's rule (n even)."""
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    total = func(a) + func(b)
    for j in range(1, n):
        x = a + j * h
        total += (4 if j % 2 == 1 else 2) * func(x)
    return total * (h / 3.0)


def check_exact_algebra():
    """Verify characteristic roots and exact rational coefficients."""
    # Characteristic polynomial r^3 - r^2 - 2r = 0
    for r in [0, 2, -1]:
        assert r ** 3 - r ** 2 - 2 * r == 0

    # Green's function kernel h(u) = alpha + beta*e^(2u) + gamma*e^(-u)
    alpha = Fraction(-1, 2)
    beta = Fraction(1, 6)
    gamma = Fraction(1, 3)
    assert alpha + beta + gamma == 0          # h(0) = 0
    assert 2 * beta - gamma == 0              # h'(0) = 0
    assert 4 * beta + gamma == 1              # h''(0) = 1
    # Check L[h] = 0 for each basis function (already guaranteed by roots)

    # Homogeneous solution y_h(x) = A + B*e^(2x) + C*e^(-x)
    # with y_h(0) = 1, y_h'(0) = -1, y_h''(0) = 0
    A = Fraction(1, 2)
    B = Fraction(-1, 6)
    C = Fraction(2, 3)
    assert A + B + C == 1
    assert 2 * B - C == -1
    assert 4 * B + C == 0

    # Particular integral for e^x: L[c * e^x] = (1 - 1 - 2) c e^x = -2 c e^x = e^x
    c_p = Fraction(-1, 2)
    assert (1 - 1 - 2) * c_p == 1

    print("[PASS] Exact rational algebra verified (characteristic roots, G, y_h, y_p)")


def h_kernel(u):
    return -0.5 + (1.0 / 6.0) * math.exp(2.0 * u) + (1.0 / 3.0) * math.exp(-u)


def G_func(x, xi):
    """Green's function for L y = y''' - y'' - 2 y' with G(0)=G'(0)=G''(0)=0."""
    if xi > 0.0:
        return h_kernel(x - xi) if x > xi else 0.0
    elif xi < 0.0:
        return -h_kernel(x - xi) if x < xi else 0.0
    raise ValueError("xi must be non-zero")


def y_h(x):
    return 0.5 - (1.0 / 6.0) * math.exp(2.0 * x) + (2.0 / 3.0) * math.exp(-x)


def y_p(x):
    if x <= 1.0:
        return 0.0
    return (
        -0.5 * math.exp(x)
        + 0.5 * math.e
        + (1.0 / 6.0) * math.exp(2.0 * x - 1.0)
        - (1.0 / 6.0) * math.exp(2.0 - x)
    )


def y_full(x):
    return y_h(x) + y_p(x)


def y_piecewise(x):
    """Equivalent combined piecewise representation."""
    if x <= 1.0:
        return 0.5 - (1.0 / 6.0) * math.exp(2.0 * x) + (2.0 / 3.0) * math.exp(-x)
    return (
        0.5 * (1.0 + math.e)
        - 0.5 * math.exp(x)
        + ((math.exp(-1.0) - 1.0) / 6.0) * math.exp(2.0 * x)
        + ((4.0 - math.e ** 2) / 6.0) * math.exp(-x)
    )


def check_greens_function():
    """Verify BCs, ODE away from x=xi, continuity of G and G', and unit jump of G''."""
    eps = 1e-4
    for xi in [0.4, 1.3, 2.5, -0.7]:
        # Homogeneous BCs at x = 0
        g0 = G_func(0.0, xi)
        gp0 = (G_func(eps, xi) - G_func(-eps, xi)) / (2.0 * eps)
        gpp0 = (G_func(eps, xi) - 2.0 * g0 + G_func(-eps, xi)) / (eps ** 2)
        assert abs(g0) < 1e-14 and abs(gp0) < 1e-14 and abs(gpp0) < 1e-14

        # Jump conditions at x = xi using one-sided finite differences
        g_left = G_func(xi - 1e-12, xi)
        g_right = G_func(xi + 1e-12, xi)
        assert abs(g_right - g_left) < 1e-11

        gp_left = (G_func(xi - 1e-7, xi) - G_func(xi - 2e-7, xi)) / 1e-7
        gp_right = (G_func(xi + 2e-7, xi) - G_func(xi + 1e-7, xi)) / 1e-7
        assert abs(gp_right - gp_left) < 1e-6

        h = 1e-4
        gpp_left = (
            G_func(xi - h, xi) - 2.0 * G_func(xi - 2 * h, xi) + G_func(xi - 3 * h, xi)
        ) / (h ** 2)
        gpp_right = (
            G_func(xi + 3 * h, xi) - 2.0 * G_func(xi + 2 * h, xi) + G_func(xi + h, xi)
        ) / (h ** 2)
        jump = gpp_right - gpp_left
        assert abs(jump - 1.0) < 2e-3, f"Jump mismatch at xi={xi}: jump={jump}"

        # ODE L_x G = 0 away from x = xi (use h=1e-3 to balance O(h^2) truncation vs eps/h^3 roundoff)
        h_ode = 1e-3
        for x in [xi - 0.5, xi + 0.5]:
            d1 = (G_func(x + h_ode, xi) - G_func(x - h_ode, xi)) / (2.0 * h_ode)
            d2 = (
                G_func(x + h_ode, xi) - 2.0 * G_func(x, xi) + G_func(x - h_ode, xi)
            ) / (h_ode ** 2)
            d3 = (
                G_func(x + 2 * h_ode, xi)
                - 2.0 * G_func(x + h_ode, xi)
                + 2.0 * G_func(x - h_ode, xi)
                - G_func(x - 2 * h_ode, xi)
            ) / (2.0 * h_ode ** 3)
            res = d3 - d2 - 2.0 * d1
            assert abs(res) < 1e-5, f"ODE residual for G at x={x}, xi={xi}: {res}"

    print("[PASS] Green's function BCs, ODE, continuity, and unit G'' jump verified")


def check_particular_and_full_solution():
    """Verify y_p via Green's convolution, full ODE, BCs, smoothness at x=1, and RK4."""
    # 1. Green's integral y_p(x) = int_1^x G(x;xi) e^xi dxi
    max_quad_err = 0.0
    for x in [1.1, 1.5, 2.0, 2.8, 3.5]:
        quad_val = simpson_real(
            lambda xi, x=x: h_kernel(x - xi) * math.exp(xi), 1.0, x, 1600
        )
        err = abs(quad_val - y_p(x))
        max_quad_err = max(max_quad_err, err)
        assert err < 1e-11, f"Green integral mismatch at x={x}: err={err}"

    # 2. Equivalence of y_full and y_piecewise
    for x in [-1.5, -0.5, 0.0, 0.7, 1.0, 1.0001, 1.7, 2.6]:
        assert abs(y_full(x) - y_piecewise(x)) < 1e-13

    # 3. Initial conditions at x = 0
    h = 1e-4
    y0 = y_full(0.0)
    yp0 = (y_full(h) - y_full(-h)) / (2.0 * h)
    ypp0 = (y_full(h) - 2.0 * y_full(0.0) + y_full(-h)) / (h ** 2)
    assert abs(y0 - 1.0) < 1e-14
    assert abs(yp0 - (-1.0)) < 1e-7
    assert abs(ypp0 - 0.0) < 1e-7

    # 4. Continuity of y, y', y'' at x = 1 and jump of y''' equal to e
    y1_l = y_full(1.0 - 1e-12)
    y1_r = y_full(1.0 + 1e-12)
    assert abs(y1_r - y1_l) < 1e-10

    # Analytic derivatives of y_p at x = 1+
    yp_p_1 = math.e * (-0.5 + 1.0 / 3.0 + 1.0 / 6.0)
    yp_pp_1 = math.e * (-0.5 + 2.0 / 3.0 - 1.0 / 6.0)
    yp_ppp_1 = math.e * (-0.5 + 4.0 / 3.0 + 1.0 / 6.0)
    assert abs(yp_p_1) < 1e-15
    assert abs(yp_pp_1) < 1e-15
    assert abs(yp_ppp_1 - math.e) < 1e-14

    # 5. Finite-difference ODE check for x < 1 and x > 1
    h_ode = 1e-3
    max_ode_err = 0.0
    for x in [-1.2, -0.3, 0.4, 0.85, 1.15, 1.8, 2.7]:
        d1 = (y_full(x + h_ode) - y_full(x - h_ode)) / (2.0 * h_ode)
        d2 = (y_full(x + h_ode) - 2.0 * y_full(x) + y_full(x - h_ode)) / (h_ode ** 2)
        d3 = (
            y_full(x + 2 * h_ode)
            - 2.0 * y_full(x + h_ode)
            + 2.0 * y_full(x - h_ode)
            - y_full(x - 2 * h_ode)
        ) / (2.0 * h_ode ** 3)
        rhs = math.exp(x) if x > 1.0 else 0.0
        err = abs((d3 - d2 - 2.0 * d1) - rhs)
        max_ode_err = max(max_ode_err, err)
        assert err < 1e-4, f"ODE residual at x={x}: err={err}"

    # 6. Independent RK4 integration from x=0 to x=3 (split at x=1)
    def rk4_segment(x0, x1, state, steps, forcing):
        dt = (x1 - x0) / steps
        x = x0
        y0_, y1_, y2_ = state

        def f(t, s):
            return (s[1], s[2], s[2] + 2.0 * s[1] + forcing(t))

        for _ in range(steps):
            k1 = f(x, (y0_, y1_, y2_))
            k2 = f(
                x + 0.5 * dt,
                (y0_ + 0.5 * dt * k1[0], y1_ + 0.5 * dt * k1[1], y2_ + 0.5 * dt * k1[2]),
            )
            k3 = f(
                x + 0.5 * dt,
                (y0_ + 0.5 * dt * k2[0], y1_ + 0.5 * dt * k2[1], y2_ + 0.5 * dt * k2[2]),
            )
            k4 = f(
                x + dt,
                (y0_ + dt * k3[0], y1_ + dt * k3[1], y2_ + dt * k3[2]),
            )
            y0_ += (dt / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            y1_ += (dt / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            y2_ += (dt / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
            x += dt
        return (y0_, y1_, y2_)

    state_1 = rk4_segment(0.0, 1.0, (1.0, -1.0, 0.0), 1000, lambda t: 0.0)
    assert abs(state_1[0] - y_full(1.0)) < 1e-10
    state_3 = rk4_segment(1.0, 3.0, state_1, 2000, lambda t: math.exp(t))
    rk4_err = abs(state_3[0] - y_full(3.0))
    assert rk4_err < 1e-8, f"RK4 mismatch at x=3: {rk4_err}"

    print(
        f"[PASS] Full solution verified: quad err = {max_quad_err:.3e}, "
        f"ODE residual = {max_ode_err:.3e}, RK4 err at x=3 = {rk4_err:.3e}"
    )


if __name__ == "__main__":
    check_exact_algebra()
    check_greens_function()
    check_particular_and_full_solution()
    print("ALL CHECKS PASSED FOR IB_2017_P3_q15A")
