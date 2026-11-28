#!/usr/bin/env python3
"""Computational verification for IB/methods/sheet-2/qx20 (2007-08 Sheet 2, Q4)."""

import math


def gauss_legendre_8(func, a, b, panels=40):
    """Composite 8-point Gauss-Legendre quadrature on [a, b]."""
    nodes = [
        -0.9602898564975363,
        -0.7966664774136267,
        -0.5255324099163290,
        -0.1834346424956498,
        0.1834346424956498,
        0.5255324099163290,
        0.7966664774136267,
        0.9602898564975363,
    ]
    weights = [
        0.1012285362903763,
        0.2223810344533745,
        0.3137066458778873,
        0.3626837833783620,
        0.3626837833783620,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763,
    ]
    h = (b - a) / panels
    total = 0.0
    for i in range(panels):
        x0 = a + i * h
        mid = x0 + 0.5 * h
        half = 0.5 * h
        s = 0.0
        for nd, wt in zip(nodes, weights):
            s += wt * func(mid + half * nd)
        total += half * s
    return total


def G(t, tau, k):
    """Green's function for y'''' + k^2 y'' = f(t) with zero ICs at t=0."""
    if t <= tau:
        return 0.0
    s = t - tau
    return s / (k ** 2) - math.sin(k * s) / (k ** 3)


def y_exact_exp(t, k):
    """Exact solution for f(t) = e^{-t}."""
    return (
        math.exp(-t) / (1.0 + k * k)
        + (t - 1.0) / (k * k)
        + (k * math.cos(k * t) - math.sin(k * t)) / (k ** 3 * (1.0 + k * k))
    )


def y_exact_exp_derivs(t, k):
    """Exact derivatives (y, y', y'', y''', y'''') for f(t) = e^{-t}."""
    c = 1.0 + k * k
    em = math.exp(-t)
    cs = math.cos(k * t)
    sn = math.sin(k * t)
    y0 = em / c + (t - 1.0) / (k * k) + cs / (k * k * c) - sn / (k ** 3 * c)
    y1 = -em / c + 1.0 / (k * k) - sn / (k * c) - cs / (k * k * c)
    y2 = em / c - cs / c + sn / (k * c)
    y3 = -em / c + k * sn / c + cs / c
    y4 = em / c + k * k * cs / c - k * sn / c
    return y0, y1, y2, y3, y4


def check_greens_function():
    """Check ODE, continuity, jump condition, and k -> 0 limit of G(t, tau)."""
    h = 2e-2
    for k in [0.7, 1.0, 2.3, -1.5]:
        for tau in [0.4, 1.25]:
            # Check G = 0 for t < tau
            for t in [0.05, 0.5 * tau, tau - 1e-6]:
                assert G(t, tau, k) == 0.0

            # Check homogeneous ODE for t > tau via O(h^4) central finite differences
            for t in [tau + 0.3, tau + 1.1, tau + 2.7]:
                gm3 = G(t - 3 * h, tau, k)
                gm2 = G(t - 2 * h, tau, k)
                gm1 = G(t - h, tau, k)
                g0 = G(t, tau, k)
                gp1 = G(t + h, tau, k)
                gp2 = G(t + 2 * h, tau, k)
                gp3 = G(t + 3 * h, tau, k)
                g4 = (-gm3 + 12 * gm2 - 39 * gm1 + 56 * g0 - 39 * gp1 + 12 * gp2 - gp3) / (6 * h ** 4)
                g2 = (-gp2 + 16 * gp1 - 30 * g0 + 16 * gm1 - gm2) / (12 * h ** 2)
                res = g4 + k * k * g2
                assert abs(res) < 1e-5, f"ODE residual too large: {res}"

            # Check right-hand limits at t = tau^+ via forward finite differences
            eps = 5e-4
            vals = [G(tau + j * eps, tau, k) for j in range(5)]
            g0_r = vals[0]
            g1_r = (-25 * vals[0] + 48 * vals[1] - 36 * vals[2] + 16 * vals[3] - 3 * vals[4]) / (12 * eps)
            g2_r = (35 * vals[0] - 104 * vals[1] + 114 * vals[2] - 56 * vals[3] + 11 * vals[4]) / (12 * eps ** 2)
            g3_r = (-5 * vals[0] + 18 * vals[1] - 24 * vals[2] + 14 * vals[3] - 3 * vals[4]) / (2 * eps ** 3)
            assert abs(g0_r) < 1e-14
            assert abs(g1_r) < 1e-10
            assert abs(g2_r) < 1e-7
            assert abs(g3_r - 1.0) < 1e-5, f"Jump mismatch: {g3_r}"

    # Check k -> 0 limit: G(t, tau) -> (t - tau)^3 / 6
    k_small = 1e-3
    t, tau = 2.1, 0.4
    limit_val = (t - tau) ** 3 / 6.0
    assert abs(G(t, tau, k_small) - limit_val) < 1e-6


def check_exp_solution():
    """Verify closed-form solution for f(t) = e^{-t} against quadrature and ODE/ICs."""
    max_quad_err = 0.0
    max_ode_err = 0.0
    max_ic_err = 0.0
    for k in [0.6, 1.0, 1.8, 3.5, -2.1]:
        # Initial conditions at t = 0
        y0, y1, y2, y3, y4 = y_exact_exp_derivs(0.0, k)
        ic_err = max(abs(y0), abs(y1), abs(y2), abs(y3))
        max_ic_err = max(max_ic_err, ic_err)
        assert ic_err < 1e-13, f"IC error for k={k}: {ic_err}"

        for t in [0.2, 0.75, 1.6, 3.1]:
            # Compare quadrature with closed form
            y_num = gauss_legendre_8(lambda tau: G(t, tau, k) * math.exp(-tau), 0.0, t, panels=30)
            y_ex = y_exact_exp(t, k)
            err = abs(y_num - y_ex)
            max_quad_err = max(max_quad_err, err)
            assert err < 1e-12, f"Quadrature mismatch for k={k}, t={t}: {err}"

            # Check analytic ODE residual
            d0, d1, d2, d3, d4 = y_exact_exp_derivs(t, k)
            ode_res = abs(d4 + k * k * d2 - math.exp(-t))
            max_ode_err = max(max_ode_err, ode_res)
            assert ode_res < 1e-13

            # Cross-check analytic derivatives against central finite differences
            h = 1e-3
            fd1 = (y_exact_exp(t + h, k) - y_exact_exp(t - h, k)) / (2 * h)
            fd2 = (y_exact_exp(t + h, k) - 2 * y_exact_exp(t, k) + y_exact_exp(t - h, k)) / (h ** 2)
            assert abs(fd1 - d1) < 1e-5
            assert abs(fd2 - d2) < 1e-4

    # Check k -> 0 limit of y(t): e^{-t} - 1 + t - t^2/2 + t^3/6
    k_small = 1e-2
    for t in [0.5, 1.5]:
        y_lim = math.exp(-t) - 1.0 + t - 0.5 * t ** 2 + (t ** 3) / 6.0
        assert abs(y_exact_exp(t, k_small) - y_lim) < 1e-3

    return max_quad_err, max_ode_err, max_ic_err


def check_general_forcing():
    """Verify convolution y(t) = int_0^t G(t,tau) f(tau) dtau solves ODE for another f(t)."""
    k = 1.7
    f = lambda s: s * s + math.cos(2.0 * s)

    def y_conv(t):
        return gauss_legendre_8(lambda tau: G(t, tau, k) * f(tau), 0.0, t, panels=40)

    h = 2e-2
    for t in [0.6, 1.4]:
        ym3 = y_conv(t - 3 * h)
        ym2 = y_conv(t - 2 * h)
        ym1 = y_conv(t - h)
        y0 = y_conv(t)
        yp1 = y_conv(t + h)
        yp2 = y_conv(t + 2 * h)
        yp3 = y_conv(t + 3 * h)
        y4 = (-ym3 + 12 * ym2 - 39 * ym1 + 56 * y0 - 39 * yp1 + 12 * yp2 - yp3) / (6 * h ** 4)
        y2 = (-yp2 + 16 * yp1 - 30 * y0 + 16 * ym1 - ym2) / (12 * h ** 2)
        res = abs(y4 + k * k * y2 - f(t))
        assert res < 1e-4, f"General forcing ODE residual: {res}"


def main():
    check_greens_function()
    max_quad_err, max_ode_err, max_ic_err = check_exp_solution()
    check_general_forcing()
    print("IB_methods_sheet-2_qx20_check.py: ALL CHECKS PASSED")
    print(f"  max quadrature vs closed-form error (f=e^-t): {max_quad_err:.3e}")
    print(f"  max analytic ODE residual:                    {max_ode_err:.3e}")
    print(f"  max IC error at t=0:                          {max_ic_err:.3e}")


if __name__ == "__main__":
    main()
