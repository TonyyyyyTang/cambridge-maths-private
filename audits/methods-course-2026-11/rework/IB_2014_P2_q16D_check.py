#!/usr/bin/env python3
"""Computational verification for IB/2014/P2/q16D (1D & 2D diffusion Green's functions and exponential integral)."""

import math
import cmath

GL_NODES = [
    -0.96028985649753623168,
    -0.79666647741362673959,
    -0.52553240991632898582,
    -0.18343464249564980494,
    0.18343464249564980494,
    0.52553240991632898582,
    0.79666647741362673959,
    0.96028985649753623168,
]
GL_WEIGHTS = [
    0.10122853629037625915,
    0.22238103445337447054,
    0.31370664587788728734,
    0.36268378337836198297,
    0.36268378337836198297,
    0.31370664587788728734,
    0.22238103445337447054,
    0.10122853629037625915,
]

EULER_GAMMA = 0.57721566490153286060651209


def quad(func, a, b, panels=40):
    h = (b - a) / panels
    total = 0.0
    for j in range(panels):
        x0 = a + j * h
        mid = x0 + 0.5 * h
        half = 0.5 * h
        for t, w in zip(GL_NODES, GL_WEIGHTS):
            total += w * func(mid + half * t)
    return 0.5 * h * total


def Ei_neg(eta):
    """Compute Ei(-eta) = -int_eta^infty e^{-lambda}/lambda dlambda for eta > 0."""
    assert eta > 0.0
    if eta < 10.0:
        # Exact power series: Ei(-eta) = gamma + ln(eta) + sum_{n=1}^infty (-eta)^n / (n * n!)
        s = EULER_GAMMA + math.log(eta)
        term = 1.0
        for n in range(1, 120):
            term *= (-eta) / n
            add = term / n
            s += add
            if abs(add) < 1e-16 * abs(s):
                break
        return s
    # For large eta, direct quadrature with substitution lambda = eta + u
    return -quad(lambda u: math.exp(-(eta + u)) / (eta + u), 0.0, 50.0, panels=100)


def G1D(t, tp, x, xp, D):
    if t <= tp:
        return 0.0
    dt = t - tp
    return math.exp(-((x - xp) ** 2) / (4.0 * D * dt)) / math.sqrt(4.0 * math.pi * D * dt)


def G2D(t, tp, x, xp, y, yp, D):
    if t <= tp:
        return 0.0
    dt = t - tp
    r2 = (x - xp) ** 2 + (y - yp) ** 2
    return math.exp(-r2 / (4.0 * D * dt)) / (4.0 * math.pi * D * dt)


def psi_exact_quad(t, x, y, t0, p0, D):
    """Direct time-convolution integral for t > t0."""
    assert t > t0
    r2 = x * x + y * y
    return (p0 / (4.0 * math.pi * D)) * quad(
        lambda tp: math.exp(-r2 / (4.0 * D * (t - tp))) / (t - tp),
        0.0,
        t0,
        panels=80,
    )


def psi_ei(t, x, y, t0, p0, D):
    """Closed-form solution in terms of Ei(-eta) for t > t0 and (x,y) != (0,0)."""
    r2 = x * x + y * y
    eta1 = r2 / (4.0 * D * t)
    eta2 = r2 / (4.0 * D * (t - t0))
    return (p0 / (4.0 * math.pi * D)) * (Ei_neg(eta2) - Ei_neg(eta1))


def main():
    # --- Part (a): 1D Green's function Fourier inversion and PDE checks ---
    D = 1.7
    t, tp, x, xp = 2.3, 0.8, 1.4, -0.5
    dt = t - tp

    # Numerical Fourier inversion of G_tilde(t, t', k, x')
    K_max = 12.0 / math.sqrt(D * dt)
    inv_val = (1.0 / (2.0 * math.pi)) * quad(
        lambda k: (cmath.exp(-1j * k * xp) * math.exp(-D * k * k * dt) * cmath.exp(1j * k * x)).real,
        -K_max,
        K_max,
        panels=120,
    )
    g1d_val = G1D(t, tp, x, xp, D)
    err_inv = abs(inv_val - g1d_val)
    assert err_inv < 1e-13, f"1D Fourier inversion mismatch: {err_inv}"

    # Check 1D diffusion PDE by central finite differences
    h_t, h_x = 1e-5, 1e-4
    gt = (G1D(t + h_t, tp, x, xp, D) - G1D(t - h_t, tp, x, xp, D)) / (2.0 * h_t)
    gxx = (G1D(t, tp, x + h_x, xp, D) - 2.0 * g1d_val + G1D(t, tp, x - h_x, xp, D)) / (h_x ** 2)
    pde_1d_res = abs(gt - D * gxx)
    assert pde_1d_res < 1e-7, f"1D PDE residual: {pde_1d_res}"

    # Check spatial normalisation int G1D dx = 1
    L = 12.0 * math.sqrt(D * dt)
    norm_1d = quad(lambda xx: G1D(t, tp, xx, xp, D), xp - L, xp + L, panels=100)
    assert abs(norm_1d - 1.0) < 1e-13

    # --- Part (b): 2D Green's function PDE, Ei(-eta) formula, and small-eta approximation ---
    y, yp = -0.9, 0.3
    g2d_val = G2D(t, tp, x, xp, y, yp, D)
    g2t = (G2D(t + h_t, tp, x, xp, y, yp, D) - G2D(t - h_t, tp, x, xp, y, yp, D)) / (2.0 * h_t)
    g2xx = (G2D(t, tp, x + h_x, xp, y, yp, D) - 2.0 * g2d_val + G2D(t, tp, x - h_x, xp, y, yp, D)) / (h_x ** 2)
    g2yy = (G2D(t, tp, x, xp, y + h_x, yp, D) - 2.0 * g2d_val + G2D(t, tp, x, xp, y - h_x, yp, D)) / (h_x ** 2)
    pde_2d_res = abs(g2t - D * (g2xx + g2yy))
    assert pde_2d_res < 1e-7, f"2D PDE residual: {pde_2d_res}"

    # Compare Ei formula against direct time integral across several parameter sets
    p0 = 3.5
    t0 = 1.2
    max_ei_err = 0.0
    for (tt, xx, yy) in [(1.5, 0.7, -0.4), (2.5, 2.0, 1.5), (5.0, 0.2, 0.3), (10.0, 5.0, -3.0)]:
        v_quad = psi_exact_quad(tt, xx, yy, t0, p0, D)
        v_ei = psi_ei(tt, xx, yy, t0, p0, D)
        max_ei_err = max(max_ei_err, abs(v_quad - v_ei))
    assert max_ei_err < 1e-13, f"Psi Ei vs quadrature error: {max_ei_err}"

    # Verify Psi satisfies homogeneous 2D diffusion PDE for t > t0
    tt, xx, yy = 2.4, 0.8, -0.6
    ps = psi_ei(tt, xx, yy, t0, p0, D)
    pst = (psi_ei(tt + h_t, xx, yy, t0, p0, D) - psi_ei(tt - h_t, xx, yy, t0, p0, D)) / (2.0 * h_t)
    psxx = (psi_ei(tt, xx + h_x, yy, t0, p0, D) - 2.0 * ps + psi_ei(tt, xx - h_x, yy, t0, p0, D)) / (h_x ** 2)
    psyy = (psi_ei(tt, xx, yy + h_x, t0, p0, D) - 2.0 * ps + psi_ei(tt, xx, yy - h_x, t0, p0, D)) / (h_x ** 2)
    psi_pde_res = abs(pst - D * (psxx + psyy))
    assert psi_pde_res < 1e-7, f"Psi PDE residual: {psi_pde_res}"

    # Verify small-eta approximation: Psi -> (p0 / (4 pi D)) * ln(t / (t - t0))
    # and check next-order error term -(p0 / (4 pi D)) * (eta2 - eta1)
    tt = 3.0
    approx_val = (p0 / (4.0 * math.pi * D)) * math.log(tt / (tt - t0))
    # Also check exact value at origin x=y=0 equals approx_val!
    origin_val = psi_exact_quad(tt, 0.0, 0.0, t0, p0, D)
    assert abs(origin_val - approx_val) < 1e-13

    for r in [1e-1, 1e-2, 1e-3]:
        exact_v = psi_ei(tt, r, 0.0, t0, p0, D)
        eta1 = (r * r) / (4.0 * D * tt)
        eta2 = (r * r) / (4.0 * D * (tt - t0))
        lead_err = exact_v - approx_val
        predicted_next = -(p0 / (4.0 * math.pi * D)) * (eta2 - eta1)
        rel_next = abs(lead_err - predicted_next) / abs(predicted_next)
        assert rel_next < 2.0 * eta2, f"Asymptotic check failed at r={r}: rel={rel_next}"

    # Verify initial condition Psi(t, x, y) -> 0 as t -> 0+ for (x, y) != (0, 0)
    # For 0 < t < t0, Psi = -(p0/(4 pi D)) Ei(-r^2/(4 D t))
    ic_val = -(p0 / (4.0 * math.pi * D)) * Ei_neg(1.0 / (4.0 * D * 1e-3))
    assert abs(ic_val) < 1e-30

    print("IB_2014_P2_q16D verification PASSED:")
    print(f"  1D Fourier inversion error = {err_inv:.3e}, norm = {norm_1d:.16f}")
    print(f"  1D PDE residual = {pde_1d_res:.3e}, 2D PDE residual = {pde_2d_res:.3e}")
    print(f"  max |Psi_Ei - Psi_quad| = {max_ei_err:.3e}")
    print(f"  Psi(t>t0) 2D PDE residual = {psi_pde_res:.3e}")
    print(f"  Psi(t,0,0) = {origin_val:.16f} == (p0/(4piD))*ln(t/(t-t0)) = {approx_val:.16f}")
    print(f"  Initial condition check at t=1e-3, r=1: Psi = {ic_val:.3e}")


if __name__ == "__main__":
    main()
