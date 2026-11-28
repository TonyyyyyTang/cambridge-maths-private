#!/usr/bin/env python3
"""Independent computational verification for IB/2005/P2/q5E."""

import cmath
import math


_GL8_X = [
    -0.96028985649753623168,
    -0.79666647741362673959,
    -0.52553240991632898582,
    -0.18343464249564980494,
    0.18343464249564980494,
    0.52553240991632898582,
    0.79666647741362673959,
    0.96028985649753623168,
]
_GL8_W = [
    0.10122853629037625915,
    0.22238103445337447054,
    0.31370664587788728734,
    0.36268378337836198297,
    0.36268378337836198297,
    0.31370664587788728734,
    0.22238103445337447054,
    0.10122853629037625915,
]


def quad(f, a, b, panels=40):
    if b == a:
        return 0.0
    h = (b - a) / panels
    total = 0.0
    for j in range(panels):
        mid = a + (j + 0.5) * h
        half = 0.5 * h
        for xi, wi in zip(_GL8_X, _GL8_W):
            total += wi * f(mid + half * xi)
    return 0.5 * h * total


def green(t, tp, k):
    """Causal Green's function for d^2/dt^2 - k^2 with x(0)=x'(0)=0."""
    if t < tp:
        return 0.0
    if k == 0:
        return t - tp
    return cmath.sinh(k * (t - tp)) / k


def green_dt(t, tp, k):
    """dG/dt for t > tp."""
    if t < tp:
        return 0.0
    if k == 0:
        return 1.0
    return cmath.cosh(k * (t - tp))


def solve_via_green(f, t, k):
    return quad(lambda tp: green(t, tp, k) * f(tp), 0.0, t)


def main():
    # 1. Check jump conditions, homogeneity, and causality for k real, -k, 0, and imaginary
    max_jump_err = 0.0
    max_homog_err = 0.0
    for k in [1.7, -1.7, 0.4, 0.0, 2.3j]:
        tp = 0.8
        assert green(0.0, tp, k) == 0.0
        assert green(tp - 1e-6, tp, k) == 0.0
        # Continuity at t = tp+
        g_plus = green(tp, tp, k)
        max_jump_err = max(max_jump_err, abs(g_plus - 0.0))
        # Derivative jump at t = tp+
        h = 1e-6
        dg_plus = (-3 * green(tp, tp, k) + 4 * green(tp + h, tp, k) - green(tp + 2 * h, tp, k)) / (2 * h)
        max_jump_err = max(max_jump_err, abs(dg_plus - 1.0))
        # Evenness in k
        assert abs(green(1.5, tp, k) - green(1.5, tp, -k)) < 1e-14
        # Homogeneous ODE for t > tp (optimal step for 4th-order 2nd derivative ~ 2e-3)
        for t in [1.1, 1.9, 2.5]:
            hh = 2e-3
            gdd = (
                -green(t - 2 * hh, tp, k)
                + 16 * green(t - hh, tp, k)
                - 30 * green(t, tp, k)
                + 16 * green(t + hh, tp, k)
                - green(t + 2 * hh, tp, k)
            ) / (12 * hh * hh)
            max_homog_err = max(max_homog_err, abs(gdd - (k * k) * green(t, tp, k)))

    # 2. Check ODE and ICs for multiple forcing functions and k values
    max_ode_err = 0.0
    max_ic_err = 0.0
    max_exact_err = 0.0

    for k in [1.4, -0.9, 0.0, 1.8j]:
        # Exact test 1: f(t) = 1
        f1 = lambda t: 1.0
        exact1 = (lambda t, k=k: 0.5 * t * t) if k == 0 else (lambda t, k=k: (cmath.cosh(k * t) - 1.0) / (k * k))
        for t in [0.0, 0.3, 1.1, 2.0]:
            val = solve_via_green(f1, t, k)
            max_exact_err = max(max_exact_err, abs(val - exact1(t)))

        # Exact test 2: resonant forcing f(t) = exp(k*t) (for k != 0)
        if k != 0:
            f_res = lambda t, k=k: cmath.exp(k * t)
            # Solution of x'' - k^2 x = e^{kt}, x(0)=x'(0)=0:
            # x(t) = (t/(2k)) e^{kt} - (1/(2k^2)) sinh(kt)
            exact_res = lambda t, k=k: (t / (2 * k)) * cmath.exp(k * t) - cmath.sinh(k * t) / (2 * k * k)
            for t in [0.25, 0.9, 1.6]:
                val = solve_via_green(f_res, t, k)
                max_exact_err = max(max_exact_err, abs(val - exact_res(t)))

        # General non-elementary test: f(t) = exp(-0.7 t) * cos(2.3 t) + t^2
        f_gen = lambda t: math.exp(-0.7 * t) * math.cos(2.3 * t) + t * t
        # Check ICs at t = 0
        x0 = solve_via_green(f_gen, 0.0, k)
        h0 = 1e-5
        xdot0 = (-3 * solve_via_green(f_gen, 0.0, k) + 4 * solve_via_green(f_gen, h0, k) - solve_via_green(f_gen, 2 * h0, k)) / (2 * h0)
        max_ic_err = max(max_ic_err, abs(x0), abs(xdot0))

        # Check ODE x'' - k^2 x = f(t) at interior points
        hh = 2e-3
        for t in [0.4, 1.0, 1.8]:
            x_t = solve_via_green(f_gen, t, k)
            xdd_t = (
                -solve_via_green(f_gen, t - 2 * hh, k)
                + 16 * solve_via_green(f_gen, t - hh, k)
                - 30 * x_t
                + 16 * solve_via_green(f_gen, t + hh, k)
                - solve_via_green(f_gen, t + 2 * hh, k)
            ) / (12 * hh * hh)
            max_ode_err = max(max_ode_err, abs(xdd_t - (k * k) * x_t - f_gen(t)))

            # Also verify x'(t) == int_0^t G_t(t, tp) f(tp) dtp
            xdot_fd = (
                solve_via_green(f_gen, t - 2 * hh, k)
                - 8 * solve_via_green(f_gen, t - hh, k)
                + 8 * solve_via_green(f_gen, t + hh, k)
                - solve_via_green(f_gen, t + 2 * hh, k)
            ) / (12 * hh)
            xdot_int = quad(lambda tp: green_dt(t, tp, k) * f_gen(tp), 0.0, t)
            assert abs(xdot_fd - xdot_int) < 1e-9

    # 3. Check k -> 0 limit
    lim_err = abs(green(2.3, 0.7, 1e-7) - green(2.3, 0.7, 0.0))

    assert max_jump_err < 1e-9
    assert max_homog_err < 1e-9
    assert max_ic_err < 1e-9
    assert max_ode_err < 1e-8
    assert max_exact_err < 1e-12
    assert lim_err < 1e-12

    print("IB_2005_P2_q5E_check.py: ALL CHECKS PASSED")
    print(f"  max_jump_err  = {max_jump_err:.3e}")
    print(f"  max_homog_err = {max_homog_err:.3e}")
    print(f"  max_ic_err    = {max_ic_err:.3e}")
    print(f"  max_ode_err   = {max_ode_err:.3e}")
    print(f"  max_exact_err = {max_exact_err:.3e}")
    print(f"  k->0 lim_err  = {lim_err:.3e}")


if __name__ == "__main__":
    main()
