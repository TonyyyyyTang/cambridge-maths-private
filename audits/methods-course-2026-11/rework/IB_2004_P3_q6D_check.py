#!/usr/bin/env python3
"""Independent computational verification for IB/2004/P3/q6D."""

import math


# 8-point Gauss-Legendre nodes and weights on [-1, 1]
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
    """Composite 8-point Gauss-Legendre quadrature on [a, b]."""
    h = (b - a) / panels
    total = 0.0
    for k in range(panels):
        x0 = a + k * h
        mid = x0 + 0.5 * h
        half = 0.5 * h
        for xi, wi in zip(_GL8_X, _GL8_W):
            total += wi * f(mid + half * xi)
    return 0.5 * h * total


def make_xc(omega, T, a, b):
    s = math.sin(omega * T)
    assert abs(s) > 1e-12, "Need sin(omega*T) != 0"

    def x(t):
        return (a * math.sin(omega * (T - t)) + b * math.sin(omega * t)) / s

    def xdot(t):
        return omega * (-a * math.cos(omega * (T - t)) + b * math.cos(omega * t)) / s

    def xddot(t):
        return -(omega ** 2) * x(t)

    return x, xdot, xddot


def action(x_fn, xdot_fn, omega, T):
    return quad(lambda t: 0.5 * (xdot_fn(t) ** 2 - (omega ** 2) * (x_fn(t) ** 2)), 0.0, T)


def main():
    cases = [
        (1.3, 1.7, 2.5, -1.1),   # omega*T in (0, pi)
        (2.1, 2.2, -0.8, 3.4),   # omega*T in (pi, 2*pi)
        (0.7, 3.0, 1.0, 1.0),    # symmetric endpoints
        (3.5, 0.4, 0.0, 2.0),    # one zero endpoint
        (0.5, math.pi, 1.0, 1.0),  # omega*T = pi/2
    ]

    max_bc_err = 0.0
    max_ode_err = 0.0
    max_var_err = 0.0
    max_S_ibp_err = 0.0
    max_S_formula_err = 0.0

    for omega, T, a, b in cases:
        xc, xcdot, _ = make_xc(omega, T, a, b)

        # 1. Check boundary conditions
        max_bc_err = max(max_bc_err, abs(xc(0.0) - a), abs(xc(T) - b))

        # 2. Check ODE x'' + omega^2 x = 0 via 4th-order central finite difference
        # Optimal step for 4th-order 2nd derivative in float64 is ~eps^(1/6) ~ 2e-3
        h = 2e-3
        for j in range(1, 20):
            t = T * j / 20.0
            xdd_fd = (-xc(t - 2 * h) + 16 * xc(t - h) - 30 * xc(t) + 16 * xc(t + h) - xc(t + 2 * h)) / (12 * h * h)
            res = abs(xdd_fd + (omega ** 2) * xc(t))
            max_ode_err = max(max_ode_err, res)

        # 3. Check first variation delta S = 0 for multiple variations with eta(0) = eta(T) = 0
        S_xc = action(xc, xcdot, omega, T)
        variations = [
            (
                lambda t, T=T: math.sin(math.pi * t / T),
                lambda t, T=T: (math.pi / T) * math.cos(math.pi * t / T),
            ),
            (
                lambda t, T=T: math.sin(3.0 * math.pi * t / T),
                lambda t, T=T: (3.0 * math.pi / T) * math.cos(3.0 * math.pi * t / T),
            ),
            (
                lambda t, T=T: t * (T - t) * math.exp(0.5 * t),
                lambda t, T=T: ((T - 2.0 * t) + 0.5 * t * (T - t)) * math.exp(0.5 * t),
            ),
        ]
        eps = 1e-2
        for eta, etadot in variations:
            assert abs(eta(0.0)) < 1e-14 and abs(eta(T)) < 1e-14
            S_plus = action(lambda t: xc(t) + eps * eta(t), lambda t: xcdot(t) + eps * etadot(t), omega, T)
            S_minus = action(lambda t: xc(t) - eps * eta(t), lambda t: xcdot(t) - eps * etadot(t), omega, T)
            dS_deps = (S_plus - S_minus) / (2.0 * eps)
            dS_int = quad(lambda t: xcdot(t) * etadot(t) - (omega ** 2) * xc(t) * eta(t), 0.0, T)
            max_var_err = max(max_var_err, abs(dS_deps), abs(dS_int))

            # Exact quadratic identity: S[xc + eps*eta] - S[xc] == eps^2 * S[eta]
            S_eta = action(eta, etadot, omega, T)
            assert abs((S_plus - S_xc) - (eps ** 2) * S_eta) < 1e-12

        # 4. Check non-extremal path has nonzero first variation
        x_lin = lambda t, T=T, a=a, b=b: a + (b - a) * (t / T)
        x_lin_dot = lambda t, T=T, a=a, b=b: (b - a) / T
        eta0, etadot0 = variations[0]
        dS_non = (
            action(lambda t: x_lin(t) + eps * eta0(t), lambda t: x_lin_dot(t) + eps * etadot0(t), omega, T)
            - action(lambda t: x_lin(t) - eps * eta0(t), lambda t: x_lin_dot(t) - eps * etadot0(t), omega, T)
        ) / (2.0 * eps)
        if abs(a + b) > 1e-6:
            assert abs(dS_non) > 1e-3, "Expected nonzero first variation on non-extremal path"

        # 5. Check S[xc] == [1/2 x_c x_c']_0^T == formula
        S_ibp = 0.5 * (xc(T) * xcdot(T) - xc(0.0) * xcdot(0.0))
        S_formula = (omega / (2.0 * math.sin(omega * T))) * ((a * a + b * b) * math.cos(omega * T) - 2.0 * a * b)
        max_S_ibp_err = max(max_S_ibp_err, abs(S_xc - S_ibp))
        max_S_formula_err = max(max_S_formula_err, abs(S_xc - S_formula))

    # 6. Check omega -> 0 free-particle limit
    T0, a0, b0 = 2.3, -1.2, 3.7
    omega_small = 1e-6
    S_small = (omega_small / (2.0 * math.sin(omega_small * T0))) * (
        (a0 * a0 + b0 * b0) * math.cos(omega_small * T0) - 2.0 * a0 * b0
    )
    S_free = ((b0 - a0) ** 2) / (2.0 * T0)
    lim_err = abs(S_small - S_free)

    assert max_bc_err < 1e-14
    assert max_ode_err < 1e-9
    assert max_var_err < 1e-11
    assert max_S_ibp_err < 1e-12
    assert max_S_formula_err < 1e-12
    assert lim_err < 1e-10

    print("IB_2004_P3_q6D_check.py: ALL CHECKS PASSED")
    print(f"  max_bc_err        = {max_bc_err:.3e}")
    print(f"  max_ode_err       = {max_ode_err:.3e}")
    print(f"  max_var_err       = {max_var_err:.3e}")
    print(f"  max_S_ibp_err     = {max_S_ibp_err:.3e}")
    print(f"  max_S_formula_err = {max_S_formula_err:.3e}")
    print(f"  omega->0 lim_err  = {lim_err:.3e}")


if __name__ == "__main__":
    main()
