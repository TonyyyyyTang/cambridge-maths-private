#!/usr/bin/env python3
"""Independent computational verification for IB/2005/P4/q5H."""

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


def quad(f, a, b, panels=64):
    h = (b - a) / panels
    total = 0.0
    for j in range(panels):
        mid = a + (j + 0.5) * h
        half = 0.5 * h
        for xi, wi in zip(_GL8_X, _GL8_W):
            total += wi * f(mid + half * xi)
    return 0.5 * h * total


def main():
    c = 2.7
    L = 1.9

    # Test 1: general asymmetric 2L-periodic Fourier series
    coeffs = [
        (1, 0.7, -1.1),
        (2, -0.4, 0.8),
        (3, 0.35, 0.25),
        (4, -0.15, -0.5),
    ]

    def g1(u):
        s = 0.2  # nonzero constant offset
        for n, an, bn in coeffs:
            k = n * math.pi / L
            s += an * math.cos(k * u) + bn * math.sin(k * u)
        return s

    def g1p(u):
        s = 0.0
        for n, an, bn in coeffs:
            k = n * math.pi / L
            s += k * (-an * math.sin(k * u) + bn * math.cos(k * u))
        return s

    # Exact Parseval integral int_{-L}^{L} g1'(u)^2 du = L * sum (n*pi/L)^2 (an^2 + bn^2)
    parseval_g1 = L * sum(((n * math.pi / L) ** 2) * (an * an + bn * bn) for n, an, bn in coeffs)

    # Test 2: non-polynomial smooth 2L-periodic function
    def g2(u):
        w = math.pi * u / L
        return math.exp(0.6 * math.sin(w)) + 0.4 * math.cos(2 * w - 0.7)

    def g2p(u):
        w = math.pi * u / L
        dw = math.pi / L
        return dw * (0.6 * math.cos(w) * math.exp(0.6 * math.sin(w)) - 0.8 * math.sin(2 * w - 0.7))

    max_pde_err = 0.0
    max_bc_err = 0.0
    max_identity_err = 0.0
    max_const_err = 0.0

    for idx, (g, gp, exact_E) in enumerate([(g1, g1p, parseval_g1), (g2, g2p, None)], start=1):
        # f(u) = -g(u), y(x, t) = f(c*t - x) + g(c*t + x)
        y = lambda x, t, g=g: g(c * t + x) - g(c * t - x)
        yt = lambda x, t, gp=gp: c * (gp(c * t + x) - gp(c * t - x))
        yx = lambda x, t, gp=gp: gp(c * t + x) + gp(c * t - x)

        # Check periodicity of g and g'
        for u in [-3.1, -0.4, 0.7, 2.5]:
            assert abs(g(u + 2 * L) - g(u)) < 1e-13
            assert abs(gp(u + 2 * L) - gp(u)) < 1e-13

        # Check BCs y(0, t) = y(L, t) = 0
        times = [-2.3, -0.7, 0.0, 0.31, 1.25, 3.77]
        for t in times:
            max_bc_err = max(max_bc_err, abs(y(0.0, t)), abs(y(L, t)))

        # Check wave equation y_tt / c^2 - y_xx = 0 via 4th-order central FD
        hx = 1e-3
        ht = hx / c
        for t in [-0.5, 0.4, 1.3]:
            for x in [0.3, 0.9, 1.5]:
                ytt = (-y(x, t - 2 * ht) + 16 * y(x, t - ht) - 30 * y(x, t) + 16 * y(x, t + ht) - y(x, t + 2 * ht)) / (12 * ht * ht)
                yxx = (-y(x - 2 * hx, t) + 16 * y(x - hx, t) - 30 * y(x, t) + 16 * y(x + hx, t) - y(x + 2 * hx, t)) / (12 * hx * hx)
                max_pde_err = max(max_pde_err, abs(ytt / (c * c) - yxx))

        # Check energy identity and constancy in t
        ref_E = quad(lambda u: gp(u) ** 2, -L, L)
        if exact_E is not None:
            assert abs(ref_E - exact_E) < 1e-12

        for t in times:
            lhs = 0.5 * quad(lambda x: (yt(x, t) / c) ** 2 + yx(x, t) ** 2, 0.0, L)
            rhs = quad(lambda x: gp(c * t + x) ** 2, -L, L)
            max_identity_err = max(max_identity_err, abs(lhs - rhs))
            max_const_err = max(max_const_err, abs(rhs - ref_E))

    # Test 3 (non-periodic counterexample): f = -g enforces y(0,t)=0 but NOT y(L,t)=0
    # Identity LHS == RHS still holds, but E(t) is NOT constant in t!
    g_np = lambda u: math.exp(-0.2 * u * u)
    gp_np = lambda u: -0.4 * u * math.exp(-0.2 * u * u)
    e_vals = []
    for t in [0.0, 1.0]:
        lhs = 0.5 * quad(
            lambda x: (gp_np(c * t + x) - gp_np(c * t - x)) ** 2 + (gp_np(c * t + x) + gp_np(c * t - x)) ** 2,
            0.0,
            L,
        )
        rhs = quad(lambda x: gp_np(c * t + x) ** 2, -L, L)
        assert abs(lhs - rhs) < 1e-13
        e_vals.append(rhs)
    assert abs(e_vals[1] - e_vals[0]) > 1e-2, "Expected time variation when g is not 2L-periodic"

    assert max_bc_err < 1e-13
    assert max_pde_err < 1e-8
    assert max_identity_err < 1e-12
    assert max_const_err < 1e-12

    print("IB_2005_P4_q5H_check.py: ALL CHECKS PASSED")
    print(f"  max_bc_err       = {max_bc_err:.3e}")
    print(f"  max_pde_err      = {max_pde_err:.3e}")
    print(f"  max_identity_err = {max_identity_err:.3e}")
    print(f"  max_const_err    = {max_const_err:.3e}")


if __name__ == "__main__":
    main()
