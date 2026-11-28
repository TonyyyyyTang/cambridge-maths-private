#!/usr/bin/env python3
"""Verification script for IB/2014/P1/q14D (Legendre polynomials & axisymmetric Laplace)."""

import math
from fractions import Fraction


def poly_deriv(p):
    if len(p) <= 1:
        return [Fraction(0)]
    return [p[k] * k for k in range(1, len(p))]


def poly_mul(a, b):
    res = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        for j, cb in enumerate(b):
            res[i + j] += ca * cb
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res


def poly_add(a, b):
    n = max(len(a), len(b))
    res = [Fraction(0)] * n
    for i in range(len(a)):
        res[i] += a[i]
    for i in range(len(b)):
        res[i] += b[i]
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res


def poly_scale(a, s):
    return [c * s for c in a]


def poly_eval(p, x):
    val = Fraction(0) if isinstance(x, Fraction) else 0.0
    for c in reversed(p):
        val = val * x + (c if isinstance(x, Fraction) else float(c))
    return val


def legendre_from_x0_recurrence(n):
    """Build polynomial solution of degree n from power series around x=0 and normalise y(1)=1."""
    coeffs = [Fraction(0)] * (n + 1)
    start = n % 2
    coeffs[start] = Fraction(1)
    for k in range(start, n - 1, 2):
        coeffs[k + 2] = -Fraction((n - k) * (n + k + 1), (k + 2) * (k + 1)) * coeffs[k]
    val1 = poly_eval(coeffs, Fraction(1))
    assert val1 != 0, f"y_{n}(1) vanished!"
    return poly_scale(coeffs, Fraction(1) / val1)


def legendre_from_x1_recurrence(n):
    """Build polynomial solution from Taylor expansion in t = x - 1 with c_0 = y(1) = 1."""
    c = [Fraction(1)]
    for k in range(n):
        c.append(Fraction((n - k) * (n + k + 1), 2 * (k + 1) ** 2) * c[-1])
    # Convert sum c_k (x - 1)^k to powers of x
    res = [Fraction(0)]
    t_pow = [Fraction(1)]  # (x - 1)^k
    xm1 = [Fraction(-1), Fraction(1)]
    for k in range(n + 1):
        res = poly_add(res, poly_scale(t_pow, c[k]))
        t_pow = poly_mul(t_pow, xm1)
    return res


def spherical_laplacian_fd(U, r, theta, hr=1e-4, ht=1e-4):
    """4th-order central finite-difference evaluation of axisymmetric spherical Laplacian."""
    U0 = U(r, theta)
    Ur = (-U(r + 2 * hr, theta) + 8 * U(r + hr, theta) - 8 * U(r - hr, theta) + U(r - 2 * hr, theta)) / (12 * hr)
    Urr = (-U(r + 2 * hr, theta) + 16 * U(r + hr, theta) - 30 * U0 + 16 * U(r - hr, theta) - U(r - 2 * hr, theta)) / (12 * hr * hr)
    Ut = (-U(r, theta + 2 * ht) + 8 * U(r, theta + ht) - 8 * U(r, theta - ht) + U(r, theta - 2 * ht)) / (12 * ht)
    Utt = (-U(r, theta + 2 * ht) + 16 * U(r, theta + ht) - 30 * U0 + 16 * U(r, theta - ht) - U(r, theta - 2 * ht)) / (12 * ht * ht)
    return Urr + (2.0 / r) * Ur + (1.0 / (r * r)) * (Utt + (math.cos(theta) / math.sin(theta)) * Ut)


def main():
    # Part (a): Verify P_n(x) for n = 0..8
    P = [legendre_from_x0_recurrence(n) for n in range(9)]
    for n, pn in enumerate(P):
        assert len(pn) == n + 1 and pn[-1] != 0, f"P_{n} does not have exact degree {n}"
        assert poly_eval(pn, Fraction(1)) == Fraction(1)
        # Check agreement with x=1 expansion
        pn_x1 = legendre_from_x1_recurrence(n)
        assert pn == pn_x1, f"Mismatch between x=0 and x=1 recurrences for n={n}"
        # Check Legendre ODE (1 - x^2) y'' - 2x y' + n(n+1) y = 0
        one_m_x2 = [Fraction(1), Fraction(0), Fraction(-1)]
        m2x = [Fraction(0), Fraction(-2)]
        ode = poly_add(
            poly_add(
                poly_mul(one_m_x2, poly_deriv(poly_deriv(pn))),
                poly_mul(m2x, poly_deriv(pn)),
            ),
            poly_scale(pn, Fraction(n * (n + 1))),
        )
        assert ode == [Fraction(0)], f"Legendre ODE failed for n={n}"

    assert P[0] == [Fraction(1)]
    assert P[1] == [Fraction(0), Fraction(1)]
    assert P[2] == [Fraction(-1, 2), Fraction(0), Fraction(3, 2)]
    print("[PASS] Part (a): P_0 = 1, P_1 = x, P_2 = (3x^2 - 1)/2 and exact Legendre ODE verified for n = 0..8")

    # Part (b): Verify separated modes (A_n r^n + B_n r^{-(n+1)}) P_n(cos theta) satisfy Laplace's eq
    max_mode_err = 0.0
    for n in range(7):
        An = 0.7 - 0.1 * n
        Bn = 1.3 + 0.2 * n
        mode = lambda r, th, n=n, An=An, Bn=Bn: (An * (r ** n) + Bn * (r ** (-(n + 1)))) * poly_eval(P[n], math.cos(th))
        for r in [0.8, 1.5, 2.7]:
            for th in [0.4, 1.1, 2.2]:
                err = abs(spherical_laplacian_fd(mode, r, th))
                scale = max(1.0, abs(mode(r, th)))
                assert err / scale < 1e-5, f"Mode n={n} failed Laplace check: {err}"
                max_mode_err = max(max_mode_err, err / scale)
    print(f"[PASS] Part (b): Separated modes n=0..6 satisfy spherical Laplace (max rel FD err = {max_mode_err:.3e})")

    # Part (b): Verify specific BVP solution U(r, theta) = v0 * (r + r0^3 / (2 r^2)) * cos(theta)
    max_lap_err = 0.0
    max_bc_err = 0.0
    max_cart_err = 0.0
    for v0, r0 in [(1.0, 1.0), (2.5, 0.6), (-1.7, 3.2)]:
        def U(r, th, v0=v0, r0=r0):
            return v0 * (r + (r0 ** 3) / (2.0 * r * r)) * math.cos(th)

        # Check spherical Laplace
        for r in [r0, 1.5 * r0, 4.0 * r0]:
            for th in [0.3, 1.0, 1.9, 2.7]:
                err = abs(spherical_laplacian_fd(U, r, th))
                max_lap_err = max(max_lap_err, err)
                assert err < 1e-6

        # Check Neumann BC dU/dr = 0 at r = r0
        hr = 1e-6 * r0
        for th in [0.0, 0.5, 1.2, math.pi / 2, 2.4, math.pi]:
            dUdr = (U(r0 + hr, th) - U(r0 - hr, th)) / (2 * hr)
            max_bc_err = max(max_bc_err, abs(dUdr))
            assert abs(dUdr) < 1e-8

        # Check far-field U(r, th) - v0*r*cos(th) -> 0 as r -> inf
        for r_big in [1e2 * r0, 1e3 * r0]:
            for th in [0.0, 0.7, 2.1, math.pi]:
                diff = abs(U(r_big, th) - v0 * r_big * math.cos(th))
                expected = abs(v0) * (r0 ** 3) / (2.0 * r_big ** 2) * abs(math.cos(th))
                assert abs(diff - expected) < 1e-10
                assert diff <= abs(v0) * (r0 ** 3) / (2.0 * r_big ** 2) + 1e-10

        # Cross-check in 3D Cartesian coordinates: U(x,y,z) = v0 * z * (1 + r0^3 / (2*(x^2+y^2+z^2)^(3/2)))
        def U_cart(x, y, z, v0=v0, r0=r0):
            r = math.sqrt(x * x + y * y + z * z)
            return v0 * z * (1.0 + (r0 ** 3) / (2.0 * r ** 3))

        h = 1e-3 * r0
        for x, y, z in [(r0, 0.5 * r0, 1.2 * r0), (-0.8 * r0, 1.4 * r0, -0.7 * r0)]:
            U0 = U_cart(x, y, z)
            uxx = (-U_cart(x + 2 * h, y, z) + 16 * U_cart(x + h, y, z) - 30 * U0 + 16 * U_cart(x - h, y, z) - U_cart(x - 2 * h, y, z)) / (12 * h * h)
            uyy = (-U_cart(x, y + 2 * h, z) + 16 * U_cart(x, y + h, z) - 30 * U0 + 16 * U_cart(x, y - h, z) - U_cart(x, y - 2 * h, z)) / (12 * h * h)
            uzz = (-U_cart(x, y, z + 2 * h) + 16 * U_cart(x, y, z + h) - 30 * U0 + 16 * U_cart(x, y, z - h) - U_cart(x, y, z - 2 * h)) / (12 * h * h)
            lap = uxx + uyy + uzz
            max_cart_err = max(max_cart_err, abs(lap))
            assert abs(lap) < 1e-8

    print(
        f"[PASS] Part (b): BVP solution U = v0(r + r0^3/(2r^2))cos(theta) verified: "
        f"max sph Lap err = {max_lap_err:.3e}, max BC err = {max_bc_err:.3e}, max Cart Lap err = {max_cart_err:.3e}"
    )


if __name__ == "__main__":
    main()
