#!/usr/bin/env python3
"""Verification script for IB/2011/P2/q5A (Legendre equation as S-L operator)."""

import cmath
import math
from fractions import Fraction


def poly_deriv(p):
    """Differentiate polynomial represented as list of Fraction coefficients [c0, c1, ...]."""
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
    val = Fraction(0)
    for c in reversed(p):
        val = val * x + c
    return val


def poly_integrate_m1_to_1(p):
    """Exact integral of polynomial p(x) over [-1, 1]."""
    total = Fraction(0)
    for k, c in enumerate(p):
        if k % 2 == 0:
            total += c * Fraction(2, k + 1)
    return total


def apply_L_poly(p):
    """Apply L y = -d/dx((1 - x^2) dy/dx) to polynomial p."""
    one_minus_x2 = [Fraction(1), Fraction(0), Fraction(-1)]
    dp = poly_deriv(p)
    prod = poly_mul(one_minus_x2, dp)
    dprod = poly_deriv(prod)
    return poly_scale(dprod, Fraction(-1))


def legendre_rodrigues(n):
    """Compute P_n(x) = 1/(2^n n!) d^n/dx^n (x^2 - 1)^n exactly."""
    base = [Fraction(-1), Fraction(0), Fraction(1)]
    poly = [Fraction(1)]
    for _ in range(n):
        poly = poly_mul(poly, base)
    for _ in range(n):
        poly = poly_deriv(poly)
    scale = Fraction(1, (2 ** n) * math.factorial(n))
    return poly_scale(poly, scale)


def composite_simpson(f, a, b, N=2000):
    assert N % 2 == 0
    h = (b - a) / N
    s = f(a) + f(b)
    for k in range(1, N):
        x = a + k * h
        s += (4 if k % 2 == 1 else 2) * f(x)
    return s * h / 3


def main():
    max_n = 8
    P = [legendre_rodrigues(n) for n in range(max_n + 1)]

    # 1. Check P_n(1) = 1 and L P_n = n(n+1) P_n exactly
    for n, pn in enumerate(P):
        assert poly_eval(pn, Fraction(1)) == Fraction(1), f"P_{n}(1) != 1"
        Lpn = apply_L_poly(pn)
        rhs = poly_scale(pn, Fraction(n * (n + 1)))
        diff = poly_add(Lpn, poly_scale(rhs, Fraction(-1)))
        assert diff == [Fraction(0)], f"Eigenvalue eq failed for n={n}: {diff}"

        # Also check expanded form (1-x^2)y'' - 2xy' + n(n+1)y = 0
        one_minus_x2 = [Fraction(1), Fraction(0), Fraction(-1)]
        minus_2x = [Fraction(0), Fraction(-2)]
        lhs = poly_add(
            poly_add(
                poly_mul(one_minus_x2, poly_deriv(poly_deriv(pn))),
                poly_mul(minus_2x, poly_deriv(pn)),
            ),
            poly_scale(pn, Fraction(n * (n + 1))),
        )
        assert lhs == [Fraction(0)], f"Expanded Legendre ODE failed for n={n}"

    print(f"[PASS] Exact eigenvalue equation L P_n = n(n+1) P_n verified for n = 0..{max_n}")

    # 2. Exact orthogonality and normalisation <P_m, P_n> = 2/(2n+1) delta_{mn}
    for m in range(max_n + 1):
        for n in range(max_n + 1):
            ip = poly_integrate_m1_to_1(poly_mul(P[m], P[n]))
            expected = Fraction(2, 2 * n + 1) if m == n else Fraction(0)
            assert ip == expected, f"<P_{m}, P_{n}> = {ip}, expected {expected}"

    print(f"[PASS] Exact orthogonality <P_m, P_n> = 2/(2n+1) delta_mn verified for 0 <= m,n <= {max_n}")

    # 3. Exact self-adjointness <u, Lv> = <Lu, v> for arbitrary polynomials u, v
    test_polys = [
        [Fraction(3), Fraction(-2), Fraction(5), Fraction(1, 7)],
        [Fraction(-1, 3), Fraction(4), Fraction(0), Fraction(-2), Fraction(9, 5)],
        [Fraction(2), Fraction(1), Fraction(-3), Fraction(7), Fraction(-1), Fraction(4)],
    ]
    for u in test_polys:
        for v in test_polys:
            lhs = poly_integrate_m1_to_1(poly_mul(u, apply_L_poly(v)))
            rhs = poly_integrate_m1_to_1(poly_mul(apply_L_poly(u), v))
            assert lhs == rhs, f"Self-adjointness failed: {lhs} != {rhs}"
    print("[PASS] Exact self-adjointness <u, Lv> = <Lu, v> verified for test polynomials")

    # 4. Numerical self-adjointness for smooth complex non-polynomial functions
    def u(x):
        return cmath.exp(1j * x) * math.cos(2 * x) + 0.5j * math.cosh(x)

    def up(x):
        return (
            1j * cmath.exp(1j * x) * math.cos(2 * x)
            - 2 * cmath.exp(1j * x) * math.sin(2 * x)
            + 0.5j * math.sinh(x)
        )

    def upp(x):
        return (
            -5 * cmath.exp(1j * x) * math.cos(2 * x)
            - 4j * cmath.exp(1j * x) * math.sin(2 * x)
            + 0.5j * math.cosh(x)
        )

    def v(x):
        return (1 + 2j * x) * math.sin(3 * x) + math.exp(-x)

    def vp(x):
        return 2j * math.sin(3 * x) + 3 * (1 + 2j * x) * math.cos(3 * x) - math.exp(-x)

    def vpp(x):
        return (
            12j * math.cos(3 * x)
            - 9 * (1 + 2j * x) * math.sin(3 * x)
            + math.exp(-x)
        )

    def Lu(x):
        return -(1 - x * x) * upp(x) + 2 * x * up(x)

    def Lv(x):
        return -(1 - x * x) * vpp(x) + 2 * x * vp(x)

    ip_u_Lv = composite_simpson(lambda x: u(x).conjugate() * Lv(x), -1.0, 1.0, N=4000)
    ip_Lu_v = composite_simpson(lambda x: Lu(x).conjugate() * v(x), -1.0, 1.0, N=4000)
    err = abs(ip_u_Lv - ip_Lu_v)
    assert err < 1e-11, f"Complex self-adjointness discrepancy: {err}"
    print(f"[PASS] Complex smooth self-adjointness |<u,Lv> - <Lu,v>| = {err:.3e}")

    # 5. Verify boundary term is NON-ZERO for singular Q_0(x) = 0.5 ln((1+x)/(1-x)) with u = P_1(x) = x
    # Boundary term B(x) = (1-x^2)(u'(x) v(x) - u(x) v'(x))
    # For u=x, v=Q_0: u'=1, v'=1/(1-x^2), so B(x) = (1-x^2)Q_0(x) - x -> -1 at x->1-, +1 at x->-1+
    eps = 1e-10
    def B_sing(x):
        Q0 = 0.5 * math.log((1 + x) / (1 - x))
        return (1 - x * x) * (1.0 * Q0) - x

    b_jump = B_sing(1 - eps) - B_sing(-1 + eps)
    assert abs(b_jump - (-2.0)) < 1e-7, f"Unexpected singular boundary term: {b_jump}"
    print(f"[PASS] Singular function Q_0 boundary term [p(u'v - uv')]_(-1+eps)^(1-eps) = {b_jump:.8f} (expected -2)")


if __name__ == "__main__":
    main()
