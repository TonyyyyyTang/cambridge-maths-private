#!/usr/bin/env python3
"""Verification script for IB/2012/P1/q14C (Regular Sturm-Liouville system)."""

import cmath
import math


def composite_simpson(f, a, b, N=2000):
    assert N % 2 == 0
    h = (b - a) / N
    s = f(a) + f(b)
    for k in range(1, N):
        x = a + k * h
        s += (4 if k % 2 == 1 else 2) * f(x)
    return s * h / 3


def bisect(f, lo, hi, tol=1e-14, max_iter=200):
    flo = f(lo)
    fhi = f(hi)
    assert flo * fhi <= 0, f"Root not bracketed: f({lo})={flo}, f({hi})={fhi}"
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if abs(fmid) < 1e-15 or (hi - lo) < tol * max(1.0, abs(mid)):
            return mid
        if flo * fmid <= 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return 0.5 * (lo + hi)


def make_robin_function(a, b, A1, A2, B1, B2, c1, c2, c3):
    """
    Construct a smooth (complex) function f(x) = g(x) + alpha + beta*x satisfying
    A1*f(a) + A2*f'(a) = 0 and B1*f(b) + B2*f'(b) = 0,
    where g(x) = c1*exp(1j*x) + c2*cos(2.3*x) + c3*x**3.
    """
    def g(x):
        return c1 * cmath.exp(1j * x) + c2 * math.cos(2.3 * x) + c3 * x ** 3

    def gp(x):
        return 1j * c1 * cmath.exp(1j * x) - 2.3 * c2 * math.sin(2.3 * x) + 3 * c3 * x ** 2

    def gpp(x):
        return -c1 * cmath.exp(1j * x) - (2.3 ** 2) * c2 * math.cos(2.3 * x) + 6 * c3 * x

    # Solve:
    # A1*alpha + (A1*a + A2)*beta = -(A1*g(a) + A2*gp(a))
    # B1*alpha + (B1*b + B2)*beta = -(B1*g(b) + B2*gp(b))
    m11, m12 = A1, A1 * a + A2
    m21, m22 = B1, B1 * b + B2
    det = m11 * m22 - m12 * m21
    assert abs(det) > 1e-8
    r1 = -(A1 * g(a) + A2 * gp(a))
    r2 = -(B1 * g(b) + B2 * gp(b))
    alpha = (r1 * m22 - m12 * r2) / det
    beta = (m11 * r2 - r1 * m21) / det

    def f(x):
        return g(x) + alpha + beta * x

    def fp(x):
        return gp(x) + beta

    def fpp(x):
        return gpp(x)

    assert abs(A1 * f(a) + A2 * fp(a)) < 1e-12
    assert abs(B1 * f(b) + B2 * fp(b)) < 1e-12
    return f, fp, fpp


def test_general_self_adjointness():
    a, b = 0.3, 1.7
    A1, A2 = 1.5, -0.7
    B1, B2 = 2.1, 0.9

    def p(x):
        return 2.0 + math.cos(x)

    def pp(x):
        return -math.sin(x)

    def q(x):
        return x * x + 1.0

    def L_apply(f, fp, fpp):
        return lambda x: -(pp(x) * fp(x) + p(x) * fpp(x)) + q(x) * f(x)

    u, up, upp = make_robin_function(a, b, A1, A2, B1, B2, 1.0 + 0.4j, -0.7 + 1.2j, 0.3 - 0.5j)
    v, vp, vpp = make_robin_function(a, b, A1, A2, B1, B2, -0.6 + 0.9j, 1.1 - 0.3j, -0.8 + 0.2j)
    Lu = L_apply(u, up, upp)
    Lv = L_apply(v, vp, vpp)

    ip1 = composite_simpson(lambda x: u(x).conjugate() * Lv(x), a, b, N=2000)
    ip2 = composite_simpson(lambda x: Lu(x).conjugate() * v(x), a, b, N=2000)
    err = abs(ip1 - ip2)
    assert err < 1e-11, f"General S-L self-adjointness failed: {err}"
    print(f"[PASS] General regular S-L self-adjointness |<u,Lv> - <Lu,v>| = {err:.3e}")


def test_explicit_q0_system():
    configs = [(1.0, 1.0), (0.5, 2.0), (5.0, 0.7), (0.2, 3.1)]
    for k, ell in configs:
        # Check no root in (0, pi/(2*ell)] for F(mu) = k*sin(mu*ell) + mu*cos(mu*ell)
        # Note for mu in (0, pi/(2*ell)], sin > 0 and cos >= 0, so F(mu) > 0 strictly.
        for j in range(1, 200):
            mu = j * (0.5 * math.pi / ell) / 200.0
            assert k * math.sin(mu * ell) + mu * math.cos(mu * ell) > 0

        # Find first 50 roots in ((n - 0.5)*pi/ell, n*pi/ell)
        mus = []
        lams = []
        for n in range(1, 51):
            lo = ((n - 0.5) * math.pi + 1e-12) / ell
            hi = (n * math.pi - 1e-12) / ell
            # Check F(mu) = tan(mu*ell) + mu/k changes sign from - to +
            assert math.tan(lo * ell) + lo / k < 0
            assert math.tan(hi * ell) + hi / k > 0
            mu_n = bisect(lambda mu: math.tan(mu * ell) + mu / k, lo, hi, tol=1e-15)
            # Residual of k*sin(mu*ell) + mu*cos(mu*ell) = 0 (relative to scale k + mu)
            res = abs(k * math.sin(mu_n * ell) + mu_n * math.cos(mu_n * ell)) / (k + mu_n)
            assert res < 1e-12, f"Residual too large for n={n}: {res}"
            mus.append(mu_n)
            lams.append(mu_n ** 2)

        # Verify strict ordering
        for i in range(len(lams) - 1):
            assert lams[i] < lams[i + 1]

        # Verify Rayleigh quotient and orthogonality for first 6 eigenfunctions
        for i in range(6):
            mu_i = mus[i]
            yi = lambda x, m=mu_i: math.sin(m * x)
            yip = lambda x, m=mu_i: m * math.cos(m * x)

            # Check BCs
            assert abs(yi(0.0)) < 1e-14
            assert abs(k * yi(ell) + yip(ell)) < 1e-11

            # Check Rayleigh quotient both via Simpson (N=8000) and exact closed-form integrals:
            # lambda = (int_0^ell (y'^2) dx - [y y']_0^ell) / int_0^ell y^2 dx
            #        = (int_0^ell (y'^2) dx + k y(ell)^2) / int_0^ell y^2 dx
            int_yp2 = composite_simpson(lambda x: yip(x) ** 2, 0.0, ell, N=8000)
            num1 = int_yp2 - (yi(ell) * yip(ell) - yi(0.0) * yip(0.0))
            num2 = int_yp2 + k * yi(ell) ** 2
            den = composite_simpson(lambda x: yi(x) ** 2, 0.0, ell, N=8000)
            assert abs(num1 - num2) / max(1.0, abs(num1)) < 1e-11
            assert abs(num1 / den - lams[i]) / lams[i] < 1e-11

            # Also exact analytic check of the Rayleigh quotient
            exact_yp2 = (mu_i ** 2) * (0.5 * ell + math.sin(2 * mu_i * ell) / (4 * mu_i))
            exact_y2 = 0.5 * ell - math.sin(2 * mu_i * ell) / (4 * mu_i)
            exact_rq = (exact_yp2 + k * yi(ell) ** 2) / exact_y2
            assert abs(exact_rq - lams[i]) / lams[i] < 1e-12

            for j in range(i + 1, 6):
                mu_j = mus[j]
                ip = composite_simpson(
                    lambda x, mi=mu_i, mj=mu_j: math.sin(mi * x) * math.sin(mj * x),
                    0.0,
                    ell,
                    N=2000,
                )
                assert abs(ip) < 1e-10, f"Orthogonality failed for ({i},{j}): {ip}"

        # Check asymptotics for large n = 50, 200, 1000
        for n in [50, 200, 1000]:
            X = (n - 0.5) * math.pi
            # Solve tan(delta) = k*ell / (X + delta) for delta in (0, pi/2)
            delta = bisect(lambda d: math.tan(d) - k * ell / (X + d), 1e-14, 0.49 * math.pi)
            xi_n = X + delta
            lam_n = (xi_n / ell) ** 2
            # Two-term asymptotic formula: X^2/ell^2 + 2k/ell
            asymp2 = (X / ell) ** 2 + 2.0 * k / ell
            # Three-term asymptotic formula: X^2/ell^2 + 2k/ell - ((k*ell)^2 + (2/3)*(k*ell)^3)/(ell^2 * X^2)
            asymp3 = asymp2 - ((k * ell) ** 2 + (2.0 / 3.0) * (k * ell) ** 3) / (ell ** 2 * X ** 2)
            err2 = abs(lam_n - asymp2)
            err3 = abs(lam_n - asymp3)
            # err2 should be O(1/n^2), err3 should be O(1/n^4)
            assert err2 < 50.0 * max(1.0, (k * ell) ** 3) / (ell ** 2 * n ** 2)
            assert err3 < err2 * 0.05 or err3 < 1e-9

    print("[PASS] Explicit q=0 system: roots, ordering, BCs, Rayleigh quotient, orthogonality, and O(1/n^2) asymptotics verified")


def test_nonzero_q_positive_eigenvalue():
    """Verify strict positivity & Rayleigh quotient for variable q(x) = x^2 >= 0 via RK4 shooting."""
    ell = 1.5
    k = 0.8

    def shoot(lam, steps=1000):
        # Solve -y'' + x^2 y = lam y, y(0)=0, y'(0)=1
        h = ell / steps
        y, yp = 0.0, 1.0
        ys = [(0.0, y, yp)]
        for i in range(steps):
            x = i * h
            def f(xv, yv, ypv):
                return ypv, (xv * xv - lam) * yv
            k1y, k1p = f(x, y, yp)
            k2y, k2p = f(x + 0.5 * h, y + 0.5 * h * k1y, yp + 0.5 * h * k1p)
            k3y, k3p = f(x + 0.5 * h, y + 0.5 * h * k2y, yp + 0.5 * h * k2p)
            k4y, k4p = f(x + h, y + h * k3y, yp + h * k3p)
            y += (h / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y)
            yp += (h / 6.0) * (k1p + 2 * k2p + 2 * k3p + k4p)
            ys.append(((i + 1) * h, y, yp))
        return k * y + yp, ys

    # Confirm no eigenvalue for lam <= 0 (shooting residual stays strictly positive)
    for test_lam in [-10.0, -2.0, -0.5, 0.0]:
        res, _ = shoot(test_lam)
        assert res > 0, f"Unexpected non-positive eigenvalue near lam={test_lam}"

    # Find first positive eigenvalue
    lam1 = bisect(lambda l: shoot(l, steps=600)[0], 0.0, 5.0, tol=1e-12)
    assert lam1 > 0
    res, ys = shoot(lam1, steps=2000)
    assert abs(res) < 1e-10

    # Compute Rayleigh quotient via Simpson's rule on grid
    h = ell / 2000
    num = 0.0
    den = 0.0
    for i, (x, y, yp) in enumerate(ys):
        w = 1 if (i == 0 or i == len(ys) - 1) else (4 if i % 2 == 1 else 2)
        num += w * (yp * yp + x * x * y * y)
        den += w * (y * y)
    num = num * h / 3.0 + k * ys[-1][1] ** 2
    den = den * h / 3.0
    rq = num / den
    err = abs(rq - lam1)
    assert err < 1e-10, f"Rayleigh quotient mismatch for variable q(x)=x^2: {err}"
    print(f"[PASS] Variable q(x)=x^2 >= 0: lambda_1 = {lam1:.10f} > 0, |RQ - lambda_1| = {err:.3e}")


def main():
    test_general_self_adjointness()
    test_explicit_q0_system()
    test_nonzero_q_positive_eigenvalue()


if __name__ == "__main__":
    main()
