#!/usr/bin/env python3
"""Computational verification for IB/methods/sheet-2/qx24 (2007-08 Sheet 2, Q8)."""

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


def fd1(func, x, h=1e-5):
    """4th-order central finite difference first derivative."""
    return (-func(x + 2 * h) + 8 * func(x + h) - 8 * func(x - h) + func(x - 2 * h)) / (12 * h)


def check_part_i():
    """Check Sturm-Liouville form of Legendre's equation."""
    # Arbitrary smooth test function y(x) = exp(0.7x) * cos(1.3x)
    y = lambda x: math.exp(0.7 * x) * math.cos(1.3 * x)
    yp = lambda x: math.exp(0.7 * x) * (0.7 * math.cos(1.3 * x) - 1.3 * math.sin(1.3 * x))
    ypp = lambda x: math.exp(0.7 * x) * ((0.49 - 1.69) * math.cos(1.3 * x) - 2 * 0.7 * 1.3 * math.sin(1.3 * x))

    max_err = 0.0
    for n in [0, 1, 2, 3.5]:
        for x in [-0.7, -0.2, 0.35, 0.85]:
            orig = (1.0 - x * x) * ypp(x) - 2.0 * x * yp(x) + n * (n + 1.0) * y(x)
            sl = fd1(lambda s: (1.0 - s * s) * yp(s), x) + n * (n + 1.0) * y(x)
            err = abs(orig - sl)
            max_err = max(max_err, err)
            assert err < 1e-9

    # Check Legendre polynomials P_2(x) = (3x^2 - 1)/2 and P_3(x) = (5x^3 - 3x)/2
    P2 = lambda x: 0.5 * (3 * x * x - 1)
    P2p = lambda x: 3 * x
    P3 = lambda x: 0.5 * (5 * x ** 3 - 3 * x)
    P3p = lambda x: 0.5 * (15 * x * x - 3)
    for x in [-0.6, 0.1, 0.8]:
        assert abs(fd1(lambda s: (1 - s * s) * P2p(s), x) + 6 * P2(x)) < 1e-10
        assert abs(fd1(lambda s: (1 - s * s) * P3p(s), x) + 12 * P3(x)) < 1e-10
    return max_err


def check_part_ii():
    """Check Sturm-Liouville form of the hypergeometric equation on all three intervals."""
    y = lambda x: math.exp(0.4 * x) * math.sin(1.1 * x + 0.3)
    yp = lambda x: math.exp(0.4 * x) * (0.4 * math.sin(1.1 * x + 0.3) + 1.1 * math.cos(1.1 * x + 0.3))
    ypp = lambda x: math.exp(0.4 * x) * ((0.16 - 1.21) * math.sin(1.1 * x + 0.3) + 0.88 * math.cos(1.1 * x + 0.3))

    max_err = 0.0
    params = [(0.4, 1.3, 0.7), (-0.6, 2.1, 1.8), (2.0, -3.0, 1.5)]
    for a, b, c in params:
        # Interval (0, 1): d/dx [x^c (1-x)^{1+a+b-c} y'] - ab x^{c-1} (1-x)^{a+b-c} y = 0
        for x in [0.15, 0.45, 0.82]:
            orig = x * (x - 1) * ypp(x) + ((1 + a + b) * x - c) * yp(x) + a * b * y(x)
            mu = -(x ** (c - 1)) * ((1 - x) ** (a + b - c))
            sl = fd1(lambda s: (s ** c) * ((1 - s) ** (1 + a + b - c)) * yp(s), x) - a * b * (x ** (c - 1)) * ((1 - x) ** (a + b - c)) * y(x)
            err = abs(mu * orig - sl)
            max_err = max(max_err, err)
            assert err < 1e-8, f"(0,1) mismatch: {err}"

        # Interval (1, inf): d/dx [x^c (x-1)^{1+a+b-c} y'] + ab x^{c-1} (x-1)^{a+b-c} y = 0
        for x in [1.3, 2.2, 3.5]:
            orig = x * (x - 1) * ypp(x) + ((1 + a + b) * x - c) * yp(x) + a * b * y(x)
            mu = (x ** (c - 1)) * ((x - 1) ** (a + b - c))
            sl = fd1(lambda s: (s ** c) * ((s - 1) ** (1 + a + b - c)) * yp(s), x) + a * b * (x ** (c - 1)) * ((x - 1) ** (a + b - c)) * y(x)
            err = abs(mu * orig - sl)
            max_err = max(max_err, err)
            assert err < 1e-7, f"(1,inf) mismatch: {err}"

        # Interval (-inf, 0): d/dx [(-x)^c (1-x)^{1+a+b-c} y'] + ab (-x)^{c-1} (1-x)^{a+b-c} y = 0
        for x in [-2.3, -0.7, -0.25]:
            orig = x * (x - 1) * ypp(x) + ((1 + a + b) * x - c) * yp(x) + a * b * y(x)
            mu = ((-x) ** (c - 1)) * ((1 - x) ** (a + b - c))
            sl = fd1(lambda s: ((-s) ** c) * ((1 - s) ** (1 + a + b - c)) * yp(s), x) + a * b * ((-x) ** (c - 1)) * ((1 - x) ** (a + b - c)) * y(x)
            err = abs(mu * orig - sl)
            max_err = max(max_err, err)
            assert err < 1e-7, f"(-inf,0) mismatch: {err}"

    # Verify explicit polynomial solution for a=-2, b=3, c=1: y(x) = 1 - 6x + 6x^2 = P_2(1-2x)
    a, b, c = -2, 3, 1
    ypoly = lambda x: 1.0 - 6.0 * x + 6.0 * x * x
    ypolyp = lambda x: -6.0 + 12.0 * x
    for x in [0.2, 0.5, 0.85]:
        sl_val = fd1(lambda s: s * (1 - s) * ypolyp(s), x) - a * b * ypoly(x)
        assert abs(sl_val) < 1e-10

    return max_err


def shoot_part_iii(lam, steps=1000):
    """Solve y'' + 4y' + (4+lam)y = 0, y(0)=0, y'(0)=1 via RK4 and return y(1)."""
    h = 1.0 / steps
    y, v = 0.0, 1.0
    for _ in range(steps):
        def f(yy, vv):
            return vv, -4.0 * vv - (4.0 + lam) * yy

        k1y, k1v = f(y, v)
        k2y, k2v = f(y + 0.5 * h * k1y, v + 0.5 * h * k1v)
        k3y, k3v = f(y + 0.5 * h * k2y, v + 0.5 * h * k2v)
        k4y, k4v = f(y + h * k3y, v + h * k3v)
        y += (h / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y)
        v += (h / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
    return y


def check_part_iii():
    """Check SL recasting, eigenvalues, eigenfunctions, shooting, and orthogonality."""
    # Check SL recasting for arbitrary test function
    y = lambda x: math.sin(2.3 * x) + 0.5 * math.exp(-x)
    yp = lambda x: 2.3 * math.cos(2.3 * x) - 0.5 * math.exp(-x)
    ypp = lambda x: -5.29 * math.sin(2.3 * x) + 0.5 * math.exp(-x)
    for lam in [-3.0, 0.0, 7.5]:
        for x in [0.1, 0.5, 0.9]:
            orig = ypp(x) + 4.0 * yp(x) + (4.0 + lam) * y(x)
            sl = fd1(lambda s: math.exp(4.0 * s) * yp(s), x) + (4.0 + lam) * math.exp(4.0 * x) * y(x)
            assert abs(math.exp(4.0 * x) * orig - sl) < 1e-8

    # Check eigenfunctions y_n(x) = e^{-2x} sin(n pi x) with lam_n = n^2 pi^2
    for n in range(1, 7):
        lam = (n * math.pi) ** 2
        yn = lambda x, n=n: math.exp(-2.0 * x) * math.sin(n * math.pi * x)
        ynp = lambda x, n=n: math.exp(-2.0 * x) * (n * math.pi * math.cos(n * math.pi * x) - 2.0 * math.sin(n * math.pi * x))
        assert abs(yn(0.0)) < 1e-15
        assert abs(yn(1.0)) < 1e-14
        for x in [0.15, 0.43, 0.77]:
            sl = fd1(lambda s: math.exp(4.0 * s) * ynp(s), x) + (4.0 + lam) * math.exp(4.0 * x) * yn(x)
            assert abs(sl) < 1e-7

    # Shooting scan on [-20, 100] to confirm eigenvalues are ONLY n^2 pi^2 (n=1,2,3)
    grid = [-20.0 + 0.5 * k for k in range(241)]
    roots = []
    prev_l = grid[0]
    prev_v = shoot_part_iii(prev_l)
    for cur_l in grid[1:]:
        cur_v = shoot_part_iii(cur_l)
        if prev_v * cur_v <= 0.0:
            # Bisection
            lo, hi = prev_l, cur_l
            flo = prev_v
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                fmid = shoot_part_iii(mid)
                if flo * fmid <= 0.0:
                    hi = mid
                else:
                    lo, flo = mid, fmid
            roots.append(0.5 * (lo + hi))
        prev_l, prev_v = cur_l, cur_v

    expected = [(n * math.pi) ** 2 for n in [1, 2, 3]]
    assert len(roots) == 3, f"Expected 3 eigenvalues in [-20,100], got {roots}"
    max_eig_err = max(abs(r - e) for r, e in zip(roots, expected))
    assert max_eig_err < 1e-7, f"Eigenvalue shooting mismatch: {max_eig_err}"

    # Check weighted orthogonality: int_0^1 e^{4x} y_m(x) y_n(x) dx = (1/2) delta_{mn}
    max_orth_err = 0.0
    for m in range(1, 7):
        for n in range(1, 7):
            val = gauss_legendre_8(
                lambda x: math.exp(4.0 * x)
                * (math.exp(-2.0 * x) * math.sin(m * math.pi * x))
                * (math.exp(-2.0 * x) * math.sin(n * math.pi * x)),
                0.0,
                1.0,
                panels=40,
            )
            target = 0.5 if m == n else 0.0
            err = abs(val - target)
            max_orth_err = max(max_orth_err, err)
            assert err < 1e-12, f"Orthogonality error m={m}, n={n}: {err}"

    return max_eig_err, max_orth_err


def main():
    err_i = check_part_i()
    err_ii = check_part_ii()
    eig_err, orth_err = check_part_iii()
    print("IB_methods_sheet-2_qx24_check.py: ALL CHECKS PASSED")
    print(f"  part (i) max SL recasting error:         {err_i:.3e}")
    print(f"  part (ii) max SL recasting error:        {err_ii:.3e}")
    print(f"  part (iii) max shooting eigenvalue err:  {eig_err:.3e}")
    print(f"  part (iii) max orthogonality error:      {orth_err:.3e}")


if __name__ == "__main__":
    main()
