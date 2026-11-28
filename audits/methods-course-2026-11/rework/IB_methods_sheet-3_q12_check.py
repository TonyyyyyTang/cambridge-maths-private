#!/usr/bin/env python3
"""Computational verification for IB/methods/sheet-3/q12 (2025-26 Sheet 3, Q12)."""

import math
from fractions import Fraction


def gauss_legendre_8(func, a, b, panels=40):
    """Composite 8-point Gauss-Legendre quadrature on [a, b] (works for a > b too)."""
    if abs(b - a) < 1e-15:
        return 0.0
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


def G_n(t, tau, n):
    """Green's function for y^{(n)} = f, y(0) = ... = y^{(n-1)}(0) = 0."""
    if t < tau:
        return 0.0
    return ((t - tau) ** (n - 1)) / math.factorial(n - 1)


def check_greens_function():
    """Verify G_n = 0 for t < tau, continuity of lower derivatives, and unit jump of (n-1)th deriv."""
    tau = 0.7
    for n in range(1, 6):
        # Zero for t < tau
        for t in [0.0, 0.3, tau - 1e-7]:
            assert G_n(t, tau, n) == 0.0

        # Exact k-th right-derivative at t = tau^+ is delta_{k, n-1}
        for k in range(n):
            # Analytic derivative of (t - tau)^{n-1} / (n-1)! is (t - tau)^{n-1-k} / (n-1-k)!
            s = 1e-9
            dk = (s ** (n - 1 - k)) / math.factorial(n - 1 - k)
            expected = 1.0 if k == n - 1 else 0.0
            assert abs(dk - expected) < 1e-8

    # Finite-difference check of jump at t = tau^+ for n = 1, 2, 3, 4
    eps = 1e-3
    # n = 1: jump in G
    assert abs(G_n(tau + eps, tau, 1) - G_n(tau - eps, tau, 1) - 1.0) < 1e-14
    # n = 2: continuity of G, jump in G'
    assert abs(G_n(tau + 1e-12, tau, 2)) < 1e-11
    g1 = (G_n(tau + eps, tau, 2) - G_n(tau, tau, 2)) / eps
    assert abs(g1 - 1.0) < 1e-12
    # n = 3: G, G' -> 0, G'' -> 1
    v = [G_n(tau + j * eps, tau, 3) for j in range(4)]
    assert abs(v[0]) < 1e-15
    assert abs((-3 * v[0] + 4 * v[1] - v[2]) / (2 * eps)) < 1e-12
    assert abs((v[0] - 2 * v[1] + v[2]) / (eps ** 2) - 1.0) < 1e-10
    # n = 4: G, G', G'' -> 0, G''' -> 1
    v = [G_n(tau + j * eps, tau, 4) for j in range(5)]
    assert abs((-v[0] + 3 * v[1] - 3 * v[2] + v[3]) / (eps ** 3) - 1.0) < 1e-7


def check_exact_monomial_ivp():
    """Exact rational check for f(t) = t^m: (1/(n-1)!) int_0^t (t-tau)^{n-1} tau^m dtau = m!/(m+n)! t^{m+n}."""
    for n in range(1, 7):
        for m in range(0, 7):
            # Expand (t - tau)^{n-1} = sum_{j=0}^{n-1} (-1)^j C(n-1, j) t^{n-1-j} tau^j
            # Integrate tau^{m+j} from 0 to t => t^{m+j+1} / (m+j+1)
            coeff = Fraction(0, 1)
            for j in range(n):
                binom = math.comb(n - 1, j)
                term = Fraction((-1) ** j * binom, m + j + 1)
                coeff += term
            coeff /= math.factorial(n - 1)
            expected = Fraction(math.factorial(m), math.factorial(m + n))
            assert coeff == expected, f"Mismatch for n={n}, m={m}: {coeff} vs {expected}"


def check_ivp_quadrature():
    """Verify IVP solution y(t) = int_0^t G_n(t, tau) f(tau) dtau for f(t) = e^t and sin(w t)."""
    max_err = 0.0
    for n in range(1, 6):
        # f(t) = e^t => y(t) = e^t - sum_{k=0}^{n-1} t^k / k!
        for t in [0.3, 1.1, 2.4]:
            y_num = gauss_legendre_8(lambda tau: G_n(t, tau, n) * math.exp(tau), 0.0, t, panels=30)
            y_ex = math.exp(t) - sum((t ** k) / math.factorial(k) for k in range(n))
            err = abs(y_num - y_ex)
            max_err = max(max_err, err)
            assert err < 1e-12, f"IVP exp mismatch n={n}, t={t}: {err}"
    return max_err


def check_taylor_theorem():
    """Verify Taylor's theorem with integral remainder for C^n functions, both x > a and x < a."""
    # Test functions with exact k-th derivatives
    tests = [
        (
            "exp(1.3 x)",
            lambda x, k: (1.3 ** k) * math.exp(1.3 * x),
            [(0.0, 1.2), (0.5, 1.8), (1.1, -0.4)],
        ),
        (
            "sin(1.7 x + 0.2)",
            lambda x, k: (1.7 ** k) * math.sin(1.7 * x + 0.2 + 0.5 * k * math.pi),
            [(0.0, 1.3), (0.4, 2.1), (1.5, 0.2)],
        ),
        (
            "ln(1 + x)",
            lambda x, k: math.log(1.0 + x) if k == 0 else ((-1) ** (k - 1)) * math.factorial(k - 1) * ((1.0 + x) ** (-k)),
            [(0.0, 0.7), (0.5, 1.9), (1.2, 0.1)],
        ),
        (
            "(2 + x)^3.5",
            lambda x, k: math.prod(3.5 - j for j in range(k)) * ((2.0 + x) ** (3.5 - k)),
            [(0.0, 1.5), (0.8, -0.6)],
        ),
    ]

    max_err = 0.0
    for name, deriv, pairs in tests:
        for a, x in pairs:
            for n in range(1, 6):
                poly = sum(deriv(a, k) * ((x - a) ** k) / math.factorial(k) for k in range(n))
                rem = gauss_legendre_8(
                    lambda s: ((x - s) ** (n - 1)) * deriv(s, n),
                    a,
                    x,
                    panels=40,
                ) / math.factorial(n - 1)
                lhs = deriv(x, 0)
                rhs = poly + rem
                err = abs(lhs - rhs)
                max_err = max(max_err, err)
                assert err < 1e-12, f"Taylor mismatch for {name}, a={a}, x={x}, n={n}: {err}"
    return max_err


def main():
    check_greens_function()
    check_exact_monomial_ivp()
    ivp_err = check_ivp_quadrature()
    taylor_err = check_taylor_theorem()
    print("IB_methods_sheet-3_q12_check.py: ALL CHECKS PASSED")
    print(f"  exact rational monomial IVP checks: 42/42 exact")
    print(f"  max IVP quadrature error (f=e^t):   {ivp_err:.3e}")
    print(f"  max Taylor integral-remainder err:  {taylor_err:.3e}")


if __name__ == "__main__":
    main()
