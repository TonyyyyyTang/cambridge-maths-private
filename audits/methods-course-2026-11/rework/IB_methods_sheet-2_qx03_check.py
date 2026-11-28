#!/usr/bin/env python3
"""Verification script for IB/methods/sheet-2/qx03 (2016-17 Sheet 2, Q3).

Checks:
1. nabla^2 (1 / |r - r'|) = 0 for r != r' in R^3 via 4th-order finite differences.
2. Exact rational verification (using Fraction) of Legendre polynomial orthogonality
   and the exact coefficients c_l for the split-hemisphere Dirichlet problem:
   c_{2k} = 0, c_{2k+1} = V * (-1)^k * (4k+3) * (2k)! / (2^{2k+1} k! (k+1)!).
3. Numerical verification of:
   - Spherical-polar Laplace equation nabla^2 Phi = 0 inside r < a,
   - Exact closed form on the z-axis (theta = 0, pi) and vanishing on equator (theta = pi/2),
   - Independent 3D Poisson integral formula at general interior points (r, theta),
   - Convergence to +/- V as r -> a^-.
"""

from fractions import Fraction
import math


def d2_central_4th(func, t, h=1e-3):
    return (-func(t + 2 * h) + 16 * func(t + h) - 30 * func(t)
            + 16 * func(t - h) - func(t - 2 * h)) / (12.0 * h * h)


def d1_central_4th(func, t, h=1e-3):
    return (-func(t + 2 * h) + 8 * func(t + h)
            - 8 * func(t - h) + func(t - 2 * h)) / (12.0 * h)


def poly_mul(p, q):
    res = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            res[i + j] += a * b
    return res


def poly_int(p, a, b):
    """Exact integral of polynomial p(x) = sum p[i] x^i from a to b."""
    def eval_anti(t):
        total = Fraction(0)
        power = Fraction(t)
        for i, c in enumerate(p):
            total += c * power / (i + 1)
            power *= t
        return total
    return eval_anti(b) - eval_anti(a)


def poly_eval(p, x):
    total = Fraction(0)
    power = Fraction(1)
    for c in p:
        total += c * power
        power *= x
    return total


def legendre_polys_exact(max_l):
    """Return list of Legendre polynomials P_0 .. P_{max_l} as Fraction coefficient lists."""
    P = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    for l in range(1, max_l):
        # (l+1) P_{l+1}(x) = (2l+1) x P_l(x) - l P_{l-1}(x)
        prev = P[l - 1]
        curr = P[l]
        nxt = [Fraction(0)] * (l + 2)
        for i, c in enumerate(curr):
            nxt[i + 1] += Fraction(2 * l + 1, l + 1) * c
        for i, c in enumerate(prev):
            nxt[i] -= Fraction(l, l + 1) * c
        P.append(nxt)
    return P


def c_odd_formula_exact(k):
    """Exact rational coefficient c_{2k+1} / V."""
    num = ((-1) ** k) * (4 * k + 3) * math.factorial(2 * k)
    den = (2 ** (2 * k + 1)) * math.factorial(k) * math.factorial(k + 1)
    return Fraction(num, den)


def double_factorial(n):
    if n <= 0:
        return 1
    res = 1
    while n > 1:
        res *= n
        n -= 2
    return res


def legendre_eval_float(n, x):
    """Evaluate P_n(x) in floating point via recurrence."""
    if n == 0:
        return 1.0
    if n == 1:
        return x
    p0, p1 = 1.0, x
    for l in range(1, n):
        p0, p1 = p1, ((2 * l + 1) * x * p1 - l * p0) / (l + 1)
    return p1


def c_odd_float(k):
    # Recurrence: c_{2k+1}/V = (4k+3)/(2(k+1)) * P_{2k}(0),
    # with P_0(0) = 1, P_{2k}(0) = -((2k-1)/(2k)) P_{2k-2}(0).
    p2k0 = 1.0
    for m in range(1, k + 1):
        p2k0 *= -(2 * m - 1) / (2 * m)
    return ((4 * k + 3) / (2 * (k + 1))) * p2k0


def phi_series(r, theta, a=1.0, V=1.0, n_terms=60):
    t = r / a
    x = math.cos(theta)
    # Build P_0(x), P_1(x), ... on the fly
    total = 0.0
    p_prev = 1.0
    p_curr = x
    p2k0 = 1.0
    t_pow = t
    for k in range(n_terms):
        # At step k, p_curr = P_{2k+1}(x), p2k0 = P_{2k}(0)
        ck = ((4 * k + 3) / (2 * (k + 1))) * p2k0
        total += ck * t_pow * p_curr
        # Advance P_{2k+1} -> P_{2k+2} -> P_{2k+3}
        l1 = 2 * k + 1
        p_even = ((2 * l1 + 1) * x * p_curr - l1 * p_prev) / (l1 + 1)
        l2 = 2 * k + 2
        p_next = ((2 * l2 + 1) * x * p_even - l2 * p_curr) / (l2 + 1)
        p_prev, p_curr = p_even, p_next
        p2k0 *= -(2 * k + 1) / (2 * k + 2)
        t_pow *= t * t
    return V * total


def simpson_real(func, a, b, n):
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    total = func(a) + func(b)
    for j in range(1, n):
        total += (4.0 if j % 2 == 1 else 2.0) * func(a + j * h)
    return total * (h / 3.0)


def phi_poisson_3d(r, theta, a=1.0, V=1.0, n_mu=120, n_phi=120):
    """Independent 3D Poisson integral for the ball of radius a with +/-V hemispheres.

    Phi(r, theta) = (a(a^2 - r^2) / (4pi)) * int_{S^2} Phi(a, theta') / |r - r'|^3 dOmega'
    Exploiting symmetry phi' -> -phi' (integrate phi' in [0, pi] and double),
    and mu' = cos(theta') in [0, 1] minus [-1, 0]:
    |r - r'|^2 = r^2 + a^2 - 2 a r (cos theta cos theta' + sin theta sin theta' cos phi').
    """
    ct, st = math.cos(theta), math.sin(theta)
    pref = a * (a * a - r * r) / (2.0 * math.pi)  # doubled for phi' in [0, pi]

    def integrand_th(th_p):
        mu = math.cos(th_p)
        s_p = math.sin(th_p)
        # Upper hemisphere (+mu) minus lower hemisphere (-mu), with Jacobian sin(th_p)
        def phi_int(phi_p):
            c_p = math.cos(phi_p)
            d_plus = r * r + a * a - 2.0 * a * r * (ct * mu + st * s_p * c_p)
            d_minus = r * r + a * a - 2.0 * a * r * (-ct * mu + st * s_p * c_p)
            return (d_plus ** -1.5) - (d_minus ** -1.5)
        return s_p * simpson_real(phi_int, 0.0, math.pi, n_phi)

    return V * pref * simpson_real(integrand_th, 0.0, math.pi / 2.0, n_mu)


def main():
    # 1. Part 1: Check nabla^2 (1 / |r - r'|) = 0 for r != r' in R^3
    max_part1_err = 0.0
    test_pairs = [
        ((1.2, -0.5, 0.8), (0.0, 0.0, 0.0)),
        ((2.5, 1.1, -3.2), (-0.7, 2.3, 0.4)),
        ((0.3, 0.4, 0.5), (0.1, -0.2, 0.9)),
    ]
    for r, rp in test_pairs:
        def inv_dist(x, y, z):
            dx, dy, dz = x - rp[0], y - rp[1], z - rp[2]
            return 1.0 / math.sqrt(dx * dx + dy * dy + dz * dz)
        lap = (d2_central_4th(lambda x: inv_dist(x, r[1], r[2]), r[0], h=3e-3)
               + d2_central_4th(lambda y: inv_dist(r[0], y, r[2]), r[1], h=3e-3)
               + d2_central_4th(lambda z: inv_dist(r[0], r[1], z), r[2], h=3e-3))
        max_part1_err = max(max_part1_err, abs(lap))
        assert abs(lap) < 1e-10, f"Part 1 Laplacian residual too large: {lap}"

    # 2. Part 2: Exact rational check of orthogonality and coefficients c_l (l = 0..15)
    max_l = 15
    P = legendre_polys_exact(max_l)
    for l in range(max_l + 1):
        for m in range(max_l + 1):
            val = poly_int(poly_mul(P[l], P[m]), -1, 1)
            expected = Fraction(2, 2 * l + 1) if l == m else Fraction(0)
            assert val == expected, f"Orthogonality failed for ({l}, {m}): {val} != {expected}"

    for l in range(max_l + 1):
        c_l = Fraction(2 * l + 1, 2) * (poly_int(P[l], 0, 1) - poly_int(P[l], -1, 0))
        if l % 2 == 0:
            assert c_l == 0, f"Even coefficient c_{l} != 0: {c_l}"
        else:
            k = (l - 1) // 2
            c_formula = c_odd_formula_exact(k)
            c_double_fact = Fraction(((-1) ** k) * (4 * k + 3) * double_factorial(2 * k - 1),
                                     double_factorial(2 * k + 2))
            assert c_l == c_formula == c_double_fact, f"Mismatch at l={l}: {c_l}, {c_formula}, {c_double_fact}"
            # Also check P_{l-1}(0) - P_{l+1}(0) when l+1 <= max_l
            if l + 1 <= max_l:
                assert c_l == poly_eval(P[l - 1], 0) - poly_eval(P[l + 1], 0)

    # Check first 4 explicit odd coefficients
    assert [c_odd_formula_exact(k) for k in range(4)] == [
        Fraction(3, 2), Fraction(-7, 8), Fraction(11, 16), Fraction(-75, 128)
    ]

    # 3. Spherical-polar Laplace equation check inside r < a
    a, V = 2.5, 3.7
    max_pde_err = 0.0
    for r in (0.4, 1.0, 1.6):
        for theta in (0.35, 0.9, 1.4, 2.2, 2.75):
            # nabla^2 Phi = Phi_rr + (2/r) Phi_r + (1/r^2)(Phi_tt + cot(theta) Phi_t)
            phi_rr = d2_central_4th(lambda u: phi_series(u, theta, a, V), r, h=2e-3)
            phi_r = d1_central_4th(lambda u: phi_series(u, theta, a, V), r, h=2e-3)
            phi_tt = d2_central_4th(lambda th: phi_series(r, th, a, V), theta, h=2e-3)
            phi_t = d1_central_4th(lambda th: phi_series(r, th, a, V), theta, h=2e-3)
            lap = phi_rr + (2.0 / r) * phi_r + (phi_tt + (math.cos(theta) / math.sin(theta)) * phi_t) / (r * r)
            max_pde_err = max(max_pde_err, abs(lap))
            assert abs(lap) < 5e-9, f"Spherical PDE residual too large at ({r}, {theta}): {lap}"

    # 4. Exact closed-form check on z-axis (theta = 0, pi) and equator (theta = pi/2)
    max_axis_err = 0.0
    for t in (0.05, 0.25, 0.5, 0.75, 0.9):
        r = t * a
        exact_north = V * (1.0 / t - (1.0 - t * t) / (t * math.sqrt(1.0 + t * t)))
        ser_north = phi_series(r, 0.0, a, V, n_terms=150)
        ser_south = phi_series(r, math.pi, a, V, n_terms=150)
        ser_eq = phi_series(r, math.pi / 2.0, a, V, n_terms=150)
        max_axis_err = max(max_axis_err, abs(ser_north - exact_north), abs(ser_south + exact_north))
        assert abs(ser_north - exact_north) < 1e-12
        assert abs(ser_south + exact_north) < 1e-12
        assert abs(ser_eq) < 1e-14  # cos(pi/2) ~ 6.12e-17 in float64

    # 5. Independent 3D Poisson integral check at general off-axis interior points
    max_poisson_err = 0.0
    for r_frac, theta in ((0.2, 0.4), (0.35, 1.1), (0.5, 0.7), (0.5, 2.3), (0.6, 1.9)):
        r = r_frac * a
        val_ser = phi_series(r, theta, a, V, n_terms=80)
        val_poi = phi_poisson_3d(r, theta, a, V, n_mu=600, n_phi=100)
        err = abs(val_ser - val_poi)
        max_poisson_err = max(max_poisson_err, err)
        assert err < 1e-10, f"Poisson mismatch at ({r_frac}, {theta}): {val_ser} vs {val_poi}"

    # 6. Boundary convergence check as r -> a^-
    val_upper_bc = phi_series(0.999 * a, 0.6, a, V, n_terms=4000)
    val_lower_bc = phi_series(0.999 * a, 2.5, a, V, n_terms=4000)
    assert abs(val_upper_bc - V) < 1e-2 and abs(val_lower_bc + V) < 1e-2

    print("IB_methods_sheet-2_qx03_check.py: ALL CHECKS PASSED")
    print(f"  Part 1 max 4th-order FD |nabla^2 (1/|r-r'|)|:       {max_part1_err:.3e}")
    print(f"  Part 2 exact rational orthogonality & c_l (l<=15):  EXACT")
    print(f"  Part 2 first 4 odd coefficients c_{{2k+1}}/V:        3/2, -7/8, 11/16, -75/128")
    print(f"  Part 2 max spherical-polar Laplace residual:        {max_pde_err:.3e}")
    print(f"  Part 2 max z-axis series-vs-closed-form error:      {max_axis_err:.3e}")
    print(f"  Part 2 max off-axis series-vs-3D-Poisson error:     {max_poisson_err:.3e}")
    print(f"  Part 2 near-boundary r=0.999a (theta=0.6, 2.5):     {val_upper_bc:.5f}, {val_lower_bc:.5f} (target +/-{V})")


if __name__ == "__main__":
    main()
