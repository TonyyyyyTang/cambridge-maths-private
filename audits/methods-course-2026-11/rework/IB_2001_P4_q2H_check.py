#!/usr/bin/env python3
"""Independent verification for IB/2001/P4/q2H (Sturm-Liouville form & orthogonality of R_n = P_n')."""

from fractions import Fraction
import math


# --- Exact rational polynomial arithmetic: list of Fraction coeffs [c0, c1, ..., cd] ---

def poly_trim(p):
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q if q else [Fraction(0)]


def poly_add(p, q):
    n = max(len(p), len(q))
    res = [
        (p[i] if i < len(p) else Fraction(0)) + (q[i] if i < len(q) else Fraction(0))
        for i in range(n)
    ]
    return poly_trim(res)


def poly_sub(p, q):
    n = max(len(p), len(q))
    res = [
        (p[i] if i < len(p) else Fraction(0)) - (q[i] if i < len(q) else Fraction(0))
        for i in range(n)
    ]
    return poly_trim(res)


def poly_scale(c, p):
    return poly_trim([c * x for x in p])


def poly_mul(p, q):
    res = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            res[i + j] += a * b
    return poly_trim(res)


def poly_deriv(p):
    if len(p) <= 1:
        return [Fraction(0)]
    return poly_trim([Fraction(i) * p[i] for i in range(1, len(p))])


def poly_eval(p, x):
    acc = Fraction(0) if isinstance(x, Fraction) else 0.0
    for c in reversed(p):
        acc = acc * x + (c if isinstance(x, Fraction) else float(c))
    return acc


def poly_integrate_m1_to_1(p):
    """Exact integral of polynomial p(x) from -1 to 1."""
    total = Fraction(0)
    for k, c in enumerate(p):
        if k % 2 == 0:
            total += c * Fraction(2, k + 1)
    return total


def legendre_polys(max_n):
    """Generate exact Legendre polynomials P_0, ..., P_{max_n} via Bonnet's recurrence."""
    P = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    x_poly = [Fraction(0), Fraction(1)]
    for n in range(1, max_n):
        # (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
        term1 = poly_scale(Fraction(2 * n + 1, n + 1), poly_mul(x_poly, P[n]))
        term2 = poly_scale(Fraction(n, n + 1), P[n - 1])
        P.append(poly_sub(term1, term2))
    return P


def gauss_legendre(n):
    nodes = []
    weights = []
    for k in range(1, n + 1):
        x = math.cos(math.pi * (k - 0.25) / (n + 0.5))
        for _ in range(30):
            p0, p1 = 1.0, x
            for j in range(2, n + 1):
                p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-15:
                break
        p0, p1 = 1.0, x
        for j in range(2, n + 1):
            p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
        dp = n * (x * p1 - p0) / (x * x - 1.0)
        w = 2.0 / ((1.0 - x * x) * dp * dp)
        nodes.append(x)
        weights.append(w)
    return nodes, weights


def gl_integrate(f, nodes, weights):
    return sum(w * f(x) for x, w in zip(nodes, weights))


def main():
    max_n = 10
    P = legendre_polys(max_n)
    R = [poly_deriv(p) for p in P]

    one_minus_x2 = [Fraction(1), Fraction(0), Fraction(-1)]        # 1 - x^2
    one_minus_x2_sq = poly_mul(one_minus_x2, one_minus_x2)          # (1 - x^2)^2
    minus_2x = [Fraction(0), Fraction(-2)]                          # -2x
    minus_4x = [Fraction(0), Fraction(-4)]                          # -4x
    zero_poly = [Fraction(0)]

    # 1. Verify Legendre ODE, differentiated ODE for R_n, and Sturm-Liouville form identically
    for n in range(max_n + 1):
        # Legendre ODE: (1-x^2) P_n'' - 2x P_n' + n(n+1) P_n = 0
        Pn = P[n]
        Pnp = poly_deriv(Pn)
        Pnpp = poly_deriv(Pnp)
        res_leg = poly_add(
            poly_add(poly_mul(one_minus_x2, Pnpp), poly_mul(minus_2x, Pnp)),
            poly_scale(Fraction(n * (n + 1)), Pn),
        )
        assert res_leg == zero_poly, f"Legendre ODE failed for n={n}"

        # Differentiated ODE for R_n = P_n': (1-x^2) R_n'' - 4x R_n' + (n-1)(n+2) R_n = 0
        Rn = R[n]
        Rnp = poly_deriv(Rn)
        Rnpp = poly_deriv(Rnp)
        lam_n = Fraction((n - 1) * (n + 2))
        res_diff = poly_add(
            poly_add(poly_mul(one_minus_x2, Rnpp), poly_mul(minus_4x, Rnp)),
            poly_scale(lam_n, Rn),
        )
        assert res_diff == zero_poly, f"Differentiated ODE failed for n={n}"

        # Sturm-Liouville form: d/dx[(1-x^2)^2 R_n'] + (n-1)(n+2)(1-x^2) R_n = 0
        sl_term1 = poly_deriv(poly_mul(one_minus_x2_sq, Rnp))
        sl_term2 = poly_scale(lam_n, poly_mul(one_minus_x2, Rn))
        res_sl = poly_add(sl_term1, sl_term2)
        assert res_sl == zero_poly, f"Sturm-Liouville equation failed for n={n}"

        # Check R_0 == 0 and deg(R_n) == n-1 for n >= 1
        if n == 0:
            assert Rn == zero_poly
        else:
            assert len(Rn) - 1 == n - 1 and Rn[-1] != 0

    print(
        f"OK: Legendre ODE, differentiated ODE, and Sturm-Liouville form "
        f"d/dx[(1-x^2)^2 R_n'] + (n-1)(n+2)(1-x^2) R_n = 0 hold identically for n=0..{max_n}."
    )

    # 2. Check distinctness of eigenvalues and exact orthogonality + normalisation for 0 <= n, m <= 10
    for n in range(max_n + 1):
        for m in range(max_n + 1):
            lam_n = (n - 1) * (n + 2)
            lam_m = (m - 1) * (m + 2)
            if n != m:
                assert lam_n != lam_m
                assert lam_n - lam_m == (n - m) * (n + m + 1)

            # Boundary term [(1-x^2)^2 (R_m R_n' - R_n R_m')]_{-1}^{1}
            bt_poly = poly_mul(
                one_minus_x2_sq,
                poly_sub(poly_mul(R[m], poly_deriv(R[n])), poly_mul(R[n], poly_deriv(R[m]))),
            )
            bt_val = poly_eval(bt_poly, Fraction(1)) - poly_eval(bt_poly, Fraction(-1))
            assert bt_val == Fraction(0)

            # Exact orthogonality integral int_{-1}^1 (1-x^2) R_n(x) R_m(x) dx
            integrand = poly_mul(one_minus_x2, poly_mul(R[n], R[m]))
            val = poly_integrate_m1_to_1(integrand)
            if n != m:
                assert val == Fraction(0), f"Orthogonality failed for n={n}, m={m}: got {val}"
            else:
                expected = Fraction(2 * n * (n + 1), 2 * n + 1)
                assert val == expected, f"Norm failed for n={n}: got {val}, expected {expected}"

    print(
        f"OK: Exact rational integrals int_{{-1}}^1 (1-x^2) R_n(x) R_m(x) dx verified for all "
        f"0 <= n, m <= {max_n}: 0 for n != m, and 2n(n+1)/(2n+1) for n = m."
    )

    # 3. Independent numerical quadrature cross-check (16-point Gauss-Legendre) for all 0 <= n, m <= 10
    nodes, weights = gauss_legendre(16)
    max_num_err = 0.0
    for n in range(max_n + 1):
        for m in range(max_n + 1):
            num_val = gl_integrate(
                lambda x, Rn=R[n], Rm=R[m]: (1.0 - x * x) * poly_eval(Rn, x) * poly_eval(Rm, x),
                nodes,
                weights,
            )
            expected = 0.0 if n != m else (2.0 * n * (n + 1)) / (2.0 * n + 1)
            max_num_err = max(max_num_err, abs(num_val - expected))
    assert max_num_err < 1e-12
    print(
        f"OK: Independent 16-point Gauss-Legendre quadrature (0 <= n,m <= {max_n}) passed: "
        f"max err={max_num_err:.2e}."
    )


if __name__ == "__main__":
    main()
