#!/usr/bin/env python3
"""Independent verification for IB/methods/sheet-4/qx10 (Symmetric tensors & invariants)."""

from fractions import Fraction
import math
import random


def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def mat_transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def mat_det3(A):
    return (
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
    )


def invariants_index(T):
    """Compute alpha = T_ii, beta = T_ik T_ki, gamma = T_ik T_km T_mi directly by summation."""
    zero = T[0][0] * 0
    alpha = sum((T[i][i] for i in range(3)), zero)
    beta = sum((T[i][k] * T[k][i] for i in range(3) for k in range(3)), zero)
    gamma = sum(
        (T[i][k] * T[k][m] * T[m][i] for i in range(3) for k in range(3) for m in range(3)),
        zero,
    )
    return alpha, beta, gamma


def char_poly_from_invariants(lam, alpha, beta, gamma):
    """Evaluate lambda^3 - alpha*lambda^2 + (1/2)(alpha^2-beta)*lambda - (1/6)(alpha^3-3*alpha*beta+2*gamma)."""
    if isinstance(alpha, Fraction):
        c2 = -alpha
        c1 = Fraction(1, 2) * (alpha * alpha - beta)
        c0 = -Fraction(1, 6) * (alpha**3 - 3 * alpha * beta + 2 * gamma)
    else:
        c2 = -alpha
        c1 = 0.5 * (alpha * alpha - beta)
        c0 = -(1.0 / 6.0) * (alpha**3 - 3.0 * alpha * beta + 2.0 * gamma)
    return lam**3 + c2 * lam**2 + c1 * lam + c0


def random_orthogonal(rng, det_sign=1):
    while True:
        cols = [[rng.gauss(0.0, 1.0) for _ in range(3)] for _ in range(3)]
        ortho = []
        ok = True
        for v in cols:
            w = v[:]
            for u in ortho:
                dot = sum(w[k] * u[k] for k in range(3))
                w = [w[k] - dot * u[k] for k in range(3)]
            norm = math.sqrt(sum(x * x for x in w))
            if norm < 1e-8:
                ok = False
                break
            ortho.append([x / norm for x in w])
        if not ok:
            continue
        L = [[ortho[j][i] for j in range(3)] for i in range(3)]
        if (mat_det3(L) > 0) != (det_sign > 0):
            L[2] = [-x for x in L[2]]
        return L


def main():
    rng = random.Random(20261127)

    # Part (i): invariance of alpha, beta, gamma for arbitrary (non-symmetric) tensors under SO(3) and O(3)
    max_da = max_db = max_dg = 0.0
    for trial in range(60):
        T = [[rng.uniform(-5.0, 5.0) for _ in range(3)] for _ in range(3)]
        det_sign = 1 if trial < 50 else -1
        L = random_orthogonal(rng, det_sign=det_sign)
        Tp = mat_mul(mat_mul(L, T), mat_transpose(L))

        a, b, g = invariants_index(T)
        ap, bp, gp = invariants_index(Tp)

        max_da = max(max_da, abs(ap - a))
        max_db = max(max_db, abs(bp - b))
        max_dg = max(max_dg, abs(gp - g))

    assert max_da < 1e-11
    assert max_db < 1e-10
    assert max_dg < 1e-9
    print(
        f"OK: Part (i) verified on 60 arbitrary non-symmetric tensors (50 SO(3), 10 O(3)\\SO(3)) — "
        f"max |da|={max_da:.2e}, |db|={max_db:.2e}, |dg|={max_dg:.2e}."
    )

    # Part (ii): exact rational checks for eigenvalue power sums and cubic equation
    test_triples = [
        (Fraction(1), Fraction(2), Fraction(3)),
        (Fraction(-4), Fraction(0), Fraction(7)),
        (Fraction(3), Fraction(3), Fraction(-6)),
        (Fraction(5), Fraction(5), Fraction(5)),
        (Fraction(1, 2), Fraction(-3, 4), Fraction(5, 7)),
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(-1, 3), Fraction(-2, 5), Fraction(-7, 11)),
    ]
    for l1, l2, l3 in test_triples:
        D = [
            [l1, Fraction(0), Fraction(0)],
            [Fraction(0), l2, Fraction(0)],
            [Fraction(0), Fraction(0), l3],
        ]
        a, b, g = invariants_index(D)
        assert a == l1 + l2 + l3
        assert b == l1**2 + l2**2 + l3**2
        assert g == l1**3 + l2**3 + l3**3

        # Exact coefficient check against (lambda - l1)(lambda - l2)(lambda - l3)
        e1 = l1 + l2 + l3
        e2 = l1 * l2 + l2 * l3 + l3 * l1
        e3 = l1 * l2 * l3
        assert a == e1
        assert Fraction(1, 2) * (a**2 - b) == e2
        assert Fraction(1, 6) * (a**3 - 3 * a * b + 2 * g) == e3
        assert e3 == mat_det3(D)

        for r in (l1, l2, l3):
            assert char_poly_from_invariants(r, a, b, g) == Fraction(0)
        for test_x in (Fraction(-2), Fraction(0), Fraction(1, 3), Fraction(5)):
            lhs = char_poly_from_invariants(test_x, a, b, g)
            rhs = (test_x - l1) * (test_x - l2) * (test_x - l3)
            assert lhs == rhs
    print("OK: Part (ii) exact rational identities verified on 7 eigenvalue triples (incl. repeated/zero).")

    # Part (ii): rotated dense symmetric tensors, det(lambda I - T), and Cayley-Hamilton check
    max_root_err = 0.0
    max_det_err = 0.0
    max_ch_err = 0.0
    for _ in range(30):
        eigs = [rng.uniform(-4.0, 4.0) for _ in range(3)]
        D = [[eigs[i] if i == j else 0.0 for j in range(3)] for i in range(3)]
        Q = random_orthogonal(rng, det_sign=1)
        T = mat_mul(mat_mul(Q, D), mat_transpose(Q))

        # Check symmetry
        for i in range(3):
            for j in range(3):
                assert abs(T[i][j] - T[j][i]) < 1e-13

        a, b, g = invariants_index(T)
        assert abs(a - sum(eigs)) < 1e-12
        assert abs(b - sum(x**2 for x in eigs)) < 1e-11
        assert abs(g - sum(x**3 for x in eigs)) < 1e-10

        for r in eigs:
            max_root_err = max(max_root_err, abs(char_poly_from_invariants(r, a, b, g)))

        for x in (-3.5, 0.0, 1.25, 4.7):
            M = [[(x if i == j else 0.0) - T[i][j] for j in range(3)] for i in range(3)]
            max_det_err = max(max_det_err, abs(mat_det3(M) - char_poly_from_invariants(x, a, b, g)))

        # Cayley-Hamilton: T^3 - a T^2 + e2 T - e3 I = 0
        T2 = mat_mul(T, T)
        T3 = mat_mul(T2, T)
        e2 = 0.5 * (a * a - b)
        e3 = (1.0 / 6.0) * (a**3 - 3.0 * a * b + 2.0 * g)
        for i in range(3):
            for j in range(3):
                val = T3[i][j] - a * T2[i][j] + e2 * T[i][j] - (e3 if i == j else 0.0)
                max_ch_err = max(max_ch_err, abs(val))

    assert max_root_err < 1e-10
    assert max_det_err < 1e-10
    assert max_ch_err < 1e-10
    print(
        f"OK: 30 dense rotated symmetric tensors verified — max root err={max_root_err:.2e}, "
        f"max det err={max_det_err:.2e}, max Cayley-Hamilton err={max_ch_err:.2e}."
    )


if __name__ == "__main__":
    main()
