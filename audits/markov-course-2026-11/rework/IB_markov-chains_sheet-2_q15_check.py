#!/usr/bin/env python3
"""Exact rational verification for IB/markov-chains/sheet-2/q15 (Birth-and-death chain on {0,...,b})."""

from fractions import Fraction
import itertools


def mat_mul(A, B):
    n = len(A)
    C = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            aik = A[i][k]
            if aik != 0:
                for j in range(n):
                    C[i][j] += aik * B[k][j]
    return C


def mat_pow(A, n):
    size = len(A)
    res = [[Fraction(1 if i == j else 0) for j in range(size)] for i in range(size)]
    base = [row[:] for row in A]
    e = n
    while e > 0:
        if e & 1:
            res = mat_mul(res, base)
        base = mat_mul(base, base)
        e >>= 1
    return res


def rref_nullspace_exact(M):
    A = [row[:] for row in M]
    rows = len(A)
    cols = len(A[0])
    pivot_cols = []
    r = 0
    for c in range(cols):
        pivot = None
        for i in range(r, rows):
            if A[i][c] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        pv = A[r][c]
        for j in range(c, cols):
            A[r][j] /= pv
        for i in range(rows):
            if i != r and A[i][c] != 0:
                factor = A[i][c]
                for j in range(c, cols):
                    A[i][j] -= factor * A[r][j]
        pivot_cols.append(c)
        r += 1
    free_cols = [c for c in range(cols) if c not in pivot_cols]
    basis = []
    for fc in free_cols:
        vec = [Fraction(0)] * cols
        vec[fc] = Fraction(1)
        for idx, pc in enumerate(pivot_cols):
            vec[pc] = -A[idx][fc]
        basis.append(vec)
    return basis


def build_chain(b, lambdas, mu_b):
    """Build P on {0,...,b}.
    lambdas has length b: lambda_0, ..., lambda_{b-1} in (0,1).
    mu_i = 1 - lambda_i for 1 <= i < b, and mu_b in (0,1).
    """
    assert b >= 1
    assert len(lambdas) == b
    assert all(0 < x < 1 for x in lambdas)
    assert 0 < mu_b < 1
    mus = [None] + [1 - lambdas[i] for i in range(1, b)] + [mu_b]
    assert all(0 < mus[i] < 1 for i in range(1, b + 1))

    n = b + 1
    P = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    P[0][0] = 1 - lambdas[0]
    P[b][b] = 1 - mu_b
    for i in range(b):
        P[i][i + 1] = lambdas[i]
        P[i + 1][i] = mus[i + 1]

    # Verify stochastic matrix
    for i in range(n):
        assert sum(P[i]) == 1
        assert all(0 <= x <= 1 for x in P[i])
        if 1 <= i < b:
            assert P[i][i] == 0
    return P, mus


def formula_pi(b, lambdas, mus):
    w = [Fraction(1)]
    for i in range(b):
        w.append(w[-1] * lambdas[i] / mus[i + 1])
    Z = sum(w)
    return [x / Z for x in w]


def verify_case(name, b, lambdas, mu_b, check_paths=True):
    P, mus = build_chain(b, lambdas, mu_b)
    n = b + 1

    # 1. Irreducibility & aperiodicity: P^b has all entries strictly positive
    Pb = mat_pow(P, b)
    assert all(Pb[i][j] > 0 for i in range(n) for j in range(n))

    # 2. Exact closed-form invariant distribution
    pi = formula_pi(b, lambdas, mus)
    assert sum(pi) == 1
    assert all(x > 0 for x in pi)

    # 3. Independently verify uniqueness via exact nullspace of P^T - I
    PT_I = [[P[j][i] - (1 if i == j else 0) for j in range(n)] for i in range(n)]
    ns = rref_nullspace_exact(PT_I)
    assert len(ns) == 1, f"Expected unique invariant distribution, got dim {len(ns)}"
    v = ns[0]
    s = sum(v)
    pi_solved = [x / s for x in v]
    assert pi_solved == pi, f"Formula mismatch: {pi} != {pi_solved}"

    # 4. Detailed balance for all i, j
    for i in range(n):
        for j in range(n):
            assert pi[i] * P[i][j] == pi[j] * P[j][i]

    # 5. Full joint path time-reversal check for all paths up to N=3 (or 4 for small b)
    if check_paths:
        max_N = 4 if b <= 2 else 3
        for N in range(1, max_N + 1):
            for path in itertools.product(range(n), repeat=N + 1):
                pf = pi[path[0]]
                for t in range(N):
                    pf *= P[path[t]][path[t + 1]]
                rev = path[::-1]
                pr = pi[rev[0]]
                for t in range(N):
                    pr *= P[rev[t]][rev[t + 1]]
                assert pf == pr
        path_msg = f"all paths N<={max_N} reversible"
    else:
        path_msg = "detailed balance verified (full path check skipped for large b)"

    print(
        f"[PASS] {name} (b={b}): P^{b} > 0 (irreducible & aperiodic), unique pi matches formula, "
        f"detailed balance holds, {path_msg}."
    )


def main():
    # Edge case b = 1 (interior condition 1 <= i < b is vacuous)
    verify_case("Boundary b=1", 1, [Fraction(2, 5)], Fraction(3, 7), check_paths=True)
    # b = 2 asymmetric
    verify_case(
        "Asymmetric b=2",
        2,
        [Fraction(1, 3), Fraction(3, 4)],
        Fraction(2, 5),
        check_paths=True,
    )
    # b = 3 heterogeneous
    verify_case(
        "Heterogeneous b=3",
        3,
        [Fraction(1, 2), Fraction(2, 5), Fraction(4, 7)],
        Fraction(1, 6),
        check_paths=True,
    )
    # b = 4 strongly biased
    verify_case(
        "Biased b=4",
        4,
        [Fraction(4, 5), Fraction(3, 4), Fraction(2, 3), Fraction(5, 6)],
        Fraction(1, 7),
        check_paths=True,
    )
    # b = 7 larger chain
    verify_case(
        "Larger b=7",
        7,
        [
            Fraction(1, 4),
            Fraction(2, 3),
            Fraction(1, 5),
            Fraction(3, 7),
            Fraction(5, 8),
            Fraction(2, 9),
            Fraction(4, 5),
        ],
        Fraction(3, 11),
        check_paths=False,
    )
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
