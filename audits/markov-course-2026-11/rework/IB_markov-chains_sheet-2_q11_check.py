#!/usr/bin/env python3
"""Exact rational verification for IB/markov-chains/sheet-2/q11 (Sum of fair die throws mod 13)."""

import cmath
from fractions import Fraction
import itertools
import math


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


def matrix_rank_exact(M):
    A = [row[:] for row in M]
    rows = len(A)
    cols = len(A[0])
    rank = 0
    for c in range(cols):
        pivot = None
        for r in range(rank, rows):
            if A[r][c] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        A[rank], A[pivot] = A[pivot], A[rank]
        pv = A[rank][c]
        for j in range(c, cols):
            A[rank][j] /= pv
        for r in range(rows):
            if r != rank and A[r][c] != 0:
                factor = A[r][c]
                for j in range(c, cols):
                    A[r][j] -= factor * A[rank][j]
        rank += 1
    return rank


def main():
    m = 13
    P = [[Fraction(0) for _ in range(m)] for _ in range(m)]
    for i in range(m):
        for r in range(1, 7):
            P[i][(i + r) % m] += Fraction(1, 6)

    # 1. Verify doubly stochastic
    for i in range(m):
        assert sum(P[i]) == 1
        assert sum(P[r][i] for r in range(m)) == 1
    print("[PASS] 13x13 transition matrix P is doubly stochastic.")

    # 2. Verify P^1, P^2, P^3, P^4 diagonal entries and positivity
    P2 = mat_pow(P, 2)
    P3 = mat_pow(P, 3)
    P4 = mat_pow(P, 4)
    assert P[0][0] == 0
    assert P2[0][0] == 0
    assert P3[0][0] > 0
    assert P4[0][0] > 0
    # In fact every entry of P^3 and P^4 is strictly positive!
    assert all(P3[i][j] > 0 for i in range(m) for j in range(m))
    assert all(P4[i][j] > 0 for i in range(m) for j in range(m))
    print(
        f"[PASS] Irreducibility & aperiodicity verified: "
        f"p_00^(1)={P[0][0]}, p_00^(2)={P2[0][0]}, p_00^(3)={P3[0][0]}, p_00^(4)={P4[0][0]} "
        f"(gcd(3,4)=1; all entries of P^3 and P^4 are strictly positive)."
    )

    # 3. Verify uniform distribution pi is the unique invariant distribution
    pi = [Fraction(1, m)] * m
    piP = [sum(pi[i] * P[i][j] for i in range(m)) for j in range(m)]
    assert piP == pi
    PT_minus_I = [
        [P[j][i] - (1 if i == j else 0) for j in range(m)] for i in range(m)
    ]
    rk = matrix_rank_exact(PT_minus_I)
    assert rk == m - 1
    print(
        f"[PASS] pi = (1/13,...,1/13) satisfies pi*P = pi and rank(P^T - I) = {rk} = 12 (unique)."
    )

    # 4. Brute-force check for n = 1..5 against exact matrix powers
    for n in range(1, 6):
        Pn = mat_pow(P, n)
        cnt = sum(
            1 for rolls in itertools.product(range(1, 7), repeat=n) if sum(rolls) % 13 == 0
        )
        bf_prob = Fraction(cnt, 6 ** n)
        assert Pn[0][0] == bf_prob, f"Mismatch at n={n}: {Pn[0][0]} != {bf_prob}"
    print("[PASS] Brute-force die enumeration matches p_00^(n) exactly for n = 1..5.")

    # 5. Track exact convergence of p_00^(n) to 1/13
    target = Fraction(1, 13)
    print("       Convergence of exact p_00^(n) to 1/13 = 0.076923076923...:")
    prev_err = None
    for n in [3, 4, 5, 10, 20, 30, 50, 80]:
        Pn = mat_pow(P, n)
        err = abs(float(Pn[0][0] - target))
        print(f"         n={n:2d}: p_00^(n) = {float(Pn[0][0]):.15f}, |err| = {err:.3e}")
        if n >= 10 and prev_err is not None:
            assert err < prev_err
        prev_err = err
    assert prev_err < 1e-12

    # 6. Check Fourier/eigenvalue moduli directly
    omega = cmath.exp(2j * math.pi / 13)
    eigs = [sum(omega ** (k * r) for r in range(1, 7)) / 6.0 for k in range(13)]
    max_non_trivial = max(abs(eigs[k]) for k in range(1, 13))
    assert max_non_trivial < 1.0
    # Cross-check DFT formula at n=20 against exact P^20
    dft_20 = sum(lam ** 20 for lam in eigs).real / 13.0
    exact_20 = float(mat_pow(P, 20)[0][0])
    assert abs(dft_20 - exact_20) < 1e-14
    print(
        f"[PASS] Fourier cross-check: max_{{k!=0}} |lambda_k| = {max_non_trivial:.6f} < 1; "
        f"DFT matches exact P^20 to {abs(dft_20 - exact_20):.2e}."
    )
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
