"""Exact rational verification for IB/2016/P4/q9H (Ehrenfest urn model).

Verifies for k = 1, 2, ..., 12:
  1. Transition matrix P is stochastic and irreducible.
  2. pi_i = C(k, i) * 2^{-k} satisfies detailed balance pi_i p_{ij} = pi_j p_{ji}
     for all i, j in {0, ..., k}.
  3. pi is the unique invariant distribution (pi P = pi, sum(pi) = 1, rank(P^T - I) = k).
  4. Direct first-step linear solve for h_i = E_i[inf{n >= 1 : X_n = 0}] gives
     E_0[T] = 1 + h_1 = 2^k = 1 / pi_0.
"""

from fractions import Fraction
import math


def solve_linear(A, b):
    n = len(b)
    M = [list(map(Fraction, row)) + [Fraction(val)] for row, val in zip(A, b)]
    for col in range(n):
        pivot = None
        for r in range(col, n):
            if M[r][col] != 0:
                pivot = r
                break
        assert pivot is not None, "Singular matrix"
        M[col], M[pivot] = M[pivot], M[col]
        pv = M[col][col]
        for j in range(col, n + 1):
            M[col][j] /= pv
        for r in range(n):
            if r != col and M[r][col] != 0:
                factor = M[r][col]
                for j in range(col, n + 1):
                    M[r][j] -= factor * M[col][j]
    return [M[r][n] for r in range(n)]


def matrix_rank(A):
    m = len(A)
    n = len(A[0])
    M = [list(map(Fraction, row)) for row in A]
    rank = 0
    for col in range(n):
        pivot = None
        for r in range(rank, m):
            if M[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        M[rank], M[pivot] = M[pivot], M[rank]
        pv = M[rank][col]
        for j in range(col, n):
            M[rank][j] /= pv
        for r in range(m):
            if r != rank and M[r][col] != 0:
                factor = M[r][col]
                for j in range(col, n):
                    M[r][j] -= factor * M[rank][j]
        rank += 1
    return rank


def ehrenfest_matrix(k):
    n = k + 1
    P = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        if i + 1 <= k:
            P[i][i + 1] = Fraction(k - i, k)
        if i - 1 >= 0:
            P[i][i - 1] = Fraction(i, k)
    return P


def main():
    print("=== IB/2016/P4/q9H Exact Verification ===")
    for k in range(1, 13):
        P = ehrenfest_matrix(k)
        n = k + 1

        # 1. Stochasticity & irreducibility (adjacent transitions strictly positive)
        for i in range(n):
            assert all(x >= 0 for x in P[i])
            assert sum(P[i]) == 1
            if i < k:
                assert P[i][i + 1] > 0
            if i > 0:
                assert P[i][i - 1] > 0

        # 2. Candidate equilibrium distribution pi_i = C(k, i) / 2^k
        pi = [Fraction(math.comb(k, i), 2**k) for i in range(n)]
        assert sum(pi) == 1

        # Detailed balance: pi_i P_{ij} = pi_j P_{ji} for all i, j
        for i in range(n):
            for j in range(n):
                assert pi[i] * P[i][j] == pi[j] * P[j][i]

        # 3. Invariance pi P = pi and uniqueness via rank(P^T - I) == k
        piP = [sum(pi[i] * P[i][j] for i in range(n)) for j in range(n)]
        assert piP == pi
        PT_minus_I = [
            [P[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)]
            for i in range(n)
        ]
        assert matrix_rank(PT_minus_I) == k

        # 4. Direct first-step analysis for hitting time of 0 from states 1..k:
        #    h_i = 1 + sum_{j=1}^k P[i][j] h_j  (since h_0 = 0 for hitting 0)
        #    => (I - Q) h = 1, then E_0[T] = 1 + P[0][1] h_1 = 1 + h_1.
        A = [
            [
                (Fraction(1) if a == b else Fraction(0)) - P[a][b]
                for b in range(1, n)
            ]
            for a in range(1, n)
        ]
        b_vec = [Fraction(1)] * k
        h = solve_linear(A, b_vec)  # h[0] is h_1, ..., h[k-1] is h_k
        E_T = Fraction(1) + h[0]
        assert E_T == Fraction(2**k)
        assert E_T == Fraction(1) / pi[0]
        print(
            f"  k = {k:2d}: detailed balance OK, unique pi OK, "
            f"h_1 = {h[0]}, E_0[T] = {E_T} == 2^{k} == 1/pi_0"
        )

    print("\nALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
