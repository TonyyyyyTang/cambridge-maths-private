#!/usr/bin/env python3
"""Exact rational verification for IB/2021/P2/q18H."""

from fractions import Fraction


def solve_linear_system(A, b):
    """Exact Gaussian elimination over Fraction for A x = b."""
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        pivot = None
        for r in range(col, n):
            if M[r][col] != 0:
                pivot = r
                break
        assert pivot is not None, "Singular matrix"
        M[col], M[pivot] = M[pivot], M[col]
        pv = M[col][col]
        for c in range(col, n + 1):
            M[col][c] /= pv
        for r in range(n):
            if r != col and M[r][col] != 0:
                factor = M[r][col]
                for c in range(col, n + 1):
                    M[r][c] -= factor * M[col][c]
    return [M[i][n] for i in range(n)]


def p_trans(i, j, q, p):
    if i == 0:
        if j == 1:
            return p
        if j == 0:
            return 1 - p
        return Fraction(0)
    if j == i + 1:
        return q ** (-(i + 2))
    if j == i - 1:
        return q ** (-i)
    if j == i:
        return 1 - q ** (-(i + 2)) - q ** (-i)
    return Fraction(0)


def check_infinite_identities(q, p, max_i=30):
    # Formula for pi
    denom = p * q * q + q - 1
    pi_0 = (q - 1) / denom
    pi = [pi_0] + [p * q * (q - 1) / denom * (q ** (-(i - 1))) for i in range(1, max_i + 2)]

    # 1. Check valid transition probabilities and row sums for i = 0..max_i
    for i in range(max_i + 1):
        nbrs = [j for j in (i - 1, i, i + 1) if j >= 0]
        probs = [p_trans(i, j, q, p) for j in nbrs]
        assert all(val >= 0 for val in probs)
        assert sum(probs) == 1

    # 2. Check detailed balance pi_i p_{i,j} = pi_j p_{j,i}
    for i in range(max_i + 1):
        for j in (i - 1, i, i + 1):
            if 0 <= j <= max_i:
                assert pi[i] * p_trans(i, j, q, p) == pi[j] * p_trans(j, i, q, p)

    # 3. Check invariance (pi P)_j = pi_j for j = 0..max_i
    for j in range(max_i + 1):
        lhs = sum(pi[i] * p_trans(i, j, q, p) for i in (j - 1, j, j + 1) if i >= 0)
        assert lhs == pi[j]

    # 4. Check exact hitting times k_i = E_i[H_0] = q^2 (q^i - 1) / (q - 1)^2
    k = [q * q * (q ** i - 1) / ((q - 1) ** 2) for i in range(max_i + 2)]
    assert k[0] == 0
    assert k[1] == q * q / (q - 1)
    for i in range(1, max_i + 1):
        rhs = (
            1
            + p_trans(i, i - 1, q, p) * k[i - 1]
            + p_trans(i, i, q, p) * k[i]
            + p_trans(i, i + 1, q, p) * k[i + 1]
        )
        assert k[i] == rhs

    # 5. Check m_0 = 1/pi_0 = 1 + p * k_1
    m0 = Fraction(1) / pi[0]
    assert m0 == 1 + p * k[1]

    # 6. Check closed-form answers for (b)(i), (b)(ii), (b)(iii)
    ans_bi = q / (q * q - 1)
    ans_bii = (p * q * q + q - 1) / (p * q * (q - 1))
    ans_biii = q * q / (q - 1)
    assert ans_bii == Fraction(1) / pi[1]
    assert ans_biii == k[1]

    # Check partial sums of gamma_{2m}^1 = pi_{2m}/pi_1
    partial_even = sum(pi[2 * m] / pi[1] for m in range(1, max_i // 2 + 1))
    exact_partial_even = q * (1 - q ** (-2 * (max_i // 2))) / (q * q - 1)
    assert partial_even == exact_partial_even
    assert 0 < ans_bi - partial_even < Fraction(1, 10**4)


def check_truncated_chain(q, p, N):
    """Reflect at state N (set p_{N,N} = 1 - q^{-N}) and solve linear systems directly."""
    P = [[Fraction(0) for _ in range(N + 1)] for _ in range(N + 1)]
    for i in range(N + 1):
        for j in (i - 1, i, i + 1):
            if 0 <= j <= N:
                P[i][j] = p_trans(i, j, q, p)
        if i == N:
            # Absorb the missing upward jump into self-loop so detailed balance is preserved
            P[N][N] = 1 - P[N][N - 1]
        assert sum(P[i]) == 1

    # Direct calculation of E_i[T_1] via first-step equations for hitting 1 from i != 1:
    # h_i = 1 + sum_{j != 1} P[i][j] h_j for i != 1; then m_1 = 1 + sum_{j != 1} P[1][j] h_j.
    others = [i for i in range(N + 1) if i != 1]
    idx = {s: r for r, s in enumerate(others)}
    A = [[Fraction(0) for _ in others] for _ in others]
    b = [Fraction(1) for _ in others]
    for r, s in enumerate(others):
        A[r][r] = Fraction(1)
        for j in others:
            A[r][idx[j]] -= P[s][j]
    h_vals = solve_linear_system(A, b)
    h_to_1 = {s: h_vals[idx[s]] for s in others}
    m1_N = 1 + sum(P[1][j] * h_to_1[j] for j in others)

    # Direct calculation of expected time in E = {2, 4, ..., <= N} before returning to 1:
    # u_i = E_i[time in E before hitting 1] for i != 1:
    # u_i = 1_{i in E} + sum_{j != 1} P[i][j] u_j; then answer = sum_{j != 1} P[1][j] u_j
    b_even = [Fraction(1 if (s >= 2 and s % 2 == 0) else 0) for s in others]
    u_vals = solve_linear_system(A, b_even)
    u_even = {s: u_vals[idx[s]] for s in others}
    even_N = sum(P[1][j] * u_even[j] for j in others)

    # Direct calculation of k_i^{(N)} = E_i[H_0] for i = 1..N:
    A0 = [[Fraction(0) for _ in range(N)] for _ in range(N)]
    b0 = [Fraction(1) for _ in range(N)]
    for i in range(1, N + 1):
        A0[i - 1][i - 1] = Fraction(1)
        for j in range(1, N + 1):
            A0[i - 1][j - 1] -= P[i][j]
    k_N = solve_linear_system(A0, b0)
    k1_N = k_N[0]

    # Exact closed-form formulas for the reflected chain of size N:
    # pi_i^{(N)} proportional to pi_i for 0 <= i <= N, so:
    # even_N = sum_{2<=2m<=N} q^{-(2m-1)} = q (1 - q^{-2*(N//2)}) / (q^2 - 1)
    # m1_N   = 1/(pq) + sum_{i=1}^N q^{-(i-1)} = 1/(pq) + q (1 - q^{-N}) / (q - 1)
    # k1_N   = sum_{j=1}^N q^{2-j} = q^2 (1 - q^{-N}) / (q - 1)
    assert even_N == q * (1 - q ** (-2 * (N // 2))) / (q * q - 1)
    assert m1_N == Fraction(1) / (p * q) + q * (1 - q ** (-N)) / (q - 1)
    assert k1_N == q * q * (1 - q ** (-N)) / (q - 1)
    return even_N, m1_N, k1_N


def main():
    test_cases = [
        (Fraction(2), Fraction(1)),
        (Fraction(2), Fraction(1, 2)),
        (Fraction(2), Fraction(1, 3)),
        (Fraction(3), Fraction(2, 5)),
        (Fraction(5), Fraction(3, 7)),
        (Fraction(5, 2), Fraction(1, 4)),
    ]
    for q, p in test_cases:
        check_infinite_identities(q, p, max_i=30)
        e20, m20, k20 = check_truncated_chain(q, p, N=20)
        ans_bi = q / (q * q - 1)
        ans_bii = (p * q * q + q - 1) / (p * q * (q - 1))
        ans_biii = q * q / (q - 1)
        print(
            f"(q={q}, p={p}): "
            f"(b)(i)={ans_bi} [N=20 err={float(ans_bi - e20):.2e}], "
            f"(b)(ii)={ans_bii} [N=20 err={float(ans_bii - m20):.2e}], "
            f"(b)(iii)={ans_biii} [N=20 err={float(ans_biii - k20):.2e}]"
        )
    print("ALL CHECKS PASSED for IB/2021/P2/q18H.")


if __name__ == "__main__":
    main()
