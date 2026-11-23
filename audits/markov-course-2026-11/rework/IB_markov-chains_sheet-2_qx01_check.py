#!/usr/bin/env python3
"""Exact rational verification for IB/markov-chains/sheet-2/qx01 (2012-13 Sheet 2, Q6)."""

from fractions import Fraction
from itertools import product


def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def mat_pow(A, k):
    n = len(A)
    res = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    base = [row[:] for row in A]
    while k > 0:
        if k & 1:
            res = mat_mul(res, base)
        base = mat_mul(base, base)
        k >>= 1
    return res


def solve_linear(A, b):
    """Exact Gaussian elimination over Q for square non-singular A x = b."""
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for col in range(n):
        pivot = None
        for r in range(col, n):
            if M[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            raise ValueError("Singular matrix")
        M[col], M[pivot] = M[pivot], M[col]
        pv = M[col][col]
        M[col] = [val / pv for val in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                factor = M[r][col]
                M[r] = [M[r][c] - factor * M[col][c] for c in range(n + 1)]
    return [M[i][n] for i in range(n)]


def matrix_rank(A):
    """Exact rank over Q."""
    m, n = len(A), len(A[0])
    M = [row[:] for row in A]
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
        M[rank] = [val / pv for val in M[rank]]
        for r in range(m):
            if r != rank and M[r][col] != 0:
                factor = M[r][col]
                M[r] = [M[r][c] - factor * M[rank][c] for c in range(n)]
        rank += 1
    return rank


def hamming(u, v):
    return sum(a != b for a, b in zip(u, v))


def main():
    vertices = list(product([0, 1], repeat=3))
    n = len(vertices)
    assert n == 8
    idx_i = vertices.index((0, 0, 0))
    idx_o = vertices.index((1, 1, 1))

    # Build 8x8 transition matrix P
    P = [[Fraction(1, 3) if hamming(u, v) == 1 else Fraction(0) for v in vertices] for u in vertices]

    # 1. Verify stochasticity, symmetry, irreducibility, periodicity
    for row in P:
        assert all(x >= 0 for x in row)
        assert sum(row) == 1
    for a in range(n):
        for b in range(n):
            assert P[a][b] == P[b][a]

    P2 = mat_pow(P, 2)
    P3 = mat_pow(P, 3)
    # Check every pair communicates within 3 steps
    for a in range(n):
        for b in range(n):
            assert P[a][b] + P2[a][b] + P3[a][b] > 0

    # Check bipartite / period 2
    for step in range(1, 10):
        Pk = mat_pow(P, step)
        if step % 2 == 1:
            assert Pk[idx_i][idx_i] == 0
        else:
            assert Pk[idx_i][idx_i] > 0

    # 2. Verify unique invariant distribution pi = (1/8, ..., 1/8)
    pi = [Fraction(1, 8)] * n
    pi_P = [sum(pi[a] * P[a][b] for a in range(n)) for b in range(n)]
    assert pi_P == pi
    I_minus_P = [[Fraction(1 if a == b else 0) - P[a][b] for b in range(n)] for a in range(n)]
    assert matrix_rank(I_minus_P) == n - 1  # unique invariant distribution

    # 3. Part (i): Expected return time to i from full 8-vertex first-step analysis
    other_from_i = [v for v in range(n) if v != idx_i]
    A_i = [
        [Fraction(1 if u == w else 0) - P[u][w] for w in other_from_i]
        for u in other_from_i
    ]
    b_i = [Fraction(1) for _ in other_from_i]
    k_vec = solve_linear(A_i, b_i)
    k_map = {idx_i: Fraction(0)}
    for u, val in zip(other_from_i, k_vec):
        k_map[u] = val
    m_i = Fraction(1) + sum(P[idx_i][w] * k_map[w] for w in range(n))
    assert m_i == Fraction(8)
    assert m_i == Fraction(1) / pi[idx_i]

    # 4. Part (ii): Expected visits to o before first return to i (full 8-vertex solve)
    # Also check gamma_v^i = pi_v / pi_i = 1 for every v != i
    gamma_from_i = {idx_i: Fraction(1)}
    for target in other_from_i:
        b_v = [Fraction(1 if u == target else 0) for u in other_from_i]
        v_sol = solve_linear(A_i, b_v)
        v_map = {u: val for u, val in zip(other_from_i, v_sol)}
        gamma_target = sum(P[idx_i][w] * v_map[w] for w in other_from_i)
        gamma_from_i[target] = gamma_target
        assert gamma_target == pi[target] / pi[idx_i] == Fraction(1)
    assert gamma_from_i[idx_o] == Fraction(1)

    # 5. Part (iii): Expected hitting time to o from full 8-vertex first-step analysis
    other_from_o = [v for v in range(n) if v != idx_o]
    A_o = [
        [Fraction(1 if u == w else 0) - P[u][w] for w in other_from_o]
        for u in other_from_o
    ]
    b_o = [Fraction(1) for _ in other_from_o]
    h_vec = solve_linear(A_o, b_o)
    h_map = {idx_o: Fraction(0)}
    for u, val in zip(other_from_o, h_vec):
        h_map[u] = val
    for u_idx, u in enumerate(vertices):
        d_uo = hamming(u, (1, 1, 1))
        expected_by_dist = {0: Fraction(0), 1: Fraction(7), 2: Fraction(9), 3: Fraction(10)}[d_uo]
        assert h_map[u_idx] == expected_by_dist
    assert h_map[idx_i] == Fraction(10)

    # 6. Verify lumpability to distance chain Y_n = d(X_n, i) on {0, 1, 2, 3}
    P_tilde = [
        [Fraction(0), Fraction(1), Fraction(0), Fraction(0)],
        [Fraction(1, 3), Fraction(0), Fraction(2, 3), Fraction(0)],
        [Fraction(0), Fraction(2, 3), Fraction(0), Fraction(1, 3)],
        [Fraction(0), Fraction(0), Fraction(1), Fraction(0)],
    ]
    for u_idx, u in enumerate(vertices):
        du = hamming(u, (0, 0, 0))
        for dv in range(4):
            prob = sum(P[u_idx][v_idx] for v_idx, v in enumerate(vertices) if hamming(v, (0, 0, 0)) == dv)
            assert prob == P_tilde[du][dv]

    # Verify hitting probability P_1(H_3 < H_0) = 2/5, P_2(H_3 < H_0) = 3/5
    # h1 = (2/3) h2, h2 = (2/3) h1 + 1/3
    h_prob = solve_linear(
        [[Fraction(1), Fraction(-2, 3)], [Fraction(-2, 3), Fraction(1)]],
        [Fraction(0), Fraction(1, 3)],
    )
    assert h_prob == [Fraction(2, 5), Fraction(3, 5)]
    # Expected visits via geometric distribution: h1 * 1/(1 - h2) = (2/5) * (5/2) = 1
    assert h_prob[0] / (1 - h_prob[1]) == Fraction(1)

    # 7. Independent time-evolution / exact partial-sum + exact Markov tail check up to N = 200 steps
    N_steps = 200
    # Track distribution absorbed at 0 (for T_i and visits to o before T_i)
    # dist[s] = P_i(Y_step = s, T_i > step) for s in {1, 2, 3}
    dist = [Fraction(0), Fraction(1), Fraction(0), Fraction(0)]  # at step = 1
    partial_E_Ti = Fraction(0)
    partial_gamma_o = Fraction(0)
    for step in range(1, N_steps):
        partial_gamma_o += dist[3]
        next_dist = [Fraction(0)] * 4
        for s in (1, 2, 3):
            for t in range(4):
                next_dist[t] += dist[s] * P_tilde[s][t]
        prob_return_at_next = next_dist[0]
        partial_E_Ti += (step + 1) * prob_return_at_next
        next_dist[0] = Fraction(0)
        dist = next_dist
    # Now dist is at step = N_steps; verify exact Markov tail identities!
    k_by_dist = [Fraction(0), Fraction(7), Fraction(9), Fraction(10)]
    v_by_dist = [Fraction(0), Fraction(1), Fraction(3, 2), Fraction(5, 2)]
    exact_E_Ti_with_tail = partial_E_Ti + sum(dist[s] * (N_steps + k_by_dist[s]) for s in (1, 2, 3))
    exact_gamma_o_with_tail = partial_gamma_o + sum(dist[s] * v_by_dist[s] for s in (1, 2, 3))
    assert exact_E_Ti_with_tail == Fraction(8)
    assert exact_gamma_o_with_tail == Fraction(1)

    # Track distribution absorbed at 3 (for H_o starting from 0)
    dist_o = [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]  # at step = 0
    partial_E_Ho = Fraction(0)
    for step in range(0, N_steps):
        next_dist = [Fraction(0)] * 4
        for s in (0, 1, 2):
            for t in range(4):
                next_dist[t] += dist_o[s] * P_tilde[s][t]
        partial_E_Ho += (step + 1) * next_dist[3]
        next_dist[3] = Fraction(0)
        dist_o = next_dist
    # Now dist_o is at step = N_steps; distance to o is 3 - s
    exact_E_Ho_with_tail = partial_E_Ho + sum(dist_o[s] * (N_steps + k_by_dist[3 - s]) for s in (0, 1, 2))
    assert exact_E_Ho_with_tail == Fraction(10)

    assert abs(float(partial_E_Ti) - 8.0) < 1e-8
    assert abs(float(partial_gamma_o) - 1.0) < 1e-9
    assert abs(float(partial_E_Ho) - 10.0) < 1e-8

    print("IB/markov-chains/sheet-2/qx01 verification PASSED")
    print(f"Unique invariant distribution: pi_v = {pi[0]} for all 8 vertices (rank(I-P) = {matrix_rank(I_minus_P)})")
    print(f"(i)   E_i[T_i] (exact 8-vertex solve) = {m_i} (1/pi_i = {Fraction(1)/pi[idx_i]})")
    print(f"(ii)  gamma_o^i (exact 8-vertex solve) = {gamma_from_i[idx_o]} (pi_o/pi_i = {pi[idx_o]/pi[idx_i]})")
    print(f"(iii) E_i[H_o] (exact 8-vertex solve) = {h_map[idx_i]}")
    print(f"      E_v[H_o] by distance to o: d=0 -> {h_map[idx_o]}, d=1 -> 7, d=2 -> 9, d=3 -> 10")
    print(f"Hitting probabilities on distance chain: P_1(H_3 < H_0) = {h_prob[0]}, P_2(H_3 < H_0) = {h_prob[1]}")
    print(f"200-step exact partial sums: E_i[T_i] ~ {float(partial_E_Ti):.10f}, gamma_o^i ~ {float(partial_gamma_o):.10f}, E_i[H_o] ~ {float(partial_E_Ho):.10f}")
    print(f"Exact partial + Markov tail: E_i[T_i] = {exact_E_Ti_with_tail}, gamma_o^i = {exact_gamma_o_with_tail}, E_i[H_o] = {exact_E_Ho_with_tail}")


if __name__ == "__main__":
    main()
