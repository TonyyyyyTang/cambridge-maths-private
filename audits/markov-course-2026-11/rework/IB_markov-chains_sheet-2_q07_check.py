#!/usr/bin/env python3
"""Exact rational verification for IB/markov-chains/sheet-2/q07 (Random walk on rooted binary tree)."""

from fractions import Fraction
import math


def solve_linear_exact(A, b):
    """Exact rational Gaussian elimination solving Ax = b."""
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            raise ValueError("Singular matrix")
        M[col], M[pivot] = M[pivot], M[col]
        pv = M[col][col]
        for j in range(col, n + 1):
            M[col][j] /= pv
        for row in range(n):
            if row != col and M[row][col] != 0:
                factor = M[row][col]
                for j in range(col, n + 1):
                    M[row][j] -= factor * M[col][j]
    return [M[i][n] for i in range(n)]


def build_truncated_tree(max_depth):
    """Build rooted binary tree up to max_depth.
    Vertex 0 is R (depth 0, degree 1); vertex 1 has depth 1;
    every vertex v >= 1 with depth < max_depth has 1 parent and 2 children.
    """
    depths = {0: 0}
    adj = {0: []}
    next_id = 1
    # R has single child 1
    if max_depth >= 1:
        adj[0].append(1)
        adj[1] = [0]
        depths[1] = 1
        next_id = 2
        frontier = [1]
        for d in range(1, max_depth):
            new_frontier = []
            for v in frontier:
                for _ in range(2):
                    c = next_id
                    next_id += 1
                    depths[c] = d + 1
                    adj[v].append(c)
                    adj[c] = [v]
                    new_frontier.append(c)
            frontier = new_frontier
    return depths, adj


def check_tree_lumpability_and_hitting(max_depth=6):
    depths, adj = build_truncated_tree(max_depth)
    n = len(depths)
    # Verify vertex count: 1 + (1 + 2 + ... + 2^(max_depth-1)) = 2^max_depth
    assert n == 2 ** max_depth

    # Verify degrees and lumpability for all interior vertices (depth < max_depth)
    for v in range(n):
        d = depths[v]
        if d == 0:
            assert len(adj[v]) == 1
            assert depths[adj[v][0]] == 1
        elif d < max_depth:
            assert len(adj[v]) == 3
            p_down = sum(Fraction(1, 3) for u in adj[v] if depths[u] == d - 1)
            p_up = sum(Fraction(1, 3) for u in adj[v] if depths[u] == d + 1)
            p_same = sum(Fraction(1, 3) for u in adj[v] if depths[u] == d)
            assert p_down == Fraction(1, 3)
            assert p_up == Fraction(2, 3)
            assert p_same == 0

    # Solve exact hitting probability h_v = P_v(hit R before hitting depth max_depth)
    A = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    b = [Fraction(0) for _ in range(n)]
    for v in range(n):
        d = depths[v]
        A[v][v] = Fraction(1)
        if d == 0:
            b[v] = Fraction(1)
        elif d == max_depth:
            b[v] = Fraction(0)
        else:
            deg = len(adj[v])
            for u in adj[v]:
                A[v][u] -= Fraction(1, deg)
    h = solve_linear_exact(A, b)

    # Verify h_v depends only on depth d(v) and equals (2^-d - 2^-L)/(1 - 2^-L)
    for v in range(n):
        d = depths[v]
        expected = (Fraction(1, 2 ** d) - Fraction(1, 2 ** max_depth)) / (
            1 - Fraction(1, 2 ** max_depth)
        )
        assert h[v] == expected, f"Mismatch at v={v}, d={d}: {h[v]} != {expected}"
    print(
        f"[PASS] Truncated tree T_{max_depth} ({n} vertices): lumpability verified and "
        f"exact hitting prob h_v = (2^-d - 2^-{max_depth})/(1 - 2^-{max_depth}) for all vertices."
    )
    print(f"       h_1^(L={max_depth}) = {h[1]} = {float(h[1]):.8f} (limit as L->inf: 1/2)")


def check_distance_chain_returns(max_steps=30):
    # Compute exact n-step distribution of Y_n starting from Y_0 = 0 (horizon max_steps)
    M = max_steps + 2
    dist = [Fraction(0)] * M
    dist[0] = Fraction(1)
    q_n_00 = [dist[0]]
    for _ in range(max_steps):
        ndist = [Fraction(0)] * M
        ndist[1] += dist[0]
        for i in range(1, M - 1):
            if dist[i] != 0:
                ndist[i - 1] += dist[i] * Fraction(1, 3)
                ndist[i + 1] += dist[i] * Fraction(2, 3)
        dist = ndist
        q_n_00.append(dist[0])

    # Compute first-return probabilities f_00^(n) via renewal equation:
    # q^(n)(0,0) = sum_{k=1}^n f_00^(k) q^(n-k)(0,0) for n >= 1
    f_00 = [Fraction(0)] * (max_steps + 1)
    for n in range(1, max_steps + 1):
        s = q_n_00[n]
        for k in range(1, n):
            s -= f_00[k] * q_n_00[n - k]
        f_00[n] = s

    # Verify odd steps are 0 and even steps match Catalan formula:
    # f_00^(2m) = (1/3) * C_{m-1} * (2/9)^(m-1)
    for n in range(1, max_steps + 1):
        if n % 2 == 1:
            assert f_00[n] == 0
        else:
            m = n // 2
            cat = math.comb(2 * (m - 1), m - 1) // m
            expected = Fraction(1, 3) * cat * Fraction(2, 9) ** (m - 1)
            assert f_00[n] == expected, f"f_00^({n}) mismatch: {f_00[n]} != {expected}"

    cum_return = sum(f_00)
    cum_visits = sum(q_n_00)
    assert cum_return < Fraction(1, 2)
    assert cum_visits < Fraction(2)
    print(
        f"[PASS] Exact Catalan formula f_00^(2m) = (1/3)*C_(m-1)*(2/9)^(m-1) verified for m=1..{max_steps//2}."
    )
    print(
        f"       Partial sum P_0(T_0 <= {max_steps}) = {float(cum_return):.10f} (< 0.5, exact limit 1/2)"
    )
    print(
        f"       Partial sum sum_{{n=0}}^{max_steps} q^(n)(0,0) = {float(cum_visits):.10f} (< 2.0, exact limit 2)"
    )


def check_monotone_bound(max_n=30, max_k=30):
    # Verify h_k^(n) = P_k(H_0 <= n) <= 2^-k by exact DP recursion
    h = [Fraction(1) if k == 0 else Fraction(0) for k in range(max_k + max_n + 2)]
    for n in range(max_n):
        for k in range(max_k + 1):
            assert h[k] <= Fraction(1, 2 ** k)
        nh = [Fraction(0)] * len(h)
        nh[0] = Fraction(1)
        for k in range(1, len(h) - 1):
            nh[k] = Fraction(2, 3) * h[k + 1] + Fraction(1, 3) * h[k - 1]
        h = nh
    for k in range(max_k + 1):
        assert h[k] <= Fraction(1, 2 ** k)
    print(
        f"[PASS] Direct induction bound h_k^(n) <= 2^(-k) verified for 0 <= n,k <= {max_n}; "
        f"h_1^({max_n}) = {float(h[1]):.10f}."
    )


if __name__ == "__main__":
    check_tree_lumpability_and_hitting(max_depth=6)
    check_distance_chain_returns(max_steps=30)
    check_monotone_bound(max_n=30, max_k=30)
    print("ALL CHECKS PASSED.")
