"""Exact rational verification for IB/2016/P3/q9H.

Verifies the identity:
    sum_{n=0}^infty P_i(X_n = i) = 1 / P_i(X_n != i for all n >= 1)
and the intermediate identity:
    P_i(T_i^(r) < infty) = f_i^r
on concrete finite Markov chains covering all three regimes:
    1. Transient state with 0 < f_i < 1
    2. Degenerate transient state with f_i = 0
    3. Recurrent state with f_i = 1
"""

from fractions import Fraction


def solve_linear(A, b):
    """Exact Gaussian elimination over Q for A x = b."""
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


def mat_mul(A, B):
    n = len(A)
    return [
        [sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def check_stochastic(P):
    n = len(P)
    for row in P:
        assert len(row) == n
        assert all(x >= 0 for x in row)
        assert sum(row) == 1


def exact_hitting_probs(P, target, absorbing_set):
    """Solve h_x = P_x(hit target before never hitting) for chain with absorbing_set.

    Here target is treated as absorbing (h_target = 1), and states in absorbing_set
    (disjoint from {target}) have h_s = 0. Assumes all remaining states can reach
    target or absorbing_set so (I - Q) is invertible.
    """
    n = len(P)
    unknowns = [s for s in range(n) if s != target and s not in absorbing_set]
    idx = {s: k for k, s in enumerate(unknowns)}
    m = len(unknowns)
    A = [[Fraction(0)] * m for _ in range(m)]
    b = [Fraction(0)] * m
    for s in unknowns:
        r = idx[s]
        A[r][r] = Fraction(1)
        for t in range(n):
            if t == target:
                b[r] += P[s][t]
            elif t in idx:
                A[r][idx[t]] -= P[s][t]
    sol = solve_linear(A, b)
    h = {}
    for s in range(n):
        if s == target:
            h[s] = Fraction(1)
        elif s in absorbing_set:
            h[s] = Fraction(0)
        else:
            h[s] = sol[idx[s]]
    return h


def exact_green_function(P, transient_states):
    """Compute G = (I - Q)^{-1} on transient_states, so G[i][j] = sum_{n>=0} P^n[i][j]."""
    m = len(transient_states)
    A = [
        [
            (Fraction(1) if a == b else Fraction(0)) - P[transient_states[a]][transient_states[b]]
            for b in range(m)
        ]
        for a in range(m)
    ]
    G = []
    for col in range(m):
        e = [Fraction(1) if r == col else Fraction(0) for r in range(m)]
        # Solve (I - Q) x = e -> col of (I - Q)^{-1}
        x = solve_linear(A, e)
        G.append(x)
    # Transpose so G_mat[a][b] = (I - Q)^{-1}_{a,b}
    return [[G[b][a] for b in range(m)] for a in range(m)]


def exact_rth_return_prob(P, state_i, absorbing_set, max_r):
    """Compute P_i(T_i^(r) < infty) for r = 0..max_r on lifted state space (state, visits)."""
    # Lifted chain tracks (s, k) where k in {0, 1, ..., max_r} is number of returns to state_i
    # made at times >= 1. Entering (state_i, r) is absorbed at (state_i, r) when computing
    # P(T_i^(r) < infty).
    n = len(P)
    results = [Fraction(1)]  # r = 0: T_i^(0) = 0 < infty always
    for r in range(1, max_r + 1):
        # States: (s, k) with s in 0..n-1, k in 0..r-1; goal is hitting (state_i, r)
        # Absorbing failure states: (s, k) with s in absorbing_set
        non_abs = [s for s in range(n) if s not in absorbing_set]
        nodes = [(s, k) for k in range(r) for s in non_abs]
        idx = {node: pos for pos, node in enumerate(nodes)}
        m = len(nodes)
        A = [[Fraction(0)] * m for _ in range(m)]
        b = [Fraction(0)] * m
        for s, k in nodes:
            row = idx[(s, k)]
            A[row][row] = Fraction(1)
            for t in range(n):
                p = P[s][t]
                if p == 0 or t in absorbing_set:
                    continue
                nk = k + 1 if t == state_i else k
                if nk == r:
                    b[row] += p
                else:
                    A[row][idx[(t, nk)]] -= p
        sol = solve_linear(A, b)
        results.append(sol[idx[(state_i, 0)]])
    return results


def main():
    print("=== IB/2016/P3/q9H Exact Verification ===")

    # Case 1: Transient states {0, 1, 2}, absorbing state {3}
    P1 = [
        [Fraction(1, 5), Fraction(2, 5), Fraction(1, 5), Fraction(1, 5)],
        [Fraction(1, 3), Fraction(1, 6), Fraction(1, 3), Fraction(1, 6)],
        [Fraction(1, 4), Fraction(1, 4), Fraction(1, 4), Fraction(1, 4)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(1)],
    ]
    check_stochastic(P1)
    G1 = exact_green_function(P1, [0, 1, 2])
    print("\nCase 1: 4-state chain with transient states {0,1,2} and absorbing state {3}")
    for i in [0, 1, 2]:
        # First-step return prob f_i = P(i,i) + sum_{j != i} P(i,j) h_{j -> i}
        h = exact_hitting_probs(P1, target=i, absorbing_set={3})
        f_i = sum(P1[i][j] * h[j] for j in range(4))
        never_return = 1 - f_i
        lhs = G1[i][i]
        rhs = Fraction(1) / never_return
        assert lhs == rhs, f"Mismatch at state {i}: {lhs} != {rhs}"

        # Also check P_i(T_i^(r) < infty) == f_i^r for r = 0..5 via lifted chain
        ret_probs = exact_rth_return_prob(P1, i, {3}, max_r=5)
        for r, prob in enumerate(ret_probs):
            assert prob == f_i**r, f"Mismatch at r={r}: {prob} != {f_i**r}"

        print(
            f"  state {i}: f_{i} = {f_i}, P_{i}(never return) = {never_return}, "
            f"sum P^n({i},{i}) = {lhs} == 1/(1-f_{i}) = {rhs}; "
            f"P(T^(r)<inf)=f^r verified for r=0..5"
        )

    # Verify partial sums converge monotonically to G1[0][0]
    P_pow = [[Fraction(1) if a == b else Fraction(0) for b in range(4)] for a in range(4)]
    partial = Fraction(0)
    for _ in range(80):
        partial += P_pow[0][0]
        P_pow = mat_mul(P_pow, P1)
    assert partial < G1[0][0]
    assert float(G1[0][0] - partial) < 1e-6
    print(f"  80-term partial sum for state 0: {float(partial):.8f} vs exact {float(G1[0][0]):.8f}")

    # Case 2: Degenerate transient state with f_0 = 0 (strictly acyclic 0 -> 1 -> 2 -> 2)
    P2 = [
        [Fraction(0), Fraction(1, 2), Fraction(1, 2)],
        [Fraction(0), Fraction(0), Fraction(1)],
        [Fraction(0), Fraction(0), Fraction(1)],
    ]
    check_stochastic(P2)
    G2 = exact_green_function(P2, [0, 1])
    h2 = exact_hitting_probs(P2, target=0, absorbing_set={2})
    f_0 = sum(P2[0][j] * h2[j] for j in range(3))
    assert f_0 == 0
    assert G2[0][0] == 1 == Fraction(1) / (1 - f_0)
    ret_probs_0 = exact_rth_return_prob(P2, 0, {2}, max_r=4)
    assert ret_probs_0 == [Fraction(1), Fraction(0), Fraction(0), Fraction(0), Fraction(0)]
    print(
        f"\nCase 2: Degenerate f_0 = 0 chain: f_0 = {f_0}, P_0(never return) = {1 - f_0}, "
        f"sum P^n(0,0) = {G2[0][0]} == 1/1; P(T^(r)<inf) = {ret_probs_0}"
    )

    # Case 3: Recurrent chain (3-state irreducible), f_i = 1, sum diverges to +infty
    P3 = [
        [Fraction(1, 2), Fraction(1, 2), Fraction(0)],
        [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)],
        [Fraction(0), Fraction(1, 2), Fraction(1, 2)],
    ]
    check_stochastic(P3)
    h3 = exact_hitting_probs(P3, target=0, absorbing_set=set())
    f_0_rec = sum(P3[0][j] * h3[j] for j in range(3))
    assert f_0_rec == 1
    assert 1 - f_0_rec == 0
    # Compute partial sum up to N=100 to verify linear growth -> +infty
    P_pow = [[Fraction(1) if a == b else Fraction(0) for b in range(3)] for a in range(3)]
    s50 = Fraction(0)
    s100 = Fraction(0)
    for n in range(101):
        if n <= 50:
            s50 += P_pow[0][0]
        s100 += P_pow[0][0]
        P_pow = mat_mul(P_pow, P3)
    assert s100 > s50 > 15
    print(
        f"\nCase 3: Recurrent irreducible chain: f_0 = {f_0_rec}, P_0(never return) = {1 - f_0_rec}, "
        f"partial sums grow linearly: S_50 = {float(s50):.4f}, S_100 = {float(s100):.4f} -> +infty = 1/0"
    )

    print("\nALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
