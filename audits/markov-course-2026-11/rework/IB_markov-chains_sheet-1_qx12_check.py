"""Exact rational verification for IB/markov-chains/sheet-1/qx12 (referencing Example 8)."""

from fractions import Fraction


def solve_linear(A, b):
    """Exact Gaussian elimination over Q."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] != 0:
                piv = r
                break
        assert piv is not None, "Singular matrix"
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [val / pv for val in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                fac = M[r][col]
                M[r] = [M[r][c] - fac * M[col][c] for c in range(n + 1)]
    return [M[i][n] for i in range(n)]


def mat_mul(A, B):
    n = len(A)
    return [
        [sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def main():
    # Transition matrix P from Example 8 (2010-11/sheet-1.pdf, page 2), states {1, 2, 3, 4, 5}
    P = [
        [Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0), Fraction(1, 2)],
        [Fraction(0), Fraction(1, 2), Fraction(0), Fraction(1, 2), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(1), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(1, 4), Fraction(1, 4), Fraction(1, 4), Fraction(1, 4)],
        [Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0), Fraction(1, 2)],
    ]
    n = 5
    for i in range(n):
        assert all(x >= 0 for x in P[i])
        assert sum(P[i]) == 1

    # 1. Reachability matrix via Warshall transitive closure (i -> j in >= 0 steps)
    reach = [[i == j or P[i][j] > 0 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])

    # Communicating classes (1-indexed labels)
    seen = set()
    classes = []
    for i in range(n):
        if i in seen:
            continue
        cls = frozenset(j + 1 for j in range(n) if reach[i][j] and reach[j][i])
        for s in cls:
            seen.add(s - 1)
        classes.append(cls)

    closed_classes = []
    open_classes = []
    for cls in classes:
        is_closed = all(
            P[i - 1][j - 1] == 0
            for i in cls
            for j in range(1, n + 1)
            if j not in cls
        )
        if is_closed:
            closed_classes.append(cls)
        else:
            open_classes.append(cls)

    assert set(classes) == {frozenset({1, 5}), frozenset({3}), frozenset({2, 4})}
    assert set(closed_classes) == {frozenset({1, 5}), frozenset({3})}
    assert set(open_classes) == {frozenset({2, 4})}

    # 2. Fundamental matrix (I - Q)^{-1} on transient states {2, 4} (indices 1, 3)
    tr_idx = [1, 3]
    Q = [[P[tr_idx[a]][tr_idx[b]] for b in range(2)] for a in range(2)]
    I_minus_Q = [
        [Fraction(1 if a == b else 0) - Q[a][b] for b in range(2)]
        for a in range(2)
    ]
    col0 = solve_linear(I_minus_Q, [Fraction(1), Fraction(0)])
    col1 = solve_linear(I_minus_Q, [Fraction(0), Fraction(1)])
    N_mat = [[col0[0], col1[0]], [col0[1], col1[1]]]
    assert N_mat == [[Fraction(3), Fraction(2)], [Fraction(1), Fraction(2)]]

    # Exact rank-1 identity: Q^2 = (3/4) Q, so Q^m = (3/4)^{m-1} Q for m >= 1
    Q2 = mat_mul(Q, Q)
    assert Q2 == [[Fraction(3, 4) * Q[a][b] for b in range(2)] for a in range(2)]

    # Exact return probabilities from N_ii = 1 / (1 - f_i)
    f_2 = 1 - Fraction(1) / N_mat[0][0]
    f_4 = 1 - Fraction(1) / N_mat[1][1]
    assert f_2 == Fraction(2, 3)
    assert f_4 == Fraction(1, 2)

    # 3. Direct first-step analysis for h_i = P_i(H_4 < inf) and g_i = P_i(H_2 < inf)
    # h_4 = 1, h_1 = h_3 = h_5 = 0, h_2 = (1/2) h_2 + (1/2) h_4
    h2 = solve_linear([[Fraction(1, 2)]], [Fraction(1, 2)])[0]
    assert h2 == Fraction(1)
    assert P[3][1] * h2 + P[3][3] == f_4

    # g_2 = 1, g_1 = g_3 = g_5 = 0, g_4 = (1/4) g_2 + (1/4) g_4
    g4 = solve_linear([[Fraction(3, 4)]], [Fraction(1, 4)])[0]
    assert g4 == Fraction(1, 3)
    assert P[1][1] + P[1][3] * g4 == f_2

    # 4. Exact verification of P^m for m = 1..30 (exact closed forms, zero approximation)
    Pm = [row[:] for row in P]
    sum_Q_pow = [[Fraction(1 if a == b else 0) for b in range(2)] for a in range(2)]
    M = 30
    for m in range(1, M + 1):
        # Closed class {1,5} entries are identically 1/2 for all m >= 1
        assert Pm[0][0] == Fraction(1, 2) and Pm[0][4] == Fraction(1, 2)
        assert Pm[4][0] == Fraction(1, 2) and Pm[4][4] == Fraction(1, 2)
        # Absorbing state 3
        assert Pm[2][2] == Fraction(1)
        # Exact rank-1 closed form on transient block {2, 4}
        scale = Fraction(3, 4) ** (m - 1)
        for a in range(2):
            for b in range(2):
                assert Pm[tr_idx[a]][tr_idx[b]] == scale * Q[a][b]
                sum_Q_pow[a][b] += Pm[tr_idx[a]][tr_idx[b]]
        Pm = mat_mul(Pm, P)

    # Exact partial sum + geometric tail == (I - Q)^{-1}
    tail_scale = 4 * (Fraction(3, 4) ** M)
    for a in range(2):
        for b in range(2):
            assert sum_Q_pow[a][b] + tail_scale * Q[a][b] == N_mat[a][b]

    # 5. Exact verification of first-return probabilities f_i^{(m)} = P_i(T_i = m) for m = 1..30
    def exact_fm(state, m):
        if state in (1, 5):
            return Fraction(1, 2) ** m
        if state == 3:
            return Fraction(1) if m == 1 else Fraction(0)
        if state == 2:
            return Fraction(1, 2) if m == 1 else Fraction(1, 2) * (Fraction(1, 4) ** (m - 1))
        if state == 4:
            return Fraction(1, 4) if m == 1 else Fraction(1, 4) * (Fraction(1, 2) ** (m - 1))
        raise ValueError(state)

    for i in range(n):
        state = i + 1
        dist = [Fraction(0)] * n
        for j in range(n):
            if j == i:
                assert P[i][j] == exact_fm(state, 1)
            else:
                dist[j] = P[i][j]
        for m in range(2, M + 1):
            new_dist = [Fraction(0)] * n
            fm = Fraction(0)
            for u in range(n):
                if dist[u] == 0:
                    continue
                for v in range(n):
                    if v == i:
                        fm += dist[u] * P[u][v]
                    else:
                        new_dist[v] += dist[u] * P[u][v]
            assert fm == exact_fm(state, m), f"State {state}, m={m}: {fm} != {exact_fm(state, m)}"
            dist = new_dist

    # 6. Verify invariant distributions on the two closed recurrent classes
    pi1 = [Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0), Fraction(1, 2)]
    pi2 = [Fraction(0), Fraction(0), Fraction(1), Fraction(0), Fraction(0)]
    for pi in (pi1, pi2):
        assert sum(pi) == 1
        assert [sum(pi[i] * P[i][j] for i in range(n)) for j in range(n)] == pi

    print("PASS: Communicating classes = {1, 5} (closed), {3} (closed), {2, 4} (open).")
    print("PASS: Exact rank-1 identity Q^m = (3/4)^{m-1} Q verified for m = 1..30.")
    print("PASS: Exact first-return probabilities f_i^{(m)} verified for all states i in {1..5}, m = 1..30.")
    print("PASS: Recurrent states = {1, 3, 5} (f_1 = f_3 = f_5 = 1); transient states = {2, 4} (f_2 = 2/3, f_4 = 1/2).")
    print(f"PASS: Fundamental matrix on {{2, 4}} is (I - Q)^(-1) = {[[str(x) for x in row] for row in N_mat]}.")


if __name__ == "__main__":
    main()
