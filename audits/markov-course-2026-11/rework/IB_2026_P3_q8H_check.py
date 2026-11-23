#!/usr/bin/env python3
"""Exact rational verification for IB/2026/P3/q8H."""

from fractions import Fraction


def solve_linear_system(A, b):
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        pivot = next(r for r in range(col, n) if M[r][col] != 0)
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


def mat_mul(A, B):
    n = len(A)
    return [
        [sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def main():
    P = [
        [Fraction(0), Fraction(1, 3), Fraction(1, 3), Fraction(1, 6), Fraction(0), Fraction(1, 6)],
        [Fraction(1, 3), Fraction(2, 3), Fraction(0), Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(1, 3), Fraction(1, 3), Fraction(0), Fraction(1, 3), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(1, 2), Fraction(1, 2), Fraction(0)],
        [Fraction(1, 2), Fraction(0), Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0)],
    ]
    n = 6

    # 1. Verify stochasticity
    for i, row in enumerate(P):
        assert all(0 <= x <= 1 for x in row), f"Entry out of [0,1] in row {i+1}"
        assert sum(row) == 1, f"Row {i+1} sums to {sum(row)}"

    # 2. Compute exact reachability (transitive closure)
    reach = [[i == j or P[i][j] > 0 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if reach[i][k] and reach[k][j]:
                    reach[i][j] = True

    visited = set()
    classes = []
    for i in range(n):
        if i in visited:
            continue
        cls = {j + 1 for j in range(n) if reach[i][j] and reach[j][i]}
        for s in cls:
            visited.add(s - 1)
        is_closed = all(
            (j + 1 in cls) for s in cls for j in range(n) if P[s - 1][j] > 0
        )
        classes.append((sorted(cls), is_closed))

    assert classes == [([1, 2, 3, 6], False), ([4, 5], True)]
    print("Communicating classes and closedness:")
    for cls, is_closed in classes:
        print(f"  {cls}: {'closed' if is_closed else 'not closed'}")

    # 3. Solve hitting probabilities h_i = P_i(H_2 < inf)
    # Closed class {4, 5} does not contain 2 => h_4 = h_5 = 0; h_2 = 1.
    # Solve (I - Q) h_T = r for T = {1, 3, 6} (indices 0, 2, 5)
    T = [0, 2, 5]
    A = [[Fraction(1 if a == b else 0) - P[T[a]][T[b]] for b in range(3)] for a in range(3)]
    b_hit2 = [P[T[a]][1] for a in range(3)]
    h_T = solve_linear_system(A, b_hit2)
    h = [h_T[0], Fraction(1), h_T[1], Fraction(0), Fraction(0), h_T[2]]
    assert h == [Fraction(13, 22), Fraction(1), Fraction(1, 2), Fraction(0), Fraction(0), Fraction(6, 11)]

    # Verify full harmonic equations h_i = sum_j P[i][j] h_j for all i != 2
    for i in range(n):
        if i == 1:
            assert h[i] == 1
        else:
            assert h[i] == sum(P[i][j] * h[j] for j in range(n))

    # Also compute probability of hitting {4, 5} before 2 from T = {1, 3, 6}
    b_hit45 = [P[T[a]][3] + P[T[a]][4] for a in range(3)]
    g_T = solve_linear_system(A, b_hit45)
    assert [h_T[a] + g_T[a] for a in range(3)] == [Fraction(1), Fraction(1), Fraction(1)]

    # Check general non-negative solution family h_4 = h_5 = c >= 0
    for c in [Fraction(0), Fraction(1, 3), Fraction(1)]:
        hc = [
            (13 + 9 * c) / 22,
            Fraction(1),
            (1 + c) / 2,
            c,
            c,
            (12 + 10 * c) / 22,
        ]
        for i in range(n):
            if i != 1:
                assert hc[i] == sum(P[i][j] * hc[j] for j in range(n))

    # 4. Independent check: powers of modified chain with state 2 absorbing
    P_abs = [list(row) for row in P]
    P_abs[1] = [Fraction(1 if j == 1 else 0) for j in range(n)]
    M = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    for step in range(1, 61):
        M = mat_mul(M, P_abs)
        if step in (5, 10, 20, 40, 60):
            err = h[0] - M[0][1]
            assert err > 0
            print(f"  After {step:2d} steps: P_1(H_2 <= {step}) = {float(M[0][1]):.12f} (gap to 13/22: {float(err):.2e})")

    assert h[0] - M[0][1] < Fraction(1, 10**10)
    print(f"Exact hitting probabilities: h = {[str(x) for x in h]}")
    print(f"Exact P_1(ever reaches 2) = h_1 = {h[0]} = {float(h[0]):.12f}")
    print("ALL CHECKS PASSED for IB/2026/P3/q8H.")


if __name__ == "__main__":
    main()
