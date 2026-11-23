"""Exact rational verification for IB/2015/P3/q9H."""
from fractions import Fraction

# States 1, 2, 3, 4 indexed 0, 1, 2, 3
P = [
    [Fraction(1, 2), Fraction(0), Fraction(1, 2), Fraction(0)],
    [Fraction(0), Fraction(1, 2), Fraction(0), Fraction(1, 2)],
    [Fraction(1, 2), Fraction(0), Fraction(1, 2), Fraction(0)],
    [Fraction(1, 2), Fraction(0), Fraction(0), Fraction(1, 2)],
]


def mat_mul(A, B):
    n = len(A)
    return [
        [sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def main():
    # 1. Verify P is stochastic
    for i, row in enumerate(P):
        assert sum(row) == 1
        assert all(x >= 0 for x in row)

    # 2. Reachability and communicating classes
    reach = [[i == j or P[i][j] > 0 for j in range(4)] for i in range(4)]
    for k in range(4):
        for i in range(4):
            for j in range(4):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])

    visited = set()
    classes = []
    for i in range(4):
        if i in visited:
            continue
        cls = {j + 1 for j in range(4) if reach[i][j] and reach[j][i]}
        for s in cls:
            visited.add(s - 1)
        is_closed = all(
            not reach[a - 1][b - 1]
            for a in cls
            for b in range(1, 5)
            if b not in cls
        )
        classes.append((sorted(cls), is_closed))

    classes.sort()
    assert classes == [([1, 3], True), ([2], False), ([4], False)]
    print("Communicating classes (states, is_closed):", classes)

    # 3. Distribution of N = inf{n >= 1 : X_n != 2} given X_0 = 2
    p22 = P[1][1]
    for n in range(1, 21):
        prob_N_n = (p22 ** (n - 1)) * (1 - p22)
        assert prob_N_n == Fraction(1, 2 ** n)
    # Exact expectation via first-step equation E(N) = 1 + p22 * E(N)
    E_N = Fraction(1) / (Fraction(1) - p22)
    assert E_N == 2
    print("P(N=n) = 2^{-n} verified for n=1..20; E(N) =", E_N)

    # 4. Exact n-step distribution from X_0 = 2 via exact matrix powers P^n
    Pn = [[Fraction(1 if i == j else 0) for j in range(4)] for i in range(4)]
    assert Pn[1] == [Fraction(0), Fraction(1), Fraction(0), Fraction(0)]
    for n in range(1, 21):
        Pn = mat_mul(Pn, P)
        row2 = Pn[1]
        expected = [
            Fraction(1, 2) - Fraction(1, 2 ** n),
            Fraction(1, 2 ** n),
            Fraction(1, 2) - Fraction(n, 2 ** n),
            Fraction(n, 2 ** n),
        ]
        assert row2 == expected, f"Mismatch at n={n}: {row2} != {expected}"
    print("Exact row 2 of P^n verified for n=1..20: (1/2 - 2^{-n}, 2^{-n}, 1/2 - n*2^{-n}, n*2^{-n}).")

    # 5. Verify joint law of (N, M, X_{N+M}, X_{N+M+1})
    for k in range(1, 7):
        for m in range(1, 7):
            for x_after in [1, 3]:
                # Path: 2 (k times: indices 0..k-1), 4 (m times: indices k..k+m-1), 1 (at k+m), x_after (at k+m+1)
                prob = (
                    (P[1][1] ** (k - 1))
                    * P[1][3]
                    * (P[3][3] ** (m - 1))
                    * P[3][0]
                    * P[0][x_after - 1]
                )
                assert prob == Fraction(1, 2 ** k) * Fraction(1, 2 ** m) * Fraction(1, 2)
    print("Verified independent Geom(1/2) holding times N (at 2), M (at 4), entry at 1, then i.i.d. Unif{1,3}.")

    # 6. Verify unique invariant distribution pi = (1/2, 0, 1/2, 0)
    pi = [Fraction(1, 2), Fraction(0), Fraction(1, 2), Fraction(0)]
    for j in range(4):
        assert sum(pi[i] * P[i][j] for i in range(4)) == pi[j]
    print("Verified invariant distribution pi = (1/2, 0, 1/2, 0).")
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
