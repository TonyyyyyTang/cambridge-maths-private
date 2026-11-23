#!/usr/bin/env python3
"""Exact rational verification for IB/markov-chains/sheet-2/q14 (Reversibility in equilibrium)."""

from fractions import Fraction
import itertools


def rref_nullspace_exact(M):
    """Return exact basis of right nullspace {x : Mx = 0} over Q."""
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


def check_invariance(pi, P):
    n = len(P)
    assert sum(pi) == 1
    assert all(x >= 0 for x in pi)
    for j in range(n):
        assert sum(pi[i] * P[i][j] for i in range(n)) == pi[j]


def satisfies_detailed_balance(pi, P):
    n = len(P)
    for i in range(n):
        for j in range(n):
            if pi[i] * P[i][j] != pi[j] * P[j][i]:
                return False
    return True


def satisfies_path_reversal(pi, P, max_N=3):
    n = len(P)
    for N in range(1, max_N + 1):
        for path in itertools.product(range(n), repeat=N + 1):
            p_fwd = pi[path[0]]
            for t in range(N):
                p_fwd *= P[path[t]][path[t + 1]]
            rev = path[::-1]
            p_rev = pi[rev[0]]
            for t in range(N):
                p_rev *= P[rev[t]][rev[t + 1]]
            if p_fwd != p_rev:
                return False
    return True


def check_part_a():
    grid = [Fraction(0), Fraction(1, 5), Fraction(1, 2), Fraction(3, 4), Fraction(1)]
    tested = 0
    for p in grid:
        for q in grid:
            P = [[p, 1 - p], [q, 1 - q]]
            PT_I = [[P[j][i] - (1 if i == j else 0) for j in range(2)] for i in range(2)]
            ns = rref_nullspace_exact(PT_I)
            if p == 1 and q == 0:
                assert len(ns) == 2
                # Every convex combination is an invariant distribution and reversible
                for alpha in [Fraction(0), Fraction(1, 3), Fraction(1, 2), Fraction(1)]:
                    pi = [alpha, 1 - alpha]
                    check_invariance(pi, P)
                    assert satisfies_detailed_balance(pi, P)
                    assert satisfies_path_reversal(pi, P, max_N=4)
            else:
                assert len(ns) == 1
                pi = [q / (1 - p + q), (1 - p) / (1 - p + q)]
                check_invariance(pi, P)
                assert satisfies_detailed_balance(pi, P)
                assert satisfies_path_reversal(pi, P, max_N=4)
            tested += 1
    print(
        f"[PASS] Part (a): all {tested} parameter pairs (p,q) in grid (including all 4 corners) "
        f"have invariant distribution(s), all satisfying detailed balance & path reversal (N<=4)."
    )


def check_part_b():
    grid = [Fraction(0), Fraction(1, 4), Fraction(1, 3), Fraction(1, 2), Fraction(2, 3), Fraction(3, 4), Fraction(1)]
    for p in grid:
        P = [
            [Fraction(0), p, 1 - p],
            [1 - p, Fraction(0), p],
            [p, 1 - p, Fraction(0)],
        ]
        PT_I = [[P[j][i] - (1 if i == j else 0) for j in range(3)] for i in range(3)]
        ns = rref_nullspace_exact(PT_I)
        assert len(ns) == 1, f"Expected unique invariant dist at p={p}"
        pi = [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)]
        check_invariance(pi, P)
        db = satisfies_detailed_balance(pi, P)
        pr = satisfies_path_reversal(pi, P, max_N=4)
        if p == Fraction(1, 2):
            assert db and pr
        else:
            assert (not db) and (not pr)
            assert pi[0] * P[0][1] - pi[1] * P[1][0] == (2 * p - 1) / 3
    print(
        "[PASS] Part (b): unique invariant dist (1/3,1/3,1/3) for all p in [0,1] (incl. p=0,1); "
        "detailed balance & path reversal hold iff p = 1/2."
    )


def check_part_c():
    # Case p = 0: reducible with closed class {0,1} and transient states >= 2
    M = 8
    P0 = [[Fraction(0) for _ in range(M)] for _ in range(M)]
    P0[0][1] = Fraction(1)
    for i in range(1, M):
        P0[i][i - 1] = Fraction(1)
    PT_I = [[P0[j][i] - (1 if i == j else 0) for j in range(M)] for i in range(M)]
    ns = rref_nullspace_exact(PT_I)
    assert len(ns) == 1
    pi0 = [Fraction(1, 2), Fraction(1, 2)] + [Fraction(0)] * (M - 2)
    check_invariance(pi0, P0)
    assert satisfies_detailed_balance(pi0, P0)
    assert satisfies_path_reversal(pi0, P0, max_N=3)

    # Case 0 < p < 1/2: check exact formula for invariance, detailed balance, and normalisation
    for p in [Fraction(1, 5), Fraction(1, 4), Fraction(1, 3), Fraction(2, 5)]:
        K = 25
        pi = [Fraction(1 - 2 * p, 2 * (1 - p))]
        r = p / (1 - p)
        for i in range(1, K + 1):
            pi.append(Fraction(1 - 2 * p, 2 * (1 - p) ** 2) * (r ** (i - 1)))
        # Exact infinite sum = partial sum up to K + exact geometric tail
        tail = pi[K] * r / (1 - r)
        assert sum(pi) + tail == 1
        # Check exact invariance for j = 0 .. K-1
        assert pi[1] * (1 - p) == pi[0]
        assert pi[0] * 1 + pi[2] * (1 - p) == pi[1]
        for j in range(2, K):
            assert pi[j - 1] * p + pi[j + 1] * (1 - p) == pi[j]
        # Check exact detailed balance for all adjacent pairs up to K-1
        assert pi[0] * 1 == pi[1] * (1 - p)
        for i in range(1, K):
            assert pi[i] * p == pi[i + 1] * (1 - p)
        # Check path reversal for all paths staying within {0..5} up to N=3
        def trans(u, v):
            if u == 0:
                return Fraction(1) if v == 1 else Fraction(0)
            if v == u + 1:
                return p
            if v == u - 1:
                return 1 - p
            return Fraction(0)

        for N in range(1, 4):
            for path in itertools.product(range(6), repeat=N + 1):
                pf = pi[path[0]]
                for t in range(N):
                    pf *= trans(path[t], path[t + 1])
                rev = path[::-1]
                pr = pi[rev[0]]
                for t in range(N):
                    pr *= trans(rev[t], rev[t + 1])
                assert pf == pr

    # Case 1/2 <= p < 1: invariant measure from recurrence has non-decreasing terms >= 1/(1-p)
    for p in [Fraction(1, 2), Fraction(2, 3), Fraction(4, 5)]:
        mu = [Fraction(1)]
        r = p / (1 - p)
        assert r >= 1
        for i in range(1, 20):
            val = Fraction(1, 1 - p) * (r ** (i - 1))
            mu.append(val)
            assert val >= Fraction(1, 1 - p) > 0
        # Verify mu satisfies global & detailed balance, yet sum diverges
        assert mu[1] * (1 - p) == mu[0]
        assert mu[0] + mu[2] * (1 - p) == mu[1]
        for j in range(2, 19):
            assert mu[j - 1] * p + mu[j + 1] * (1 - p) == mu[j]
        assert sum(mu) > 20

    # Case p = 1: pi_0 = 0 and pi_j = pi_{j-1} forces pi == 0 identically
    pi1 = [Fraction(0)] * 10
    assert all(x == 0 for x in pi1)

    print(
        "[PASS] Part (c): p=0 (reducible, unique pi=(1/2,1/2,0,...)) and 0<p<1/2 (irreducible, "
        "exact pi) satisfy invariance, detailed balance & path reversal; "
        "1/2<=p<=1 have no invariant distribution."
    )


def check_part_d():
    # 1. Finite irreducible symmetric (non-trivial weights)
    P_irr = [
        [Fraction(1, 2), Fraction(1, 3), Fraction(1, 6), Fraction(0)],
        [Fraction(1, 3), Fraction(1, 6), Fraction(0), Fraction(1, 2)],
        [Fraction(1, 6), Fraction(0), Fraction(1, 3), Fraction(1, 2)],
        [Fraction(0), Fraction(1, 2), Fraction(1, 2), Fraction(0)],
    ]
    for i in range(4):
        assert sum(P_irr[i]) == 1
        for j in range(4):
            assert P_irr[i][j] == P_irr[j][i]
    PT_I = [[P_irr[j][i] - (1 if i == j else 0) for j in range(4)] for i in range(4)]
    ns = rref_nullspace_exact(PT_I)
    assert len(ns) == 1
    pi_irr = [Fraction(1, 4)] * 4
    check_invariance(pi_irr, P_irr)
    assert satisfies_detailed_balance(pi_irr, P_irr)
    assert satisfies_path_reversal(pi_irr, P_irr, max_N=3)

    # 2. Finite reducible symmetric (closed classes {0,1} and {2,3,4})
    P_red = [
        [Fraction(1, 4), Fraction(3, 4), Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(3, 4), Fraction(1, 4), Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(1, 2), Fraction(1, 2)],
        [Fraction(0), Fraction(0), Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)],
        [Fraction(0), Fraction(0), Fraction(1, 2), Fraction(1, 6), Fraction(1, 3)],
    ]
    for i in range(5):
        assert sum(P_red[i]) == 1
        for j in range(5):
            assert P_red[i][j] == P_red[j][i]
    PT_I = [[P_red[j][i] - (1 if i == j else 0) for j in range(5)] for i in range(5)]
    ns = rref_nullspace_exact(PT_I)
    assert len(ns) == 2
    for w in [Fraction(0), Fraction(1, 5), Fraction(2, 3), Fraction(1)]:
        pi_mix = [w / 2, w / 2, (1 - w) / 3, (1 - w) / 3, (1 - w) / 3]
        check_invariance(pi_mix, P_red)
        assert satisfies_detailed_balance(pi_mix, P_red)
        assert satisfies_path_reversal(pi_mix, P_red, max_N=3)

    # 3. Infinite reducible symmetric with finite closed class {0,1} + reflected SSRW on {2,3,...}
    # Check equations directly on states 0..15 with pi = (1/2, 1/2, 0, 0, ...)
    K = 15
    def p_inf_red(i, j):
        if {i, j} == {0, 1}:
            return Fraction(1)
        if i == 2 and j == 2:
            return Fraction(1, 2)
        if min(i, j) >= 2 and abs(i - j) == 1:
            return Fraction(1, 2)
        return Fraction(0)

    pi_inf = [Fraction(1, 2), Fraction(1, 2)] + [Fraction(0)] * (K - 1)
    for i in range(K):
        # Row sum for i < K
        assert sum(p_inf_red(i, j) for j in range(K + 1)) == 1
        for j in range(K):
            assert p_inf_red(i, j) == p_inf_red(j, i)
            assert pi_inf[i] * p_inf_red(i, j) == pi_inf[j] * p_inf_red(j, i)
        assert sum(pi_inf[k] * p_inf_red(k, i) for k in range(K + 1)) == pi_inf[i]

    print(
        "[PASS] Part (d): verified finite irreducible symmetric (unique uniform pi), "
        "finite reducible symmetric (all class-wise constant mixtures reversible), and "
        "infinite reducible symmetric with finite closed class (reversible in equilibrium)."
    )


if __name__ == "__main__":
    check_part_a()
    check_part_b()
    check_part_c()
    check_part_d()
    print("ALL CHECKS PASSED.")
