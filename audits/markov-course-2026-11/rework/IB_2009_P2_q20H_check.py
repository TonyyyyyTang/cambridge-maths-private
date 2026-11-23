"""Exact rational verification for IB/2009/P2/q20H (superharmonic functions & recurrence)."""
from fractions import Fraction


def mat_vec(P, v):
    return [sum(P[i][j] * v[j] for j in range(len(v))) for i in range(len(P))]


def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def solve_linear(A, b):
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        pivot = None
        for r in range(col, n):
            if M[r][col] != 0:
                pivot = r
                break
        assert pivot is not None
        M[col], M[pivot] = M[pivot], M[col]
        piv = M[col][col]
        for j in range(col, n + 1):
            M[col][j] /= piv
        for r in range(n):
            if r != col and M[r][col] != 0:
                factor = M[r][col]
                for j in range(col, n + 1):
                    M[r][j] -= factor * M[col][j]
    return [M[i][n] for i in range(n)]


def make_tilde_P(P, B):
    n = len(P)
    P_tilde = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        if i in B:
            P_tilde[i][i] = Fraction(1)
        else:
            for j in range(n):
                P_tilde[i][j] = P[i][j]
    return P_tilde


def test_reducible_chain():
    """Chain on {0,1,2,3,4}:
    - 0 is in B (can jump to 1 or 4, showing strict inequality h_0 > (Ph)_0!)
    - 1, 2 are transient states communicating with B={0} and absorbing sink {4}
    - 3 is in B (also communicates with 2)
    - 4 is an absorbing state outside B (so h_4 = 0)
    """
    F = Fraction
    P = [
        [F(1, 2), F(1, 4), F(0), F(0), F(1, 4)],  # 0 in B
        [F(1, 2), F(0), F(1, 2), F(0), F(0)],     # 1
        [F(0), F(1, 3), F(0), F(1, 3), F(1, 3)],  # 2
        [F(0), F(0), F(1, 2), F(1, 2), F(0)],     # 3 in B
        [F(0), F(0), F(0), F(0), F(1)],           # 4 absorbing sink
    ]
    n = len(P)
    for row in P:
        assert sum(row) == 1
    B = {0, 3}
    P_tilde = make_tilde_P(P, B)

    # Exact h_i = P_i(tau < inf): h_0 = h_3 = 1, h_4 = 0, solve for h_1, h_2
    # h_1 = 1/2 * 1 + 1/2 * h_2
    # h_2 = 1/3 * h_1 + 1/3 * 1 + 1/3 * 0
    A = [[F(1), F(-1, 2)], [F(-1, 3), F(1)]]
    rhs = [F(1, 2), F(1, 3)]
    h1, h2 = solve_linear(A, rhs)
    h = [F(1), h1, h2, F(1), F(0)]
    assert h == [F(1), F(4, 5), F(3, 5), F(1), F(0)]

    # Check (a) and (b) for h
    Ph = mat_vec(P, h)
    for i in range(n):
        assert h[i] >= Ph[i] >= 0
        if i in B:
            assert h[i] == 1
        else:
            assert h[i] == Ph[i]
    # Note strict inequality at i=0 and i=3!
    assert h[0] > Ph[0] and h[3] > Ph[3]

    # Check P_tilde^n * 1_B increases monotonically toward h
    ind_B = [F(1) if i in B else F(0) for i in range(n)]
    Pn = [[F(1) if i == j else F(0) for j in range(n)] for i in range(n)]
    prev = [F(0)] * n
    for step in range(25):
        prob_le_step = mat_vec(Pn, ind_B)
        for i in range(n):
            assert prev[i] <= prob_le_step[i] <= h[i]
        prev = prob_le_step
        Pn = mat_mul(Pn, P_tilde)
    # Check close convergence by step 24
    for i in range(n):
        assert h[i] - prev[i] < F(1, 10**6)

    # Test several superharmonic solutions g satisfying (a) and (b)
    # g1 = [1, 1, 1, 1, 1] (constant 1)
    # g2 = h + [0, 0, 0, 0, 1] = [1, 4/5, 3/5, 1, 1]? Wait: check (Pg2)_2 = 1/3*(4/5)+1/3*1+1/3*1 = 14/15 > 3/5 (not superharmonic!)
    # Construct genuine superharmonic g > h: take g = h + c * (1 - h) + d, e.g. convex combination / adding positive constant on sink?
    # For any alpha >= 0, g = h + alpha * (1 - h) satisfies:
    #   for i not in B: (Pg)_i = (Ph)_i + alpha*(1 - (Ph)_i) = h_i + alpha*(1 - h_i) = g_i
    #   for i in B: g_i = 1 >= (Pg)_i (provided alpha <= 1, so g_j <= 1!)
    # Also g = h + [0, 1/10, 1/20, 0, 0]? Check!
    test_gs = [
        [F(1)] * n,
        [h[i] + F(1, 2) * (F(1) - h[i]) for i in range(n)],
        [h[i] + F(1, 5) * (F(1) - h[i]) for i in range(n)],
    ]
    # Also test strictly superharmonic solutions (including one exceeding 1 outside B!):
    g_strict_1 = [F(1), F(19, 20), F(17, 20), F(1), F(1, 2)]
    g_exceed_1 = [F(1), F(5, 4), F(1), F(1), F(1, 2)]
    test_gs.extend([g_strict_1, g_exceed_1])

    for idx, g in enumerate(test_gs):
        Pg = mat_vec(P, g)
        Ptg = mat_vec(P_tilde, g)
        for i in range(n):
            assert g[i] >= Pg[i] >= 0, f"g[{idx}] fails (a) at {i}"
            if i in B:
                assert g[i] == 1, f"g[{idx}] fails (b) at {i}"
            assert g[i] >= Ptg[i], f"g[{idx}] fails (c) at {i}"
            assert g[i] >= h[i], f"Minimality fails for g[{idx}] at {i}"
        # Check g >= P_tilde^k g >= P_tilde^k 1_B for k=0..10
        cur = list(g)
        cur_B = list(ind_B)
        for _ in range(10):
            for i in range(n):
                assert g[i] >= cur[i] >= cur_B[i]
            cur = mat_vec(P_tilde, cur)
            cur_B = mat_vec(P_tilde, cur_B)

    print(
        f"[PASS] Reducible chain: h={h}, Ph={Ph} (strict at B: h_0={h[0]}>{Ph[0]}, h_3={h[3]}>{Ph[3]}), "
        f"all {len(test_gs)} superharmonic g satisfy (c) and g >= h."
    )


def test_irreducible_finite_chain():
    """On a finite irreducible (hence recurrent) chain, verify that (P - I)g <= 0 forces g constant.
    Indeed, if pi > 0 satisfies pi P = pi, then pi (g - Pg) = 0 with g - Pg >= 0 forces g = Pg,
    and ker(I - P) has dimension 1 (spanned by 1).
    """
    F = Fraction
    P = [
        [F(0), F(1, 2), F(1, 2), F(0)],
        [F(1, 3), F(0), F(1, 3), F(1, 3)],
        [F(0), F(1, 2), F(0), F(1, 2)],
        [F(1, 4), F(0), F(3, 4), F(0)],
    ]
    n = len(P)
    # Find unique invariant distribution pi
    # (P^T - I) pi = 0 with sum pi = 1
    A = [[P[j][i] - (F(1) if i == j else F(0)) for j in range(n)] for i in range(n - 1)]
    A.append([F(1)] * n)
    b = [F(0)] * (n - 1) + [F(1)]
    pi = solve_linear(A, b)
    assert all(x > 0 for x in pi)
    # Check pi P = pi
    for j in range(n):
        assert sum(pi[i] * P[i][j] for i in range(n)) == pi[j]

    # Verify rank(I - P) == n - 1 by solving (I - P)g = 0 with g_0 = 1 -> unique solution g = [1,1,1,1]
    A_harm = [[(F(1) if i == j else F(0)) - P[i][j] for j in range(n)] for i in range(1, n)]
    A_harm.insert(0, [F(1)] + [F(0)] * (n - 1))
    g_sol = solve_linear(A_harm, [F(1)] + [F(0)] * (n - 1))
    assert g_sol == [F(1)] * n
    print(f"[PASS] Finite irreducible chain: pi={pi}, unique harmonic solution with g_0=1 is {g_sol}.")


def test_transient_infinite_chain():
    """Verify non-constant solution h_i = (q/p)^i on transient biased random walk on Z+ (p=2/3, q=1/3)."""
    F = Fraction
    p, q = F(2, 3), F(1, 3)
    # p_{00} = q, p_{01} = p; for i >= 1: p_{i,i+1} = p, p_{i,i-1} = q
    # B = {0}, h_i = (q/p)^i = (1/2)^i
    h = [F(1, 2**i) for i in range(30)]
    Ph0 = q * h[0] + p * h[1]
    assert h[0] == 1 and Ph0 == F(2, 3) < h[0] and Ph0 >= 0  # P_0(T_0 < inf) = 2/3 < 1
    for i in range(1, 29):
        Phi = q * h[i - 1] + p * h[i + 1]
        assert h[i] == Phi >= 0
    assert h[0] != h[1]  # non-constant solution to (a)!
    print(
        f"[PASS] Transient biased RW (p=2/3, q=1/3): non-constant solution h_i=(1/2)^i satisfies (a), "
        f"with (Ph)_0 = P_0(T_0 < inf) = {Ph0} < 1 = h_0."
    )


def main():
    test_reducible_chain()
    test_irreducible_finite_chain()
    test_transient_infinite_chain()
    print("All IB/2009/P2/q20H checks passed.")


if __name__ == "__main__":
    main()
