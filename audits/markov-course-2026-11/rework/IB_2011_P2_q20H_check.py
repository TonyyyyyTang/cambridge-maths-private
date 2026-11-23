"""Exact rational verification for IB/2011/P2/q20H (flea on a polyhedron)."""
from fractions import Fraction


def poly_mul(a, b):
    res = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        for j, cb in enumerate(b):
            res[i + j] += ca * cb
    return res


def poly_eval(a, x):
    val = Fraction(0)
    xp = Fraction(1)
    for c in a:
        val += c * xp
        xp *= x
    return val


def poly_deriv_at_1(a):
    return sum(Fraction(k) * c for k, c in enumerate(a))


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


def verify_part_i_general():
    """Verify part (i) on a general finite Markov chain, including states where H=inf has positive prob."""
    F = Fraction
    # States {0, 1, 2, 3}, A = {0}, state 3 is absorbing outside A (so H=inf from 3)
    P = [
        [F(1, 2), F(1, 2), F(0), F(0)],      # 0 in A
        [F(1, 3), F(1, 3), F(1, 3), F(0)],   # 1
        [F(0), F(1, 2), F(0), F(1, 2)],      # 2
        [F(0), F(0), F(0), F(1)],            # 3 (never hits A)
    ]
    A_set = {0}
    # Case 1: g(i) < 1 on S \ A -> on {H=inf}, infinite product is 0, so V_3 = 0
    g1 = [F(1, 2), F(3, 4), F(2, 5), F(1, 3)]
    # Solve linear system: V_0 = 1, V_i = g(i) sum_j P_{ij} V_j for i not in A
    M = [[F(0)] * 4 for _ in range(4)]
    rhs = [F(0)] * 4
    for i in range(4):
        if i in A_set:
            M[i][i] = F(1)
            rhs[i] = F(1)
        else:
            M[i][i] = F(1)
            for j in range(4):
                M[i][j] -= g1[i] * P[i][j]
    V1 = solve_linear(M, rhs)
    assert V1[0] == 1 and V1[3] == 0
    # Verify against truncated expectation E[prod_{n=0}^{min(H,T)-1} g(X_n)] as T grows
    # Recursion: W_i^{(0)} = 1 (upper bound since g<=1), W_i^{(T+1)} = 1 if i in A else g(i) sum_j P_{ij} W_j^{(T)}
    W = [F(1)] * 4
    for _ in range(40):
        W = [F(1) if i in A_set else g1[i] * sum(P[i][j] * W[j] for j in range(4)) for i in range(4)]
    for i in range(4):
        assert abs(W[i] - V1[i]) < F(1, 10**6)

    # Case 2: g(3) = 1 at absorbing state 3 outside A -> on {H=inf} at state 3, prod = 1!
    # The equation V_3 = 1 * 1 * V_3 is satisfied by V_3 = 1 (and note non-uniqueness if g=1 on closed class outside A!).
    g2 = [F(1, 2), F(3, 4), F(2, 5), F(1)]
    W2 = [F(1)] * 4
    for _ in range(40):
        W2 = [F(1) if i in A_set else g2[i] * sum(P[i][j] * W2[j] for j in range(4)) for i in range(4)]
    # Exact V2 with V2[3] = 1:
    M2 = [[F(1), F(0), F(0), F(0)],
          [-g2[1] * P[1][0], F(1) - g2[1] * P[1][1], -g2[1] * P[1][2], -g2[1] * P[1][3]],
          [-g2[2] * P[2][0], -g2[2] * P[2][1], F(1) - g2[2] * P[2][2], -g2[2] * P[2][3]],
          [F(0), F(0), F(0), F(1)]]
    V2 = solve_linear(M2, [F(1), F(0), F(0), F(1)])
    for i in range(4):
        assert abs(W2[i] - V2[i]) < F(1, 10**6)
        if i in A_set:
            assert V2[i] == 1
        else:
            assert V2[i] == g2[i] * sum(P[i][j] * V2[j] for j in range(4))
    print(f"[PASS] Part (i) general chain: V1={V1} (g(3)<1), V2={V2} (g(3)=1).")


def verify_flea(N):
    """Exact polynomial and transition-matrix verification for the flea on N vertices."""
    F = Fraction
    # Transition matrix P (0-indexed for vertices 1..N)
    P = [[F(0)] * N for _ in range(N)]
    P[0][0] = F(1)  # vertex 1 is absorbing
    for i in range(2, N + 1):
        for j in range(1, i):
            P[i - 1][j - 1] = F(1, i - 1)
    for row in P:
        assert sum(row) == 1

    # Exact PGF G_i(s) = E(s^H | X_0 = i) as polynomial in s via direct DP over transition probabilities:
    # G_1(s) = 1
    # G_i(s) = s * sum_{j=1}^{i-1} P_{ij} G_j(s)  for i >= 2
    G = [None] * (N + 1)
    G[1] = [F(1)]
    for i in range(2, N + 1):
        # Multiply each G[j] by s / (i-1) and sum
        poly = [F(0)] * i
        for j in range(1, i):
            for deg, coeff in enumerate(G[j]):
                poly[deg + 1] += P[i - 1][j - 1] * coeff
        G[i] = poly

    # Check every G[i] is a valid PGF (non-negative coeffs summing to 1)
    for i in range(1, N + 1):
        assert all(c >= 0 for c in G[i])
        assert sum(G[i]) == 1

    # Check recurrence (i/s) V_{i+1} = V_i + ((i-1)/s) V_i  <=>  i * G_{i+1}(s) = (s + i - 1) * G_i(s)
    for i in range(1, N):
        lhs = [F(i) * c for c in G[i + 1]]
        rhs = poly_mul([F(i - 1), F(1)], G[i])
        assert lhs == rhs, f"Recurrence failed at N={N}, i={i}: {lhs} vs {rhs}"

    # Check closed-form product: prod_{i=1}^{N-1} (1 + (s-1)/i) = prod_{i=1}^{N-1} ((i-1)/i + (1/i) s)
    prod_poly = [F(1)]
    for i in range(1, N):
        prod_poly = poly_mul(prod_poly, [F(i - 1, i), F(1, i)])
    assert G[N] == prod_poly, f"Product mismatch at N={N}: {G[N]} vs {prod_poly}"

    # Check at concrete rational values of s in (0, 1] via direct linear system from part (i)
    for s in [F(1, 5), F(1, 2), F(3, 4), F(1)]:
        # Solve V_1 = 1, V_i = s * sum_j P_{ij} V_j
        M = [[F(0)] * N for _ in range(N)]
        rhs = [F(0)] * N
        M[0][0] = F(1)
        rhs[0] = F(1)
        for i in range(1, N):
            M[i][i] = F(1)
            for j in range(N):
                M[i][j] -= s * P[i][j]
        V = solve_linear(M, rhs)
        assert V[0] == 1
        for i in range(1, N):
            assert (F(i) / s) * V[i] == V[i - 1] + (F(i - 1) / s) * V[i - 1]
        expected = F(1)
        for i in range(1, N):
            expected *= F(1) + (s - F(1)) / F(i)
        assert V[N - 1] == expected == poly_eval(G[N], s)

    # Check part (iv): E(H | X_0 = N) = sum_{i=1}^{N-1} 1/i
    # Compare three ways: derivative of PGF at 1, direct hitting-time linear solve, harmonic sum
    E_from_pgf = poly_deriv_at_1(G[N])
    harmonic = sum((F(1, i) for i in range(1, N)), F(0))
    # Direct hitting-time solve: h_1 = 0, h_i = 1 + sum_{j=1}^{i-1} P_{ij} h_j
    h_direct = [F(0)] * (N + 1)
    for i in range(2, N + 1):
        h_direct[i] = F(1) + sum(P[i - 1][j - 1] * h_direct[j] for j in range(1, i))
    assert E_from_pgf == harmonic == h_direct[N]
    return G[N], harmonic


def main():
    verify_part_i_general()
    for N in range(1, 13):
        poly, harm = verify_flea(N)
        if N <= 6 or N == 12:
            print(f"[PASS] Flea N={N:2d}: E(H|X_0=N) = {harm} (~{float(harm):.6f}), deg(PGF)={len(poly)-1}")
    print("All IB/2011/P2/q20H checks passed.")


if __name__ == "__main__":
    main()
