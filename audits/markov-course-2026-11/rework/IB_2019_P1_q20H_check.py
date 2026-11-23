"""Exact rational verification for IB/2019/P1/q20H.

Verifies parts (a), (b), (c) across multiple exact finite irreducible Markov chains
(N = 1, 2, 3, 4, 5, including reversible, non-reversible aperiodic, and periodic chains):
  (a) P*(x, y) = P[X_0 = y | X_1 = x] = pi(y) P(y, x) / pi(x), P* stochastic, pi P* = pi.
  (b) E[T] = N where T = min{n >= 1 : X_n = X_0} and X_0 ~ pi.
  (c) E[tau] = N - 1 where tau = min{n >= 0 : C_n = R}, verified BOTH via exact hitting-time
      linear solves AND via direct distribution check P(tau = k) = P(T = k + 1) for k = 0..15.
  Also checks contrasts:
    - If strict n >= 1 were required (tau_+ = min{n >= 1 : C_n = R}), E[tau_+] = N - 1 + tr(P).
    - If the robber stepped via P instead of P* on a non-reversible chain, E[tau] != N - 1.
"""

from fractions import Fraction


def solve_linear(A, b):
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


def invariant_distribution(P):
    n = len(P)
    # Replace last equation of (P^T - I) pi = 0 by sum(pi) = 1
    A = [
        [P[j][i] - (Fraction(1) if i == j else Fraction(0)) for j in range(n)]
        for i in range(n)
    ]
    A[-1] = [Fraction(1)] * n
    b = [Fraction(0)] * (n - 1) + [Fraction(1)]
    pi = solve_linear(A, b)
    assert sum(pi) == 1
    assert all(x > 0 for x in pi)
    for j in range(n):
        assert sum(pi[i] * P[i][j] for i in range(n)) == pi[j]
    return pi


def hitting_times_to_target(P, target):
    """Exact h(x, target) = E_x[min{n >= 0 : X_n = target}] for all x."""
    n = len(P)
    others = [x for x in range(n) if x != target]
    if not others:
        return [Fraction(0)]
    m = len(others)
    idx = {x: k for k, x in enumerate(others)}
    A = [[Fraction(0)] * m for _ in range(m)]
    b = [Fraction(1)] * m
    for x in others:
        r = idx[x]
        A[r][r] = Fraction(1)
        for y in others:
            A[r][idx[y]] -= P[x][y]
    sol = solve_linear(A, b)
    h = [Fraction(0)] * n
    for x in others:
        h[x] = sol[idx[x]]
    return h


def exact_return_and_catch_pmfs(P, pi, P_star, max_k=15):
    """Independently compute P(T = k+1) and P(tau = k) for k = 0..max_k."""
    n = len(P)
    pmf_T = []    # pmf_T[k] = P(T = k + 1)
    pmf_tau = []  # pmf_tau[k] = P(tau = k)

    # For T: condition on X_0 = y, track state distribution at time m >= 1
    # on trajectories that have not yet returned to y.
    for k in range(max_k + 1):
        # Compute P(T = k + 1)
        p_T = Fraction(0)
        for y in range(n):
            # Distribution after 1 step from y
            dist = [P[y][x] for x in range(n)]
            if k == 0:
                p_T += pi[y] * dist[y]
            else:
                # Zero out y (haven't returned at time 1), evolve k-1 more steps avoiding y,
                # then step into y at time k+1
                dist[y] = Fraction(0)
                for _ in range(k - 1):
                    nxt = [Fraction(0)] * n
                    for u in range(n):
                        if dist[u] != 0:
                            for v in range(n):
                                if v != y:
                                    nxt[v] += dist[u] * P[u][v]
                    dist = nxt
                hit = sum(dist[u] * P[u][y] for u in range(n))
                p_T += pi[y] * hit
        pmf_T.append(p_T)

        # Compute P(tau = k) directly from (Z ~ pi, R | Z ~ P*, cop C_0 = Z moves via P)
        p_tau = Fraction(0)
        for z in range(n):
            for r in range(n):
                joint = pi[z] * P_star[z][r]
                if joint == 0:
                    continue
                if k == 0:
                    if z == r:
                        p_tau += joint
                else:
                    if z == r:
                        continue
                    dist = [Fraction(0)] * n
                    dist[z] = Fraction(1)
                    for _ in range(k - 1):
                        nxt = [Fraction(0)] * n
                        for u in range(n):
                            if dist[u] != 0:
                                for v in range(n):
                                    if v != r:
                                        nxt[v] += dist[u] * P[u][v]
                        dist = nxt
                    hit = sum(dist[u] * P[u][r] for u in range(n))
                    p_tau += joint * hit
        pmf_tau.append(p_tau)

    return pmf_T, pmf_tau


def verify_chain(name, P, check_non_rev_contrast=False):
    n = len(P)
    for row in P:
        assert len(row) == n
        assert all(x >= 0 for x in row)
        assert sum(row) == 1

    pi = invariant_distribution(P)

    # Part (a): Direct Bayes calculation vs formula pi(y) P(y, x) / pi(x)
    # Joint table J[y][x] = P(X_0 = y, X_1 = x) = pi[y] * P[y][x]
    J = [[pi[y] * P[y][x] for x in range(n)] for y in range(n)]
    P_X1 = [sum(J[y][x] for y in range(n)) for x in range(n)]
    assert P_X1 == pi

    P_star = [[J[y][x] / P_X1[x] for y in range(n)] for x in range(n)]
    for x in range(n):
        for y in range(n):
            assert P_star[x][y] == pi[y] * P[y][x] / pi[x]
        assert sum(P_star[x]) == 1
    for y in range(n):
        assert sum(pi[x] * P_star[x][y] for x in range(n)) == pi[y]

    # Exact hitting times h[x][y] = E_x[H_y]
    h = [[Fraction(0)] * n for _ in range(n)]
    for y in range(n):
        hy = hitting_times_to_target(P, y)
        for x in range(n):
            h[x][y] = hy[x]

    # Part (b): E_y[T_y] = 1 + sum_x P(y,x) h(x,y) == 1/pi(y), and E[T] = N
    E_T_given = []
    for y in range(n):
        e_ty = Fraction(1) + sum(P[y][x] * h[x][y] for x in range(n))
        assert e_ty == Fraction(1) / pi[y]
        E_T_given.append(e_ty)
    E_T = sum(pi[y] * E_T_given[y] for y in range(n))
    assert E_T == Fraction(n)

    # Part (c): E[tau] = sum_{x,y} pi(x) P*(x,y) h(x,y) == N - 1
    E_tau = sum(pi[x] * P_star[x][y] * h[x][y] for x in range(n) for y in range(n))
    assert E_tau == Fraction(n - 1)

    # Direct distribution identity P(tau = k) == P(T = k + 1) for k = 0..15
    pmf_T, pmf_tau = exact_return_and_catch_pmfs(P, pi, P_star, max_k=15)
    assert pmf_T == pmf_tau

    # Contrast 1: strict positive catch time tau_+ = min{n >= 1 : C_n = R}
    tr_P = sum(P[x][x] for x in range(n))
    E_tau_plus = sum(
        pi[x] * P_star[x][y] * (E_T_given[y] if x == y else h[x][y])
        for x in range(n)
        for y in range(n)
    )
    assert E_tau_plus == Fraction(n - 1) + tr_P

    # Contrast 2: if robber stepped via P instead of P*
    E_tau_wrong_P = sum(pi[x] * P[x][y] * h[x][y] for x in range(n) for y in range(n))
    if check_non_rev_contrast:
        assert E_tau_wrong_P != Fraction(n - 1)

    print(
        f"  [{name}] N={n}: E[T]={E_T}, E[tau]={E_tau} (==N-1), "
        f"PMF match k=0..15 OK, tr(P)={tr_P}, E[tau_+]={E_tau_plus}, "
        f"E[tau if robber uses P]={E_tau_wrong_P}"
    )


def main():
    print("=== IB/2019/P1/q20H Exact Verification ===")

    # 1. N = 1 degenerate case
    verify_chain("N=1 trivial", [[Fraction(1)]])

    # 2. N = 2 reversible aperiodic
    verify_chain(
        "N=2 aperiodic",
        [
            [Fraction(1, 3), Fraction(2, 3)],
            [Fraction(1, 4), Fraction(3, 4)],
        ],
    )

    # 3. N = 3 non-reversible aperiodic
    verify_chain(
        "N=3 non-rev aperiodic",
        [
            [Fraction(1, 5), Fraction(3, 5), Fraction(1, 5)],
            [Fraction(1, 6), Fraction(1, 3), Fraction(1, 2)],
            [Fraction(3, 4), Fraction(1, 8), Fraction(1, 8)],
        ],
        check_non_rev_contrast=True,
    )

    # 4. N = 4 non-reversible sparse aperiodic
    verify_chain(
        "N=4 sparse non-rev",
        [
            [Fraction(1, 4), Fraction(3, 4), Fraction(0), Fraction(0)],
            [Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(0)],
            [Fraction(1, 5), Fraction(0), Fraction(0), Fraction(4, 5)],
            [Fraction(1, 2), Fraction(1, 6), Fraction(1, 6), Fraction(1, 6)],
        ],
        check_non_rev_contrast=True,
    )

    # 5. N = 5 non-reversible aperiodic
    verify_chain(
        "N=5 non-rev aperiodic",
        [
            [Fraction(1, 10), Fraction(2, 5), Fraction(1, 5), Fraction(1, 5), Fraction(1, 10)],
            [Fraction(1, 4), Fraction(1, 4), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)],
            [Fraction(0), Fraction(1, 3), Fraction(1, 6), Fraction(1, 2), Fraction(0)],
            [Fraction(1, 2), Fraction(0), Fraction(1, 4), Fraction(0), Fraction(1, 4)],
            [Fraction(1, 3), Fraction(1, 3), Fraction(0), Fraction(1, 6), Fraction(1, 6)],
        ],
        check_non_rev_contrast=True,
    )

    # 6. Bonus: N = 3 periodic non-reversible 3-cycle (aperiodicity not needed!)
    verify_chain(
        "N=3 periodic 3-cycle",
        [
            [Fraction(0), Fraction(1), Fraction(0)],
            [Fraction(0), Fraction(0), Fraction(1)],
            [Fraction(1), Fraction(0), Fraction(0)],
        ],
        check_non_rev_contrast=True,
    )

    print("\nALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
