"""Exact rational verification for IB/2013/P1/q20H."""
from fractions import Fraction
import itertools

STATES = ["a", "b", "c"]
IDX = {s: i for i, s in enumerate(STATES)}

P = [
    [Fraction(0), Fraction(3, 5), Fraction(2, 5)],
    [Fraction(3, 4), Fraction(0), Fraction(1, 4)],
    [Fraction(2, 3), Fraction(1, 3), Fraction(0)],
]


def solve_linear(A, b):
    """Exact Gaussian elimination over Q for square matrix A and vector b."""
    n = len(b)
    M = [[Fraction(A[i][j]) for j in range(n)] + [Fraction(b[i])] for i in range(n)]
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
    return [M[i][n] for i in range(n)]


def main():
    # 1. Stochasticity
    for i, row in enumerate(P):
        assert sum(row) == 1, f"Row {i} does not sum to 1"
        assert all(x >= 0 for x in row)

    # 2. Equivalent 6-edge multigraph on {a, b, c}
    edges = {("a", "b"): 3, ("a", "c"): 2, ("b", "c"): 1}
    total_edges = sum(edges.values())
    assert total_edges == 6
    deg = {}
    for u in STATES:
        deg[u] = sum(
            edges.get((u, v), edges.get((v, u), 0))
            for v in STATES
            if v != u
        )
    assert deg == {"a": 5, "b": 4, "c": 3}
    assert sum(deg.values()) == 2 * total_edges == 12
    for u in STATES:
        for v in STATES:
            e_uv = 0 if u == v else edges.get((u, v), edges.get((v, u), 0))
            assert P[IDX[u]][IDX[v]] == Fraction(e_uv, deg[u])

    # 3. Stationary distribution and detailed balance
    pi = [Fraction(5, 12), Fraction(4, 12), Fraction(3, 12)]
    assert sum(pi) == 1
    for i in range(3):
        for j in range(3):
            assert pi[i] * P[i][j] == pi[j] * P[j][i], "Detailed balance failed"
    for j in range(3):
        assert sum(pi[i] * P[i][j] for i in range(3)) == pi[j], "Invariance failed"

    # 4. Mean hitting times h(i, j) = E_i[H_j] (0 when i=j) and first passage times k(i, j) = E_i[T_j]
    h = [[Fraction(0) for _ in range(3)] for _ in range(3)]
    for target in range(3):
        others = [i for i in range(3) if i != target]
        A = []
        b = []
        for i in others:
            row = [(Fraction(1) if i == j else Fraction(0)) - P[i][j] for j in others]
            A.append(row)
            b.append(Fraction(1))
        sol = solve_linear(A, b)
        for idx_o, i in enumerate(others):
            h[i][target] = sol[idx_o]

    k = [[Fraction(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            if i != j:
                k[i][j] = h[i][j]
            else:
                # First return time: 1 + sum_l p_{il} h(l, i)
                k[i][i] = Fraction(1) + sum(P[i][l] * h[l][i] for l in range(3))
                assert k[i][i] == Fraction(1) / pi[i]

    assert k[0][0] == Fraction(12, 5)
    assert k[1][1] == Fraction(3, 1)
    assert k[2][2] == Fraction(4, 1)
    assert k[0][1] == Fraction(21, 11)
    assert k[1][0] == Fraction(15, 11)
    assert k[2][0] == Fraction(16, 11)
    assert k[2][1] == Fraction(25, 11)

    print("Part (0): 6-edge multigraph edges:", edges, "degrees:", deg)
    print("Part (i): k(a,a) =", k[0][0], ", k(a,b) =", k[0][1])
    print("All k(i,j):")
    for i, si in enumerate(STATES):
        print(f"  from {si}: " + ", ".join(f"k({si},{sj})={k[i][j]}" for j, sj in enumerate(STATES)))

    # 5. Part (ii): Expected steps until first completing a step b -> c
    # Transient substochastic matrix Q on {a, b, c} removing the b -> c transition
    Q = [[P[i][j] for j in range(3)] for i in range(3)]
    Q[IDX["b"]][IDX["c"]] = Fraction(0)
    I_minus_Q = [
        [(Fraction(1) if i == j else Fraction(0)) - Q[i][j] for j in range(3)]
        for i in range(3)
    ]
    g = solve_linear(I_minus_Q, [Fraction(1), Fraction(1), Fraction(1)])
    assert g == [Fraction(128, 11), Fraction(107, 11), Fraction(12, 1)]
    # Cross-check with strong Markov / edge return time identities
    assert g[0] == k[0][1] + g[1]
    assert g[1] == Fraction(4) + Fraction(3) * k[0][1]
    assert g[2] == Fraction(1) / (pi[1] * P[1][2])
    print("Part (ii): g(a) =", g[0], ", g(b) =", g[1], ", g(c) =", g[2])

    # 6. Part (iii): Reversibility over all paths of length m = 1..5
    for m in range(1, 6):
        for path in itertools.product(range(3), repeat=m + 1):
            p_fwd = pi[path[0]]
            for t in range(m):
                p_fwd *= P[path[t]][path[t + 1]]
            rev = path[::-1]
            p_rev = pi[rev[0]]
            for t in range(m):
                p_rev *= P[rev[t]][rev[t + 1]]
            assert p_fwd == p_rev
    print("Part (iii): Exact path-reversal distribution equality verified for all 3^(m+1) paths, m=1..5.")

    # Exact distribution and expectation of tau(u, v) via augmented absorbing chain
    # Transient states: (vertex, stage) where stage=0 means u not yet seen, stage=1 means u seen but v not yet seen after u.
    def analyze_tau(u, v):
        t_states = []
        for s in range(3):
            if s != u:
                t_states.append((s, 0))
            if s != v:
                t_states.append((s, 1))
        t_idx = {st: idx for idx, st in enumerate(t_states)}
        n_t = len(t_states)
        Q_tau = [[Fraction(0) for _ in range(n_t)] for _ in range(n_t)]
        for st, i_st in t_idx.items():
            s, stage = st
            for nxt in range(3):
                prob = P[s][nxt]
                if prob == 0:
                    continue
                if stage == 0:
                    nxt_stage = 1 if nxt == u else 0
                    Q_tau[i_st][t_idx[(nxt, nxt_stage)]] += prob
                else:
                    if nxt == v:
                        # absorbed!
                        pass
                    else:
                        Q_tau[i_st][t_idx[(nxt, 1)]] += prob
        # Initial distribution over transient states at time 0
        init_t = [Fraction(0) for _ in range(n_t)]
        for s in range(3):
            st0 = (s, 1) if s == u else (s, 0)
            init_t[t_idx[st0]] += pi[s]
        # Exact P(tau > m) for m = 0..25
        vec = list(init_t)
        tail = [sum(vec)]
        for _ in range(25):
            vec = [sum(vec[r] * Q_tau[r][c] for r in range(n_t)) for c in range(n_t)]
            tail.append(sum(vec))
        # Exact expected hitting times from each transient state
        A_tau = [
            [(Fraction(1) if r == c else Fraction(0)) - Q_tau[r][c] for c in range(n_t)]
            for r in range(n_t)
        ]
        exp_t = solve_linear(A_tau, [Fraction(1)] * n_t)
        exp_pi = sum(init_t[r] * exp_t[r] for r in range(n_t))
        return tail, exp_pi

    tail_ab, exp_ab = analyze_tau(0, 1)
    tail_ba, exp_ba = analyze_tau(1, 0)
    assert tail_ab == tail_ba, "Distributions of tau(a,b) and tau(b,a) differ!"
    assert exp_ab == exp_ba == Fraction(30, 11)
    print("Part (iii): P_pi(tau(a,b) > m) == P_pi(tau(b,a) > m) verified exactly for m=0..25.")
    print("Part (iii): E_pi[tau(a,b)] =", exp_ab, ", E_pi[tau(b,a)] =", exp_ba)

    # Verify the identity for all pairs (u, v) under BOTH conventions k(i,i)=1/pi_i and h(i,i)=0
    for u in range(3):
        for v in range(3):
            if u == v:
                continue
            lhs = k[v][u] - k[u][v]
            rhs_k = sum(pi[i] * (k[i][u] - k[i][v]) for i in range(3))
            rhs_h = sum(pi[i] * (h[i][u] - h[i][v]) for i in range(3))
            assert lhs == rhs_k == rhs_h
            tail_uv, exp_uv = analyze_tau(u, v)
            tail_vu, exp_vu = analyze_tau(v, u)
            assert tail_uv == tail_vu and exp_uv == exp_vu
    print("Part (iii): k(b,a) - k(a,b) =", k[1][0] - k[0][1],
          "== sum pi_i [k(i,a) - k(i,b)] =", sum(pi[i] * (k[i][0] - k[i][1]) for i in range(3)))
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
