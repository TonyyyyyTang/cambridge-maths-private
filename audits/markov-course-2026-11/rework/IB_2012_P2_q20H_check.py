"""Exact rational verification for IB/2012/P2/q20H (symmetric random walk on a graph / 10-vertex tree)."""
from collections import deque
from fractions import Fraction


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


def build_transition_matrix(n, edges):
    adj = {i: [] for i in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    P = [[Fraction(0)] * n for _ in range(n)]
    for u in range(1, n + 1):
        d = len(adj[u])
        assert d >= 1
        for v in adj[u]:
            P[u - 1][v - 1] = Fraction(1, d)
    return adj, P


def compute_all_q_and_m(n, P):
    """Compute exact q_{ij} = P_i(T_j < T_i) (i != j) and m_{ij} = E_i(T_j) (all i,j)."""
    F = Fraction
    q = {}
    m = {}
    # Compute m_{ij} for each target j
    for j in range(1, n + 1):
        others = [k for k in range(1, n + 1) if k != j]
        idx_map = {k: idx for idx, k in enumerate(others)}
        sz = len(others)
        A = [[F(0)] * sz for _ in range(sz)]
        b = [F(1)] * sz
        for k in others:
            r = idx_map[k]
            A[r][r] = F(1)
            for ell in others:
                A[r][idx_map[ell]] -= P[k - 1][ell - 1]
        sol = solve_linear(A, b)
        for k in others:
            m[(k, j)] = sol[idx_map[k]]
        # Mean return time m_{jj} = 1 + sum_{ell != j} P_{j,ell} m_{ell,j}
        m[(j, j)] = F(1) + sum(P[j - 1][ell - 1] * m[(ell, j)] for ell in others)

    # Compute q_{ij} = P_i(T_j < T_i) for all i != j
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if i == j:
                continue
            # Solve h_k = P_k(hit j before i) for k not in {i, j}, with h_j=1, h_i=0
            others = [k for k in range(1, n + 1) if k not in (i, j)]
            idx_map = {k: idx for idx, k in enumerate(others)}
            sz = len(others)
            h = {i: F(0), j: F(1)}
            if sz > 0:
                A = [[F(0)] * sz for _ in range(sz)]
                b = [F(0)] * sz
                for k in others:
                    r = idx_map[k]
                    A[r][r] = F(1)
                    b[r] = P[k - 1][j - 1]  # h_j = 1
                    for ell in others:
                        A[r][idx_map[ell]] -= P[k - 1][ell - 1]
                sol = solve_linear(A, b)
                for k in others:
                    h[k] = sol[idx_map[k]]
            q[(i, j)] = sum(P[i - 1][ell - 1] * h[ell] for ell in range(1, n + 1))
    return q, m


def bfs_distances(adj, start):
    dist = {start: 0}
    dq = deque([start])
    while dq:
        u = dq.popleft()
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                dq.append(v)
    return dist


def component_edges_without(edges, u_rem, v_rem, start):
    """Count edges in the connected component of `start` after removing edge (u_rem, v_rem)."""
    rem_set = {u_rem, v_rem}
    adj = {}
    kept_edges = []
    for a, b in edges:
        if {a, b} == rem_set:
            continue
        kept_edges.append((a, b))
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    visited = {start}
    dq = deque([start])
    while dq:
        x = dq.popleft()
        for y in adj.get(x, []):
            if y not in visited:
                visited.add(y)
                dq.append(y)
    return sum(1 for a, b in kept_edges if a in visited and b in visited), visited


def verify_general_graph():
    """Verify parts (i) and (ii) on a connected graph WITH cycles."""
    n = 6
    edges = [(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 6), (6, 4), (2, 5)]
    adj, P = build_transition_matrix(n, edges)
    E = len(edges)
    pi = [Fraction(len(adj[i]), 2 * E) for i in range(1, n + 1)]
    assert sum(pi) == 1
    # Check detailed balance & invariance
    for i in range(n):
        for j in range(n):
            assert pi[i] * P[i][j] == pi[j] * P[j][i]
        assert sum(pi[k] * P[k][i] for k in range(n)) == pi[i]
    q, m = compute_all_q_and_m(n, P)
    for i in range(1, n + 1):
        assert m[(i, i)] == Fraction(1) / pi[i - 1]
        for j in range(1, n + 1):
            if i != j:
                assert pi[i - 1] * q[(i, j)] == pi[j - 1] * q[(j, i)]
                # Commute time identity holds on ANY connected graph!
                assert m[(i, j)] + m[(j, i)] == Fraction(1) / (pi[i - 1] * q[(i, j)])
    print(f"[PASS] General cyclic graph (n={n}, |E|={E}): detailed balance, pi_i q_ij = pi_j q_ji, commute times verified.")


def verify_printed_tree():
    """Verify parts (i)-(v) on the exact 10-vertex tree printed in IB/2012/P2/q20H."""
    n = 10
    edges = [
        (1, 3),
        (2, 3),
        (3, 5),
        (4, 5),
        (5, 6),
        (5, 8),
        (7, 8),
        (8, 9),
        (8, 10),
    ]
    assert len(edges) == 9
    adj, P = build_transition_matrix(n, edges)
    E = len(edges)
    degrees = {i: len(adj[i]) for i in range(1, n + 1)}
    assert degrees == {1: 1, 2: 1, 3: 3, 4: 1, 5: 4, 6: 1, 7: 1, 8: 4, 9: 1, 10: 1}
    assert sum(degrees.values()) == 2 * E == 18

    pi = [Fraction(degrees[i], 2 * E) for i in range(1, n + 1)]
    assert sum(pi) == 1
    for i in range(n):
        for j in range(n):
            assert pi[i] * P[i][j] == pi[j] * P[j][i]
        assert sum(pi[k] * P[k][i] for k in range(n)) == pi[i]

    q, m = compute_all_q_and_m(n, P)

    # Part (ii): pi_i q_ij = pi_j q_ji for all i != j
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if i != j:
                assert pi[i - 1] * q[(i, j)] == pi[j - 1] * q[(j, i)]

    # Part (iii): m_33 = 6
    assert m[(3, 3)] == Fraction(6)
    for i in range(1, n + 1):
        assert m[(i, i)] == Fraction(1) / pi[i - 1]

    # Part (iv): m_ij = 1 + 2 e_ij(i) for every directed edge (i, j)
    edge_results = []
    for u, v in edges:
        for i, j in [(u, v), (v, u)]:
            e_i, comp_i = component_edges_without(edges, u, v, i)
            e_j, comp_j = component_edges_without(edges, u, v, j)
            assert comp_i.isdisjoint(comp_j) and len(comp_i | comp_j) == n
            assert e_i + e_j == E - 1
            assert m[(i, j)] == 1 + 2 * e_i
            edge_results.append((i, j, e_i, m[(i, j)]))

    # Part (v): m_ij + m_ji = 18 * d(i, j) in {18, 36, 54, 72} for all i != j
    commute_vals = set()
    for i in range(1, n + 1):
        dists = bfs_distances(adj, i)
        for j in range(1, n + 1):
            if i == j:
                continue
            c_ij = m[(i, j)] + m[(j, i)]
            assert c_ij == 18 * dists[j]
            assert c_ij == Fraction(1) / (pi[i - 1] * q[(i, j)])
            commute_vals.add(c_ij)
    assert commute_vals == {Fraction(18), Fraction(36), Fraction(54), Fraction(72)}

    print(f"[PASS] Printed 10-vertex tree (|E|={E}, deg(3)={degrees[3]}):")
    print(f"       pi_3 = {pi[2]}, m_33 = {m[(3, 3)]}")
    print(f"       All 18 directed edges satisfy m_ij = 1 + 2*e_ij(i), e.g. m_35={m[(3,5)]} (e=2), m_53={m[(5,3)]} (e=6), m_58={m[(5,8)]} (e=5), m_85={m[(8,5)]} (e=3).")
    print(f"       Commute times {{m_ij + m_ji : i != j}} = {sorted(int(x) for x in commute_vals)}.")


def main():
    verify_general_graph()
    verify_printed_tree()
    print("All IB/2012/P2/q20H checks passed.")


if __name__ == "__main__":
    main()
