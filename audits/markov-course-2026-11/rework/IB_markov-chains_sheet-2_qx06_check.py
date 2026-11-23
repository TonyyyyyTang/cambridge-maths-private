#!/usr/bin/env python3
"""Exact rational verification for IB/markov-chains/sheet-2/qx06 (2012-13 Sheet 2, Q14)."""

from collections import deque
from fractions import Fraction
from math import gcd


def build_Q(P, alpha):
    """Construct exact Metropolis-Hastings transition matrix Q from P and positive weights alpha."""
    n = len(P)
    assert len(alpha) == n and all(a > 0 for a in alpha)
    Q = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        off_sum = Fraction(0)
        for j in range(n):
            if i != j:
                q_ij = min(P[i][j], alpha[j] * P[j][i] / alpha[i])
                Q[i][j] = q_ij
                off_sum += q_ij
        Q[i][i] = Fraction(1) - off_sum
    return Q


def check_stochastic(M):
    n = len(M)
    for i in range(n):
        assert len(M[i]) == n
        assert all(val >= 0 for val in M[i]), f"Negative entry in row {i}: {M[i]}"
        assert sum(M[i]) == 1, f"Row {i} sums to {sum(M[i])}, not 1"


def check_detailed_balance_and_invariance(Q, alpha):
    n = len(Q)
    Z = sum(alpha)
    pi = [a / Z for a in alpha]
    assert sum(pi) == 1
    for i in range(n):
        for j in range(n):
            assert pi[i] * Q[i][j] == pi[j] * Q[j][i], f"Detailed balance failed at ({i},{j})"
    pi_Q = [sum(pi[i] * Q[i][j] for i in range(n)) for j in range(n)]
    assert pi_Q == pi, "Invariance pi Q = pi failed"
    return pi


def check_proposal_accept_reject(P, alpha, Q):
    """Verify two-step proposal + accept/reject reproduces Q exactly."""
    n = len(P)
    for i in range(n):
        stay = P[i][i]
        for j in range(n):
            if i != j:
                if P[i][j] == 0:
                    move = Fraction(0)
                else:
                    a_ij = min(Fraction(1), (alpha[j] * P[j][i]) / (alpha[i] * P[i][j]))
                    move = P[i][j] * a_ij
                    stay += P[i][j] * (1 - a_ij)
                assert move == Q[i][j]
        assert stay == Q[i][i]


def reachable_from(M, start):
    n = len(M)
    seen = {start}
    q = deque([start])
    while q:
        u = q.popleft()
        for v in range(n):
            if M[u][v] > 0 and v not in seen:
                seen.add(v)
                q.append(v)
    return seen


def is_irreducible(M):
    n = len(M)
    return all(len(reachable_from(M, i)) == n for i in range(n))


def period_of_irreducible(M):
    """Exact period of an irreducible chain via gcd of (depth[u] + 1 - depth[v]) over all edges u->v."""
    assert is_irreducible(M)
    n = len(M)
    depth = {0: 0}
    q = deque([0])
    while q:
        u = q.popleft()
        for v in range(n):
            if M[u][v] > 0 and v not in depth:
                depth[v] = depth[u] + 1
                q.append(v)
    d = 0
    for u in range(n):
        for v in range(n):
            if M[u][v] > 0:
                diff = abs(depth[u] + 1 - depth[v])
                d = gcd(d, diff)
    return d


def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def mat_pow(A, k):
    n = len(A)
    res = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    base = [row[:] for row in A]
    while k > 0:
        if k & 1:
            res = mat_mul(res, base)
        base = mat_mul(base, base)
        k >>= 1
    return res


def max_tv_distance(M, pi):
    """max_i (1/2) sum_j |M_{ij} - pi_j|."""
    n = len(M)
    return max(sum(abs(M[i][j] - pi[j]) for j in range(n)) / 2 for i in range(n))


def main():
    # Example 1: Asymmetric proposal P with symmetric support, p_ii = 0, non-uniform alpha
    P1 = [
        [Fraction(0), Fraction(2, 3), Fraction(1, 3), Fraction(0)],
        [Fraction(1, 4), Fraction(0), Fraction(1, 2), Fraction(1, 4)],
        [Fraction(1, 5), Fraction(2, 5), Fraction(0), Fraction(2, 5)],
        [Fraction(0), Fraction(3, 5), Fraction(2, 5), Fraction(0)],
    ]
    alpha1 = [Fraction(1), Fraction(3), Fraction(2), Fraction(5)]
    check_stochastic(P1)
    assert is_irreducible(P1)
    Q1 = build_Q(P1, alpha1)
    check_stochastic(Q1)
    pi1 = check_detailed_balance_and_invariance(Q1, alpha1)
    check_proposal_accept_reject(P1, alpha1, Q1)
    assert is_irreducible(Q1)
    assert period_of_irreducible(Q1) == 1
    tv1_10 = max_tv_distance(mat_pow(Q1, 10), pi1)
    tv1_30 = max_tv_distance(mat_pow(Q1, 30), pi1)
    assert tv1_30 < tv1_10 < Fraction(1, 10)

    # Example 2 (Counterexample to irreducibility of Q): directed 3-cycle P
    P2 = [
        [Fraction(0), Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(1)],
        [Fraction(1), Fraction(0), Fraction(0)],
    ]
    alpha2 = [Fraction(1), Fraction(2), Fraction(3)]
    check_stochastic(P2)
    assert is_irreducible(P2)
    Q2 = build_Q(P2, alpha2)
    check_stochastic(Q2)
    pi2 = check_detailed_balance_and_invariance(Q2, alpha2)
    check_proposal_accept_reject(P2, alpha2, Q2)
    assert Q2 == [[Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]
    assert not is_irreducible(Q2)
    assert max_tv_distance(mat_pow(Q2, 50), pi2) == Fraction(5, 6)  # no convergence

    # Example 3 (Counterexample 2 to irreducibility of Q): two symmetric pairs joined by directed bridge
    P3 = [
        [Fraction(0), Fraction(1), Fraction(0), Fraction(0)],
        [Fraction(1, 2), Fraction(0), Fraction(1, 2), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0), Fraction(1)],
        [Fraction(1, 2), Fraction(0), Fraction(1, 2), Fraction(0)],
    ]
    alpha3 = [Fraction(1), Fraction(2), Fraction(3), Fraction(4)]
    check_stochastic(P3)
    assert is_irreducible(P3)
    Q3 = build_Q(P3, alpha3)
    check_stochastic(Q3)
    check_detailed_balance_and_invariance(Q3, alpha3)
    check_proposal_accept_reject(P3, alpha3, Q3)
    assert not is_irreducible(Q3)
    assert reachable_from(Q3, 0) == {0, 1}
    assert reachable_from(Q3, 2) == {2, 3}

    # Example 4 (Counterexample to aperiodicity of Q): symmetric 2-cycle with uniform alpha
    P4 = [
        [Fraction(0), Fraction(1)],
        [Fraction(1), Fraction(0)],
    ]
    alpha4 = [Fraction(1), Fraction(1)]
    check_stochastic(P4)
    assert is_irreducible(P4)
    Q4 = build_Q(P4, alpha4)
    check_stochastic(Q4)
    pi4 = check_detailed_balance_and_invariance(Q4, alpha4)
    check_proposal_accept_reject(P4, alpha4, Q4)
    assert Q4 == P4
    assert is_irreducible(Q4)
    assert period_of_irreducible(Q4) == 2
    assert max_tv_distance(mat_pow(Q4, 50), pi4) == Fraction(1, 2)  # oscillates forever

    # Example 5 (Contrast: periodic symmetric P becomes APERIODIC Q when alpha is non-uniform)
    alpha5 = [Fraction(1), Fraction(2)]
    Q5 = build_Q(P4, alpha5)
    check_stochastic(Q5)
    pi5 = check_detailed_balance_and_invariance(Q5, alpha5)
    check_proposal_accept_reject(P4, alpha5, Q5)
    assert Q5 == [[Fraction(0), Fraction(1)], [Fraction(1, 2), Fraction(1, 2)]]
    assert is_irreducible(Q5)
    assert period_of_irreducible(Q5) == 1
    tv5_20 = max_tv_distance(mat_pow(Q5, 20), pi5)
    assert tv5_20 == Fraction(1, 3 * 2**19)

    # Example 6 (Lazy fix on Example 4): P_lazy = (I + P4)/2 makes Q aperiodic
    P6 = [[(Fraction(1 if i == j else 0) + P4[i][j]) / 2 for j in range(2)] for i in range(2)]
    Q6 = build_Q(P6, alpha4)
    check_stochastic(Q6)
    pi6 = check_detailed_balance_and_invariance(Q6, alpha4)
    assert is_irreducible(Q6) and period_of_irreducible(Q6) == 1
    assert max_tv_distance(mat_pow(Q6, 1), pi6) == Fraction(0)

    # Example 7 (Truncated geometric target on {0, ..., 6} with reflected SSRW proposal)
    m = 7
    P7 = [[Fraction(0)] * m for _ in range(m)]
    P7[0][1] = Fraction(1)
    P7[m - 1][m - 2] = Fraction(1)
    for i in range(1, m - 1):
        P7[i][i - 1] = Fraction(1, 2)
        P7[i][i + 1] = Fraction(1, 2)
    alpha7 = [Fraction(2, 3) ** k for k in range(m)]
    check_stochastic(P7)
    assert is_irreducible(P7)
    assert period_of_irreducible(P7) == 2
    Q7 = build_Q(P7, alpha7)
    check_stochastic(Q7)
    pi7 = check_detailed_balance_and_invariance(Q7, alpha7)
    check_proposal_accept_reject(P7, alpha7, Q7)
    assert is_irreducible(Q7)
    assert period_of_irreducible(Q7) == 1  # self-loops introduced by rejection!
    tv7_50 = max_tv_distance(mat_pow(Q7, 50), pi7)
    assert tv7_50 < Fraction(1, 100)

    print("IB/markov-chains/sheet-2/qx06 verification PASSED")
    print(f"Example 1 (4-state asymmetric P, sym support): irreducible=True, period=1, TV(10)={float(tv1_10):.6e}, TV(30)={float(tv1_30):.6e}")
    print(f"Example 2 (directed 3-cycle P): P irreducible=True, Q=I_3 irreducible=False, TV(50)={max_tv_distance(mat_pow(Q2, 50), pi2)}")
    print(f"Example 3 (2 sym pairs + directed bridge): P irreducible=True, Q irreducible=False (classes {{0,1}}, {{2,3}})")
    print(f"Example 4 (2-cycle, alpha=(1,1)): Q=P irreducible=True, period=2, TV(50)={max_tv_distance(mat_pow(Q4, 50), pi4)}")
    print(f"Example 5 (2-cycle, alpha=(1,2)): Q irreducible=True, period=1, TV(20)={tv5_20} ({float(tv5_20):.6e})")
    print(f"Example 6 (lazy 2-cycle, alpha=(1,1)): Q irreducible=True, period=1, TV(1)={max_tv_distance(mat_pow(Q6, 1), pi6)}")
    print(f"Example 7 (truncated Geom(2/3), m=7, periodic P): Q irreducible=True, period=1, TV(50)={float(tv7_50):.6e}")


if __name__ == "__main__":
    main()
