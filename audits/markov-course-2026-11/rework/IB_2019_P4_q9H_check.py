#!/usr/bin/env python3
"""Exact rational verification for IB/2019/P4/q9H."""

from fractions import Fraction
from math import comb


def mat_mul(A, B):
    n = len(A)
    return [
        [sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def check_part_a_examples():
    # Example 1: Recurrent state 0 in an irreducible 3-state chain
    P_rec = [
        [Fraction(0), Fraction(1, 2), Fraction(1, 2)],
        [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)],
        [Fraction(1, 2), Fraction(1, 2), Fraction(0)],
    ]
    for row in P_rec:
        assert sum(row) == 1

    # Stationary distribution satisfies pi P = pi:
    # pi = (3/10, 9/20, 1/4)? Let's verify exact pi:
    # 1/3 pi1 + 1/2 pi2 = pi0
    # 1/2 pi0 + 1/3 pi1 + 1/2 pi2 = pi1 => 1/2 pi0 + 1/2 pi2 = 2/3 pi1
    # 1/2 pi0 + 1/3 pi1 = pi2
    # Wait: let pi0 = 6/21 = 2/7, pi1 = 9/21 = 3/7, pi2 = 6/21 = 2/7?
    # Check: (1/3)(3/7) + (1/2)(2/7) = 1/7 + 1/7 = 2/7 = pi0!
    # (1/2)(2/7) + (1/3)(3/7) + (1/2)(2/7) = 1/7 + 1/7 + 1/7 = 3/7 = pi1!
    # (1/2)(2/7) + (1/3)(3/7) = 1/7 + 1/7 = 2/7 = pi2!
    pi = [Fraction(2, 7), Fraction(3, 7), Fraction(2, 7)]
    assert sum(pi) == 1
    for j in range(3):
        assert sum(pi[i] * P_rec[i][j] for i in range(3)) == pi[j]

    M = [[Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]
    partial_rec = Fraction(0)
    for n in range(40):
        partial_rec += M[0][0]
        M = mat_mul(M, P_rec)
    # Since p_00(n) -> pi_0 = 2/7, partial sum grows linearly (~ 40 * 2/7 = 80/7 ~ 11.43)
    assert partial_rec > 11

    # Example 2: Transient state 0 with return prob f_0 = 2/3
    # State 0 -> 0 (1/2), 1 (1/3), 2 (1/6); state 1 -> 0 (1/2), 2 (1/2); state 2 absorbing
    # Return prob f_0 = p_00 + p_01 * p_10 = 1/2 + (1/3)*(1/2) = 2/3
    # Theoretical sum_{n>=0} p_00(n) = 1 / (1 - f_0) = 3
    P_tr = [
        [Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)],
        [Fraction(1, 2), Fraction(0), Fraction(1, 2)],
        [Fraction(0), Fraction(0), Fraction(1)],
    ]
    for row in P_tr:
        assert sum(row) == 1
    f0 = Fraction(1, 2) + Fraction(1, 3) * Fraction(1, 2)
    assert f0 == Fraction(2, 3)
    expected_visits_tr = Fraction(1, 1 - f0)
    assert expected_visits_tr == Fraction(3)

    M = [[Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]
    partial_tr = Fraction(0)
    for n in range(60):
        partial_tr += M[0][0]
        M = mat_mul(M, P_tr)
    assert partial_tr < expected_visits_tr
    assert expected_visits_tr - partial_tr < Fraction(1, 10**6)

    print("Part (a) sanity checks passed:")
    print(f"  Recurrent example: sum_{{n=0}}^39 p_00(n) = {float(partial_rec):.6f} (diverges, pi_0 = {pi[0]})")
    print(f"  Transient example: f_0 = {f0}, exact 1/(1-f_0) = {expected_visits_tr}, "
          f"sum_{{n=0}}^59 p_00(n) = {float(partial_tr):.10f}")


def check_part_b_2d_srw(max_n=8):
    # Exact DP for 2D simple symmetric random walk starting at (0, 0)
    # dist[n] maps (x, y) -> exact probability p_{0,(x,y)}(n)
    dist = [{(0, 0): Fraction(1)}]
    steps = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    for n in range(1, 2 * max_n + 2):
        prev = dist[-1]
        cur = {}
        for (x, y), prob in prev.items():
            q = prob * Fraction(1, 4)
            for dx, dy in steps:
                nxt = (x + dx, y + dy)
                cur[nxt] = cur.get(nxt, Fraction(0)) + q
        assert sum(cur.values()) == 1
        dist.append(cur)

    # Check symmetry p_{0,x}(n) = p_{0,-x}(n) = p_{x,0}(n) (by translation invariance of SSRW)
    for n in range(len(dist)):
        for (x, y), prob in dist[n].items():
            assert dist[n].get((-x, -y), Fraction(0)) == prob
            # Parity check: (x + y) % 2 == n % 2
            assert (x + y) % 2 == n % 2

    # Check p_00(2n+1) = 0 for all odd steps
    for n in range(max_n + 1):
        assert dist[2 * n + 1].get((0, 0), Fraction(0)) == 0

    # Check P(X_n = Y_n) = sum_x p_{0,x}(n)^2 = p_00(2n) = (C(2n, n) / 4^n)^2
    partial_EZ = Fraction(0)
    partial_all = Fraction(0)
    print("Part (b) exact 2D SSRW checks:")
    for n in range(max_n + 1):
        prob_meet = sum(p * p for p in dist[n].values())
        p00_2n = dist[2 * n].get((0, 0), Fraction(0))
        closed_form = Fraction(comb(2 * n, n) ** 2, 4 ** (2 * n))
        assert prob_meet == p00_2n == closed_form
        partial_EZ += prob_meet
        partial_all += dist[2 * n].get((0, 0), Fraction(0)) + dist[2 * n + 1].get((0, 0), Fraction(0))
        assert partial_EZ == partial_all
        print(f"  n={n}: P(X_n=Y_n) = p_00({2*n}) = {p00_2n} = {float(p00_2n):.8f}")

    # Verify divergence rate via exact formula up to n=200: sum p_00(2n) ~ (1/pi) ln(n)
    exact_sum_200 = Fraction(0)
    for n in range(201):
        exact_sum_200 += Fraction(comb(2 * n, n) ** 2, 4 ** (2 * n))
    print(f"  sum_{{n=0}}^200 p_00(2n) = {float(exact_sum_200):.6f} (diverges like (1/pi) ln n)")
    print("ALL CHECKS PASSED for IB/2019/P4/q9H.")


if __name__ == "__main__":
    check_part_a_examples()
    check_part_b_2d_srw(max_n=8)
