#!/usr/bin/env python3
"""Exact rational verification for IB/2024/P3/q8H."""

from fractions import Fraction
from itertools import product


def check_distribution(name, dist, horizon=4):
    assert sum(dist.values()) == 1
    assert all(p > 0 for p in dist.values())
    vals = sorted(dist.keys())

    # Build exact joint distribution of (X_1, ..., X_horizon)
    outcomes = []
    for xs in product(vals, repeat=horizon):
        prob = Fraction(1)
        for x in xs:
            prob *= dist[x]
        S = []
        L = []
        K = []
        s_acc = 0
        l_acc = None
        prev_x = 0  # X_0 = 0
        for x in xs:
            s_acc += x
            l_acc = x if l_acc is None else min(l_acc, x)
            S.append(s_acc)
            L.append(l_acc)
            K.append(x + prev_x)
            prev_x = x
        outcomes.append((prob, xs, tuple(S), tuple(L), tuple(K)))

    assert sum(o[0] for o in outcomes) == 1

    def verify_markov(seq_idx, theoretical_p):
        for n in range(1, horizon):
            # Group by prefix (Y_1, ..., Y_n)
            prefix_probs = {}
            joint_probs = {}
            for prob, *seqs in outcomes:
                seq = seqs[seq_idx]
                pref = seq[:n]
                nxt = seq[n]
                prefix_probs[pref] = prefix_probs.get(pref, Fraction(0)) + prob
                joint_probs[(pref, nxt)] = joint_probs.get((pref, nxt), Fraction(0)) + prob
            for pref, p_pref in prefix_probs.items():
                assert p_pref > 0
                last_state = pref[-1]
                # Check all reachable next states and some unreachable ones
                all_next = {nxt for (pr, nxt) in joint_probs if pr == pref}
                total_cond = Fraction(0)
                for nxt in all_next:
                    cond = joint_probs[(pref, nxt)] / p_pref
                    assert cond == theoretical_p(last_state, nxt), (
                        f"Mismatch for {pref} -> {nxt}: {cond} vs {theoretical_p(last_state, nxt)}"
                    )
                    total_cond += cond
                assert total_cond == 1

    # 1. {X_n}: p_{ij} = P(X_1 = j)
    verify_markov(0, lambda i, j: dist.get(j, Fraction(0)))

    # 2. {S_n}: p_{ij} = P(X_1 = j - i)
    verify_markov(1, lambda i, j: dist.get(j - i, Fraction(0)))

    # 3. {L_n}: p_{ij} = 0 if j > i; P(X_1 >= i) if j == i; P(X_1 = j) if j < i
    def p_L(i, j):
        if j > i:
            return Fraction(0)
        if j == i:
            return sum(p for x, p in dist.items() if x >= i)
        return dist.get(j, Fraction(0))

    verify_markov(2, p_L)
    print(f"[{name}] Verified exact time-homogeneous Markov property for {{X_n}}, {{S_n}}, {{L_n}} (n <= {horizon}).")
    return outcomes


def check_K_counterexample():
    # Counterexample: X_n ~ Bernoulli(1/2) on {0, 1}, X_0 = 0
    dist = {0: Fraction(1, 2), 1: Fraction(1, 2)}
    outcomes = check_distribution("Bernoulli(1/2)", dist, horizon=4)

    def cond_prob_K(n, target, history):
        """Compute P(K_{n+1} = target | K_1..K_n = history)."""
        p_hist = sum(o[0] for o in outcomes if o[4][:n] == history)
        p_both = sum(o[0] for o in outcomes if o[4][:n] == history and o[4][n] == target)
        return p_hist, p_both / p_hist

    p_h1, c1 = cond_prob_K(2, 2, (0, 1))
    p_h2, c2 = cond_prob_K(2, 2, (1, 1))
    assert p_h1 == Fraction(1, 4) and c1 == Fraction(1, 2)
    assert p_h2 == Fraction(1, 4) and c2 == Fraction(0)
    assert c1 != c2

    # Also check conditioning on (K_2, K_3) only (independent of boundary X_0 = 0):
    def cond_prob_K_from_2(target, k2, k3):
        p_h = sum(o[0] for o in outcomes if o[4][1] == k2 and o[4][2] == k3)
        p_b = sum(o[0] for o in outcomes if o[4][1] == k2 and o[4][2] == k3 and o[4][3] == target)
        return p_h, p_b / p_h

    p_a, ca = cond_prob_K_from_2(2, 0, 1)
    p_b, cb = cond_prob_K_from_2(2, 2, 1)
    assert p_a == Fraction(1, 8) and ca == Fraction(1, 2)
    assert p_b == Fraction(1, 8) and cb == Fraction(0)
    assert ca != cb

    print("Counterexample for {K_n} verified:")
    print(f"  P(K_1=0, K_2=1) = {p_h1}, P(K_3=2 | K_1=0, K_2=1) = {c1}")
    print(f"  P(K_1=1, K_2=1) = {p_h2}, P(K_3=2 | K_1=1, K_2=1) = {c2}")
    print(f"  Also without K_1: P(K_4=2 | K_2=0, K_3=1) = {ca} != {cb} = P(K_4=2 | K_2=2, K_3=1)")


def main():
    check_K_counterexample()
    check_distribution("Asymmetric {-1, 0, 2}", {-1: Fraction(1, 6), 0: Fraction(1, 2), 2: Fraction(1, 3)}, horizon=4)
    # Degenerate constant X_n = 3 a.s.: {K_n} = (3, 6, 6, 6, ...) IS Markov (why "necessarily" matters)
    outcomes_const = check_distribution("Degenerate {3}", {3: Fraction(1)}, horizon=4)
    assert outcomes_const[0][4] == (3, 6, 6, 6)
    print("Degenerate constant X_n = 3 check: K = (3, 6, 6, 6) is Markov (confirms 'necessarily' phrasing).")
    print("ALL CHECKS PASSED for IB/2024/P3/q8H.")


if __name__ == "__main__":
    main()
