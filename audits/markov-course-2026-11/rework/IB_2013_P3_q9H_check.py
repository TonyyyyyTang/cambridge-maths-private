"""Exact rational verification for IB/2013/P3/q9H."""
from fractions import Fraction
import math


def transition_prob(s1, s2, alpha, beta, gamma):
    """Exact transition probability p_{(i,j),(i',j')} for the 2-urn model."""
    i, j = s1
    i2, j2 = s2
    n = i + j
    if n == 0:
        if (i2, j2) == (0, 0):
            return beta + gamma
        if (i2, j2) == (1, 0) or (i2, j2) == (0, 1):
            return alpha / 2
        return Fraction(0)
    # n >= 1
    if (i2, j2) == (i + 1, j) or (i2, j2) == (i, j + 1):
        return alpha / 2
    if (i2, j2) == (i - 1, j) and i >= 1:
        return gamma * Fraction(i, n)
    if (i2, j2) == (i, j - 1) and j >= 1:
        return gamma * Fraction(j, n)
    if (i2, j2) == (i - 1, j + 1) and i >= 1:
        return beta * Fraction(i, n)
    if (i2, j2) == (i + 1, j - 1) and j >= 1:
        return beta * Fraction(j, n)
    return Fraction(0)


def neighbors(s):
    """All states reachable from s in one step (plus s itself)."""
    i, j = s
    cand = [(i, j), (i + 1, j), (i, j + 1)]
    if i >= 1:
        cand.append((i - 1, j))
        cand.append((i - 1, j + 1))
    if j >= 1:
        cand.append((i, j - 1))
        cand.append((i + 1, j - 1))
    return cand


def lam(i, j, alpha, gamma):
    n = i + j
    binom = Fraction(math.factorial(n), math.factorial(i) * math.factorial(j))
    ratio = alpha / (2 * gamma)
    # Note: 0**0 == 1 in Python
    return binom * (ratio ** n)


def verify_params(alpha, beta, gamma, max_n=10):
    assert alpha >= 0 and beta >= 0 and gamma > 0
    assert alpha + beta + gamma == 1
    assert alpha < gamma

    states = [(i, n - i) for n in range(max_n + 1) for i in range(n + 1)]

    # 1. Row sums = 1
    for s in states:
        row_sum = sum(transition_prob(s, s2, alpha, beta, gamma) for s2 in neighbors(s))
        assert row_sum == 1, f"Row sum failed at {s}: {row_sum}"

    # 2. Detailed balance for all pairs with i+j <= max_n
    for s1 in states:
        l1 = lam(s1[0], s1[1], alpha, gamma)
        for s2 in neighbors(s1):
            l2 = lam(s2[0], s2[1], alpha, gamma)
            lhs = l1 * transition_prob(s1, s2, alpha, beta, gamma)
            rhs = l2 * transition_prob(s2, s1, alpha, beta, gamma)
            assert lhs == rhs, f"Detailed balance failed for {s1} <-> {s2}: {lhs} != {rhs}"

    # 3. Direct global balance (lambda P)_{k,l} = lambda_{k,l} for k+l <= max_n - 1
    for n in range(max_n):
        for k in range(n + 1):
            l = n - k
            s_target = (k, l)
            # Predecessors can only be among neighbors(s_target) because transitions are symmetric in support
            inflow = sum(
                lam(s_prev[0], s_prev[1], alpha, gamma)
                * transition_prob(s_prev, s_target, alpha, beta, gamma)
                for s_prev in neighbors(s_target)
            )
            assert inflow == lam(k, l, alpha, gamma), f"Global balance failed at {s_target}"

    # 4. Normalization and n-ball mass formula
    Z = Fraction(1) / (Fraction(1) - alpha / gamma)
    for n in range(max_n + 5):
        mass_n = sum(lam(i, n - i, alpha, gamma) for i in range(n + 1))
        assert mass_n == (alpha / gamma) ** n
        prob_n = mass_n / Z
        expected_prob_n = (Fraction(1) - alpha / gamma) * ((alpha / gamma) ** n)
        assert prob_n == expected_prob_n

    # 5. Exact invariance on reflecting finite truncation i+j <= M (M=6)
    M = 6
    trunc_states = [(i, n - i) for n in range(M + 1) for i in range(n + 1)]
    trunc_set = set(trunc_states)
    Z_M = sum(lam(i, j, alpha, gamma) for i, j in trunc_states)
    pi_M = {s: lam(s[0], s[1], alpha, gamma) / Z_M for s in trunc_states}
    assert sum(pi_M.values()) == 1
    for s_t in trunc_states:
        inflow = Fraction(0)
        for s_p in trunc_states:
            if s_p == s_t:
                # self-loop plus redirected outward transitions
                out_outside = sum(
                    transition_prob(s_p, nxt, alpha, beta, gamma)
                    for nxt in neighbors(s_p)
                    if nxt not in trunc_set
                )
                p_val = transition_prob(s_p, s_t, alpha, beta, gamma) + out_outside
            else:
                p_val = transition_prob(s_p, s_t, alpha, beta, gamma)
            inflow += pi_M[s_p] * p_val
        assert inflow == pi_M[s_t]


def main():
    test_cases = [
        ("generic_1", Fraction(1, 6), Fraction(1, 2), Fraction(1, 3)),
        ("generic_2", Fraction(1, 5), Fraction(1, 5), Fraction(3, 5)),
        ("beta_zero", Fraction(1, 4), Fraction(0), Fraction(3, 4)),
        ("alpha_zero", Fraction(0), Fraction(1, 3), Fraction(2, 3)),
    ]
    for name, a, b, g in test_cases:
        verify_params(a, b, g, max_n=10)
        print(f"Verified {name} (alpha={a}, beta={b}, gamma={g}): "
              f"row sums, detailed balance, global balance, Z={1/(1-a/g)}, "
              f"P(n balls)=(1-{a/g})*({a/g})^n, truncated exact invariance OK.")
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
