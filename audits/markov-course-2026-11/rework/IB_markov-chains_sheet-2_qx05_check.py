#!/usr/bin/env python3
"""Exact rational verification for IB/markov-chains/sheet-2/qx05 (2012-13 Sheet 2, Q13)."""

from fractions import Fraction


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


def g(u):
    """Long-run average profit (in thousands of pounds) as a function of u = sqrt(x)."""
    return (Fraction(3, 2) * u - u * u) / (1 + 2 * u)


def g_prime(u):
    """Exact derivative g'(u)."""
    return (3 + 2 * u) * (1 - 2 * u) / (2 * (1 + 2 * u) ** 2)


def main():
    # 1. Edge case x = 0 (u = 0): reducible chain, state 0 absorbing, state 1 transient
    P0 = [
        [Fraction(1), Fraction(0)],
        [Fraction(1, 2), Fraction(1, 2)],
    ]
    for n in range(1, 21):
        Pn = mat_pow(P0, n)
        assert Pn == [
            [Fraction(1), Fraction(0)],
            [1 - Fraction(1, 2**n), Fraction(1, 2**n)],
        ]
    # Unique invariant distribution is (1, 0), giving long-run average profit 0
    pi0 = [Fraction(1), Fraction(0)]
    assert [sum(pi0[i] * P0[i][j] for i in range(2)) for j in range(2)] == pi0
    # Expected average profit over first N nights starting from state 1:
    # (1/N) * sum_{k=0}^{N-1} (3/4) * (1/2)^k = (3/(4N)) * (2 - 2^{1-N}) -> 0
    for N in (1, 5, 10, 20):
        avg_from_1 = sum(Fraction(3, 4) * mat_pow(P0, k)[1][1] for k in range(N)) / N
        closed_form = Fraction(3, 4 * N) * (2 - Fraction(1, 2 ** (N - 1)))
        assert avg_from_1 == closed_form
    assert g(Fraction(0)) == Fraction(0)

    # 2. Irreducible case 0 < x <= 1 (u = sqrt(x) in (0, 1])
    grid = [Fraction(k, 100) for k in range(1, 101)] + [Fraction(1, 3), Fraction(2, 3), Fraction(1, 4), Fraction(3, 4)]
    for u in grid:
        x = u * u
        P = [
            [1 - u, u],
            [Fraction(1, 2), Fraction(1, 2)],
        ]
        # Verify stochasticity and irreducibility (off-diagonal > 0) and aperiodicity (p_11 > 0)
        assert all(val >= 0 for row in P for val in row)
        assert all(sum(row) == 1 for row in P)
        assert P[0][1] > 0 and P[1][0] > 0 and P[1][1] > 0

        # Exact invariant distribution
        pi = [Fraction(1) / (1 + 2 * u), (2 * u) / (1 + 2 * u)]
        assert sum(pi) == 1
        assert [sum(pi[i] * P[i][j] for i in range(2)) for j in range(2)] == pi

        # Long-run profit from invariant distribution vs formula g(u)
        profit = Fraction(3, 4) * pi[1] - x * pi[0]
        assert profit == g(u)

        # Renewal-cycle check: E[K] = 2 concerts, E[M] = 1/u flower days per cycle
        E_K = Fraction(2)
        E_M = Fraction(1) / u
        renewal_profit = (Fraction(3, 4) * E_K - x * E_M) / (E_K + E_M)
        assert renewal_profit == g(u)

        # Exact algebraic global-maximum identity: 1/4 - g(u) = (2u - 1)^2 / (4(1 + 2u))
        diff = Fraction(1, 4) - g(u)
        assert diff == (2 * u - 1) ** 2 / (4 * (1 + 2 * u))
        if u == Fraction(1, 2):
            assert diff == 0
            assert g_prime(u) == 0
        else:
            assert diff > 0
            if u < Fraction(1, 2):
                assert g_prime(u) > 0
            else:
                assert g_prime(u) < 0

        # Verify derivative formula via exact symbolic/rational secant identity:
        # (g(u + h) - g(u)) / h at u with symbolic-style expansion check for tiny h
        h = Fraction(1, 10**6)
        secant = (g(u + h) - g(u)) / h
        assert abs(secant - g_prime(u)) < Fraction(1, 10**4)

        # Verify Bellman optimality inequality for dynamic policies with V(0)=0, V(1)=1
        V = [Fraction(0), Fraction(1)]
        bellman_1 = Fraction(3, 4) + P[1][0] * V[0] + P[1][1] * V[1] - V[1]
        bellman_0 = -x + P[0][0] * V[0] + P[0][1] * V[1] - V[0]
        assert bellman_1 == Fraction(1, 4)
        assert bellman_0 == Fraction(1, 4) - (u - Fraction(1, 2)) ** 2 <= Fraction(1, 4)

    # 3. Verify exact convergence of Cesaro averages at optimum x = 1/4 (u = 1/2)
    u_opt = Fraction(1, 2)
    x_opt = u_opt * u_opt
    P_opt = [
        [1 - u_opt, u_opt],
        [Fraction(1, 2), Fraction(1, 2)],
    ]
    # Note: at u = 1/2, P_opt has identical rows (1/2, 1/2), so P_opt^n = P_opt for all n >= 1!
    for n in range(1, 10):
        assert mat_pow(P_opt, n) == [[Fraction(1, 2), Fraction(1, 2)], [Fraction(1, 2), Fraction(1, 2)]]

    r_opt = [-x_opt, Fraction(3, 4)]
    for init in (0, 1):
        avg_N100 = sum(sum(mat_pow(P_opt, k)[init][j] * r_opt[j] for j in range(2)) for k in range(100)) / 100
        # Exact value: (r_opt[init] + 99 * (1/4)) / 100
        assert avg_N100 == (r_opt[init] + 99 * Fraction(1, 4)) / 100

    print("IB/markov-chains/sheet-2/qx05 verification PASSED")
    print(f"Optimal u = sqrt(x) = {u_opt}, optimal x = {x_opt} thousand pounds = £{int(x_opt * 1000)}")
    print(f"At optimum: pi = (1/2, 1/2), long-run profit = {g(u_opt)} thousand pounds/night = £{int(g(u_opt) * 1000)}/night")
    print(f"Endpoint values: g(0) = {g(Fraction(0))} (£0/night), g(1) = {g(Fraction(1))} (£{float(g(Fraction(1))*1000):.2f}/night)")
    print("Global-maximum identity 1/4 - g(u) = (2u-1)^2 / (4(1+2u)) and Bellman inequality verified on 104 rational points.")


if __name__ == "__main__":
    main()
