"""Exact rational verification for IB/markov-chains/sheet-1/q04."""

from fractions import Fraction
import itertools


def expected_X_trans(p, k, j):
    q = 1 - p
    if k == 0:
        return Fraction(1) if j == 1 else Fraction(0)
    denom = p ** k + q ** k
    if j == k + 1:
        return (p ** (k + 1) + q ** (k + 1)) / denom
    if j == k - 1:
        return (p * q * (p ** (k - 1) + q ** (k - 1))) / denom
    return Fraction(0)


def expected_Y_trans(p, i, j):
    q = 1 - p
    if i == 0:
        if j == 0:
            return p
        if j == 1:
            return q
        return Fraction(0)
    if j == i - 1:
        return p
    if j == i + 1:
        return q
    return Fraction(0)


def verify_X_and_Y():
    N = 8  # 2^8 = 256 paths per parameter p -> exact enumeration
    test_ps = [Fraction(0), Fraction(1, 5), Fraction(1, 3), Fraction(1, 2), Fraction(2, 3), Fraction(1)]

    total_X_histories = 0
    total_Y_histories = 0

    for p in test_ps:
        q = 1 - p
        # Check row sums of theoretical transition matrices for k = 0..10
        for k in range(11):
            row_sum_X = sum(expected_X_trans(p, k, j) for j in range(12))
            row_sum_Y = sum(expected_Y_trans(p, k, j) for j in range(12))
            assert row_sum_X == 1, f"Row sum X != 1 at p={p}, k={k}"
            assert row_sum_Y == 1, f"Row sum Y != 1 at p={p}, k={k}"

        # Enumerate all paths of S_0, ..., S_N with S_0 = 0
        paths_data = []
        for steps in itertools.product((-1, 1), repeat=N):
            prob = Fraction(1)
            S = [0]
            M = [0]
            for xi in steps:
                prob *= p if xi == 1 else q
                s_new = S[-1] + xi
                S.append(s_new)
                M.append(max(M[-1], s_new))
            X = tuple(abs(s) for s in S)
            Y = tuple(M[t] - S[t] for t in range(N + 1))

            # Pointwise check of Lindley recursion Y_{n+1} = max(Y_n - xi_{n+1}, 0)
            for t in range(N):
                assert Y[t + 1] == max(Y[t] - steps[t], 0)

            if prob > 0:
                paths_data.append((prob, S, X, Y))

        assert sum(item[0] for item in paths_data) == 1

        # Verify Markov property and exact transition probabilities for X_n = |S_n|
        for n in range(N):
            hist_prob = {}
            hist_sign_plus = {}
            hist_next = {}
            for prob, S, X, _ in paths_data:
                hist = X[: n + 1]
                hist_prob[hist] = hist_prob.get(hist, Fraction(0)) + prob
                if S[n] > 0:
                    hist_sign_plus[hist] = hist_sign_plus.get(hist, Fraction(0)) + prob
                nxt = (hist, X[n + 1])
                hist_next[nxt] = hist_next.get(nxt, Fraction(0)) + prob

            for hist, pr in hist_prob.items():
                total_X_histories += 1
                k = hist[-1]
                if k >= 1:
                    cond_plus = hist_sign_plus.get(hist, Fraction(0)) / pr
                    assert cond_plus == (p ** k) / (p ** k + q ** k), (
                        f"Sign dist mismatch at p={p}, hist={hist}: {cond_plus}"
                    )
                for j in range(N + 2):
                    cond = hist_next.get((hist, j), Fraction(0)) / pr
                    exp = expected_X_trans(p, k, j)
                    assert cond == exp, f"X trans mismatch at p={p}, hist={hist}, j={j}: {cond} != {exp}"

        # Verify Markov property and exact transition probabilities for Y_n = M_n - S_n
        for n in range(N):
            hist_prob = {}
            hist_next = {}
            for prob, _, _, Y in paths_data:
                hist = Y[: n + 1]
                hist_prob[hist] = hist_prob.get(hist, Fraction(0)) + prob
                nxt = (hist, Y[n + 1])
                hist_next[nxt] = hist_next.get(nxt, Fraction(0)) + prob

            for hist, pr in hist_prob.items():
                total_Y_histories += 1
                i = hist[-1]
                for j in range(N + 2):
                    cond = hist_next.get((hist, j), Fraction(0)) / pr
                    exp = expected_Y_trans(p, i, j)
                    assert cond == exp, f"Y trans mismatch at p={p}, hist={hist}, j={j}: {cond} != {exp}"

    print(f"PASS: X_n = |S_n| Markov property & transition probabilities verified across {total_X_histories} histories (N={N}, 6 values of p).")
    print(f"PASS: Y_n = M_n - S_n Markov property, Lindley recursion & transition probabilities verified across {total_Y_histories} histories.")


def verify_S0_nonzero_counterexample():
    """Verify why S_0 = 0 is essential for X_n = |S_n| when p != 1/2."""
    p = Fraction(1, 3)
    q = Fraction(2, 3)
    # Start at S_0 = 1 (so X_0 = 1)
    # P(X_1 = 2 | X_0 = 1) = P(S_1 = 2 | S_0 = 1) = p = 1/3
    # P(X_3 = 2 | X_0 = 1, X_1 = 0, X_2 = 1):
    #   S_1 = 0 (prob q); then S_2 = +1 (prob p) or -1 (prob q)
    #   so given X_0=1, X_1=0, X_2=1, S_2 = +1 w.p. p and -1 w.p. q -> X_3 = 2 w.p. p^2 + q^2 = 5/9 != 1/3
    p1 = Fraction(0)
    num_3 = Fraction(0)
    den_3 = Fraction(0)
    for steps in itertools.product((-1, 1), repeat=3):
        prob = Fraction(1)
        S = [1]
        for xi in steps:
            prob *= p if xi == 1 else q
            S.append(S[-1] + xi)
        X = tuple(abs(s) for s in S)
        if X[1] == 2:
            p1 += prob
        if X[:3] == (1, 0, 1):
            den_3 += prob
            if X[3] == 2:
                num_3 += prob
    cond1 = p1
    cond3 = num_3 / den_3
    assert cond1 == Fraction(1, 3)
    assert cond3 == Fraction(5, 9)
    assert cond1 != cond3
    print(f"PASS: Counterexample with S_0=1, p=1/3 confirmed: P(X_1=2|X_0=1)={cond1} != P(X_3=2|X_0=1,X_1=0,X_2=1)={cond3}.")


if __name__ == "__main__":
    verify_X_and_Y()
    verify_S0_nonzero_counterexample()
