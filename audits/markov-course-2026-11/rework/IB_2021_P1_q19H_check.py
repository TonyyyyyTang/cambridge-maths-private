"""Exact rational check for IB/2021/P1/q19H (Unopoly)."""
from fractions import Fraction


def solve_linear_fraction(A, b):
    """Exact Gaussian elimination over Q."""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        div = M[col][col]
        M[col] = [x / div for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                factor = M[r][col]
                M[r] = [M[r][c] - factor * M[col][c] for c in range(n + 1)]
    return [M[i][n] for i in range(n)]


# 1. Verify cubic equation for alpha = P_{(2,m)}(hit (2, m-1))
alpha = Fraction(2, 3)
assert alpha == Fraction(1, 2) + Fraction(1, 8) * alpha**2 + Fraction(3, 8) * alpha**3
# Verify factorization 3 a^3 + a^2 - 8 a + 4 = (a - 1)(3a - 2)(a + 2)
for a in [Fraction(-3), Fraction(-2), Fraction(0), Fraction(2, 3), Fraction(1), Fraction(5, 2)]:
    lhs = 3 * a**3 + a**2 - 8 * a + 4
    rhs = (a - 1) * (3 * a - 2) * (a + 2)
    assert lhs == rhs
print("Cubic check: alpha = 2/3 satisfies alpha = 1/2 + (1/8)alpha^2 + (3/8)alpha^3, roots are {-2, 2/3, 1}.")


# 2. Verify full harmonic equations on S = {1,2,4} x {0,1,...,K}
def g(s, k):
    if s == 2:
        return Fraction(2, 3) ** k
    if s == 4:
        return Fraction(1, 2) * (Fraction(2, 3) ** k)
    if s == 1:
        return Fraction(0) if k == 0 else Fraction(5, 4) * (Fraction(2, 3) ** k)
    raise ValueError(s)


for k in range(1, 15):
    assert g(2, k) == Fraction(1, 2) * g(2, k - 1) + Fraction(1, 2) * g(4, k)
    assert g(1, k) == Fraction(1, 2) * g(2, k) + Fraction(1, 2) * g(2, k - 1)
    assert g(4, k) == Fraction(1, 2) * g(1, k + 2) + Fraction(1, 2) * g(2, k + 2)
    assert g(1, k) == Fraction(5, 6) * (Fraction(2, 3) ** (k - 1))
assert g(2, 0) == 1
assert g(4, 0) == Fraction(1, 2) * g(1, 2) + Fraction(1, 2) * g(2, 2)
print("Harmonic equations verified identically for k=1..14.")


# 3. Verify by exact finite-truncation hitting probabilities (absorbing at wealth >= M as upper/lower bounds)
# Lower bound: absorb at 0 if wealth >= M (so hitting (2,0) before wealth >= M is <= true hitting prob)
# Upper bound: reflect / set h=1 at wealth >= M? Even better: compute exact hitting probability before wealth >= M
# and check monotone convergence as M -> infinity to (2/3)^m and (5/6)(2/3)^{m-1}!
def truncated_hit_prob(M_cap):
    # States: (s, k) for s in (1,2,4), 0 <= k < M_cap; states with k >= M_cap have value 0 (absorbing win)
    states = [(s, k) for k in range(M_cap) for s in (1, 2, 4)]
    idx = {st: i for i, st in enumerate(states)}
    n = len(states)
    A = [[Fraction(0)] * n for _ in range(n)]
    b = [Fraction(0)] * n
    for i, (s, k) in enumerate(states):
        A[i][i] = Fraction(1)
        if (s, k) == (2, 0):
            b[i] = Fraction(1)
        elif (s, k) == (1, 0):
            b[i] = Fraction(0)
        elif s == 1:
            # -> (2, k) w.p. 1/2, -> (2, k-1) w.p. 1/2
            A[i][idx[(2, k)]] -= Fraction(1, 2)
            A[i][idx[(2, k - 1)]] -= Fraction(1, 2)
        elif s == 2:
            # -> (2, k-1) w.p. 1/2, -> (4, k) w.p. 1/2
            A[i][idx[(2, k - 1)]] -= Fraction(1, 2)
            A[i][idx[(4, k)]] -= Fraction(1, 2)
        elif s == 4:
            # -> (1, k+2) w.p. 1/2, -> (2, k+2) w.p. 1/2
            if k + 2 < M_cap:
                A[i][idx[(1, k + 2)]] -= Fraction(1, 2)
                A[i][idx[(2, k + 2)]] -= Fraction(1, 2)
    sol = solve_linear_fraction(A, b)
    return {st: sol[idx[st]] for st in states}


sol25 = truncated_hit_prob(25)
sol35 = truncated_hit_prob(35)
for m in range(1, 6):
    exact_2 = g(2, m)
    exact_1 = g(1, m)
    err2_25 = float(exact_2 - sol25[(2, m)])
    err2_35 = float(exact_2 - sol35[(2, m)])
    err1_35 = float(exact_1 - sol35[(1, m)])
    assert sol25[(2, m)] < sol35[(2, m)] < exact_2
    assert sol25[(1, m)] < sol35[(1, m)] < exact_1
    assert err2_35 < 1e-5 and err1_35 < 1e-5
    print(
        f"m={m}: h(2,{m}) exact={exact_2} ({float(exact_2):.8f}), "
        f"trunc35={float(sol35[(2,m)]):.8f}; "
        f"h(1,{m}) exact={exact_1} ({float(exact_1):.8f}), "
        f"trunc35={float(sol35[(1,m)]):.8f}"
    )

print("ALL CHECKS PASSED.")
