"""Exact rational check for IB/2004/P2/q11H (reversibility <=> detailed balance)."""
from fractions import Fraction


def mat_vec_left(pi, P):
    n = len(pi)
    return [sum(pi[i] * P[i][j] for i in range(n)) for j in range(n)]


def check_joint_distribution(pi, P, N):
    """Return True iff P(X_0=i0,...,X_N=iN) == P(X_N=i0,...,X_0=iN) for all paths."""
    n = len(pi)
    from itertools import product

    for path in product(range(n), repeat=N + 1):
        p_fwd = pi[path[0]]
        for t in range(N):
            p_fwd *= P[path[t]][path[t + 1]]
        rev = path[::-1]
        p_rev = pi[rev[0]]
        for t in range(N):
            p_rev *= P[rev[t]][rev[t + 1]]
        if p_fwd != p_rev:
            return False, path, p_fwd, p_rev
    return True, None, None, None


def detailed_balance(pi, P):
    n = len(pi)
    return all(pi[i] * P[i][j] == pi[j] * P[j][i] for i in range(n) for j in range(n))


# Example 1: Irreducible, positive-recurrent, reversible 3-state chain
P_rev = [
    [Fraction(1, 2), Fraction(1, 2), Fraction(0)],
    [Fraction(1, 4), Fraction(1, 4), Fraction(1, 2)],
    [Fraction(0), Fraction(1, 3), Fraction(2, 3)],
]
# pi_0 * (1/2) = pi_1 * (1/4) => pi_1 = 2 pi_0
# pi_1 * (1/2) = pi_2 * (1/3) => pi_2 = (3/2) pi_1 = 3 pi_0
pi_rev = [Fraction(1, 6), Fraction(2, 6), Fraction(3, 6)]
assert sum(pi_rev) == 1
assert mat_vec_left(pi_rev, P_rev) == pi_rev
assert detailed_balance(pi_rev, P_rev)
for N in range(1, 6):
    ok, *_ = check_joint_distribution(pi_rev, P_rev, N)
    assert ok, f"Failed reversibility at N={N}"
print("Example 1 (reversible 3-state chain): invariant, detailed balance=True, joint laws equal for N=1..5.")

# Example 2: Irreducible, positive-recurrent, NON-reversible 3-state cycle
P_nonrev = [
    [Fraction(1, 6), Fraction(2, 3), Fraction(1, 6)],
    [Fraction(1, 6), Fraction(1, 6), Fraction(2, 3)],
    [Fraction(2, 3), Fraction(1, 6), Fraction(1, 6)],
]
pi_nonrev = [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)]
assert sum(pi_nonrev) == 1
assert mat_vec_left(pi_nonrev, P_nonrev) == pi_nonrev
assert not detailed_balance(pi_nonrev, P_nonrev)
ok, witness, p_fwd, p_rev = check_joint_distribution(pi_nonrev, P_nonrev, 1)
assert not ok
print(f"Example 2 (non-reversible 3-state chain): detailed balance=False, N=1 witness {witness}: fwd={p_fwd} != rev={p_rev}.")

# Check time-reversal transition matrix P_hat_{ij} = pi_j P_{ji} / pi_i on Example 2
P_hat = [
    [pi_nonrev[j] * P_nonrev[j][i] / pi_nonrev[i] for j in range(3)]
    for i in range(3)
]
assert all(sum(row) == 1 for row in P_hat)
assert mat_vec_left(pi_nonrev, P_hat) == pi_nonrev
# Verify reversed joint law equals Markov(pi, P_hat) for N=1..4
from itertools import product
for N in range(1, 5):
    for path in product(range(3), repeat=N + 1):
        rev = path[::-1]
        p_rev_orig = pi_nonrev[rev[0]]
        for t in range(N):
            p_rev_orig *= P_nonrev[rev[t]][rev[t + 1]]
        p_hat_fwd = pi_nonrev[path[0]]
        for t in range(N):
            p_hat_fwd *= P_hat[path[t]][path[t + 1]]
        assert p_rev_orig == p_hat_fwd
print("Time-reversed law matches Markov(pi, P_hat) identically for N=1..4; P_hat != P.")
print("ALL CHECKS PASSED.")
