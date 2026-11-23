"""Exact rational verification for IB/markov-chains/sheet-1/q03."""

from fractions import Fraction
import itertools


def mat_mul(A, B):
    n = len(A)
    m = len(B[0])
    k = len(B)
    return [
        [sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)]
        for i in range(n)
    ]


def mat_pow(P, e):
    n = len(P)
    res = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    base = [row[:] for row in P]
    while e > 0:
        if e & 1:
            res = mat_mul(res, base)
        base = mat_mul(base, base)
        e >>= 1
    return res


def vec_mat(v, P):
    m = len(P[0])
    return [sum(v[i] * P[i][j] for i in range(len(v))) for j in range(m)]


def verify_subsampled_homogeneous():
    """Verify Y_r = X_{n_r} is an inhomogeneous Markov chain with Q(r) = P^{n_{r+1}-n_r}."""
    P = [
        [Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)],
        [Fraction(1, 4), Fraction(0), Fraction(3, 4)],
        [Fraction(2, 5), Fraction(2, 5), Fraction(1, 5)],
    ]
    lam = [Fraction(1, 6), Fraction(1, 2), Fraction(1, 3)]
    for row in P:
        assert sum(row) == 1
    assert sum(lam) == 1

    # Irregular increasing sequence of positive integers (and test starting at 0 too)
    for seq in ([1, 2, 5, 7], [0, 1, 3, 6], [1, 3, 3, 5]):
        max_t = seq[-1]
        # Build exact joint distribution of (Y_0, ..., Y_R) directly by summing over all X paths
        # via forward propagation on full X paths of length max_t (3^8 = 6561 max, instant)
        y_joint = {}
        for path in itertools.product(range(3), repeat=max_t + 1):
            prob = lam[path[0]]
            for t in range(max_t):
                prob *= P[path[t]][path[t + 1]]
            if prob == 0:
                continue
            y_tuple = tuple(path[idx] for idx in seq)
            y_joint[y_tuple] = y_joint.get(y_tuple, Fraction(0)) + prob

        assert sum(y_joint.values()) == 1

        # Theoretical initial distribution and transition matrices
        mu = vec_mat(lam, mat_pow(P, seq[0]))
        Qs = [mat_pow(P, seq[r + 1] - seq[r]) for r in range(len(seq) - 1)]

        # Check joint distribution formula for ALL tuples
        for y_tuple in itertools.product(range(3), repeat=len(seq)):
            expected = mu[y_tuple[0]]
            for r in range(len(seq) - 1):
                expected *= Qs[r][y_tuple[r]][y_tuple[r + 1]]
            actual = y_joint.get(y_tuple, Fraction(0))
            assert actual == expected, f"Mismatch for {seq}, {y_tuple}: {actual} != {expected}"

        # Check conditional Markov property for every prefix with positive probability
        for r in range(len(seq) - 1):
            prefix_probs = {}
            next_probs = {}
            for y_tuple, pr in y_joint.items():
                pref = y_tuple[: r + 1]
                prefix_probs[pref] = prefix_probs.get(pref, Fraction(0)) + pr
                nxt = y_tuple[: r + 2]
                next_probs[nxt] = next_probs.get(nxt, Fraction(0)) + pr
            for pref, p_pref in prefix_probs.items():
                assert p_pref > 0
                for j in range(3):
                    cond = next_probs.get(pref + (j,), Fraction(0)) / p_pref
                    assert cond == Qs[r][pref[-1]][j]

        if seq == [1, 2, 5, 7]:
            assert Qs[0] != Qs[1] and Qs[1] != Qs[2], "Expected genuinely inhomogeneous Q(r)"

    print("PASS: Subsampled homogeneous chain Y_r = X_{n_r} verified for sequences [1,2,5,7], [0,1,3,6], [1,3,3,5].")


def verify_subsampled_inhomogeneous():
    """Verify Y_r = X_{n_r} when X itself is inhomogeneous."""
    Ps = [
        [[Fraction(1, 2), Fraction(1, 2)], [Fraction(1, 3), Fraction(2, 3)]],
        [[Fraction(1, 4), Fraction(3, 4)], [Fraction(2, 5), Fraction(3, 5)]],
        [[Fraction(3, 5), Fraction(2, 5)], [Fraction(1, 6), Fraction(5, 6)]],
        [[Fraction(1, 3), Fraction(2, 3)], [Fraction(3, 4), Fraction(1, 4)]],
        [[Fraction(4, 5), Fraction(1, 5)], [Fraction(1, 2), Fraction(1, 2)]],
    ]
    lam = [Fraction(2, 7), Fraction(5, 7)]
    seq = [1, 2, 5]
    y_joint = {}
    for path in itertools.product(range(2), repeat=6):
        prob = lam[path[0]]
        for t in range(5):
            prob *= Ps[t][path[t]][path[t + 1]]
        y_tuple = tuple(path[idx] for idx in seq)
        y_joint[y_tuple] = y_joint.get(y_tuple, Fraction(0)) + prob

    mu = vec_mat(lam, Ps[0])
    Q0 = Ps[1]
    Q1 = mat_mul(mat_mul(Ps[2], Ps[3]), Ps[4])
    for y_tuple in itertools.product(range(2), repeat=3):
        expected = mu[y_tuple[0]] * Q0[y_tuple[0]][y_tuple[1]] * Q1[y_tuple[1]][y_tuple[2]]
        assert y_joint.get(y_tuple, Fraction(0)) == expected

    print("PASS: Subsampled inhomogeneous chain verified with Q(r) = P(n_r)...P(n_{r+1}-1).")


def verify_srw_even_times():
    """Verify Y_r = X_{2r} for simple random walk on Z with general p in [0,1]."""
    test_ps = [Fraction(0), Fraction(1, 5), Fraction(1, 3), Fraction(1, 2), Fraction(3, 4), Fraction(1)]
    R = 3  # steps 0, 2, 4, 6 -> 2^6 = 64 paths per starting state
    for p in test_ps:
        q = 1 - p
        # Expected 2-step transition law from any state i
        expected_step = {2: p * p, 0: 2 * p * q, -2: q * q}
        assert sum(expected_step.values()) == 1

        # Exact enumeration from starting state i0 = 0 (and i0 = 3 for odd parity)
        for i0 in (0, 3):
            y_joint = {}
            for steps in itertools.product((-1, 1), repeat=2 * R):
                prob = Fraction(1)
                pos = [i0]
                for s in steps:
                    prob *= p if s == 1 else q
                    pos.append(pos[-1] + s)
                if prob == 0:
                    continue
                y_tuple = tuple(pos[2 * r] for r in range(R + 1))
                y_joint[y_tuple] = y_joint.get(y_tuple, Fraction(0)) + prob

            assert sum(y_joint.values()) == 1
            for y_tuple, prob in y_joint.items():
                theo = Fraction(1)
                for r in range(R):
                    diff = y_tuple[r + 1] - y_tuple[r]
                    theo *= expected_step.get(diff, Fraction(0))
                assert prob == theo, f"SRW mismatch for p={p}, {y_tuple}: {prob} != {theo}"

            # Check conditional probabilities directly
            for r in range(R):
                pref_probs = {}
                next_probs = {}
                for y_tuple, prob in y_joint.items():
                    pref = y_tuple[: r + 1]
                    pref_probs[pref] = pref_probs.get(pref, Fraction(0)) + prob
                    nxt = y_tuple[: r + 2]
                    next_probs[nxt] = next_probs.get(nxt, Fraction(0)) + prob
                for pref, p_pref in pref_probs.items():
                    cur = pref[-1]
                    for diff in range(-4, 5):
                        actual = next_probs.get(pref + (cur + diff,), Fraction(0)) / p_pref
                        assert actual == expected_step.get(diff, Fraction(0))

    print("PASS: Simple random walk Y_r = X_{2r} transition probabilities verified for p in {0, 1/5, 1/3, 1/2, 3/4, 1}.")
    print("  General p: q_{i,i+2} = p^2, q_{i,i} = 2p(1-p), q_{i,i-2} = (1-p)^2.")
    print("  Symmetric p = 1/2: q_{i,i+2} = 1/4, q_{i,i} = 1/2, q_{i,i-2} = 1/4.")


if __name__ == "__main__":
    verify_subsampled_homogeneous()
    verify_subsampled_inhomogeneous()
    verify_srw_even_times()
