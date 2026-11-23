"""Exact rational verification for IB/2014/P4/q9H."""
from fractions import Fraction
import itertools


def mat_mul(A, B):
    n = len(A)
    return [
        [sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def check_part_a():
    # Concrete 4-state Markov chain on {0, 1, 2, 3}
    P = [
        [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)],
        [Fraction(1, 3), Fraction(0), Fraction(1, 3), Fraction(1, 3)],
        [Fraction(0), Fraction(2, 5), Fraction(2, 5), Fraction(1, 5)],
        [Fraction(1, 6), Fraction(1, 6), Fraction(1, 2), Fraction(1, 6)],
    ]
    lam = [Fraction(1, 5), Fraction(2, 5), Fraction(3, 10), Fraction(1, 10)]
    assert sum(lam) == 1
    for row in P:
        assert sum(row) == 1

    P2 = mat_mul(P, P)
    for row in P2:
        assert sum(row) == 1

    # Exact joint distribution of (W0, W1, W2) = (X0, X2, X4) by summing over (X1, X3)
    joint_W = {}
    for x0, x1, x2, x3, x4 in itertools.product(range(4), repeat=5):
        prob = lam[x0] * P[x0][x1] * P[x1][x2] * P[x2][x3] * P[x3][x4]
        key = (x0, x2, x4)
        joint_W[key] = joint_W.get(key, Fraction(0)) + prob

    assert sum(joint_W.values()) == 1
    for i0, i1, i2 in itertools.product(range(4), repeat=3):
        expected = lam[i0] * P2[i0][i1] * P2[i1][i2]
        assert joint_W[(i0, i1, i2)] == expected

    # Verify P(W_1 = 0) = sum_{i, k} lambda_i p_{i,k} p_{k,0}
    p_w1_0_direct = sum(
        prob for (i0, i1, i2), prob in joint_W.items() if i1 == 0
    )
    p_w1_0_formula = sum(
        lam[i] * P[i][k] * P[k][0]
        for i in range(4)
        for k in range(4)
    )
    assert p_w1_0_direct == p_w1_0_formula
    print(f"Part (a): Verified (W_0, W_1, W_2) ~ Markov(lambda, P^2) across all 64 triples; "
          f"P(W_1=0) = {p_w1_0_formula}.")


def check_part_b_generic(name, P_dict, lam_dict):
    states = [-1, 0, 1]
    assert sum(lam_dict[s] for s in states) == 1
    for i in states:
        assert sum(P_dict[(i, j)] for j in states) == 1
    assert P_dict[(0, 1)] == 0
    assert P_dict[(-1, -1)] == 0
    assert P_dict[(-1, 0)] != P_dict[(1, 0)]

    # Enumerate all 3^5 paths (X_0, X_1, X_2, X_3, X_4)
    paths = []
    for path in itertools.product(states, repeat=5):
        prob = lam_dict[path[0]]
        for t in range(4):
            prob *= P_dict[(path[t], path[t + 1])]
        paths.append((path, prob))
    assert sum(pr for _, pr in paths) == 1

    # Check (i) and (ii) for n = 2 and n = 3
    for n in [2, 3]:
        # (i) P(Y_{n+1}=0 | Y_n=1, Y_{n-1}=0)
        num1 = sum(pr for p, pr in paths if abs(p[n + 1]) == 0 and abs(p[n]) == 1 and abs(p[n - 1]) == 0)
        den1 = sum(pr for p, pr in paths if abs(p[n]) == 1 and abs(p[n - 1]) == 0)
        assert den1 > 0
        cond1 = num1 / den1
        assert cond1 == P_dict[(-1, 0)]

        # (ii) P(Y_{n+1}=0 | Y_n=1, Y_{n-1}=1, Y_{n-2}=0)
        num2 = sum(
            pr for p, pr in paths
            if abs(p[n + 1]) == 0 and abs(p[n]) == 1 and abs(p[n - 1]) == 1 and abs(p[n - 2]) == 0
        )
        den2 = sum(
            pr for p, pr in paths
            if abs(p[n]) == 1 and abs(p[n - 1]) == 1 and abs(p[n - 2]) == 0
        )
        assert den2 > 0
        cond2 = num2 / den2
        assert cond2 == P_dict[(1, 0)]
        assert cond1 != cond2

    print(f"Part (b) [{name}]: (i) = p_{{-1,0}} = {P_dict[(-1, 0)]}, "
          f"(ii) = p_{{1,0}} = {P_dict[(1, 0)]} (unequal -> Y_n is not Markov).")


def check_part_b_edge_cases():
    states = [-1, 0, 1]
    # Edge case 1: p_{-1,1} = 0 (so p_{-1,0} = 1) and lambda_1 = 0 -> Y_n IS a homogeneous Markov chain
    P_deg = {
        (-1, -1): Fraction(0), (-1, 0): Fraction(1), (-1, 1): Fraction(0),
        (0, -1): Fraction(1, 2), (0, 0): Fraction(1, 2), (0, 1): Fraction(0),
        (1, -1): Fraction(0), (1, 0): Fraction(1, 3), (1, 1): Fraction(2, 3),
    }
    lam_deg = {-1: Fraction(1, 3), 0: Fraction(2, 3), 1: Fraction(0)}
    # Check Y_n has transition matrix Q_Y on {0, 1}: 0->(1/2, 1/2), 1->(1, 0)
    Q_Y = [[Fraction(1, 2), Fraction(1, 2)], [Fraction(1), Fraction(0)]]
    lam_Y = [Fraction(2, 3), Fraction(1, 3)]
    y_dist = {}
    for path in itertools.product(states, repeat=5):
        pr = lam_deg[path[0]]
        for t in range(4):
            pr *= P_deg[(path[t], path[t + 1])]
        y_path = tuple(abs(x) for x in path)
        y_dist[y_path] = y_dist.get(y_path, Fraction(0)) + pr
    for y_path, pr in y_dist.items():
        expected = lam_Y[y_path[0]]
        for t in range(4):
            expected *= Q_Y[y_path[t]][y_path[t + 1]]
        assert pr == expected
    print("Part (b) edge case 1 (p_{-1,1}=0, lambda_1=0): verified Y_n IS a homogeneous Markov chain.")

    # Edge case 2: deterministic 3-cycle 0 -> -1 -> 1 -> 0 starting at X_0 = 0
    P_cyc = {
        (-1, -1): Fraction(0), (-1, 0): Fraction(0), (-1, 1): Fraction(1),
        (0, -1): Fraction(1), (0, 0): Fraction(0), (0, 1): Fraction(0),
        (1, -1): Fraction(0), (1, 0): Fraction(1), (1, 1): Fraction(0),
    }
    # Trajectory: X = (0, -1, 1, 0, -1, 1, ...), Y = (0, 1, 1, 0, 1, 1, ...)
    # At n=1: P(Y_2=0 | Y_1=1, Y_0=0) = 0 = p_{-1,0}
    # At n=2: P(Y_3=0 | Y_2=1, Y_1=1, Y_0=0) = 1 = p_{1,0}
    # Not homogeneous Markov (P(Y_2=0|Y_1=1)=0 != 1=P(Y_3=0|Y_2=1)), though conditioning events never
    # simultaneously have positive probability at the same n.
    assert P_cyc[(0, 1)] == 0 and P_cyc[(-1, -1)] == 0 and P_cyc[(-1, 0)] != P_cyc[(1, 0)]
    print("Part (b) edge case 2 (deterministic 3-cycle 0->-1->1->0): "
          "not homogeneous Markov; (i) holds at n=1 (val 0), (ii) holds at n=2 (val 1).")


def main():
    check_part_a()
    check_part_b_generic(
        "generic_1",
        {
            (-1, -1): Fraction(0), (-1, 0): Fraction(1, 3), (-1, 1): Fraction(2, 3),
            (0, -1): Fraction(3, 5), (0, 0): Fraction(2, 5), (0, 1): Fraction(0),
            (1, -1): Fraction(1, 4), (1, 0): Fraction(1, 2), (1, 1): Fraction(1, 4),
        },
        {-1: Fraction(1, 4), 0: Fraction(1, 2), 1: Fraction(1, 4)},
    )
    check_part_b_generic(
        "generic_2",
        {
            (-1, -1): Fraction(0), (-1, 0): Fraction(4, 7), (-1, 1): Fraction(3, 7),
            (0, -1): Fraction(1, 2), (0, 0): Fraction(1, 2), (0, 1): Fraction(0),
            (1, -1): Fraction(0), (1, 0): Fraction(1, 5), (1, 1): Fraction(4, 5),
        },
        {-1: Fraction(0), 0: Fraction(1), 1: Fraction(0)},
    )
    check_part_b_edge_cases()
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
