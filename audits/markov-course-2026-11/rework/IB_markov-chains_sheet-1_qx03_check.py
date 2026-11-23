"""Exact rational verification for IB/markov-chains/sheet-1/qx03."""

from fractions import Fraction


def verify_part1_finite_and_infinite():
    """Verify P(A | B_n) = p for all n (with P(B_n) > 0) => P(A) = p."""
    # Case 1: Finite partition with a zero-probability piece included
    p = Fraction(3, 7)
    pB = [Fraction(1, 6), Fraction(0), Fraction(1, 2), Fraction(1, 3)]
    assert sum(pB) == 1
    pA_and_B = [p * b for b in pB]
    assert sum(pA_and_B) == p

    # Case 2: Countably infinite partition P(B_n) = (1-r) r^{n-1} with r = 2/5,
    # truncated at N with exact remainder tail term
    r = Fraction(2, 5)
    N = 30
    partial_B = sum((1 - r) * (r ** (n - 1)) for n in range(1, N + 1))
    tail_B = r ** N
    assert partial_B + tail_B == 1
    partial_AB = sum(p * (1 - r) * (r ** (n - 1)) for n in range(1, N + 1))
    tail_AB = p * tail_B
    assert partial_AB + tail_AB == p
    print("PASS: Part 1 verified on finite (with zero-measure event) and geometric countably infinite partitions.")


def check_equivalence(joint):
    """Given joint[(x, y)] = P(X=x, Y=y), independently test (a) and (b) and assert (a) <=> (b)."""
    xs = sorted({x for x, _ in joint})
    ys = sorted({y for _, y in joint})
    assert sum(joint.values()) == 1

    marg_X = {x: sum(joint.get((x, y), Fraction(0)) for y in ys) for x in xs}
    marg_Y = {y: sum(joint.get((x, y), Fraction(0)) for x in xs) for y in ys}
    assert sum(marg_X.values()) == 1
    assert sum(marg_Y.values()) == 1

    # Condition (a): X and Y are independent
    cond_a = all(
        joint.get((x, y), Fraction(0)) == marg_X[x] * marg_Y[y]
        for x in xs
        for y in ys
    )

    # Condition (b): conditional distribution of X given Y=y is independent of y (over P(Y=y) > 0)
    pos_ys = [y for y in ys if marg_Y[y] > 0]
    assert len(pos_ys) > 0
    cond_dists = {
        y: tuple(joint.get((x, y), Fraction(0)) / marg_Y[y] for x in xs)
        for y in pos_ys
    }
    first_dist = cond_dists[pos_ys[0]]
    cond_b = all(dist == first_dist for dist in cond_dists.values())

    assert cond_a == cond_b
    if cond_b:
        # Also verify the common conditional distribution equals the marginal distribution of X
        assert first_dist == tuple(marg_X[x] for x in xs)
    return cond_a


def verify_part2():
    # Example 1: Independent discrete RVs, including a value y=2 with P(Y=2) = 0
    pX = {0: Fraction(1, 5), 1: Fraction(3, 10), 2: Fraction(1, 2)}
    pY = {0: Fraction(1, 4), 1: Fraction(3, 4), 2: Fraction(0)}
    joint_indep = {(x, y): pX[x] * pY[y] for x in pX for y in pY}
    assert check_equivalence(joint_indep) is True

    # Example 2: Degenerate Y (takes single value with prob 1) -> always independent
    joint_degen = {(0, 0): Fraction(2, 7), (1, 0): Fraction(5, 7), (0, 1): Fraction(0), (1, 1): Fraction(0)}
    assert check_equivalence(joint_degen) is True

    # Example 3: Dependent discrete RVs
    joint_dep = {
        (0, 0): Fraction(1, 4),
        (0, 1): Fraction(1, 8),
        (1, 0): Fraction(1, 8),
        (1, 1): Fraction(1, 2),
    }
    assert check_equivalence(joint_dep) is False

    # Example 4: Dependent RVs with matching marginals to an independent pair
    joint_dep2 = {
        (0, 0): Fraction(1, 6),
        (0, 1): Fraction(1, 3),
        (1, 0): Fraction(1, 3),
        (1, 1): Fraction(1, 6),
    }
    assert check_equivalence(joint_dep2) is False

    print("PASS: Part 2 equivalence (a) <=> (b) verified on independent (with zero-prob y, degenerate Y) and dependent distributions.")


if __name__ == "__main__":
    verify_part1_finite_and_infinite()
    verify_part2()
