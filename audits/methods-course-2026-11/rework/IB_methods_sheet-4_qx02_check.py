#!/usr/bin/env python3
"""Independent verification for IB/methods/sheet-4/qx02 (Tensor transformations)."""

from fractions import Fraction
import math
import random


def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def mat_transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def mat_det3(A):
    return (
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
    )


def mat_trace(A):
    return sum(A[i][i] for i in range(len(A)))


def transform_rank2(L, T):
    """Compute T'_{ij} = l_{ip} l_{jq} T_{pq} = (L T L^T)_{ij}."""
    return mat_mul(mat_mul(L, T), mat_transpose(L))


def build_frame_rotation(n):
    """
    Construct L^{(n+1)} where the i-th axis of frame 1 coincides with
    the (i+n)-th axis (mod 3, 1-indexed) of frame (n+1).
    l_{r, i} = e^{(n+1)}_r . e^{(1)}_i = 1 iff r == ((i - 1 + n) % 3) + 1.
    """
    L = [[Fraction(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):  # 0-indexed axis of frame 1
        r = (i + n) % 3  # 0-indexed axis of frame (n+1)
        L[r][i] = Fraction(1)
    return L


def random_so3(rng):
    """Generate a random proper orthogonal 3x3 matrix via Gram-Schmidt."""
    while True:
        cols = [[rng.gauss(0.0, 1.0) for _ in range(3)] for _ in range(3)]
        ortho = []
        ok = True
        for v in cols:
            w = v[:]
            for u in ortho:
                dot = sum(w[k] * u[k] for k in range(3))
                w = [w[k] - dot * u[k] for k in range(3)]
            norm = math.sqrt(sum(x * x for x in w))
            if norm < 1e-8:
                ok = False
                break
            ortho.append([x / norm for x in w])
        if not ok:
            continue
        L = [[ortho[j][i] for j in range(3)] for i in range(3)]
        if mat_det3(L) < 0:
            # Flip last row to ensure det(L) = +1
            L[2] = [-x for x in L[2]]
        return L


def main():
    I3 = [[Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]

    # Build L^{(1)}, L^{(2)}, L^{(3)}
    Ls = [build_frame_rotation(n) for n in range(3)]
    expected_L2 = [
        [Fraction(0), Fraction(0), Fraction(1)],
        [Fraction(1), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(1), Fraction(0)],
    ]
    expected_L3 = [
        [Fraction(0), Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(1)],
        [Fraction(1), Fraction(0), Fraction(0)],
    ]
    assert Ls[0] == I3
    assert Ls[1] == expected_L2
    assert Ls[2] == expected_L3

    for idx, L in enumerate(Ls, start=1):
        assert mat_mul(L, mat_transpose(L)) == I3, f"L^({idx}) not orthogonal"
        assert mat_mul(mat_transpose(L), L) == I3, f"L^({idx}) not orthogonal"
        assert mat_det3(L) == Fraction(1), f"L^({idx}) not proper rotation"
    print("OK: Frame rotation matrices L^(1), L^(2), L^(3) are exact SO(3) cyclic permutations.")

    # Part (i): given arrays
    A1 = [[Fraction(1 if i == j else 0) for j in range(3)] for i in range(3)]
    A2 = [[Fraction(1 if (i == j and i < 2) else 0) for j in range(3)] for i in range(3)]
    A3 = [[Fraction(1 if (i == j and i < 1) else 0) for j in range(3)] for i in range(3)]

    A2_pred = transform_rank2(Ls[1], A1)
    A3_pred = transform_rank2(Ls[2], A1)
    assert A2_pred == I3 and A2_pred != A2
    assert A3_pred == I3 and A3_pred != A3
    assert (mat_trace(A1), mat_trace(A2), mat_trace(A3)) == (Fraction(3), Fraction(2), Fraction(1))
    print("OK: Part (i) verified — L A^(1) L^T = I != A^(2), A^(3); traces 3, 2, 1 disagree.")

    # Part (ii): given arrays
    T1 = [[Fraction(1 if (i == 0 and j == 0) else 0) for j in range(3)] for i in range(3)]
    T2 = [[Fraction(1 if (i == 1 and j == 1) else 0) for j in range(3)] for i in range(3)]
    T3 = [[Fraction(1 if (i == 2 and j == 2) else 0) for j in range(3)] for i in range(3)]

    assert transform_rank2(Ls[1], T1) == T2
    assert transform_rank2(Ls[2], T1) == T3
    # Also check Frame 2 -> Frame 3 directly: L^{(3<-2)} = L^{(3)} (L^{(2)})^T
    L32 = mat_mul(Ls[2], mat_transpose(Ls[1]))
    assert L32 == Ls[1]
    assert transform_rank2(L32, T2) == T3
    print("OK: Part (ii) consistency across frames 1, 2, 3 verified exactly.")

    # Part (ii): arbitrary Cartesian frame with direction cosines (lambda, mu, nu) to 1-axis of frame 1
    T1_f = [[float(x) for x in row] for row in T1]
    rng = random.Random(20261127)
    max_err = 0.0
    max_proj_err = 0.0
    max_tr_err = 0.0
    for _ in range(50):
        L = random_so3(rng)
        # Direction cosines of new axes (e'_1, e'_2, e'_3) with respect to 1-axis of frame 1 (e^{(1)}_1)
        lam, mu, nu = L[0][0], L[1][0], L[2][0]
        assert abs(lam * lam + mu * mu + nu * nu - 1.0) < 1e-13

        T_prime = transform_rank2(L, T1_f)
        T_formula = [
            [lam * lam, lam * mu, lam * nu],
            [lam * mu, mu * mu, mu * nu],
            [lam * nu, mu * nu, nu * nu],
        ]
        for i in range(3):
            for j in range(3):
                err = abs(T_prime[i][j] - T_formula[i][j])
                assert err < 1e-13
                max_err = max(max_err, err)
                # Index formula l_{i1} l_{j1}
                assert abs(T_prime[i][j] - L[i][0] * L[j][0]) < 1e-13

        # Trace = 1 and idempotence T'^2 = T'
        max_tr_err = max(max_tr_err, abs(mat_trace(T_prime) - 1.0))
        T_sq = mat_mul(T_prime, T_prime)
        for i in range(3):
            for j in range(3):
                max_proj_err = max(max_proj_err, abs(T_sq[i][j] - T_prime[i][j]))

    assert max_err < 1e-13
    assert max_tr_err < 1e-13
    assert max_proj_err < 1e-13
    print(
        f"OK: 50 random SO(3) frames verified — max formula err={max_err:.2e}, "
        f"max trace err={max_tr_err:.2e}, max idempotence err={max_proj_err:.2e}."
    )


if __name__ == "__main__":
    main()
