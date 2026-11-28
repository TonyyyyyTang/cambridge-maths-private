#!/usr/bin/env python3
"""Independent verification for IB/2004/P2/q17B."""

from fractions import Fraction
import math
import random


def delta(i, j):
    return 1.0 if i == j else 0.0


def inertia_discrete(masses, positions, P):
    I = [[0.0] * 3 for _ in range(3)]
    for m, r in zip(masses, positions):
        rp = [r[i] - P[i] for i in range(3)]
        rp2 = sum(c * c for c in rp)
        for i in range(3):
            for j in range(3):
                I[i][j] += m * (rp2 * delta(i, j) - rp[i] * rp[j])
    return I


def check_parallel_axis_theorem():
    rng = random.Random(200417)
    max_err = 0.0
    for _ in range(20):
        n = 15
        masses = [rng.uniform(0.5, 4.0) for _ in range(n)]
        M = sum(masses)
        positions = [[rng.uniform(-3.0, 3.0) for _ in range(3)] for _ in range(n)]
        G = [sum(masses[k] * positions[k][i] for k in range(n)) / M for i in range(3)]
        P = [rng.uniform(-4.0, 4.0) for _ in range(3)]
        X = [P[i] - G[i] for i in range(3)]  # vector GP
        X2 = sum(c * c for c in X)

        I_P = inertia_discrete(masses, positions, P)
        I_G = inertia_discrete(masses, positions, G)
        for i in range(3):
            for j in range(3):
                rhs = I_G[i][j] + M * (X2 * delta(i, j) - X[i] * X[j])
                err = abs(I_P[i][j] - rhs)
                max_err = max(max_err, err)
                assert err < 1e-11
    print(f"[PASS] Parallel-axis theorem verified on random rigid bodies (max err = {max_err:.2e})")


def gauss_legendre_3pt():
    # Exact for polynomials up to degree 5 on [-1, 1]
    s = math.sqrt(3.0 / 5.0)
    return [-s, 0.0, s], [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]


def check_cube_inertia_quadrature():
    nodes, weights = gauss_legendre_3pt()
    max_err_G = 0.0
    max_err_P = 0.0
    for a, rho in ((1.0, 1.0), (2.3, 0.7), (0.5, 4.2)):
        M = 8.0 * (a ** 3) * rho
        I_G_num = [[0.0] * 3 for _ in range(3)]
        I_P_num = [[0.0] * 3 for _ in range(3)]
        P = [a, a, a]
        for xi, wx in zip(nodes, weights):
            x = a * xi
            for yj, wy in zip(nodes, weights):
                y = a * yj
                for zk, wz in zip(nodes, weights):
                    z = a * zk
                    dV = (a ** 3) * wx * wy * wz
                    r = [x, y, z]
                    r2 = x * x + y * y + z * z
                    rp = [x - P[0], y - P[1], z - P[2]]
                    rp2 = sum(c * c for c in rp)
                    for i in range(3):
                        for j in range(3):
                            I_G_num[i][j] += rho * (r2 * delta(i, j) - r[i] * r[j]) * dV
                            I_P_num[i][j] += rho * (rp2 * delta(i, j) - rp[i] * rp[j]) * dV

        for i in range(3):
            for j in range(3):
                exp_G = (2.0 / 3.0) * M * a * a if i == j else 0.0
                exp_P = (8.0 / 3.0) * M * a * a if i == j else -M * a * a
                max_err_G = max(max_err_G, abs(I_G_num[i][j] - exp_G))
                max_err_P = max(max_err_P, abs(I_P_num[i][j] - exp_P))
                assert abs(I_G_num[i][j] - exp_G) < 1e-12
                assert abs(I_P_num[i][j] - exp_P) < 1e-12
    print(
        f"[PASS] Cube I_ij(G) and I_ij(P) verified by exact 3D Gauss quadrature "
        f"(err_G={max_err_G:.2e}, err_P={max_err_P:.2e})"
    )


def check_eigen_exact():
    # Exact rational matrix I(P) / (M a^2)
    A = [
        [Fraction(8, 3), Fraction(-1), Fraction(-1)],
        [Fraction(-1), Fraction(8, 3), Fraction(-1)],
        [Fraction(-1), Fraction(-1), Fraction(8, 3)],
    ]

    def mv(M, v):
        return [sum(M[i][j] * v[j] for j in range(3)) for i in range(3)]

    def dot_f(u, v):
        return sum(u[i] * v[i] for i in range(3))

    v1 = [Fraction(1), Fraction(1), Fraction(1)]
    v2 = [Fraction(1), Fraction(-1), Fraction(0)]
    v3 = [Fraction(1), Fraction(1), Fraction(-2)]

    lam1 = Fraction(2, 3)
    lam2 = Fraction(11, 3)

    assert mv(A, v1) == [lam1 * c for c in v1]
    assert mv(A, v2) == [lam2 * c for c in v2]
    assert mv(A, v3) == [lam2 * c for c in v3]

    # Mutual orthogonality
    assert dot_f(v1, v2) == 0
    assert dot_f(v1, v3) == 0
    assert dot_f(v2, v3) == 0

    # Check trace and determinant
    tr = sum(A[i][i] for i in range(3))
    assert tr == lam1 + 2 * lam2 == Fraction(8)
    det = (
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
    )
    assert det == lam1 * lam2 * lam2 == Fraction(242, 27)

    # Check characteristic polynomial det(A - mu I) == (2/3 - mu)(11/3 - mu)^2 at several rational points
    for mu in (Fraction(0), Fraction(1), Fraction(2, 3), Fraction(11, 3), Fraction(-5, 2), Fraction(7)):
        B = [[A[i][j] - (mu if i == j else 0) for j in range(3)] for i in range(3)]
        dB = (
            B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
            - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0])
            + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0])
        )
        assert dB == (Fraction(2, 3) - mu) * ((Fraction(11, 3) - mu) ** 2)

    # Check orthonormal diagonalisation Q^T A Q = diag(2/3, 11/3, 11/3)
    e1 = [1.0 / math.sqrt(3.0)] * 3
    e2 = [1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0), 0.0]
    e3 = [1.0 / math.sqrt(6.0), 1.0 / math.sqrt(6.0), -2.0 / math.sqrt(6.0)]
    Q = [[e1[i], e2[i], e3[i]] for i in range(3)]
    Af = [[float(A[i][j]) for j in range(3)] for i in range(3)]
    AQ = [[sum(Af[i][k] * Q[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    D = [[sum(Q[k][i] * AQ[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    expected_D = [2.0 / 3.0, 11.0 / 3.0, 11.0 / 3.0]
    for i in range(3):
        for j in range(3):
            target = expected_D[i] if i == j else 0.0
            assert abs(D[i][j] - target) < 1e-14

    print("[PASS] Exact rational eigenvalues (2/3, 11/3, 11/3)*Ma^2, eigenvectors, charpoly & diagonalisation verified")


if __name__ == "__main__":
    check_parallel_axis_theorem()
    check_cube_inertia_quadrature()
    check_eigen_exact()
    print("ALL CHECKS PASSED FOR IB/2004/P2/q17B")
