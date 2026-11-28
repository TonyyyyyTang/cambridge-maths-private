#!/usr/bin/env python3
"""Independent verification for IB/2003/P2/q11C."""

import math
import random


def dot(u, v):
    return sum(u[i] * v[i] for i in range(3))


def det3(M):
    return (
        M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
        - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
        + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
    )


def random_orthogonal(rng, proper=True):
    cols = []
    for _ in range(3):
        v = [rng.gauss(0.0, 1.0) for _ in range(3)]
        for u in cols:
            proj = dot(v, u)
            v = [v[i] - proj * u[i] for i in range(3)]
        nrm = math.sqrt(dot(v, v))
        cols.append([x / nrm for x in v])
    R = [[cols[j][i] for j in range(3)] for i in range(3)]
    d = det3(R)
    if (proper and d < 0) or ((not proper) and d > 0):
        R[0] = [-x for x in R[0]]
    return R


def delta(i, j):
    return 1.0 if i == j else 0.0


def build_c(alpha, beta, gamma):
    return [
        [
            [
                [
                    alpha * delta(i, j) * delta(k, l)
                    + beta * delta(i, k) * delta(j, l)
                    + gamma * delta(i, l) * delta(j, k)
                    for l in range(3)
                ]
                for k in range(3)
            ]
            for j in range(3)
        ]
        for i in range(3)
    ]


def check_isotropy():
    rng = random.Random(200311)
    max_err = 0.0
    for proper in (True, False):
        for _ in range(10):
            alpha = rng.uniform(-4.0, 4.0)
            beta = rng.uniform(-4.0, 4.0)
            gamma = rng.uniform(-4.0, 4.0)
            c = build_c(alpha, beta, gamma)
            R = random_orthogonal(rng, proper=proper)
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            val = 0.0
                            for p in range(3):
                                rip = R[i][p]
                                for q in range(3):
                                    rip_rjq = rip * R[j][q]
                                    for r in range(3):
                                        rip_rjq_rkr = rip_rjq * R[k][r]
                                        for s in range(3):
                                            val += rip_rjq_rkr * R[l][s] * c[p][q][r][s]
                            err = abs(val - c[i][j][k][l])
                            max_err = max(max_err, err)
                            assert err < 1e-12
    print(f"[PASS] Isotropy of c_ijkl verified over O(3) (max err = {max_err:.2e})")


def check_stress_strain_and_decomposition():
    rng = random.Random(200312)
    max_sigma_err = 0.0
    max_E_err = 0.0
    for _ in range(50):
        alpha = rng.uniform(-5.0, 5.0)
        beta = rng.uniform(-5.0, 5.0)
        gamma = rng.uniform(-5.0, 5.0)
        lam = alpha
        mu = 0.5 * (beta + gamma)
        c = build_c(alpha, beta, gamma)

        # Random symmetric strain tensor e_ij
        A = [[rng.uniform(-3.0, 3.0) for _ in range(3)] for _ in range(3)]
        e = [[0.5 * (A[i][j] + A[j][i]) for j in range(3)] for i in range(3)]

        # sigma via c_ijkl e_kl
        sigma1 = [
            [sum(c[i][j][k][l] * e[k][l] for k in range(3) for l in range(3)) for j in range(3)]
            for i in range(3)
        ]
        # sigma via lambda e_kk delta_ij + 2 mu e_ij
        tr_e = sum(e[k][k] for k in range(3))
        sigma2 = [[lam * tr_e * delta(i, j) + 2.0 * mu * e[i][j] for j in range(3)] for i in range(3)]

        for i in range(3):
            for j in range(3):
                err = abs(sigma1[i][j] - sigma2[i][j])
                max_sigma_err = max(max_sigma_err, err)
                assert err < 1e-12

        # Decomposition e_ij = p delta_ij + d_ij
        p = tr_e / 3.0
        d = [[e[i][j] - p * delta(i, j) for j in range(3)] for i in range(3)]
        assert abs(sum(d[i][i] for i in range(3))) < 1e-13
        for i in range(3):
            for j in range(3):
                assert abs(d[i][j] - d[j][i]) < 1e-14
                assert abs(p * delta(i, j) + d[i][j] - e[i][j]) < 1e-14

        # Energy density E = (1/2) sigma_ij e_ij
        E_direct = 0.5 * sum(sigma2[i][j] * e[i][j] for i in range(3) for j in range(3))
        d_sq = sum(d[i][j] * d[i][j] for i in range(3) for j in range(3))
        E_decomp = 1.5 * (3.0 * lam + 2.0 * mu) * (p * p) + mu * d_sq
        E_trace = 0.5 * (lam + (2.0 / 3.0) * mu) * (tr_e * tr_e) + mu * d_sq
        max_E_err = max(max_E_err, abs(E_direct - E_decomp), abs(E_direct - E_trace))
        assert abs(E_direct - E_decomp) < 1e-12
        assert abs(E_direct - E_trace) < 1e-12

    print(
        f"[PASS] Stress-strain form & energy decomposition verified "
        f"(max sigma err={max_sigma_err:.2e}, max E err={max_E_err:.2e})"
    )


def energy(e, lam, mu):
    tr_e = sum(e[k][k] for k in range(3))
    return 0.5 * sum((lam * tr_e * delta(i, j) + 2.0 * mu * e[i][j]) * e[i][j] for i in range(3) for j in range(3))


def check_positivity_conditions():
    rng = random.Random(200313)

    # 1. Sufficiency: test interior and boundary points of {mu >= 0, lam >= -2/3 mu}
    valid_pairs = [
        (0.0, 0.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (-2.0 / 3.0, 1.0),
        (-2.0, 3.0),
        (2.5, 1.7),
        (-0.5, 1.0),
    ]
    for lam, mu in valid_pairs:
        for _ in range(100):
            A = [[rng.uniform(-3.0, 3.0) for _ in range(3)] for _ in range(3)]
            e = [[0.5 * (A[i][j] + A[j][i]) for j in range(3)] for i in range(3)]
            assert energy(e, lam, mu) >= -1e-13

    # 2. Necessity: whenever mu < 0 or lam < -2/3 mu, explicit test strain gives E < 0
    pure_shear = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    pure_dilation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

    invalid_pairs = [
        (5.0, -0.01),      # mu < 0, lam + 2/3 mu > 0
        (-1.0, -0.5),      # both violated
        (-0.67, 1.0),      # mu > 0, lam < -2/3 mu
        (-1e-6, 0.0),      # mu = 0, lam < 0
    ]
    for lam, mu in invalid_pairs:
        e_shear = energy(pure_shear, lam, mu)
        e_dil = energy(pure_dilation, lam, mu)
        assert abs(e_shear - 2.0 * mu) < 1e-13
        assert abs(e_dil - 1.5 * (3.0 * lam + 2.0 * mu)) < 1e-13
        assert min(e_shear, e_dil) < 0.0

    print("[PASS] Necessity & sufficiency of mu >= 0 and lambda >= -2/3 mu verified")


if __name__ == "__main__":
    check_isotropy()
    check_stress_strain_and_decomposition()
    check_positivity_conditions()
    print("ALL CHECKS PASSED FOR IB/2003/P2/q11C")
