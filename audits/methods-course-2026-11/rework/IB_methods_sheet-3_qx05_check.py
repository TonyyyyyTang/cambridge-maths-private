#!/usr/bin/env python3
"""Verification script for IB/methods/sheet-3/qx05 (Stationary points)."""

import math
import random


def phi(x1, x2, x3):
    return 0.25 * (x1**4 + x2**4 + x3**4) - x2 * x3 - x3 * x1 - x1 * x2


def phi_sos(x1, x2, x3):
    return (
        -3.0
        + 0.25 * ((x1**2 - 2.0) ** 2 + (x2**2 - 2.0) ** 2 + (x3**2 - 2.0) ** 2)
        + 0.5 * ((x1 - x2) ** 2 + (x2 - x3) ** 2 + (x3 - x1) ** 2)
    )


def grad_phi(x1, x2, x3):
    return (x1**3 - x2 - x3, x2**3 - x3 - x1, x3**3 - x1 - x2)


def hess_phi(x1, x2, x3):
    return [
        [3.0 * x1**2, -1.0, -1.0],
        [-1.0, 3.0 * x2**2, -1.0],
        [-1.0, -1.0, 3.0 * x3**2],
    ]


def matvec(M, v):
    return tuple(sum(M[i][j] * v[j] for j in range(3)) for i in range(3))


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def norm(v):
    return math.sqrt(dot(v, v))


def main():
    random.seed(42)

    # 1. Verify sum-of-squares identity and global lower bound -3
    max_id_err = 0.0
    for _ in range(5000):
        x1 = random.uniform(-5.0, 5.0)
        x2 = random.uniform(-5.0, 5.0)
        x3 = random.uniform(-5.0, 5.0)
        val = phi(x1, x2, x3)
        sos = phi_sos(x1, x2, x3)
        max_id_err = max(max_id_err, abs(val - sos))
        assert val >= -3.0 - 1e-12
    assert max_id_err < 1e-11, f"SOS identity error {max_id_err}"
    print(f"[PASS] Sum-of-squares identity verified on 5000 random points (max err = {max_id_err:.2e})")

    # 2. Verify all three stationary points and finite-difference gradients/Hessians
    s2 = math.sqrt(2.0)
    stat_pts = [
        ("min_+", (s2, s2, s2), -3.0),
        ("min_-", (-s2, -s2, -s2), -3.0),
        ("saddle", (0.0, 0.0, 0.0), 0.0),
    ]
    for name, pt, expected_val in stat_pts:
        val = phi(*pt)
        g = grad_phi(*pt)
        assert abs(val - expected_val) < 1e-13, f"{name}: val={val}"
        assert norm(g) < 1e-13, f"{name}: grad={g}"

        # Finite-difference check on gradient and Hessian
        h = 1e-5
        H_exact = hess_phi(*pt)
        for i in range(3):
            e_i = [0.0, 0.0, 0.0]
            e_i[i] = h
            fd_g = (phi(*(pt[k] + e_i[k] for k in range(3))) - phi(*(pt[k] - e_i[k] for k in range(3)))) / (2 * h)
            assert abs(fd_g - g[i]) < 1e-8
            for j in range(3):
                e_j = [0.0, 0.0, 0.0]
                e_j[j] = h
                f_pp = phi(*(pt[k] + e_i[k] + e_j[k] for k in range(3)))
                f_pm = phi(*(pt[k] + e_i[k] - e_j[k] for k in range(3)))
                f_mp = phi(*(pt[k] - e_i[k] + e_j[k] for k in range(3)))
                f_mm = phi(*(pt[k] - e_i[k] - e_j[k] for k in range(3)))
                fd_H = (f_pp - f_pm - f_mp + f_mm) / (4 * h * h)
                assert abs(fd_H - H_exact[i][j]) < 1e-5
    print("[PASS] Stationary points, values, gradients, and Hessians verified")

    # 3. Verify no other stationary points on a dense 3D grid in [-3, 3]^3
    # (Outside [-3, 3]^3, if max|x_i| = M > 3, then M^3 >= 27 > 6 >= |x_j + x_k|?
    # More precisely, for M = max|x_i| > sqrt(2), |x_i^3| = M^3 > 2M >= |x_j + x_k|, so grad != 0!)
    n_grid = 61
    min_grad_outside = float("inf")
    for i in range(n_grid):
        x1 = -2.0 + 4.0 * i / (n_grid - 1)
        for j in range(n_grid):
            x2 = -2.0 + 4.0 * j / (n_grid - 1)
            for k in range(n_grid):
                x3 = -2.0 + 4.0 * k / (n_grid - 1)
                if min(norm((x1 - p[0], x2 - p[1], x3 - p[2])) for _, p, _ in stat_pts) < 0.15:
                    continue
                min_grad_outside = min(min_grad_outside, norm(grad_phi(x1, x2, x3)))
    assert min_grad_outside > 0.05, f"Unexpected near-stationary point: {min_grad_outside}"
    print(f"[PASS] No extraneous stationary points on [-2,2]^3 grid (min |grad| outside 0.15-balls = {min_grad_outside:.4f})")

    # 4. Verify eigenpairs of Hessians at minima and saddle
    e1 = (1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0))
    e2 = (1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0), 0.0)
    e3 = (1.0 / math.sqrt(6.0), 1.0 / math.sqrt(6.0), -2.0 / math.sqrt(6.0))
    basis = [e1, e2, e3]
    for i in range(3):
        for j in range(3):
            assert abs(dot(basis[i], basis[j]) - (1.0 if i == j else 0.0)) < 1e-14

    H_min = hess_phi(s2, s2, s2)
    for vec, lam in zip(basis, [4.0, 7.0, 7.0]):
        res = norm(tuple(matvec(H_min, vec)[k] - lam * vec[k] for k in range(3)))
        assert res < 1e-13
    H_sad = hess_phi(0.0, 0.0, 0.0)
    for vec, lam in zip(basis, [-2.0, 1.0, 1.0]):
        res = norm(tuple(matvec(H_sad, vec)[k] - lam * vec[k] for k in range(3)))
        assert res < 1e-13
    print("[PASS] Hessian eigenpairs verified: minima eigenvalues (4, 7, 7), saddle eigenvalues (-2, 1, 1)")

    # 5. Verify tangent double cone at (0,0,0) with semi-angle alpha = atan(sqrt(2))
    alpha = math.atan(math.sqrt(2.0))
    # Coordinate axes lie on cone at angle alpha from e1
    for axis in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]:
        theta_axis = math.acos(dot(axis, e1))
        assert abs(theta_axis - alpha) < 1e-14
        Q_axis = 0.5 * dot(axis, matvec(H_sad, axis))
        assert abs(Q_axis) < 1e-14

    # Points on surface phi = 0 approaching the origin along both nappes
    max_phi_err = 0.0
    max_normal_angle_ratio = 0.0
    for nappe_sign in (+1, -1):
        for eps in [1e-2, 1e-3, 1e-4, 1e-5]:
            theta = (alpha - eps) if nappe_sign == +1 else (math.pi - alpha + eps)
            for m in range(24):
                psi = 2.0 * math.pi * m / 24.0
                n_hat = tuple(
                    math.cos(theta) * e1[k]
                    + math.sin(theta) * (math.cos(psi) * e2[k] + math.sin(psi) * e3[k])
                    for k in range(3)
                )
                Q_n = 0.5 * dot(n_hat, matvec(H_sad, n_hat))
                assert Q_n < 0.0
                quart = sum(c**4 for c in n_hat)
                r = math.sqrt(-4.0 * Q_n / quart)
                x = tuple(r * c for c in n_hat)
                max_phi_err = max(max_phi_err, abs(phi(*x)))

                # Compare surface normal grad_phi(x) with cone normal grad_Q(x) = H_sad * x
                g_phi = grad_phi(*x)
                g_Q = matvec(H_sad, x)
                cos_ang = dot(g_phi, g_Q) / (norm(g_phi) * norm(g_Q))
                cos_ang = max(-1.0, min(1.0, cos_ang))
                ang = math.acos(cos_ang)
                max_normal_angle_ratio = max(max_normal_angle_ratio, ang / (r * r))
    assert max_phi_err < 1e-15
    assert max_normal_angle_ratio < 2.0
    print(
        f"[PASS] Surface phi=0 approaches both nappes of cone (semi-angle atan(sqrt(2))={alpha:.6f} rad): "
        f"max |phi|={max_phi_err:.2e}, max normal angle/r^2={max_normal_angle_ratio:.4f}"
    )
    print("ALL CHECKS PASSED for IB_methods_sheet-3_qx05.")


if __name__ == "__main__":
    main()
