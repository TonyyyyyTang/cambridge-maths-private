#!/usr/bin/env python3
"""Independent verification for IB/2002/P2/q11C."""

import math
import random


def mat_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mat_vec(A, v):
    return [sum(A[i][j] * v[j] for j in range(3)) for i in range(3)]


def dot(u, v):
    return sum(u[i] * v[i] for i in range(3))


def cross(u, v):
    return [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]


def T_tensor(B, alpha, beta):
    B2 = dot(B, B)
    return [
        [alpha * B[i] * B[j] + (beta * B2 if i == j else 0.0) for j in range(3)]
        for i in range(3)
    ]


def det3(M):
    return (
        M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
        - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
        + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
    )


def random_orthogonal(rng, proper=True):
    # Gram-Schmidt on 3 random vectors
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


def check_tensor_and_symmetry():
    rng = random.Random(200211)
    max_err = 0.0
    for proper in (True, False):
        for _ in range(25):
            R = random_orthogonal(rng, proper=proper)
            B = [rng.uniform(-3.0, 3.0) for _ in range(3)]
            alpha = rng.uniform(-5.0, 5.0)
            beta = rng.uniform(-5.0, 5.0)
            T = T_tensor(B, alpha, beta)
            # Symmetry
            for i in range(3):
                for j in range(3):
                    assert abs(T[i][j] - T[j][i]) < 1e-14
            # Transformation as polar vector and as axial vector (pseudovector)
            for sign in (1.0, det3(R)):
                B_prime = [sign * x for x in mat_vec(R, B)]
                T_prime_direct = T_tensor(B_prime, alpha, beta)
                # R T R^T
                RT = mat_mul(R, T)
                R_T = [[R[j][i] for j in range(3)] for i in range(3)]
                T_prime_law = mat_mul(RT, R_T)
                for i in range(3):
                    for j in range(3):
                        err = abs(T_prime_direct[i][j] - T_prime_law[i][j])
                        max_err = max(max_err, err)
                        assert err < 1e-12
    print(f"[PASS] Tensor transformation & symmetry checked (max err = {max_err:.2e})")


def check_eigenvalues_eigenvectors():
    rng = random.Random(200212)
    max_err = 0.0
    for _ in range(30):
        B = [rng.uniform(-3.0, 3.0) for _ in range(3)]
        if dot(B, B) < 1e-4:
            B = [1.0, 2.0, -1.5]
        alpha = rng.uniform(-5.0, 5.0)
        beta = rng.uniform(-5.0, 5.0)
        B2 = dot(B, B)
        T = T_tensor(B, alpha, beta)

        # Eigenvalue 1: (alpha + beta) B^2, eigenvector B
        lam1 = (alpha + beta) * B2
        Tv1 = mat_vec(T, B)
        for i in range(3):
            err = abs(Tv1[i] - lam1 * B[i])
            max_err = max(max_err, err)
            assert err < 1e-12

        # Build two orthonormal vectors perpendicular to B
        ref = [1.0, 0.0, 0.0] if abs(B[0]) < 0.9 * math.sqrt(B2) else [0.0, 1.0, 0.0]
        e2 = cross(B, ref)
        n2 = math.sqrt(dot(e2, e2))
        e2 = [x / n2 for x in e2]
        e3 = cross(B, e2)
        n3 = math.sqrt(dot(e3, e3))
        e3 = [x / n3 for x in e3]

        lam2 = beta * B2
        for v in (e2, e3):
            assert abs(dot(B, v)) < 1e-13
            Tv = mat_vec(T, v)
            for i in range(3):
                err = abs(Tv[i] - lam2 * v[i])
                max_err = max(max_err, err)
                assert err < 1e-12

        # Trace and determinant checks
        tr = sum(T[i][i] for i in range(3))
        assert abs(tr - (lam1 + 2.0 * lam2)) < 1e-12
        assert abs(det3(T) - lam1 * lam2 * lam2) < 1e-10

    # Degenerate cases: alpha = 0 and B = 0
    T_alpha0 = T_tensor([1.0, -2.0, 3.0], 0.0, 2.5)
    for v in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.3, -0.7, 1.1]):
        Tv = mat_vec(T_alpha0, v)
        for i in range(3):
            assert abs(Tv[i] - 2.5 * 14.0 * v[i]) < 1e-12
    T_B0 = T_tensor([0.0, 0.0, 0.0], 3.0, -4.0)
    assert all(T_B0[i][j] == 0.0 for i in range(3) for j in range(3))
    print(f"[PASS] Eigenvalues, eigenvectors & degenerate cases checked (max err = {max_err:.2e})")


def curl_field(F, x, y, z, h=1e-5):
    def df(comp, axis):
        p = [x, y, z]
        p_p = list(p)
        p_m = list(p)
        p_p2 = list(p)
        p_m2 = list(p)
        p_p[axis] += h
        p_m[axis] -= h
        p_p2[axis] += 2 * h
        p_m2[axis] -= 2 * h
        return (-F(*p_p2)[comp] + 8 * F(*p_p)[comp] - 8 * F(*p_m)[comp] + F(*p_m2)[comp]) / (12 * h)

    return [
        df(2, 1) - df(1, 2),
        df(0, 2) - df(2, 0),
        df(1, 0) - df(0, 1),
    ]


def div_field(F, x, y, z, h=1e-5):
    s = 0.0
    for axis in range(3):
        p_p = [x, y, z]
        p_m = [x, y, z]
        p_p2 = [x, y, z]
        p_m2 = [x, y, z]
        p_p[axis] += h
        p_m[axis] -= h
        p_p2[axis] += 2 * h
        p_m2[axis] -= 2 * h
        s += (-F(*p_p2)[axis] + 8 * F(*p_p)[axis] - 8 * F(*p_m)[axis] + F(*p_m2)[axis]) / (12 * h)
    return s


def check_divergence_identity():
    # Construct a genuinely 3D, non-symmetric divergence-free field B = curl(A)
    def A(x, y, z):
        return [
            math.sin(x * y) + 0.3 * z * z * x,
            math.cos(y + 2 * z) + 0.5 * x * x * y,
            math.exp(0.2 * x) * math.sin(y) + x * y * z,
        ]

    # Analytic B = curl(A)
    def B(x, y, z):
        dAz_dy = math.exp(0.2 * x) * math.cos(y) + x * z
        dAy_dz = -2.0 * math.sin(y + 2 * z)
        dAx_dz = 0.6 * z * x
        dAz_dx = 0.2 * math.exp(0.2 * x) * math.sin(y) + y * z
        dAy_dx = 1.0 * x * y
        dAx_dy = x * math.cos(x * y)
        return [dAz_dy - dAy_dz, dAx_dz - dAz_dx, dAy_dx - dAx_dy]

    alpha, beta = 1.0, -0.5
    rng = random.Random(200213)
    max_divB = 0.0
    max_err = 0.0
    h = 1e-4

    for _ in range(20):
        x, y, z = [rng.uniform(-1.5, 1.5) for _ in range(3)]
        max_divB = max(max_divB, abs(div_field(B, x, y, z, h)))
        assert max_divB < 1e-10

        # LHS: d/dx_j T_{ij} using 4th-order central finite differences
        lhs = [0.0, 0.0, 0.0]
        for j in range(3):
            p_p = [x, y, z]
            p_m = [x, y, z]
            p_p2 = [x, y, z]
            p_m2 = [x, y, z]
            p_p[j] += h
            p_m[j] -= h
            p_p2[j] += 2 * h
            p_m2[j] -= 2 * h
            Tp = T_tensor(B(*p_p), alpha, beta)
            Tm = T_tensor(B(*p_m), alpha, beta)
            Tp2 = T_tensor(B(*p_p2), alpha, beta)
            Tm2 = T_tensor(B(*p_m2), alpha, beta)
            for i in range(3):
                lhs[i] += (-Tp2[i][j] + 8 * Tp[i][j] - 8 * Tm[i][j] + Tm2[i][j]) / (12 * h)

        # RHS: (curl B) x B
        rhs = cross(curl_field(B, x, y, z, h), B(x, y, z))
        for i in range(3):
            err = abs(lhs[i] - rhs[i])
            max_err = max(max_err, err)
            assert err < 1e-8

    print(
        f"[PASS] d_j T_ij = [(curl B) x B]_i for alpha=1, beta=-1/2 "
        f"(max |div B|={max_divB:.2e}, max err={max_err:.2e})"
    )


def gauss_legendre_nodes_weights(n):
    # Compute n-point Gauss-Legendre nodes & weights on [-1, 1] via Newton's method
    nodes = []
    weights = []
    for k in range(1, n + 1):
        x = math.cos(math.pi * (k - 0.25) / (n + 0.5))
        for _ in range(30):
            p0 = 1.0
            p1 = x
            for j in range(2, n + 1):
                p2 = ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
                p0, p1 = p1, p2
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-15:
                break
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * dp * dp))
    return nodes, weights


def check_volume_integral():
    # Smooth divergence-free B = curl(A) vanishing on boundary of [-1, 1]^3:
    # A(x,y,z) = phi(x,y,z) * g(x,y,z) where phi = (1-x^2)^2 (1-y^2)^2 (1-z^2)^2
    # Since phi and grad(phi) both vanish on boundary, B = curl(A) = 0 identically on d[-1,1]^3.
    # Choose asymmetric g so integrand (curl B) x B is not trivially odd!
    def A(x, y, z):
        phi = ((1.0 - x * x) * (1.0 - y * y) * (1.0 - z * z)) ** 2
        return [
            phi * (1.0 + 0.7 * x - 0.4 * y + 0.5 * z),
            phi * (0.6 - 0.3 * x + 0.8 * y * z),
            phi * (-0.5 + 0.9 * x * y + 0.2 * z),
        ]

    # Analytic B = curl A via exact derivatives of phi * g
    def phi_and_grad(x, y, z):
        u = 1.0 - x * x
        v = 1.0 - y * y
        w = 1.0 - z * z
        p = (u * v * w) ** 2
        px = -4.0 * x * u * (v * w) ** 2
        py = -4.0 * y * v * (u * w) ** 2
        pz = -4.0 * z * w * (u * v) ** 2
        return p, px, py, pz

    def B(x, y, z):
        p, px, py, pz = phi_and_grad(x, y, z)
        gx = 1.0 + 0.7 * x - 0.4 * y + 0.5 * z
        gy = 0.6 - 0.3 * x + 0.8 * y * z
        gz = -0.5 + 0.9 * x * y + 0.2 * z
        # dAz/dy - dAy/dz
        Bx = (py * gz + p * 0.9 * x) - (pz * gy + p * 0.8 * y)
        # dAx/dz - dAz/dx
        By = (pz * gx + p * 0.5) - (px * gz + p * 0.9 * y)
        # dAy/dx - dAx/dy
        Bz = (px * gy + p * (-0.3)) - (py * gx + p * (-0.4))
        return [Bx, By, Bz]

    # Verify B vanishes on all 6 faces
    for s in (-1.0, 1.0):
        for u in (-0.3, 0.5):
            for v in (-0.7, 0.2):
                for pt in ((s, u, v), (u, s, v), (u, v, s)):
                    assert max(abs(c) for c in B(*pt)) == 0.0

    # Verify integrand is non-trivial inside
    test_val = cross(curl_field(B, 0.2, -0.3, 0.1, h=1e-4), B(0.2, -0.3, 0.1))
    assert math.sqrt(dot(test_val, test_val)) > 1e-2

    # Integrate (curl B) x B over [-1, 1]^3 using 12-point Gauss-Legendre in each axis
    # (polynomial of degree <= 20 in each variable: 12-pt GL is exact up to degree 23!)
    nodes, weights = gauss_legendre_nodes_weights(12)
    integral = [0.0, 0.0, 0.0]
    abs_norm = 0.0
    for xi, wx in zip(nodes, weights):
        for yj, wy in zip(nodes, weights):
            wxy = wx * wy
            for zk, wz in zip(nodes, weights):
                w = wxy * wz
                bval = B(xi, yj, zk)
                cb = curl_field(B, xi, yj, zk, h=1e-4)
                f = cross(cb, bval)
                for c in range(3):
                    integral[c] += w * f[c]
                abs_norm += w * math.sqrt(dot(f, f))

    assert abs_norm > 0.1
    assert max(abs(c) for c in integral) < 1e-9
    print(
        f"[PASS] Volume integral of (curl B) x B over [-1,1]^3 = "
        f"[{integral[0]:.2e}, {integral[1]:.2e}, {integral[2]:.2e}] (L1 norm = {abs_norm:.3f})"
    )

    # Also verify counterexample when div(B) != 0: B = (1 - r^2)(1, x, 0) on unit ball
    # Exact integral of (curl B) x B is (0, 16*pi/105, 0) != 0
    n_r, n_t = 10, 10
    rn, rw = gauss_legendre_nodes_weights(n_r)
    un, uw = gauss_legendre_nodes_weights(n_t)  # u = cos(theta)
    # phi integral of (1-r^2)(1-r^2-2x^2+2y): since x = r sin(t) cos(p), y = r sin(t) sin(p),
    # integrate over phi with 8-point trapezoidal rule (exact for trig deg <= 7)
    iy = 0.0
    for r_raw, wr in zip(rn, rw):
        r = 0.5 * (r_raw + 1.0)
        wr_s = 0.5 * wr
        for u, wu in zip(un, uw):
            st = math.sqrt(max(0.0, 1.0 - u * u))
            for k in range(8):
                phi = 2.0 * math.pi * k / 8
                x = r * st * math.cos(phi)
                y = r * st * math.sin(phi)
                # (curl B)_z * B_x = (1 - r^2 - 2x^2 + 2y) * (1 - r^2)
                val = (1.0 - r * r) * (1.0 - r * r - 2.0 * x * x + 2.0 * y)
                iy += wr_s * wu * (2.0 * math.pi / 8.0) * (r * r) * val
    expected = 16.0 * math.pi / 105.0
    assert abs(iy - expected) < 1e-12
    print(f"[PASS] Counterexample (div B != 0) confirmed: int_V = {iy:.10f} == 16*pi/105 = {expected:.10f}")


if __name__ == "__main__":
    check_tensor_and_symmetry()
    check_eigenvalues_eigenvectors()
    check_divergence_identity()
    check_volume_integral()
    print("ALL CHECKS PASSED FOR IB/2002/P2/q11C")
