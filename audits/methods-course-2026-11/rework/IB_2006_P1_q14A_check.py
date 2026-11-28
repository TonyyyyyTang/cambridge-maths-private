#!/usr/bin/env python3
"""Independent computational verification for IB/2006/P1/q14A."""

import math


_GL8_X = [
    -0.96028985649753623168,
    -0.79666647741362673959,
    -0.52553240991632898582,
    -0.18343464249564980494,
    0.18343464249564980494,
    0.52553240991632898582,
    0.79666647741362673959,
    0.96028985649753623168,
]
_GL8_W = [
    0.10122853629037625915,
    0.22238103445337447054,
    0.31370664587788728734,
    0.36268378337836198297,
    0.36268378337836198297,
    0.31370664587788728734,
    0.22238103445337447054,
    0.10122853629037625915,
]


def quad1d(f, a, b, panels=12):
    h = (b - a) / panels
    total = 0.0
    for j in range(panels):
        mid = a + (j + 0.5) * h
        half = 0.5 * h
        for xi, wi in zip(_GL8_X, _GL8_W):
            total += wi * f(mid + half * xi)
    return 0.5 * h * total


def rodrigues(axis, theta, reflect=False):
    """Build 3x3 orthogonal matrix in O(3)."""
    ux, uy, uz = axis
    norm = math.sqrt(ux * ux + uy * uy + uz * uz)
    ux, uy, uz = ux / norm, uy / norm, uz / norm
    c = math.cos(theta)
    s = math.sin(theta)
    omc = 1.0 - c
    R = [
        [c + ux * ux * omc, ux * uy * omc - uz * s, ux * uz * omc + uy * s],
        [uy * ux * omc + uz * s, c + uy * uy * omc, uy * uz * omc - ux * s],
        [uz * ux * omc - uy * s, uz * uy * omc + ux * s, c + uz * uz * omc],
    ]
    if reflect:
        R = [[-v for v in row] for row in R]
    return R


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def transpose(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]


def trace(A):
    return sum(A[i][i] for i in range(3))


def inertia_from_2d_integral(int_fn):
    """Compute 3x3 inertia tensor for a planar plate (z=0) given an integrator over D."""
    M = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            delta = 1.0 if i == j else 0.0

            def integrand(x, y, i=i, j=j, delta=delta):
                coords = (x, y, 0.0)
                r2 = x * x + y * y
                return r2 * delta - coords[i] * coords[j]

            M[i][j] = int_fn(integrand)
    return M


def main():
    # 1. Verify trace invariance M'_{ii} = M_{ii} for arbitrary rank-2 tensors under O(3)
    T = [
        [2.3, -1.7, 0.5],
        [0.9, -4.1, 3.2],
        [-2.0, 1.4, 5.6],
    ]  # Deliberately non-symmetric rank-2 tensor!
    max_trace_err = 0.0
    for axis, theta, refl in [
        ((1, 2, -3), 0.87, False),
        ((-4, 1, 2), 2.13, True),
        ((0, 0, 1), 1.23, False),
        ((1, 1, 1), math.pi / 3, True),
    ]:
        R = rodrigues(axis, theta, reflect=refl)
        # Check orthogonality R^T R = I
        RtR = matmul(transpose(R), R)
        for i in range(3):
            for j in range(3):
                assert abs(RtR[i][j] - (1.0 if i == j else 0.0)) < 1e-14
        # Transform T'_{ij} = R_{ip} R_{jq} T_{pq} = (R T R^T)_{ij}
        Tp = matmul(matmul(R, T), transpose(R))
        max_trace_err = max(max_trace_err, abs(trace(Tp) - trace(T)))

    # 2. Verify eigenvector e_3, eigenvalue M_perp, and perpendicular axis theorem M_perp = M_1 + M_2
    #    on an asymmetric off-centre plate with non-uniform density rho(x, y)
    #    D = rectangle [-0.7, 1.8] x [-1.2, 0.9] with asymmetric density rho(x, y) > 0
    def rho_nonunif(x, y):
        return 2.0 + 0.7 * x * x - 0.4 * x * y + 0.5 * math.cos(y) + 0.3 * math.sin(x)

    def int_rect(f):
        return quad1d(lambda x: quad1d(lambda y: f(x, y) * rho_nonunif(x, y), -1.2, 0.9, panels=8), -0.7, 1.8, panels=8)

    M_asym = inertia_from_2d_integral(int_rect)
    M_perp_int = int_rect(lambda x, y: x * x + y * y)

    # Check M e_3 = M_perp e_3
    e3 = [0.0, 0.0, 1.0]
    Me3 = [sum(M_asym[i][j] * e3[j] for j in range(3)) for i in range(3)]
    eig_z_err = max(abs(Me3[0]), abs(Me3[1]), abs(Me3[2] - M_perp_int))

    # Check symmetry and non-zero off-diagonal M_12 (so plate is genuinely rotated/asymmetric!)
    assert abs(M_asym[0][1] - M_asym[1][0]) < 1e-14
    assert abs(M_asym[0][1]) > 0.1, "Expected non-diagonal 2x2 block for asymmetric plate"

    # Eigenvalues M_1, M_2 of upper-left 2x2 symmetric block
    a11, a12, a22 = M_asym[0][0], M_asym[0][1], M_asym[1][1]
    tr2 = a11 + a22
    disc = math.sqrt((a11 - a22) ** 2 + 4.0 * a12 * a12)
    M1 = 0.5 * (tr2 - disc)
    M2 = 0.5 * (tr2 + disc)

    perp_thm_err = abs(M_perp_int - (M1 + M2))
    trace_2m_err = abs(trace(M_asym) - 2.0 * M_perp_int)

    # Check all three eigenvalues satisfy det(M - lambda I) = 0
    for lam in (M1, M2, M_perp_int):
        # Since M is block diagonal, det(M - lam I) = ((a11-lam)(a22-lam) - a12^2)*(M33 - lam)
        det_val = ((a11 - lam) * (a22 - lam) - a12 * a12) * (M_asym[2][2] - lam)
        assert abs(det_val) < 1e-12

    # 3. Verify circular disc of radius a and uniform density rho (polar AND Cartesian quadrature)
    max_disc_polar_err = 0.0
    max_disc_cart_err = 0.0
    for a, rho in [(1.0, 1.0), (2.3, 0.85), (0.7, 3.4)]:
        expected = [
            [0.25 * math.pi * rho * (a ** 4), 0.0, 0.0],
            [0.0, 0.25 * math.pi * rho * (a ** 4), 0.0],
            [0.0, 0.0, 0.50 * math.pi * rho * (a ** 4)],
        ]

        # Polar integration
        def int_disc_polar(f, a=a, rho=rho):
            return rho * quad1d(
                lambda r: r * quad1d(lambda th: f(r * math.cos(th), r * math.sin(th)), 0.0, 2.0 * math.pi, panels=12),
                0.0,
                a,
                panels=8,
            )

        M_polar = inertia_from_2d_integral(int_disc_polar)

        # Independent Cartesian integration: x in [-a, a], y in [-sqrt(a^2-x^2), sqrt(a^2-x^2)]
        # Substitute x = a*sin(u) to remove endpoint square-root singularity in quadrature!
        def int_disc_cart(f, a=a, rho=rho):
            def inner(u):
                x = a * math.sin(u)
                ym = a * math.cos(u)
                return (a * math.cos(u)) * quad1d(lambda y: f(x, y), -ym, ym, panels=10)

            return rho * quad1d(inner, -0.5 * math.pi, 0.5 * math.pi, panels=12)

        M_cart = inertia_from_2d_integral(int_disc_cart)

        for i in range(3):
            for j in range(3):
                # Also check index formula (pi*rho*a^4/4)*(delta_ij + delta_i3*delta_j3)
                idx_val = 0.25 * math.pi * rho * (a ** 4) * ((1.0 if i == j else 0.0) + (1.0 if i == 2 and j == 2 else 0.0))
                assert abs(idx_val - expected[i][j]) < 1e-14
                max_disc_polar_err = max(max_disc_polar_err, abs(M_polar[i][j] - expected[i][j]))
                max_disc_cart_err = max(max_disc_cart_err, abs(M_cart[i][j] - expected[i][j]))

    assert max_trace_err < 1e-13
    assert eig_z_err < 1e-13
    assert perp_thm_err < 1e-13
    assert trace_2m_err < 1e-13
    assert max_disc_polar_err < 1e-12
    assert max_disc_cart_err < 1e-12

    print("IB_2006_P1_q14A_check.py: ALL CHECKS PASSED")
    print(f"  max_trace_err      = {max_trace_err:.3e}")
    print(f"  eig_z_err          = {eig_z_err:.3e}")
    print(f"  perp_thm_err       = {perp_thm_err:.3e}")
    print(f"  trace_2m_err       = {trace_2m_err:.3e}")
    print(f"  max_disc_polar_err = {max_disc_polar_err:.3e}")
    print(f"  max_disc_cart_err  = {max_disc_cart_err:.3e}")


if __name__ == "__main__":
    main()
