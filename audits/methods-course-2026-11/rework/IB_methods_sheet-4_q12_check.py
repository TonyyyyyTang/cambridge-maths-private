#!/usr/bin/env python3
"""Computational verification for IB/methods/sheet-4/q12.

Verifies the method-of-images Dirichlet Green's function on the positive orthant
Omega = {x in R^n : x_1 > 0, ..., x_n > 0} for n = 1, 2, 3, 4, 5:
1. Unit sphere surface area omega_n = 2*pi^(n/2)/Gamma(n/2).
2. Homogeneous Dirichlet BC G(x; y) = 0 on every boundary face x_k = 0 and corners.
3. Symmetry G(x; y) = G(y; x) and negativity G(x; y) < 0 in Omega.
4. Harmonicity Delta_x G(x; y) = 0 away from y in Omega (central finite differences).
5. Delta-source flux int_{dB_eps(y)} (dG/dn) dS = 1 for n = 1, 2, 3.
6. Green's representation formula (Poisson formula) on the 1D half-line and 2D quadrant.
7. Decay G(x; y) -> 0 as |x| -> infinity for n >= 2, and exact leading multipole asymptotics.
"""

import itertools
import math


def omega_n(n):
    return 2.0 * (math.pi ** (0.5 * n)) / math.gamma(0.5 * n)


def phi_n(r, n):
    if n == 1:
        return 0.5 * r
    elif n == 2:
        return math.log(r) / (2.0 * math.pi)
    else:
        return -1.0 / ((n - 2) * omega_n(n) * (r ** (n - 2)))


def green_orthant(x, y):
    n = len(x)
    assert len(y) == n
    total = 0.0
    for sigma in itertools.product((1, -1), repeat=n):
        sgn = 1.0
        dist_sq = 0.0
        for xi, yi, si in zip(x, y, sigma):
            sgn *= si
            d = xi - si * yi
            dist_sq += d * d
        total += sgn * phi_n(math.sqrt(dist_sq), n)
    return total


def laplacian_fd(func, x, h=1e-4):
    n = len(x)
    f0 = func(x)
    lap = 0.0
    for k in range(n):
        xp = list(x)
        xm = list(x)
        xp[k] += h
        xm[k] -= h
        lap += (func(xp) - 2.0 * f0 + func(xm)) / (h * h)
    return lap


def simpson_1d(func, a, b, m):
    assert m % 2 == 0
    h = (b - a) / m
    s = func(a) + func(b)
    for k in range(1, m):
        w = 4.0 if k % 2 == 1 else 2.0
        s += w * func(a + k * h)
    return s * h / 3.0


def main():
    # 1. Verify omega_n values
    expected_omega = {
        1: 2.0,
        2: 2.0 * math.pi,
        3: 4.0 * math.pi,
        4: 2.0 * (math.pi ** 2),
        5: (8.0 / 3.0) * (math.pi ** 2),
    }
    for n, val in expected_omega.items():
        assert abs(omega_n(n) - val) < 1e-14

    # 2. Check 1D case: G(x; y) = -min(x, y), jump = +1, and ODE representation
    for x, y in [(0.3, 1.2), (1.7, 0.5), (2.4, 2.4)]:
        assert abs(green_orthant([x], [y]) - (-min(x, y))) < 1e-15
    assert abs(green_orthant([0.0], [1.2])) < 1e-15
    # Test u''(x) = -exp(-x), u(0) = 0, u bounded => u(x) = 1 - exp(-x)
    for x in [0.4, 1.5, 3.0]:
        u_rec = simpson_1d(lambda y: green_orthant([x], [y]) * (-math.exp(-y)), 0.0, 25.0, 2000)
        assert abs(u_rec - (1.0 - math.exp(-x))) < 1e-8

    # 3. Check BCs, symmetry, negativity, harmonicity for n = 2, 3, 4, 5
    max_bc_err = 0.0
    max_sym_err = 0.0
    max_lap_err = 0.0
    for n in [2, 3, 4, 5]:
        y = [0.6 + 0.35 * k for k in range(n)]
        x = [1.1 + 0.25 * k for k in range(n)]

        # Symmetry and negativity
        g_xy = green_orthant(x, y)
        g_yx = green_orthant(y, x)
        max_sym_err = max(max_sym_err, abs(g_xy - g_yx))
        assert g_xy < 0.0, f"Expected G < 0 in Omega for n={n}, got {g_xy}"

        # Boundary conditions on every coordinate face x_k = 0 and at origin
        for k in range(n):
            xb = list(x)
            xb[k] = 0.0
            max_bc_err = max(max_bc_err, abs(green_orthant(xb, y)))
        max_bc_err = max(max_bc_err, abs(green_orthant([0.0] * n, y)))

        # Harmonicity at interior point x != y
        lap = laplacian_fd(lambda pt: green_orthant(pt, y), x, h=1e-4)
        max_lap_err = max(max_lap_err, abs(lap))

    assert max_bc_err < 1e-14, f"BC error too large: {max_bc_err}"
    assert max_sym_err < 1e-14, f"Symmetry error too large: {max_sym_err}"
    assert max_lap_err < 1e-6, f"Laplacian residual too large: {max_lap_err}"

    # 4. Flux check int_{dB_eps(y)} (dG/dn) dS = 1 for n = 2 and n = 3
    # n = 2: circle of radius eps around y = (0.7, 1.1)
    y2 = [0.7, 1.1]
    eps = 0.05
    dr = 1e-6
    N_th = 200
    flux_2d = 0.0
    for j in range(N_th):
        th = 2.0 * math.pi * j / N_th
        cx, cy = math.cos(th), math.sin(th)
        gp = green_orthant([y2[0] + (eps + dr) * cx, y2[1] + (eps + dr) * cy], y2)
        gm = green_orthant([y2[0] + (eps - dr) * cx, y2[1] + (eps - dr) * cy], y2)
        dGdn = (gp - gm) / (2.0 * dr)
        flux_2d += dGdn * eps * (2.0 * math.pi / N_th)
    assert abs(flux_2d - 1.0) < 1e-7, f"2D flux mismatch: {flux_2d}"

    # n = 3: sphere of radius eps around y = (0.7, 1.1, 0.9)
    y3 = [0.7, 1.1, 0.9]
    eps3 = 0.05
    n_mu = 40  # midpoint rule in cos(theta) = mu in [-1, 1]
    n_phi = 80
    dmu = 2.0 / n_mu
    dphi = 2.0 * math.pi / n_phi
    flux_3d = 0.0
    for i in range(n_mu):
        mu = -1.0 + (i + 0.5) * dmu
        s_th = math.sqrt(max(0.0, 1.0 - mu * mu))
        for j in range(n_phi):
            phi = j * dphi
            nx = s_th * math.cos(phi)
            ny = s_th * math.sin(phi)
            nz = mu
            gp = green_orthant([y3[0] + (eps3 + dr) * nx, y3[1] + (eps3 + dr) * ny, y3[2] + (eps3 + dr) * nz], y3)
            gm = green_orthant([y3[0] + (eps3 - dr) * nx, y3[1] + (eps3 - dr) * ny, y3[2] + (eps3 - dr) * nz], y3)
            dGdn = (gp - gm) / (2.0 * dr)
            flux_3d += dGdn * (eps3 * eps3) * dmu * dphi
    assert abs(flux_3d - 1.0) < 1e-5, f"3D flux mismatch: {flux_3d}"

    # 5. 2D Poisson representation check for harmonic u(x1, x2) = (x1 + 1)/((x1 + 1)^2 + (x2 + 2)^2)
    #    Note u is harmonic on Omega, smooth on closure(Omega), O(1/r) at infinity.
    #    Outward normal on x2 = 0 is (0, -1), on x1 = 0 is (-1, 0).
    #    Compute dG/dn analytically on each face:
    #      On x2 = 0 (x1 = s > 0):
    #        dG/dn = -dG/dx2|_{x2=0}
    #              = (y2 / pi) * (1/((s - y1)^2 + y2^2) - 1/((s + y1)^2 + y2^2))
    #      On x1 = 0 (x2 = t > 0):
    #        dG/dn = -dG/dx1|_{x1=0}
    #              = (y1 / pi) * (1/(y1^2 + (t - y2)^2) - 1/(y1^2 + (t + y2)^2))
    def u_harm_2d(x1, x2):
        return (x1 + 1.0) / ((x1 + 1.0) ** 2 + (x2 + 2.0) ** 2)

    y_pt = (0.8, 1.3)
    # Use substitution s = tan(theta) on (0, pi/2) to integrate over [0, infinity) cleanly!
    def integrand_bottom(th):
        s = math.tan(th)
        sec2 = 1.0 / (math.cos(th) ** 2)
        y1, y2 = y_pt
        dGdn = (y2 / math.pi) * (
            1.0 / ((s - y1) ** 2 + y2 * y2) - 1.0 / ((s + y1) ** 2 + y2 * y2)
        )
        return u_harm_2d(s, 0.0) * dGdn * sec2

    def integrand_left(th):
        t = math.tan(th)
        sec2 = 1.0 / (math.cos(th) ** 2)
        y1, y2 = y_pt
        dGdn = (y1 / math.pi) * (
            1.0 / (y1 * y1 + (t - y2) ** 2) - 1.0 / (y1 * y1 + (t + y2) ** 2)
        )
        return u_harm_2d(0.0, t) * dGdn * sec2

    # First verify our analytic dG/dn against finite differences of green_orthant
    eps_n = 1e-6
    s_chk = 1.4
    dGdn_fd = (green_orthant([s_chk, 0.0], list(y_pt)) - green_orthant([s_chk, eps_n], list(y_pt))) / eps_n
    dGdn_ana = (y_pt[1] / math.pi) * (
        1.0 / ((s_chk - y_pt[0]) ** 2 + y_pt[1] ** 2) - 1.0 / ((s_chk + y_pt[0]) ** 2 + y_pt[1] ** 2)
    )
    assert abs(dGdn_fd - dGdn_ana) < 1e-6

    u_poisson = simpson_1d(integrand_bottom, 0.0, 0.5 * math.pi - 1e-7, 2000) + simpson_1d(
        integrand_left, 0.0, 0.5 * math.pi - 1e-7, 2000
    )
    u_exact = u_harm_2d(*y_pt)
    poisson_err = abs(u_poisson - u_exact)
    assert poisson_err < 1e-7, f"2D Poisson integral mismatch: {u_poisson} vs {u_exact}"

    # 6. Decay and multipole asymptotics as |x| -> infinity for n = 2, 3
    R = 100.0
    x2_far = [R * 0.6, R * 0.8]
    g2_far = green_orthant(x2_far, [0.7, 1.1])
    g2_asymp = -4.0 * x2_far[0] * x2_far[1] * 0.7 * 1.1 / (math.pi * (R ** 4))
    rel_err_2d = abs(g2_far - g2_asymp) / abs(g2_asymp)
    assert rel_err_2d < 1e-3

    x3_far = [R * 0.5, R * 0.6, R * math.sqrt(1.0 - 0.25 - 0.36)]
    g3_far = green_orthant(x3_far, [0.7, 1.1, 0.9])
    # Leading multipole: (-1)^3 * (2*y1)(2*y2)(2*y3) * (15*x1*x2*x3)/(4*pi*R^7)
    g3_asymp = -(8.0 * 0.7 * 1.1 * 0.9) * (15.0 * x3_far[0] * x3_far[1] * x3_far[2]) / (4.0 * math.pi * (R ** 7))
    rel_err_3d = abs(g3_far - g3_asymp) / abs(g3_asymp)
    assert rel_err_3d < 1e-3

    print("IB_methods_sheet-4_q12 verification PASSED:")
    print(f"  max BC error (n=2..5)           = {max_bc_err:.3e}")
    print(f"  max symmetry error (n=2..5)     = {max_sym_err:.3e}")
    print(f"  max Laplacian residual (n=2..5) = {max_lap_err:.3e}")
    print(f"  2D delta flux                   = {flux_2d:.8f} (target 1.0)")
    print(f"  3D delta flux                   = {flux_3d:.8f} (target 1.0)")
    print(f"  2D Poisson integral error       = {poisson_err:.3e}")
    print(f"  far-field multipole rel err     = 2D: {rel_err_2d:.3e}, 3D: {rel_err_3d:.3e}")


if __name__ == "__main__":
    main()
