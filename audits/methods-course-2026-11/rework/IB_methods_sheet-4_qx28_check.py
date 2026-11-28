#!/usr/bin/env python3
"""Independent verification for IB/methods/sheet-4/qx28 (Symmetry of 3D Dirichlet Green's function)."""

import math
import random


def v_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def v_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def v_scale(c, a):
    return (c * a[0], c * a[1], c * a[2])


def v_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def v_norm(a):
    return math.sqrt(v_dot(a, a))


# --- Exact 3D Dirichlet Green's functions (convention: laplacian G = +delta) ---

def G_ball(r, r0, a):
    """Dirichlet Green's function for the ball |r| < a via Kelvin inversion."""
    d = v_norm(v_sub(r, r0))
    val = -1.0 / (4.0 * math.pi * d)
    r0_norm = v_norm(r0)
    if r0_norm > 0.0:
        r0_star = v_scale((a * a) / (r0_norm * r0_norm), r0)
        d_star = v_norm(v_sub(r, r0_star))
        val += (a / r0_norm) / (4.0 * math.pi * d_star)
    else:
        val += 1.0 / (4.0 * math.pi * a)
    return val


def grad_G_ball(r, r0, a):
    """Exact gradient with respect to r of G_ball(r; r0)."""
    diff = v_sub(r, r0)
    d = v_norm(diff)
    g = v_scale(1.0 / (4.0 * math.pi * d**3), diff)
    r0_norm = v_norm(r0)
    if r0_norm > 0.0:
        r0_star = v_scale((a * a) / (r0_norm * r0_norm), r0)
        diff_s = v_sub(r, r0_star)
        ds = v_norm(diff_s)
        g_img = v_scale(-(a / r0_norm) / (4.0 * math.pi * ds**3), diff_s)
        g = v_add(g, g_img)
    return g


def G_halfspace(r, r0):
    """Dirichlet Green's function for the half-space z > 0 via method of images."""
    r0_img = (r0[0], r0[1], -r0[2])
    d = v_norm(v_sub(r, r0))
    d_img = v_norm(v_sub(r, r0_img))
    return -1.0 / (4.0 * math.pi * d) + 1.0 / (4.0 * math.pi * d_img)


def grad_G_halfspace(r, r0):
    r0_img = (r0[0], r0[1], -r0[2])
    diff = v_sub(r, r0)
    diff_img = v_sub(r, r0_img)
    d = v_norm(diff)
    d_img = v_norm(diff_img)
    g1 = v_scale(1.0 / (4.0 * math.pi * d**3), diff)
    g2 = v_scale(-1.0 / (4.0 * math.pi * d_img**3), diff_img)
    return v_add(g1, g2)


# --- Gauss-Legendre nodes/weights on [-1, 1] (dependency-free Newton method) ---

def gauss_legendre(n):
    nodes = []
    weights = []
    for k in range(1, n + 1):
        # Initial guess
        x = math.cos(math.pi * (k - 0.25) / (n + 0.5))
        for _ in range(30):
            p0, p1 = 1.0, x
            for j in range(2, n + 1):
                p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-15:
                break
        # Recompute derivative at converged root
        p0, p1 = 1.0, x
        for j in range(2, n + 1):
            p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
        dp = n * (x * p1 - p0) / (x * x - 1.0)
        w = 2.0 / ((1.0 - x * x) * dp * dp)
        nodes.append(x)
        weights.append(w)
    return nodes, weights


def sphere_green_integral(center, radius, outward_from_ball, u_fn, gu_fn, v_fn, gv_fn, n_mu=20, n_phi=40):
    """
    Compute int_S (u dv/dn - v du/dn) dS on sphere |r - center| = radius.
    If outward_from_ball is False, normal n points INWARD toward center (outward from excised domain D_eps).
    """
    mus, w_mus = gauss_legendre(n_mu)
    dphi = 2.0 * math.pi / n_phi
    sgn = 1.0 if outward_from_ball else -1.0
    total = 0.0
    for mu, w_mu in zip(mus, w_mus):
        sin_t = math.sqrt(max(0.0, 1.0 - mu * mu))
        for j in range(n_phi):
            phi = j * dphi
            er = (sin_t * math.cos(phi), sin_t * math.sin(phi), mu)
            r = v_add(center, v_scale(radius, er))
            n_vec = v_scale(sgn, er)
            u = u_fn(r)
            v = v_fn(r)
            du_dn = v_dot(gu_fn(r), n_vec)
            dv_dn = v_dot(gv_fn(r), n_vec)
            total += w_mu * dphi * (u * dv_dn - v * du_dn)
    return (radius * radius) * total


def hemisphere_green_integral(R, u_fn, gu_fn, v_fn, gv_fn, n_mu=24, n_phi=48):
    """Integrate (u dv/dn - v du/dn) over upper hemisphere |r| = R, z > 0 with outward normal +e_r."""
    mus, w_mus = gauss_legendre(n_mu)
    # Map [-1, 1] -> [0, 1] for mu = cos(theta) in [0, 1]
    dphi = 2.0 * math.pi / n_phi
    total = 0.0
    for xi, w_xi in zip(mus, w_mus):
        mu = 0.5 * (xi + 1.0)
        w_mu = 0.5 * w_xi
        sin_t = math.sqrt(max(0.0, 1.0 - mu * mu))
        for j in range(n_phi):
            phi = j * dphi
            er = (sin_t * math.cos(phi), sin_t * math.sin(phi), mu)
            r = v_scale(R, er)
            u = u_fn(r)
            v = v_fn(r)
            du_dn = v_dot(gu_fn(r), er)
            dv_dn = v_dot(gv_fn(r), er)
            total += w_mu * dphi * (u * dv_dn - v * du_dn)
    return (R * R) * total


def fd_laplacian(fn, r, h=1e-4):
    val0 = fn(r)
    s = 0.0
    for k in range(3):
        rp = list(r)
        rm = list(r)
        rp[k] += h
        rm[k] -= h
        s += (fn(tuple(rp)) - 2.0 * val0 + fn(tuple(rm))) / (h * h)
    return s


def main():
    rng = random.Random(20261127)
    a = 2.5

    # 1. Symmetry, Dirichlet BC, and harmonicity on the ball |r| < a
    max_sym_ball = 0.0
    max_bc_ball = 0.0
    for _ in range(25):
        while True:
            r1 = tuple(rng.uniform(-1.5, 1.5) for _ in range(3))
            r2 = tuple(rng.uniform(-1.5, 1.5) for _ in range(3))
            if v_norm(r1) < a - 0.3 and v_norm(r2) < a - 0.3 and v_norm(v_sub(r1, r2)) > 0.25:
                break
        g12 = G_ball(r1, r2, a)
        g21 = G_ball(r2, r1, a)
        max_sym_ball = max(max_sym_ball, abs(g12 - g21))

        # Check BC on sphere |r| = a
        dir_vec = tuple(rng.gauss(0.0, 1.0) for _ in range(3))
        rb = v_scale(a / v_norm(dir_vec), dir_vec)
        max_bc_ball = max(max_bc_ball, abs(G_ball(rb, r1, a)))

    assert max_sym_ball < 1e-14
    assert max_bc_ball < 1e-14

    r1 = (0.7, -0.4, 0.5)
    r2 = (-0.6, 0.9, -0.3)
    lap_ball = abs(fd_laplacian(lambda r: G_ball(r, r1, a), r2))
    assert lap_ball < 1e-7
    print(
        f"OK: Ball (a={a}) — max symmetry err={max_sym_ball:.2e}, "
        f"max BC err={max_bc_ball:.2e}, |laplacian G|={lap_ball:.2e}."
    )

    # 2. Direct numerical Green's second identity on excised ball D_eps = B_a \ (B_eps(r1) U B_eps(r2))
    u_b = lambda r: G_ball(r, r1, a)
    gu_b = lambda r: grad_G_ball(r, r1, a)
    v_b = lambda r: G_ball(r, r2, a)
    gv_b = lambda r: grad_G_ball(r, r2, a)

    G12 = G_ball(r1, r2, a)
    G21 = G_ball(r2, r1, a)
    assert abs(G12 - G21) < 1e-15

    I_outer_ball = sphere_green_integral((0.0, 0.0, 0.0), a, True, u_b, gu_b, v_b, gv_b, n_mu=20, n_phi=40)
    assert abs(I_outer_ball) < 1e-14

    def sphere_abs_u_dvdn(center, radius, u_fn, gv_fn, n_mu=20, n_phi=40):
        mus, w_mus = gauss_legendre(n_mu)
        dphi = 2.0 * math.pi / n_phi
        total = 0.0
        for mu, w_mu in zip(mus, w_mus):
            sin_t = math.sqrt(max(0.0, 1.0 - mu * mu))
            for j in range(n_phi):
                phi = j * dphi
                er = (sin_t * math.cos(phi), sin_t * math.sin(phi), mu)
                r = v_add(center, v_scale(radius, er))
                n_vec = v_scale(-1.0, er)
                total += w_mu * dphi * abs(u_fn(r) * v_dot(gv_fn(r), n_vec))
        return (radius * radius) * total

    prev_l1 = None
    for eps in (0.1, 0.05, 0.02, 0.01, 0.005):
        I1 = sphere_green_integral(r1, eps, False, u_b, gu_b, v_b, gv_b, n_mu=20, n_phi=40)
        I2 = sphere_green_integral(r2, eps, False, u_b, gu_b, v_b, gv_b, n_mu=20, n_phi=40)
        # Green's second identity holds EXACTLY at every finite eps, and by harmonicity
        # I1(eps) == G(r1; r2) and I2(eps) == -G(r2; r1) identically for every finite eps!
        assert abs(I1 + I2 + I_outer_ball) < 1e-11, f"Green's identity failed at eps={eps}: {I1+I2}"
        assert abs(I1 - G12) < 1e-12, f"I1 mismatch at eps={eps}: {I1} vs {G12}"
        assert abs(I2 - (-G21)) < 1e-12, f"I2 mismatch at eps={eps}: {I2} vs {-G21}"
        l1_term = sphere_abs_u_dvdn(r1, eps, u_b, gv_b, n_mu=20, n_phi=40)
        if prev_l1 is not None:
            assert l1_term < prev_l1
        prev_l1 = l1_term
        if eps == 0.005:
            print(
                f"OK: Excised-sphere Green's 2nd identity on ball (eps=0.1..0.005): "
                f"I1={I1:.10f} (G(r1;r2)={G12:.10f}), I2={I2:.10f} (-G(r2;r1)={-G21:.10f}), "
                f"|I1+I2|={abs(I1+I2):.2e}, int|u dv/dn|={l1_term:.2e}."
            )

    # 3. Unbounded domain: Half-space z > 0
    r1_h = (0.3, -0.8, 1.2)
    r2_h = (-0.9, 0.5, 0.7)
    u_h = lambda r: G_halfspace(r, r1_h)
    gu_h = lambda r: grad_G_halfspace(r, r1_h)
    v_h = lambda r: G_halfspace(r, r2_h)
    gv_h = lambda r: grad_G_halfspace(r, r2_h)

    assert abs(G_halfspace(r1_h, r2_h) - G_halfspace(r2_h, r1_h)) < 1e-15
    assert abs(G_halfspace((1.7, -2.3, 0.0), r1_h)) < 1e-15
    assert abs(fd_laplacian(u_h, r2_h)) < 1e-7

    I1_h = sphere_green_integral(r1_h, 0.005, False, u_h, gu_h, v_h, gv_h, n_mu=20, n_phi=40)
    I2_h = sphere_green_integral(r2_h, 0.005, False, u_h, gu_h, v_h, gv_h, n_mu=20, n_phi=40)
    assert abs(I1_h - G_halfspace(r1_h, r2_h)) < 1e-4
    assert abs(I2_h + G_halfspace(r2_h, r1_h)) < 1e-4
    assert abs(I1_h + I2_h) < 1e-11

    # Outer hemisphere integral decay as R -> infinity
    outer_vals = []
    for R in (10.0, 50.0, 200.0):
        I_R = hemisphere_green_integral(R, u_h, gu_h, v_h, gv_h, n_mu=24, n_phi=48)
        outer_vals.append((R, abs(I_R)))
    assert outer_vals[0][1] > outer_vals[1][1] > outer_vals[2][1]
    assert outer_vals[-1][1] < 1e-10
    print(
        f"OK: Half-space z>0 verified — |I1+I2|={abs(I1_h+I2_h):.2e}, "
        f"outer hemisphere |I_R|: R=10 -> {outer_vals[0][1]:.2e}, "
        f"R=50 -> {outer_vals[1][1]:.2e}, R=200 -> {outer_vals[2][1]:.2e}."
    )


if __name__ == "__main__":
    main()
