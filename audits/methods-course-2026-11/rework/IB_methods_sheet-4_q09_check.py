#!/usr/bin/env python3
"""Computational verification for IB/methods/sheet-4/q09.

Verifies the Dirichlet Green's function for the Laplacian in R^3 for both:
  (i)  the interior of a sphere of radius a (B_a = {|x| < a}), including x0 -> 0,
  (ii) the exterior of a sphere of radius a (E_a = {|x| > a}), including decay at infinity.
Checks:
1. Image location (|x0*| > a for interior, |x0*| < a for exterior).
2. Homogeneous Dirichlet BC G(x; x0) = 0 on |x| = a.
3. Harmonicity away from x0 via central finite differences.
4. Delta-function flux int_{dB_eps(x0)} (dG/dn) dS = 1.
5. Symmetry G(x; x0) = G(x0; x) and regularised form valid at x0 = 0.
6. Poisson integral formula on both interior and exterior domains for exact harmonic test functions.
7. Far-field decay G_ext(x; x0) ~ -(1 - a/|x0|)/(4*pi*|x|) as |x| -> infinity.
"""

import math


def vnorm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def vsub(u, v):
    return (u[0] - v[0], u[1] - v[1], u[2] - v[2])


def vdot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def invert_sphere(x0, a):
    r0_sq = vdot(x0, x0)
    scale = (a * a) / r0_sq
    return (scale * x0[0], scale * x0[1], scale * x0[2])


def green_sphere(x, x0, a):
    """Formula valid for x0 != 0 (both interior |x0|<a and exterior |x0|>a)."""
    r0 = vnorm(x0)
    x0_star = invert_sphere(x0, a)
    d1 = vnorm(vsub(x, x0))
    d2 = vnorm(vsub(x, x0_star))
    return -1.0 / (4.0 * math.pi * d1) + (a / r0) / (4.0 * math.pi * d2)


def green_sphere_reg(x, x0, a):
    """Regularised form valid for ALL x0 (including x0 = 0)."""
    d1 = vnorm(vsub(x, x0))
    denom = math.sqrt(vdot(x, x) * vdot(x0, x0) - 2.0 * (a * a) * vdot(x, x0) + a ** 4)
    return -1.0 / (4.0 * math.pi * d1) + a / (4.0 * math.pi * denom)


def laplacian_fd(func, x, h=1e-4):
    f0 = func(x)
    lap = 0.0
    for k in range(3):
        xp = list(x)
        xm = list(x)
        xp[k] += h
        xm[k] -= h
        lap += (func(tuple(xp)) - 2.0 * f0 + func(tuple(xm))) / (h * h)
    return lap


def flux_around_source(x0, a, eps=0.04):
    dr = 1e-6
    n_mu = 40
    n_phi = 80
    dmu = 2.0 / n_mu
    dphi = 2.0 * math.pi / n_phi
    flux = 0.0
    for i in range(n_mu):
        mu = -1.0 + (i + 0.5) * dmu
        s_th = math.sqrt(max(0.0, 1.0 - mu * mu))
        for j in range(n_phi):
            phi = j * dphi
            nvec = (s_th * math.cos(phi), s_th * math.sin(phi), mu)
            xp = (x0[0] + (eps + dr) * nvec[0], x0[1] + (eps + dr) * nvec[1], x0[2] + (eps + dr) * nvec[2])
            xm = (x0[0] + (eps - dr) * nvec[0], x0[1] + (eps - dr) * nvec[1], x0[2] + (eps - dr) * nvec[2])
            dGdn = (green_sphere(xp, x0, a) - green_sphere(xm, x0, a)) / (2.0 * dr)
            flux += dGdn * (eps * eps) * dmu * dphi
    return flux


def integrate_sphere(func, a, n_mu=80, n_phi=160):
    """Integrate over sphere |x|=a using composite Simpson's rule in mu=cos(theta)
    and periodic trapezoidal rule in phi."""
    assert n_mu % 2 == 0
    dmu = 2.0 / n_mu
    dphi = 2.0 * math.pi / n_phi
    total = 0.0
    for i in range(n_mu + 1):
        mu = -1.0 + i * dmu
        s_th = math.sqrt(max(0.0, 1.0 - mu * mu))
        if i == 0 or i == n_mu:
            w_mu = 1.0
        elif i % 2 == 1:
            w_mu = 4.0
        else:
            w_mu = 2.0
        phi_sum = 0.0
        for j in range(n_phi):
            phi = j * dphi
            pt = (a * s_th * math.cos(phi), a * s_th * math.sin(phi), a * mu)
            phi_sum += func(pt)
        total += w_mu * phi_sum * dphi
    return total * (a * a) * (dmu / 3.0)


def main():
    a = 2.5

    # ------------------------------------------------------------------
    # Part (i): Interior of sphere |x| < a
    # ------------------------------------------------------------------
    x0_int = (0.7, -0.5, 1.1)  # |x0| ~ 1.396 < 2.5
    assert vnorm(x0_int) < a
    x0_int_star = invert_sphere(x0_int, a)
    assert vnorm(x0_int_star) > a

    # Check BC on |x| = a
    max_bc_int = 0.0
    for mu in [-0.9, -0.2, 0.4, 0.95]:
        s_th = math.sqrt(1.0 - mu * mu)
        for phi in [0.0, 1.1, 2.7, 4.5]:
            xb = (a * s_th * math.cos(phi), a * s_th * math.sin(phi), a * mu)
            max_bc_int = max(max_bc_int, abs(green_sphere(xb, x0_int, a)))
    assert max_bc_int < 1e-14, f"Interior BC error too large: {max_bc_int}"

    # Check harmonicity, symmetry, negativity, equivalence with regularised form
    x_int = (-0.6, 0.9, 0.4)
    lap_int = laplacian_fd(lambda p: green_sphere(p, x0_int, a), x_int)
    assert abs(lap_int) < 1e-6
    assert abs(green_sphere(x_int, x0_int, a) - green_sphere(x0_int, x_int, a)) < 1e-15
    assert abs(green_sphere(x_int, x0_int, a) - green_sphere_reg(x_int, x0_int, a)) < 1e-15
    assert green_sphere(x_int, x0_int, a) < 0.0

    # Check limit x0 -> 0 matches -1/(4*pi*|x|) + 1/(4*pi*a)
    x0_tiny = (1e-8, -2e-8, 1.5e-8)
    g_zero_exact = -1.0 / (4.0 * math.pi * vnorm(x_int)) + 1.0 / (4.0 * math.pi * a)
    assert abs(green_sphere_reg(x_int, (0.0, 0.0, 0.0), a) - g_zero_exact) < 1e-15
    assert abs(green_sphere(x_int, x0_tiny, a) - g_zero_exact) < 1e-8

    # Check delta flux = 1
    flux_int = flux_around_source(x0_int, a)
    assert abs(flux_int - 1.0) < 1e-5, f"Interior flux error: {flux_int}"

    # Check interior Poisson formula: dG/dn = (a^2 - |x0|^2) / (4*pi*a*|x - x0|^3)
    def dGdn_int_fd(xb, x0):
        eps_n = 1e-6
        r = vnorm(xb)
        xin = tuple(c * (r - eps_n) / r for c in xb)
        return (green_sphere(xb, x0, a) - green_sphere(xin, x0, a)) / eps_n

    def dGdn_int_ana(xb, x0):
        return (a * a - vdot(x0, x0)) / (4.0 * math.pi * a * (vnorm(vsub(xb, x0)) ** 3))

    xb_chk = (a * 0.6, a * 0.8, 0.0)
    assert abs(dGdn_int_fd(xb_chk, x0_int) - dGdn_int_ana(xb_chk, x0_int)) < 1e-6

    # Test harmonic u1(x) = 1 and u2(x) = x^2 - y^2 + 3*x*z - 2*y + 5
    def u_int_harm(p):
        return p[0] ** 2 - p[1] ** 2 + 3.0 * p[0] * p[2] - 2.0 * p[1] + 5.0

    int_const = integrate_sphere(lambda xb: dGdn_int_ana(xb, x0_int), a)
    assert abs(int_const - 1.0) < 1e-5

    int_harm = integrate_sphere(lambda xb: u_int_harm(xb) * dGdn_int_ana(xb, x0_int), a)
    int_harm_err = abs(int_harm - u_int_harm(x0_int))
    assert int_harm_err < 1e-4, f"Interior Poisson mismatch: {int_harm} vs {u_int_harm(x0_int)}"

    # ------------------------------------------------------------------
    # Part (ii): Exterior of sphere |x| > a
    # ------------------------------------------------------------------
    x0_ext = (2.1, -1.8, 3.2)  # |x0| ~ 4.230 > 2.5
    assert vnorm(x0_ext) > a
    x0_ext_star = invert_sphere(x0_ext, a)
    assert vnorm(x0_ext_star) < a

    # Check BC on |x| = a
    max_bc_ext = 0.0
    for mu in [-0.9, -0.2, 0.4, 0.95]:
        s_th = math.sqrt(1.0 - mu * mu)
        for phi in [0.0, 1.1, 2.7, 4.5]:
            xb = (a * s_th * math.cos(phi), a * s_th * math.sin(phi), a * mu)
            max_bc_ext = max(max_bc_ext, abs(green_sphere(xb, x0_ext, a)))
    assert max_bc_ext < 1e-14, f"Exterior BC error too large: {max_bc_ext}"

    # Check harmonicity, symmetry, negativity
    x_ext = (-3.0, 1.5, 2.2)
    lap_ext = laplacian_fd(lambda p: green_sphere(p, x0_ext, a), x_ext)
    assert abs(lap_ext) < 1e-6
    assert abs(green_sphere(x_ext, x0_ext, a) - green_sphere(x0_ext, x_ext, a)) < 1e-15
    assert green_sphere(x_ext, x0_ext, a) < 0.0

    # Check delta flux = 1
    flux_ext = flux_around_source(x0_ext, a)
    assert abs(flux_ext - 1.0) < 1e-5, f"Exterior flux error: {flux_ext}"

    # Check far-field decay G_ext(x; x0) ~ -(1 - a/|x0|) / (4*pi*|x|)
    R_far = 1000.0
    x_far = (R_far * 0.3, R_far * 0.4, R_far * math.sqrt(1.0 - 0.09 - 0.16))
    g_far = green_sphere(x_far, x0_ext, a)
    g_far_asymp = -(1.0 - a / vnorm(x0_ext)) / (4.0 * math.pi * R_far)
    far_rel_err = abs(g_far - g_far_asymp) / abs(g_far_asymp)
    assert far_rel_err < 1e-2

    # Check exterior Poisson formula: outward normal from E_a is -x_hat!
    # So dG_ext/dn = -dG_ext/dr|_{r=a} = (|x0|^2 - a^2) / (4*pi*a*|x - x0|^3) > 0.
    def dGdn_ext_fd(xb, x0):
        eps_n = 1e-6
        r = vnorm(xb)
        xout = tuple(c * (r + eps_n) / r for c in xb)
        return (green_sphere(xb, x0, a) - green_sphere(xout, x0, a)) / eps_n

    def dGdn_ext_ana(xb, x0):
        return (vdot(x0, x0) - a * a) / (4.0 * math.pi * a * (vnorm(vsub(xb, x0)) ** 3))

    assert abs(dGdn_ext_fd(xb_chk, x0_ext) - dGdn_ext_ana(xb_chk, x0_ext)) < 1e-6

    # 1) u(x) = a/|x| (= 1 on |x|=a, decays to 0 at infinity)
    ext_const = integrate_sphere(lambda xb: dGdn_ext_ana(xb, x0_ext), a)
    ext_const_err = abs(ext_const - a / vnorm(x0_ext))
    assert ext_const_err < 1e-5

    # 2) u(x) = 1/|x - p| with |p| < a (harmonic on E_a, decays at infinity)
    p_in = (0.5, 0.8, -0.6)
    assert vnorm(p_in) < a

    def u_ext_harm(pt):
        return 1.0 / vnorm(vsub(pt, p_in))

    ext_harm = integrate_sphere(lambda xb: u_ext_harm(xb) * dGdn_ext_ana(xb, x0_ext), a)
    ext_harm_err = abs(ext_harm - u_ext_harm(x0_ext))
    assert ext_harm_err < 1e-5, f"Exterior Poisson mismatch: {ext_harm} vs {u_ext_harm(x0_ext)}"

    print("IB_methods_sheet-4_q09 verification PASSED:")
    print(f"  Interior max BC error          = {max_bc_int:.3e}")
    print(f"  Interior Laplacian residual    = {abs(lap_int):.3e}")
    print(f"  Interior delta flux            = {flux_int:.8f} (target 1.0)")
    print(f"  Interior Poisson err (poly)    = {int_harm_err:.3e}")
    print(f"  Exterior max BC error          = {max_bc_ext:.3e}")
    print(f"  Exterior Laplacian residual    = {abs(lap_ext):.3e}")
    print(f"  Exterior delta flux            = {flux_ext:.8f} (target 1.0)")
    print(f"  Exterior far-field rel err     = {far_rel_err:.3e}")
    print(f"  Exterior Poisson err (a/|x|)   = {ext_const_err:.3e}")
    print(f"  Exterior Poisson err (1/|x-p|) = {ext_harm_err:.3e}")


if __name__ == "__main__":
    main()
