"""Computational verification for IB/2006/P2/q15G (Green's function ODE)."""
import math


def gl_quad(func, a, b, n_panels=64):
    """Composite 4-point Gauss-Legendre quadrature on [a, b]."""
    nodes = (
        -0.8611363115940526,
        -0.33998104358485626,
        0.33998104358485626,
        0.8611363115940526,
    )
    weights = (
        0.34785484513745385,
        0.6521451548625461,
        0.6521451548625461,
        0.34785484513745385,
    )
    h = (b - a) / n_panels
    total = 0.0
    for i in range(n_panels):
        mid = a + (i + 0.5) * h
        half = 0.5 * h
        for nd, wt in zip(nodes, weights):
            total += wt * func(mid + half * nd)
    return total * (0.5 * h)


def L_fd(func, x, h=1e-4):
    """Finite-difference evaluation of L[y] = y'' + (x+1)/(x+2) y' - 1/(x+2) y."""
    y0 = func(x)
    yp = (func(x + h) - func(x - h)) / (2.0 * h)
    ypp = (func(x + h) - 2.0 * y0 + func(x - h)) / (h * h)
    return ypp + ((x + 1.0) / (x + 2.0)) * yp - (1.0 / (x + 2.0)) * y0


def G(x, xi):
    if x < xi:
        return -(x + 1.0) / (xi + 2.0)
    return -((xi + 1.0) / (xi + 2.0)) * math.exp(xi - x)


def y_exact(x):
    return (0.5 * x * x + 2.0 * x + 1.0) * math.exp(-x)


def main():
    # 1. Verify homogeneous solutions y1 = e^{-x} and y2 = x + 1
    max_hom1 = 0.0
    max_hom2 = 0.0
    for i in range(1, 50):
        x = i * 0.2
        y1 = math.exp(-x)
        y1p = -math.exp(-x)
        y1pp = math.exp(-x)
        res1 = (x + 2.0) * y1pp + (x + 1.0) * y1p - y1
        max_hom1 = max(max_hom1, abs(res1))

        y2 = x + 1.0
        y2p = 1.0
        y2pp = 0.0
        res2 = (x + 2.0) * y2pp + (x + 1.0) * y2p - y2
        max_hom2 = max(max_hom2, abs(res2))
    assert max_hom1 < 1e-14, f"y1 residual too large: {max_hom1}"
    assert max_hom2 < 1e-14, f"y2 residual too large: {max_hom2}"

    # 2. Verify Green's function properties for several xi > 0
    max_ode_err = 0.0
    max_cont_err = 0.0
    max_jump_err = 0.0
    max_bc0_err = 0.0
    max_sym_err = 0.0

    xis = [0.25, 0.7, 1.5, 3.0, 6.2]
    for xi in xis:
        # L_x[G] = 0 for x < xi and x > xi
        for x in [0.05 * xi, 0.5 * xi, 0.9 * xi, xi + 0.2, xi + 1.5, xi + 4.0]:
            err = abs(L_fd(lambda t: G(t, xi), x, h=1e-4))
            max_ode_err = max(max_ode_err, err)

        # Boundary condition at x = 0: G_x(0, xi) = G(0, xi)
        g0 = G(0.0, xi)
        # Forward 2nd-order derivative at 0
        h = 1e-5
        gp0 = (-3.0 * G(0.0, xi) + 4.0 * G(h, xi) - G(2.0 * h, xi)) / (2.0 * h)
        max_bc0_err = max(max_bc0_err, abs(gp0 - g0))

        # Decay as x -> infinity
        assert abs(G(xi + 50.0, xi)) < 1e-20

        # Continuity at x = xi
        eps = 1e-12
        max_cont_err = max(max_cont_err, abs(G(xi + eps, xi) - G(xi - eps, xi)))

        # Derivative jump [G_x]_{xi-}^{xi+} = 1
        h = 1e-6
        gp_plus = (-3.0 * G(xi, xi) + 4.0 * G(xi + h, xi) - G(xi + 2.0 * h, xi)) / (2.0 * h)
        gp_minus = (3.0 * G(xi, xi) - 4.0 * G(xi - h, xi) + G(xi - 2.0 * h, xi)) / (2.0 * h)
        max_jump_err = max(max_jump_err, abs((gp_plus - gp_minus) - 1.0))

    # Self-adjoint symmetry: (xi + 2) e^{-xi} G(x, xi) is symmetric in (x, xi)
    for x in [0.1, 0.8, 2.3, 5.0]:
        for xi in [0.3, 1.2, 3.7]:
            if abs(x - xi) < 1e-9:
                continue
            gt1 = (xi + 2.0) * math.exp(-xi) * G(x, xi)
            gt2 = (x + 2.0) * math.exp(-x) * G(xi, x)
            max_sym_err = max(max_sym_err, abs(gt1 - gt2))

    assert max_ode_err < 1e-7, f"G ODE err: {max_ode_err}"
    assert max_bc0_err < 1e-10, f"G BC0 err: {max_bc0_err}"
    assert max_cont_err < 1e-10, f"G continuity err: {max_cont_err}"
    assert max_jump_err < 1e-8, f"G jump err: {max_jump_err}"
    assert max_sym_err < 1e-14, f"Gtilde symmetry err: {max_sym_err}"

    # 3. Verify particular solution y(x) = (x^2/2 + 2x + 1) e^{-x}
    #    both by Green's convolution and by direct ODE + BC substitution
    max_conv_err = 0.0
    max_sol_ode_err = 0.0
    f_prob = lambda xi: -(xi + 2.0) * math.exp(-xi)

    for x in [0.0, 0.25, 0.8, 1.5, 3.0, 5.5, 10.0]:
        part1 = gl_quad(lambda xi: G(x, xi) * f_prob(xi), 0.0, x, 32) if x > 0 else 0.0
        # Map [x, inf) via xi = x + t/(1-t), t in [0, 1)
        part2 = gl_quad(
            lambda t: G(x, x + t / (1.0 - t)) * f_prob(x + t / (1.0 - t)) / ((1.0 - t) ** 2),
            0.0,
            1.0 - 1e-12,
            128,
        )
        y_num = part1 + part2
        max_conv_err = max(max_conv_err, abs(y_num - y_exact(x)))

        # Direct analytic derivatives of y_exact
        yp = (-0.5 * x * x - x + 1.0) * math.exp(-x)
        ypp = (0.5 * x * x - 2.0) * math.exp(-x)
        lhs = ypp + ((x + 1.0) / (x + 2.0)) * yp - (1.0 / (x + 2.0)) * y_exact(x)
        max_sol_ode_err = max(max_sol_ode_err, abs(lhs - f_prob(x)))

    # BCs of y_exact
    y0 = y_exact(0.0)
    yp0 = 1.0  # analytic y'(0)
    assert abs(y0 - yp0) < 1e-15
    assert abs(y_exact(60.0)) < 1e-20
    assert max_conv_err < 1e-11, f"Convolution err: {max_conv_err}"
    assert max_sol_ode_err < 1e-14, f"Solution ODE err: {max_sol_ode_err}"

    # 4. Cross-check Green's function on an independent forcing f2(x) = e^{-2x} cos(x)
    f2 = lambda xi: math.exp(-2.0 * xi) * math.cos(xi)

    def y_from_f2(x):
        p1 = gl_quad(lambda xi: G(x, xi) * f2(xi), 0.0, x, 64) if x > 0 else 0.0
        p2 = gl_quad(
            lambda t: G(x, x + t / (1.0 - t)) * f2(x + t / (1.0 - t)) / ((1.0 - t) ** 2),
            0.0,
            1.0 - 1e-12,
            128,
        )
        return p1 + p2

    max_f2_err = 0.0
    for x in [0.3, 1.1, 2.7]:
        res = abs(L_fd(y_from_f2, x, h=1e-4) - f2(x))
        max_f2_err = max(max_f2_err, res)
    # Check BC at x = 0 for y_from_f2
    h = 1e-5
    y2_0 = y_from_f2(0.0)
    y2p_0 = (-3.0 * y2_0 + 4.0 * y_from_f2(h) - y_from_f2(2.0 * h)) / (2.0 * h)
    bc_f2_err = abs(y2p_0 - y2_0)
    assert max_f2_err < 1e-6, f"Independent forcing ODE err: {max_f2_err}"
    assert bc_f2_err < 1e-8, f"Independent forcing BC0 err: {bc_f2_err}"

    print("IB_2006_P2_q15G verification PASSED")
    print(f"  Homogeneous residuals: y1={max_hom1:.2e}, y2={max_hom2:.2e}")
    print(
        f"  Green's function checks: ODE={max_ode_err:.2e}, BC0={max_bc0_err:.2e}, "
        f"cont={max_cont_err:.2e}, jump={max_jump_err:.2e}, sym={max_sym_err:.2e}"
    )
    print(
        f"  Solution y(x)=(x^2/2+2x+1)e^(-x): conv_err={max_conv_err:.2e}, "
        f"analytic_ODE_err={max_sol_ode_err:.2e}"
    )
    print(f"  Independent forcing f2=e^(-2x)cos(x): ODE_err={max_f2_err:.2e}, BC0_err={bc_f2_err:.2e}")


if __name__ == "__main__":
    main()
