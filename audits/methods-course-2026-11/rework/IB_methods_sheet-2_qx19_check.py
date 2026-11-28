"""Computational verification for IB/methods/sheet-2/qx19."""
import math


def gauss_legendre_8():
    nodes = [
        -0.96028985649753623,
        -0.79666647741362674,
        -0.52553240991632899,
        -0.18343464249564980,
         0.18343464249564980,
         0.52553240991632899,
         0.79666647741362674,
         0.96028985649753623,
    ]
    weights = [
        0.10122853629037626,
        0.22238103445337447,
        0.31370664587788729,
        0.36268378337836198,
        0.36268378337836198,
        0.31370664587788729,
        0.22238103445337447,
        0.10122853629037626,
    ]
    return nodes, weights


def integrate_gl(func, a, b, n_panels=32):
    nodes, weights = gauss_legendre_8()
    h = (b - a) / n_panels
    total = 0.0
    for i in range(n_panels):
        mid = a + (i + 0.5) * h
        half = 0.5 * h
        s = 0.0
        for xi, wi in zip(nodes, weights):
            s += wi * func(mid + half * xi)
        total += half * s
    return total


def L_fd(func, x, h=1e-5):
    """Finite-difference evaluation of L_x y = -(1/x^2) d/dx(x^2 dy/dx) + y."""
    def x2_yp(t):
        return (t * t) * (func(t + h) - func(t - h)) / (2.0 * h)
    d_x2_yp = (x2_yp(x + h) - x2_yp(x - h)) / (2.0 * h)
    return -d_x2_yp / (x * x) + func(x)


def G(x, a):
    """Green's function for L_x G(x, a) = delta(x - a) bounded as x -> 0 and x -> inf."""
    if x == 0.0:
        return a * math.exp(-a)
    if x <= a:
        return a * math.exp(-a) * math.sinh(x) / x
    else:
        return a * math.sinh(a) * math.exp(-x) / x


def y_exact(x, R):
    """Analytical solution of L_x y = 1_{0<=x<=R}."""
    A = -(R + 1.0) * math.exp(-R)
    B = R * math.cosh(R) - math.sinh(R)
    if x <= R:
        if x == 0.0:
            return 1.0 + A
        return 1.0 + A * math.sinh(x) / x
    else:
        return B * math.exp(-x) / x


def main():
    # 1. Verify homogeneous solutions y1 = sinh(x)/x and y2 = exp(-x)/x
    y1 = lambda x: math.sinh(x) / x
    y2 = lambda x: math.exp(-x) / x
    max_hom_err = 0.0
    for x in [0.3, 1.1, 2.7, 5.0]:
        max_hom_err = max(max_hom_err, abs(L_fd(y1, x)), abs(L_fd(y2, x)))
    assert max_hom_err < 1e-6, f"Homogeneous ODE error: {max_hom_err}"

    # 2. Verify Green's function ODE, continuity, jump condition, and weighted symmetry
    max_g_ode_err = 0.0
    max_cont_err = 0.0
    max_jump_err = 0.0
    max_sym_err = 0.0
    eps = 1e-6
    for a in [0.4, 1.3, 2.8]:
        for x in [0.2, 0.9, 1.9, 4.1]:
            if abs(x - a) > 0.05:
                max_g_ode_err = max(max_g_ode_err, abs(L_fd(lambda t, aa=a: G(t, aa), x)))
            # Weighted symmetry x^2 G(x,a) = a^2 G(a,x)
            sym_err = abs((x * x) * G(x, a) - (a * a) * G(a, x))
            max_sym_err = max(max_sym_err, sym_err)
        # Continuity at x = a
        cont_err = abs(G(a - 1e-12, a) - G(a + 1e-12, a))
        max_cont_err = max(max_cont_err, cont_err)
        # Jump in dG/dx at x = a must be -1 (and jump in x^2 dG/dx must be -a^2)
        dG_left = (G(a, a) - G(a - eps, a)) / eps
        dG_right = (G(a + eps, a) - G(a, a)) / eps
        jump_err = abs((dG_right - dG_left) - (-1.0))
        max_jump_err = max(max_jump_err, jump_err)

    assert max_g_ode_err < 1e-6
    assert max_cont_err < 1e-10
    assert max_jump_err < 1e-4
    assert max_sym_err < 1e-14

    # 3. Verify Green's reproduction: for smooth psi(x) = exp(-x^2) (bounded at 0, fast decay at inf),
    # int_0^inf G(x, a) (L_a psi)(a) da = psi(x)
    # Note: psi'(a) = -2a e^{-a^2}, (a^2 psi')' = (-2a^3 e^{-a^2})' = (-6a^2 + 4a^4) e^{-a^2}
    # so L_a psi(a) = (6 - 4a^2 + 1) e^{-a^2} = (7 - 4a^2) e^{-a^2}
    def psi(x):
        return math.exp(-x * x)
    def L_psi(a):
        return (7.0 - 4.0 * a * a) * math.exp(-a * a)

    max_reprod_err = 0.0
    for x in [0.25, 1.0, 2.3]:
        # Split integral at a = x to integrate smoothly across the corner
        val = (integrate_gl(lambda a, xx=x: G(xx, a) * L_psi(a), 0.0, x, 32)
               + integrate_gl(lambda a, xx=x: G(xx, a) * L_psi(a), x, 10.0, 64))
        err = abs(val - psi(x))
        max_reprod_err = max(max_reprod_err, err)
    assert max_reprod_err < 1e-12, f"Green reproduction error: {max_reprod_err}"

    # 4. Verify step-RHS solution against numerical quadrature int_0^R G(x, a) da, ODE, and C^1 matching at R
    max_quad_err = 0.0
    max_step_ode_err = 0.0
    max_c0_R_err = 0.0
    max_c1_R_err = 0.0
    max_j2_R_err = 0.0
    h = 1e-5

    for R in [0.7, 1.5, 3.2]:
        # Quadrature vs exact formula
        for x in [0.0, 0.25 * R, 0.7 * R, R, 1.4 * R, 2.5 * R]:
            if 0.0 < x < R:
                num_y = (integrate_gl(lambda a, xx=x: G(xx, a), 0.0, x, 32)
                         + integrate_gl(lambda a, xx=x: G(xx, a), x, R, 32))
            else:
                num_y = integrate_gl(lambda a, xx=x: G(xx, a), 0.0, R, 32)
            ex_y = y_exact(x, R)
            err = abs(num_y - ex_y)
            max_quad_err = max(max_quad_err, err)

        # ODE interior checks
        for x in [0.3 * R, 0.8 * R]:
            max_step_ode_err = max(max_step_ode_err, abs(L_fd(lambda t, RR=R: y_exact(t, RR), x) - 1.0))
        for x in [1.3 * R, 2.5 * R]:
            max_step_ode_err = max(max_step_ode_err, abs(L_fd(lambda t, RR=R: y_exact(t, RR), x) - 0.0))

        # C^0 and C^1 continuity at x = R, and jump of +1 in y'' at x = R
        c0_err = abs(y_exact(R, R) - y_exact(R + 1e-13, R))
        yp_left = (y_exact(R, R) - y_exact(R - eps, R)) / eps
        yp_right = (y_exact(R + eps, R) - y_exact(R, R)) / eps
        c1_err = abs(yp_right - yp_left)
        ypp_left = (y_exact(R, R) - 2.0 * y_exact(R - h, R) + y_exact(R - 2.0 * h, R)) / (h * h)
        ypp_right = (y_exact(R + 2.0 * h, R) - 2.0 * y_exact(R + h, R) + y_exact(R, R)) / (h * h)
        j2_err = abs((ypp_right - ypp_left) - 1.0)

        max_c0_R_err = max(max_c0_R_err, c0_err)
        max_c1_R_err = max(max_c1_R_err, c1_err)
        max_j2_R_err = max(max_j2_R_err, j2_err)

    assert max_quad_err < 1e-13, f"Quadrature error: {max_quad_err}"
    assert max_step_ode_err < 1e-6, f"Step ODE error: {max_step_ode_err}"
    assert max_c0_R_err < 1e-11
    assert max_c1_R_err < 1e-5
    assert max_j2_R_err < 1e-3

    print("IB_methods_sheet-2_qx19_check: ALL PASSED")
    print(f"  homogeneous ODE err: {max_hom_err:.3e}")
    print(f"  Green's fn: ODE err={max_g_ode_err:.3e}, cont err={max_cont_err:.3e}, jump err={max_jump_err:.3e}, sym err={max_sym_err:.3e}")
    print(f"  Green reproduction err (psi=exp(-x^2)): {max_reprod_err:.3e}")
    print(f"  step-RHS: quad vs closed-form err={max_quad_err:.3e}, ODE err={max_step_ode_err:.3e}")
    print(f"  matching at x=R: C^0 err={max_c0_R_err:.3e}, C^1 err={max_c1_R_err:.3e}, y'' jump(+1) err={max_j2_R_err:.3e}")


if __name__ == "__main__":
    main()
