"""Computational verification for IB/2008/P1/q14D (Fermat's principle / light ray)."""
import math


def gl_quad(func, a, b, n_panels=64):
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


def main():
    # 1. General c(x) test: c(x) = 1 + 0.3 x + 0.1 x^2 on [x1, x2] = [0.5, 2.0]
    c_gen = lambda x: 1.0 + 0.3 * x + 0.1 * x * x
    x1, x2, y1 = 0.5, 2.0, 1.0
    k = 3.0  # > c(2.0) = 2.0
    yp_gen = lambda x: c_gen(x) / math.sqrt(k * k - c_gen(x) ** 2)
    y2 = y1 + gl_quad(yp_gen, x1, x2, 128)

    # Check first integral d f / d y' = 1 / k across [x1, x2]
    max_fi_err = 0.0
    for i in range(1, 40):
        x = x1 + (x2 - x1) * i / 40.0
        yp = yp_gen(x)
        df_dyp = yp / (c_gen(x) * math.sqrt(1.0 + yp * yp))
        max_fi_err = max(max_fi_err, abs(df_dyp - 1.0 / k))
    assert max_fi_err < 1e-15, f"First integral err: {max_fi_err}"

    # Check that T is strictly minimised among perturbations vanishing at x1, x2
    T_opt_gen = gl_quad(lambda x: math.sqrt(1.0 + yp_gen(x) ** 2) / c_gen(x), x1, x2, 128)
    for m in (1, 2, 3):
        etap = lambda x, mm=m: (mm * math.pi / (x2 - x1)) * math.cos(mm * math.pi * (x - x1) / (x2 - x1))
        for amp in (-0.2, -0.05, 0.05, 0.2):
            T_pert = gl_quad(
                lambda x: math.sqrt(1.0 + (yp_gen(x) + amp * etap(x)) ** 2) / c_gen(x),
                x1,
                x2,
                128,
            )
            assert T_pert > T_opt_gen + 1e-10
    print(f"General c(x) test: y2={y2:.10f}, T_opt={T_opt_gen:.10f}, first_integral_err={max_fi_err:.2e}")

    # 2. Specific case: c(x) = a * x, k = a * x2 -> (y2 - y)^2 = x2^2 - x^2
    cases = [
        (1.0, 0.6, 1.0, 2.0),
        (2.5, 1.0, 3.0, 5.0),
        (0.7, 0.0, 2.0, 3.5),
    ]
    for a, x1, x2, y2 in cases:
        k_val = a * x2
        y1 = y2 - math.sqrt(x2 * x2 - x1 * x1)
        assert y1 > 0.0 and y2 > y1

        y_path = lambda x: y2 - math.sqrt(x2 * x2 - x * x)
        yp_path = lambda x: x / math.sqrt(x2 * x2 - x * x)

        # Check equation (y2 - y)^2 = x2^2 - x^2, ODE, first integral, and BCs
        assert abs(y_path(x1) - y1) < 1e-15
        assert abs(y_path(x2) - y2) < 1e-15
        max_circle_err = 0.0
        max_ode_err = 0.0
        max_fi2_err = 0.0
        for i in range(1, 50):
            x = x1 + (x2 - x1) * i / 50.0
            yy = y_path(x)
            yp = yp_path(x)
            cx = a * x
            max_circle_err = max(max_circle_err, abs((y2 - yy) ** 2 - (x2 * x2 - x * x)))
            rhs = cx / math.sqrt(k_val * k_val - cx * cx)
            max_ode_err = max(max_ode_err, abs(yp - rhs) / (1.0 + abs(yp)))
            if x > 0:
                df_dyp = yp / (cx * math.sqrt(1.0 + yp * yp))
                max_fi2_err = max(max_fi2_err, abs(df_dyp - 1.0 / k_val))

        # Check endpoint integral y2 - y1 = int_{x1}^{x2} c(x)/sqrt(k^2 - c(x)^2) dx
        # Substitute x = x2 * sin(phi) to remove integrable endpoint singularity at x = x2
        phi1 = math.asin(x1 / x2)
        integral_val = gl_quad(
            lambda phi: (a * x2 * math.sin(phi)) / (a * x2 * math.cos(phi)) * (x2 * math.cos(phi)),
            phi1,
            0.5 * math.pi,
            64,
        )
        # Also test substitution u = sqrt(x2 - x) => x = x2 - u^2, dx = -2 u du
        # Using k - c(x) = a * u^2 avoids floating-point cancellation near u = 0.
        u1 = math.sqrt(x2 - x1)
        integral_u = gl_quad(
            lambda u: (a * (x2 - u * u))
            / math.sqrt((a * u * u) * (k_val + a * (x2 - u * u)))
            * (2.0 * u),
            0.0,
            u1,
            128,
        )
        assert abs(integral_val - (y2 - y1)) < 1e-14
        assert abs(integral_u - (y2 - y1)) < 1e-13
        assert max_circle_err < 1e-14
        assert max_ode_err < 1e-14
        assert max_fi2_err < 1e-14

        # Check travel time minimality when x1 > 0 (so T is finite; exact T = (1/a) arccosh(x2/x1))
        if x1 > 0:
            T_exact = math.acosh(x2 / x1) / a
            # Compare against straight-line path between (x1, y1) and (x2, y2)
            slope_line = (y2 - y1) / (x2 - x1)
            T_line = (math.sqrt(1.0 + slope_line * slope_line) / a) * math.log(x2 / x1)
            assert T_line > T_exact + 1e-6
            # Numerical quadrature of T along optimal path using u = sqrt(x2 - x)
            T_num = gl_quad(
                lambda u: (x2 / (a * (x2 - u * u) * math.sqrt(2.0 * x2 - u * u))) * 2.0,
                0.0,
                u1,
                128,
            )
            assert abs(T_num - T_exact) < 1e-13
            print(
                f"c(x)={a}x on [{x1},{x2}]: y1={y1:.6f}, int_err={abs(integral_u-(y2-y1)):.2e}, "
                f"T_opt={T_exact:.10f} < T_line={T_line:.10f}"
            )
        else:
            print(
                f"c(x)={a}x on [{x1},{x2}]: y1={y1:.6f}, int_err={abs(integral_u-(y2-y1)):.2e} (x1=0 boundary)"
            )

    print("IB_2008_P1_q14D verification PASSED")


if __name__ == "__main__":
    main()
