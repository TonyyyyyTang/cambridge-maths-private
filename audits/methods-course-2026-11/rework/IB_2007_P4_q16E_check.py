"""Computational verification for IB/2007/P4/q16E (Variational road travel time)."""
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


def F(r, rp, k):
    return math.sqrt(rp * rp + r * r) / (k * r * r)


def travel_time(r_fn, rp_fn, k):
    return gl_quad(lambda th: F(r_fn(th), rp_fn(th), k), 0.0, 0.5 * math.pi, 128)


def main():
    params = [(1.0, 1.0), (2.5, 0.7), (0.4, 3.2)]
    for a, k in params:
        # 1. Circular road r = a
        T_circ_num = travel_time(lambda th: a, lambda th: 0.0, k)
        T_circ_exact = math.pi / (2.0 * k * a)
        assert abs(T_circ_num - T_circ_exact) < 1e-14

        # 2. Optimal road r*(theta) = a*(cos(theta) + sin(theta)) = sqrt(2)*a*cos(theta - pi/4)
        r_opt = lambda th: a * (math.cos(th) + math.sin(th))
        rp_opt = lambda th: a * (-math.sin(th) + math.cos(th))

        # Check equivalence with sqrt(2)*a*cos(theta - pi/4) and BCs
        assert abs(r_opt(0.0) - a) < 1e-15
        assert abs(r_opt(0.5 * math.pi) - a) < 1e-15
        for i in range(21):
            th = 0.5 * math.pi * i / 20.0
            alt = math.sqrt(2.0) * a * math.cos(th - 0.25 * math.pi)
            assert abs(r_opt(th) - alt) < 1e-14

        # Check Beltrami first integral F - r' dF/dr' = 1/(k*sqrt(2)*a)
        C_expected = 1.0 / (k * math.sqrt(2.0) * a)
        max_beltrami_err = 0.0
        max_el_err = 0.0
        h = 1e-5
        for i in range(1, 20):
            th = 0.5 * math.pi * i / 20.0
            r = r_opt(th)
            rp = rp_opt(th)
            dF_drp = rp / (k * r * r * math.sqrt(rp * rp + r * r))
            beltrami = F(r, rp, k) - rp * dF_drp
            max_beltrami_err = max(max_beltrami_err, abs(beltrami - C_expected))

            # Full Euler-Lagrange: d/dtheta (dF/dr') - dF/dr = 0
            def g(t):
                rt = r_opt(t)
                rpt = rp_opt(t)
                return rpt / (k * rt * rt * math.sqrt(rpt * rpt + rt * rt))

            d_dF_drp = (g(th + h) - g(th - h)) / (2.0 * h)
            dF_dr = -(2.0 * rp * rp + r * r) / (k * (r ** 3) * math.sqrt(rp * rp + r * r))
            max_el_err = max(max_el_err, abs(d_dF_drp - dF_dr))

        assert max_beltrami_err < 1e-14, f"Beltrami err: {max_beltrami_err}"
        assert max_el_err < 1e-8, f"Euler-Lagrange err: {max_el_err}"

        # Check optimal travel time T_min = sqrt(2) / (k * a)
        T_opt_num = travel_time(r_opt, rp_opt, k)
        T_opt_exact = math.sqrt(2.0) / (k * a)
        assert abs(T_opt_num - T_opt_exact) < 1e-13
        assert T_opt_exact < T_circ_exact

        # 3. Check strict minimality against diverse admissible perturbations vanishing at 0 and pi/2
        perturbations = [
            (lambda th: math.sin(2.0 * th), lambda th: 2.0 * math.cos(2.0 * th)),
            (lambda th: math.sin(4.0 * th), lambda th: 4.0 * math.cos(4.0 * th)),
            (lambda th: math.sin(6.0 * th), lambda th: 6.0 * math.cos(6.0 * th)),
            (
                lambda th: th * (0.5 * math.pi - th) * (1.0 + th),
                lambda th: (0.5 * math.pi - 2.0 * th) * (1.0 + th) + th * (0.5 * math.pi - th),
            ),
            (
                lambda th: math.sin(2.0 * th) - 0.5 * math.sin(4.0 * th) + 0.3 * math.sin(6.0 * th),
                lambda th: 2.0 * math.cos(2.0 * th) - 2.0 * math.cos(4.0 * th) + 1.8 * math.cos(6.0 * th),
            ),
        ]
        for eta, etap in perturbations:
            # First variation = 0 (use eps1 = 1e-6*a), second variation > 0 (use eps2 = 1e-4*a)
            eps1 = 1e-6 * a
            T_p1 = travel_time(lambda th: r_opt(th) + eps1 * eta(th), lambda th: rp_opt(th) + eps1 * etap(th), k)
            T_m1 = travel_time(lambda th: r_opt(th) - eps1 * eta(th), lambda th: rp_opt(th) - eps1 * etap(th), k)
            d1 = (T_p1 - T_m1) / (2.0 * eps1)
            eps2 = 1e-4 * a
            T_p2 = travel_time(lambda th: r_opt(th) + eps2 * eta(th), lambda th: rp_opt(th) + eps2 * etap(th), k)
            T_m2 = travel_time(lambda th: r_opt(th) - eps2 * eta(th), lambda th: rp_opt(th) - eps2 * etap(th), k)
            d2 = (T_p2 - 2.0 * T_opt_num + T_m2) / (eps2 * eps2)
            assert abs(d1) < 1e-9 / (k * a * a), f"First variation nonzero: {d1}"
            assert d2 > 0.0, f"Second variation not positive: {d2}"

            # Finite perturbations
            for amp in [-0.3 * a, -0.1 * a, 0.1 * a, 0.3 * a, 0.6 * a]:
                # Ensure r > 0 on [0, pi/2]
                if all(r_opt(0.5 * math.pi * j / 40.0) + amp * eta(0.5 * math.pi * j / 40.0) > 0.05 * a for j in range(41)):
                    T_pert = travel_time(
                        lambda th: r_opt(th) + amp * eta(th),
                        lambda th: rp_opt(th) + amp * etap(th),
                        k,
                    )
                    assert T_pert > T_opt_num + 1e-12

        print(
            f"Params (a={a}, k={k}): T_circ={T_circ_exact:.10f}, T_min={T_opt_exact:.10f} "
            f"(quad_err={abs(T_opt_num - T_opt_exact):.2e}), Beltrami_err={max_beltrami_err:.2e}, EL_err={max_el_err:.2e}"
        )

    print("IB_2007_P4_q16E verification PASSED")


if __name__ == "__main__":
    main()
