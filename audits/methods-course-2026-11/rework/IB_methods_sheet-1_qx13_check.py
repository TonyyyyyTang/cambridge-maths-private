"""Computational verification for IB/methods/sheet-1/qx13 (2001-02 Sheet 1 Q7)."""

import math


def rk4_second_order(rhs_acc, q0, dq0, t_end, n_steps):
    """Integrate q'' = rhs_acc(t, q, dq) from t=0 to t_end using RK4."""
    dt = t_end / n_steps
    t = 0.0
    q = q0
    v = dq0
    history = [(t, q, v)]
    for _ in range(n_steps):
        k1_q = v
        k1_v = rhs_acc(t, q, v)

        k2_q = v + 0.5 * dt * k1_v
        k2_v = rhs_acc(t + 0.5 * dt, q + 0.5 * dt * k1_q, v + 0.5 * dt * k1_v)

        k3_q = v + 0.5 * dt * k2_v
        k3_v = rhs_acc(t + 0.5 * dt, q + 0.5 * dt * k2_q, v + 0.5 * dt * k2_v)

        k4_q = v + dt * k3_v
        k4_v = rhs_acc(t + dt, q + dt * k3_q, v + dt * k3_v)

        q += (dt / 6.0) * (k1_q + 2.0 * k2_q + 2.0 * k3_q + k4_q)
        v += (dt / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v)
        t += dt
        history.append((t, q, v))
    return history


def verify_all():
    pi = math.pi
    param_sets = [
        # (ell, mu, c, A, F0)
        (1.0, 1.0, 1.0, 1.0, 1.0),
        (2.5, 1.7, 3.2, 0.8, 4.5),
        (0.7, 0.4, 1.9, -1.3, 2.1),
    ]

    for ell, mu, c, A, F0 in param_sets:
        T = mu * c * c
        k = pi * c / ell

        # --- Part (i): Free damped motion ---
        def y_free(x, t):
            return A * (1.0 + k * t) * math.exp(-k * t) * math.sin(pi * x / ell)

        # Check BCs and ICs
        for t in [0.0, 0.25, 1.1]:
            assert abs(y_free(0.0, t)) < 1e-15
            assert abs(y_free(ell, t)) < 1e-15

        dt_fd = 1e-4 / k
        dx_fd = 1e-4 * ell
        for frac in [0.15, 0.5, 0.77]:
            x = frac * ell
            assert abs(y_free(x, 0.0) - A * math.sin(pi * x / ell)) < 1e-14
            yt0_fd = (y_free(x, dt_fd) - y_free(x, -dt_fd)) / (2.0 * dt_fd)
            assert abs(yt0_fd) < 1e-7

        # Check PDE c^2 y_xx = y_tt + 2k y_t via finite differences
        max_pde1_res = 0.0
        for x in [0.2 * ell, 0.5 * ell, 0.81 * ell]:
            for t in [0.1 / k, 0.7 / k, 2.3 / k]:
                yxx = (y_free(x + dx_fd, t) - 2.0 * y_free(x, t) + y_free(x - dx_fd, t)) / (dx_fd ** 2)
                ytt = (y_free(x, t + dt_fd) - 2.0 * y_free(x, t) + y_free(x, t - dt_fd)) / (dt_fd ** 2)
                yt = (y_free(x, t + dt_fd) - y_free(x, t - dt_fd)) / (2.0 * dt_fd)
                max_pde1_res = max(max_pde1_res, abs(c * c * yxx - (ytt + 2.0 * k * yt)))
        assert max_pde1_res < 1e-5 * abs(A) * (k ** 2), f"max_pde1_res={max_pde1_res}"

        # Independent RK4 check of mode ODE q'' + 2k q' + k^2 q = 0
        hist1 = rk4_second_order(lambda t, q, v: -2.0 * k * v - (k ** 2) * q, A, 0.0, 4.0 / k, 4000)
        max_rk4_err1 = max(abs(q - A * (1.0 + k * t) * math.exp(-k * t)) for t, q, _ in hist1)
        assert max_rk4_err1 < 1e-10, f"max_rk4_err1={max_rk4_err1}"

        # --- Part (ii): Forced oscillation ---
        amp1 = F0 / (2.0 * mu * k * k)
        amp2 = F0 * ell * ell / (2.0 * (pi ** 2) * T)
        assert abs(amp1 - amp2) < 1e-14

        def y_forced(x, t):
            return amp1 * math.sin(pi * x / ell) * math.sin(k * t)

        # Check BCs and forced PDE: mu y_tt + 2k mu y_t - T y_xx = F0 sin(pi x/l) cos(pi c t/l)
        max_pde2_res = 0.0
        for x in [0.23 * ell, 0.5 * ell, 0.79 * ell]:
            for t in [0.0, 0.4 / k, 1.7 / k, 3.1 / k]:
                assert abs(y_forced(0.0, t)) < 1e-15
                assert abs(y_forced(ell, t)) < 1e-15
                yxx = (y_forced(x + dx_fd, t) - 2.0 * y_forced(x, t) + y_forced(x - dx_fd, t)) / (dx_fd ** 2)
                ytt = (y_forced(x, t + dt_fd) - 2.0 * y_forced(x, t) + y_forced(x, t - dt_fd)) / (dt_fd ** 2)
                yt = (y_forced(x, t + dt_fd) - y_forced(x, t - dt_fd)) / (2.0 * dt_fd)
                lhs = mu * ytt + 2.0 * k * mu * yt - T * yxx
                rhs = F0 * math.sin(pi * x / ell) * math.cos(k * t)
                max_pde2_res = max(max_pde2_res, abs(lhs - rhs))
        assert max_pde2_res < 1e-5 * abs(F0), f"max_pde2_res={max_pde2_res}"

        # Independent RK4 check of forced mode ODE with Part (i) ICs (q(0)=A, q'(0)=0)
        # Exact solution: (A + (k*A - F0/(2*mu*k))*t)*exp(-k*t) + (F0/(2*mu*k^2))*sin(k*t)
        hist2 = rk4_second_order(
            lambda t, q, v: (F0 / mu) * math.cos(k * t) - 2.0 * k * v - (k ** 2) * q,
            A, 0.0, 25.0 / k, 25000
        )
        max_rk4_err2 = 0.0
        for t, q, _ in hist2:
            q_exact = (A + (k * A - F0 / (2.0 * mu * k)) * t) * math.exp(-k * t) + amp1 * math.sin(k * t)
            max_rk4_err2 = max(max_rk4_err2, abs(q - q_exact))
        assert max_rk4_err2 < 1e-9, f"max_rk4_err2={max_rk4_err2}"

        # Check convergence to steady-state forced oscillation at large t (t in [20/k, 25/k])
        max_steady_diff = max(abs(q - amp1 * math.sin(k * t)) for t, q, _ in hist2 if t >= 20.0 / k)
        assert max_steady_diff < 1e-6 * (abs(A) + abs(amp1)), f"max_steady_diff={max_steady_diff}"

        print(
            f"Params (l={ell}, mu={mu}, c={c}, A={A}, F0={F0}): "
            f"PDE1_res={max_pde1_res:.2e}, RK4_1={max_rk4_err1:.2e}, "
            f"PDE2_res={max_pde2_res:.2e}, RK4_2={max_rk4_err2:.2e}, "
            f"steady_diff(t>=20/k)={max_steady_diff:.2e}"
        )


if __name__ == "__main__":
    verify_all()
    print("ALL CHECKS PASSED for IB/methods/sheet-1/qx13.")
