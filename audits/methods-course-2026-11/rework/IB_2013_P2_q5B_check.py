#!/usr/bin/env python3
"""Verification script for IB/2013/P2/q5B (Method of characteristics)."""

import math


def u_sol(x, y):
    return y / x


def main():
    # 1. Verify PDE x*u_x + (x+y)*u_y = 1 analytically and by central finite differences
    h = 1e-6
    max_pde_err = 0.0
    for x in [0.1, 0.5, 1.0, 2.3, 7.5, -0.4, -3.2]:
        for y in [-5.0, -1.2, 0.0, 0.7, 4.4]:
            ux = (u_sol(x + h, y) - u_sol(x - h, y)) / (2 * h)
            uy = (u_sol(x, y + h) - u_sol(x, y - h)) / (2 * h)
            lhs = x * ux + (x + y) * uy
            err = abs(lhs - 1.0)
            max_pde_err = max(max_pde_err, err)
            assert err < 1e-8, f"PDE error {err} at ({x}, {y})"
    print(f"[PASS] Finite-difference PDE check max |x u_x + (x+y) u_y - 1| = {max_pde_err:.3e}")

    # 2. Verify Cauchy data u(1, y) = y
    for y in [-10.0, -2.5, 0.0, 1.0, 3.14159, 8.7]:
        assert abs(u_sol(1.0, y) - y) < 1e-15
    print("[PASS] Cauchy data u(1, y) = y verified exactly")

    # 3. Integrate characteristic ODEs dx/dt = x, dy/dt = x+y, du/dt = 1 via RK4 from (1, s, s)
    max_traj_err = 0.0
    max_u_err = 0.0
    for s in [-4.0, -1.5, 0.0, 0.8, 3.5]:
        for direction in [1.0, -1.0]:
            T = 2.5 * direction
            steps = 2500
            dt = T / steps
            x, y, u = 1.0, s, s
            for step in range(1, steps + 1):
                t = step * dt
                # RK4 step for (x, y, u)' = (x, x+y, 1)
                k1x, k1y, k1u = x, x + y, 1.0
                k2x, k2y, k2u = x + 0.5 * dt * k1x, (x + 0.5 * dt * k1x) + (y + 0.5 * dt * k1y), 1.0
                k3x, k3y, k3u = x + 0.5 * dt * k2x, (x + 0.5 * dt * k2x) + (y + 0.5 * dt * k2y), 1.0
                k4x, k4y, k4u = x + dt * k3x, (x + dt * k3x) + (y + dt * k3y), 1.0
                x += (dt / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
                y += (dt / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y)
                u += (dt / 6.0) * (k1u + 2 * k2u + 2 * k3u + k4u)

                if step % 250 == 0:
                    x_ex = math.exp(t)
                    y_ex = (t + s) * math.exp(t)
                    u_ex = t + s
                    max_traj_err = max(max_traj_err, abs(x - x_ex), abs(y - y_ex), abs(u - u_ex))
                    max_u_err = max(max_u_err, abs(u - u_sol(x, y)))

    assert max_traj_err < 1e-11, f"Characteristic trajectory error: {max_traj_err}"
    assert max_u_err < 1e-11, f"Solution along characteristic error: {max_u_err}"
    print(f"[PASS] RK4 characteristic ODE check: max traj err = {max_traj_err:.3e}, max |u - y/x| = {max_u_err:.3e}")

    # 4. Verify non-uniqueness remark on x < 0: u_alt(x, y) = y/x + F(y/x - ln(-x)) also satisfies PDE on x < 0
    def u_alt(x, y):
        eta = y / x - math.log(-x)
        return y / x + math.sin(eta) + 0.3 * eta * eta

    max_alt_err = 0.0
    for x in [-0.2, -1.0, -3.5]:
        for y in [-2.0, 0.5, 4.0]:
            ux = (u_alt(x + h, y) - u_alt(x - h, y)) / (2 * h)
            uy = (u_alt(x, y + h) - u_alt(x, y - h)) / (2 * h)
            err = abs(x * ux + (x + y) * uy - 1.0)
            max_alt_err = max(max_alt_err, err)
            assert err < 1e-7
    print(f"[PASS] Non-uniqueness on x < 0 verified: max PDE err for u + F(y/x - ln(-x)) = {max_alt_err:.3e}")


if __name__ == "__main__":
    main()
