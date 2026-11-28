#!/usr/bin/env python3
"""Computational verification for IB/2014/P3/q7D (method of characteristics for -y u_x + x u_y = 0)."""

import math


def u(x, y):
    return math.cos(x * x + y * y)


def main():
    # 1. Check characteristic curves x(t,s) = -s sin t, y(t,s) = s cos t
    h = 1e-6
    max_ode_err = 0.0
    max_jac_err = 0.0
    for s in [0.0, 0.3, 1.2, 2.7]:
        for t in [0.0, 0.5, 1.7, 3.1, 5.4]:
            x = -s * math.sin(t)
            y = s * math.cos(t)
            dxdt = (-s * math.sin(t + h) - (-s * math.sin(t - h))) / (2.0 * h)
            dydt = (s * math.cos(t + h) - s * math.cos(t - h)) / (2.0 * h)
            max_ode_err = max(max_ode_err, abs(dxdt - (-y)), abs(dydt - x))

            dxds = (-(s + h) * math.sin(t) - (-(s - h) * math.sin(t))) / (2.0 * h)
            dyds = ((s + h) * math.cos(t) - (s - h) * math.cos(t)) / (2.0 * h)
            jac = dxdt * dyds - dxds * dydt
            max_jac_err = max(max_jac_err, abs(jac - (-s)))

            # Check u is constant along characteristics and matches initial data
            assert abs(u(x, y) - math.cos(s * s)) < 1e-14
    assert max_ode_err < 1e-9
    assert max_jac_err < 1e-9

    # 2. Check PDE -y u_x + x u_y = 0 across all four quadrants, axes, and origin
    h_fd = 1e-6
    max_pde_res = 0.0
    grid = [-2.3, -1.0, -0.4, 0.0, 0.4, 1.0, 2.3]
    for x in grid:
        for y in grid:
            ux = (u(x + h_fd, y) - u(x - h_fd, y)) / (2.0 * h_fd)
            uy = (u(x, y + h_fd) - u(x, y - h_fd)) / (2.0 * h_fd)
            res = abs(-y * ux + x * uy)
            max_pde_res = max(max_pde_res, res)
    assert max_pde_res < 1e-8, f"PDE residual too large: {max_pde_res}"

    # 3. Check boundary condition u(0, y) = cos(y^2) for y >= 0
    for y in [0.0, 0.25, 1.0, 1.73, 3.5]:
        assert abs(u(0.0, y) - math.cos(y * y)) < 1e-15

    # 4. Independent RK4 integration of characteristic ODE dx/dt = -y, dy/dt = x from (0, s)
    max_rk4_err = 0.0
    for s in [0.5, 1.5, 2.2]:
        x, y = 0.0, s
        dt = 1e-3
        steps = int(round(2.0 * math.pi / dt))
        for _ in range(steps):
            k1x, k1y = -y, x
            k2x, k2y = -(y + 0.5 * dt * k1y), (x + 0.5 * dt * k1x)
            k3x, k3y = -(y + 0.5 * dt * k2y), (x + 0.5 * dt * k2x)
            k4x, k4y = -(y + dt * k3y), (x + dt * k3x)
            x += (dt / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
            y += (dt / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y)
            max_rk4_err = max(max_rk4_err, abs(u(x, y) - math.cos(s * s)))
    assert max_rk4_err < 1e-11, f"RK4 invariance error: {max_rk4_err}"

    print("IB_2014_P3_q7D verification PASSED:")
    print(f"  Characteristic ODE error = {max_ode_err:.3e}, Jacobian error = {max_jac_err:.3e}")
    print(f"  PDE residual across 49 grid points = {max_pde_res:.3e}")
    print(f"  RK4 full-orbit invariance error = {max_rk4_err:.3e}")


if __name__ == "__main__":
    main()
