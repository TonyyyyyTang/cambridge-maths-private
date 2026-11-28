#!/usr/bin/env python3
"""Computational verification for IB/methods/sheet-4/q03.

Verifies:
1. Characteristic coordinates xi = x + 2*sqrt(-y), eta = x - 2*sqrt(-y) in y < 0
   satisfy the characteristic equation, have positive Jacobian 2/sqrt(-y), and
   invert smoothly via x = (xi+eta)/2, y = -(xi-eta)^2/16.
2. Exact operator identity u_xx + y*u_yy + (1/2)*u_y = 4*U_{eta,xi} for a
   generic non-solution smooth function U(xi, eta).
3. General solution u(x,y) = f(x + 2*sqrt(-y)) + g(x - 2*sqrt(-y)) solves the
   original PDE in the hyperbolic region y < 0.
"""

import math


def xi_coord(x, y):
    return x + 2.0 * math.sqrt(-y)


def eta_coord(x, y):
    return x - 2.0 * math.sqrt(-y)


def inv_coords(xi, eta):
    x = 0.5 * (xi + eta)
    y = -((xi - eta) ** 2) / 16.0
    return x, y


def cd1_x(func, x, y, h=1e-5):
    return (-func(x + 2.0 * h, y) + 8.0 * func(x + h, y) - 8.0 * func(x - h, y) + func(x - 2.0 * h, y)) / (12.0 * h)


def cd1_y(func, x, y, h=1e-5):
    return (-func(x, y + 2.0 * h) + 8.0 * func(x, y + h) - 8.0 * func(x, y - h) + func(x, y - 2.0 * h)) / (12.0 * h)


def cd2_xx(func, x, y, h=1e-3):
    return (
        -func(x + 2.0 * h, y)
        + 16.0 * func(x + h, y)
        - 30.0 * func(x, y)
        + 16.0 * func(x - h, y)
        - func(x - 2.0 * h, y)
    ) / (12.0 * h * h)


def cd2_yy(func, x, y, h=1e-3):
    return (
        -func(x, y + 2.0 * h)
        + 16.0 * func(x, y + h)
        - 30.0 * func(x, y)
        + 16.0 * func(x, y - h)
        - func(x, y - 2.0 * h)
    ) / (12.0 * h * h)


def cd2_xy(func, x, y, h=1e-3):
    # Richardson extrapolation of 2nd-order mixed difference: O(h^4)
    def d2(step):
        return (
            func(x + step, y + step)
            - func(x + step, y - step)
            - func(x - step, y + step)
            + func(x - step, y - step)
        ) / (4.0 * step * step)

    return (4.0 * d2(0.5 * h) - d2(h)) / 3.0


def main():
    test_points = [
        (-2.3, -0.05),
        (0.0, -0.5),
        (1.4, -1.0),
        (-0.7, -3.7),
        (4.2, -10.5),
    ]

    # 1. Verify characteristic equation, Jacobian, and inverse map
    max_char_err = 0.0
    max_jac_err = 0.0
    max_inv_err = 0.0
    for x, y in test_points:
        xi = xi_coord(x, y)
        eta = eta_coord(x, y)
        assert xi > eta, "Expected xi > eta in y < 0"

        x_back, y_back = inv_coords(xi, eta)
        max_inv_err = max(max_inv_err, abs(x_back - x), abs(y_back - y))

        h = 1e-4 if abs(y) > 0.1 else 1e-5
        xi_x = cd1_x(xi_coord, x, y, h)
        xi_y = cd1_y(xi_coord, x, y, h)
        eta_x = cd1_x(eta_coord, x, y, h)
        eta_y = cd1_y(eta_coord, x, y, h)

        char_xi = xi_x ** 2 + y * (xi_y ** 2)
        char_eta = eta_x ** 2 + y * (eta_y ** 2)
        max_char_err = max(max_char_err, abs(char_xi), abs(char_eta))

        jac_num = xi_x * eta_y - xi_y * eta_x
        jac_exact = 2.0 / math.sqrt(-y)
        max_jac_err = max(max_jac_err, abs(jac_num - jac_exact))

    assert max_inv_err < 1e-14, f"Inverse map error too large: {max_inv_err}"
    assert max_char_err < 1e-9, f"Characteristic eqn error too large: {max_char_err}"
    assert max_jac_err < 1e-8, f"Jacobian error too large: {max_jac_err}"

    # 2. Verify operator identity u_xx + y*u_yy + (1/2)*u_y = 4*U_{eta,xi}
    def U_generic(xi, eta):
        return (
            math.sin(0.8 * xi) * math.cos(1.1 * eta)
            + 0.05 * (xi ** 2) * (eta ** 3)
            + math.exp(0.15 * xi - 0.1 * eta)
        )

    def U_etaxi_exact(xi, eta):
        return (
            -0.88 * math.cos(0.8 * xi) * math.sin(1.1 * eta)
            + 0.05 * 6.0 * xi * (eta ** 2)
            - 0.015 * math.exp(0.15 * xi - 0.1 * eta)
        )

    def u_from_U(x, y):
        return U_generic(xi_coord(x, y), eta_coord(x, y))

    max_op_err = 0.0
    max_op_rel_err = 0.0
    for x, y in test_points:
        h = 3e-3 if abs(y) > 0.2 else 5e-4
        lhs = (
            cd2_xx(u_from_U, x, y, h)
            + y * cd2_yy(u_from_U, x, y, h)
            + 0.5 * cd1_y(u_from_U, x, y, h)
        )
        xi = xi_coord(x, y)
        eta = eta_coord(x, y)
        rhs = 4.0 * U_etaxi_exact(xi, eta)
        rhs_fd = 4.0 * cd2_xy(U_generic, xi, eta, 5e-3)
        assert abs(rhs - rhs_fd) < 1e-8, f"Analytic vs FD U_etaxi mismatch: {rhs} vs {rhs_fd}"
        err = abs(lhs - rhs)
        max_op_err = max(max_op_err, err)
        max_op_rel_err = max(max_op_rel_err, err / max(1.0, abs(rhs)))

    assert max_op_err < 1e-6, f"Operator identity error too large: {max_op_err}"
    assert max_op_rel_err < 1e-8, f"Operator identity relative error too large: {max_op_rel_err}"

    # 3. Verify general solution u(x,y) = f(xi) + g(eta) satisfies the PDE
    def f(s):
        return math.sin(1.3 * s) + math.exp(-0.1 * s * s)

    def g(s):
        return math.cos(0.7 * s) + (s ** 3) / (1.0 + s * s)

    def u_sol(x, y):
        return f(xi_coord(x, y)) + g(eta_coord(x, y))

    max_sol_err = 0.0
    for x, y in test_points:
        h = 1e-3 if abs(y) > 0.2 else 2e-4
        res = (
            cd2_xx(u_sol, x, y, h)
            + y * cd2_yy(u_sol, x, y, h)
            + 0.5 * cd1_y(u_sol, x, y, h)
        )
        max_sol_err = max(max_sol_err, abs(res))

    assert max_sol_err < 1e-7, f"General solution PDE residual too large: {max_sol_err}"

    print("IB_methods_sheet-4_q03 verification PASSED:")
    print(f"  max inverse map error      = {max_inv_err:.3e}")
    print(f"  max characteristic eq err  = {max_char_err:.3e}")
    print(f"  max Jacobian error         = {max_jac_err:.3e}")
    print(f"  max operator identity err  = {max_op_err:.3e} (rel {max_op_rel_err:.3e})  [u_xx + y*u_yy + 0.5*u_y == 4*U_{{eta,xi}}]")
    print(f"  max general sol PDE resid  = {max_sol_err:.3e}")


if __name__ == "__main__":
    main()
