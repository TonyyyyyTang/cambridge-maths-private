"""Independent numerical verification for IB/2018/P4/q5A.

Verifies the separated-variables solution
    u(x, y) = 2 * sin(3*pi*x) * sinh(3*pi*y) / sinh(3*pi)
to Laplace's equation u_xx + u_yy = 0 on (0,1)x(0,1) with
    u(0,y) = u(1,y) = u(x,0) = 0,  u(x,1) = 2*sin(3*pi*x).
"""

import math


def u_exact(x: float, y: float) -> float:
    return 2.0 * math.sin(3.0 * math.pi * x) * math.sinh(3.0 * math.pi * y) / math.sinh(3.0 * math.pi)


def simpson(f, a: float, b: float, n: int) -> float:
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * f(a + i * h)
    return s * h / 3.0


def main() -> None:
    # 1. Check 4th-order finite-difference Laplacian on an interior grid
    h = 1e-3
    max_lap_err = 0.0
    for ix in range(1, 10):
        x = ix / 10.0
        for iy in range(1, 10):
            y = iy / 10.0
            uxx = (-u_exact(x - 2 * h, y) + 16 * u_exact(x - h, y) - 30 * u_exact(x, y)
                   + 16 * u_exact(x + h, y) - u_exact(x + 2 * h, y)) / (12 * h * h)
            uyy = (-u_exact(x, y - 2 * h) + 16 * u_exact(x, y - h) - 30 * u_exact(x, y)
                   + 16 * u_exact(x, y + h) - u_exact(x, y + 2 * h)) / (12 * h * h)
            max_lap_err = max(max_lap_err, abs(uxx + uyy))
    print(f"Max 4th-order FD Laplacian residual: {max_lap_err:.3e}")
    assert max_lap_err < 1e-8

    # 2. Check all four boundary conditions on a fine grid
    max_bc_err = 0.0
    for i in range(101):
        s = i / 100.0
        max_bc_err = max(max_bc_err, abs(u_exact(0.0, s)))
        max_bc_err = max(max_bc_err, abs(u_exact(1.0, s)))
        max_bc_err = max(max_bc_err, abs(u_exact(s, 0.0)))
        max_bc_err = max(max_bc_err, abs(u_exact(s, 1.0) - 2.0 * math.sin(3.0 * math.pi * s)))
    print(f"Max boundary condition error:        {max_bc_err:.3e}")
    assert max_bc_err < 1e-14

    # 3. Check Fourier sine mode projection at y = 1
    max_mode_err = 0.0
    for m in range(1, 11):
        proj = 2.0 * simpson(lambda x: 2.0 * math.sin(3.0 * math.pi * x) * math.sin(m * math.pi * x), 0.0, 1.0, 400)
        expected = 2.0 if m == 3 else 0.0
        max_mode_err = max(max_mode_err, abs(proj - expected))
    print(f"Max Fourier sine coefficient error:  {max_mode_err:.3e}")
    assert max_mode_err < 1e-12

    # 4. Compare against independent 2D finite-difference Laplace solver (SOR)
    #    on two grid sizes N=20 and N=40 to confirm O(h^2) convergence (~4x reduction)
    fd_errors = []
    for N in (20, 40):
        grid = [[0.0] * (N + 1) for _ in range(N + 1)]
        for i in range(N + 1):
            grid[i][N] = 2.0 * math.sin(3.0 * math.pi * i / N)
        omega = 2.0 / (1.0 + math.sin(math.pi / N))
        for _ in range(5000):
            max_diff = 0.0
            for i in range(1, N):
                for j in range(1, N):
                    val = 0.25 * (grid[i - 1][j] + grid[i + 1][j] + grid[i][j - 1] + grid[i][j + 1])
                    new_val = (1.0 - omega) * grid[i][j] + omega * val
                    diff = abs(new_val - grid[i][j])
                    if diff > max_diff:
                        max_diff = diff
                    grid[i][j] = new_val
            if max_diff < 1e-13:
                break
        err = max(
            abs(grid[i][j] - u_exact(i / N, j / N))
            for i in range(1, N) for j in range(1, N)
        )
        fd_errors.append(err)
        print(f"2D FD Laplace solver (N={N:2d}) max error vs exact: {err:.3e}")

    ratio = fd_errors[0] / fd_errors[1]
    print(f"Error reduction ratio (N=20 -> N=40, expected ~4): {ratio:.3f}")
    assert fd_errors[1] < 2e-2
    assert 3.8 < ratio < 4.3
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
