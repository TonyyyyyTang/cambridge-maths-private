#!/usr/bin/env python3
"""Independent verification script for IB/2018/P2/q16A (Methods)."""

from fractions import Fraction
import math


def simpson_real(func, a, b, n):
    """Composite Simpson's rule (n even)."""
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    total = func(a) + func(b)
    for j in range(1, n):
        x = a + j * h
        total += (4 if j % 2 == 1 else 2) * func(x)
    return total * (h / 3.0)


def f_periodic(x):
    """2pi-periodic extension of |x| on [-pi, pi]."""
    t = (x + math.pi) % (2.0 * math.pi) - math.pi
    return abs(t)


def check_part_a():
    """Verify Fourier coefficients of f(x) = |x| on [-pi, pi]."""
    # Integrate on [-pi, 0] and [0, pi] separately so kink at 0 is on a panel boundary
    a0_num = (
        simpson_real(lambda x: -x, -math.pi, 0.0, 200)
        + simpson_real(lambda x: x, 0.0, math.pi, 200)
    ) / math.pi
    assert abs(a0_num - math.pi) < 1e-13

    max_coeff_err = 0.0
    for n in range(1, 15):
        an_num = (
            simpson_real(lambda x, n=n: -x * math.cos(n * x), -math.pi, 0.0, 4000)
            + simpson_real(lambda x, n=n: x * math.cos(n * x), 0.0, math.pi, 4000)
        ) / math.pi
        bn_num = (
            simpson_real(lambda x, n=n: -x * math.sin(n * x), -math.pi, 0.0, 4000)
            + simpson_real(lambda x, n=n: x * math.sin(n * x), 0.0, math.pi, 4000)
        ) / math.pi

        an_exact = 0.0 if n % 2 == 0 else -4.0 / (math.pi * n * n)
        bn_exact = 0.0
        err_a = abs(an_num - an_exact)
        err_b = abs(bn_num - bn_exact)
        max_coeff_err = max(max_coeff_err, err_a, err_b)
        assert err_a < 1e-9 and err_b < 1e-9, (
            f"Coeff mismatch at n={n}: err_a={err_a}, err_b={err_b}"
        )

    # Verify Fourier series convergence to |x|
    def f_series(x, n_max=500):
        s = 0.5 * math.pi
        for n in range(1, n_max + 1, 2):
            s -= (4.0 / (math.pi * n * n)) * math.cos(n * x)
        return s

    for x in [-math.pi, -2.1, -0.5, 0.0, 0.7, 1.9, math.pi]:
        assert abs(f_series(x, 2000) - f_periodic(x)) < 1e-3

    print(f"[PASS] Part (a): Fourier coefficients verified (max coeff err = {max_coeff_err:.3e})")


def y_fourier(x, n_max=2001):
    """Fourier series solution of y'' + 2y' + y = f(x) (2pi-periodic)."""
    s = 0.5 * math.pi
    for n in range(1, n_max + 1, 2):
        denom = (n * n) * ((n * n + 1) ** 2)
        s += (4.0 / math.pi) * (
            ((n * n - 1) * math.cos(n * x) - 2.0 * n * math.sin(n * x)) / denom
        )
    return s


def y_closed_form(x):
    """Independent exact closed-form 2pi-periodic piecewise solution."""
    t = (x + math.pi) % (2.0 * math.pi) - math.pi
    ep = math.exp(math.pi)
    d2 = -2.0 / (ep + 1.0)
    c2 = 2.0 * ep / (ep + 1.0)
    d1 = -4.0 / (ep + 1.0) - (2.0 * math.pi * ep) / ((ep + 1.0) ** 2)
    c1 = d1 + 4.0
    if t >= 0.0:
        return t - 2.0 + (c1 + c2 * t) * math.exp(-t)
    else:
        return -t + 2.0 + (d1 + d2 * t) * math.exp(-t)


def check_part_b():
    """Verify exact linear algebra, Fourier series vs closed form, ODE, BCs, and uniqueness."""
    # 1. Exact rational check of the 2x2 Fourier mode linear system (factoring out 1/pi)
    for n in range(1, 25):
        a_n_scaled = Fraction(0) if n % 2 == 0 else Fraction(-4, n * n)
        A_n_scaled = a_n_scaled * Fraction(1 - n * n, (n * n + 1) ** 2)
        B_n_scaled = a_n_scaled * Fraction(2 * n, (n * n + 1) ** 2)
        assert (1 - n * n) * A_n_scaled + 2 * n * B_n_scaled == a_n_scaled
        assert -2 * n * A_n_scaled + (1 - n * n) * B_n_scaled == 0
        if n % 2 == 1:
            assert A_n_scaled == Fraction(4 * (n * n - 1), n * n * (n * n + 1) ** 2)
            assert B_n_scaled == Fraction(-8, n * (n * n + 1) ** 2)

    # 2. Cross-check Fourier series against independent closed-form piecewise solution
    max_diff = 0.0
    for k in range(-25, 26):
        x = (k / 25.0) * math.pi
        diff = abs(y_fourier(x, 2001) - y_closed_form(x))
        max_diff = max(max_diff, diff)
        assert diff < 1e-9, f"FS vs closed form mismatch at x={x}: diff={diff}"

    # 3. Verify closed-form solution satisfies ODE, C^2 smoothness at 0, and 2pi-periodicity
    h = 1e-4
    max_ode_err = 0.0
    for x in [-2.8, -1.7, -0.4, 0.3, 1.5, 2.6, 3.8, -4.1]:
        yp = (y_closed_form(x + h) - y_closed_form(x - h)) / (2.0 * h)
        ypp = (
            y_closed_form(x + h) - 2.0 * y_closed_form(x) + y_closed_form(x - h)
        ) / (h ** 2)
        res = (ypp + 2.0 * yp + y_closed_form(x)) - f_periodic(x)
        max_ode_err = max(max_ode_err, abs(res))
        assert abs(res) < 1e-6, f"ODE residual at x={x}: {res}"

    # Check C^2 continuity across x = 0 and x = pi (periodicity)
    for x0 in [0.0, math.pi]:
        val_l = y_closed_form(x0 - 1e-12)
        val_r = y_closed_form(x0 + 1e-12)
        assert abs(val_r - val_l) < 1e-11
        d1_l = (y_closed_form(x0 - h) - y_closed_form(x0 - 2 * h)) / h
        d1_r = (y_closed_form(x0 + 2 * h) - y_closed_form(x0 + h)) / h
        assert abs(d1_r - d1_l) < 2e-3
        d2_l = (
            y_closed_form(x0 - h)
            - 2.0 * y_closed_form(x0 - 2 * h)
            + y_closed_form(x0 - 3 * h)
        ) / (h ** 2)
        d2_r = (
            y_closed_form(x0 + 3 * h)
            - 2.0 * y_closed_form(x0 + 2 * h)
            + y_closed_form(x0 + h)
        ) / (h ** 2)
        assert abs(d2_r - d2_l) < 2e-3

    # 4. Verify ODE directly on the Fourier series too
    max_fs_ode_err = 0.0
    for x in [-2.3, -0.8, 0.6, 2.1]:
        yp = (y_fourier(x + h, 1501) - y_fourier(x - h, 1501)) / (2.0 * h)
        ypp = (
            y_fourier(x + h, 1501)
            - 2.0 * y_fourier(x, 1501)
            + y_fourier(x - h, 1501)
        ) / (h ** 2)
        res = (ypp + 2.0 * yp + y_fourier(x, 1501)) - f_periodic(x)
        max_fs_ode_err = max(max_fs_ode_err, abs(res))
        assert abs(res) < 1e-3, f"FS ODE residual at x={x}: {res}"

    # 5. Verify no non-trivial 2pi-periodic homogeneous solution (A + Bx)e^{-x}
    for A, B in [(1.0, 0.0), (0.0, 1.0), (1.0, -2.0)]:
        yh = lambda x, A=A, B=B: (A + B * x) * math.exp(-x)
        assert max(abs(yh(2.0 * math.pi) - yh(0.0)), abs(yh(-2.0 * math.pi) - yh(0.0))) > 1.0

    print(
        f"[PASS] Part (b): FS vs closed-form max diff = {max_diff:.3e}, "
        f"closed-form ODE residual = {max_ode_err:.3e}, "
        f"FS ODE residual = {max_fs_ode_err:.3e}, uniqueness verified"
    )


if __name__ == "__main__":
    check_part_a()
    check_part_b()
    print("ALL CHECKS PASSED FOR IB_2018_P2_q16A")
