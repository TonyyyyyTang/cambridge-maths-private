"""Computational verification for IB/methods/sheet-1/q08."""
import math


def simpson(f, a, b, n=2000):
    assert n % 2 == 0
    h = (b - a) / n
    s = f(a) + f(b)
    for k in range(1, n):
        s += (4 if k % 2 == 1 else 2) * f(a + k * h)
    return s * h / 3.0


def test_part_i():
    # Test function y(x) = exp(0.7*x) + 0.3*sin(2.1*x)
    def y(x):
        return math.exp(0.7 * x) + 0.3 * math.sin(2.1 * x)

    def yp(x):
        return 0.7 * math.exp(0.7 * x) + 0.3 * 2.1 * math.cos(2.1 * x)

    def ypp(x):
        return 0.49 * math.exp(0.7 * x) - 0.3 * 4.41 * math.sin(2.1 * x)

    def p(x):
        return 1.0 - x * x

    # Verify p(-1) = p(1) = 0 (singular at both endpoints)
    assert p(-1.0) == 0.0 and p(1.0) == 0.0

    # Verify (p y')' + lambda y == (1-x^2)y'' - 2x y' + lambda y
    h = 1e-6
    max_err = 0.0
    for lam in (-2.5, 0.0, 6.0):
        for i in range(1, 20):
            x = -1.0 + 2.0 * i / 20.0
            sl_deriv = (p(x + h) * yp(x + h) - p(x - h) * yp(x - h)) / (2.0 * h)
            lhs_sl = sl_deriv + lam * y(x)
            lhs_orig = (1.0 - x * x) * ypp(x) - 2.0 * x * yp(x) + lam * y(x)
            max_err = max(max_err, abs(lhs_sl - lhs_orig))
    assert max_err < 1e-8, f"Part (i) mismatch: {max_err}"
    print(f"Part (i) SL form check passed (max_err = {max_err:.2e}), p(-1)=p(1)=0.")


def test_part_ii():
    def y(x):
        return math.cos(1.7 * x) + 0.5 * math.exp(0.9 * x)

    def yp(x):
        return -1.7 * math.sin(1.7 * x) + 0.45 * math.exp(0.9 * x)

    def ypp(x):
        return -2.89 * math.cos(1.7 * x) + 0.405 * math.exp(0.9 * x)

    params = [(2.5, 0.4), (3.0, 1.0), (4.2, 2.7), (1.3, 0.8), (2.0, 1.0)]
    h = 1e-6
    max_err = 0.0
    for a, b in params:
        assert a > b > 0

        def p(x):
            return (x ** b) * ((1.0 - x) ** (a - b))

        def w(x):
            return (x ** (b - 1.0)) * ((1.0 - x) ** (a - b - 1.0))

        # Singular endpoints: p(0) = p(1) = 0
        assert p(0.0) == 0.0 and p(1.0) == 0.0

        for lam in (-1.3, 0.0, 5.2):
            for i in range(1, 20):
                x = i / 20.0
                sl_deriv = (p(x + h) * yp(x + h) - p(x - h) * yp(x - h)) / (2.0 * h)
                lhs_sl = sl_deriv + lam * w(x) * y(x)
                lhs_orig = w(x) * (x * (1.0 - x) * ypp(x) - (a * x - b) * yp(x) + lam * y(x))
                max_err = max(max_err, abs(lhs_sl - lhs_orig))

        # Check weight integrability via smooth endpoint-desingularising substitutions:
        # On [0, 0.5], set x = u^m1 with m1*b >= 5 -> x^(b-1) dx = m1 u^(m1*b - 1) du (smooth at u=0)
        # On [0.5, 1], set 1-x = v^m2 with m2*(a-b) >= 5 -> (1-x)^(a-b-1) dx = -m2 v^(m2*(a-b)-1) dv
        beta_exact = math.gamma(b) * math.gamma(a - b) / math.gamma(a)
        m1 = max(5, math.ceil(5.0 / b))
        m2 = max(5, math.ceil(5.0 / (a - b)))
        int_left = simpson(
            lambda u: m1 * (u ** (m1 * b - 1.0)) * ((1.0 - u ** m1) ** (a - b - 1.0)),
            0.0,
            0.5 ** (1.0 / m1),
            n=2000,
        )
        int_right = simpson(
            lambda v: m2 * (v ** (m2 * (a - b) - 1.0)) * ((1.0 - v ** m2) ** (b - 1.0)),
            0.0,
            0.5 ** (1.0 / m2),
            n=2000,
        )
        err_beta = abs((int_left + int_right) - beta_exact)
        assert err_beta < 1e-10, f"Beta mismatch for (a={a}, b={b}): err={err_beta:.2e}"

    assert max_err < 1e-7, f"Part (ii) mismatch: {max_err}"
    print(f"Part (ii) SL form & Beta(b,a-b) weight checks passed for {len(params)} (a,b) pairs (max_err = {max_err:.2e}).")


def test_part_iii():
    # Verify SL identity (e^{4x} y')' + (4 + lam) e^{4x} y = e^{4x} (y'' + 4y' + (4+lam)y)
    # and verify eigenvalues lam_n = (n*pi)^2 and eigenfunctions y_n(x) = e^{-2x} sin(n*pi*x).
    max_ode_err = 0.0
    max_bc_err = 0.0
    for n in range(1, 7):
        lam = (n * math.pi) ** 2

        def yn(x, n=n):
            return math.exp(-2.0 * x) * math.sin(n * math.pi * x)

        def ynp(x, n=n):
            return math.exp(-2.0 * x) * (n * math.pi * math.cos(n * math.pi * x) - 2.0 * math.sin(n * math.pi * x))

        def ynpp(x, n=n):
            k = n * math.pi
            return math.exp(-2.0 * x) * ((4.0 - k * k) * math.sin(k * x) - 4.0 * k * math.cos(k * x))

        max_bc_err = max(max_bc_err, abs(yn(0.0)), abs(yn(1.0)))
        for i in range(1, 20):
            x = i / 20.0
            res = ynpp(x) + 4.0 * ynp(x) + (4.0 + lam) * yn(x)
            max_ode_err = max(max_ode_err, abs(res))

    assert max_bc_err < 1e-14 and max_ode_err < 1e-11

    # Verify orthogonality relation: int_0^1 e^{4x} y_n(x) y_m(x) dx = (1/2) delta_{nm}
    max_orth_err = 0.0
    for n in range(1, 7):
        for m in range(1, 7):
            val = simpson(
                lambda x, n=n, m=m: math.exp(4.0 * x)
                * (math.exp(-2.0 * x) * math.sin(n * math.pi * x))
                * (math.exp(-2.0 * x) * math.sin(m * math.pi * x)),
                0.0,
                1.0,
                n=2000,
            )
            expected = 0.5 if n == m else 0.0
            max_orth_err = max(max_orth_err, abs(val - expected))
    assert max_orth_err < 1e-11

    # Shooting check: solve y'' + 4y' + (4+lam)y = 0, y(0)=0, y'(0)=1 via exact formula & RK4
    def shoot_rk4(lam, steps=400):
        h = 1.0 / steps
        y, v = 0.0, 1.0
        for _ in range(steps):
            def f(yy, vv):
                return vv, -4.0 * vv - (4.0 + lam) * yy
            k1y, k1v = f(y, v)
            k2y, k2v = f(y + 0.5 * h * k1y, v + 0.5 * h * k1v)
            k3y, k3v = f(y + 0.5 * h * k2y, v + 0.5 * h * k2v)
            k4y, k4v = f(y + h * k3y, v + h * k3v)
            y += (h / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y)
            v += (h / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
        return y

    # Verify no zero of y(1; lam) for lam <= 0
    for lam in (-20.0, -5.0, -1.0, -0.01, 0.0):
        y1 = shoot_rk4(lam)
        assert y1 > 0.1, f"Unexpected small y(1) at lam={lam}: {y1}"

    # Verify zeros of y(1; lam) at lam = (n*pi)^2 for n=1..4
    for n in range(1, 5):
        y1 = shoot_rk4((n * math.pi) ** 2, steps=800)
        assert abs(y1) < 1e-8, f"Shooting error at n={n}: {y1}"

    print(
        f"Part (iii) checks passed: max_ode_err={max_ode_err:.2e}, "
        f"max_bc_err={max_bc_err:.2e}, max_orth_err={max_orth_err:.2e}, RK4 shooting verified."
    )


if __name__ == "__main__":
    test_part_i()
    test_part_ii()
    test_part_iii()
    print("ALL CHECKS PASSED FOR IB/methods/sheet-1/q08.")
