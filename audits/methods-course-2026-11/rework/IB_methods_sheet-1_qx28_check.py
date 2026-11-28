"""Computational verification for IB/methods/sheet-1/qx28."""
import math


def simpson(f, a, b, n=2000):
    assert n % 2 == 0
    h = (b - a) / n
    s = f(a) + f(b)
    for k in range(1, n):
        s += (4 if k % 2 == 1 else 2) * f(a + k * h)
    return s * h / 3.0


def test_part_i():
    def y(x):
        return math.exp(0.6 * x) + 0.4 * math.cos(2.3 * x)

    def yp(x):
        return 0.6 * math.exp(0.6 * x) - 0.4 * 2.3 * math.sin(2.3 * x)

    def ypp(x):
        return 0.36 * math.exp(0.6 * x) - 0.4 * 5.29 * math.cos(2.3 * x)

    def p(x):
        return 1.0 - x * x

    h = 1e-6
    max_err = 0.0
    for n in (0, 1, 2.5, 4):
        lam = n * (n + 1.0)
        for i in range(1, 20):
            x = -1.0 + 2.0 * i / 20.0
            sl = (p(x + h) * yp(x + h) - p(x - h) * yp(x - h)) / (2.0 * h) + lam * y(x)
            orig = (1.0 - x * x) * ypp(x) - 2.0 * x * yp(x) + lam * y(x)
            max_err = max(max_err, abs(sl - orig))
    assert max_err < 1e-8
    print(f"Part (i) Legendre SL check passed (max_err = {max_err:.2e}).")


def test_part_ii():
    def y(x):
        return math.sin(1.4 * x) + 0.3 * math.exp(0.5 * x)

    def yp(x):
        return 1.4 * math.cos(1.4 * x) + 0.15 * math.exp(0.5 * x)

    def ypp(x):
        return -1.96 * math.sin(1.4 * x) + 0.075 * math.exp(0.5 * x)

    h = 1e-6
    triples = [(1.7, -0.4, 0.6), (2.0, 3.0, 1.5), (-1.2, 2.8, 0.9), (0.5, 0.5, 1.0)]
    max_err_01 = 0.0
    max_err_gt1 = 0.0
    max_err_lt0 = 0.0

    for a, b, c in triples:
        # 1) Interval (0, 1): (x^c (1-x)^(1+a+b-c) y')' - ab x^(c-1) (1-x)^(a+b-c) y
        #    = - x^(c-1) (1-x)^(a+b-c) [x(x-1)y'' + ((1+a+b)x - c)y' + ab y]
        for i in range(1, 20):
            x = i / 20.0

            def p01(t):
                return (t ** c) * ((1.0 - t) ** (1.0 + a + b - c))

            w01 = (x ** (c - 1.0)) * ((1.0 - x) ** (a + b - c))
            sl = (p01(x + h) * yp(x + h) - p01(x - h) * yp(x - h)) / (2.0 * h) - a * b * w01 * y(x)
            orig = -w01 * (x * (x - 1.0) * ypp(x) + ((1.0 + a + b) * x - c) * yp(x) + a * b * y(x))
            max_err_01 = max(max_err_01, abs(sl - orig))

        # 2) Interval (1, inf): (x^c (x-1)^(1+a+b-c) y')' + ab x^(c-1) (x-1)^(a+b-c) y
        #    = + x^(c-1) (x-1)^(a+b-c) [x(x-1)y'' + ((1+a+b)x - c)y' + ab y]
        for i in range(1, 15):
            x = 1.0 + 0.15 * i

            def pgt1(t):
                return (t ** c) * ((t - 1.0) ** (1.0 + a + b - c))

            wgt1 = (x ** (c - 1.0)) * ((x - 1.0) ** (a + b - c))
            sl = (pgt1(x + h) * yp(x + h) - pgt1(x - h) * yp(x - h)) / (2.0 * h) + a * b * wgt1 * y(x)
            orig = wgt1 * (x * (x - 1.0) * ypp(x) + ((1.0 + a + b) * x - c) * yp(x) + a * b * y(x))
            max_err_gt1 = max(max_err_gt1, abs(sl - orig))

        # 3) Interval (-inf, 0): ((-x)^c (1-x)^(1+a+b-c) y')' + ab (-x)^(c-1) (1-x)^(a+b-c) y
        #    = + (-x)^(c-1) (1-x)^(a+b-c) [x(x-1)y'' + ((1+a+b)x - c)y' + ab y]
        for i in range(1, 15):
            x = -0.15 * i

            def plt0(t):
                return ((-t) ** c) * ((1.0 - t) ** (1.0 + a + b - c))

            wlt0 = ((-x) ** (c - 1.0)) * ((1.0 - x) ** (a + b - c))
            sl = (plt0(x + h) * yp(x + h) - plt0(x - h) * yp(x - h)) / (2.0 * h) + a * b * wlt0 * y(x)
            orig = wlt0 * (x * (x - 1.0) * ypp(x) + ((1.0 + a + b) * x - c) * yp(x) + a * b * y(x))
            max_err_lt0 = max(max_err_lt0, abs(sl - orig))

    assert max_err_01 < 1e-6 and max_err_gt1 < 1e-5 and max_err_lt0 < 1e-5

    # Bonus: Verify orthogonality of hypergeometric polynomials 2F1(-n, n+gamma; c; x) on (0,1)
    # with fixed gamma = a + b = 3.0, c = 1.5 -> w(x) = x^(0.5) (1-x)^(1.5)
    gamma = 3.0
    c = 1.5

    def hyp2f1_poly(n, x):
        # a = -n, b = n + gamma
        a, b = -n, n + gamma
        s = 1.0
        term = 1.0
        for k in range(n):
            term *= (a + k) * (b + k) / ((c + k) * (k + 1.0)) * x
            s += term
        return s

    max_poly_orth = 0.0
    for n in range(4):
        for m in range(n + 1, 4):
            # Substitute x = sin^2(t), dx = 2 sin(t) cos(t) dt -> x^(0.5) (1-x)^(1.5) dx = 2 sin^2(t) cos^4(t) dt
            val = simpson(
                lambda t, n=n, m=m: 2.0
                * (math.sin(t) ** 2)
                * (math.cos(t) ** 4)
                * hyp2f1_poly(n, math.sin(t) ** 2)
                * hyp2f1_poly(m, math.sin(t) ** 2),
                0.0,
                0.5 * math.pi,
                n=1000,
            )
            max_poly_orth = max(max_poly_orth, abs(val))
    assert max_poly_orth < 1e-13
    print(
        f"Part (ii) SL checks passed on (0,1) ({max_err_01:.2e}), (1,inf) ({max_err_gt1:.2e}), "
        f"(-inf,0) ({max_err_lt0:.2e}); 2F1 polynomial orthogonality err={max_poly_orth:.2e}."
    )


def test_part_iii():
    max_ode_err = 0.0
    max_orth_err = 0.0
    for n in range(1, 7):
        lam = (n * math.pi) ** 2
        for i in range(1, 20):
            x = i / 20.0
            k = n * math.pi
            y = math.exp(-2.0 * x) * math.sin(k * x)
            yp = math.exp(-2.0 * x) * (k * math.cos(k * x) - 2.0 * math.sin(k * x))
            ypp = math.exp(-2.0 * x) * ((4.0 - k * k) * math.sin(k * x) - 4.0 * k * math.cos(k * x))
            max_ode_err = max(max_ode_err, abs(ypp + 4.0 * yp + (4.0 + lam) * y))
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

    assert max_ode_err < 1e-11 and max_orth_err < 1e-11
    print(f"Part (iii) checks passed: max_ode_err={max_ode_err:.2e}, max_orth_err={max_orth_err:.2e}.")


if __name__ == "__main__":
    test_part_i()
    test_part_ii()
    test_part_iii()
    print("ALL CHECKS PASSED FOR IB/methods/sheet-1/qx28.")
