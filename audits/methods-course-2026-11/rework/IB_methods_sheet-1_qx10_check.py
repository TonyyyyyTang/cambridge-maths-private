"""Computational verification for IB/methods/sheet-1/qx10 (2001-02 Sheet 1 Q1)."""

from fractions import Fraction
import math


def simpson(func, a, b, n):
    """Composite Simpson's rule on [a, b] with even integer n subintervals."""
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    total = func(a) + func(b)
    for k in range(1, n):
        weight = 4.0 if k % 2 == 1 else 2.0
        total += weight * func(a + k * h)
    return total * h / 3.0


def verify_part_i():
    print("=== Part (i): f(x) = (x^2 - 1)^2 on [-1, 1) ===")
    f = lambda x: (x * x - 1.0) ** 2
    fp = lambda x: 4.0 * x * (x * x - 1.0)
    fpp = lambda x: 12.0 * x * x - 4.0
    fppp = lambda x: 24.0 * x

    # Exact rational check for a_0 / 2 = int_0^1 (x^4 - 2x^2 + 1) dx
    a0_half_exact = Fraction(1, 5) - Fraction(2, 3) + Fraction(1, 1)
    assert a0_half_exact == Fraction(8, 15), f"Got {a0_half_exact}"
    print(f"Exact a_0/2 = {a0_half_exact}")

    # Numerical Fourier coefficients on [-1, 1]
    a0_num = simpson(f, -1.0, 1.0, 4000)
    assert abs(a0_num / 2.0 - 8.0 / 15.0) < 1e-13, f"a0_num/2 err = {abs(a0_num / 2.0 - 8.0 / 15.0)}"

    max_a_err = 0.0
    max_b_err = 0.0
    for n in range(1, 13):
        an_num = simpson(lambda x, n=n: f(x) * math.cos(n * math.pi * x), -1.0, 1.0, 4000)
        bn_num = simpson(lambda x, n=n: f(x) * math.sin(n * math.pi * x), -1.0, 1.0, 4000)
        an_exact = 48.0 * ((-1) ** (n + 1)) / ((math.pi ** 4) * (n ** 4))
        max_a_err = max(max_a_err, abs(an_num - an_exact))
        max_b_err = max(max_b_err, abs(bn_num))
    assert max_a_err < 1e-11, f"max_a_err = {max_a_err}"
    assert max_b_err < 1e-13, f"max_b_err = {max_b_err}"
    print(f"Fourier coeffs n=1..12 verified: max |a_n err| = {max_a_err:.2e}, max |b_n| = {max_b_err:.2e}")

    # Continuity of derivatives at the periodic join x = +-1
    assert f(1.0) == f(-1.0) == 0.0
    assert fp(1.0) == fp(-1.0) == 0.0
    assert fpp(1.0) == fpp(-1.0) == 8.0
    assert fppp(1.0) == 24.0 and fppp(-1.0) == -24.0
    jump_fppp = fppp(-1.0) - fppp(1.0)
    assert jump_fppp == -48.0
    print(f"Join checks at x=+-1: f=0, f'=0, f''=8 (C^2 continuous), jump in f''' = {jump_fppp}")

    # Series convergence to f, f', f'' on [-1, 1] (and 2-periodic extension)
    N = 2000
    def series_f(x):
        s = 8.0 / 15.0
        for n in range(1, N + 1):
            s += (48.0 * ((-1) ** (n + 1)) / ((math.pi ** 4) * (n ** 4))) * math.cos(n * math.pi * x)
        return s

    def series_fp(x):
        s = 0.0
        for n in range(1, N + 1):
            s -= (48.0 * ((-1) ** (n + 1)) / ((math.pi ** 3) * (n ** 3))) * math.sin(n * math.pi * x)
        return s

    def series_fpp(x):
        s = 0.0
        for n in range(1, N + 1):
            s -= (48.0 * ((-1) ** (n + 1)) / ((math.pi ** 2) * (n ** 2))) * math.cos(n * math.pi * x)
        return s

    test_xs = [-1.0, -0.75, -1.0 / math.sqrt(3.0), -0.25, 0.0, 0.37, 1.0 / math.sqrt(3.0), 0.9, 1.0]
    err_f = max(abs(series_f(x) - f(x)) for x in test_xs)
    err_fp = max(abs(series_fp(x) - fp(x)) for x in test_xs)
    err_fpp = max(abs(series_fpp(x) - fpp(x)) for x in test_xs)
    assert err_f < 1e-10, f"err_f = {err_f}"
    assert err_fp < 1e-6, f"err_fp = {err_fp}"
    assert err_fpp < 3e-3, f"err_fpp = {err_fpp}"
    # Check 2-periodicity at x = 2.37 -> 0.37
    assert abs(series_f(2.37) - f(0.37)) < 1e-10
    assert abs(series_f(-1.63) - f(0.37)) < 1e-10
    print(f"Series convergence (N={N}): max err f={err_f:.2e}, f'={err_fp:.2e}, f''={err_fpp:.2e}")


def verify_part_ii():
    print("=== Part (ii): f(x) = x^2 on [0, pi], sine/cosine series & differentiation ===")
    pi = math.pi
    f = lambda x: x * x

    # (a) Sine series coefficients b_n
    def bn_exact(n):
        return (2.0 * pi * ((-1) ** (n + 1)) / n) - (4.0 * (1 - (-1) ** n) / (pi * (n ** 3)))

    # (b) Cosine series coefficients a_0/2 and a_n
    a0_half_exact = (pi ** 2) / 3.0
    def an_exact(n):
        return 4.0 * ((-1) ** n) / (n * n)

    max_err_a_coeff = 0.0
    max_err_b_coeff = 0.0
    a0_num = (2.0 / pi) * simpson(f, 0.0, pi, 2000)
    assert abs(a0_num / 2.0 - a0_half_exact) < 1e-13

    for n in range(1, 16):
        bn_num = (2.0 / pi) * simpson(lambda x, n=n: f(x) * math.sin(n * x), 0.0, pi, 4000)
        an_num = (2.0 / pi) * simpson(lambda x, n=n: f(x) * math.cos(n * x), 0.0, pi, 4000)
        max_err_b_coeff = max(max_err_b_coeff, abs(bn_num - bn_exact(n)))
        max_err_a_coeff = max(max_err_a_coeff, abs(an_num - an_exact(n)))
    assert max_err_b_coeff < 1e-10, f"max_err_b_coeff = {max_err_b_coeff}"
    assert max_err_a_coeff < 1e-10, f"max_err_a_coeff = {max_err_a_coeff}"
    print(f"Half-range coeffs n=1..15 verified: max |b_n err|={max_err_b_coeff:.2e}, max |a_n err|={max_err_a_coeff:.2e}")

    # Periodic extensions on [-6*pi, 6*pi]
    def reduce_to_minus_pi_pi(x):
        return ((x + pi) % (2.0 * pi)) - pi

    def f_odd_ext(x):
        u = reduce_to_minus_pi_pi(x)
        return u * abs(u)

    def f_even_ext(x):
        u = reduce_to_minus_pi_pi(x)
        return u * u

    N = 1000
    def S_a(x):
        return sum(bn_exact(n) * math.sin(n * x) for n in range(1, N + 1))

    def S_b(x):
        return a0_half_exact + sum(an_exact(n) * math.cos(n * x) for n in range(1, N + 1))

    # Check interior points across [-6*pi, 6*pi] away from jumps of f_odd at odd multiples of pi
    xs_wide = [-5.4 * pi, -3.7 * pi, -2.0 * pi, -0.6 * pi, 0.0, 0.4 * pi, 2.5 * pi, 4.3 * pi, 5.8 * pi]
    err_Sa = max(abs(S_a(x) - f_odd_ext(x)) for x in xs_wide)
    err_Sb = max(abs(S_b(x) - f_even_ext(x)) for x in xs_wide + [-5 * pi, -3 * pi, -pi, pi, 3 * pi, 5 * pi])
    assert err_Sa < 2e-2, f"err_Sa = {err_Sa}"
    assert err_Sb < 5e-3, f"err_Sb = {err_Sb}"
    # At odd multiples of pi, S_a evaluates strictly to 0 (up to float sin(n*k*pi) roundoff)
    for k in [-5, -3, -1, 1, 3, 5]:
        assert abs(S_a(k * pi)) < 1e-10, f"S_a({k}*pi) = {S_a(k * pi)}"
        assert abs(S_b(k * pi) - pi * pi) < 5e-3
    print(f"Extensions on [-6pi, 6pi] (N={N}): max err S_a={err_Sa:.2e} (0 at odd pi), max err S_b={err_Sb:.2e}")

    # Fourier series of g(x) = x and h(x) = |x| on (-pi, pi)
    def bn_g_exact(n):
        return 2.0 * ((-1) ** (n + 1)) / n

    a0_h_half = pi / 2.0
    def an_h_exact(n):
        return -2.0 * (1 - (-1) ** n) / (pi * n * n)

    for n in range(1, 12):
        bn_g_num = (1.0 / pi) * simpson(lambda x, n=n: x * math.sin(n * x), -pi, pi, 4000)
        an_h_num = (1.0 / pi) * simpson(lambda x, n=n: abs(x) * math.cos(n * x), -pi, pi, 4000)
        assert abs(bn_g_num - bn_g_exact(n)) < 1e-10
        assert abs(an_h_num - an_h_exact(n)) < 1e-10
    assert abs((1.0 / (2.0 * pi)) * simpson(abs, -pi, pi, 2000) - a0_h_half) < 1e-13

    # Term-by-term derivative of S_b: coeff of sin(nx) is -n * an_exact(n)
    for n in range(1, 50):
        sb_prime_sin_coeff = -n * an_exact(n)
        assert abs(sb_prime_sin_coeff - 2.0 * bn_g_exact(n)) < 1e-14
    print("Term-by-term S_b'(x) coeffs match 2*g(x) identically for n=1..49.")

    # Term-by-term derivative of S_a: coeff of cos(nx) is n * bn_exact(n)
    # Check: (1) terms do NOT tend to 0 (classical divergence), and
    # (2) exact distributional identity: S_a'(x) = 2*h(x) - 2*pi^2 * sum_k delta(x - (2k+1)pi)
    for n in range(1, 50):
        sa_prime_cos_coeff = n * bn_exact(n)
        expected_dist_coeff = 2.0 * an_h_exact(n) - 2.0 * pi * ((-1) ** n)
        assert abs(sa_prime_cos_coeff - expected_dist_coeff) < 1e-13
    # Constant term check: 0 = 2*(pi/2) - 2*pi^2 * (1/(2*pi))
    assert abs(0.0 - (2.0 * a0_h_half - 2.0 * pi * pi * (1.0 / (2.0 * pi)))) < 1e-14

    # Check non-decay of terms of S_a'(x) at x = 0.7 (not a multiple of pi)
    x_test = 0.7
    tail_max = max(abs(n * bn_exact(n) * math.cos(n * x_test)) for n in range(100, 120))
    assert tail_max > 5.0, f"Expected non-vanishing terms, got {tail_max}"
    print(f"S_a'(0.7) terms for n=100..119 reach {tail_max:.3f} (does not tend to 0 -> classically divergent).")

    # Distributional pairing test against smooth 2pi-periodic test functions psi(x):
    # <S_a', psi> = sum_{n=1}^\infty (n b_n) * int_{-pi}^{pi} cos(nx) psi(x) dx
    # should equal int_{-pi}^{pi} 2|x| psi(x) dx - 2*pi^2 * psi(pi).
    test_funcs = [
        ("exp(cos x)", lambda x: math.exp(math.cos(x))),
        ("1/(2 + cos x)", lambda x: 1.0 / (2.0 + math.cos(x))),
        ("cos(3x) + 0.5*sin(2x)", lambda x: math.cos(3.0 * x) + 0.5 * math.sin(2.0 * x)),
    ]
    for name, psi in test_funcs:
        lhs = 0.0
        for n in range(1, 60):
            psi_cos_n = simpson(lambda x, n=n: psi(x) * math.cos(n * x), -pi, pi, 4000)
            lhs += (n * bn_exact(n)) * psi_cos_n
        rhs = simpson(lambda x: 2.0 * abs(x) * psi(x), -pi, pi, 4000) - 2.0 * (pi ** 2) * psi(pi)
        assert abs(lhs - rhs) < 1e-9, f"Distributional test failed for {name}: lhs={lhs}, rhs={rhs}"
        print(f"Distributional pairing <S_a', {name}>: |LHS - RHS| = {abs(lhs - rhs):.2e}")


if __name__ == "__main__":
    verify_part_i()
    verify_part_ii()
    print("ALL CHECKS PASSED for IB/methods/sheet-1/qx10.")
