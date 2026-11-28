"""Computational verification for IB/methods/sheet-1/qx16 (2001-02 Sheet 1 Q10)."""

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


def delta_gaussian(u, eps):
    """Symmetric Gaussian nascent delta function."""
    s = u / eps
    if abs(s) > 40.0:
        return 0.0
    return math.exp(-0.5 * s * s) / (math.sqrt(2.0 * math.pi) * eps)


def delta_asymmetric(u, eps):
    """Strongly asymmetric nascent delta function with total integral 1 and int_0^inf > 0.7."""
    s = u / eps
    val = 0.0
    s1 = s - 1.0
    if abs(s1) <= 40.0:
        val += (0.75 / math.sqrt(math.pi)) * math.exp(-s1 * s1)
    s2 = s + 0.5
    if abs(s2) <= 40.0:
        val += (0.25 / math.sqrt(math.pi)) * math.exp(-s2 * s2)
    return val / eps


def verify_delta_kernels():
    for name, kernel in [("Gaussian", delta_gaussian), ("Asymmetric", delta_asymmetric)]:
        total = simpson(lambda u: kernel(u, 0.2), -5.0, 5.0, 4000)
        assert abs(total - 1.0) < 1e-12, f"{name} integral = {total}"
    print("Both symmetric and asymmetric nascent delta kernels integrate to 1.")


def integrate_near_zeros(integrand, a, b, zeros, width, n_local=400):
    """Integrate accurately over [a, b] by refining tightly around each zero."""
    pts = [a, b]
    for z in zeros:
        if a < z - width:
            pts.append(z - width)
        if z + width < b:
            pts.append(z + width)
    pts = sorted(set(pts))
    total = 0.0
    for left, right in zip(pts[:-1], pts[1:]):
        total += simpson(integrand, left, right, n_local)
    return total


def verify_monotone_and_general():
    eps = 1e-5

    # 1. Monotone increasing: phi(x) = exp(x) - 2 on [0, 2], zero c = ln(2), phi'(c) = 2
    phi_inc = lambda x: math.exp(x) - 2.0
    f_inc = lambda x: x * x + math.cos(x)
    c_inc = math.log(2.0)
    expected_inc = f_inc(c_inc) / abs(math.exp(c_inc))

    for kname, kernel in [("Gaussian", delta_gaussian), ("Asymmetric", delta_asymmetric)]:
        val = integrate_near_zeros(
            lambda x: f_inc(x) * kernel(phi_inc(x), eps),
            0.0, 2.0, [c_inc], width=20.0 * eps
        )
        assert abs(val - expected_inc) < 1e-5, f"Inc {kname}: {val} vs {expected_inc}"
        print(f"Monotone increasing ({kname}, eps={eps}): num={val:.10f}, exact={expected_inc:.10f}, err={abs(val-expected_inc):.2e}")

    # 2. Monotone decreasing: phi(x) = 1 - x^3 on [0.2, 2.0], zero c = 1, phi'(c) = -3
    phi_dec = lambda x: 1.0 - x ** 3
    f_dec = lambda x: math.exp(0.5 * x)
    c_dec = 1.0
    expected_dec = f_dec(c_dec) / abs(-3.0)

    for kname, kernel in [("Gaussian", delta_gaussian), ("Asymmetric", delta_asymmetric)]:
        val = integrate_near_zeros(
            lambda x: f_dec(x) * kernel(phi_dec(x), eps),
            0.2, 2.0, [c_dec], width=20.0 * eps
        )
        assert abs(val - expected_dec) < 1e-5, f"Dec {kname}: {val} vs {expected_dec}"
        print(f"Monotone decreasing ({kname}, eps={eps}): num={val:.10f}, exact={expected_dec:.10f}, err={abs(val-expected_dec):.2e}")

    # 3. General phi(x) with multiple simple zeros of alternating slope signs:
    # phi(x) = x^3 - x on [-2, 2], zeros at -1, 0, 1; phi'(x) = 3x^2 - 1 -> 2, -1, 2
    phi_gen = lambda x: x ** 3 - x
    f_gen = lambda x: 1.0 + 0.7 * x + 0.3 * math.sin(x)
    zeros_gen = [-1.0, 0.0, 1.0]
    expected_gen = sum(f_gen(z) / abs(3.0 * z * z - 1.0) for z in zeros_gen)

    for kname, kernel in [("Gaussian", delta_gaussian), ("Asymmetric", delta_asymmetric)]:
        val = integrate_near_zeros(
            lambda x: f_gen(x) * kernel(phi_gen(x), eps),
            -2.0, 2.0, zeros_gen, width=25.0 * eps
        )
        assert abs(val - expected_gen) < 1e-4, f"Gen {kname}: {val} vs {expected_gen}"
        print(f"General multi-zero ({kname}, eps={eps}): num={val:.10f}, exact={expected_gen:.10f}, err={abs(val-expected_gen):.2e}")


def verify_target_integral():
    eps = 1e-4
    # Target: int_{-infty}^{+infty} |x| delta(x^2 - a^2) dx = 1 for a != 0 (both a > 0 and a < 0)
    for a in [1.7, -0.8, 0.35]:
        R = abs(a) + 2.0
        zeros = [-abs(a), abs(a)]
        for kname, kernel in [("Gaussian", delta_gaussian), ("Asymmetric", delta_asymmetric)]:
            val = integrate_near_zeros(
                lambda x, a=a: abs(x) * kernel(x * x - a * a, eps),
                -R, R, zeros, width=30.0 * eps / (2.0 * abs(a))
            )
            assert abs(val - 1.0) < 1e-5, f"Target a={a}, {kname}: val={val}"
            print(f"Target integral a={a:+.2f} ({kname}, eps={eps}): val={val:.10f}, err={abs(val-1.0):.2e}")

    # Degenerate a = 0 check: int |x| delta_eps(x^2) dx = int_0^inf eta(s) ds (kernel-dependent!)
    eps0 = 1e-3
    w0 = math.sqrt(30.0 * eps0)
    val_gauss_0 = simpson(lambda x: abs(x) * delta_gaussian(x * x, eps0), -w0, w0, 4000)
    val_asym_0 = simpson(lambda x: abs(x) * delta_asymmetric(x * x, eps0), -w0, w0, 4000)
    # Exact int_0^inf eta_asym(s) ds = 0.75*(1 + erf(1.0))/2 + 0.25*(1 - erf(0.5))/2
    expected_asym_0 = 0.375 * (1.0 + math.erf(1.0)) + 0.125 * (1.0 - math.erf(0.5))
    assert abs(val_gauss_0 - 0.5) < 1e-8
    assert abs(val_asym_0 - expected_asym_0) < 1e-8
    assert abs(val_gauss_0 - val_asym_0) > 0.2
    print(
        f"Degenerate a=0 check: Gaussian -> {val_gauss_0:.6f} (=0.5), "
        f"Asymmetric -> {val_asym_0:.6f} (={expected_asym_0:.6f}) => ill-defined at a=0!"
    )


if __name__ == "__main__":
    verify_delta_kernels()
    verify_monotone_and_general()
    verify_target_integral()
    print("ALL CHECKS PASSED for IB/methods/sheet-1/qx16.")
