"""Computational verification for IB/methods/sheet-2/qx25."""
import math


def J0_series(z):
    """Bessel function J_0(z) via its Maclaurin series (accurate for |z| <= 12)."""
    term = 1.0
    total = 1.0
    z2_over_4 = 0.25 * z * z
    for m in range(1, 120):
        term *= -z2_over_4 / (m * m)
        total += term
        if abs(term) <= 1e-16 * max(1.0, abs(total)):
            break
    return total


def J0_prime_series(z):
    """Derivative J_0'(z) = -J_1(z) via termwise derivative of Maclaurin series."""
    term = -0.5 * z
    total = term
    z2_over_4 = 0.25 * z * z
    for m in range(2, 120):
        term *= -z2_over_4 / (m * (m - 1))
        total += term
        if abs(term) <= 1e-16 * max(1.0, abs(total)):
            break
    return total


def J0(z):
    """Cancellation-free evaluation of J_0(z) via Poisson's integral (trapezoidal rule on pi-periodic integrand)."""
    if abs(z) <= 10.0:
        return J0_series(z)
    N = max(128, int(4 * abs(z)) + 64)
    step = math.pi / N
    return sum(math.cos(z * math.sin(k * step)) for k in range(N)) / N


def J0_prime(z):
    """Cancellation-free evaluation of J_0'(z) = -(1/pi) int_0^pi sin(theta) sin(z sin(theta)) dtheta."""
    if abs(z) <= 10.0:
        return J0_prime_series(z)
    N = max(128, int(4 * abs(z)) + 64)
    step = math.pi / N
    return -sum(math.sin(k * step) * math.sin(z * math.sin(k * step)) for k in range(N)) / N


def I0(z):
    """Modified Bessel function I_0(z) = J_0(i*z) for real z."""
    term = 1.0
    total = 1.0
    z2_over_4 = 0.25 * z * z
    for m in range(1, 120):
        term *= z2_over_4 / (m * m)
        total += term
        if abs(term) <= 1e-16 * abs(total):
            break
    return total


def find_J0_zeros(count=20):
    """Locate first `count` positive zeros of J_0(z) by bracketing + Newton."""
    zeros = []
    step = 0.25
    x = 0.1
    prev_val = J0(x)
    while len(zeros) < count:
        x_next = x + step
        val = J0(x_next)
        if prev_val * val <= 0.0:
            # Bisect a little then Newton
            lo, hi = x, x_next
            flo = prev_val
            for _ in range(20):
                mid = 0.5 * (lo + hi)
                fmid = J0(mid)
                if flo * fmid <= 0.0:
                    hi = mid
                else:
                    lo = mid
                    flo = fmid
            z = 0.5 * (lo + hi)
            for _ in range(10):
                z -= J0(z) / J0_prime(z)
            zeros.append(z)
        x = x_next
        prev_val = val
    return zeros


def gauss_legendre_8():
    nodes = [
        -0.96028985649753623,
        -0.79666647741362674,
        -0.52553240991632899,
        -0.18343464249564980,
         0.18343464249564980,
         0.52553240991632899,
         0.79666647741362674,
         0.96028985649753623,
    ]
    weights = [
        0.10122853629037626,
        0.22238103445337447,
        0.31370664587788729,
        0.36268378337836198,
        0.36268378337836198,
        0.31370664587788729,
        0.22238103445337447,
        0.10122853629037626,
    ]
    return nodes, weights


def integrate_gl(func, a, b, n_panels=32):
    nodes, weights = gauss_legendre_8()
    h = (b - a) / n_panels
    total = 0.0
    for i in range(n_panels):
        mid = a + (i + 0.5) * h
        half = 0.5 * h
        s = 0.0
        for xi, wi in zip(nodes, weights):
            s += wi * func(mid + half * xi)
        total += half * s
    return total


def main():
    # Part (i): Verify Bessel ODE, zeros, and absence of lambda <= 0 eigenvalues
    h = 1e-5
    max_ode_err = 0.0
    for z in [0.5, 1.7, 4.2, 9.1]:
        # (z J_0'(z))' + z J_0(z) = 0
        deriv = ((z + h) * J0_prime(z + h) - (z - h) * J0_prime(z - h)) / (2.0 * h)
        err = abs(deriv + z * J0(z))
        max_ode_err = max(max_ode_err, err)
        assert err < 1e-8, f"Bessel ODE error: {err}"

    # Verify I_0(mu) > 1 for mu > 0 (no negative eigenvalues)
    for mu in [0.1, 1.0, 3.5, 10.0]:
        assert I0(mu) > 1.0

    zeros = find_J0_zeros(22)
    # Check standard known values for j_1, j_2, j_3
    assert abs(zeros[0] - 2.4048255576957728) < 1e-12
    assert abs(zeros[1] - 5.5200781102863106) < 1e-12
    assert abs(zeros[2] - 8.6537279129110122) < 1e-12
    for jn in zeros:
        assert abs(J0(jn)) < 1e-14
        assert abs(J0_prime(jn)) > 0.05

    # Part (ii): Verify Lommel identity for arbitrary alpha != beta
    max_lommel_err = 0.0
    for alpha, beta in [(1.3, 2.9), (2.4, 5.1), (3.7, 7.8), (0.8, 4.5), (zeros[0], 3.1)]:
        num_int = integrate_gl(lambda x, a=alpha, b=beta: J0(a * x) * J0(b * x) * x, 0.0, 1.0, 32)
        formula = (beta * J0(alpha) * J0_prime(beta) - alpha * J0(beta) * J0_prime(alpha)) / (alpha * alpha - beta * beta)
        err = abs(num_int - formula)
        max_lommel_err = max(max_lommel_err, err)
        assert err < 1e-13, f"Lommel mismatch for ({alpha},{beta}): {err}"

    # Verify orthogonality and normalisation for first 8 eigenfunctions
    max_orth_err = 0.0
    max_norm_err = 0.0
    for n in range(8):
        for m in range(8):
            jn, jm = zeros[n], zeros[m]
            val = integrate_gl(lambda x, a=jn, b=jm: J0(a * x) * J0(b * x) * x, 0.0, 1.0, 48)
            if n != m:
                max_orth_err = max(max_orth_err, abs(val))
                assert abs(val) < 1e-13
            else:
                expected = 0.5 * (J0_prime(jn) ** 2)
                err = abs(val - expected)
                max_norm_err = max(max_norm_err, err)
                assert err < 1e-13

    # Cross-check Maclaurin series vs Poisson integral representation at z = 5.0, 9.5
    for zc in [5.0, 9.5]:
        N = 128
        step = math.pi / N
        j0_int = sum(math.cos(zc * math.sin(k * step)) for k in range(N)) / N
        j0p_int = -sum(math.sin(k * step) * math.sin(zc * math.sin(k * step)) for k in range(N)) / N
        assert abs(j0_int - J0_series(zc)) < 1e-14
        assert abs(j0p_int - J0_prime_series(zc)) < 1e-14

    # Verify the epsilon -> 0 limit from the hint explicitly
    eps = 1e-4
    jn = zeros[2]
    def hint_expr(e):
        beta = jn + e
        return (beta * J0(jn) * J0_prime(beta) - jn * J0(beta) * J0_prime(jn)) / (jn * jn - beta * beta)
    # One-sided O(eps) convergence and symmetric O(eps^2) convergence
    hint_err_1 = abs(hint_expr(eps) - 0.5 * (J0_prime(jn) ** 2))
    hint_err = abs(0.5 * (hint_expr(eps) + hint_expr(-eps)) - 0.5 * (J0_prime(jn) ** 2))
    assert hint_err_1 < 1e-6, f"One-sided hint limit error: {hint_err_1}"
    assert hint_err < 1e-9, f"Symmetric hint limit error: {hint_err}"

    # Part (iii): Verify eigenfunction expansion against 3 independent exact closed-form solutions
    def solve_eigen_expansion(f_func, lam_tilde, n_terms=20):
        coeffs = []
        for n in range(n_terms):
            jn = zeros[n]
            proj = integrate_gl(lambda xi, a=jn: xi * f_func(xi) * J0(a * xi), 0.0, 1.0, 64)
            cn = 2.0 * proj / ((lam_tilde - jn * jn) * (J0_prime(jn) ** 2))
            coeffs.append((jn, cn))
        return lambda x: sum(cn * J0(jn * x) for jn, cn in coeffs)

    # Case A: lam_tilde = 0, f(x) = 1 => u(x) = (x^2 - 1)/4
    u_exp_A = solve_eigen_expansion(lambda x: 1.0, 0.0, n_terms=22)
    max_err_A = 0.0
    for x in [0.0, 0.2, 0.5, 0.8, 1.0]:
        u_ex = 0.25 * (x * x - 1.0)
        err = abs(u_exp_A(x) - u_ex)
        max_err_A = max(max_err_A, err)
    assert max_err_A < 2e-4, f"Case A error too large: {max_err_A}"

    # Case B: lam_tilde = 10.0 (between j_1^2 and j_2^2), u(x) = (1 - x^2)^3 (so f(1) = 0, fast decay)
    # u'(x) = -6x(1-x^2)^2, (x u')' = (-6x^2(1-x^2)^2)' = x * 12(1-x^2)(3x^2 - 1)
    lam_B = 10.0
    f_B = lambda x: 12.0 * (1.0 - x * x) * (3.0 * x * x - 1.0) + lam_B * ((1.0 - x * x) ** 3)
    u_exp_B_10 = solve_eigen_expansion(f_B, lam_B, n_terms=10)
    u_exp_B = solve_eigen_expansion(f_B, lam_B, n_terms=22)
    max_err_B_10 = max(abs(u_exp_B_10(x) - (1.0 - x * x) ** 3) for x in [0.0, 0.25, 0.5, 0.75, 1.0])
    max_err_B = max(abs(u_exp_B(x) - (1.0 - x * x) ** 3) for x in [0.0, 0.25, 0.5, 0.75, 1.0])
    assert max_err_B < 5e-6 and max_err_B < 0.05 * max_err_B_10, f"Case B error: {max_err_B} (10-term: {max_err_B_10})"

    # Case C: lam_tilde = -3.0 < 0, u(x) = 1 - x^4 => (x u')' = (-4x^4)' = -16x^3 = x(-16x^2)
    lam_C = -3.0
    f_C = lambda x: -16.0 * x * x + lam_C * (1.0 - x ** 4)
    u_exp_C = solve_eigen_expansion(f_C, lam_C, n_terms=22)
    max_err_C = 0.0
    for x in [0.0, 0.3, 0.6, 0.9, 1.0]:
        u_ex = 1.0 - x ** 4
        err = abs(u_exp_C(x) - u_ex)
        max_err_C = max(max_err_C, err)
    assert max_err_C < 1e-3, f"Case C error too large: {max_err_C}"

    # Case D: Single eigenmode exactness f(x) = J_0(j_2 x), lam_tilde = 12.5
    lam_D = 12.5
    j2 = zeros[1]
    u_exp_D = solve_eigen_expansion(lambda x: J0(j2 * x), lam_D, n_terms=10)
    max_err_D = 0.0
    for x in [0.0, 0.33, 0.71, 1.0]:
        u_ex = J0(j2 * x) / (lam_D - j2 * j2)
        err = abs(u_exp_D(x) - u_ex)
        max_err_D = max(max_err_D, err)
    assert max_err_D < 1e-13, f"Case D exact mode error: {max_err_D}"

    print("IB_methods_sheet-2_qx25_check: ALL PASSED")
    print(f"  first 3 zeros: j_1={zeros[0]:.10f}, j_2={zeros[1]:.10f}, j_3={zeros[2]:.10f}")
    print(f"  part (i) Bessel ODE err: {max_ode_err:.3e}")
    print(f"  part (ii) Lommel err: {max_lommel_err:.3e}, orth err: {max_orth_err:.3e}, norm err: {max_norm_err:.3e}, hint-limit err: {hint_err:.3e}")
    print(f"  part (iii) errors: Case A (lam=0, u=(x^2-1)/4)={max_err_A:.3e}, Case B (lam=10, u=(1-x^2)^3)={max_err_B:.3e}, Case C (lam=-3, u=1-x^4)={max_err_C:.3e}, Case D (single mode)={max_err_D:.3e}")


if __name__ == "__main__":
    main()
