"""Computational verification for IB/methods/sheet-1/qx12 (2001-02 Sheet 1 Q6)."""

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


def verify_exact_rationals():
    # int_0^1 u^2 (1 - u)^2 du = int_0^1 (u^2 - 2u^3 + u^4) du
    I = Fraction(1, 3) - Fraction(2, 4) + Fraction(1, 5)
    assert I == Fraction(1, 30)
    # E(0) / (mu * l * V^2) = (1/2) * 16 * I = 8 * I
    E0_coeff = Fraction(8, 1) * I
    assert E0_coeff == Fraction(4, 15)
    # Mode sum prefactor: (1/4) * 32^2 = 256
    mode_prefactor = Fraction(32 * 32, 4)
    assert mode_prefactor == Fraction(256, 1)
    # Ratio gives sum_{n odd} 1/n^6 / pi^6
    series_coeff = E0_coeff / mode_prefactor
    assert series_coeff == Fraction(1, 960)
    print(f"Exact rational check: I={I}, E0_coeff={E0_coeff}, prefactor={mode_prefactor}, ratio={series_coeff}")

    # Numerical check of sum_{n odd} 1/n^6 = pi^6 / 960
    s_odd = sum(1.0 / (n ** 6) for n in range(1, 2001, 2))
    exact_val = (math.pi ** 6) / 960.0
    err = abs(s_odd - exact_val)
    assert err < 1e-14, f"Series sum error = {err}"
    print(f"Numerical sum_{{n odd}} 1/n^6 = {s_odd:.15f}, pi^6/960 = {exact_val:.15f}, err = {err:.2e}")


def verify_string_modes_and_energies():
    pi = math.pi
    param_sets = [
        (1.0, 1.0, 1.0, 1.0),
        (2.3, 7.5, 1.8, 3.1),
        (0.4, 12.0, 2.5, 0.7),
    ]

    for mu, T, ell, V in param_sets:
        c = math.sqrt(T / mu)
        yt0 = lambda x: (4.0 * V / (ell * ell)) * x * (ell - x)

        # Verify modal amplitudes via numerical Fourier sine projection
        for n in range(1, 15):
            wn = n * pi * c / ell
            wn_Bn_num = (2.0 / ell) * simpson(lambda x, n=n: yt0(x) * math.sin(n * pi * x / ell), 0.0, ell, 2000)
            Bn_num = wn_Bn_num / wn
            if n % 2 == 0:
                Bn_exact = 0.0
            else:
                Bn_exact = 32.0 * V * ell / ((n ** 4) * (pi ** 4) * c)
            assert abs(Bn_num - Bn_exact) < 1e-11, f"n={n}: num={Bn_num}, exact={Bn_exact}"

        # Build truncated solution y(x, t) and derivatives
        N_max = 199  # odd modes up to 199
        odd_ns = list(range(1, N_max + 1, 2))
        B = {n: 32.0 * V * ell / ((n ** 4) * (pi ** 4) * c) for n in odd_ns}
        w = {n: n * pi * c / ell for n in odd_ns}

        def y(x, t):
            return sum(B[n] * math.sin(w[n] * t) * math.sin(n * pi * x / ell) for n in odd_ns)

        def y_t(x, t):
            return sum(w[n] * B[n] * math.cos(w[n] * t) * math.sin(n * pi * x / ell) for n in odd_ns)

        def y_x(x, t):
            return sum((n * pi / ell) * B[n] * math.sin(w[n] * t) * math.cos(n * pi * x / ell) for n in odd_ns)

        # Check BCs and ICs
        for t in [0.0, 0.13, 0.71, 2.4]:
            assert abs(y(0.0, t)) < 1e-14
            assert abs(y(ell, t)) < 1e-14
        for frac in [0.05, 0.2, 0.5, 0.83]:
            x = frac * ell
            assert abs(y(x, 0.0)) < 1e-15
            assert abs(y_t(x, 0.0) - yt0(x)) < 1e-5 * V

        # Check wave equation y_tt = c^2 y_xx via central finite differences
        # (use smaller N_max=25 for smooth FD check without high-frequency truncation aliasing)
        smooth_ns = list(range(1, 26, 2))
        def y_smooth(x, t):
            return sum(B[n] * math.sin(w[n] * t) * math.sin(n * pi * x / ell) for n in smooth_ns)

        h_x = 1e-5 * ell
        h_t = 1e-5 * (ell / c)
        max_pde_res = 0.0
        for x in [0.23 * ell, 0.5 * ell, 0.71 * ell]:
            for t in [0.19, 0.83]:
                ytt_fd = (y_smooth(x, t + h_t) - 2.0 * y_smooth(x, t) + y_smooth(x, t - h_t)) / (h_t * h_t)
                yxx_fd = (y_smooth(x + h_x, t) - 2.0 * y_smooth(x, t) + y_smooth(x - h_x, t)) / (h_x * h_x)
                max_pde_res = max(max_pde_res, abs(ytt_fd - c * c * yxx_fd))
        assert max_pde_res < 1e-4, f"max_pde_res = {max_pde_res}"

        # Verify K(t), U(t) by spatial integration vs analytical mode sums, and E(t) = E(0)
        E0_exact = (4.0 / 15.0) * mu * ell * V * V
        E0_num = 0.5 * mu * simpson(lambda x: yt0(x) ** 2, 0.0, ell, 4000)
        assert abs(E0_num - E0_exact) < 1e-12 * E0_exact, f"E0 err = {abs(E0_num - E0_exact)}"

        for t in [0.0, 0.17, 0.63, 1.41, 2.9]:
            K_num = 0.5 * mu * simpson(lambda x, t=t: y_t(x, t) ** 2, 0.0, ell, 600)
            U_num = 0.5 * T * simpson(lambda x, t=t: y_x(x, t) ** 2, 0.0, ell, 600)
            pref = 256.0 * mu * ell * V * V / (pi ** 6)
            K_formula = pref * sum((math.cos(w[n] * t) ** 2) / (n ** 6) for n in odd_ns)
            U_formula = pref * sum((math.sin(w[n] * t) ** 2) / (n ** 6) for n in odd_ns)
            assert abs(K_num - K_formula) < 1e-7 * E0_exact, f"K mismatch at t={t}: {K_num} vs {K_formula}"
            assert abs(U_num - U_formula) < 1e-7 * E0_exact, f"U mismatch at t={t}: {U_num} vs {U_formula}"
            assert abs((K_formula + U_formula) - E0_exact) < 1e-12 * E0_exact
        print(f"Params (mu={mu}, T={T}, l={ell}, V={V}): modes, PDE (res={max_pde_res:.2e}), K(t), U(t), E(t)=E(0)={E0_exact:.6f} verified.")


if __name__ == "__main__":
    verify_exact_rationals()
    verify_string_modes_and_energies()
    print("ALL CHECKS PASSED for IB/methods/sheet-1/qx12.")
