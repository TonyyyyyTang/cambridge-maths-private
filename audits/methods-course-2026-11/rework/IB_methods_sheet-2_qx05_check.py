"""Computational verification for IB/methods/sheet-2/qx05."""
import math


def gauss_legendre_8():
    # 8-point Gauss-Legendre nodes and weights on [-1, 1]
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


def integrate_gl(func, a, b, n_panels=64):
    nodes, weights = gauss_legendre_8()
    h = (b - a) / n_panels
    total = 0.0
    for i in range(n_panels):
        x0 = a + i * h
        x1 = x0 + h
        mid = 0.5 * (x0 + x1)
        half = 0.5 * h
        s = 0.0
        for xi, wi in zip(nodes, weights):
            s += wi * func(mid + half * xi)
        total += half * s
    return total


def y_series(x, t, A, L, c, lam, K=1200):
    """Fourier sine series solution for t > 0."""
    val = A * math.sin(math.pi * x / L) * math.sin(math.pi * c * t / L)
    s = 0.0
    for k in range(K):
        m = 2 * k + 1
        sign = 1.0 if (k % 2 == 0) else -1.0
        s += (sign / m) * math.sin(m * math.pi * x / L) * math.sin(m * math.pi * c * t / L)
    return val + (2.0 * lam / (math.pi * c)) * s


def square_wave_S(u, L):
    """2L-periodic square wave: sgn(cos(pi*u/L))."""
    val = math.cos(math.pi * u / L)
    assert abs(val) > 1e-12, f"Point u={u} is on a discontinuity"
    return 1.0 if val > 0.0 else -1.0


def y_dalembert(x, t, A, L, c, lam):
    """Exact closed-form (fundamental mode + d'Alembert piecewise constant impulse response)."""
    y0 = A * math.sin(math.pi * x / L) * math.sin(math.pi * c * t / L)
    y1 = (lam / (4.0 * c)) * (square_wave_S(x - c * t, L) - square_wave_S(x + c * t, L))
    return y0 + y1


def main():
    A = 1.3
    L = 2.7
    c = 1.9
    lam = 2.4

    # 1. Verify wave PDE and BCs on truncated Fourier series (exact for each mode)
    K_pde = 30
    h = 2e-4
    max_pde_err = 0.0
    for x_frac in [0.17, 0.41, 0.63, 0.88]:
        x = x_frac * L
        for t in [0.23, 0.71, 1.54]:
            # 4th-order central finite difference for second derivatives
            ytt = (-y_series(x, t + 2 * h, A, L, c, lam, K_pde)
                   + 16.0 * y_series(x, t + h, A, L, c, lam, K_pde)
                   - 30.0 * y_series(x, t, A, L, c, lam, K_pde)
                   + 16.0 * y_series(x, t - h, A, L, c, lam, K_pde)
                   - y_series(x, t - 2 * h, A, L, c, lam, K_pde)) / (12.0 * h * h)
            yxx = (-y_series(x + 2 * h, t, A, L, c, lam, K_pde)
                   + 16.0 * y_series(x + h, t, A, L, c, lam, K_pde)
                   - 30.0 * y_series(x, t, A, L, c, lam, K_pde)
                   + 16.0 * y_series(x - h, t, A, L, c, lam, K_pde)
                   - y_series(x - 2 * h, t, A, L, c, lam, K_pde)) / (12.0 * h * h)
            err = abs(ytt - c * c * yxx)
            max_pde_err = max(max_pde_err, err)
            assert err < 1e-5, f"PDE mismatch: {err}"

    # BCs at x=0 and x=L
    for t in [0.1, 0.7, 2.3]:
        assert abs(y_series(0.0, t, A, L, c, lam, K_pde)) < 1e-15
        assert abs(y_series(L, t, A, L, c, lam, K_pde)) < 1e-14
        assert abs(y_dalembert(0.0, t, A, L, c, lam)) < 1e-15
        assert abs(y_dalembert(L, t, A, L, c, lam)) < 1e-15

    # Initial displacement y(x, 0+) = 0
    for x_frac in [0.1, 0.35, 0.5, 0.77]:
        assert abs(y_series(x_frac * L, 0.0, A, L, c, lam, K_pde)) < 1e-15

    # 2. Compare Fourier series against independent d'Alembert closed form across multiple reflections
    # Use Cesaro / large-K sum away from discontinuities
    max_dal_err = 0.0
    test_points = [
        (0.15, 0.12), (0.48, 0.12), (0.82, 0.12),  # ct/L = 0.12 (before first reflection)
        (0.25, 0.38), (0.50, 0.38), (0.75, 0.38),  # ct/L = 0.38
        (0.18, 0.72), (0.52, 0.72), (0.87, 0.72),  # ct/L = 0.72 (after reflection at ends!)
        (0.22, 1.33), (0.61, 1.33),                # ct/L = 1.33 (negative plateau phase)
        (0.34, 1.81), (0.50, 1.81),                # ct/L = 1.81
    ]
    for x_frac, ct_frac in test_points:
        x = x_frac * L
        t = ct_frac * L / c
        val_ser = y_series(x, t, A, L, c, lam, K=6000)
        val_dal = y_dalembert(x, t, A, L, c, lam)
        err = abs(val_ser - val_dal)
        max_dal_err = max(max_dal_err, err)
        assert err < 5e-4, f"Series vs d'Alembert mismatch at ({x_frac}, {ct_frac}): {val_ser} vs {val_dal}"

    # Also verify Dirichlet midpoint convergence right ON a characteristic jump: x/L=0.31, ct/L=1.81 (x-ct = -1.5 L)
    x_j, t_j = 0.31 * L, 1.81 * L / c
    mid_dal = (A * math.sin(math.pi * x_j / L) * math.sin(math.pi * c * t_j / L)
               + (lam / (4.0 * c)) * (0.0 - square_wave_S(x_j + c * t_j, L)))
    assert abs(y_series(x_j, t_j, A, L, c, lam, K=6000) - mid_dal) < 5e-4

    # 3. Verify distributional initial velocity using smooth test function phi(x)
    # phi(x) = sin(pi*x/L) + 0.7*sin(2*pi*x/L) - 0.4*sin(3*pi*x/L)
    def phi(x):
        return (math.sin(math.pi * x / L)
                + 0.7 * math.sin(2.0 * math.pi * x / L)
                - 0.4 * math.sin(3.0 * math.pi * x / L))

    expected_pairing = (integrate_gl(lambda x: phi(x) * (math.pi * c * A / L) * math.sin(math.pi * x / L), 0.0, L)
                        + lam * phi(0.5 * L))
    dt = 1e-5
    # Exact pairing using d'Alembert expression at small t = dt:
    # For 0 < c*dt < L/2, y1(x, dt) = lam/(2c) on (L/2 - c*dt, L/2 + c*dt) and 0 elsewhere
    y0_part = integrate_gl(lambda x: phi(x) * A * math.sin(math.pi * x / L) * math.sin(math.pi * c * dt / L) / dt, 0.0, L)
    y1_part = integrate_gl(lambda x: phi(x) * (lam / (2.0 * c * dt)), 0.5 * L - c * dt, 0.5 * L + c * dt, n_panels=16)
    actual_pairing = y0_part + y1_part
    pairing_err = abs(actual_pairing - expected_pairing)
    assert pairing_err < 1e-8, f"Distributional velocity mismatch: {pairing_err}"

    # 4. Verify Gaussian-smoothed impulse approaches exact Fourier coefficients b_n
    sigma = 1e-4 * L
    max_coeff_err = 0.0
    for n in range(1, 10):
        # Exact coefficient b_n = A*delta_{n1} + 2*lam/(n*pi*c)*sin(n*pi/2)
        b_exact = (A if n == 1 else 0.0) + (2.0 * lam / (n * math.pi * c)) * math.sin(0.5 * n * math.pi)
        # Numerical projection of smoothed initial velocity v_sigma(x)
        def v_sigma(x):
            base = (math.pi * c * A / L) * math.sin(math.pi * x / L)
            gauss = lam * math.exp(-0.5 * ((x - 0.5 * L) / sigma) ** 2) / (math.sqrt(2.0 * math.pi) * sigma)
            return base + gauss

        # Integrate tightly around midpoint plus rest of interval
        a_mid = 0.5 * L - 12.0 * sigma
        b_mid = 0.5 * L + 12.0 * sigma
        proj = (integrate_gl(lambda x: v_sigma(x) * math.sin(n * math.pi * x / L), 0.0, a_mid, 32)
                + integrate_gl(lambda x: v_sigma(x) * math.sin(n * math.pi * x / L), a_mid, b_mid, 32)
                + integrate_gl(lambda x: v_sigma(x) * math.sin(n * math.pi * x / L), b_mid, L, 32))
        b_num = (2.0 / (n * math.pi * c)) * proj
        err = abs(b_num - b_exact)
        max_coeff_err = max(max_coeff_err, err)
        assert err < 1e-6, f"Coefficient mismatch for n={n}: {b_num} vs {b_exact}"

    # 5. Fundamental-mode annihilation check when lam = -pi*c*A/2
    lam_cancel = -0.5 * math.pi * c * A
    b1_cancel = A + 2.0 * lam_cancel / (math.pi * c)
    assert abs(b1_cancel) < 1e-15

    print("IB_methods_sheet-2_qx05_check: ALL PASSED")
    print(f"  max PDE finite-diff error (K={K_pde}): {max_pde_err:.3e}")
    print(f"  max series vs d'Alembert error (K=6000): {max_dal_err:.3e}")
    print(f"  distributional initial velocity error: {pairing_err:.3e}")
    print(f"  max smoothed-delta coefficient error (n=1..9): {max_coeff_err:.3e}")


if __name__ == "__main__":
    main()
