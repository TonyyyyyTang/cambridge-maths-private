"""Computational verification for IB/2009/P3/q15A (Euler-Lagrange & radial electrostatics)."""

import math


def simpson(func, a, b, n=2000):
    """Composite Simpson's rule with even number of subintervals n."""
    assert n % 2 == 0
    h = (b - a) / n
    total = func(a) + func(b)
    for i in range(1, n):
        x = a + i * h
        total += (4.0 if i % 2 == 1 else 2.0) * func(x)
    return total * h / 3.0


def phi_3d(r, R1, R2, lam, phi1):
    """Exact solution of (1/r) d^2(r phi)/dr^2 - phi/lam^2 = 0 with phi(R1)=phi1, phi(R2)=0."""
    return phi1 * (R1 / r) * math.sinh((R2 - r) / lam) / math.sinh((R2 - R1) / lam)


def dphi_3d(r, R1, R2, lam, phi1):
    """Exact analytic derivative d(phi_3d)/dr."""
    c = phi1 * R1 / math.sinh((R2 - R1) / lam)
    arg = (R2 - r) / lam
    return c * (-math.sinh(arg) / (r * r) - math.cosh(arg) / (r * lam))


def energy_3d(phi_fn, dphi_fn, R1, R2, lam, n=2000):
    """Compute E[phi] = 2*pi * int_{R1}^{R2} [0.5*(phi')^2 + 0.5*phi^2/lam^2] r^2 dr."""
    integrand = lambda r: 2.0 * math.pi * (
        0.5 * dphi_fn(r) ** 2 + 0.5 * (phi_fn(r) ** 2) / (lam * lam)
    ) * (r * r)
    return simpson(integrand, R1, R2, n=n)


def bessel_i0(x):
    """Modified Bessel function I_0(x) via power series."""
    term = 1.0
    s = 1.0
    k = 1
    while abs(term) > 1e-16 * abs(s):
        term *= (x * x) / (4.0 * k * k)
        s += term
        k += 1
    return s


def dbessel_i0(x):
    """Derivative I_0'(x) = I_1(x) via power series."""
    term = x / 2.0
    s = term
    k = 1
    while abs(term) > 1e-16 * (abs(s) + 1e-30):
        term *= (x * x) / (4.0 * k * (k + 1))
        s += term
        k += 1
    return s


def test_3d_solution():
    cases = [
        (1.0, 3.0, 1.5, 2.5),
        (0.5, 2.0, 0.7, -1.2),
        (2.0, 5.0, -1.3, 4.0),  # invariance under lambda -> -lambda
    ]
    h = 5e-3
    for R1, R2, lam, phi1 in cases:
        # 1. Boundary conditions
        assert abs(phi_3d(R1, R1, R2, lam, phi1) - phi1) < 1e-14
        assert abs(phi_3d(R2, R1, R2, lam, phi1)) < 1e-14

        # 2. Check both forms of the ODE at interior points (4th-order central differences)
        for frac in (0.15, 0.37, 0.62, 0.88):
            r = R1 + frac * (R2 - R1)
            # Form (1): (1/r) d^2(r phi)/dr^2 - phi/lam^2 = 0
            u = lambda x: x * phi_3d(x, R1, R2, lam, phi1)
            d2u = (
                -u(r + 2 * h) + 16 * u(r + h) - 30 * u(r) + 16 * u(r - h) - u(r - 2 * h)
            ) / (12.0 * h * h)
            res1 = d2u / r - phi_3d(r, R1, R2, lam, phi1) / (lam * lam)
            assert abs(res1) < 1e-10, f"ODE form (1) residual too large: {res1}"

            # Direct Euler-Lagrange form: d/dr(r^2 phi') - (r^2/lam^2) phi = 0
            fl = lambda x: x * x * dphi_3d(x, R1, R2, lam, phi1)
            dfl = (-fl(r + 2 * h) + 8 * fl(r + h) - 8 * fl(r - h) + fl(r - 2 * h)) / (12.0 * h)
            res2 = dfl - (r * r / (lam * lam)) * phi_3d(r, R1, R2, lam, phi1)
            assert abs(res2) < 1e-9, f"EL residual too large: {res2}"

        # 3. Check sign invariance under lambda -> -lambda
        r_mid = 0.5 * (R1 + R2)
        assert abs(
            phi_3d(r_mid, R1, R2, lam, phi1) - phi_3d(r_mid, R1, R2, -lam, phi1)
        ) < 1e-14

    # 4. Variational stationarity & strict minimum check
    R1, R2, lam, phi1 = 1.0, 3.0, 1.2, 2.0
    phi0 = lambda r: phi_3d(r, R1, R2, lam, phi1)
    dphi0 = lambda r: dphi_3d(r, R1, R2, lam, phi1)
    E0 = energy_3d(phi0, dphi0, R1, R2, lam)

    # Analytic minimum energy from integration by parts: E[phi0] = -pi * R1^2 * phi1 * phi0'(R1)
    E0_exact = -math.pi * R1 * R1 * phi1 * dphi0(R1)
    assert abs(E0 - E0_exact) < 1e-11, f"E0={E0}, E0_exact={E0_exact}"

    perturbations = [
        (
            lambda r: math.sin(math.pi * (r - R1) / (R2 - R1)),
            lambda r: (math.pi / (R2 - R1)) * math.cos(math.pi * (r - R1) / (R2 - R1)),
        ),
        (
            lambda r: (r - R1) * (R2 - r) * r,
            lambda r: (R2 - r) * r - (r - R1) * r + (r - R1) * (R2 - r),
        ),
    ]
    for eta, deta in perturbations:
        E_eta = energy_3d(eta, deta, R1, R2, lam)
        assert E_eta > 0.0
        for eps in (-0.3, -1e-4, 1e-4, 0.3):
            E_pert = energy_3d(
                lambda r, e=eps: phi0(r) + e * eta(r),
                lambda r, e=eps: dphi0(r) + e * deta(r),
                R1,
                R2,
                lam,
            )
            # Exact quadratic expansion: E[phi0 + eps*eta] - E[phi0] = eps^2 * E[eta]
            assert E_pert > E0
            assert abs((E_pert - E0) - eps * eps * E_eta) < 1e-11

        # First variation vanishes: check both direct integral and central difference
        direct_first_var = simpson(
            lambda r: 2.0
            * math.pi
            * (dphi0(r) * deta(r) + phi0(r) * eta(r) / (lam * lam))
            * r
            * r,
            R1,
            R2,
            n=4000,
        )
        assert abs(direct_first_var) < 1e-12, f"Direct first variation nonzero: {direct_first_var}"

        eps = 1e-4
        Ep = energy_3d(
            lambda r: phi0(r) + eps * eta(r),
            lambda r: dphi0(r) + eps * deta(r),
            R1,
            R2,
            lam,
        )
        Em = energy_3d(
            lambda r: phi0(r) - eps * eta(r),
            lambda r: dphi0(r) - eps * deta(r),
            R1,
            R2,
            lam,
        )
        first_var = (Ep - Em) / (2.0 * eps)
        assert abs(first_var) < 1e-9, f"FD first variation nonzero: {first_var}"

    # 5. Limiting cases: R2 -> inf (Yukawa) and lam -> inf (Laplace)
    r_test = 1.7
    yukawa = phi1 * (R1 / r_test) * math.exp(-(r_test - R1) / lam)
    assert abs(phi_3d(r_test, R1, 50.0, lam, phi1) - yukawa) < 1e-12
    laplace = phi1 * (R1 / r_test) * (R2 - r_test) / (R2 - R1)
    assert abs(phi_3d(r_test, R1, R2, 1e5, phi1) - laplace) < 1e-9
    print(f"3D checks passed: E[phi0] = {E0:.12f} (exact {E0_exact:.12f})")


def test_2d_annular():
    R1, R2, lam = 1.0, 3.0, 1.4
    # Solution of (1/r) d/dr(r dphi/dr) - phi/lam^2 = 0: phi(r) = I_0(r/lam)
    phi_2d = lambda r: bessel_i0(r / lam)
    dphi_2d = lambda r: dbessel_i0(r / lam) / lam

    # Verify 2D ODE directly
    h = 1e-5
    for r in (1.2, 1.8, 2.5):
        fl_p = (r + h) * dphi_2d(r + h)
        fl_m = (r - h) * dphi_2d(r - h)
        res = (fl_p - fl_m) / (2.0 * h * r) - phi_2d(r) / (lam * lam)
        assert abs(res) < 1e-8, f"2D ODE residual too large: {res}"

    # Verify stationarity of 2D energy functional E_2D[phi] = 2*pi * int [0.5*(phi')^2 + 0.5*phi^2/lam^2] r dr
    def energy_2d(f, df):
        return simpson(
            lambda r: 2.0 * math.pi * (0.5 * df(r) ** 2 + 0.5 * f(r) ** 2 / (lam * lam)) * r,
            R1,
            R2,
            n=2000,
        )

    E0 = energy_2d(phi_2d, dphi_2d)
    eta = lambda r: math.sin(math.pi * (r - R1) / (R2 - R1)) * (r + 0.5)
    deta = lambda r: (
        (math.pi / (R2 - R1)) * math.cos(math.pi * (r - R1) / (R2 - R1)) * (r + 0.5)
        + math.sin(math.pi * (r - R1) / (R2 - R1))
    )
    E_eta = energy_2d(eta, deta)
    for eps in (-0.2, 0.2):
        Ep = energy_2d(
            lambda r, e=eps: phi_2d(r) + e * eta(r),
            lambda r, e=eps: dphi_2d(r) + e * deta(r),
        )
        assert abs((Ep - E0) - eps * eps * E_eta) < 1e-11
    print(f"2D annular checks passed: E_2D[I_0(r/lam)] = {E0:.12f}")


if __name__ == "__main__":
    test_3d_solution()
    test_2d_annular()
    print("ALL CHECKS PASSED for IB_2009_P3_q15A")
