"""Independent numerical verification for IB/2019/P3/q15D.

Verifies:
1. Distributional ODE psi'' + alpha^2 psi = delta(t) for psi(t) = H(t) sin(alpha*t)/alpha
   tested against smooth rapidly decaying test functions phi(t).
2. The given 3D->1D Fourier integral identity (with exponential damping e^{-eps*k} for
   absolute convergence) and the resulting retarded Green's function
       G(x, t; y, tau) = H(t - tau) / (4*pi*c*|x - y|) * delta(|x - y| - c*(t - tau)).
3. Retarded Green's function convolution reproduces an exact inhomogeneous causal solution.
4. Kirchhoff's formula u(x, t) = <f> * t reproduces three independent exact solutions of
   u_tt - c^2 Laplacian(u) = 0 with u(x, 0) = 0, u_t(x, 0) = f(x):
   - plane wave
   - polynomial
   - radial Gaussian pulse (compared against d'Alembert radial solution)
"""

import cmath
import math


def simpson(f, a: float, b: float, n: int) -> float:
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * f(a + i * h)
    return s * h / 3.0


def gauss_legendre_nodes_weights(n: int) -> tuple[list[float], list[float]]:
    """Compute n-point Gauss-Legendre nodes & weights on [-1, 1] via Newton's method."""
    nodes = []
    weights = []
    for k in range(1, n + 1):
        x = math.cos(math.pi * (k - 0.25) / (n + 0.5))
        for _ in range(30):
            p0, p1 = 1.0, x
            for j in range(2, n + 1):
                p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = -p1 / dp
            x += dx
            if abs(dx) < 1e-15:
                break
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * dp * dp))
    return nodes, weights


def sphere_average(f, center: tuple[float, float, float], radius: float, n_mu: int = 24, n_phi: int = 48) -> float:
    """High-order Gauss-Legendre (in mu=cos theta) x trapezoidal (in phi) spherical average."""
    mus, ws = gauss_legendre_nodes_weights(n_mu)
    cx, cy, cz = center
    total = 0.0
    dphi = 2.0 * math.pi / n_phi
    for mu, w in zip(mus, ws):
        sin_th = math.sqrt(max(0.0, 1.0 - mu * mu))
        phi_sum = 0.0
        for j in range(n_phi):
            phi = j * dphi
            x = cx + radius * sin_th * math.cos(phi)
            y = cy + radius * sin_th * math.sin(phi)
            z = cz + radius * mu
            phi_sum += f(x, y, z)
        total += w * (phi_sum / n_phi)
    return 0.5 * total


def main() -> None:
    # 1. Check distributional ODE: int_0^inf (sin(a*t)/a) * (phi''(t) + a^2 phi(t)) dt == phi(0)
    #    Use phi(t) = (1 + 0.7*t - 0.3*t^2) * exp(-0.5*(t - 0.4)^2)
    def phi(t: float) -> float:
        return (1.0 + 0.7 * t - 0.3 * t * t) * math.exp(-0.5 * (t - 0.4) ** 2)

    def phi_dd(t: float) -> float:
        # Exact second derivative of P(t) * exp(g(t)), g(t) = -0.5*(t-0.4)^2
        p = 1.0 + 0.7 * t - 0.3 * t * t
        dp = 0.7 - 0.6 * t
        d2p = -0.6
        dg = -(t - 0.4)
        d2g = -1.0
        return (d2p + 2.0 * dp * dg + p * (d2g + dg * dg)) * math.exp(-0.5 * (t - 0.4) ** 2)

    max_dist_err = 0.0
    for alpha in (0.5, 1.7, 3.3):
        lhs = simpson(
            lambda t, a=alpha: (math.sin(a * t) / a) * (phi_dd(t) + a * a * phi(t)),
            0.0,
            14.0,
            12000,
        )
        max_dist_err = max(max_dist_err, abs(lhs - phi(0.0)))
    print(f"1. Distributional ODE max error:                {max_dist_err:.3e}")
    assert max_dist_err < 1e-11

    # 2. Check regularized 3D->1D Fourier integral identity (with e^{-eps*|k|} damping)
    #    LHS: (1/(2*pi)) * 2*pi * int_0^inf k^2 * (2*sin(k*r)/(k*r)) * (sin(k*c*t)/(k*c)) * e^{-eps*k} dk
    #    RHS: (-i/(c*r)) * int_{-inf}^{inf} e^{i*k*r} * sin(k*c*t) * e^{-eps*|k|} dk
    r, c, t, eps = 1.35, 1.8, 0.65, 0.25
    lhs_reg = simpson(
        lambda k: (2.0 * math.sin(k * r) / (c * r)) * math.sin(k * c * t) * math.exp(-eps * k),
        0.0,
        80.0,
        8000,
    )
    rhs_real = simpson(
        lambda k: ((-1j / (c * r)) * cmath.exp(1j * k * r) * math.sin(k * c * t) * math.exp(-eps * abs(k))).real,
        -80.0,
        80.0,
        16000,
    )
    rhs_imag = simpson(
        lambda k: ((-1j / (c * r)) * cmath.exp(1j * k * r) * math.sin(k * c * t) * math.exp(-eps * abs(k))).imag,
        -80.0,
        80.0,
        16000,
    )
    id_err = abs(lhs_reg - rhs_real) + abs(rhs_imag)
    print(f"2. Regularized 3D->1D hint identity error:      {id_err:.3e}")
    assert id_err < 1e-10

    # 3. Check retarded Green's function on inhomogeneous causal PDE:
    #    u(x, t) = t^3 * exp(-|x|^2) for t >= 0 (0 for t < 0)
    #    F = u_tt - c^2 * Laplacian(u) = (6*t - c^2 * t^3 * (4*|x|^2 - 6)) * exp(-|x|^2)
    #    Retarded solution: u(x, t) = int_0^t dtau (t - tau) * <F(., tau)>_{S(x, c*(t-tau))}
    c_val = 1.4
    x0 = (0.3, -0.4, 0.25)
    t_val = 0.85

    def source_F(x: float, y: float, z: float, tau: float) -> float:
        r2 = x * x + y * y + z * z
        return (6.0 * tau - c_val * c_val * (tau ** 3) * (4.0 * r2 - 6.0)) * math.exp(-r2)

    u_from_G = simpson(
        lambda tau: (t_val - tau) * sphere_average(
            lambda x, y, z: source_F(x, y, z, tau), x0, c_val * (t_val - tau), n_mu=16, n_phi=32
        ),
        0.0,
        t_val,
        40,
    )
    u_exact_inhom = (t_val ** 3) * math.exp(-(x0[0] ** 2 + x0[1] ** 2 + x0[2] ** 2))
    inhom_err = abs(u_from_G - u_exact_inhom)
    print(f"3. Inhomogeneous retarded Green's fn error:     {inhom_err:.3e}")
    assert inhom_err < 1e-9

    # 4. Check Kirchhoff's formula u(x, t) = <f> * t on 3 exact homogeneous IVP solutions
    # 4a. Plane wave f(x) = cos(k . x + phi0)
    kx, ky, kz, phi0 = 1.2, -0.7, 0.9, 0.4
    kmag = math.sqrt(kx * kx + ky * ky + kz * kz)
    t1 = 0.72
    u_pw_kirch = t1 * sphere_average(
        lambda x, y, z: math.cos(kx * x + ky * y + kz * z + phi0), x0, c_val * t1
    )
    u_pw_exact = (math.sin(kmag * c_val * t1) / (kmag * c_val)) * math.cos(
        kx * x0[0] + ky * x0[1] + kz * x0[2] + phi0
    )
    err_pw = abs(u_pw_kirch - u_pw_exact)
    print(f"4a. Kirchhoff plane-wave error:                 {err_pw:.3e}")
    assert err_pw < 1e-12

    # 4b. Polynomial f(x) = |x|^2 + x*y - 2*z
    u_poly_kirch = t1 * sphere_average(
        lambda x, y, z: x * x + y * y + z * z + x * y - 2.0 * z, x0, c_val * t1
    )
    u_poly_exact = t1 * (x0[0] ** 2 + x0[1] ** 2 + x0[2] ** 2 + (c_val * t1) ** 2 + x0[0] * x0[1] - 2.0 * x0[2])
    err_poly = abs(u_poly_kirch - u_poly_exact)
    print(f"4b. Kirchhoff polynomial error:                 {err_poly:.3e}")
    assert err_poly < 1e-12

    # 4c. Radial Gaussian f(x) = exp(-a*|x|^2), exact radial d'Alembert solution:
    #     u(r, t) = (exp(-a*(r - c*t)^2) - exp(-a*(r + c*t)^2)) / (4*a*c*r)
    a_val = 1.5
    r0 = math.sqrt(x0[0] ** 2 + x0[1] ** 2 + x0[2] ** 2)
    u_gauss_kirch = t1 * sphere_average(
        lambda x, y, z: math.exp(-a_val * (x * x + y * y + z * z)), x0, c_val * t1
    )
    u_gauss_exact = (
        math.exp(-a_val * (r0 - c_val * t1) ** 2) - math.exp(-a_val * (r0 + c_val * t1) ** 2)
    ) / (4.0 * a_val * c_val * r0)
    err_gauss = abs(u_gauss_kirch - u_gauss_exact)
    print(f"4c. Kirchhoff radial Gaussian error:            {err_gauss:.3e}")
    assert err_gauss < 1e-12

    # 4d. Check strong Huygens' principle (sharp trailing edge -> exact 0 after wave passes)
    #     Compactly supported f(y) = (1 - |y|^2)^2 for |y| < 1, 0 otherwise; observe at |x| = 2.5
    #     Wave must be strictly 0 for c*t < 1.5 and for c*t > 3.5!
    x_obs = (2.5, 0.0, 0.0)

    def f_compact(x: float, y: float, z: float) -> float:
        r2 = x * x + y * y + z * z
        return (1.0 - r2) ** 2 if r2 < 1.0 else 0.0

    u_before = (1.0 / c_val) * sphere_average(f_compact, x_obs, 1.0)
    u_during = (2.5 / c_val) * sphere_average(f_compact, x_obs, 2.5)
    u_after = (4.0 / c_val) * sphere_average(f_compact, x_obs, 4.0)
    print(f"4d. Huygens compact source: before={u_before:.1e}, during={u_during:.4f}, after={u_after:.1e}")
    assert u_before == 0.0 and u_during > 0.0 and u_after == 0.0
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
