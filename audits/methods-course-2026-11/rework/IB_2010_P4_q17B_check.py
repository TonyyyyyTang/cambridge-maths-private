"""Computational verification for IB/2010/P4/q17B (Green's 3rd identity & half-space Neumann problem)."""

import math


def gauss_legendre_nodes_weights(n):
    """Compute n-point Gauss-Legendre nodes and weights on [-1, 1] using Newton's method."""
    nodes = [0.0] * n
    weights = [0.0] * n
    m = (n + 1) // 2
    for i in range(m):
        # Initial guess
        x = math.cos(math.pi * (i + 0.75) / (n + 0.5))
        for _ in range(50):
            p0 = 1.0
            p1 = x
            for k in range(2, n + 1):
                pk = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                p0, p1 = p1, pk
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = -p1 / dp
            x += dx
            if abs(dx) < 1e-15:
                break
        w = 2.0 / ((1.0 - x * x) * dp * dp)
        nodes[i] = -x
        nodes[n - 1 - i] = x
        weights[i] = w
        weights[n - 1 - i] = w
    return nodes, weights


def test_greens_third_identity():
    """Verify u(r0) = int_S (u dG_{f3}/dn - (du/dn) G_{f3}) dS on a sphere S of radius R."""
    R = 2.5
    r0 = (0.6, -0.4, 0.9)  # interior point (|r0| < R)

    # Test harmonic functions:
    # 1. Harmonic polynomial: u(x,y,z) = x^2 - y^2 + 3*x*z + 2*y - 5
    def u1(x, y, z):
        return x * x - y * y + 3.0 * x * z + 2.0 * y - 5.0

    def grad_u1(x, y, z):
        return (2.0 * x + 3.0 * z, -2.0 * y + 2.0, 3.0 * x)

    # 2. Exterior point source: u(x,y,z) = 1 / |r - r_ext| with |r_ext| > R
    rext = (4.0, -2.0, 3.5)

    def u2(x, y, z):
        dx, dy, dz = x - rext[0], y - rext[1], z - rext[2]
        return 1.0 / math.sqrt(dx * dx + dy * dy + dz * dz)

    def grad_u2(x, y, z):
        dx, dy, dz = x - rext[0], y - rext[1], z - rext[2]
        d3 = (dx * dx + dy * dy + dz * dz) ** 1.5
        return (-dx / d3, -dy / d3, -dz / d3)

    mu_nodes, mu_weights = gauss_legendre_nodes_weights(40)  # mu = cos(theta) in [-1, 1]
    n_phi = 80
    dphi = 2.0 * math.pi / n_phi

    for name, u, gu in (("harmonic_poly", u1, grad_u1), ("ext_source", u2, grad_u2)):
        integral = 0.0
        for mu, wmu in zip(mu_nodes, mu_weights):
            sin_t = math.sqrt(max(0.0, 1.0 - mu * mu))
            for j in range(n_phi):
                phi = (j + 0.5) * dphi
                nx = sin_t * math.cos(phi)
                ny = sin_t * math.sin(phi)
                nz = mu
                x, y, z = R * nx, R * ny, R * nz
                rx, ry, rz = x - r0[0], y - r0[1], z - r0[2]
                dist2 = rx * rx + ry * ry + rz * rz
                dist = math.sqrt(dist2)
                G = -1.0 / (4.0 * math.pi * dist)
                # grad G = (r - r0) / (4*pi*dist^3)
                dG_dn = (rx * nx + ry * ny + rz * nz) / (4.0 * math.pi * dist2 * dist)
                gx, gy, gz = gu(x, y, z)
                du_dn = gx * nx + gy * ny + gz * nz
                integral += (u(x, y, z) * dG_dn - du_dn * G) * (R * R) * wmu * dphi

        exact = u(*r0)
        err = abs(integral - exact)
        assert err < 1e-11, f"Green's 3rd identity failed for {name}: err={err}"
        print(f"Green's 3rd identity ({name}): integral={integral:.12f}, exact={exact:.12f}, err={err:.2e}")


def test_part2_flux_condition_and_counterexample():
    """Verify Divergence Theorem flux balance and the a=1 counterexample vs a=2 dipole."""
    # Counterexample for a = 1 (monopole at (0,0,-1)): u = 1 / sqrt(x^2 + y^2 + (z+1)^2)
    # Analytically:
    #   int_{D_R} p dx dy = -2*pi * (1 - 1/sqrt(R^2 + 1)) -> -2*pi != 0
    #   int_{H_R} (du/dr) dS = -2*pi * (1 - 1/sqrt(R^2 + 1)) -> -2*pi
    # Check both integrals numerically via Gauss-Legendre in mu = cos(theta) in [0, 1]
    mu_nodes, mu_weights = gauss_legendre_nodes_weights(40)
    for R in (10.0, 100.0, 1000.0):
        # Disk integral: int_0^R (-1/(rho^2+1)^(3/2)) * 2*pi*rho drho
        disk_exact = -2.0 * math.pi * (1.0 - 1.0 / math.sqrt(R * R + 1.0))
        # Hemisphere integral: int_0^1 (du/dr at r=R) * 2*pi*R^2 dmu
        # u(r, mu) = (r^2 + 2*r*mu + 1)^(-1/2) => du/dr = -(r + mu) / (r^2 + 2*r*mu + 1)^(3/2)
        hemi = 0.0
        for t, wt in zip(mu_nodes, mu_weights):
            mu = 0.5 * (t + 1.0)
            wmu = 0.5 * wt
            du_dr = -(R + mu) / ((R * R + 2.0 * R * mu + 1.0) ** 1.5)
            hemi += du_dr * (2.0 * math.pi * R * R) * wmu
        assert abs(hemi - disk_exact) < 1e-11
        assert abs(disk_exact - (-2.0 * math.pi)) < 2.0 * math.pi / R

    # Contrast with a = 2 > 1 (vertical dipole u = (z+1) / (x^2 + y^2 + (z+1)^2)^(3/2)):
    # Wait: even better, horizontal dipole u = x / (x^2 + y^2 + (z+1)^2)^(3/2) (antisymmetric => flux = 0)
    # AND vertical dipole u = (z+1) / (x^2 + y^2 + (z+1)^2)^(3/2), which is radially symmetric in (x,y)!
    # For vertical dipole: p(x,y) = du/dz|_{z=0} = (rho^2 - 2) / (rho^2 + 1)^(5/2)
    # int_0^R p * 2*pi*rho drho = -2*pi * R^2 / (R^2 + 1)^(3/2) = O(1/R) -> 0 as R -> inf!
    for R in (10.0, 100.0, 1000.0):
        disk_dip = -2.0 * math.pi * R * R / ((R * R + 1.0) ** 1.5)
        hemi_dip = 0.0
        for t, wt in zip(mu_nodes, mu_weights):
            mu = 0.5 * (t + 1.0)
            wmu = 0.5 * wt
            # u(r, mu) = (r*mu + 1) / (r^2 + 2*r*mu + 1)^(3/2)
            D = R * R + 2.0 * R * mu + 1.0
            du_dr = mu / (D ** 1.5) - 3.0 * (R * mu + 1.0) * (R + mu) / (D ** 2.5)
            hemi_dip += du_dr * (2.0 * math.pi * R * R) * wmu
        assert abs(hemi_dip - disk_dip) < 1e-11
        assert abs(disk_dip) < 2.0 * math.pi / R

    print("Part 2 flux checks passed: a=1 monopole has flux -2*pi != 0; a=2 dipole has flux O(1/R) -> 0")


def test_part3_green_function_and_formula():
    """Verify Neumann Green's function G and integral formula u = int p * f dx dy."""
    x0, y0, z0 = 0.7, -0.5, 1.3

    def G_neumann(x, y, z):
        r1 = math.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)
        r2 = math.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z + z0) ** 2)
        return -1.0 / (4.0 * math.pi * r1) - 1.0 / (4.0 * math.pi * r2)

    # Check dG/dz = 0 on z = 0
    h = 1e-6
    for x, y in ((0.0, 0.0), (1.2, -0.8), (-2.0, 3.0)):
        dG_dz = (G_neumann(x, y, h) - G_neumann(x, y, -h)) / (2.0 * h)
        assert abs(dG_dz) < 1e-12

    # Check harmonicity of G away from (x0, y0, z0) in z > 0
    xt, yt, zt = -0.4, 0.9, 2.1
    h2 = 1e-5
    lap = (
        G_neumann(xt + h2, yt, zt)
        + G_neumann(xt - h2, yt, zt)
        + G_neumann(xt, yt + h2, zt)
        + G_neumann(xt, yt - h2, zt)
        + G_neumann(xt, yt, zt + h2)
        + G_neumann(xt, yt, zt - h2)
        - 6.0 * G_neumann(xt, yt, zt)
    ) / (h2 * h2)
    assert abs(lap) < 1e-7

    # Test integral representation on the vertical dipole u(x,y,z) = (z+1)/(x^2+y^2+(z+1)^2)^(3/2)
    # Evaluate at (0, 0, z0) using polar coordinates with substitution rho = tan(theta), theta in [0, pi/2)
    # p(rho) = (rho^2 - 2) / (rho^2 + 1)^(5/2), f(rho, z0) = -1 / (2*pi*sqrt(rho^2 + z0^2))
    nodes, weights = gauss_legendre_nodes_weights(80)
    for z_test in (0.5, 1.3, 3.0):
        val = 0.0
        for s, ws in zip(nodes, weights):
            theta = 0.25 * math.pi * (s + 1.0)
            wtheta = 0.25 * math.pi * ws
            rho = math.tan(theta)
            drho = 1.0 / (math.cos(theta) ** 2)
            p_val = (rho * rho - 2.0) / ((rho * rho + 1.0) ** 2.5)
            f_val = -1.0 / (2.0 * math.pi * math.sqrt(rho * rho + z_test * z_test))
            val += 2.0 * math.pi * rho * p_val * f_val * drho * wtheta
        exact = 1.0 / ((z_test + 1.0) ** 2)
        assert abs(val - exact) < 1e-11, f"Integral formula mismatch at z0={z_test}: {val} vs {exact}"

    # Also test on monopole u(x,y,z) = 1/sqrt(x^2+y^2+(z+1)^2) at (0,0,z0):
    for z_test in (0.5, 1.3, 3.0):
        val = 0.0
        for s, ws in zip(nodes, weights):
            theta = 0.25 * math.pi * (s + 1.0)
            wtheta = 0.25 * math.pi * ws
            rho = math.tan(theta)
            drho = 1.0 / (math.cos(theta) ** 2)
            p_val = -1.0 / ((rho * rho + 1.0) ** 1.5)
            f_val = -1.0 / (2.0 * math.pi * math.sqrt(rho * rho + z_test * z_test))
            val += 2.0 * math.pi * rho * p_val * f_val * drho * wtheta
        exact = 1.0 / (z_test + 1.0)
        assert abs(val - exact) < 1e-11, f"Monopole integral mismatch at z0={z_test}: {val} vs {exact}"

    print("Part 3 Neumann Green's function & integral formula checks passed!")


def test_part4_asymptotics():
    """Verify u(x0,y0,z0) -> -2*a^4*x0 / (3*pi*z0^3) as z0 -> inf."""
    nodes, weights = gauss_legendre_nodes_weights(24)

    def u_exact_num(x0, y0, z0, a):
        total = 0.0
        for sx, wx in zip(nodes, weights):
            x = a * sx
            for sy, wy in zip(nodes, weights):
                y = a * sy
                f_val = -1.0 / (
                    2.0 * math.pi * math.sqrt((x - x0) ** 2 + (y - y0) ** 2 + z0 * z0)
                )
                total += (a * wx) * (a * wy) * x * f_val
        return total

    for a, x0, y0 in ((1.0, 1.5, -0.7), (2.0, -0.8, 1.2), (1.3, 0.0, 0.5)):
        asymp_coeff = -2.0 * (a ** 4) * x0 / (3.0 * math.pi)
        prev_scaled_err = None
        for z0 in (20.0, 50.0, 100.0, 200.0, 500.0):
            u_val = u_exact_num(x0, y0, z0, a)
            if x0 == 0.0:
                assert abs(u_val) < 1e-15
                continue
            asymp = asymp_coeff / (z0 ** 3)
            rel_err = abs(u_val - asymp) / abs(asymp)
            # Remainder should scale like O(z0^-2) relative to leading O(z0^-3)
            scaled_err = rel_err * (z0 ** 2)
            if prev_scaled_err is not None:
                # scaled_err should stabilize to a constant as z0 grows
                assert abs(scaled_err - prev_scaled_err) / max(scaled_err, 1e-12) < 0.05
            prev_scaled_err = scaled_err
        if x0 != 0.0:
            u500 = u_exact_num(x0, y0, 500.0, a)
            ratio = u500 / (asymp_coeff / (500.0 ** 3))
            assert abs(ratio - 1.0) < 1e-4
            print(
                f"Part 4 (a={a}, x0={x0}, y0={y0}): z0=500 ratio u/u_asymp = {ratio:.10f}, "
                f"z0^2 * rel_err -> {prev_scaled_err:.6f}"
            )


if __name__ == "__main__":
    test_greens_third_identity()
    test_part2_flux_condition_and_counterexample()
    test_part3_green_function_and_formula()
    test_part4_asymptotics()
    print("ALL CHECKS PASSED for IB_2010_P4_q17B")
