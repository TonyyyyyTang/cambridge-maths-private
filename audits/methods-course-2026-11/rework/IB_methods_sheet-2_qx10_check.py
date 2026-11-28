"""Computational verification for IB/methods/sheet-2/qx10."""
import math


def gauss_legendre_6():
    nodes = [
        -0.93246951420315203,
        -0.66120938646626451,
        -0.23861918608319691,
         0.23861918608319691,
         0.66120938646626451,
         0.93246951420315203,
    ]
    weights = [
        0.17132449237917035,
        0.36076157304813861,
        0.46791393457269105,
        0.46791393457269105,
        0.36076157304813861,
        0.17132449237917035,
    ]
    return nodes, weights


def quad_1d_nodes_weights(a, b, n_panels=6):
    nodes, weights = gauss_legendre_6()
    h = (b - a) / n_panels
    xs, ws = [], []
    for i in range(n_panels):
        mid = a + (i + 0.5) * h
        half = 0.5 * h
        for xi, wi in zip(nodes, weights):
            xs.append(mid + half * xi)
            ws.append(half * wi)
    return xs, ws


def integrate_3d_separable(fx, fy, fz, a, b, d, n_panels=6):
    xs, wxs = quad_1d_nodes_weights(0.0, a, n_panels)
    ys, wys = quad_1d_nodes_weights(0.0, b, n_panels)
    zs, wzs = quad_1d_nodes_weights(0.0, d, n_panels)
    ix = sum(w * fx(x) for x, w in zip(xs, wxs))
    iy = sum(w * fy(y) for y, w in zip(ys, wys))
    iz = sum(w * fz(z) for z, w in zip(zs, wzs))
    return ix * iy * iz


def integrate_3d_general(func, a, b, d, n_panels=3):
    xs, wxs = quad_1d_nodes_weights(0.0, a, n_panels)
    ys, wys = quad_1d_nodes_weights(0.0, b, n_panels)
    zs, wzs = quad_1d_nodes_weights(0.0, d, n_panels)
    total = 0.0
    for x, wx in zip(xs, wxs):
        for y, wy in zip(ys, wys):
            wxy = wx * wy
            for z, wz in zip(zs, wzs):
                total += wxy * wz * func(x, y, z)
    return total


def main():
    # Domain Omega = (0, a) x (0, b) x (0, d)
    a, b, d = 1.3, 1.7, 2.1

    def psi(mode, x, y, z):
        p, q, r = mode
        return (math.sin(p * math.pi * x / a)
                * math.sin(q * math.pi * y / b)
                * math.sin(r * math.pi * z / d))

    def grad_sq_psi(mode, x, y, z):
        p, q, r = mode
        sx, cx = math.sin(p * math.pi * x / a), math.cos(p * math.pi * x / a)
        sy, cy = math.sin(q * math.pi * y / b), math.cos(q * math.pi * y / b)
        sz, cz = math.sin(r * math.pi * z / d), math.cos(r * math.pi * z / d)
        gx = (p * math.pi / a) * cx * sy * sz
        gy = (q * math.pi / b) * sx * cy * sz
        gz = (r * math.pi / d) * sx * sy * cz
        return gx * gx + gy * gy + gz * gz

    def eigval(mode):
        p, q, r = mode
        return math.pi ** 2 * ((p / a) ** 2 + (q / b) ** 2 + (r / d) ** 2)

    # Part (i): Verify lambda_n = int |grad psi_n|^2 / int psi_n^2 > 0 and orthogonality for distinct eigenvalues
    modes = [(1, 1, 1), (2, 1, 1), (1, 2, 1), (1, 1, 2), (2, 2, 1)]
    max_rayleigh_err = 0.0
    max_orth_err = 0.0
    for m in modes:
        lam = eigval(m)
        assert lam > 0.0
        num = integrate_3d_general(lambda x, y, z, mm=m: grad_sq_psi(mm, x, y, z), a, b, d, n_panels=4)
        den = integrate_3d_general(lambda x, y, z, mm=m: psi(mm, x, y, z) ** 2, a, b, d, n_panels=4)
        err = abs(num / den - lam)
        max_rayleigh_err = max(max_rayleigh_err, err)
        assert err < 1e-11, f"Rayleigh quotient mismatch: {err}"

    for i in range(len(modes)):
        for j in range(i + 1, len(modes)):
            assert abs(eigval(modes[i]) - eigval(modes[j])) > 1e-6
            ip = integrate_3d_general(
                lambda x, y, z, mi=modes[i], mj=modes[j]: psi(mi, x, y, z) * psi(mj, x, y, z),
                a, b, d, n_panels=4,
            )
            max_orth_err = max(max_orth_err, abs(ip))
            assert abs(ip) < 1e-12, f"Orthogonality failed: {ip}"

    # Part (ii): Degenerate eigenspace on cube (0, pi)^3 for lambda = 1^2 + 1^2 + 2^2 = 6
    cube_modes = [(1, 1, 2), (1, 2, 1), (2, 1, 1)]
    def cpsi(mode, x, y, z):
        p, q, r = mode
        return math.sin(p * x) * math.sin(q * y) * math.sin(r * z)

    # Non-orthogonal linear combinations in E_6
    coeffs_raw = [
        [1.0, 2.0, -1.0],
        [3.0, -1.0, 2.0],
        [1.0, 1.0, 1.0],
    ]
    # Inner product of cpsi(m1) and cpsi(m2) on (0,pi)^3 is (pi/2)^3 * delta_{m1, m2}
    # Perform Gram-Schmidt via numerical 3D integration
    def make_lincomb(c):
        return lambda x, y, z: sum(ci * cpsi(m, x, y, z) for ci, m in zip(c, cube_modes))

    def ip_cube(f1, f2):
        return integrate_3d_general(lambda x, y, z: f1(x, y, z) * f2(x, y, z), math.pi, math.pi, math.pi, n_panels=4)

    # Confirm raw basis is NOT orthogonal
    raw_funcs = [make_lincomb(c) for c in coeffs_raw]
    assert abs(ip_cube(raw_funcs[0], raw_funcs[1])) > 1.0

    # Gram-Schmidt on coefficient vectors using numerical inner products
    gs_coeffs = []
    for c in coeffs_raw:
        v = list(c)
        for u in gs_coeffs:
            fu = make_lincomb(u)
            fv = make_lincomb(v)
            proj = ip_cube(fv, fu) / ip_cube(fu, fu)
            v = [vi - proj * ui for vi, ui in zip(v, u)]
        gs_coeffs.append(v)

    gs_funcs = [make_lincomb(c) for c in gs_coeffs]
    max_gs_orth_err = 0.0
    for i in range(3):
        assert ip_cube(gs_funcs[i], gs_funcs[i]) > 1e-6
        for j in range(i + 1, 3):
            err = abs(ip_cube(gs_funcs[i], gs_funcs[j]))
            max_gs_orth_err = max(max_gs_orth_err, err)
            assert err < 1e-12

    # Check each GS function still satisfies -Delta psi = 6 psi
    h = 1e-4
    max_gs_eig_err = 0.0
    for f in gs_funcs:
        x0, y0, z0 = 0.8, 1.4, 2.2
        lap = ((f(x0 + h, y0, z0) - 2 * f(x0, y0, z0) + f(x0 - h, y0, z0))
               + (f(x0, y0 + h, z0) - 2 * f(x0, y0, z0) + f(x0, y0 - h, z0))
               + (f(x0, y0, z0 + h) - 2 * f(x0, y0, z0) + f(x0, y0, z0 - h))) / (h * h)
        err = abs(-lap - 6.0 * f(x0, y0, z0))
        max_gs_eig_err = max(max_gs_eig_err, err)
        assert err < 1e-6

    # Part (iii): Heat equation (A) with f(x,y,z) = x(a-x) y(b-y) z(d-z) (vanishes on boundary)
    kappa = 0.35
    Nmax = 9  # odd modes 1, 3, 5, 7, 9
    # Compute alpha_{p,q,r} via separable 1D quadrature
    alpha_modes = []
    for p in range(1, Nmax + 1):
        for q in range(1, Nmax + 1):
            for r in range(1, Nmax + 1):
                num = integrate_3d_separable(
                    lambda x, pp=p: x * (a - x) * math.sin(pp * math.pi * x / a),
                    lambda y, qq=q: y * (b - y) * math.sin(qq * math.pi * y / b),
                    lambda z, rr=r: z * (d - z) * math.sin(rr * math.pi * z / d),
                    a, b, d, n_panels=8,
                )
                den = (a / 2.0) * (b / 2.0) * (d / 2.0)
                an = num / den
                if abs(an) > 1e-14:
                    alpha_modes.append(((p, q, r), an, eigval((p, q, r))))

    def phi_heat(x, y, z, t):
        return sum(an * math.exp(-lam * kappa * t) * psi(m, x, y, z) for m, an, lam in alpha_modes)

    # Check IC at t=0
    x0, y0, z0 = 0.47, 0.91, 1.35
    f_exact = x0 * (a - x0) * y0 * (b - y0) * z0 * (d - z0)
    heat_ic_err = abs(phi_heat(x0, y0, z0, 0.0) - f_exact)
    assert heat_ic_err < 1e-4, f"Heat IC mismatch: {heat_ic_err}"

    # Check heat PDE at t=0.25
    t0 = 0.25
    ht = 1e-6
    hx = 1e-5
    phit = (phi_heat(x0, y0, z0, t0 + ht) - phi_heat(x0, y0, z0, t0 - ht)) / (2.0 * ht)
    lap_heat = ((phi_heat(x0 + hx, y0, z0, t0) - 2 * phi_heat(x0, y0, z0, t0) + phi_heat(x0 - hx, y0, z0, t0))
                + (phi_heat(x0, y0 + hx, z0, t0) - 2 * phi_heat(x0, y0, z0, t0) + phi_heat(x0, y0 - hx, z0, t0))
                + (phi_heat(x0, y0, z0 + hx, t0) - 2 * phi_heat(x0, y0, z0, t0) + phi_heat(x0, y0, z0 - hx, t0))) / (hx * hx)
    heat_pde_err = abs(phit - kappa * lap_heat)
    assert heat_pde_err < 1e-5, f"Heat PDE mismatch: {heat_pde_err}"

    # Check formal operator series exp(kappa*t*Delta) on a finite eigenfunction combination
    finite_modes = [((1, 1, 1), 0.8), ((2, 1, 1), -0.5), ((1, 2, 2), 0.3)]
    t_op = 0.08
    exact_heat_op = sum(c0 * math.exp(-eigval(m) * kappa * t_op) * psi(m, x0, y0, z0) for m, c0 in finite_modes)
    series_heat_op = 0.0
    for p_ord in range(40):
        # Delta^p applied to sum c0 psi_m is sum c0 (-lam_m)^p psi_m
        term = sum(c0 * ((-eigval(m)) ** p_ord) * psi(m, x0, y0, z0) for m, c0 in finite_modes)
        series_heat_op += ((kappa * t_op) ** p_ord / math.factorial(p_ord)) * term
    heat_op_err = abs(exact_heat_op - series_heat_op)
    assert heat_op_err < 1e-13, f"Heat operator series mismatch: {heat_op_err}"

    # Part (iv): Wave equation (B) with non-trivial g and h
    c_wave = 1.4
    g_modes = [((1, 1, 1), 0.9), ((2, 1, 1), -0.4), ((1, 3, 2), 0.25)]
    h_modes = [((1, 1, 1), -0.6), ((1, 2, 1), 0.7), ((2, 2, 1), 0.35)]

    all_wave_modes = {}
    for m, beta in g_modes:
        all_wave_modes.setdefault(m, [0.0, 0.0])[0] = beta
    for m, gamma in h_modes:
        all_wave_modes.setdefault(m, [0.0, 0.0])[1] = gamma

    def phi_wave(x, y, z, t):
        s = 0.0
        for m, (beta, gamma) in all_wave_modes.items():
            w = c_wave * math.sqrt(eigval(m))
            s += (beta * math.cos(w * t) + gamma * (math.sin(w * t) / w)) * psi(m, x, y, z)
        return s

    # Check wave ICs
    g_exact = sum(beta * psi(m, x0, y0, z0) for m, beta in g_modes)
    h_exact = sum(gamma * psi(m, x0, y0, z0) for m, gamma in h_modes)
    wave_ic_g_err = abs(phi_wave(x0, y0, z0, 0.0) - g_exact)
    wave_ic_h_err = abs((phi_wave(x0, y0, z0, ht) - phi_wave(x0, y0, z0, -ht)) / (2.0 * ht) - h_exact)
    assert wave_ic_g_err < 1e-15
    assert wave_ic_h_err < 1e-7

    # Check wave PDE at t = 0.31
    tw = 0.31
    phitt = (phi_wave(x0, y0, z0, tw + hx) - 2 * phi_wave(x0, y0, z0, tw) + phi_wave(x0, y0, z0, tw - hx)) / (hx * hx)
    lap_wave = ((phi_wave(x0 + hx, y0, z0, tw) - 2 * phi_wave(x0, y0, z0, tw) + phi_wave(x0 - hx, y0, z0, tw))
                + (phi_wave(x0, y0 + hx, z0, tw) - 2 * phi_wave(x0, y0, z0, tw) + phi_wave(x0, y0 - hx, z0, tw))
                + (phi_wave(x0, y0, z0 + hx, tw) - 2 * phi_wave(x0, y0, z0, tw) + phi_wave(x0, y0, z0 - hx, tw))) / (hx * hx)
    wave_pde_err = abs(phitt - (c_wave ** 2) * lap_wave)
    assert wave_pde_err < 1e-4, f"Wave PDE mismatch: {wave_pde_err}"

    # Check formal wave operator series:
    # sin(c*t*sqrt(-Delta))/(c*sqrt(-Delta)) h + cos(c*t*sqrt(-Delta)) g
    # = sum_{p>=0} (-1)^p c^{2p} t^{2p+1} / (2p+1)! (-Delta)^p h + sum_{p>=0} (-1)^p (ct)^{2p} / (2p)! (-Delta)^p g
    series_wave_op = 0.0
    for p_ord in range(30):
        sign = -1.0 if (p_ord % 2 == 1) else 1.0
        h_term = sum(gamma * (eigval(m) ** p_ord) * psi(m, x0, y0, z0) for m, gamma in h_modes)
        g_term = sum(beta * (eigval(m) ** p_ord) * psi(m, x0, y0, z0) for m, beta in g_modes)
        series_wave_op += (sign * (c_wave ** (2 * p_ord)) * (tw ** (2 * p_ord + 1)) / math.factorial(2 * p_ord + 1)) * h_term
        series_wave_op += (sign * ((c_wave * tw) ** (2 * p_ord)) / math.factorial(2 * p_ord)) * g_term
    wave_op_err = abs(phi_wave(x0, y0, z0, tw) - series_wave_op)
    assert wave_op_err < 1e-13, f"Wave operator series mismatch: {wave_op_err}"

    # Check boundary conditions on all 6 faces
    for bx, by, bz in [(0.0, y0, z0), (a, y0, z0), (x0, 0.0, z0), (x0, b, z0), (x0, y0, 0.0), (x0, y0, d)]:
        assert abs(phi_heat(bx, by, bz, t0)) < 1e-14
        assert abs(phi_wave(bx, by, bz, tw)) < 1e-14

    print("IB_methods_sheet-2_qx10_check: ALL PASSED")
    print(f"  part (i) max Rayleigh quotient error: {max_rayleigh_err:.3e}")
    print(f"  part (i) max distinct-eigenvalue orthogonality error: {max_orth_err:.3e}")
    print(f"  part (ii) degenerate E_6 Gram-Schmidt orth error: {max_gs_orth_err:.3e}, eig eq error: {max_gs_eig_err:.3e}")
    print(f"  part (iii) heat IC err: {heat_ic_err:.3e}, PDE err: {heat_pde_err:.3e}, operator series err: {heat_op_err:.3e}")
    print(f"  part (iv) wave IC err (g,h): ({wave_ic_g_err:.3e}, {wave_ic_h_err:.3e}), PDE err: {wave_pde_err:.3e}, operator series err: {wave_op_err:.3e}")


if __name__ == "__main__":
    main()
