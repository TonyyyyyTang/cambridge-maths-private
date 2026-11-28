"""Computational verification for IB/2007/P2/q15D (Sturm-Liouville solvability & bifurcation)."""
import math


def simpson_grid(vals, h):
    n = len(vals) - 1
    assert n % 2 == 0
    s = vals[0] + vals[-1]
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * vals[i]
    return s * h / 3.0


def rk4_shoot(p, dp, q, w, lam, slope0, nonlinear=False, n_steps=400):
    """Integrate (p y')' + (q + lam w) y = (y^3 if nonlinear else 0) on [0, 1] with y(0)=0, y'(0)=slope0."""
    h = 1.0 / n_steps
    y = 0.0
    v = slope0  # v = y'
    ys = [y]

    def rhs(x, yy, vv):
        px = p(x)
        dpx = dp(x)
        src = (yy ** 3) if nonlinear else 0.0
        dv = (src - (q(x) + lam * w(x)) * yy - dpx * vv) / px
        return vv, dv

    for i in range(n_steps):
        x = i * h
        k1y, k1v = rhs(x, y, v)
        k2y, k2v = rhs(x + 0.5 * h, y + 0.5 * h * k1y, v + 0.5 * h * k1v)
        k3y, k3v = rhs(x + 0.5 * h, y + 0.5 * h * k2y, v + 0.5 * h * k2v)
        k4y, k4v = rhs(x + h, y + h * k3y, v + h * k3v)
        y += (h / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
        v += (h / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)
        ys.append(y)
    return ys


def secant_root(func, x0, x1, tol=1e-14, max_iter=60):
    f0 = func(x0)
    f1 = func(x1)
    for _ in range(max_iter):
        if abs(f1) < tol:
            return x1
        denom = f1 - f0
        assert abs(denom) > 0
        x2 = x1 - f1 * (x1 - x0) / denom
        x0, f0 = x1, f1
        x1, f1 = x2, func(x2)
        if abs(x1 - x0) < tol:
            return x1
    raise RuntimeError("Secant method failed to converge")


def run_case(name, p, dp, q, w, lam_guess0, lam_guess1, n_steps=800):
    h = 1.0 / n_steps
    xs = [i * h for i in range(n_steps + 1)]

    # Find ground-state eigenvalue lam0 of linear problem
    lam0 = secant_root(
        lambda lam: rk4_shoot(p, dp, q, w, lam, 1.0, nonlinear=False, n_steps=n_steps)[-1],
        lam_guess0,
        lam_guess1,
    )
    y0_raw = rk4_shoot(p, dp, q, w, lam0, 1.0, nonlinear=False, n_steps=n_steps)
    norm_sq = simpson_grid([w(x) * y * y for x, y in zip(xs, y0_raw)], h)
    scale = 1.0 / math.sqrt(norm_sq)
    y0 = [scale * y for y in y0_raw]
    y0_slope = scale  # y0'(0)

    # Check normalization and y0(1) = 0
    assert abs(simpson_grid([w(x) * y * y for x, y in zip(xs, y0)], h) - 1.0) < 1e-13
    assert abs(y0[-1]) < 1e-12

    # Test solvability condition (*) on a smooth test function y(x) vanishing at 0 and 1
    def y_test(x):
        return x * (1.0 - x) * (1.0 + 0.7 * math.sin(3.0 * x) + 0.4 * x * x)

    def dy_test(x):
        u = x * (1.0 - x)
        du = 1.0 - 2.0 * x
        v = 1.0 + 0.7 * math.sin(3.0 * x) + 0.4 * x * x
        dv = 2.1 * math.cos(3.0 * x) + 0.8 * x
        return du * v + u * dv

    def d2y_test(x):
        u = x * (1.0 - x)
        du = 1.0 - 2.0 * x
        d2u = -2.0
        v = 1.0 + 0.7 * math.sin(3.0 * x) + 0.4 * x * x
        dv = 2.1 * math.cos(3.0 * x) + 0.8 * x
        d2v = -6.3 * math.sin(3.0 * x) + 0.8
        return d2u * v + 2.0 * du * dv + u * d2v

    f_vals = [
        p(x) * d2y_test(x) + dp(x) * dy_test(x) + (q(x) + lam0 * w(x)) * y_test(x)
        for x in xs
    ]
    solv_integral = simpson_grid([y0_i * f_i for y0_i, f_i in zip(y0, f_vals)], h)
    assert abs(solv_integral) < 1e-10, f"Solvability integral failed: {solv_integral}"

    # Compute predicted mu0 = int_0^1 y0^4 dx
    mu0 = simpson_grid([y ** 4 for y in y0], h)

    # Solve nonlinear BVP L(y; lam) = y^3, y(0) = y(1) = 0 for decreasing epsilon
    epsilons = [0.08, 0.04, 0.02, 0.01]
    errs = []
    for eps in epsilons:
        # Shoot with y'(0) = eps * y0'(0) so y = eps * y0 + O(eps^3)
        lam_nl = secant_root(
            lambda lam: rk4_shoot(p, dp, q, w, lam, eps * y0_slope, nonlinear=True, n_steps=n_steps)[-1],
            lam0 + eps * eps * mu0 * 0.9,
            lam0 + eps * eps * mu0 * 1.1,
        )
        y_nl = rk4_shoot(p, dp, q, w, lam_nl, eps * y0_slope, nonlinear=True, n_steps=n_steps)
        # Check exact identity (lam - lam0) * int w y0 y = int y0 y^3
        lhs_exact = (lam_nl - lam0) * simpson_grid(
            [w(x) * a * b for x, a, b in zip(xs, y0, y_nl)], h
        )
        rhs_exact = simpson_grid([a * (b ** 3) for a, b in zip(y0, y_nl)], h)
        assert abs(lhs_exact - rhs_exact) < 1e-11

        mu_num = (lam_nl - lam0) / (eps * eps)
        err = abs(mu_num - mu0)
        errs.append((eps, mu_num, err))

        # Also test non-orthogonal re-parameterisation eps_tilde = eps / (1 + c * eps)
        # so y = eps_tilde * y0 + eps_tilde^2 * (c * y0) + O(eps_tilde^3) -> O(eps_tilde) error!
        c = 0.5
        eps_tilde = eps / (1.0 + c * eps)
        mu_tilde = (lam_nl - lam0) / (eps_tilde * eps_tilde)
        expected_shift = 2.0 * c * eps_tilde * mu0
        actual_shift = mu_tilde - mu0
        assert abs(actual_shift - expected_shift) < 5.0 * (eps_tilde ** 2) * mu0 + 1e-6

    # Check convergence to mu0 and O(eps^2) scaling when y1 = O(eps)
    assert errs[-1][2] < 1e-4, f"mu did not converge to mu0: {errs}"
    ratio = errs[1][2] / errs[2][2]
    assert 3.5 < ratio < 4.5, f"Expected O(eps^2) ratio ~4, got {ratio}"

    print(f"Case [{name}]: lam0={lam0:.10f}, int(y0^4)={mu0:.10f}, solv_int={solv_integral:.2e}")
    for eps, mu_num, err in errs:
        print(f"  eps={eps:.3f}: mu_num={mu_num:.10f}, |mu_num - mu0|={err:.3e}")
    print(f"  error halving ratio (eps 0.04->0.02): {ratio:.3f} (O(eps^2) when y1=O(eps))")
    return lam0, mu0


def main():
    # Case 1: constant coefficients y'' + lam y = y^3 -> lam0 = pi^2, int y0^4 = 3/2
    lam0, mu0 = run_case(
        "p=1, q=0, w=1",
        p=lambda x: 1.0,
        dp=lambda x: 0.0,
        q=lambda x: 0.0,
        w=lambda x: 1.0,
        lam_guess0=9.5,
        lam_guess1=10.2,
    )
    assert abs(lam0 - math.pi ** 2) < 1e-10
    assert abs(mu0 - 1.5) < 1e-10

    # Case 2: variable coefficients p(x)=1+x, q(x)=-x^2, w(x)=1+2x
    run_case(
        "p=1+x, q=-x^2, w=1+2x",
        p=lambda x: 1.0 + x,
        dp=lambda x: 1.0,
        q=lambda x: -x * x,
        w=lambda x: 1.0 + 2.0 * x,
        lam_guess0=6.0,
        lam_guess1=8.5,
    )
    print("IB_2007_P2_q15D verification PASSED")


if __name__ == "__main__":
    main()
