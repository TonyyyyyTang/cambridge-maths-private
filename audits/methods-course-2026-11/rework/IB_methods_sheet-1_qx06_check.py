"""Computational verification for IB/methods/sheet-1/qx06."""
import math


def j0_all(z):
    """Compute J_0(z), J_0'(z), J_0''(z) via Bessel's integral (spectrally accurate for all z >= 0)."""
    # J_0(z) = (1/pi) int_0^pi cos(z sin t) dt (even about pi/2 as well -> integrate on [0, pi/2])
    N = max(48, int(1.5 * abs(z)) + 32)
    h = (0.5 * math.pi) / N
    # Composite trapezoidal rule on [0, pi/2] (f'(0)=0, f'(pi/2)=0, symmetric periodic -> exponential accuracy)
    s0 = 0.5 * (1.0 + math.cos(z))
    s1 = 0.5 * (-math.sin(z))
    s2 = 0.5 * (-math.cos(z))
    for k in range(1, N):
        t = k * h
        st = math.sin(t)
        ang = z * st
        ca = math.cos(ang)
        sa = math.sin(ang)
        s0 += ca
        s1 -= st * sa
        s2 -= st * st * ca
    scale = (2.0 / math.pi) * h
    return s0 * scale, s1 * scale, s2 * scale


def i0(z):
    """Modified Bessel function I_0(z) for checking lambda < 0."""
    half_z = 0.5 * z
    term = 1.0
    val = 1.0
    for k in range(1, 120):
        term *= (half_z * half_z) / (k * k)
        val += term
        if term < 1e-17 * val:
            break
    return val


def simpson(f, a, b, n=1200):
    assert n % 2 == 0
    h = (b - a) / n
    s = f(a) + f(b)
    for k in range(1, n):
        s += (4 if k % 2 == 1 else 2) * f(a + k * h)
    return s * h / 3.0


def find_j0_zeros(count=20):
    zeros = []
    step = 0.25
    x = 0.1
    prev = j0_all(x)[0]
    while len(zeros) < count:
        nx = x + step
        cur = j0_all(nx)[0]
        if prev * cur <= 0.0:
            lo, hi = x, nx
            flo = prev
            for _ in range(70):
                mid = 0.5 * (lo + hi)
                fmid = j0_all(mid)[0]
                if flo * fmid <= 0.0:
                    hi = mid
                else:
                    lo = mid
                    flo = fmid
            root = 0.5 * (lo + hi)
            for _ in range(5):
                fv, fp, _ = j0_all(root)
                root -= fv / fp
            zeros.append(root)
        x = nx
        prev = cur
    return zeros


def test_part_i(zeros):
    # Verify ODE (x u')' + j_n^2 x u = 0, boundedness at 0, u(1) = 0
    max_ode_err = 0.0
    max_bc_err = 0.0
    for jn in zeros[:8]:
        u0 = j0_all(0.0)[0]
        u1 = j0_all(jn)[0]
        assert u0 == 1.0
        max_bc_err = max(max_bc_err, abs(u1))
        for i in range(1, 20):
            x = i / 20.0
            v, dp, dpp = j0_all(jn * x)
            # (x u')' + jn^2 x u = u' + x u'' + jn^2 x u
            res = jn * dp + x * jn * jn * dpp + jn * jn * x * v
            max_ode_err = max(max_ode_err, abs(res))
    assert max_bc_err < 1e-13 and max_ode_err < 1e-10

    # For lambda = -kappa^2 < 0, bounded solution is I_0(kappa x) which satisfies I_0(kappa) > 1 > 0
    for kappa in (0.1, 1.0, 3.0, 10.0):
        assert i0(kappa) > 1.0

    print(f"Part (i) checks passed: first 5 zeros = {[round(z, 6) for z in zeros[:5]]}, max_ode_err={max_ode_err:.2e}.")


def test_part_ii(zeros):
    # Check general alpha != beta formula
    max_gen_err = 0.0
    for alpha, beta in ((1.2, 4.3), (2.7, 7.1), (5.0, 2.3)):
        num = simpson(lambda x: j0_all(alpha * x)[0] * j0_all(beta * x)[0] * x, 0.0, 1.0, n=1200)
        ja, jpa, _ = j0_all(alpha)
        jb, jpb, _ = j0_all(beta)
        exact = (beta * ja * jpb - alpha * jb * jpa) / (alpha * alpha - beta * beta)
        max_gen_err = max(max_gen_err, abs(num - exact))
    assert max_gen_err < 1e-11

    # Check orthogonality & norm for n, m = 1..6
    max_orth_err = 0.0
    for n in range(6):
        jn = zeros[n]
        for m in range(6):
            jm = zeros[m]
            num = simpson(lambda x, jn=jn, jm=jm: j0_all(jn * x)[0] * j0_all(jm * x)[0] * x, 0.0, 1.0, n=1600)
            if n == m:
                jp = j0_all(jn)[1]
                exact = 0.5 * jp * jp
            else:
                exact = 0.0
            max_orth_err = max(max_orth_err, abs(num - exact))
    assert max_orth_err < 1e-10

    # Verify beta = j_n + eps limit explicitly
    jn = zeros[0]
    jp = j0_all(jn)[1]
    exact_norm = 0.5 * jp * jp
    eps = 1e-6
    beta = jn + eps
    ja, jpa, _ = j0_all(jn)
    jb, jpb, _ = j0_all(beta)
    lim_approx = (beta * ja * jpb - jn * jb * jpa) / (jn * jn - beta * beta)
    assert abs(lim_approx - exact_norm) < 1e-6

    print(f"Part (ii) checks passed: max_gen_err={max_gen_err:.2e}, max_orth_err={max_orth_err:.2e}, eps-limit error={abs(lim_approx - exact_norm):.2e}.")


def test_part_iii(zeros):
    # Exact test solution u(x) = (1 - x^2)^3 (bounded at 0, u(1) = 0).
    # (x u')' = x(1 - x^2)(-12 + 36 x^2).
    # So f(x) = (1 - x^2)(-12 + 36 x^2) + lam_tilde * (1 - x^2)^3 satisfies f bounded at 0 and f(1) = 0!
    def u_exact(x):
        return (1.0 - x * x) ** 3

    for lam_tilde in (3.5, 15.0, -4.0):
        # Verify lam_tilde is not an eigenvalue
        for jn in zeros:
            assert abs(lam_tilde - jn * jn) > 1.0

        def f(x, lt=lam_tilde):
            t = 1.0 - x * x
            return t * (-12.0 + 36.0 * x * x) + lt * (t ** 3)

        assert abs(f(1.0)) < 1e-15

        # Compute Fourier-Bessel coefficients c_n and a_n = c_n / (lam_tilde - j_n^2)
        a_coeffs = []
        for jn in zeros:
            jp = j0_all(jn)[1]
            cn = (2.0 / (jp * jp)) * simpson(lambda x, jn=jn: f(x) * j0_all(jn * x)[0] * x, 0.0, 1.0, n=2000)
            a_coeffs.append(cn / (lam_tilde - jn * jn))

        # Check reconstruction error on [0, 1]
        max_rec_err = 0.0
        for i in range(0, 41):
            x = i / 40.0
            u_approx = sum(a * j0_all(jn * x)[0] for a, jn in zip(a_coeffs, zeros))
            max_rec_err = max(max_rec_err, abs(u_approx - u_exact(x)))

        assert max_rec_err < 1e-5, f"Expansion error at lam_tilde={lam_tilde}: {max_rec_err}"
        print(f"  Part (iii) test (lam_tilde={lam_tilde:5.1f}, 20 modes): max |u_20(x) - u_exact(x)| = {max_rec_err:.2e}")


if __name__ == "__main__":
    zeros = find_j0_zeros(count=20)
    test_part_i(zeros)
    test_part_ii(zeros)
    test_part_iii(zeros)
    print("ALL CHECKS PASSED FOR IB/methods/sheet-1/qx06.")
