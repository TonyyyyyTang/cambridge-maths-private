"""Computational verification for IB/methods/sheet-1/q10."""
import math


def bessel_j_and_derivs(m, z):
    """Compute J_m(z), J_m'(z), J_m''(z).
    Uses Bessel's integral (periodic trapezoidal rule, exponentially accurate without
    alternating-series cancellation) for integer m >= 0, and exact/series for non-integer m > -1.
    """
    if isinstance(m, int) and m >= 0:
        # J_m(z) = (1/pi) int_0^pi cos(m*t - z*sin(t)) dt
        # Even 2*pi-periodic integrand -> composite trapezoidal rule on [0, pi] is spectrally accurate.
        N = max(64, int(2.5 * z) + 48)
        h = math.pi / N
        s0 = 0.5 * (math.cos(0.0) + math.cos(m * math.pi))
        s1 = 0.0  # sin(0)=sin(pi)=0
        s2 = 0.0
        for k in range(1, N):
            t = k * h
            st = math.sin(t)
            ang = m * t - z * st
            ca = math.cos(ang)
            sa = math.sin(ang)
            s0 += ca
            s1 += st * sa
            s2 -= st * st * ca
        inv_pi = h / math.pi
        return s0 * inv_pi, s1 * inv_pi, s2 * inv_pi

    # Non-integer m > -1 (used only for moderate z <= 10)
    if z == 0.0:
        raise ValueError("Derivative singular at z=0 for non-integer m < 2")
    half_z = 0.5 * z
    term = (half_z ** m) / math.gamma(m + 1.0)
    val = term
    dval = term * (m / z)
    ddval = term * (m * (m - 1.0) / (z * z))
    for k in range(1, 120):
        term *= -(half_z * half_z) / (k * (k + m))
        p = 2 * k + m
        val += term
        dval += term * (p / z)
        ddval += term * (p * (p - 1.0) / (z * z))
        if abs(term) < 1e-17 * max(1.0, abs(val)) and abs(term * p / z) < 1e-17 * max(1.0, abs(dval)):
            break
    return val, dval, ddval


def simpson(f, a, b, n=1000):
    assert n % 2 == 0
    h = (b - a) / n
    s = f(a) + f(b)
    for k in range(1, n):
        s += (4 if k % 2 == 1 else 2) * f(a + k * h)
    return s * h / 3.0


def find_positive_zeros(m, count=5):
    zeros = []
    step = 0.25
    x = 0.1
    prev_val = bessel_j_and_derivs(m, x)[0]
    while len(zeros) < count:
        nx = x + step
        cur_val = bessel_j_and_derivs(m, nx)[0]
        if prev_val * cur_val <= 0.0:
            lo, hi = x, nx
            flo = prev_val
            for _ in range(70):
                mid = 0.5 * (lo + hi)
                fmid = bessel_j_and_derivs(m, mid)[0]
                if flo * fmid <= 0.0:
                    hi = mid
                else:
                    lo = mid
                    flo = fmid
            root = 0.5 * (lo + hi)
            # Refine with Newton steps
            for _ in range(5):
                fv, fp, _ = bessel_j_and_derivs(m, root)
                root -= fv / fp
            zeros.append(root)
        x = nx
        prev_val = cur_val
    return zeros


def test_ode_and_wronskian():
    max_ode_err = 0.0
    max_wron_err = 0.0
    h = 1e-6
    for m in (0, 1, 2, 3, 0.5):
        for alpha, beta in ((1.7, 3.4), (2.5, 5.1), (4.2, 1.3)):
            for i in range(1, 20):
                r = i / 20.0
                ja, jpa, jppa = bessel_j_and_derivs(m, alpha * r)
                ya = ja
                yap = alpha * jpa
                yapp = alpha * alpha * jppa
                # -(r y')' + (m^2/r) y = alpha^2 r y
                lhs = -(yap + r * yapp) + (m * m / r) * ya
                rhs = alpha * alpha * r * ya
                max_ode_err = max(max_ode_err, abs(lhs - rhs))

                def bracket(rr):
                    va, dpa, _ = bessel_j_and_derivs(m, alpha * rr)
                    vb, dpb, _ = bessel_j_and_derivs(m, beta * rr)
                    return rr * (va * beta * dpb - vb * alpha * dpa)

                d_bracket = (bracket(r + h) - bracket(r - h)) / (2.0 * h)
                vb = bessel_j_and_derivs(m, beta * r)[0]
                rhs_w = (alpha * alpha - beta * beta) * r * ya * vb
                max_wron_err = max(max_wron_err, abs(d_bracket - rhs_w))

    assert max_ode_err < 1e-11, f"ODE error: {max_ode_err}"
    assert max_wron_err < 1e-7, f"Wronskian derivative error: {max_wron_err}"
    print(f"ODE & Wronskian identity passed: max_ode_err={max_ode_err:.2e}, max_wron_err={max_wron_err:.2e}.")


def test_lommel_integral():
    max_err = 0.0
    for m in (0, 1, 2, 3, 0.5):
        for alpha, beta in ((1.3, 4.7), (2.8, 6.2), (5.5, 2.1)):
            def integrand(r):
                if r == 0.0:
                    return 0.0
                return bessel_j_and_derivs(m, alpha * r)[0] * bessel_j_and_derivs(m, beta * r)[0] * r

            num = simpson(integrand, 0.0, 1.0, n=1200)
            ja, jpa, _ = bessel_j_and_derivs(m, alpha)
            jb, jpb, _ = bessel_j_and_derivs(m, beta)
            exact = (beta * ja * jpb - alpha * jb * jpa) / (alpha * alpha - beta * beta)
            max_err = max(max_err, abs(num - exact))

    assert max_err < 1e-10, f"Lommel integral mismatch: {max_err}"
    print(f"General alpha != beta integral check passed (max_err = {max_err:.2e}).")


def test_orthogonality_zeros():
    max_err = 0.0
    for m in (0, 1, 2, 3):
        zeros = find_positive_zeros(m, count=5)
        for k, jk in enumerate(zeros):
            assert abs(bessel_j_and_derivs(m, jk)[0]) < 1e-13
            for l, jl in enumerate(zeros):
                def integrand(r, jk=jk, jl=jl, m=m):
                    if r == 0.0:
                        return 0.0
                    return bessel_j_and_derivs(m, jk * r)[0] * bessel_j_and_derivs(m, jl * r)[0] * r

                num = simpson(integrand, 0.0, 1.0, n=1600)
                if k == l:
                    jp = bessel_j_and_derivs(m, jk)[1]
                    exact = 0.5 * jp * jp
                else:
                    exact = 0.0
                max_err = max(max_err, abs(num - exact))
        print(f"  m={m}: first 5 positive zeros = {[round(z, 6) for z in zeros]}")

    assert max_err < 1e-10, f"Orthogonality mismatch: {max_err}"
    print(f"Zero-orthogonality & norm checks passed for m=0,1,2,3 (max_err = {max_err:.2e}).")

    # Sanity check on edge case m=1, alpha=0: integral is 0 whereas (1/2)[J_1'(0)]^2 = 1/8 != 0
    j1_0, j1p_0, _ = bessel_j_and_derivs(1, 0.0)
    assert abs(j1_0) < 1e-15 and abs(0.5 * j1p_0 * j1p_0 - 0.125) < 1e-14
    print("Confirmed edge case: j_{mk} must be positive zeros (at z=0 for m=1, J_1(0)=0 but (1/2)[J_1'(0)]^2 = 1/8 != 0).")


if __name__ == "__main__":
    test_ode_and_wronskian()
    test_lommel_integral()
    test_orthogonality_zeros()
    print("ALL CHECKS PASSED FOR IB/methods/sheet-1/q10.")
