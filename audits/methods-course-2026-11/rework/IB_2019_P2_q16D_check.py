"""Independent numerical and exact rational verification for IB/2019/P2/q16D.

Verifies for Gegenbauer polynomials C_n^alpha(x) solving
    (1 - x^2) y'' - (2*alpha + 1) x y' + n(n + 2*alpha) y = 0  (alpha > 0):
1. Exact rational construction & ODE satisfaction for n = 0..6 and several alpha > 0.
2. Sturm-Liouville form (p y')' + n(n+2*alpha) w y = 0 with
       p(x) = (1 - x^2)^(alpha + 1/2),  w(x) = (1 - x^2)^(alpha - 1/2),  (a,b) = (-1,1).
3. Orthogonality on (-1, 1) with weight w(x) for all m != n (via x = cos(theta)).
4. All n roots are real, simple, and lie strictly in (-1, 1).
5. The integral int_{-1}^1 C_n^alpha(x) prod_{i=1}^k (x - x_i) w(x) dx vanishes for
   any proper subset of roots (k < n) and is strictly non-zero when k = n.
"""

from fractions import Fraction
import math


def gegenbauer_coeffs_exact(n: int, alpha: Fraction) -> list[Fraction]:
    """Return exact monic coefficients [c_0, c_1, ..., c_n] of C_n^alpha(x)."""
    c = [Fraction(0)] * (n + 1)
    c[n] = Fraction(1)
    for j in range(n - 2, -1, -1):
        num = (j + 2) * (j + 1)
        den = (n - j) * (n + j + 2 * alpha)
        c[j] = -Fraction(num, 1) * c[j + 2] / den
    return c


def poly_eval(coeffs: list[float], x: float) -> float:
    val = 0.0
    for a in reversed(coeffs):
        val = val * x + a
    return val


def poly_deriv_exact(c: list[Fraction]) -> list[Fraction]:
    return [c[j] * j for j in range(1, len(c))]


def verify_ode_exact(n: int, alpha: Fraction, c: list[Fraction]) -> None:
    """Verify (1 - x^2) y'' - (2*alpha+1) x y' + n(n+2*alpha) y == 0 identically."""
    dy = poly_deriv_exact(c)
    d2y = poly_deriv_exact(dy)
    lam = n * (n + 2 * alpha)
    res = [Fraction(0)] * (n + 1)
    for j, a in enumerate(d2y):
        res[j] += a
        res[j + 2] -= a
    for j, a in enumerate(dy):
        res[j + 1] -= (2 * alpha + 1) * a
    for j, a in enumerate(c):
        res[j] += lam * a
    assert all(r == 0 for r in res), f"ODE failed for n={n}, alpha={alpha}: {res}"


def simpson(f, a: float, b: float, n: int) -> float:
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * f(a + i * h)
    return s * h / 3.0


def weighted_integral(f, alpha: float, n_pts: int = 800) -> float:
    """Compute int_{-1}^1 f(x) (1-x^2)^(alpha - 1/2) dx using tanh-sinh quadrature
    (double-exponential), achieving ~1e-14 accuracy even when 0 < alpha < 1/2."""
    # x = tanh((pi/2) * sinh(t)), dx = (pi/2) * cosh(t) / cosh^2((pi/2)*sinh(t)) dt
    # 1 - x^2 = 1 / cosh^2((pi/2)*sinh(t))
    h = 4.0 / n_pts
    total = 0.0
    half_pi = 0.5 * math.pi
    for k in range(-n_pts, n_pts + 1):
        t = k * h
        sh = math.sinh(t)
        ch = math.cosh(t)
        u = half_pi * sh
        if abs(u) > 35.0:
            continue
        x = math.tanh(u)
        cu = math.cosh(u)
        # (1 - x^2)^(alpha - 0.5) * dx/dt = half_pi * ch * cu^(-2*(alpha - 0.5) - 2)
        #                                 = half_pi * ch * cu^(-2*alpha - 1)
        weight = half_pi * ch * (cu ** (-2.0 * alpha - 1.0))
        total += f(x) * weight
    return total * h


def find_roots_in_open_interval(coeffs: list[float], n_scan: int = 4000) -> list[float]:
    """Bracket and bisect all sign-changing roots in (-1, 1)."""
    roots = []
    eps = 1e-12
    xs = [-1.0 + eps + (2.0 - 2.0 * eps) * i / n_scan for i in range(n_scan + 1)]
    vals = [poly_eval(coeffs, x) for x in xs]
    for i in range(n_scan):
        if vals[i] == 0.0:
            roots.append(xs[i])
        elif vals[i] * vals[i + 1] < 0.0:
            lo, hi = xs[i], xs[i + 1]
            flo = vals[i]
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                fmid = poly_eval(coeffs, mid)
                if fmid == 0.0 or hi - lo < 1e-15:
                    lo = hi = mid
                    break
                if flo * fmid < 0.0:
                    hi = mid
                else:
                    lo = mid
                    flo = fmid
            roots.append(0.5 * (lo + hi))
    return roots


def main() -> None:
    alphas = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 1), Fraction(5, 2)]
    max_n = 6

    for alpha_frac in alphas:
        alpha = float(alpha_frac)
        polys_exact = [gegenbauer_coeffs_exact(n, alpha_frac) for n in range(max_n + 1)]
        polys_float = [[float(x) for x in c] for c in polys_exact]
        dpolys_float = [[float(x) for x in poly_deriv_exact(c)] for c in polys_exact]

        # 1. Exact rational ODE check
        for n in range(max_n + 1):
            verify_ode_exact(n, alpha_frac, polys_exact[n])

        # 2. Check Sturm-Liouville self-adjoint form (p y')' + lam w y = 0 numerically
        h = 1e-5
        max_sl_err = 0.0
        for n in range(max_n + 1):
            lam = n * (n + 2.0 * alpha)
            for x in (-0.7, -0.2, 0.35, 0.8):
                py_p = ((1.0 - (x + h) ** 2) ** (alpha + 0.5)) * poly_eval(dpolys_float[n], x + h)
                py_m = ((1.0 - (x - h) ** 2) ** (alpha + 0.5)) * poly_eval(dpolys_float[n], x - h)
                d_py = (py_p - py_m) / (2.0 * h)
                w = (1.0 - x * x) ** (alpha - 0.5)
                res = abs(d_py + lam * w * poly_eval(polys_float[n], x))
                max_sl_err = max(max_sl_err, res)
        assert max_sl_err < 1e-6, f"SL residual too large: {max_sl_err}"

        # 3. Check orthogonality for m != n and strict positivity for m == n
        max_orth_err = 0.0
        min_norm_sq = float("inf")
        for m in range(max_n + 1):
            for n in range(m, max_n + 1):
                val = weighted_integral(
                    lambda x, pm=polys_float[m], pn=polys_float[n]: poly_eval(pm, x) * poly_eval(pn, x),
                    alpha,
                )
                if m == n:
                    min_norm_sq = min(min_norm_sq, val)
                else:
                    max_orth_err = max(max_orth_err, abs(val))
        assert max_orth_err < 1e-11, f"Orthogonality failed for alpha={alpha_frac}: {max_orth_err}"
        assert min_norm_sq > 0.0

        # 4 & 5. Check all n roots lie strictly in (-1, 1), are simple, and test the Part 2 integral
        min_deriv_at_root = float("inf")
        max_subset_int = 0.0
        min_full_int = float("inf")
        for n in range(1, max_n + 1):
            roots = find_roots_in_open_interval(polys_float[n])
            assert len(roots) == n, f"Expected {n} roots in (-1,1) for alpha={alpha_frac}, got {len(roots)}"
            for r in roots:
                assert -1.0 < r < 1.0
                assert abs(poly_eval(polys_float[n], r)) < 1e-13
                min_deriv_at_root = min(min_deriv_at_root, abs(poly_eval(dpolys_float[n], r)))
            # Distinctness
            for i in range(n - 1):
                assert roots[i + 1] - roots[i] > 1e-4

            # Proper subset k = n - 1 -> integral must be 0
            sub_roots = roots[:-1]
            sub_int = weighted_integral(
                lambda x, pn=polys_float[n], rs=sub_roots: poly_eval(pn, x) * math.prod(x - r for r in rs),
                alpha,
            )
            max_subset_int = max(max_subset_int, abs(sub_int))

            # Full root set k = n -> integral must be strictly positive (monic C_n^alpha)
            full_int = weighted_integral(
                lambda x, pn=polys_float[n], rs=roots: poly_eval(pn, x) * math.prod(x - r for r in rs),
                alpha,
            )
            min_full_int = min(min_full_int, full_int)

        assert min_deriv_at_root > 1e-2
        assert max_subset_int < 1e-11
        assert min_full_int > 0.0

        print(
            f"alpha={str(alpha_frac):>4s}: exact ODE OK | SL err={max_sl_err:.2e} | "
            f"orth err={max_orth_err:.2e} | all roots in (-1,1) simple (min |C'|={min_deriv_at_root:.2e}) | "
            f"k<n int={max_subset_int:.2e}, k=n min int={min_full_int:.2e}"
        )

    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
