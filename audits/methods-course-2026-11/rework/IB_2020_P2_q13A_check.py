"""Independent numerical verification for IB/2020/P2/q13A.

Verifies:
1. Power-series evaluation of J_0(x) and J_0'(x) satisfies Bessel's ODE (x F')' + alpha^2 x F = 0.
2. Lommel integral identity for arbitrary alpha != beta > 0.
3. Orthogonality and norm-squared identity for zeros gamma_k of J_0.
4. Integral int_0^b r J_0(gamma_k r) dr = -(b / gamma_k) J_0'(gamma_k b).
5. Circular drum series solution z(r, t) = sum C_k J_0(gamma_k r) sin(gamma_k t) with
       C_k = -2 * b * U * J_0'(gamma_k b) / (gamma_k^2 * (J_0'(gamma_k))^2)
   satisfies the PDE, BCs, initial conditions (including pointwise convergence of z_t(r, 0)
   to the piecewise constant profile and Parseval's identity).
"""

import math


def _hankel_jn(nu: int, ax: float) -> float:
    """Hankel asymptotic expansion for J_nu(ax), accurate to < 1e-15 for ax >= 24."""
    chi = ax - (0.5 * nu + 0.25) * math.pi
    inv8 = 1.0 / (8.0 * ax)
    p = 1.0
    q = 0.0
    ak = 1.0
    nu4 = 4.0 * nu * nu
    for k in range(1, 15):
        ak *= (nu4 - (2 * k - 1) ** 2) * inv8 / k
        if k % 2 == 0:
            r = k // 2
            p += ak if r % 2 == 0 else -ak
        else:
            r = (k - 1) // 2
            q += ak if r % 2 == 0 else -ak
    return math.sqrt(2.0 / (math.pi * ax)) * (p * math.cos(chi) - q * math.sin(chi))


# Precompute 64-point half-period midpoint angles for Bessel's integral representation
# (exact to < 1e-15 for |x| <= 24 since cos(x sin th) is even pi-periodic and analytic)
_N_INT = 64
_SIN_TH = [math.sin((i + 0.5) * math.pi / _N_INT) for i in range(_N_INT)]


def j0(x: float) -> float:
    """Evaluate J_0(x) to ~1e-15 accuracy via Bessel's integral (|x| <= 24) and
    Hankel asymptotic expansion (|x| > 24)."""
    ax = abs(x)
    if ax > 24.0:
        return _hankel_jn(0, ax)
    return sum(math.cos(ax * s) for s in _SIN_TH) / _N_INT


def j0_prime(x: float) -> float:
    """Evaluate J_0'(x) = -J_1(x) to ~1e-15 accuracy via Bessel's integral (|x| <= 24)
    and Hankel asymptotic expansion (|x| > 24)."""
    ax = abs(x)
    sgn = 1.0 if x >= 0.0 else -1.0
    if ax > 24.0:
        return -sgn * _hankel_jn(1, ax)
    return -sgn * sum(s * math.sin(ax * s) for s in _SIN_TH) / _N_INT


def find_j0_zeros(count: int) -> list[float]:
    """Find first `count` positive zeros of J_0 using McMahon guess + Newton's method."""
    zeros = []
    for k in range(1, count + 1):
        # McMahon asymptotic approximation for k-th zero of J_0
        beta = (k - 0.25) * math.pi
        x = beta + 1.0 / (8.0 * beta) - 31.0 / (384.0 * beta ** 3)
        for _ in range(30):
            fx = j0(x)
            dfx = j0_prime(x)
            dx = -fx / dfx
            x += dx
            if abs(dx) < 1e-14:
                break
        assert abs(j0(x)) < 1e-12, f"Zero {k} did not converge: J_0({x})={j0(x)}"
        if zeros:
            assert x > zeros[-1] + 2.0
        zeros.append(x)
    return zeros


def simpson(f, a: float, b: float, n: int) -> float:
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * f(a + i * h)
    return s * h / 3.0


def main() -> None:
    # 1. Verify Bessel ODE (x F')' + alpha^2 x F = 0 for F(x) = J_0(alpha * x)
    h = 1e-5
    max_ode_err = 0.0
    for alpha in (0.7, 2.4, 5.1):
        for x in (0.1, 0.4, 0.85):
            # (x F')' = F' + x F''
            f_val = j0(alpha * x)
            xp = (x + h) * alpha * j0_prime(alpha * (x + h))
            xm = (x - h) * alpha * j0_prime(alpha * (x - h))
            lhs = (xp - xm) / (2.0 * h) + alpha * alpha * x * f_val
            max_ode_err = max(max_ode_err, abs(lhs))
    print(f"1. Bessel ODE residual max error:               {max_ode_err:.3e}")
    assert max_ode_err < 1e-7

    # 2. Verify Lommel integral for arbitrary positive alpha != beta
    max_lommel_err = 0.0
    for alpha, beta in ((1.2, 3.7), (2.5, 4.1), (0.8, 6.3), (3.1, 3.2)):
        num = simpson(lambda x, a=alpha, b=beta: j0(a * x) * j0(b * x) * x, 0.0, 1.0, 1000)
        ana = (beta * j0(alpha) * j0_prime(beta) - alpha * j0(beta) * j0_prime(alpha)) / (alpha * alpha - beta * beta)
        max_lommel_err = max(max_lommel_err, abs(num - ana))
    print(f"2. Lommel integral max error:                   {max_lommel_err:.3e}")
    assert max_lommel_err < 1e-12

    # 3. Compute zeros gamma_k and check orthogonality + normalization
    gammas = find_j0_zeros(500)
    print(f"   First 3 zeros of J_0: {gammas[0]:.10f}, {gammas[1]:.10f}, {gammas[2]:.10f}")
    assert abs(gammas[0] - 2.4048255577) < 1e-9
    assert abs(gammas[1] - 5.5200781103) < 1e-9

    max_orth_err = 0.0
    max_norm_err = 0.0
    for k in range(8):
        gk = gammas[k]
        for l in range(k, 8):
            gl = gammas[l]
            val = simpson(lambda x, a=gk, b=gl: j0(a * x) * j0(b * x) * x, 0.0, 1.0, 4000)
            if k == l:
                expected = 0.5 * (j0_prime(gk) ** 2)
                max_norm_err = max(max_norm_err, abs(val - expected))
            else:
                max_orth_err = max(max_orth_err, abs(val))
    print(f"3. Orthogonality max error (k!=l, <=8):         {max_orth_err:.3e}")
    print(f"   Normalization max error (k=l,  <=8):         {max_norm_err:.3e}")
    assert max_orth_err < 1e-13
    assert max_norm_err < 1e-13

    # 4. Verify int_0^b r J_0(gamma_k r) dr = -(b / gamma_k) J_0'(gamma_k b)
    max_int_err = 0.0
    for b in (0.25, 0.43, 0.78):
        for k in range(10):
            gk = gammas[k]
            num = simpson(lambda r, a=gk: r * j0(a * r), 0.0, b, 4000)
            ana = -(b / gk) * j0_prime(gk * b)
            max_int_err = max(max_int_err, abs(num - ana))
    print(f"4. Integral int_0^b r J_0(gamma_k r) dr error:  {max_int_err:.3e}")
    assert max_int_err < 1e-13

    # 5. Verify Part (ii) drum solution C_k, PDE, BCs, initial velocity & Parseval
    U = 2.7
    b = 0.45
    C = [-2.0 * b * U * j0_prime(gk * b) / (gk * gk * (j0_prime(gk) ** 2)) for gk in gammas]

    # 5a. Check PDE on 20-mode truncation via finite differences
    def z_trunc(r: float, t: float, M: int = 20) -> float:
        return sum(C[k] * j0(gammas[k] * r) * math.sin(gammas[k] * t) for k in range(M))

    h = 1e-5
    max_pde_err = 0.0
    for r in (0.15, 0.4, 0.75):
        for t in (0.2, 0.9):
            # (1/r) d/dr (r z_r) = z_rr + (1/r) z_r
            z_rr = (z_trunc(r + h, t) - 2.0 * z_trunc(r, t) + z_trunc(r - h, t)) / (h * h)
            z_r = (z_trunc(r + h, t) - z_trunc(r - h, t)) / (2.0 * h)
            z_tt = (z_trunc(r, t + h) - 2.0 * z_trunc(r, t) + z_trunc(r, t - h)) / (h * h)
            max_pde_err = max(max_pde_err, abs(z_rr + z_r / r - z_tt))
    print(f"5a. Drum PDE residual (20-mode FD):             {max_pde_err:.3e}")
    assert max_pde_err < 1e-5

    # 5b. Boundary and initial displacement conditions
    max_bc_err = max(abs(z_trunc(1.0, t, 50)) for t in (0.0, 0.3, 1.1, 2.5))
    max_ic0_err = max(abs(z_trunc(r, 0.0, 50)) for r in (0.0, 0.2, 0.5, 0.9))
    print(f"5b. z(1, t) err={max_bc_err:.2e}, z(r, 0) err={max_ic0_err:.2e}")
    assert max_bc_err < 1e-12 and max_ic0_err == 0.0

    # 5c. Initial velocity z_t(r, 0) = sum C_k gamma_k J_0(gamma_k r)
    def zt0(r: float, M: int = 500) -> float:
        return sum(C[k] * gammas[k] * j0(gammas[k] * r) for k in range(M))

    err_in = max(abs(zt0(r) - U) for r in (0.05, 0.15, 0.25, 0.35))
    err_out = max(abs(zt0(r) - 0.0) for r in (0.55, 0.65, 0.8, 0.95))
    err_jump = abs(zt0(b) - 0.5 * U)
    print(f"5c. z_t(r,0) 500-term err: r<b {err_in:.3e} | r>b {err_out:.3e} | at r=b (vs U/2) {err_jump:.3e}")
    assert err_in < 1e-2 and err_out < 1e-2 and err_jump < 1e-2

    # 5d. Parseval's identity: sum 0.5 * (C_k * gamma_k)^2 * (J_0'(gamma_k))^2 == 0.5 * U^2 * b^2
    #     (step-function tail decays as O(1/M))
    parseval_200 = sum(0.5 * ((C[k] * gammas[k]) ** 2) * (j0_prime(gammas[k]) ** 2) for k in range(200))
    parseval_500 = sum(0.5 * ((C[k] * gammas[k]) ** 2) * (j0_prime(gammas[k]) ** 2) for k in range(500))
    parseval_rhs = 0.5 * U * U * b * b
    rel_200 = abs(parseval_200 - parseval_rhs) / parseval_rhs
    rel_500 = abs(parseval_500 - parseval_rhs) / parseval_rhs
    print(
        f"5d. Parseval identity: rhs={parseval_rhs:.6f} | "
        f"200 modes err={rel_200:.3e} | 500 modes err={rel_500:.3e} (ratio={rel_200/rel_500:.2f})"
    )
    assert rel_500 < 1e-3 and 2.2 < rel_200 / rel_500 < 2.8
    print("ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
