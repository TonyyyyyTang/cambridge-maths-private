#!/usr/bin/env python3
"""Computational verification for IB/2015/P1/q14C (Sturm-Liouville self-adjointness and Fredholm integral equation)."""

import math

GL_NODES = [
    -0.96028985649753623168,
    -0.79666647741362673959,
    -0.52553240991632898582,
    -0.18343464249564980494,
    0.18343464249564980494,
    0.52553240991632898582,
    0.79666647741362673959,
    0.96028985649753623168,
]
GL_WEIGHTS = [
    0.10122853629037625915,
    0.22238103445337447054,
    0.31370664587788728734,
    0.36268378337836198297,
    0.36268378337836198297,
    0.31370664587788728734,
    0.22238103445337447054,
    0.10122853629037625915,
]


def quad(func, a, b, panels=40):
    h = (b - a) / panels
    total = 0.0
    for j in range(panels):
        x0 = a + j * h
        mid = x0 + 0.5 * h
        half = 0.5 * h
        for t, w in zip(GL_NODES, GL_WEIGHTS):
            total += w * func(mid + half * t)
    return 0.5 * h * total


def test_part_ii_self_adjointness():
    """Verify Green's identity and vanishing of boundary terms for 5 BC types."""
    # General regular operator L y = (p y')' - q y on [0, 1]
    p = lambda x: 1.0 + 0.5 * math.sin(math.pi * x)  # p(0) = p(1) = 1
    dp = lambda x: 0.5 * math.pi * math.cos(math.pi * x)
    q = lambda x: 2.0 + x * x

    def L(y, dy, d2y, x):
        return p(x) * d2y(x) + dp(x) * dy(x) - q(x) * y(x)

    def check_pair(y1, dy1, d2y1, y2, dy2, d2y2, a=0.0, b=1.0, p_fn=p, dp_fn=dp, q_fn=q):
        def Lloc(y, dy, d2y, x):
            return p_fn(x) * d2y(x) + dp_fn(x) * dy(x) - q_fn(x) * y(x)

        lhs = quad(lambda x: y1(x) * Lloc(y2, dy2, d2y2, x) - y2(x) * Lloc(y1, dy1, d2y1, x), a, b, panels=40)
        bdy = (
            p_fn(b) * (y1(b) * dy2(b) - dy1(b) * y2(b))
            - p_fn(a) * (y1(a) * dy2(a) - dy1(a) * y2(a))
        )
        assert abs(lhs - bdy) < 1e-13, f"Green's identity failed: lhs={lhs}, bdy={bdy}"
        return lhs, bdy

    # 1. Dirichlet: y(0)=y(1)=0
    lhs, bdy = check_pair(
        lambda x: math.sin(math.pi * x),
        lambda x: math.pi * math.cos(math.pi * x),
        lambda x: -(math.pi ** 2) * math.sin(math.pi * x),
        lambda x: x * (1.0 - x) * math.exp(x),
        lambda x: (1.0 - x - x * x) * math.exp(x),
        lambda x: (-3.0 * x - x * x) * math.exp(x),
    )
    assert abs(bdy) < 1e-14

    # 2. Neumann: y'(0)=y'(1)=0
    lhs, bdy = check_pair(
        lambda x: math.cos(math.pi * x),
        lambda x: -math.pi * math.sin(math.pi * x),
        lambda x: -(math.pi ** 2) * math.cos(math.pi * x),
        lambda x: 3.0 * x ** 2 - 2.0 * x ** 3,
        lambda x: 6.0 * x - 6.0 * x ** 2,
        lambda x: 6.0 - 12.0 * x,
    )
    assert abs(bdy) < 1e-14

    # 3. Robin: y'(0) + 2 y(0) = 0, 3 y'(1) - y(1) = 0
    # Construct polynomials satisfying exact Robin BCs:
    # u(x) = 1 - 2x + c x^2 => u'(0)+2u(0)=0; 3(-2+2c) - (1-2+c) = 0 => 5c - 5 = 0 => c = 1
    # v(x) = 1 - 2x + d x^3 => v'(0)+2v(0)=0; 3(-2+3d) - (1-2+d) = 0 => 8d - 5 = 0 => d = 5/8
    lhs, bdy = check_pair(
        lambda x: 1.0 - 2.0 * x + x ** 2,
        lambda x: -2.0 + 2.0 * x,
        lambda x: 2.0,
        lambda x: 1.0 - 2.0 * x + 0.625 * x ** 3,
        lambda x: -2.0 + 1.875 * x ** 2,
        lambda x: 3.75 * x,
    )
    assert abs(bdy) < 1e-14

    # 4. Periodic (with p(0) = p(1) = 1): y(0)=y(1), y'(0)=y'(1)
    lhs, bdy = check_pair(
        lambda x: math.cos(2.0 * math.pi * x) + 0.3 * math.sin(4.0 * math.pi * x),
        lambda x: -2.0 * math.pi * math.sin(2.0 * math.pi * x) + 1.2 * math.pi * math.cos(4.0 * math.pi * x),
        lambda x: -4.0 * math.pi ** 2 * math.cos(2.0 * math.pi * x) - 4.8 * math.pi ** 2 * math.sin(4.0 * math.pi * x),
        lambda x: math.sin(2.0 * math.pi * x) + 0.7,
        lambda x: 2.0 * math.pi * math.cos(2.0 * math.pi * x),
        lambda x: -4.0 * math.pi ** 2 * math.sin(2.0 * math.pi * x),
    )
    assert abs(bdy) < 1e-13

    # 5. Singular (Legendre: p(x) = 1 - x^2 vanishes at a=-1, b=1)
    lhs, bdy = check_pair(
        lambda x: x ** 2,
        lambda x: 2.0 * x,
        lambda x: 2.0,
        lambda x: math.exp(x),
        lambda x: math.exp(x),
        lambda x: math.exp(x),
        a=-1.0,
        b=1.0,
        p_fn=lambda x: 1.0 - x * x,
        dp_fn=lambda x: -2.0 * x,
        q_fn=lambda x: 0.0,
    )
    assert abs(bdy) < 1e-14


def test_part_iii_general_fredholm():
    """Verify part (iii) for non-constant p(x)=1+x, q(x)=0, w(x)=1+x^2, general Robin BCs, non-zero f and lambda."""
    a, b = 0.0, 1.0
    p = lambda x: 1.0 + x
    dp = lambda x: 1.0
    w = lambda x: 1.0 + x * x
    lam = -1.7

    # Robin BCs: alpha y'(0) + beta y(0) = 0, mu y'(1) + nu y(1) = 0
    alpha, beta = 1.0, 2.0
    mu, nu = 3.0, -1.0

    # Homogeneous solutions of ((1+x) u')' = 0 are 1 and ln(1+x).
    # u1(x) = alpha * ln(1+x) - beta  => u1'(0) = alpha, u1(0) = -beta => alpha*u1'(0)+beta*u1(0) = 0? Wait:
    # u1(x) = -beta * ln(1+x) + alpha => u1'(0) = -beta, u1(0) = alpha => alpha*(-beta) + beta*alpha = 0!
    u1 = lambda x: -beta * math.log(1.0 + x) + alpha
    du1 = lambda x: -beta / (1.0 + x)

    # u2(x) satisfies mu u2'(1) + nu u2(1) = 0:
    # Let u2(x) = c * ln(1+x) + d => u2'(1) = c/2, u2(1) = c*ln2 + d
    # Choose c = -nu, then mu*(-nu/2) + nu*(-nu*ln2 + d) = 0 => d = mu/2 + nu*ln2
    ln2 = math.log(2.0)
    u2 = lambda x: -nu * math.log(1.0 + x) + (0.5 * mu + nu * ln2)
    du2 = lambda x: -nu / (1.0 + x)

    assert abs(alpha * du1(a) + beta * u1(a)) < 1e-15
    assert abs(mu * du2(b) + nu * u2(b)) < 1e-15

    # Wronskian p(x) * (u1(x) u2'(x) - u1'(x) u2(x)) is constant C_W:
    C_W = p(0.3) * (u1(0.3) * du2(0.3) - du1(0.3) * u2(0.3))
    assert abs(C_W - p(0.8) * (u1(0.8) * du2(0.8) - du1(0.8) * u2(0.8))) < 1e-14

    # Green's function for ((1+x) G')' = delta(x - xi):
    # G(x, xi) = u1(x) u2(xi) / C_W for x <= xi, u1(xi) u2(x) / C_W for x >= xi
    def G(x, xi):
        if x <= xi:
            return u1(x) * u2(xi) / C_W
        return u1(xi) * u2(x) / C_W

    # Verify symmetry and jump condition p(xi)*(G_x(xi+) - G_x(xi-)) = 1
    for xi in [0.2, 0.55, 0.85]:
        assert abs(G(0.3, xi) - G(xi, 0.3)) < 1e-14
        jump = p(xi) * (u1(xi) * du2(xi) / C_W - du1(xi) * u2(xi) / C_W)
        assert abs(jump - 1.0) < 1e-14

    def K(x, xi):
        return math.sqrt(w(xi) * w(x)) * G(x, xi)

    # Pick an exact y(x) satisfying the same Robin BCs:
    # y(x) = 1 - 2x + x^2 + 0.5*(1 - 2x + 0.625 x^3)
    y = lambda x: (1.0 - 2.0 * x + x ** 2) + 0.5 * (1.0 - 2.0 * x + 0.625 * x ** 3)
    dy = lambda x: (-2.0 + 2.0 * x) + 0.5 * (-2.0 + 1.875 * x ** 2)
    d2y = lambda x: 2.0 + 0.5 * (3.75 * x)
    assert abs(alpha * dy(a) + beta * y(a)) < 1e-14
    assert abs(mu * dy(b) + nu * y(b)) < 1e-14

    # Compute forcing f(x) = (1/w) (p y')' - lam y
    f = lambda x: (p(x) * d2y(x) + dp(x) * dy(x)) / w(x) - lam * y(x)
    psi = lambda x: math.sqrt(w(x)) * y(x)

    # Check psi(xi) == U(xi) + lam * int_a^b K(x, xi) psi(x) dx
    # splitting the integral at x = xi for exact quadrature accuracy
    max_fred_err = 0.0
    for xi in [0.1, 0.25, 0.5, 0.73, 0.9]:
        U_xi = quad(lambda x: K(x, xi) * math.sqrt(w(x)) * f(x), a, xi, panels=30) + quad(
            lambda x: K(x, xi) * math.sqrt(w(x)) * f(x), xi, b, panels=30
        )
        I_xi = quad(lambda x: K(x, xi) * psi(x), a, xi, panels=30) + quad(
            lambda x: K(x, xi) * psi(x), xi, b, panels=30
        )
        rhs = U_xi + lam * I_xi
        max_fred_err = max(max_fred_err, abs(psi(xi) - rhs))
    assert max_fred_err < 1e-13, f"Part (iii) Fredholm error: {max_fred_err}"
    return max_fred_err


def test_part_iv_dirichlet_fredholm():
    """Verify part (iv) Green's function and Fredholm integral equation for y'' - lam y = 0, y(0)=y(1)=0."""
    def G(x, xi):
        if x <= xi:
            return -x * (1.0 - xi)
        return -xi * (1.0 - x)

    # Check BCs, continuity, symmetry, jump
    for xi in [0.15, 0.4, 0.8]:
        assert abs(G(0.0, xi)) < 1e-15
        assert abs(G(1.0, xi)) < 1e-15
        assert abs(G(0.3, xi) - G(xi, 0.3)) < 1e-15
        h = 1e-7
        gx_minus = (G(xi, xi) - G(xi - h, xi)) / h
        gx_plus = (G(xi + h, xi) - G(xi, xi)) / h
        assert abs((gx_plus - gx_minus) - 1.0) < 1e-9

    # Check eigenfunctions y_n(x) = sin(n pi x) with lambda_n = -(n pi)^2
    max_eig_err = 0.0
    for n in range(1, 6):
        lam = -((n * math.pi) ** 2)
        yn = lambda x, n=n: math.sin(n * math.pi * x)
        for xi in [0.1, 0.3, 0.5, 0.7, 0.9]:
            rhs = -lam * (
                (1.0 - xi) * quad(lambda x: x * yn(x), 0.0, xi, panels=30)
                + xi * quad(lambda x: (1.0 - x) * yn(x), xi, 1.0, panels=30)
            )
            max_eig_err = max(max_eig_err, abs(yn(xi) - rhs))
    assert max_eig_err < 1e-13, f"Part (iv) eigenfunction Fredholm error: {max_eig_err}"

    # Independent check: power iteration on discretised K gives smallest |lambda_1| = pi^2
    M = 200
    dx = 1.0 / M
    xs = [(i + 0.5) * dx for i in range(M)]
    v = [math.sin(math.pi * x) + 0.2 * math.sin(3.0 * math.pi * x) for x in xs]
    mu = 0.0
    for _ in range(30):
        w_vec = [sum(G(xs[j], xs[i]) * v[j] * dx for j in range(M)) for i in range(M)]
        norm = math.sqrt(sum(val * val for val in w_vec))
        # Rayleigh quotient: v^T K v / v^T v = 1 / lambda_1
        rq = sum(v[i] * w_vec[i] for i in range(M)) / sum(v[i] * v[i] for i in range(M))
        mu = 1.0 / rq
        v = [val / norm for val in w_vec]
    assert abs(mu - (-(math.pi ** 2))) < 1e-3, f"Discretised eigenvalue mismatch: {mu}"
    return max_eig_err, mu


def main():
    test_part_ii_self_adjointness()
    err_iii = test_part_iii_general_fredholm()
    err_iv, lam1_num = test_part_iv_dirichlet_fredholm()
    print("IB_2015_P1_q14C verification PASSED:")
    print("  Part (ii): Green's identity & vanishing boundary term verified for all 5 BC types")
    print(f"  Part (iii): general Fredholm equation (non-const p, w, Robin BCs, f!=0) max error = {err_iii:.3e}")
    print(f"  Part (iv): eigenfunction check (n=1..5) max error = {err_iv:.3e}")
    print(f"  Part (iv): discretised Fredholm operator lambda_1 = {lam1_num:.6f} (exact -pi^2 = {-(math.pi**2):.6f})")


if __name__ == "__main__":
    main()
