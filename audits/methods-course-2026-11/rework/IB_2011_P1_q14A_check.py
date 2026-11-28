"""Computational verification for IB/2011/P1/q14A (Parseval's theorem & square wave Fourier series)."""

import math


def gauss_legendre_nodes_weights(n):
    """Compute n-point Gauss-Legendre nodes and weights on [-1, 1]."""
    nodes = [0.0] * n
    weights = [0.0] * n
    m = (n + 1) // 2
    for i in range(m):
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


def integrate_gl(func, a, b, n=64):
    nodes, weights = gauss_legendre_nodes_weights(n)
    mid = 0.5 * (a + b)
    half = 0.5 * (b - a)
    return half * sum(w * func(mid + half * x) for x, w in zip(nodes, weights))


def test_parseval_general_T():
    """Test Parseval's theorem on (-T, T) for f(t) = t^2 + 2t - 1 with T = 2.5."""
    T = 2.5
    f = lambda t: t * t + 2.0 * t - 1.0
    # Exact integral (1/T) * int_{-T}^T (t^2 + 2t - 1)^2 dt
    # = (2/T) * [T^5/5 + 2*T^3/3 + T] = 2*(T^4/5 + 2*T^2/3 + 1)
    lhs_exact = 2.0 * (T ** 4 / 5.0 + 2.0 * T * T / 3.0 + 1.0)
    lhs_num = integrate_gl(lambda t: f(t) ** 2, -T, T, n=20) / T
    assert abs(lhs_num - lhs_exact) < 1e-12, f"lhs_num={lhs_num}, lhs_exact={lhs_exact}"

    # Exact Fourier coefficients of f(t) = t^2 + 2t - 1 on (-T, T):
    # a_0 = 2*(T^2/3 - 1)
    # a_n = 4*T^2*(-1)^n / (n*pi)^2
    # b_n = -4*T*(-1)^n / (n*pi)
    a0_num = integrate_gl(f, -T, T, n=20) / T
    a0_exact = 2.0 * (T * T / 3.0 - 1.0)
    assert abs(a0_num - a0_exact) < 1e-13

    for n in range(1, 15):
        an_num = integrate_gl(lambda t, k=n: f(t) * math.cos(k * math.pi * t / T), -T, T, n=40) / T
        bn_num = integrate_gl(lambda t, k=n: f(t) * math.sin(k * math.pi * t / T), -T, T, n=40) / T
        an_ex = 4.0 * T * T * ((-1) ** n) / ((n * math.pi) ** 2)
        bn_ex = -4.0 * T * ((-1) ** n) / (n * math.pi)
        assert abs(an_num - an_ex) < 1e-12
        assert abs(bn_num - bn_ex) < 1e-12

    # Sum Parseval RHS using exact formula up to N = 50000 (tail ~ 16*T^2 / (pi^2 * N))
    N = 50000
    rhs = 0.5 * a0_exact ** 2 + sum(
        (4.0 * T * T / ((n * math.pi) ** 2)) ** 2 + (4.0 * T / (n * math.pi)) ** 2
        for n in range(1, N + 1)
    )
    # Add leading tail correction for 1/n^2: sum_{n>N} 1/n^2 ~ 1/(N + 0.5)
    rhs_corr = rhs + (16.0 * T * T / (math.pi ** 2)) / (N + 0.5)
    assert abs(rhs_corr - lhs_exact) < 1e-8
    print(f"Parseval general T={T}: LHS={lhs_exact:.12f}, RHS(corr)={rhs_corr:.12f}")


def test_square_wave_fourier():
    """Verify Fourier coefficients, pointwise convergence, Gibbs overshoot, and exact sums."""
    # Numerical coefficients for g(t) = 0 on (-pi, 0], 1 on (0, pi)
    a0 = integrate_gl(lambda t: 1.0, 0.0, math.pi, n=20) / math.pi
    assert abs(a0 - 1.0) < 1e-14
    for n in range(1, 20):
        an = integrate_gl(lambda t, k=n: math.cos(k * t), 0.0, math.pi, n=40) / math.pi
        bn = integrate_gl(lambda t, k=n: math.sin(k * t), 0.0, math.pi, n=40) / math.pi
        assert abs(an) < 1e-14
        bn_exact = (2.0 / (n * math.pi)) if (n % 2 == 1) else 0.0
        assert abs(bn - bn_exact) < 1e-14

    def S_N(t, N):
        return 0.5 + (2.0 / math.pi) * sum(
            math.sin((2 * k + 1) * t) / (2 * k + 1) for k in range(N)
        )

    # Check exact value at t = 0 for all N
    for N in (1, 5, 50, 500):
        assert S_N(0.0, N) == 0.5

    # Check pointwise convergence away from jump
    assert abs(S_N(0.5 * math.pi, 10000) - 1.0) < 1e-4
    assert abs(S_N(-0.5 * math.pi, 10000) - 0.0) < 1e-4
    assert abs(S_N(1.0, 10000) - 1.0) < 1e-4
    assert abs(S_N(-1.0, 10000) - 0.0) < 1e-4

    # Check Gibbs overshoot near t = 0+ (first peak at t = pi / (2N))
    N_gibbs = 500
    overshoot = S_N(math.pi / (2.0 * N_gibbs), N_gibbs)
    # Theoretical Wilbraham-Gibbs limit: 0.5 + (1/pi) * Si(pi) ~ 1.0894898722
    assert 1.089 < overshoot < 1.090

    # Check Leibniz sum: sum_{k=0}^inf (-1)^k / (2k+1) = pi / 4
    # Averaging consecutive partial sums dramatically improves accuracy
    M = 20000
    s_even = sum(((-1) ** k) / (2 * k + 1) for k in range(M))
    s_odd = s_even + ((-1) ** M) / (2 * M + 1)
    s_avg = 0.5 * (s_even + s_odd)
    assert abs(s_avg - math.pi / 4.0) < 1e-9

    # Check Parseval sum: sum_{k=0}^inf 1 / (2k+1)^2 = pi^2 / 8
    # Tail: sum_{k=M}^inf 1/(2k+1)^2 = int_{M-1/2}^inf dx/(2x+1)^2 + O(M^-3) = 1/(4M)
    M2 = 50000
    s2 = sum(1.0 / ((2 * k + 1) ** 2) for k in range(M2))
    s2_corr = s2 + 1.0 / (4.0 * M2)
    assert abs(s2_corr - (math.pi ** 2) / 8.0) < 1e-10

    # Verify Parseval identity directly for g(t)
    lhs_g = 1.0  # (1/pi) * int_0^pi 1^2 dt
    rhs_g = 0.5 * (1.0 ** 2) + (4.0 / (math.pi ** 2)) * ((math.pi ** 2) / 8.0)
    assert abs(lhs_g - rhs_g) < 1e-15

    print(
        f"Square wave checks passed: S_N(0)=0.5, Gibbs overshoot={overshoot:.6f}, "
        f"Leibniz={s_avg:.10f} (pi/4={math.pi/4:.10f}), "
        f"odd sq sum={s2_corr:.12f} (pi^2/8={math.pi**2/8:.12f})"
    )


if __name__ == "__main__":
    test_parseval_general_T()
    test_square_wave_fourier()
    print("ALL CHECKS PASSED for IB_2011_P1_q14A")
