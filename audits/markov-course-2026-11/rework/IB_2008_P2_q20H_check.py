"""Exact rational verification for IB/2008/P2/q20H (birth-death Markov chain)."""
from fractions import Fraction


def solve_linear(A, b):
    """Exact Gaussian elimination over Fraction for A x = b."""
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        pivot = None
        for r in range(col, n):
            if M[r][col] != 0:
                pivot = r
                break
        assert pivot is not None, "Singular matrix"
        M[col], M[pivot] = M[pivot], M[col]
        piv = M[col][col]
        for j in range(col, n + 1):
            M[col][j] /= piv
        for r in range(n):
            if r != col and M[r][col] != 0:
                factor = M[r][col]
                for j in range(col, n + 1):
                    M[r][j] -= factor * M[col][j]
    return [M[i][n] for i in range(n)]


def gamma_seq(p_list, q_list, N):
    """Return gamma_k = prod_{r=1}^k (q_r / p_r) for k = 0, ..., N-1."""
    gammas = [Fraction(1)]
    for r in range(1, N):
        gammas.append(gammas[-1] * q_list[r] / p_list[r])
    return gammas


def a_seq(p_list, q_list, N):
    """Return a_n = prod_{r=1}^n (p_{r-1} / q_r) for n = 0, ..., N."""
    a = [Fraction(1)]
    for r in range(1, N + 1):
        a.append(a[-1] * p_list[r - 1] / q_list[r])
    return a


def verify_truncation_hitting(name, p_fn, q_fn, N=15):
    """Verify exact u_i^(N) = P_i(H_0 < H_N) against closed-form formula."""
    p_list = [p_fn(i) for i in range(N + 1)]
    q_list = [q_fn(i) for i in range(N + 1)]
    for i in range(N + 1):
        assert p_list[i] + q_list[i] == 1
        assert p_list[i] > 0
        if i >= 1:
            assert q_list[i] > 0
        else:
            assert q_list[0] >= 0

    # Solve u_i = p_i u_{i+1} + q_i u_{i-1} for i = 1, ..., N-1 with u_0 = 1, u_N = 0
    m = N - 1
    A = [[Fraction(0)] * m for _ in range(m)]
    b = [Fraction(0)] * m
    for idx in range(m):
        i = idx + 1
        A[idx][idx] = Fraction(1)
        if i - 1 >= 1:
            A[idx][idx - 1] = -q_list[i]
        else:
            b[idx] += q_list[i] * Fraction(1)  # u_0 = 1
        if i + 1 <= N - 1:
            A[idx][idx + 1] = -p_list[i]
        # u_N = 0 contributes 0
    u_interior = solve_linear(A, b)
    u = [Fraction(1)] + u_interior + [Fraction(0)]

    gammas = gamma_seq(p_list, q_list, N)
    denom = sum(gammas)
    for i in range(N + 1):
        formula = sum(gammas[i:]) / denom
        assert u[i] == formula, f"Mismatch at i={i}: {u[i]} vs {formula}"

    # Verify a_n * gamma_n == p_0 / p_n
    a = a_seq(p_list, q_list, N)
    for n in range(1, N):
        assert a[n] * gammas[n] == p_list[0] / p_list[n]

    # Verify balance equations (pi P)_j = pi_j for j = 0, ..., N-1
    # Row 0: a_0 * q_0 + a_1 * q_1 == a_0
    assert a[0] * q_list[0] + a[1] * q_list[1] == a[0]
    for j in range(1, N):
        assert a[j - 1] * p_list[j - 1] + a[j + 1] * q_list[j + 1] == a[j]

    # Reflected finite chain on {0, ..., N} with p_{N,N} = p_N, p_{N,N-1} = q_N:
    # Its invariant distribution has pi_n = a_n / sum_{k=0}^N a_k, and E_0(T_0) = 1/pi_0 = sum_{k=0}^N a_k.
    # Verify E_0(T_0) directly from first-step hitting-time equations!
    # For i = 1, ..., N-1: m_i = 1 + p_i m_{i+1} + q_i m_{i-1} (with m_0 = 0)
    # For i = N: m_N = 1 + p_N m_N + q_N m_{N-1} => q_N m_N - q_N m_{N-1} = 1
    Am = [[Fraction(0)] * N for _ in range(N)]
    bm = [Fraction(1)] * N
    for idx in range(N):
        i = idx + 1
        if i < N:
            Am[idx][idx] = Fraction(1)
            if i - 1 >= 1:
                Am[idx][idx - 1] = -q_list[i]
            Am[idx][idx + 1] = -p_list[i]
        else:
            Am[idx][idx] = q_list[N]
            Am[idx][idx - 1] = -q_list[N]
    m_interior = solve_linear(Am, bm)  # m_1, ..., m_N
    E0_T0 = Fraction(1) + p_list[0] * m_interior[0]
    assert E0_T0 == sum(a), f"E0(T0) mismatch: {E0_T0} vs {sum(a)}"

    print(
        f"[PASS] {name} (N={N}): u_1^(N)={u[1]} (~{float(u[1]):.6f}), "
        f"reflected E_0(T_0)={E0_T0} (~{float(E0_T0):.6f})"
    )


def main():
    # 1. Constant biased walk: transient (p=2/3, q=1/3), limit u_1 = 1/2
    verify_truncation_hitting(
        "Constant transient (p=2/3, q=1/3, q0=1/3)",
        lambda i: Fraction(2, 3),
        lambda i: Fraction(1, 3),
        N=12,
    )
    # 2. Constant symmetric walk: null-recurrent (p=1/2, q=1/2)
    verify_truncation_hitting(
        "Constant null-recurrent (p=1/2, q=1/2, q0=1/2)",
        lambda i: Fraction(1, 2),
        lambda i: Fraction(1, 2),
        N=12,
    )
    # 3. Constant biased walk: positive-recurrent (p=1/3, q=2/3)
    verify_truncation_hitting(
        "Constant positive-recurrent (p=1/3, q=2/3, q0=2/3)",
        lambda i: Fraction(1, 3),
        lambda i: Fraction(2, 3),
        N=12,
    )
    # 4. Boundary case q_0 = 0 (p_0 = 1, period 2)
    verify_truncation_hitting(
        "Boundary q_0=0, p_0=1 (p_i=1/3 for i>=1)",
        lambda i: Fraction(1) if i == 0 else Fraction(1, 3),
        lambda i: Fraction(0) if i == 0 else Fraction(2, 3),
        N=12,
    )
    # 5. Inhomogeneous transient: p_r = (r+2)/(2r+2), q_r = r/(2r+2) for r>=1
    #    gamma_n = 2/((n+1)(n+2)), sum_{n>=1} gamma_n = 1, so exact h_1 = 1 - 1/(1+1) = 1/2
    verify_truncation_hitting(
        "Inhomogeneous transient (gamma_n = 2/((n+1)(n+2)))",
        lambda i: Fraction(1, 2) if i == 0 else Fraction(i + 2, 2 * i + 2),
        lambda i: Fraction(1, 2) if i == 0 else Fraction(i, 2 * i + 2),
        N=20,
    )
    # Check exact telescoping sum for example 5: sum_{k=0}^{N-1} gamma_k = 2 - 2/(N+1)
    gam5 = [Fraction(2, (k + 1) * (k + 2)) for k in range(20)]
    gam5[0] = Fraction(1)
    assert sum(gam5) == Fraction(2) - Fraction(2, 21)
    # 6. Inhomogeneous null-recurrent: p_r = (r+1)/(2r+1), q_r = r/(2r+1) for r>=1
    #    gamma_n = 1/(n+1) (harmonic -> diverges), a_n >= p_0*(n+1) (diverges)
    verify_truncation_hitting(
        "Inhomogeneous null-recurrent (gamma_n = 1/(n+1))",
        lambda i: Fraction(1, 2) if i == 0 else Fraction(i + 1, 2 * i + 1),
        lambda i: Fraction(1, 2) if i == 0 else Fraction(i, 2 * i + 1),
        N=15,
    )
    # 7. Inhomogeneous positive-recurrent: p_r = r/(2r+1), q_r = (r+1)/(2r+1) for r>=1
    #    gamma_n = n+1 (diverges), a_n = p_0 * (2n+1) / (n*(n+1)*...) ~ O(1/n) -> check p_r = r/(3r+1)
    verify_truncation_hitting(
        "Inhomogeneous positive-recurrent (p_r = r/(3r+1))",
        lambda i: Fraction(1, 4) if i == 0 else Fraction(i, 3 * i + 1),
        lambda i: Fraction(3, 4) if i == 0 else Fraction(2 * i + 1, 3 * i + 1),
        N=15,
    )
    print("All IB/2008/P2/q20H checks passed.")


if __name__ == "__main__":
    main()
