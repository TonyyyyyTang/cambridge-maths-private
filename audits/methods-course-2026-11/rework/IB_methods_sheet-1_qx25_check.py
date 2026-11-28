#!/usr/bin/env python3
"""Verification script for IB/methods/sheet-1/qx25 (2007-08 Sheet 1, Q7).

Checks:
1. Forward FT of f(x) = e^{-x} H(x) equals (1 - ik)/(1 + k^2).
2. Symmetric inverse FT integral on [-R, R] at x = 0 equals arctan(R)/pi -> 1/2.
3. Asymmetric truncation [-A, B] fails to converge in imaginary part unless A/B -> 1.
4. Parseval's identity and general-x inversion (via tail subtraction).
"""

import cmath
import math


def simpson_complex(func, a, b, n=4000):
    assert n % 2 == 0 and n > 0
    h = (b - a) / n
    total = func(a) + func(b)
    for j in range(1, n):
        x = a + j * h
        total += (4.0 if j % 2 == 1 else 2.0) * func(x)
    return total * (h / 3.0)


def f_tilde(k):
    return (1.0 - 1j * k) / (1.0 + k * k)


def main():
    # 1. Forward Fourier transform check: int_0^{inf} e^{-x} e^{-ikx} dx = (1-ik)/(1+k^2)
    # Truncate at X = 40 (tail <= e^{-40} ~ 4.2e-18).
    max_fwd_err = 0.0
    for k in (-5.0, -1.7, -0.3, 0.0, 0.5, 1.0, 2.4, 7.5):
        num = simpson_complex(lambda x: math.exp(-x) * cmath.exp(-1j * k * x), 0.0, 40.0, n=30000)
        ex = f_tilde(k)
        err = abs(num - ex)
        max_fwd_err = max(max_fwd_err, err)
        assert err < 1e-11, f"Forward FT mismatch at k={k}: {num} vs {ex}"

    # 2. Symmetric inverse FT at x = 0: (1/(2pi)) int_{-R}^{R} f_tilde(k) dk = arctan(R)/pi
    max_sym_err = 0.0
    for R in (1.0, 3.5, 10.0, 50.0, 200.0):
        n_panels = max(4000, int(100 * R))
        if n_panels % 2 == 1:
            n_panels += 1
        val = simpson_complex(f_tilde, -R, R, n=n_panels) / (2.0 * math.pi)
        expected = math.atan(R) / math.pi
        err = abs(val - expected)
        max_sym_err = max(max_sym_err, err)
        assert err < 1e-10, f"Symmetric integral mismatch at R={R}: {val} vs {expected}"
        assert abs(val.imag) < 1e-13

    # Check limit as R -> inf approaches 1/2
    val_large = math.atan(1e9) / math.pi
    assert abs(val_large - 0.5) < 1e-9

    # 3. Demonstrate why Cauchy principal value (symmetric limit) is essential:
    # (1/(2pi)) int_{-A}^{B} f_tilde(k) dk = (atan(B)+atan(A))/(2pi) - (i/(4pi)) ln((1+B^2)/(1+A^2))
    A, B = 10.0, 30.0
    asym_num = simpson_complex(f_tilde, -A, B, n=8000) / (2.0 * math.pi)
    asym_ex = ((math.atan(B) + math.atan(A)) / (2.0 * math.pi)
               - 1j * math.log((1.0 + B * B) / (1.0 + A * A)) / (4.0 * math.pi))
    asym_err = abs(asym_num - asym_ex)
    assert asym_err < 1e-11
    assert abs(asym_num.imag) > 0.1  # Non-vanishing imaginary part when A != B

    # 4. Parseval check: int_0^inf e^{-2x} dx = 1/2 vs (1/(2pi)) int_{-inf}^{inf} |f_tilde(k)|^2 dk
    # Use substitution k = tan(t), t in (-pi/2, pi/2): |f_tilde(k)|^2 dk = dt!
    parseval_num = simpson_complex(
        lambda t: (abs(f_tilde(math.tan(t))) ** 2) / (math.cos(t) ** 2),
        -math.pi / 2.0 + 1e-12, math.pi / 2.0 - 1e-12, n=2000
    ).real / (2.0 * math.pi)
    parseval_err = abs(parseval_num - 0.5)
    assert parseval_err < 1e-11

    # 5. General-x inversion check using exact tail separation:
    # f_inv(x) = (1/pi) int_0^inf cos(kx)/(1+k^2) dk + (1/pi) int_0^inf k*sin(kx)/(1+k^2) dk
    # Write k/(1+k^2) = 1/k - 1/(k(1+k^2)), and use int_0^inf sin(kx)/k dk = (pi/2)*sgn(x).
    # Integrate remaining O(1/k^2) and O(1/k^3) integrands via k = tan(t) on (0, pi/2)
    # (Wait: oscillatory near pi/2, so integrate directly to K=2000 for smooth check).
    max_gen_err = 0.0
    for x in (-2.0, -0.5, 0.0, 0.5, 1.5, 3.0):
        if x == 0.0:
            rec = 0.5
            exact = 0.5
        else:
            # First term: (1/pi) int_0^inf cos(kx)/(1+k^2) dk = (1/2) e^{-|x|}
            # Second term: (1/2)*sgn(x) - (1/pi) int_0^inf sin(kx)/(k(1+k^2)) dk
            # Integrate directly on [0, K] with K = 800 (tail < 1/(pi*800^2) ~ 5e-7 for O(1/k^3),
            # and for cos(kx)/(1+k^2) integrate by parts once:
            # int_0^inf cos(kx)/(1+k^2) dk = int_0^inf 2k*sin(kx)/(x*(1+k^2)^2) dk, which is O(1/k^3)!)
            K = 400.0
            N = 40000
            term1 = (1.0 / math.pi) * simpson_complex(
                lambda k: 2.0 * k * math.sin(k * x) / (x * (1.0 + k * k) ** 2), 0.0, K, n=N
            ).real
            def g2(k):
                if k == 0.0:
                    return x
                return math.sin(k * x) / (k * (1.0 + k * k))
            term2 = 0.5 * (1.0 if x > 0 else -1.0) - (1.0 / math.pi) * simpson_complex(g2, 0.0, K, n=N).real
            rec = term1 + term2
            exact = math.exp(-x) if x > 0 else 0.0
        err = abs(rec - exact)
        max_gen_err = max(max_gen_err, err)
        assert err < 1e-5, f"General-x inversion mismatch at x={x}: {rec} vs {exact}"

    print("IB_methods_sheet-1_qx25_check.py: ALL CHECKS PASSED")
    print(f"  Max forward FT quadrature error:                {max_fwd_err:.3e}")
    print(f"  Max symmetric finite-R integral error (x=0):    {max_sym_err:.3e}")
    print(f"  Asymmetric [-10, 30] analytical-vs-num error:   {asym_err:.3e} (Im={asym_num.imag:.6f})")
    print(f"  Parseval check error:                           {parseval_err:.3e}")
    print(f"  Max general-x inversion error (x in [-2, 3]):   {max_gen_err:.3e}")


if __name__ == "__main__":
    main()
