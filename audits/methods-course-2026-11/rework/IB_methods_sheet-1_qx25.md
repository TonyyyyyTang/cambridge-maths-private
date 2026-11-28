---
id: IB/methods/sheet-1/qx25
slug: IB_methods_sheet-1_qx25
source_archive: private/raw/sheets/IB/methods/2007-08/sheet-1.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-1_qx25_p2.png` (page 2 of `2007-08/sheet-1.pdf`, Question 7). Question 8 immediately follows on page 2, so there is no page spillover.

- **Question 7 body:** Exact verbatim match to the printed sheet, including lowercase `fourier transform` / `fourier series`, wording (`takes the value of 1/2 at x = 0`), formula `\tilde f(k)=\frac{1-ik}{1+k^2}.`, and final parenthetical remark.
- **Purely typographic differences:** Printed title is underlined italic *Functions with discontinuities.*; stored `question.tex` uses `\emph{Functions with discontinuities.}`.
- **Sheet convention context:** Page 2's section preamble above Question 6 defines $\tilde{f}(k) = \int_{-\infty}^{\infty} f(x) e^{-ikx}\,dx$, which is also uniquely confirmed here by the target formula $\tilde{f}(k) = \frac{1-ik}{1+k^2}$.

## Full Solution

We use the Fourier transform pair from the sheet's convention:
\[
\tilde{f}(k) = \int_{-\infty}^{\infty} f(x) e^{-ikx}\,dx, \qquad
f_{\mathrm{inv}}(x) = \lim_{R \to \infty} \frac{1}{2\pi} \int_{-R}^{R} \tilde{f}(k) e^{ikx}\,dk.
\]

### Part 1: Computing $\tilde{f}(k)$

Here $f(x) = e^{-x}$ for $x > 0$ and $f(x) = 0$ for $x < 0$ (the single point $x = 0$ has measure zero and does not affect the integral). Since $f \in L^1(\mathbb{R})$, its Fourier transform exists for every $k \in \mathbb{R}$:
\[
\tilde{f}(k) = \int_{0}^{\infty} e^{-x} e^{-ikx}\,dx = \int_{0}^{\infty} e^{-(1 + ik)x}\,dx = \lim_{X \to \infty} \left[ \frac{e^{-(1 + ik)x}}{-(1 + ik)} \right]_{0}^{X}.
\]
For every real $k$, $\operatorname{Re}(1 + ik) = 1 > 0$, so
\[
\left|e^{-(1 + ik)X}\right| = e^{-X} \to 0 \quad \text{as } X \to \infty.
\]
Hence
\[
\tilde{f}(k) = \frac{1}{1 + ik}.
\]
Multiplying numerator and denominator by the complex conjugate $1 - ik$ gives
\[
\tilde{f}(k) = \frac{1 - ik}{(1 + ik)(1 - ik)} = \frac{1 - ik}{1 + k^2},
\]
as required.

### Part 2: Inverse Fourier transform at $x = 0$

At $x = 0$, $e^{ikx} = 1$, so the Fourier inversion integral is
\[
\mathcal{F}^{-1}[\tilde{f}](0) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \frac{1 - ik}{1 + k^2}\,dk.
\]
**Convergence remark (Cauchy principal value):** Note that $\tilde{f}(k) \notin L^1(\mathbb{R})$ because the imaginary part $-\dfrac{k}{1 + k^2} \sim -\dfrac{1}{k}$ as $|k| \to \infty$ is not absolutely integrable. Indeed, on an asymmetric interval $[-A, B]$,
\[
\int_{-A}^{B} \frac{-ik}{1 + k^2}\,dk = -\frac{i}{2}\ln\!\left(\frac{1 + B^2}{1 + A^2}\right),
\]
which has no limit as $A, B \to \infty$ independently. As in the Fourier inversion theorem for piecewise smooth functions with jump discontinuities (analogous to symmetric partial sums $\sum_{n=-N}^{N}$ of Fourier series), the inversion integral must be interpreted as the symmetric (Cauchy principal value) limit $R \to \infty$ over $[-R, R]$:
\[
\mathcal{F}^{-1}[\tilde{f}](0) = \lim_{R \to \infty} \frac{1}{2\pi} \int_{-R}^{R} \frac{1 - ik}{1 + k^2}\,dk
= \lim_{R \to \infty} \frac{1}{2\pi}\left(\int_{-R}^{R} \frac{dk}{1 + k^2} - i \int_{-R}^{R} \frac{k\,dk}{1 + k^2}\right).
\]
For every finite $R > 0$:
1. The integrand $\dfrac{k}{1 + k^2}$ is continuous and **odd** on the symmetric interval $[-R, R]$, so
   \[
   \int_{-R}^{R} \frac{k}{1 + k^2}\,dk = \left[\frac{1}{2}\ln(1 + k^2)\right]_{-R}^{R} = 0.
   \]
2. The integrand $\dfrac{1}{1 + k^2}$ is even and integrable on $\mathbb{R}$, with
   \[
   \int_{-R}^{R} \frac{dk}{1 + k^2} = \Big[\arctan k\Big]_{-R}^{R} = 2\arctan R.
   \]

Taking the limit $R \to \infty$, since $\arctan R \to \dfrac{\pi}{2}$:
\[
\mathcal{F}^{-1}[\tilde{f}](0) = \lim_{R \to \infty} \frac{2\arctan R}{2\pi} = \lim_{R \to \infty} \frac{\arctan R}{\pi} = \frac{\pi/2}{\pi} = \frac{1}{2}.
\]
Thus the inverse Fourier transform of $\tilde{f}(k)$ equals $\frac{1}{2}$ at $x = 0$, which is $\frac{1}{2}\bigl(f(0^+) + f(0^-)\bigr) = \frac{1}{2}(1 + 0)$, the midpoint of the jump discontinuity at $x = 0$.

## Final Answers

- Forward Fourier transform:
  \[
  \tilde{f}(k) = \frac{1}{1 + ik} = \frac{1 - ik}{1 + k^2}.
  \]
- Inverse Fourier transform at $x = 0$ (in the symmetric / Cauchy principal value sense):
  \[
  \lim_{R \to \infty} \frac{1}{2\pi} \int_{-R}^{R} \frac{1 - ik}{1 + k^2}\,dk = \lim_{R \to \infty} \frac{\arctan R}{\pi} = \frac{1}{2}.
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_qx25_check.py` (executed with `python3`, exit code 0).

```text
IB_methods_sheet-1_qx25_check.py: ALL CHECKS PASSED
  Max forward FT quadrature error:                7.606e-12
  Max symmetric finite-R integral error (x=0):    1.277e-15
  Asymmetric [-10, 30] analytical-vs-num error:   2.388e-15 (Im=-0.174146)
  Parseval check error:                           3.184e-13
  Max general-x inversion error (x in [-2, 3]):   2.357e-08
```

## Confidence & Edge-Case Notes

- **Confidence:** High.
- **Phase B watchpoints:**
  1. **Principal value / symmetric cutoff rigor:** A naive write-up that simply declares $\int_{-\infty}^\infty \frac{k}{1+k^2}\,dk = 0$ "because it is odd" without noting that the integral is not absolutely convergent and must be taken as the symmetric limit $\lim_{R\to\infty}\int_{-R}^R$ is a classic rigor gap. Check whether the existing solution addresses the symmetric limit properly.
  2. **Inverse factor $\frac{1}{2\pi}$:** Verify consistency with the sheet's convention $\tilde{f}(k) = \int_{-\infty}^\infty f(x)e^{-ikx}\,dx$ and $\mathcal{F}^{-1}[\tilde{f}](x) = \frac{1}{2\pi}\int_{-\infty}^\infty \tilde{f}(k)e^{ikx}\,dk$.
  3. **Justification of the upper limit at $+\infty$:** Verify that $\lim_{X \to \infty} e^{-(1+ik)X} = 0$ is justified via $|e^{-(1+ik)X}| = e^{-X} \to 0$.
