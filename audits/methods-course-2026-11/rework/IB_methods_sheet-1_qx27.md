---
id: IB/methods/sheet-1/qx27
slug: IB_methods_sheet-1_qx27
source_archive: private/raw/sheets/IB/methods/2007-08/sheet-1.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-1_qx27_p2.png` (page 2 of `2007-08/sheet-1.pdf`, Question 11). Page 2 is the last page of the PDF (`IB_methods_sheet-1_qx27_p3.png` does not exist; the page ends with the sheet footnote directly below Question 11).

- **Question 11 body:** Exact verbatim match to the printed sheet, including:
  - lowercase `fourier transform` / `fourier transforms`,
  - exact piecewise definition punctuation (`e^k-e^{-k},&|k|\le1,\\ 0&|k|>1,`, with no comma after `0` and a comma after `|k|>1`),
  - exact target formula `f(x)=\frac{2i}{\pi(1+x^2)}(\cosh 1\sin x-x\cos x\sinh 1).`,
  - exact PDE line `\nabla^2\psi=0;\quad-\infty<x<\infty,\quad 0<y<1,`,
  - exact wording `where $\psi(x,0)=f(x)$ the function given above, and $\psi(x,1)=0$ for $-\infty<x<\infty$.`,
  - exact historical note `(This was a long tripos question (2004/4/II/15A) for Complex Methods on material now in the Methods schedule.)`.
- **Purely typographic differences:** Printed title is underlined italic *Laplace’s equation.* vs `\emph{Laplace’s equation.}`; the final parenthetical remark starts on a new line without extra vertical space on the sheet vs a paragraph break in `question.tex`.
- **Sheet convention context:** Page 2's section preamble defines $\tilde{f}(k) = \int_{-\infty}^{\infty} f(x)e^{-ikx}\,dx$, which corresponds to the inverse transform $f(x) = \frac{1}{2\pi}\int_{-\infty}^{\infty}\tilde{f}(k)e^{ikx}\,dk$ (also uniquely confirmed by the given formula for $f(x)$).

## Full Solution

### Part 1: Inverse Fourier transform of $\tilde{f}(k)$

Using the inverse Fourier transform corresponding to the sheet's convention $\tilde{f}(k) = \int_{-\infty}^{\infty} f(x)e^{-ikx}\,dx$,
\[
f(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{f}(k) e^{ikx}\,dk
= \frac{1}{2\pi} \int_{-1}^{1} \left(e^k - e^{-k}\right) e^{ikx}\,dk.
\]
Split the integral into two exponentials:
\[
f(x) = \frac{1}{2\pi} \left( \int_{-1}^{1} e^{(1 + ix)k}\,dk - \int_{-1}^{1} e^{(-1 + ix)k}\,dk \right).
\]
For all $x \in \mathbb{R}$, both $1 + ix \neq 0$ and $-1 + ix \neq 0$, so direct integration gives
\[
\int_{-1}^{1} e^{(1 + ix)k}\,dk = \left[\frac{e^{(1 + ix)k}}{1 + ix}\right]_{-1}^{1} = \frac{e^{1 + ix} - e^{-1 - ix}}{1 + ix},
\]
and (either directly or by substituting $u = -k$)
\[
\int_{-1}^{1} e^{(-1 + ix)k}\,dk = \left[\frac{e^{(-1 + ix)k}}{-1 + ix}\right]_{-1}^{1} = \frac{e^{-1 + ix} - e^{1 - ix}}{-1 + ix} = \frac{e^{1 - ix} - e^{-1 + ix}}{1 - ix}.
\]
Writing $J(x) = \dfrac{e^{1 + ix} - e^{-1 - ix}}{1 + ix}$, note that for real $x$ the second integral is simply $\overline{J(x)}$ (equivalently $J(-x)$), so
\[
f(x) = \frac{1}{2\pi}\bigl(J(x) - \overline{J(x)}\bigr) = \frac{i}{\pi}\operatorname{Im} J(x).
\]
Now simplify $J(x)$ by multiplying numerator and denominator by $1 - ix$ and using $e^{\pm(1+ix)} = e^{\pm 1}(\cos x \pm i\sin x)$:
\[
e^{1 + ix} - e^{-1 - ix} = (e - e^{-1})\cos x + i(e + e^{-1})\sin x = 2\sinh 1 \cos x + 2i\cosh 1 \sin x.
\]
Hence
\[
\begin{aligned}
J(x) &= \frac{(1 - ix)(2\sinh 1 \cos x + 2i\cosh 1 \sin x)}{1 + x^2} \\
&= \frac{2(\sinh 1 \cos x + x\cosh 1 \sin x) + 2i(\cosh 1 \sin x - x\cos x \sinh 1)}{1 + x^2}.
\end{aligned}
\]
Taking the imaginary part and multiplying by $\dfrac{i}{\pi}$ (equivalently computing $\frac{1}{2\pi}(J(x) - J(-x))$ directly) yields
\[
f(x) = \frac{2i}{\pi(1 + x^2)}\bigl(\cosh 1 \sin x - x\cos x \sinh 1\bigr),
\]
as required.

*(Alternative check via integration by parts: since $e^k - e^{-k} = 2\sinh k$ is odd, $f(x) = \frac{2i}{\pi}\int_0^1 \sinh k \sin(kx)\,dk$, and integrating by parts twice gives the same result.)*

### Part 2: Laplace's equation in the infinite strip $0 \le y \le 1$

We solve
\[
\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial^2 \psi}{\partial y^2} = 0, \qquad -\infty < x < \infty, \quad 0 < y < 1,
\]
with boundary conditions
\[
\psi(x, 0) = f(x), \qquad \psi(x, 1) = 0, \qquad -\infty < x < \infty.
\]
To ensure uniqueness (excluding exponentially growing homogeneous solutions such as $\sin(n\pi y)e^{\pm n\pi x}$ that vanish at $y = 0, 1$) and applicability of the Fourier transform in $x$, we seek the solution that decays as $|x| \to \infty$ (with $\psi(x, y), \psi_x(x, y) \to 0$ as $|x| \to \infty$ for each $y \in [0, 1]$).

Define the Fourier transform in $x$:
\[
\tilde{\psi}(k, y) = \int_{-\infty}^{\infty} \psi(x, y) e^{-ikx}\,dx.
\]
Transforming Laplace's equation with respect to $x$ and integrating by parts twice (using $\psi, \psi_x \to 0$ as $|x| \to \infty$) gives
\[
-k^2 \tilde{\psi}(k, y) + \frac{\partial^2 \tilde{\psi}}{\partial y^2}(k, y) = 0, \qquad 0 < y < 1,
\]
subject to the transformed boundary conditions
\[
\tilde{\psi}(k, 0) = \tilde{f}(k) = \begin{cases}
e^k - e^{-k} = 2\sinh k, & |k| \le 1, \\
0, & |k| > 1,
\end{cases}
\qquad
\tilde{\psi}(k, 1) = 0.
\]

We solve this second-order linear ODE in $y$ for each fixed $k \in \mathbb{R}$:
- **Case $k \neq 0$:** The general solution can be written conveniently as
  \[
  \tilde{\psi}(k, y) = A(k)\sinh(ky) + B(k)\sinh\bigl(k(1 - y)\bigr).
  \]
  Applying the boundary condition at $y = 1$:
  \[
  0 = \tilde{\psi}(k, 1) = A(k)\sinh k \implies A(k) = 0 \quad (\text{since } \sinh k \neq 0 \text{ for } k \neq 0).
  \]
  Applying the boundary condition at $y = 0$:
  \[
  \tilde{f}(k) = \tilde{\psi}(k, 0) = B(k)\sinh k \implies B(k) = \frac{\tilde{f}(k)}{\sinh k}.
  \]
  Thus, for $k \neq 0$,
  \[
  \tilde{\psi}(k, y) = \tilde{f}(k)\frac{\sinh\bigl(k(1 - y)\bigr)}{\sinh k}.
  \]
- **Case $k = 0$:** The ODE reduces to $\dfrac{\partial^2 \tilde{\psi}}{\partial y^2}(0, y) = 0$, so $\tilde{\psi}(0, y) = c_0 + c_1 y$. Since $\tilde{\psi}(0, 0) = \tilde{f}(0) = e^0 - e^0 = 0$ and $\tilde{\psi}(0, 1) = 0$, we obtain $c_0 = c_1 = 0$, i.e. $\tilde{\psi}(0, y) = 0$.

Substituting the explicit formula for $\tilde{f}(k)$ produces an exact cancellation of $\sinh k$ on $0 < |k| \le 1$:
\[
\tilde{\psi}(k, y) = \begin{cases}
2\sinh\bigl(k(1 - y)\bigr) = e^{k(1 - y)} - e^{-k(1 - y)}, & |k| \le 1, \\
0, & |k| > 1,
\end{cases}
\]
which also agrees at $k = 0$ where $2\sinh 0 = 0$.

Now invert the Fourier transform to recover $\psi(x, y)$. Set $\alpha = 1 - y \in [0, 1]$:
\[
\psi(x, y) = \frac{1}{2\pi} \int_{-1}^{1} \left(e^{\alpha k} - e^{-\alpha k}\right) e^{ikx}\,dk.
\]
For $(\alpha, x) \neq (0, 0)$ (which holds for all $x \in \mathbb{R}$ whenever $0 \le y < 1$, and for all $x \neq 0$ when $y = 1$), we evaluate the two exponential integrals exactly as in Part 1:
\[
\int_{-1}^{1} e^{(\alpha + ix)k}\,dk = \frac{e^{\alpha + ix} - e^{-\alpha - ix}}{\alpha + ix}, \qquad
\int_{-1}^{1} e^{(-\alpha + ix)k}\,dk = \frac{e^{\alpha - ix} - e^{-\alpha + ix}}{\alpha - ix}.
\]
Setting $J_\alpha(x) = \dfrac{e^{\alpha + ix} - e^{-\alpha - ix}}{\alpha + ix}$, we have
\[
e^{\alpha + ix} - e^{-\alpha - ix} = 2\sinh\alpha \cos x + 2i\cosh\alpha \sin x,
\]
and multiplying numerator and denominator by $\alpha - ix$ gives
\[
J_\alpha(x) = \frac{2(\alpha\sinh\alpha\cos x + x\cosh\alpha\sin x) + 2i(\alpha\cosh\alpha\sin x - x\cos x\sinh\alpha)}{\alpha^2 + x^2}.
\]
Therefore,
\[
\psi(x, y) = \frac{1}{2\pi}\bigl(J_\alpha(x) - \overline{J_\alpha(x)}\bigr)
= \frac{i}{\pi}\operatorname{Im} J_\alpha(x)
= \frac{2i}{\pi(\alpha^2 + x^2)}\bigl(\alpha\cosh\alpha\sin x - x\cos x\sinh\alpha\bigr).
\]
Restoring $\alpha = 1 - y$:
\[
\psi(x, y) = \frac{2i}{\pi\bigl((1 - y)^2 + x^2\bigr)}\Bigl((1 - y)\cosh(1 - y)\sin x - x\cos x\sinh(1 - y)\Bigr).
\]

### Verification of regularity, PDE, boundary conditions, and decay

1. **Smoothness and removable singularity at $(x, y) = (0, 1)$:**
   Introducing $z = (1 - y) + ix$, the integral representation evaluates identically to
   \[
   \psi(x, y) = \frac{1}{2\pi}\int_{-1}^{1}\bigl(e^{kz} - e^{-k\bar{z}}\bigr)\,dk
   = \frac{1}{\pi}\left(\frac{\sinh z}{z} - \frac{\sinh\bar{z}}{\bar{z}}\right)
   = \frac{2i}{\pi}\operatorname{Im}\!\left(\frac{\sinh z}{z}\right).
   \]
   Since $\dfrac{\sinh z}{z} = \sum_{n=0}^{\infty} \dfrac{z^{2n}}{(2n+1)!}$ is entire, $\psi(x, y)$ is smooth (real-analytic) everywhere on $\mathbb{R}^2$. In particular, at $(x, y) = (0, 1)$ ($z = 0$), the closed-form rational expression has a removable $0/0$ singularity with
   \[
   \psi(x, y) = \frac{2i}{3\pi}(1 - y)x + O\Bigl(\bigl((1 - y)^2 + x^2\bigr)^2\Bigr) \to 0 = \psi(0, 1).
   \]
2. **Laplace's equation:** Since $\dfrac{\sinh z}{z}$ is holomorphic in $z = (1 - y) + ix$, its imaginary part is harmonic everywhere, so $\nabla^2 \psi = 0$ on $-\infty < x < \infty$, $0 < y < 1$ (equivalently, differentiating under the finite integral $\int_{-1}^1$ immediately gives $\psi_{xx} + \psi_{yy} = 0$ because $(-k^2 + k^2) = 0$).
3. **Boundary conditions:**
   - At $y = 0$: $\psi(x, 0) = \dfrac{2i}{\pi(1 + x^2)}(\cosh 1 \sin x - x\cos x \sinh 1) = f(x)$.
   - At $y = 1$: $z = ix$ makes $\dfrac{\sinh(ix)}{ix} = \dfrac{\sin x}{x} \in \mathbb{R}$, so $\psi(x, 1) = 0$ for all $x \in \mathbb{R}$.
4. **Decay as $|x| \to \infty$:** Uniformly for $y \in [0, 1]$, $|\psi(x, y)| = O(1/|x|) \to 0$ and $|\psi_x(x, y)| = O(1/|x|) \to 0$ as $|x| \to \infty$, confirming consistency with the Fourier-transform derivation and uniqueness among solutions vanishing as $|x| \to \infty$.

## Final Answers

1. Inverse Fourier transform of $\tilde{f}(k) = (e^k - e^{-k})\mathbf{1}_{|k|\le 1}$:
   \[
   f(x) = \frac{2i}{\pi(1 + x^2)}(\cosh 1 \sin x - x\cos x \sinh 1).
   \]
2. Solution of $\nabla^2\psi = 0$ in the strip $-\infty < x < \infty$, $0 \le y \le 1$ with $\psi(x, 0) = f(x)$, $\psi(x, 1) = 0$, and $\psi(x, y) \to 0$ as $|x| \to \infty$:
   \[
   \psi(x, y) = \frac{2i}{\pi\bigl((1 - y)^2 + x^2\bigr)}\Bigl((1 - y)\cosh(1 - y)\sin x - x\cos x\sinh(1 - y)\Bigr)
   \]
   (with the removable singularity at $(x, y) = (0, 1)$ taking the value $\psi(0, 1) = 0$; equivalently $\psi(x, y) = \dfrac{2i}{\pi}\operatorname{Im}\!\left(\dfrac{\sinh z}{z}\right)$ where $z = (1 - y) + ix$).

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_qx27_check.py` (executed with `python3`, exit code 0).

```text
IB_methods_sheet-1_qx27_check.py: ALL CHECKS PASSED
  Part 1 max inverse FT quadrature error:      2.669e-13
  Part 2 max closed-vs-integral error:         9.580e-14
  Part 2 max closed-vs-holomorphic error:      8.327e-17
  Boundary condition errors (y=0, y=1):        0.000e+00, 0.000e+00
  Max 4th-order FD Laplace residual |nabla^2|: 2.557e-11
```

## Confidence & Edge-Case Notes

- **Confidence:** High.
- **Phase B watchpoints:**
  1. **Factor of $(1 - y)$ in numerator:** When generalizing Part 1 ($\alpha = 1$) to Part 2 ($\alpha = 1 - y$), watch carefully that the first term inside the parentheses acquires a factor of $(1 - y)$: $(1 - y)\cosh(1 - y)\sin x - x\cos x\sinh(1 - y)$, and the denominator becomes $(1 - y)^2 + x^2$. Forgetting the $(1 - y)$ in front of $\cosh(1 - y)\sin x$ is an easy scaling slip!
  2. **Special mode $k = 0$:** The ODE $\tilde{\psi}_{yy} - k^2\tilde{\psi} = 0$ should explicitly handle $k = 0$ (where the general solution is $c_0 + c_1 y$, not $\sinh(ky)$) rather than dividing by $\sinh 0 = 0$.
  3. **Removable singularity at $(x, y) = (0, 1)$ and decay as $|x| \to \infty$:** Check that $\psi(0, 1) = 0$ is noted (since the rational formula evaluates to $0/0$ there), and that the decay assumption $\psi, \psi_x \to 0$ as $|x| \to \infty$ is stated to justify transforming $\psi_{xx}$ and ensuring uniqueness.
