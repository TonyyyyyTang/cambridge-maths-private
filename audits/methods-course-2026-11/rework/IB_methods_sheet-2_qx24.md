---
id: IB/methods/sheet-2/qx24
slug: IB_methods_sheet-2_qx24
source_archive: private/raw/sheets/IB/methods/2007-08/sheet-2.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-2_qx24_p2.png` (Question 8 on page 2 of Michaelmas 2007 Methods Examples II, Dr C. P. Caulfield; cleanly followed on the same page by Questions 9 and 10).

**Verdict: `ok`.** Every word, sign, coefficient, constant, and boundary condition matches verbatim, including the exact ordering `x(x-1)y'' + [(1+a+b)x - c]y' + aby = 0` in part (ii) and `y'' + 4y' + (4+\lambda)y = 0, y(0) = y(1) = 0.` in part (iii). Purely typographic differences only: the printed title *Recasting in Sturm-Liouville form.* is italicised and underlined with the full stop inside the underline, and `(iii)` starts on a new paragraph in the LaTeX fragment rather than inline after "constants.".

## Full Solution

A second-order linear differential equation
\[
A(x)y'' + B(x)y' + \big(C(x) + \lambda D(x)\big)y = 0
\]
is recast into Sturm–Liouville form
\[
\frac{d}{dx}\!\left(p(x)\frac{dy}{dx}\right) + \big(q(x) + \lambda w(x)\big)y = 0
\]
on an interval where $A(x)$ does not vanish by multiplying by an integrating factor $\mu(x)$ such that $(\mu A)' = \mu B$, i.e.
\[
p(x) := \mu(x) A(x) = \exp\!\left(\int \frac{B(x)}{A(x)}\,dx\right) \quad (\text{up to a non-zero constant}), \qquad \mu(x) = \frac{p(x)}{A(x)}.
\]

### Part (i)

The equation is Legendre's equation:
\[
(1 - x^2)y'' - 2xy' + n(n+1)y = 0.
\]
Here $A(x) = 1 - x^2$ and $B(x) = -2x = A'(x)$, so $\mu(x) = 1$ and the equation is already in Sturm–Liouville form:
\[
\frac{d}{dx}\!\left[(1 - x^2)\frac{dy}{dx}\right] + n(n+1)y = 0.
\]
Comparing with $\frac{d}{dx}(p y') + (q + \lambda w)y = 0$, we identify
\[
p(x) = 1 - x^2, \qquad q(x) = 0, \qquad \lambda = n(n+1), \qquad w(x) = 1
\]
(on $(-1,1)$, where $p(x) > 0$ and $p(\pm 1) = 0$, this is a singular Sturm–Liouville problem).

### Part (ii)

The equation is the hypergeometric equation:
\[
x(x - 1)y'' + \big[(1 + a + b)x - c\big]y' + aby = 0.
\]
Here $A(x) = x(x-1)$ and $B(x) = (1+a+b)x - c$. Decomposing $B(x)/A(x)$ into partial fractions:
\[
\frac{(1+a+b)x - c}{x(x-1)} = \frac{c}{x} + \frac{1 + a + b - c}{x - 1}.
\]
(Check: $c(x-1) + (1+a+b-c)x = (1+a+b)x - c$.) Integrating gives
\[
\int \frac{B(x)}{A(x)}\,dx = c\ln|x| + (1 + a + b - c)\ln|x - 1| + \text{const},
\]
so a positive function $p(x)$ on any of the three intervals $(-\infty,0)$, $(0,1)$, or $(1,\infty)$ is
\[
p(x) = |x|^c |x - 1|^{1 + a + b - c}.
\]
Because $A(x) = x(x-1)$ changes sign between $(0,1)$ and the exterior intervals, care must be taken with the sign of the integrating factor $\mu(x) = p(x)/[x(x-1)]$:

- **On $(0,1)$** (where $x > 0$ and $1 - x > 0$, so $x(x-1) = -x(1-x) < 0$):
  Take $p(x) = x^c (1 - x)^{1 + a + b - c} > 0$, giving $\mu(x) = -x^{c-1}(1-x)^{a+b-c}$. Since
  \[
  \frac{d}{dx}\!\left[x^c(1-x)^{1+a+b-c}\right] = x^{c-1}(1-x)^{a+b-c}\big[c(1-x) - (1+a+b-c)x\big] = -x^{c-1}(1-x)^{a+b-c}\big[(1+a+b)x - c\big],
  \]
  multiplying the ODE by $\mu(x)$ yields
  \[
  \frac{d}{dx}\!\left[x^c(1 - x)^{1 + a + b - c}\frac{dy}{dx}\right] - ab\, x^{c-1}(1 - x)^{a + b - c}\, y = 0.
  \]
  Holding $a+b$ and $c$ fixed so that $p$ and $w$ are independent of the eigenvalue parameter, this has Sturm–Liouville form with $p(x) = x^c(1-x)^{1+a+b-c} > 0$, $q(x) = 0$, $\lambda = -ab$, and weight $w(x) = x^{c-1}(1-x)^{a+b-c} > 0$.

- **On $(1,\infty)$** (or formally without absolute values, where $x > 0$ and $x - 1 > 0$):
  Take $p(x) = x^c(x - 1)^{1 + a + b - c} > 0$, giving $\mu(x) = +x^{c-1}(x-1)^{a+b-c}$ and
  \[
  \frac{d}{dx}\!\left[x^c(x - 1)^{1 + a + b - c}\frac{dy}{dx}\right] + ab\, x^{c-1}(x - 1)^{a + b - c}\, y = 0.
  \]
  (Similarly, on $(-\infty, 0)$, where $-x > 0$ and $1-x > 0$ so $x(x-1) = (-x)(1-x) > 0$, we get $\frac{d}{dx}\!\left[(-x)^c(1-x)^{1+a+b-c}y'\right] + ab\,(-x)^{c-1}(1-x)^{a+b-c}y = 0$. All three intervals are summarised simultaneously by $\frac{d}{dx}\!\left(|x|^c |x-1|^{1+a+b-c} y'\right) + \operatorname{sgn}\!\big(x(x-1)\big)\,ab\,|x|^{c-1}|x-1|^{a+b-c} y = 0$.)

### Part (iii)

Consider the boundary value problem
\[
y'' + 4y' + (4 + \lambda)y = 0, \qquad y(0) = y(1) = 0.
\]

#### 1. Sturm–Liouville form and reality of eigenvalues
Here $A(x) = 1$ and $B(x) = 4$, so the integrating factor is
\[
\mu(x) = \exp\!\left(\int 4\,dx\right) = e^{4x}.
\]
Multiplying by $e^{4x}$ yields the regular Sturm–Liouville problem on $[0,1]$:
\[
\frac{d}{dx}\!\left(e^{4x}\frac{dy}{dx}\right) + (4 + \lambda)e^{4x} y = 0, \qquad y(0) = y(1) = 0,
\]
with $p(x) = e^{4x} > 0$, $q(x) = 4e^{4x}$, and weight function $w(x) = e^{4x} > 0$.

For any two eigenpairs $(\lambda_m, y_m)$ and $(\lambda_n, y_n)$ (allowing complex $\lambda$ a priori), integrating by parts twice and using $y_m(0)=y_m(1)=y_n(0)=y_n(1)=0$:
\[
\begin{aligned}
(\lambda_m^* - \lambda_n)\int_0^1 e^{4x} y_m^*(x) y_n(x)\,dx
&= \int_0^1 \left[ y_n \frac{d}{dx}\!\left(e^{4x}(y_m^*)'\right) - y_m^*\frac{d}{dx}\!\left(e^{4x} y_n'\right) \right]\,dx \\
&= \Big[e^{4x}\big(y_n (y_m^*)' - y_m^* y_n'\big)\Big]_0^1 = 0.
\end{aligned}
\]
Taking $m = n$, since $e^{4x} > 0$ and $y_n \not\equiv 0$, we have $\int_0^1 e^{4x}|y_n(x)|^2\,dx > 0$, hence $\lambda_n^* = \lambda_n$; all eigenvalues are real.

#### 2. Eigenvalues and eigenfunctions
Set $y(x) = e^{-2x} u(x)$. Then $y' = e^{-2x}(u' - 2u)$ and $y'' = e^{-2x}(u'' - 4u' + 4u)$, so
\[
y'' + 4y' + (4 + \lambda)y = e^{-2x}\big[(u'' - 4u' + 4u) + 4(u' - 2u) + (4 + \lambda)u\big] = e^{-2x}(u'' + \lambda u) = 0.
\]
Since $e^{-2x} \neq 0$, the problem is equivalent to
\[
u'' + \lambda u = 0, \qquad u(0) = u(1) = 0.
\]
We check all real values of $\lambda$:

- **Case $\lambda < 0$:** Write $\lambda = -\alpha^2$ with $\alpha > 0$. The general solution is $u(x) = A\cosh(\alpha x) + B\sinh(\alpha x)$. Then $u(0) = A = 0$, and $u(1) = B\sinh(\alpha) = 0 \implies B = 0$ (as $\sinh\alpha > 0$ for $\alpha > 0$). Thus only the trivial solution $u \equiv 0$ exists; there are no negative eigenvalues.

- **Case $\lambda = 0$:** The general solution is $u(x) = A + Bx$. Then $u(0) = A = 0$ and $u(1) = B = 0$, giving only $u \equiv 0$. Thus $\lambda = 0$ is not an eigenvalue.

- **Case $\lambda > 0$:** Write $\lambda = \alpha^2$ with $\alpha > 0$. The general solution is $u(x) = A\cos(\alpha x) + B\sin(\alpha x)$. Then $u(0) = A = 0$, and $u(1) = B\sin(\alpha) = 0$. For a non-trivial solution ($B \neq 0$), we require $\sin(\alpha) = 0$, so $\alpha = n\pi$ for $n = 1, 2, 3, \dots$ (since $\alpha > 0$).

Transforming back via $y(x) = e^{-2x} u(x)$, the eigenvalues and corresponding eigenfunctions (up to arbitrary non-zero constant multiples) are
\[
\lambda_n = n^2 \pi^2, \qquad y_n(x) = e^{-2x}\sin(n\pi x), \qquad n = 1, 2, 3, \dots
\]

#### 3. Orthogonality relation
From the Sturm–Liouville identity derived above, for distinct eigenvalues $\lambda_m \neq \lambda_n$ (i.e. $m \neq n$), the eigenfunctions are orthogonal on $[0,1]$ with respect to the weight function $w(x) = e^{4x}$:
\[
\int_0^1 e^{4x} y_m(x) y_n(x)\,dx = 0 \qquad (m \neq n).
\]
In fact, evaluating directly for all $m, n \in \mathbb{N}$ with $y_n(x) = e^{-2x}\sin(n\pi x)$:
\[
\int_0^1 e^{4x} y_m(x) y_n(x)\,dx = \int_0^1 \sin(m\pi x)\sin(n\pi x)\,dx = \frac{1}{2}\delta_{mn}.
\]
(Thus the orthonormal eigenfunctions with respect to $w(x) = e^{4x}$ are $\hat{y}_n(x) = \sqrt{2}\,e^{-2x}\sin(n\pi x)$.)

## Final Answers

- **(i)** $\dfrac{d}{dx}\!\left[(1-x^2)\dfrac{dy}{dx}\right] + n(n+1)y = 0$ (with $p(x) = 1-x^2$, $q(x) = 0$, $\lambda = n(n+1)$, $w(x) = 1$).
- **(ii)** On $(0,1)$: $\dfrac{d}{dx}\!\left[x^c(1-x)^{1+a+b-c}\dfrac{dy}{dx}\right] - ab\, x^{c-1}(1-x)^{a+b-c} y = 0$ (with $p(x) = x^c(1-x)^{1+a+b-c}$, $w(x) = x^{c-1}(1-x)^{a+b-c}$, $\lambda = -ab$ for fixed $a+b, c$). On $(1,\infty)$ (or formally): $\dfrac{d}{dx}\!\left[x^c(x-1)^{1+a+b-c}\dfrac{dy}{dx}\right] + ab\, x^{c-1}(x-1)^{a+b-c} y = 0$ (and on $(-\infty,0)$: $\dfrac{d}{dx}\!\left[(-x)^c(1-x)^{1+a+b-c}\dfrac{dy}{dx}\right] + ab\,(-x)^{c-1}(1-x)^{a+b-c} y = 0$).
- **(iii)** Sturm–Liouville form: $\dfrac{d}{dx}\!\left(e^{4x}\dfrac{dy}{dx}\right) + (4+\lambda)e^{4x} y = 0$.
  - **Eigenvalues:** $\lambda_n = n^2\pi^2$ for $n = 1, 2, 3, \dots$
  - **Eigenfunctions:** $y_n(x) = e^{-2x}\sin(n\pi x)$ (up to a non-zero constant multiple), $n = 1, 2, 3, \dots$
  - **Orthogonality relation:** $\displaystyle \int_0^1 e^{4x} y_m(x) y_n(x)\,dx = 0$ for $m \neq n$ (explicitly $\displaystyle \int_0^1 e^{4x} y_m(x) y_n(x)\,dx = \frac{1}{2}\delta_{mn}$ for $y_n(x) = e^{-2x}\sin(n\pi x)$).

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-2_qx24_check.py`

Verbatim stdout:
```text
IB_methods_sheet-2_qx24_check.py: ALL CHECKS PASSED
  part (i) max SL recasting error:         2.012e-11
  part (ii) max SL recasting error:        1.662e-09
  part (iii) max shooting eigenvalue err:  6.639e-09
  part (iii) max orthogonality error:      2.220e-16
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  - **Sign trap in part (ii)!** Because the question writes $x(x-1)$ rather than $x(1-x)$, recasting with $(1-x)$ on $(0,1)$ flips the sign on the $aby$ term to **$-ab\,x^{c-1}(1-x)^{a+b-c}y$**, whereas recasting formally with $(x-1)$ (valid as a real positive $p(x)$ on $(1,\infty)$) gives **$+ab\,x^{c-1}(x-1)^{a+b-c}y$**. Any solution mixing $(1-x)$ in $p(x)$ with $+ab$ in the zero-order term has a sign error!
  - In part (iii), verify that $\lambda < 0$ and $\lambda = 0$ (and reality of $\lambda$) are explicitly ruled out, not merely assumed positive, and that the weight function $w(x) = e^{4x}$ (NOT $e^{-4x}$ or $e^{2x}$!) appears in the orthogonality relation.
