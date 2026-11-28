---
id: IB/methods/sheet-1/q08
slug: IB_methods_sheet-1_q08
source_archive: private/raw/sheets/IB/methods/2025-26/sheet-1.pdf
pages_checked: [1, 2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against the 200 DPI renders `IB_methods_sheet-1_q08_p1.png` (bottom of page 1, Question 8) and `IB_methods_sheet-1_q08_p2.png` (top of page 2, continuation of Question 8, followed immediately by Question 9).

- **Result:** `ok` — exact mathematical and verbal match across all three subparts (i), (ii), (iii), including every coefficient, sign, domain ($[-1,1]$, $[0,1]$, $0 < x < 1$), condition ($a > b > 0$ constant, $\lambda$ constant), boundary condition ($y(0) = 0$, $y(1) = 0$), and punctuation mark.
- **Purely typographic differences:** In the printed PDF, equations (i) and (ii) are set side-by-side on a single display line whereas `question.tex` stacks them in an `aligned` environment; in (iii), the three lines inside the left brace are aligned on `=` in the PDF rather than left-aligned via `cases`. No mathematical or verbal content is affected.

## Full Solution

Recall the standard Sturm–Liouville eigenvalue form on a finite interval $[c, d]$:
\[
-\frac{\mathrm{d}}{\mathrm{d}x}\!\left(p(x)\frac{\mathrm{d}y}{\mathrm{d}x}\right) + q(x)\,y = \lambda\, w(x)\,y
\qquad\Longleftrightarrow\qquad
\frac{\mathrm{d}}{\mathrm{d}x}\!\left(p(x)\frac{\mathrm{d}y}{\mathrm{d}x}\right) + \big(\lambda w(x) - q(x)\big)y = 0,
\]
with $p(x) > 0$ and $w(x) > 0$ on $(c, d)$. Given a second-order linear ODE
\[
A(x)y'' + B(x)y' + C(x)y + \lambda D(x)y = 0
\]
with $A(x) > 0$ on $(c, d)$, multiplying through by the positive integrating factor
\[
\mu(x) = \frac{p(x)}{A(x)}, \qquad \text{where}\quad p(x) = \exp\!\left(\int \frac{B(x)}{A(x)}\,\mathrm{d}x\right) > 0,
\]
converts the leading two terms into $\frac{\mathrm{d}}{\mathrm{d}x}(p(x)y')$ since $p'(x) = \frac{B(x)}{A(x)}p(x) = \mu(x)B(x)$.

A Sturm–Liouville problem on a finite closed interval $[c, d]$ is **regular** if $p, p', q, w$ are continuous on $[c, d]$ and $p(x) > 0$, $w(x) > 0$ for all $x \in [c, d]$ (together with regular separated boundary conditions); it is **singular** if $p$ or $w$ vanishes or fails to be continuous/finite at either endpoint (or if the interval is unbounded).

---

### Part (i)

The equation on $[-1, 1]$ is Legendre's equation:
\[
(1-x^2)y'' - 2xy' + \lambda y = 0.
\]
Since $\frac{\mathrm{d}}{\mathrm{d}x}(1-x^2) = -2x$, the equation is already in self-adjoint Sturm–Liouville form (the integrating factor is $\mu(x) = 1$):
\[
-\frac{\mathrm{d}}{\mathrm{d}x}\!\left((1-x^2)\frac{\mathrm{d}y}{\mathrm{d}x}\right) = \lambda y
\qquad\Longleftrightarrow\qquad
\frac{\mathrm{d}}{\mathrm{d}x}\!\left((1-x^2)\frac{\mathrm{d}y}{\mathrm{d}x}\right) + \lambda y = 0,
\]
with
\[
p(x) = 1 - x^2, \qquad q(x) = 0, \qquad w(x) = 1.
\]

**Singularity on $[-1, 1]$:** Although $p(x) = 1 - x^2 > 0$ on $(-1, 1)$ and $w(x) = 1 > 0$ on $[-1, 1]$, we have $p(-1) = 0$ and $p(1) = 0$. Because $p(x)$ vanishes at both endpoints $x = \pm 1$, problem (i) is **singular** (at both endpoints).

---

### Part (ii)

The equation on $[0, 1]$ is
\[
x(1-x)y'' - (ax - b)y' + \lambda y = 0, \qquad a > b > 0.
\]
Here $A(x) = x(1-x) > 0$ on $(0, 1)$ and $B(x) = b - ax$. We seek $p(x) > 0$ on $(0, 1)$ satisfying
\[
\frac{p'(x)}{p(x)} = \frac{B(x)}{A(x)} = \frac{b - ax}{x(1-x)}.
\]
Decomposing into partial fractions:
\[
\frac{b - ax}{x(1-x)} = \frac{K}{x} + \frac{M}{1-x} = \frac{K + (M - K)x}{x(1-x)},
\]
giving $K = b$ and $M - K = -a$, so $M = -(a - b)$. Hence
\[
\frac{p'(x)}{p(x)} = \frac{b}{x} - \frac{a-b}{1-x}.
\]
Integrating on $(0, 1)$ (where $x > 0$ and $1 - x > 0$, noting $\int \frac{-1}{1-x}\,\mathrm{d}x = +\ln(1-x)$) gives
\[
\ln p(x) = b\ln x + (a - b)\ln(1 - x) + \text{const}.
\]
Choosing the constant of integration to be $0$, we obtain
\[
p(x) = x^b (1 - x)^{a - b},
\]
and the positive integrating factor on $(0, 1)$ is
\[
\mu(x) = \frac{p(x)}{x(1-x)} = x^{b-1}(1-x)^{a-b-1}.
\]
Check directly by differentiation:
\[
p'(x) = b x^{b-1}(1-x)^{a-b} - (a-b)x^b(1-x)^{a-b-1} = x^{b-1}(1-x)^{a-b-1}\big[b(1-x) - (a-b)x\big] = \mu(x)(b - ax).
\]
Multiplying the original equation by $\mu(x) = x^{b-1}(1-x)^{a-b-1}$ yields the Sturm–Liouville form on $[0, 1]$:
\[
-\frac{\mathrm{d}}{\mathrm{d}x}\!\left(x^b(1-x)^{a-b}\frac{\mathrm{d}y}{\mathrm{d}x}\right) = \lambda\, x^{b-1}(1-x)^{a-b-1} y
\qquad\Longleftrightarrow\qquad
\frac{\mathrm{d}}{\mathrm{d}x}\!\left(x^b(1-x)^{a-b}\frac{\mathrm{d}y}{\mathrm{d}x}\right) + \lambda\, x^{b-1}(1-x)^{a-b-1} y = 0,
\]
with
\[
p(x) = x^b(1-x)^{a-b}, \qquad q(x) = 0, \qquad w(x) = x^{b-1}(1-x)^{a-b-1}.
\]

**Singularity on $[0, 1]$:** Since $a > b > 0$, both exponents $b > 0$ and $a - b > 0$ are strictly positive. Consequently,
\[
p(0) = 0 \qquad\text{and}\qquad p(1) = 0.
\]
(Additionally, the weight function $w(x) = x^{b-1}(1-x)^{a-b-1}$ vanishes at $x = 0$ if $b > 1$ or diverges as $x \to 0^+$ if $0 < b < 1$, and vanishes at $x = 1$ if $a - b > 1$ or diverges as $x \to 1^-$ if $0 < a - b < 1$; note that $w$ is integrable on $(0, 1)$ with $\int_0^1 w(x)\,\mathrm{d}x = B(b, a - b) < \infty$ precisely because $b > 0$ and $a - b > 0$.) Because $p(x)$ vanishes at both endpoints $x = 0$ and $x = 1$, problem (ii) is also **singular** (at both endpoints).

Thus **both** problems (i) and (ii) are singular Sturm–Liouville problems.

---

### Part (iii)

Consider the boundary value problem
\[
y'' + 4y' + (4 + \lambda)y = 0 \quad (0 < x < 1), \qquad y(0) = 0, \quad y(1) = 0.
\]

#### 1. Sturm–Liouville form

The integrating factor is $\mu(x) = \exp\!\left(\int 4\,\mathrm{d}x\right) = e^{4x}$. Multiplying the ODE by $e^{4x}$ gives
\[
e^{4x}y'' + 4e^{4x}y' + (4 + \lambda)e^{4x}y = 0
\qquad\Longleftrightarrow\qquad
-\frac{\mathrm{d}}{\mathrm{d}x}\!\left(e^{4x}\frac{\mathrm{d}y}{\mathrm{d}x}\right) - 4e^{4x}y = \lambda e^{4x}y,
\]
(equivalently $\frac{\mathrm{d}}{\mathrm{d}x}(e^{4x}y') + (4 + \lambda)e^{4x}y = 0$). This is a **regular** Sturm–Liouville problem on $[0, 1]$ with $p(x) = e^{4x} > 0$, $q(x) = -4e^{4x}$, weight function $w(x) = e^{4x} > 0$, and homogeneous Dirichlet boundary conditions $y(0) = y(1) = 0$. In particular, all eigenvalues $\lambda$ are real.

#### 2. Eigenvalues and eigenfunctions

The ODE has constant coefficients, with auxiliary equation
\[
r^2 + 4r + (4 + \lambda) = 0 \iff (r + 2)^2 = -\lambda.
\]
(Equivalently, setting $y(x) = e^{-2x}u(x)$ reduces the BVP to $u'' + \lambda u = 0$ with $u(0) = u(1) = 0$.) We examine all real values of $\lambda$ (note also that the complex check below rules out non-real $\lambda$ independently of SL theory):

- **Case $\lambda < 0$:** Write $\lambda = -\kappa^2$ with $\kappa > 0$. The roots are $r = -2 \pm \kappa$, giving the general solution
  \[
  y(x) = e^{-2x}\big(A\cosh(\kappa x) + B\sinh(\kappa x)\big).
  \]
  The condition $y(0) = 0$ implies $A = 0$, and then $y(1) = 0$ gives $B e^{-2}\sinh\kappa = 0$. Since $\kappa > 0$, $\sinh\kappa > 0$, so $B = 0$ and $y \equiv 0$. Thus there are no negative eigenvalues.

- **Case $\lambda = 0$:** There is a repeated root $r = -2$, giving
  \[
  y(x) = (A + Bx)e^{-2x}.
  \]
  Then $y(0) = 0 \implies A = 0$, and $y(1) = 0 \implies B e^{-2} = 0 \implies B = 0$, so $y \equiv 0$. Thus $\lambda = 0$ is not an eigenvalue.

- **Case $\lambda > 0$:** Write $\lambda = k^2$ with $k > 0$. The roots are $r = -2 \pm ik$, giving the general solution
  \[
  y(x) = e^{-2x}\big(A\cos(kx) + B\sin(kx)\big).
  \]
  The condition $y(0) = 0$ gives $A = 0$. Then $y(1) = 0$ requires
  \[
  B e^{-2}\sin k = 0.
  \]
  For a non-trivial solution ($B \neq 0$) with $k > 0$, we must have $\sin k = 0$, i.e. $k = n\pi$ for $n = 1, 2, 3, \dots$.
  (Even allowing arbitrary $k \in \mathbb{C}\setminus\{0\}$ with $\lambda = k^2$, $\sin k = 0$ still forces $k \in \pi\mathbb{Z}\setminus\{0\}$, independently confirming no non-real eigenvalues exist.)

Therefore, the eigenvalues and eigenfunctions are
\[
\lambda_n = n^2\pi^2, \qquad y_n(x) = e^{-2x}\sin(n\pi x) \quad (\text{up to a non-zero scalar multiple}), \qquad n = 1, 2, 3, \dots
\]

#### 3. Orthogonality relation

From the Sturm–Liouville form with weight function $w(x) = e^{4x}$, let $y_n, y_m$ be eigenfunctions with eigenvalues $\lambda_n, \lambda_m$. Integrating by parts on $[0, 1]$ (Lagrange's identity):
\[
(\lambda_n - \lambda_m)\int_0^1 e^{4x} y_n(x) y_m(x)\,\mathrm{d}x
= \int_0^1 \Big[y_n(e^{4x}y_m')' - y_m(e^{4x}y_n')'\Big]\mathrm{d}x
= \Big[e^{4x}\big(y_n y_m' - y_m y_n'\big)\Big]_0^1 = 0,
\]
since $y_n(0) = y_n(1) = y_m(0) = y_m(1) = 0$. Hence for $n \neq m$ ($\lambda_n \neq \lambda_m$),
\[
\int_0^1 e^{4x} y_n(x) y_m(x)\,\mathrm{d}x = 0.
\]
More explicitly, evaluating the weighted inner product for all $n, m \ge 1$ with $y_n(x) = e^{-2x}\sin(n\pi x)$:
\[
\int_0^1 e^{4x} y_n(x) y_m(x)\,\mathrm{d}x
= \int_0^1 e^{4x} \big(e^{-2x}\sin(n\pi x)\big)\big(e^{-2x}\sin(m\pi x)\big)\,\mathrm{d}x
= \int_0^1 \sin(n\pi x)\sin(m\pi x)\,\mathrm{d}x
= \frac{1}{2}\delta_{nm}.
\]
(Normalised eigenfunctions with respect to $w(x) = e^{4x}$ are $\hat{y}_n(x) = \sqrt{2}\,e^{-2x}\sin(n\pi x)$, satisfying $\int_0^1 e^{4x}\hat{y}_n(x)\hat{y}_m(x)\,\mathrm{d}x = \delta_{nm}$.)

## Final Answers

- **(i) Sturm–Liouville form on $[-1, 1]$:**
  \[
  -\frac{\mathrm{d}}{\mathrm{d}x}\!\left((1-x^2)\frac{\mathrm{d}y}{\mathrm{d}x}\right) = \lambda y
  \qquad \big(p(x) = 1 - x^2,\; q(x) = 0,\; w(x) = 1\big).
  \]
- **(ii) Sturm–Liouville form on $[0, 1]$ ($a > b > 0$):**
  \[
  -\frac{\mathrm{d}}{\mathrm{d}x}\!\left(x^b(1-x)^{a-b}\frac{\mathrm{d}y}{\mathrm{d}x}\right) = \lambda\, x^{b-1}(1-x)^{a-b-1} y
  \qquad \big(p(x) = x^b(1-x)^{a-b},\; q(x) = 0,\; w(x) = x^{b-1}(1-x)^{a-b-1}\big).
  \]
- **Singularity:** **Both** problems (i) and (ii) are **singular** at both endpoints: in (i), $p(\pm 1) = 0$; in (ii), since $b > 0$ and $a - b > 0$, $p(0) = p(1) = 0$ (and $w$ also vanishes or blows up at an endpoint whenever $b \neq 1$ or $a - b \neq 1$).
- **(iii) Eigenvalues, eigenfunctions, and orthogonality relation:**
  \[
  \lambda_n = n^2\pi^2, \qquad y_n(x) = e^{-2x}\sin(n\pi x) \quad (n = 1, 2, 3, \dots),
  \]
  with weight function $w(x) = e^{4x}$ and orthogonality relation
  \[
  \int_0^1 e^{4x} y_n(x) y_m(x)\,\mathrm{d}x = \frac{1}{2}\delta_{nm} \qquad \left(\text{in particular } \int_0^1 e^{4x} y_n(x) y_m(x)\,\mathrm{d}x = 0 \text{ for } n \neq m\right).
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_q08_check.py` (exit code 0).

```text
Part (i) SL form check passed (max_err = 1.87e-10), p(-1)=p(1)=0.
Part (ii) SL form & Beta(b,a-b) weight checks passed for 5 (a,b) pairs (max_err = 9.92e-11).
Part (iii) checks passed: max_ode_err=2.84e-14, max_bc_err=9.94e-17, max_orth_err=7.77e-16, RK4 shooting verified.
ALL CHECKS PASSED FOR IB/methods/sheet-1/q08.
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  1. **Sign error when integrating $-\frac{a-b}{1-x}$ in (ii):** A very common error is writing $(1-x)^{-(a-b)}$ instead of $(1-x)^{+(a-b)}$, because $\int \frac{\mathrm{d}x}{1-x} = -\ln(1-x)$. Verify $p(x) = x^b(1-x)^{a-b}$ and $w(x) = x^{b-1}(1-x)^{a-b-1}$.
  2. **Singularity classification in (ii):** Both endpoints $0$ and $1$ are singular precisely because $a > b > 0 \implies b > 0$ and $a - b > 0$, forcing $p(0) = p(1) = 0$ (even in the special cases $b = 1$ or $a - b = 1$ where $w$ stays finite and non-zero at an endpoint!). Watch out for answers that claim only one problem is singular or forget to justify using $a > b > 0$.
  3. **Checking $\lambda \le 0$ in (iii):** Any complete solution must explicitly rule out $\lambda < 0$ and $\lambda = 0$ (not merely assume $\lambda > 0$).
  4. **Weight function in (iii):** The weight function in the orthogonality relation is $w(x) = e^{4x}$ (NOT $e^{-4x}$ or $e^{2x}$); check that $e^{4x}(e^{-2x})^2 = 1$ cancels cleanly.
