---
id: IB/methods/sheet-1/qx28
slug: IB_methods_sheet-1_qx28
source_archive: private/raw/sheets/IB/methods/2011-12/sheet-1.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against the 200 DPI render `IB_methods_sheet-1_qx28_p2.png` (upper-middle of page 2, Question 7, followed immediately by Question 8).

- **Result:** `ok` — exact mathematical and verbal match across all three subparts (i), (ii), (iii), including every coefficient, sign (`x(x-1)`, `+[(1+a+b)x-c]`, `+ab`), constant (`n, a, b, and c`), boundary condition (`y(0) = y(1) = 0`), and wording.
- **Purely typographic differences:** In the printed PDF, the heading *Recasting in Sturm-Liouville form.* is underlined as well as italicised (with the full stop inside the underline), whereas `question.tex` writes `\emph{Recasting in Sturm-Liouville form}.`. No mathematical or verbal content is affected.

## Full Solution

Recall that a second-order linear differential equation in Sturm–Liouville (self-adjoint) form has the structure
\[
\frac{\mathrm{d}}{\mathrm{d}x}\!\left(p(x)\frac{\mathrm{d}y}{\mathrm{d}x}\right) + \big(\tilde{q}(x) + \lambda w(x)\big)y = 0
\qquad\Longleftrightarrow\qquad
-\frac{\mathrm{d}}{\mathrm{d}x}\!\left(p(x)\frac{\mathrm{d}y}{\mathrm{d}x}\right) + q(x)\,y = \lambda\, w(x)\,y
\]
(where $q = -\tilde{q}$), with $p(x) > 0$ and $w(x) > 0$ on the open interval of interest. Given $A(x)y'' + B(x)y' + C(x)y = 0$ with $A(x) \neq 0$, multiplying by the integrating factor $\mu(x) = p(x)/A(x)$, where
\[
p(x) = \exp\!\left(\int \frac{B(x)}{A(x)}\,\mathrm{d}x\right) > 0,
\]
ensures $p'(x) = \frac{B(x)}{A(x)}p(x) = \mu(x)B(x)$, so that $\mu(x)(A y'' + B y') = \frac{\mathrm{d}}{\mathrm{d}x}(p y')$.

---

### Part (i)

The equation is Legendre's equation:
\[
(1-x^2)y'' - 2xy' + n(n+1)y = 0.
\]
Since $\frac{\mathrm{d}}{\mathrm{d}x}(1-x^2) = -2x$, the leading two terms are already an exact derivative (integrating factor $\mu(x) = 1$). Thus the Sturm–Liouville form is
\[
\frac{\mathrm{d}}{\mathrm{d}x}\!\left((1-x^2)\frac{\mathrm{d}y}{\mathrm{d}x}\right) + n(n+1)y = 0
\qquad\Longleftrightarrow\qquad
-\frac{\mathrm{d}}{\mathrm{d}x}\!\left((1-x^2)\frac{\mathrm{d}y}{\mathrm{d}x}\right) = n(n+1)y,
\]
with $p(x) = 1 - x^2$ (positive on the standard domain $(-1, 1)$), $q(x) = 0$, weight function $w(x) = 1$, and eigenvalue parameter $\lambda = n(n+1)$.

---

### Part (ii)

The equation is Gauss's hypergeometric differential equation (written with overall opposite sign):
\[
x(x-1)y'' + \big[(1 + a + b)x - c\big]y' + ab\, y = 0.
\]
Here $A(x) = x(x-1)$ and $B(x) = (1 + a + b)x - c$. Computing $B(x)/A(x)$ by partial fractions:
\[
\frac{(1 + a + b)x - c}{x(x-1)} = \frac{K}{x} + \frac{M}{x-1} = \frac{(K + M)x - K}{x(x-1)}.
\]
Comparing coefficients gives $-K = -c \implies K = c$ and $K + M = 1 + a + b \implies M = 1 + a + b - c$. Hence
\[
\frac{p'(x)}{p(x)} = \frac{c}{x} + \frac{1 + a + b - c}{x - 1}.
\]
Integrating on any open interval avoiding the singular points $x = 0$ and $x = 1$ gives
\[
p(x) = |x|^c |x - 1|^{1 + a + b - c} > 0,
\]
and therefore the integrating factor is
\[
\mu(x) = \frac{p(x)}{x(x-1)} = \operatorname{sgn}\!\big(x(x-1)\big)\, |x|^{c-1}|x - 1|^{a+b-c}.
\]
Because the problem does not specify an interval, we write out the explicit real-valued positive-$p$, positive-$w$ Sturm–Liouville form on each connected component of $\mathbb{R}\setminus\{0, 1\}$, paying careful attention to the sign of $x(x-1)$:

1. **On the standard bounded interval $0 < x < 1$** (where $x > 0$, $1 - x > 0$, and $x(x-1) = -x(1-x) < 0$):
   We have $p(x) = x^c(1-x)^{1+a+b-c} > 0$, whose derivative is
   \[
   p'(x) = x^{c-1}(1-x)^{a+b-c}\big[c(1-x) - (1+a+b-c)x\big] = x^{c-1}(1-x)^{a+b-c}\big[c - (1+a+b)x\big].
   \]
   Multiplying the given ODE by $\mu(x) = -x^{c-1}(1-x)^{a+b-c}$ (equivalently multiplying $x(1-x)y'' + [c - (1+a+b)x]y' - aby = 0$ by the positive factor $x^{c-1}(1-x)^{a+b-c}$) yields
   \[
   \frac{\mathrm{d}}{\mathrm{d}x}\!\left(x^c(1-x)^{1+a+b-c}\frac{\mathrm{d}y}{\mathrm{d}x}\right) - ab\, x^{c-1}(1-x)^{a+b-c}\, y = 0,
   \]
   with $p(x) = x^c(1-x)^{1+a+b-c} > 0$ and $w(x) = x^{c-1}(1-x)^{a+b-c} > 0$.
   *(Note: When $a + b$ and $c$ are held fixed while the product $ab$ varies — e.g. for shifted Jacobi polynomials where $a = -n$ and $b = n + \gamma$ with $\gamma = a+b$ fixed — both $p(x)$ and $w(x)$ are independent of $n$, and $\lambda = -ab = n(n+\gamma)$ is the Sturm–Liouville eigenvalue with $q(x) = 0$.)*

2. **On $x > 1$** (where $x > 0$, $x - 1 > 0$, and $x(x-1) > 0$):
   Here $p(x) = x^c(x-1)^{1+a+b-c} > 0$ and $\mu(x) = +x^{c-1}(x-1)^{a+b-c} > 0$, giving
   \[
   \frac{\mathrm{d}}{\mathrm{d}x}\!\left(x^c(x-1)^{1+a+b-c}\frac{\mathrm{d}y}{\mathrm{d}x}\right) + ab\, x^{c-1}(x-1)^{a+b-c}\, y = 0,
   \]
   with $p(x) = x^c(x-1)^{1+a+b-c} > 0$, $w(x) = x^{c-1}(x-1)^{a+b-c} > 0$, and $\lambda = +ab$ (when $a + b, c$ are fixed).

3. **On $x < 0$** (where $-x > 0$, $1 - x > 0$, and $x(x-1) = (-x)(1-x) > 0$):
   Here $p(x) = (-x)^c(1-x)^{1+a+b-c} > 0$ and $\mu(x) = +(-x)^{c-1}(1-x)^{a+b-c} > 0$, giving
   \[
   \frac{\mathrm{d}}{\mathrm{d}x}\!\left((-x)^c(1-x)^{1+a+b-c}\frac{\mathrm{d}y}{\mathrm{d}x}\right) + ab\, (-x)^{c-1}(1-x)^{a+b-c}\, y = 0.
   \]

---

### Part (iii)

Consider the boundary value problem
\[
y'' + 4y' + (4 + \lambda)y = 0 \quad (0 < x < 1), \qquad y(0) = y(1) = 0.
\]

#### 1. Sturm–Liouville form

The integrating factor is $\mu(x) = \exp\!\left(\int 4\,\mathrm{d}x\right) = e^{4x}$. Multiplying the ODE by $e^{4x}$ gives
\[
\frac{\mathrm{d}}{\mathrm{d}x}\!\left(e^{4x}\frac{\mathrm{d}y}{\mathrm{d}x}\right) + 4e^{4x} y + \lambda e^{4x} y = 0
\qquad\Longleftrightarrow\qquad
-\frac{\mathrm{d}}{\mathrm{d}x}\!\left(e^{4x}\frac{\mathrm{d}y}{\mathrm{d}x}\right) - 4e^{4x} y = \lambda e^{4x} y.
\]
This is a regular Sturm–Liouville problem on $[0, 1]$ with $p(x) = e^{4x} > 0$, $q(x) = -4e^{4x}$ (or $\tilde{q}(x) = 4e^{4x}$), weight function $w(x) = e^{4x} > 0$, and homogeneous Dirichlet boundary conditions $y(0) = y(1) = 0$. In particular, all eigenvalues $\lambda$ are real.

#### 2. Eigenvalues and eigenfunctions

The characteristic equation is
\[
r^2 + 4r + (4 + \lambda) = 0 \iff (r + 2)^2 = -\lambda.
\]
(Equivalently, setting $y(x) = e^{-2x}u(x)$ transforms the problem into $u'' + \lambda u = 0$, $u(0) = u(1) = 0$.) Checking all cases:

- **$\lambda < 0$:** Write $\lambda = -\kappa^2$ ($\kappa > 0$). Then $r = -2 \pm \kappa$, so $y(x) = e^{-2x}\big(A\cosh(\kappa x) + B\sinh(\kappa x)\big)$. Now $y(0) = 0 \implies A = 0$, and $y(1) = 0 \implies B e^{-2}\sinh\kappa = 0 \implies B = 0$ (as $\sinh\kappa > 0$). No negative eigenvalues.
- **$\lambda = 0$:** Repeated root $r = -2$, so $y(x) = (A + Bx)e^{-2x}$. Then $y(0) = 0 \implies A = 0$, and $y(1) = 0 \implies B e^{-2} = 0 \implies B = 0$. Hence $\lambda = 0$ is not an eigenvalue.
- **$\lambda > 0$:** Write $\lambda = k^2$ ($k > 0$). Then $r = -2 \pm ik$, so $y(x) = e^{-2x}\big(A\cos(kx) + B\sin(kx)\big)$. Now $y(0) = 0 \implies A = 0$, and $y(1) = 0 \implies B e^{-2}\sin k = 0$. For a non-trivial solution ($B \neq 0$), we require $\sin k = 0$, so $k = n\pi$ for $n = 1, 2, 3, \dots$. (Allowing $k \in \mathbb{C}\setminus\{0\}$ still gives $\sin k = 0 \iff k \in \pi\mathbb{Z}\setminus\{0\}$, confirming independently that there are no non-real eigenvalues.)

Thus the eigenvalues and eigenfunctions are
\[
\lambda_n = n^2\pi^2, \qquad y_n(x) = e^{-2x}\sin(n\pi x) \quad (\text{up to a non-zero constant multiple}), \qquad n = 1, 2, 3, \dots
\]

#### 3. Orthogonality relation

By Sturm–Liouville theory (directly verified via integration by parts: $(\lambda_n - \lambda_m)\int_0^1 e^{4x} y_n y_m\,\mathrm{d}x = \big[e^{4x}(y_n y_m' - y_m y_n')\big]_0^1 = 0$), the eigenfunctions are orthogonal on $[0, 1]$ with respect to the weight function $w(x) = e^{4x}$:
\[
\int_0^1 e^{4x} y_n(x) y_m(x)\,\mathrm{d}x = 0 \qquad (n \neq m).
\]
Evaluating explicitly for $y_n(x) = e^{-2x}\sin(n\pi x)$ for all $n, m \ge 1$:
\[
\int_0^1 e^{4x} y_n(x) y_m(x)\,\mathrm{d}x
= \int_0^1 \sin(n\pi x)\sin(m\pi x)\,\mathrm{d}x
= \frac{1}{2}\delta_{nm}.
\]

## Final Answers

- **(i)** $\displaystyle\frac{\mathrm{d}}{\mathrm{d}x}\!\left((1-x^2)\frac{\mathrm{d}y}{\mathrm{d}x}\right) + n(n+1)y = 0$ (with $p(x) = 1 - x^2$, $w(x) = 1$).
- **(ii)**
  - On $(0, 1)$ (standard bounded domain, $p, w > 0$):
    \[
    \frac{\mathrm{d}}{\mathrm{d}x}\!\left(x^c(1-x)^{1+a+b-c}\frac{\mathrm{d}y}{\mathrm{d}x}\right) - ab\, x^{c-1}(1-x)^{a+b-c}\, y = 0.
    \]
  - On $(1, \infty)$ ($p, w > 0$):
    \[
    \frac{\mathrm{d}}{\mathrm{d}x}\!\left(x^c(x-1)^{1+a+b-c}\frac{\mathrm{d}y}{\mathrm{d}x}\right) + ab\, x^{c-1}(x-1)^{a+b-c}\, y = 0.
    \]
  - (And on $(-\infty, 0)$, replace $x$ by $-x$ and $x - 1$ by $1 - x$ in the $(1, \infty)$ formula, retaining $+ab$.)
- **(iii)** Eigenvalues $\lambda_n = n^2\pi^2$, eigenfunctions $y_n(x) = e^{-2x}\sin(n\pi x)$ ($n = 1, 2, 3, \dots$), and orthogonality relation with weight $w(x) = e^{4x}$:
  \[
  \int_0^1 e^{4x} y_n(x) y_m(x)\,\mathrm{d}x = \frac{1}{2}\delta_{nm} \qquad \left(\text{in particular } \int_0^1 e^{4x} y_n(x) y_m(x)\,\mathrm{d}x = 0 \text{ for } n \neq m\right).
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_qx28_check.py` (exit code 0).

```text
Part (i) Legendre SL check passed (max_err = 1.70e-10).
Part (ii) SL checks passed on (0,1) (1.08e-10), (1,inf) (9.04e-08), (-inf,0) (1.91e-07); 2F1 polynomial orthogonality err=4.16e-17.
Part (iii) checks passed: max_ode_err=2.84e-14, max_orth_err=7.77e-16.
ALL CHECKS PASSED FOR IB/methods/sheet-1/qx28.
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  1. **Critical sign change between $(0, 1)$ and $(1, \infty)$ in part (ii):** Because the leading coefficient is printed as $x(x-1)$ (negative on $(0, 1)$ and positive on $(1, \infty)$):
     - Using $p(x) = x^c(1-x)^{1+a+b-c}$ on $(0, 1)$ requires a **minus** sign before $ab$, i.e. $\big(x^c(1-x)^{1+a+b-c}y'\big)' - ab\,x^{c-1}(1-x)^{a+b-c}y = 0$.
     - Using $p(x) = x^c(x-1)^{1+a+b-c}$ on $(1, \infty)$ keeps a **plus** sign before $ab$, i.e. $\big(x^c(x-1)^{1+a+b-c}y'\big)' + ab\,x^{c-1}(x-1)^{a+b-c}y = 0$.
     Watch carefully for any solution that writes $(1-x)^{1+a+b-c}$ with $+ab$ (wrong sign!) OR writes $(x-1)^{1+a+b-c}$ without noting domain $x > 1$ (since $(x-1)^{1+a+b-c}$ is generally non-real on $(0, 1)$ for real non-integer $a+b-c$).
  2. **Role of $a, b$ as parameters in part (ii):** Because $a+b$ appears in $p(x)$ and $w(x)$, treating $\mp ab$ as an eigenvalue parameter requires holding $a + b$ and $c$ fixed (as in Jacobi polynomials).
  3. **Part (iii) pitfalls:** Check that $\lambda < 0$ and $\lambda = 0$ are explicitly excluded, and that the orthogonality weight is $w(x) = e^{4x}$ (cancelling $e^{-2x} \cdot e^{-2x}$).
