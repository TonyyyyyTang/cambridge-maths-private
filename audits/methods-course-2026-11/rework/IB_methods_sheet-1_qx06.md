---
id: IB/methods/sheet-1/qx06
slug: IB_methods_sheet-1_qx06
source_archive: private/raw/sheets/IB/methods/2015-16/all-sheets.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against the 200 DPI render `IB_methods_sheet-1_qx06_p2.png` (middle of page 2, Question 8, followed immediately by Question 9).

- **Result:** `ok` — exact mathematical and verbal match across all three subparts (i), (ii), (iii), including every equation, boundary condition, recall note, hint (`\beta = j_n + \epsilon` as `\epsilon \to 0`), sign, factor, and condition (`\beta \ne \alpha`, `n \ne m`, `\tilde{\lambda}`).
- **Purely typographic differences:** In the printed PDF, the heading *Bessel's equation.* is underlined as well as italicised; in (ii), the second and third integrals together with the hint appear on a single display line whereas `question.tex` splits them across separate display blocks and italicises only `Hint:` rather than the entire hint sentence; and in (iii), `question.tex` uses `\widetilde\lambda` rather than `\tilde\lambda`. No mathematical or verbal content is affected.

## Full Solution

### Part (i)

Consider the singular Sturm–Liouville problem
\[
\frac{\mathrm{d}}{\mathrm{d}x}\!\left(x\frac{\mathrm{d}u}{\mathrm{d}x}\right) + \lambda x u = 0, \qquad 0 < x < 1,
\]
with $u(x)$ bounded as $x \to 0^+$ and $u(1) = 0$.

#### Step 1: Proof that every eigenvalue satisfies $\lambda \in \mathbb{R}$ and $\lambda > 0$

Let $u \not\equiv 0$ be a (possibly complex-valued) eigenfunction with eigenvalue $\lambda \in \mathbb{C}$, so $u \in C^2(0, 1]$, $|u(x)| \le M$ on $(0, 1)$ for some constant $M < \infty$, and $u(1) = 0$.

First, we show that $\lim_{x \to 0^+} x u'(x) = 0$. Integrating $(t u'(t))' = -\lambda t u(t)$ over $[x, 1]$ for $0 < x < 1$ gives
\[
u'(1) - x u'(x) = -\lambda \int_x^1 t u(t)\,\mathrm{d}t.
\]
Since $|t u(t)| \le M t$ is integrable on $(0, 1)$, the limit
\[
C := \lim_{x \to 0^+} x u'(x) = u'(1) + \lambda \int_0^1 t u(t)\,\mathrm{d}t
\]
exists and is finite. Suppose for contradiction that $C \neq 0$. Then for sufficiently small $\delta \in (0, 1)$ we have $\operatorname{Re}(\bar{C}\, t u'(t)) \ge \frac{1}{2}|C|^2$ for all $t \in (0, \delta)$, so $\operatorname{Re}(\bar{C}\, u'(t)) \ge \frac{|C|^2}{2t}$. Integrating from $x$ to $\delta$ yields
\[
\operatorname{Re}\!\big(\bar{C}(u(\delta) - u(x))\big) \ge \frac{|C|^2}{2}\ln\!\left(\frac{\delta}{x}\right) \to +\infty \qquad \text{as } x \to 0^+,
\]
contradicting the boundedness of $u$ as $x \to 0^+$. Hence $\lim_{x \to 0^+} x u'(x) = 0$.

Now multiply the ODE $(x u')' + \lambda x u = 0$ by $\overline{u(x)}$ and integrate by parts over $[\varepsilon, 1]$ ($0 < \varepsilon < 1$):
\[
\Big[x u'(x)\overline{u(x)}\Big]_\varepsilon^1 - \int_\varepsilon^1 x |u'(x)|^2\,\mathrm{d}x + \lambda \int_\varepsilon^1 x |u(x)|^2\,\mathrm{d}x = 0.
\]
At $x = 1$, $\overline{u(1)} = 0$. As $\varepsilon \to 0^+$, $|\varepsilon u'(\varepsilon)\overline{u(\varepsilon)}| \le M |\varepsilon u'(\varepsilon)| \to 0$. Letting $\varepsilon \to 0^+$ therefore yields
\[
\lambda \int_0^1 x |u(x)|^2\,\mathrm{d}x = \int_0^1 x |u'(x)|^2\,\mathrm{d}x.
\]
Since $u \not\equiv 0$ is continuous on $(0, 1]$, we have $\int_0^1 x |u(x)|^2\,\mathrm{d}x > 0$ and $\int_0^1 x |u'(x)|^2\,\mathrm{d}x \ge 0$. Thus $\lambda$ must be **real** and **non-negative** ($\lambda \ge 0$).

Moreover, if $\lambda = 0$, then $\int_0^1 x |u'(x)|^2\,\mathrm{d}x = 0$, so $u'(x) \equiv 0$ on $(0, 1)$ and $u$ is constant; the boundary condition $u(1) = 0$ then forces $u \equiv 0$. (Equivalently, directly solving $(x u')' = 0$ gives $u(x) = A\ln x + B$; boundedness as $x \to 0^+$ forces $A = 0$, and $u(1) = 0$ forces $B = 0$.) Hence $\lambda = 0$ is not an eigenvalue, and every eigenvalue must satisfy
\[
\lambda > 0.
\]

#### Step 2: Reduction to Bessel's equation of order zero ($\lambda > 0$)

Write $\lambda = \alpha^2$ with $\alpha = \sqrt{\lambda} > 0$, and make the change of variable $z = \alpha x \in (0, \alpha]$, setting $y(z) = u(z/\alpha) = u(x)$. Since $\frac{\mathrm{d}}{\mathrm{d}x} = \alpha \frac{\mathrm{d}}{\mathrm{d}z}$ and $x = z/\alpha$, the ODE becomes
\[
\alpha \frac{\mathrm{d}}{\mathrm{d}z}\!\left(\frac{z}{\alpha}\cdot \alpha \frac{\mathrm{d}y}{\mathrm{d}z}\right) + \alpha^2 \frac{z}{\alpha} y = 0
\iff \alpha\left[\frac{\mathrm{d}}{\mathrm{d}z}\!\left(z\frac{\mathrm{d}y}{\mathrm{d}z}\right) + z y\right] = 0.
\]
Dividing by $\alpha > 0$, $y(z)$ satisfies Bessel's equation of order zero:
\[
\frac{\mathrm{d}}{\mathrm{d}z}\!\left(z\frac{\mathrm{d}y}{\mathrm{d}z}\right) + z y = 0, \qquad z > 0.
\]
By the given recall, the general solution for $z > 0$ is
\[
y(z) = A J_0(z) + B Y_0(z)
\qquad\Longrightarrow\qquad
u(x) = A J_0(\alpha x) + B Y_0(\alpha x),
\]
where $J_0(z)$ is regular (hence bounded) at $z = 0$ and $Y_0(z)$ is singular (unbounded) as $z \to 0^+$.

Applying the boundary conditions:
1. $u(x)$ bounded as $x \to 0^+$ (i.e. $z = \alpha x \to 0^+$) forces $B = 0$, leaving $u(x) = A J_0(\alpha x)$.
2. $u(1) = 0$ gives $A J_0(\alpha) = 0$.

For a non-trivial solution ($A \neq 0$), we require $J_0(\alpha) = 0$ with $\alpha > 0$. Letting $0 < j_1 < j_2 < j_3 < \cdots$ denote the positive zeros of $J_0(z)$ arranged in ascending order, we obtain $\alpha = j_n$ ($n = 1, 2, \ldots$), and therefore the eigenvalues are
\[
\lambda = \lambda_n = j_n^2 \qquad (n = 1, 2, \ldots),
\]
with corresponding eigenfunctions $u_n(x) = J_0(j_n x)$ (up to a non-zero constant multiple).

---

### Part (ii)

Let $y_\alpha(x) = J_0(\alpha x)$ and $y_\beta(x) = J_0(\beta x)$ (with $\alpha, \beta > 0$, or more generally $\alpha^2 \neq \beta^2$). From part (i), they satisfy
\[
\frac{\mathrm{d}}{\mathrm{d}x}\!\left(x y_\alpha'\right) + \alpha^2 x y_\alpha = 0,
\qquad
\frac{\mathrm{d}}{\mathrm{d}x}\!\left(x y_\beta'\right) + \beta^2 x y_\beta = 0,
\qquad 0 < x < 1.
\]
Multiply the first equation by $y_\beta(x)$ and integrate over $[0, 1]$ (interpreting the lower limit as $\varepsilon \to 0^+$):
\[
\int_0^1 y_\beta (x y_\alpha')'\,\mathrm{d}x + \alpha^2 \int_0^1 y_\alpha y_\beta x\,\mathrm{d}x = 0.
\]
Integrating the first term by parts gives
\[
\Big[x y_\beta(x) y_\alpha'(x)\Big]_0^1 - \int_0^1 x y_\beta'(x) y_\alpha'(x)\,\mathrm{d}x + \alpha^2 \int_0^1 y_\alpha y_\beta x\,\mathrm{d}x = 0.
\]
Similarly, multiplying the equation for $y_\beta$ by $y_\alpha(x)$ and integrating by parts over $[0, 1]$ gives
\[
\Big[x y_\alpha(x) y_\beta'(x)\Big]_0^1 - \int_0^1 x y_\alpha'(x) y_\beta'(x)\,\mathrm{d}x + \beta^2 \int_0^1 y_\alpha y_\beta x\,\mathrm{d}x = 0.
\]
Subtracting the second identity from the first eliminates the symmetric integral $-\int_0^1 x y_\alpha' y_\beta'\,\mathrm{d}x$:
\[
\Big[x\big(y_\beta(x) y_\alpha'(x) - y_\alpha(x) y_\beta'(x)\big)\Big]_0^1 + (\alpha^2 - \beta^2)\int_0^1 y_\alpha(x) y_\beta(x) x\,\mathrm{d}x = 0,
\]
and hence
\[
(\alpha^2 - \beta^2)\int_0^1 J_0(\alpha x) J_0(\beta x) x\,\mathrm{d}x
= \Big[x\big(y_\alpha(x) y_\beta'(x) - y_\beta(x) y_\alpha'(x)\big)\Big]_0^1.
\]
Since $J_0$ is regular at $0$, $y_\alpha, y_\beta, y_\alpha', y_\beta'$ are all bounded at $x = 0$, so the boundary term at $x = 0$ vanishes. At $x = 1$, using the chain rule $y_\alpha'(1) = \alpha J_0'(\alpha)$ and $y_\beta'(1) = \beta J_0'(\beta)$, we obtain
\[
(\alpha^2 - \beta^2)\int_0^1 J_0(\alpha x) J_0(\beta x) x\,\mathrm{d}x
= \beta J_0(\alpha) J_0'(\beta) - \alpha J_0(\beta) J_0'(\alpha).
\]
Dividing by $\alpha^2 - \beta^2$ (valid for $\beta \neq \alpha$ with $\alpha, \beta > 0$, or more generally $\alpha^2 \neq \beta^2$) yields
\[
\int_0^1 J_0(\alpha x) J_0(\beta x) x\,\mathrm{d}x
= \frac{\beta J_0(\alpha) J_0'(\beta) - \alpha J_0(\beta) J_0'(\alpha)}{\alpha^2 - \beta^2} \qquad (\beta \neq \alpha).
\]

#### Orthogonality ($n \neq m$)

Set $\alpha = j_n$ and $\beta = j_m$ with $n \neq m$. Because $0 < j_1 < j_2 < \cdots$, we have $j_n \neq j_m$ and $j_n^2 \neq j_m^2$, while $J_0(j_n) = J_0(j_m) = 0$. Substituting into the formula above immediately gives
\[
\int_0^1 J_0(j_n x) J_0(j_m x) x\,\mathrm{d}x
= \frac{j_m \cdot 0 \cdot J_0'(j_m) - j_n \cdot 0 \cdot J_0'(j_n)}{j_n^2 - j_m^2} = 0 \qquad (n \neq m).
\]

#### Normalisation ($n = m$)

Following the hint, set $\alpha = j_n$ (so $J_0(j_n) = 0$) and $\beta = j_n + \epsilon$ with $0 < |\epsilon| < j_n$. Then
\[
\int_0^1 J_0(j_n x) J_0\big((j_n + \epsilon)x\big) x\,\mathrm{d}x
= \frac{0 - j_n J_0(j_n + \epsilon) J_0'(j_n)}{j_n^2 - (j_n + \epsilon)^2}
= \frac{-j_n J_0'(j_n)\, J_0(j_n + \epsilon)}{-2j_n\epsilon - \epsilon^2}
= \frac{j_n J_0'(j_n)}{2j_n + \epsilon}\cdot \frac{J_0(j_n + \epsilon) - J_0(j_n)}{\epsilon}.
\]
(Equivalently, Taylor expanding $J_0(j_n + \epsilon) = \epsilon J_0'(j_n) + O(\epsilon^2)$ gives $\frac{-j_n J_0'(j_n)[\epsilon J_0'(j_n) + O(\epsilon^2)]}{-2j_n \epsilon - \epsilon^2} = \frac{1}{2}[J_0'(j_n)]^2 + O(\epsilon)$.)

Taking $\epsilon \to 0$: on the left-hand side, continuity of $(x, \epsilon) \mapsto J_0(j_n x) J_0((j_n + \epsilon)x) x$ on $[0, 1] \times [-\delta, \delta]$ justifies passing the limit inside the integral; on the right-hand side, since $j_n > 0$ and $J_0$ is differentiable at $j_n$, $\frac{j_n}{2j_n + \epsilon} \to \frac{1}{2}$ and $\frac{J_0(j_n + \epsilon) - J_0(j_n)}{\epsilon} \to J_0'(j_n)$. Therefore,
\[
\int_0^1 \big[J_0(j_n x)\big]^2 x\,\mathrm{d}x = \frac{1}{2}\big[J_0'(j_n)\big]^2.
\]
(Note $J_0'(j_n) \neq 0$, because $j_n > 0$ is an ordinary point of Bessel's equation so $J_0(j_n) = J_0'(j_n) = 0$ would imply $J_0 \equiv 0$.)

---

### Part (iii)

Consider the inhomogeneous equation
\[
\frac{\mathrm{d}}{\mathrm{d}x}\!\left(x\frac{\mathrm{d}u}{\mathrm{d}x}\right) + \widetilde{\lambda} x u = x f(x), \qquad 0 < x < 1,
\]
where $\widetilde{\lambda}$ is not an eigenvalue (so $\widetilde{\lambda} \neq j_n^2$ for all $n \ge 1$), and both $u$ and $f$ are bounded as $x \to 0^+$ and satisfy $u(1) = f(1) = 0$.

By the assumed completeness of the orthogonal eigenfunctions $u_n(x) = J_0(j_n x)$ ($n = 1, 2, \ldots$) with respect to the weight $w(x) = x$ on $(0, 1)$, we may expand $u$ and $f$ as Fourier–Bessel series:
\[
u(x) = \sum_{n=1}^\infty a_n J_0(j_n x), \qquad f(x) = \sum_{n=1}^\infty c_n J_0(j_n x),
\]
where, by the orthogonality and normalisation formulas from part (ii), the coefficients are
\[
a_n = \frac{\int_0^1 u(x) J_0(j_n x) x\,\mathrm{d}x}{\int_0^1 [J_0(j_n x)]^2 x\,\mathrm{d}x} = \frac{2}{[J_0'(j_n)]^2}\int_0^1 u(x) J_0(j_n x) x\,\mathrm{d}x,
\]
and
\[
c_n = \frac{2}{[J_0'(j_n)]^2}\int_0^1 f(x) J_0(j_n x) x\,\mathrm{d}x.
\]

To relate $a_n$ directly to $c_n$ without assuming term-by-term differentiability of the series for $u$, multiply the inhomogeneous equation by $u_n(x) = J_0(j_n x)$ and integrate over $[\varepsilon, 1]$:
\[
\int_\varepsilon^1 u_n (x u')'\,\mathrm{d}x + \widetilde{\lambda}\int_\varepsilon^1 u u_n x\,\mathrm{d}x = \int_\varepsilon^1 f u_n x\,\mathrm{d}x.
\]
Integrating the first term by parts twice gives
\[
\int_\varepsilon^1 u_n (x u')'\,\mathrm{d}x
= \Big[x\big(u_n u' - u_n' u\big)\Big]_\varepsilon^1 + \int_\varepsilon^1 u (x u_n')'\,\mathrm{d}x.
\]
Let us verify that the boundary terms vanish:
- At $x = 1$: $u_n(1) = J_0(j_n) = 0$ and $u(1) = 0$, so $\left.x(u_n u' - u_n' u)\right|_{x=1} = 0$.
- As $\varepsilon \to 0^+$: since $(x u')' = x(f(x) - \widetilde{\lambda} u(x))$ is bounded on $(0, 1)$ (both $f$ and $u$ are bounded), the exact same argument as in part (i) shows $\lim_{\varepsilon \to 0^+} \varepsilon u'(\varepsilon) = 0$. Since $u_n$, $u_n'$, and $u$ are bounded near $0$, we have $\varepsilon\big(u_n(\varepsilon) u'(\varepsilon) - u_n'(\varepsilon) u(\varepsilon)\big) \to 0$ as $\varepsilon \to 0^+$.

Using $(x u_n')' = -j_n^2 x u_n$ and taking $\varepsilon \to 0^+$ therefore yields
\[
(\widetilde{\lambda} - j_n^2)\int_0^1 u(x) J_0(j_n x) x\,\mathrm{d}x = \int_0^1 f(x) J_0(j_n x) x\,\mathrm{d}x.
\]
Dividing both sides by $\frac{1}{2}[J_0'(j_n)]^2$ gives
\[
(\widetilde{\lambda} - j_n^2) a_n = c_n.
\]
Since $\widetilde{\lambda}$ is not an eigenvalue, $\widetilde{\lambda} - j_n^2 \neq 0$ for every $n \ge 1$, so
\[
a_n = \frac{c_n}{\widetilde{\lambda} - j_n^2} = \frac{2}{(\widetilde{\lambda} - j_n^2)[J_0'(j_n)]^2}\int_0^1 f(t) J_0(j_n t) t\,\mathrm{d}t.
\]
Hence the eigenfunction expansion of $u(x)$ is
\[
u(x) = \sum_{n=1}^\infty \frac{c_n}{\widetilde{\lambda} - j_n^2}\, J_0(j_n x),
\qquad \text{where}\quad c_n = \frac{2}{[J_0'(j_n)]^2}\int_0^1 f(t) J_0(j_n t) t\,\mathrm{d}t,
\]
or equivalently, written out explicitly,
\[
u(x) = \sum_{n=1}^\infty \frac{2}{(\widetilde{\lambda} - j_n^2)[J_0'(j_n)]^2}\left(\int_0^1 f(t) J_0(j_n t) t\,\mathrm{d}t\right) J_0(j_n x).
\]

## Final Answers

- **(i)** Every eigenvalue is strictly positive, and the eigenvalues are $\lambda_n = j_n^2$ ($n = 1, 2, \ldots$), where $0 < j_1 < j_2 < \cdots$ are the positive zeros of $J_0(z)$, with eigenfunctions $u_n(x) = J_0(j_n x)$.
- **(ii)** All three integral identities are established:
  \[
  \int_0^1 J_0(\alpha x) J_0(\beta x) x\,\mathrm{d}x = \frac{\beta J_0(\alpha) J_0'(\beta) - \alpha J_0(\beta) J_0'(\alpha)}{\alpha^2 - \beta^2} \quad (\alpha \neq \pm\beta),
  \]
  \[
  \int_0^1 J_0(j_n x) J_0(j_m x) x\,\mathrm{d}x = 0 \quad (n \neq m), \qquad
  \int_0^1 \big[J_0(j_n x)\big]^2 x\,\mathrm{d}x = \frac{1}{2}\big[J_0'(j_n)\big]^2.
  \]
- **(iii)** The eigenfunction expansion of $u(x)$ is
  \[
  u(x) = \sum_{n=1}^\infty \frac{c_n}{\widetilde{\lambda} - j_n^2}\, J_0(j_n x),
  \qquad c_n = \frac{2}{[J_0'(j_n)]^2}\int_0^1 f(t) J_0(j_n t) t\,\mathrm{d}t.
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_qx06_check.py` (exit code 0).

```text
Part (i) checks passed: first 5 zeros = [2.404826, 5.520078, 8.653728, 11.791534, 14.930918], max_ode_err=4.26e-14.
Part (ii) checks passed: max_gen_err=3.92e-13, max_orth_err=8.31e-13, eps-limit error=5.60e-08.
  Part (iii) test (lam_tilde=  3.5, 20 modes): max |u_20(x) - u_exact(x)| = 3.65e-06
  Part (iii) test (lam_tilde= 15.0, 20 modes): max |u_20(x) - u_exact(x)| = 3.65e-06
  Part (iii) test (lam_tilde= -4.0, 20 modes): max |u_20(x) - u_exact(x)| = 3.65e-06
ALL CHECKS PASSED FOR IB/methods/sheet-1/qx06.
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  1. **Ruling out $\lambda \le 0$ (and non-real $\lambda$) in part (i):** A very frequent gap is jumping straight to $z = \sqrt{\lambda}\, x$ without proving $\lambda > 0$. Checking $\lim_{x \to 0^+} x u'(x) = 0$ from boundedness of $u$ and integrating $[(x u')' + \lambda x u]\bar{u}$ by parts rigorously rules out $\lambda \notin \mathbb{R}$, $\lambda < 0$, and $\lambda = 0$ without needing any unstated facts about modified Bessel functions $I_0, K_0$.
  2. **Sign of the denominator in part (iii):** Because the ODE is $+\frac{\mathrm{d}}{\mathrm{d}x}(x\frac{\mathrm{d}u}{\mathrm{d}x}) + \widetilde{\lambda} x u = +x f(x)$ and $(x u_n')' = -j_n^2 x u_n$, the coefficient is $a_n = \frac{c_n}{\widetilde{\lambda} - j_n^2}$, NOT $\frac{c_n}{j_n^2 - \widetilde{\lambda}}$! Our numerical test with $\widetilde{\lambda} = 15.0 \in (j_1^2, j_2^2)$ explicitly verifies both signs of $\widetilde{\lambda} - j_n^2$.
  3. **Factor of $2/[J_0'(j_n)]^2$ in $c_n$:** Because $\int_0^1 [J_0(j_n x)]^2 x\,\mathrm{d}x = \frac{1}{2}[J_0'(j_n)]^2$, watch for missing the factor of $2$ (or misplacing it in the denominator) when writing $c_n$.
  4. **Avoiding unproved term-by-term differentiation in part (iii):** Multiplying the ODE by $J_0(j_n x)$ and integrating by parts twice derives $(\widetilde{\lambda} - j_n^2) a_n = c_n$ directly from the ODE and BCs rather than formally differentiating $\sum a_n J_0(j_n x)$ twice term-by-term.
