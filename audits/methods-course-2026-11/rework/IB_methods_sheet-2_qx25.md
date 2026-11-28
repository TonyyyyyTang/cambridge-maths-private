---
id: IB/methods/sheet-2/qx25
slug: IB_methods_sheet-2_qx25
source_archive: private/raw/sheets/IB/methods/2004-05/sheet-2.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-2_qx25_p2.png` (Question 8, *Bessel's equation*, on page 2 of the 2004–05 Part IB Methods Example Sheet 2). Question 9 follows immediately below Question 8 on page 2, so there is no spillover.

Every word, symbol, sign, subscript, superscript, equation, note, hint, and subpart (i)–(iii) matches the printed PDF exactly. Purely typographic differences only: underlined italics for the heading *Bessel's equation.*, inline placement of the hint beside the second display integral in (ii), and `\widetilde\lambda` vs `\tilde{\lambda}` in (iii). Verdict: `transcription: ok`.

## Full Solution

### Part (i): Eigenvalues of the Sturm–Liouville problem

We seek all $\lambda$ for which the singular Sturm–Liouville problem
\[
\frac{d}{dx}\!\left(x\frac{du}{dx}\right) + \lambda x u = 0, \qquad 0 < x < 1,
\]
has a non-trivial solution $u(x)$ bounded as $x \to 0^+$ with $u(1) = 0$.

#### Step 1: Every eigenvalue is real and strictly positive ($\lambda > 0$)

Let $u(x)$ be any (possibly complex-valued) non-trivial solution bounded as $x \to 0^+$ with $u(1) = 0$. First we verify that the singular boundary term at $x = 0$ vanishes. Integrating $(xu')' = -\lambda xu$ from $\varepsilon$ to $x$ ($0 < \varepsilon < x < 1$):
\[
x u'(x) - \varepsilon u'(\varepsilon) = -\lambda \int_\varepsilon^x t u(t)\,dt.
\]
Since $u$ is bounded on $(0,1)$, say $|u(t)| \le M$, the integral $\int_0^x t u(t)\,dt$ converges and is $O(x^2)$. Hence $\varepsilon u'(\varepsilon)$ approaches a finite limit $C = x u'(x) + \lambda\int_0^x t u(t)\,dt$ as $\varepsilon \to 0^+$, so
\[
u'(x) = \frac{C}{x} + O(x) \implies u(1) - u(x) = -C\log x + O(1).
\]
Boundedness of $u(x)$ as $x \to 0^+$ forces $C = 0$, whence $\lim_{\varepsilon \to 0^+} \varepsilon u'(\varepsilon) = 0$ (indeed $u'(x) = O(x)$ as $x \to 0^+$).

Now multiply $(xu')' + \lambda xu = 0$ by $\overline{u}(x)$, integrate over $[\varepsilon, 1]$, and integrate the first term by parts:
\[
\left[x\overline{u}(x)u'(x)\right]_\varepsilon^1 - \int_\varepsilon^1 x|u'(x)|^2\,dx + \lambda\int_\varepsilon^1 x|u(x)|^2\,dx = 0.
\]
At $x = 1$, $u(1) = 0$; as $\varepsilon \to 0^+$, $\overline{u}(\varepsilon)$ is bounded and $\varepsilon u'(\varepsilon) \to 0$, so the boundary term vanishes in the limit $\varepsilon \to 0^+$. Since $u'(x) = O(x)$ near $0$, $\int_0^1 x|u'|^2\,dx < \infty$, and since $u \not\equiv 0$ is continuous on $(0,1]$, $\int_0^1 x|u|^2\,dx > 0$. Thus
\[
\lambda = \frac{\int_0^1 x|u'(x)|^2\,dx}{\int_0^1 x|u(x)|^2\,dx},
\]
which shows that $\lambda$ is real and $\lambda \ge 0$.

If $\lambda = 0$, then $\int_0^1 x|u'(x)|^2\,dx = 0$, forcing $u'(x) \equiv 0$ on $(0,1)$, so $u(x)$ is constant; $u(1) = 0$ then forces $u \equiv 0$, contradiction. (Equivalently, solving $(xu')' = 0$ directly gives $u(x) = A\log x + B$; boundedness at $0$ gives $A = 0$ and $u(1) = 0$ gives $B = 0$.)

Hence every eigenvalue must satisfy $\lambda > 0$. (One may also check real $\lambda = -\mu^2 < 0$ directly: changing variable to $z = i\mu x$ gives the modified Bessel function $I_0(\mu x) = \sum_{m=0}^\infty \frac{(\mu x)^{2m}}{2^{2m}(m!)^2} \ge 1 > 0$ as the unique bounded-at-zero solution up to scale, which never vanishes at $x = 1$.)

#### Step 2: Reduction to Bessel's equation of order zero for $\lambda > 0$

Write $\lambda = k^2$ with $k > 0$, and change independent variable to $z = kx$ ($z > 0$). Since $\frac{d}{dx} = k\frac{d}{dz}$ and $x = z/k$, the ODE becomes
\[
k\frac{d}{dz}\!\left(\frac{z}{k}\cdot k\frac{du}{dz}\right) + k^2\frac{z}{k}u = k\left[\frac{d}{dz}\!\left(z\frac{du}{dz}\right) + zu\right] = 0.
\]
Dividing by $k > 0$, $y(z) = u(z/k)$ satisfies Bessel's equation of order zero:
\[
\frac{d}{dz}\!\left(z\frac{dy}{dz}\right) + zy = 0, \qquad z > 0.
\]
By the given note, its general solution is
\[
u(x) = A J_0(kx) + B\bigl[g(kx) + J_0(kx)\log(kx)\bigr],
\]
where $A, B$ are constants, $J_0(z) = \sum_{m=0}^\infty (-1)^m \frac{z^{2m}}{2^{2m}(m!)^2}$ is smooth everywhere with $J_0(0) = 1$, and $g(z)$ is regular (hence bounded) as $z \to 0^+$.

As $x \to 0^+$, $J_0(kx) \to 1$ and $g(kx)$ stays finite, while $\log(kx) \to -\infty$, so $J_0(kx)\log(kx) \sim \log(kx) \to -\infty$. Therefore boundedness of $u(x)$ as $x \to 0^+$ requires $B = 0$, leaving
\[
u(x) = A J_0(kx),
\]
with $A \neq 0$ for a non-trivial eigenfunction. Applying $u(1) = 0$ yields
\[
J_0(k) = 0 \qquad (k > 0).
\]
Hence $k = j_n$ ($n = 1, 2, \dots$), where $0 < j_1 < j_2 < \dots$ are the positive zeros of $J_0(z)$ arranged in ascending order (note $J_0(0) = 1 \neq 0$). Consequently, the eigenvalues are
\[
\lambda = \lambda_n = j_n^2 \qquad (n = 1, 2, \dots),
\]
with corresponding eigenfunctions $u_n(x) = J_0(j_n x)$ (up to a non-zero constant multiple).

### Part (ii): Lommel integral identity, orthogonality, and normalisation

Let $u(x) = J_0(\alpha x)$ and $v(x) = J_0(\beta x)$ (with $\alpha, \beta > 0$ and $\beta \neq \alpha$; more generally $\alpha^2 \neq \beta^2$). From part (i), $u$ and $v$ are smooth on $[0, 1]$ and satisfy
\[
\frac{d}{dx}\!\left(x u'\right) + \alpha^2 x u = 0, \qquad \frac{d}{dx}\!\left(x v'\right) + \beta^2 x v = 0.
\]
Multiply the first equation by $v(x)$, the second by $u(x)$, and integrate from $x = 0$ to $x = 1$:
\[
\int_0^1 v(x)\frac{d}{dx}\!\left(x u'(x)\right)\,dx + \alpha^2 \int_0^1 x u(x)v(x)\,dx = 0,
\]
\[
\int_0^1 u(x)\frac{d}{dx}\!\left(x v'(x)\right)\,dx + \beta^2 \int_0^1 x u(x)v(x)\,dx = 0.
\]
Integrating by parts in both equations gives
\[
\left[x v u'\right]_0^1 - \int_0^1 x u' v'\,dx + \alpha^2 \int_0^1 x u v\,dx = 0,
\]
\[
\left[x u v'\right]_0^1 - \int_0^1 x u' v'\,dx + \beta^2 \int_0^1 x u v\,dx = 0.
\]
Subtracting the second equation from the first eliminates the $-\int_0^1 x u' v'\,dx$ integral:
\[
(\alpha^2 - \beta^2)\int_0^1 x u(x)v(x)\,dx = \left[x\bigl(u(x)v'(x) - v(x)u'(x)\bigr)\right]_0^1.
\]
Since $u, v, u', v'$ are smooth (bounded) at $x = 0$, the boundary term at $x = 0$ vanishes. Using $u'(x) = \alpha J_0'(\alpha x)$ and $v'(x) = \beta J_0'(\beta x)$ at $x = 1$,
\[
(\alpha^2 - \beta^2)\int_0^1 J_0(\alpha x)J_0(\beta x)x\,dx = \beta J_0(\alpha)J_0'(\beta) - \alpha J_0(\beta)J_0'(\alpha).
\]
Dividing by $\alpha^2 - \beta^2 \neq 0$ yields the required identity:
\[
\int_0^1 J_0(\alpha x)J_0(\beta x)x\,dx = \frac{\beta J_0(\alpha)J_0'(\beta) - \alpha J_0(\beta)J_0'(\alpha)}{\alpha^2 - \beta^2} \qquad (\beta \neq \alpha).
\]

**Orthogonality ($n \neq m$):**

Set $\alpha = j_n$ and $\beta = j_m$ with $n \neq m$. Since the positive zeros are distinct ($j_n \neq j_m$, hence $j_n^2 - j_m^2 \neq 0$) and $J_0(j_n) = J_0(j_m) = 0$, the numerator vanishes identically, giving
\[
\int_0^1 J_0(j_n x)J_0(j_m x)x\,dx = 0 \qquad (n \neq m).
\]

**Normalisation ($n = m$):**

Following the hint, set $\alpha = j_n$ and $\beta = j_n + \epsilon$ with $\epsilon \neq 0$ sufficiently small (so $\alpha^2 \neq \beta^2$). Using $J_0(j_n) = 0$, the identity reduces to
\[
\int_0^1 J_0(j_n x)J_0\bigl((j_n + \epsilon)x\bigr)x\,dx
= \frac{-j_n J_0(j_n + \epsilon)J_0'(j_n)}{j_n^2 - (j_n + \epsilon)^2}
= \frac{-j_n J_0'(j_n)\bigl(J_0(j_n + \epsilon) - J_0(j_n)\bigr)}{-\epsilon(2j_n + \epsilon)}
= \frac{j_n J_0'(j_n)}{2j_n + \epsilon}\cdot\frac{J_0(j_n + \epsilon) - J_0(j_n)}{\epsilon}.
\]
Take the limit $\epsilon \to 0$. On the left-hand side, $J_0((j_n + \epsilon)x) \to J_0(j_n x)$ uniformly on $x \in [0, 1]$ by continuity of $J_0$, so the limit may be taken inside the integral. On the right-hand side, $\frac{J_0(j_n+\epsilon)-J_0(j_n)}{\epsilon} \to J_0'(j_n)$ and $\frac{j_n}{2j_n+\epsilon} \to \frac{1}{2}$ (since $j_n > 0$). Hence
\[
\int_0^1 [J_0(j_n x)]^2 x\,dx = \frac{1}{2}[J_0'(j_n)]^2.
\]
(Note that this also confirms $J_0'(j_n) \neq 0$ for every $n$, since the integral on the left is strictly positive.)

### Part (iii): Eigenfunction expansion of the inhomogeneous problem

We solve
\[
\frac{d}{dx}\!\left(x\frac{du}{dx}\right) + \widetilde\lambda x u = x f(x), \qquad 0 < x < 1,
\]
with $u(x)$ bounded as $x \to 0^+$ and $u(1) = 0$, where $\widetilde\lambda \neq j_n^2$ for all $n \ge 1$.

Assuming completeness of $\{J_0(j_n x)\}_{n=1}^\infty$, write the eigenfunction expansion
\[
u(x) = \sum_{n=1}^\infty c_n J_0(j_n x).
\]
By the orthogonality and normalisation relations in part (ii) (multiplying by $x J_0(j_m x)$ and integrating over $[0,1]$), the expansion coefficients are
\[
c_n = \frac{\int_0^1 x u(x) J_0(j_n x)\,dx}{\int_0^1 x [J_0(j_n x)]^2\,dx} = \frac{2}{[J_0'(j_n)]^2}\int_0^1 x u(x) J_0(j_n x)\,dx.
\]

To determine $\int_0^1 x u(x) J_0(j_n x)\,dx$ directly without assuming term-by-term differentiation of the series for $u$, multiply the inhomogeneous ODE by $J_0(j_n x)$ and integrate over $[\varepsilon, 1]$:
\[
\int_\varepsilon^1 J_0(j_n x)\frac{d}{dx}\!\left(x u'(x)\right)\,dx + \widetilde\lambda \int_\varepsilon^1 x u(x) J_0(j_n x)\,dx = \int_\varepsilon^1 x f(x) J_0(j_n x)\,dx.
\]
Integrating the first term by parts twice yields
\[
\int_\varepsilon^1 J_0(j_n x)\frac{d}{dx}\!\left(x u'(x)\right)\,dx
= \left[x\left(J_0(j_n x)u'(x) - u(x)\frac{d}{dx}J_0(j_n x)\right)\right]_\varepsilon^1
+ \int_\varepsilon^1 u(x)\frac{d}{dx}\!\left(x\frac{d}{dx}J_0(j_n x)\right)\,dx.
\]
Check the boundary terms:
- At $x = 1$: $J_0(j_n) = 0$ and $u(1) = 0$, so both terms vanish.
- As $\varepsilon \to 0^+$: $J_0(j_n \varepsilon) \to 1$, $\frac{d}{dx}J_0(j_n x)\big|_{x=\varepsilon} = j_n J_0'(j_n\varepsilon) \to 0$, $u(\varepsilon)$ is bounded, and by the same integration argument as in part (i) (for integrable $x f(x)$), $\varepsilon u'(\varepsilon) \to 0$. Thus the lower boundary term vanishes as $\varepsilon \to 0^+$.

In the remaining integral, $\frac{d}{dx}\!\left(x\frac{d}{dx}J_0(j_n x)\right) = -j_n^2 x J_0(j_n x)$. Letting $\varepsilon \to 0^+$ gives
\[
(\widetilde\lambda - j_n^2)\int_0^1 x u(x) J_0(j_n x)\,dx = \int_0^1 x f(x) J_0(j_n x)\,dx.
\]
Since $\widetilde\lambda$ is not an eigenvalue, $\widetilde\lambda - j_n^2 \neq 0$ for every $n \ge 1$. Dividing by $\widetilde\lambda - j_n^2$ and substituting into the formula for $c_n$ (using a dummy integration variable $\xi$ to distinguish it from $x$) gives
\[
c_n = \frac{2}{(\widetilde\lambda - j_n^2)[J_0'(j_n)]^2}\int_0^1 \xi f(\xi) J_0(j_n \xi)\,d\xi,
\]
and therefore the eigenfunction expansion of $u(x)$ is
\[
u(x) = \sum_{n=1}^\infty \frac{2 J_0(j_n x)}{(\widetilde\lambda - j_n^2)[J_0'(j_n)]^2}\int_0^1 \xi f(\xi) J_0(j_n \xi)\,d\xi.
\]
(Equivalently, writing $f(x) = \sum_{n=1}^\infty f_n J_0(j_n x)$ with $f_n = \frac{2}{[J_0'(j_n)]^2}\int_0^1 \xi f(\xi)J_0(j_n\xi)\,d\xi$, we have $c_n = \frac{f_n}{\widetilde\lambda - j_n^2}$.)

## Final Answers

- **(i)** Eigenvalues $\lambda_n = j_n^2$ ($n = 1, 2, \dots$), where $0 < j_1 < j_2 < \dots$ are the positive zeros of $J_0(z)$; eigenfunctions $u_n(x) = J_0(j_n x)$. (No eigenvalues with $\lambda \le 0$ or $\lambda \notin \mathbf{R}$.)
- **(ii)** General identity proved by integration by parts; orthogonality $\int_0^1 J_0(j_n x)J_0(j_m x)x\,dx = 0$ ($n \neq m$); normalisation $\int_0^1 [J_0(j_n x)]^2 x\,dx = \frac{1}{2}[J_0'(j_n)]^2$.
- **(iii)** Eigenfunction expansion of $u(x)$:
  \[
  u(x) = \sum_{n=1}^\infty c_n J_0(j_n x), \qquad c_n = \frac{2}{(\widetilde\lambda - j_n^2)[J_0'(j_n)]^2}\int_0^1 \xi f(\xi) J_0(j_n \xi)\,d\xi.
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-2_qx25_check.py`

Verbatim stdout (exit code 0):
```text
IB_methods_sheet-2_qx25_check: ALL PASSED
  first 3 zeros: j_1=2.4048255577, j_2=5.5200781103, j_3=8.6537279129
  part (i) Bessel ODE err: 4.805e-09
  part (ii) Lommel err: 6.956e-16, orth err: 1.653e-16, norm err: 1.041e-15, hint-limit err: 6.452e-11
  part (iii) errors: Case A (lam=0, u=(x^2-1)/4)=3.061e-05, Case B (lam=10, u=(1-x^2)^3)=2.393e-06, Case C (lam=-3, u=1-x^4)=4.894e-04, Case D (single mode)=1.249e-16
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`. Every analytical derivation is airtight and independently verified to machine precision for the Lommel/orthogonality/normalisation identities and against four exact closed-form inhomogeneous test problems (covering $\widetilde\lambda = 0$, $\widetilde\lambda > 0$, $\widetilde\lambda < 0$, and single-mode forcing).
- **Phase B watchpoints:**
  - In part (i), watch whether the existing solution explicitly excludes $\lambda = 0$ and $\lambda < 0$ (many student solutions assume $\lambda > 0$ without justification!), and verify careful handling of the singular endpoint $x = 0$ (justifying $x u'(x) \to 0$ if using an integration-by-parts energy proof).
  - In part (ii), watch the sign in $\frac{\beta J_0(\alpha)J_0'(\beta) - \alpha J_0(\beta)J_0'(\alpha)}{\alpha^2 - \beta^2}$ and check that the $\epsilon \to 0$ limit correctly uses $j_n^2 - (j_n+\epsilon)^2 = -\epsilon(2j_n+\epsilon)$ to give $+\frac{1}{2}[J_0'(j_n)]^2$.
  - In part (iii), watch the denominator sign: it MUST be $\widetilde\lambda - j_n^2$ (NOT $j_n^2 - \widetilde\lambda$!) because the PDE is $(xu')' + \widetilde\lambda xu = xf(x)$ with $+xf(x)$ on the RHS and $(x J_0(j_n x)')' = -j_n^2 x J_0(j_n x)$. Also check that the weight $\xi$ is present in $\int_0^1 \xi f(\xi) J_0(j_n \xi)\,d\xi$ and the factor $2/[J_0'(j_n)]^2$ is not inverted or missing the $2$.
