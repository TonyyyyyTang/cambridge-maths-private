---
id: IB/methods/sheet-2/qx03
slug: IB_methods_sheet-2_qx03
source_archive: private/raw/sheets/IB/methods/2016-17/sheet-2.pdf
pages_checked: [1, 2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-2_qx03_p1.png` (page 1 of `2016-17/sheet-2.pdf`, Question 3) and confirmed via `IB_methods_sheet-2_qx03_p2.png` that nothing spills over (Question 4 immediately follows Question 3 on page 1).

- **Question 3 body:** Exact verbatim match to all words, punctuation, and mathematical symbols in the printed sheet.
- **Purely typographic difference:** In the printed PDF, the entire bracketed hint (`[Hint: Recall the orthogonality relation` $\int_{-1}^1 P_\ell(x) P_m(x)\,dx = \frac{2}{2\ell+1}\delta_{\ell m}$.`]`) has italicized text (including the brackets and `Recall the orthogonality relation`), whereas `question.tex` italicizes only `\emph{Hint:}` and leaves the remaining words in roman type. No mathematical content or wording is affected.

## Full Solution

### Part 1: Showing $1/|\mathbf{r} - \mathbf{r}'|$ satisfies Laplace's equation for $\mathbf{r} \neq \mathbf{r}'$

We work in three dimensions ($\mathbb{R}^3$). Fix $\mathbf{r}' \in \mathbb{R}^3$ and set $\mathbf{R} = \mathbf{r} - \mathbf{r}'$, with components $R_i = x_i - x'_i$ ($i = 1, 2, 3$) and length $R = |\mathbf{r} - \mathbf{r}'| = (R_j R_j)^{1/2}$. Whenever $\mathbf{r} \neq \mathbf{r}'$, we have $R > 0$, so $1/R$ is smooth.

Differentiating $R^2 = R_j R_j$ with respect to $x_i$ (using the summation convention over $i, j \in \{1, 2, 3\}$):
\[
2R \frac{\partial R}{\partial x_i} = 2 R_j \frac{\partial R_j}{\partial x_i} = 2 R_j \delta_{ij} = 2 R_i
\implies \frac{\partial R}{\partial x_i} = \frac{R_i}{R}.
\]
Hence the gradient of $1/R$ is
\[
\frac{\partial}{\partial x_i}\left(\frac{1}{R}\right) = -\frac{1}{R^2}\frac{\partial R}{\partial x_i} = -\frac{R_i}{R^3}.
\]
Taking the divergence (i.e. applying $\dfrac{\partial}{\partial x_i}$ and summing over $i = 1, 2, 3$):
\[
\nabla^2 \left(\frac{1}{R}\right)
= \frac{\partial}{\partial x_i}\left(-\frac{R_i}{R^3}\right)
= -\frac{1}{R^3}\frac{\partial R_i}{\partial x_i} - R_i \left(-\frac{3}{R^4}\frac{\partial R}{\partial x_i}\right).
\]
In $\mathbb{R}^3$, $\dfrac{\partial R_i}{\partial x_i} = \delta_{ii} = 3$, and $R_i \dfrac{\partial R}{\partial x_i} = \dfrac{R_i R_i}{R} = \dfrac{R^2}{R} = R$. Therefore,
\[
\nabla^2 \left(\frac{1}{|\mathbf{r} - \mathbf{r}'|}\right) = -\frac{3}{R^3} + \frac{3R}{R^4} = -\frac{3}{R^3} + \frac{3}{R^3} = 0
\]
for all $\mathbf{r} \neq \mathbf{r}'$.

*(Equivalently, in spherical polar coordinates centred at $\mathbf{r}'$, $\nabla^2(R^{-1}) = \frac{1}{R^2}\frac{d}{dR}\bigl(R^2 \frac{d}{dR}(R^{-1})\bigr) = \frac{1}{R^2}\frac{d}{dR}(-1) = 0$. By symmetry, $1/|\mathbf{r}-\mathbf{r}'|$ also satisfies Laplace's equation with respect to $\mathbf{r}'$.)*

### Part 2: Potential inside a sphere bounded by two hemispheres at $\pm V$

Let the spherical region have radius $a > 0$ and centre at the origin. Choose spherical polar coordinates $(r, \theta, \phi)$ such that the polar axis ($\theta = 0$) is perpendicular to the equatorial plane separating the two hemispheres, with the boundary potential given by
\[
\Phi(a, \theta, \phi) = \begin{cases}
+V, & 0 \le \theta < \dfrac{\pi}{2}, \\[4pt]
-V, & \dfrac{\pi}{2} < \theta \le \pi,
\end{cases}
\]
(assuming a vanishingly thin insulating gap at the equator $\theta = \pi/2$). Inside the sphere ($0 \le r < a$), the electrostatic potential satisfies Laplace's equation $\nabla^2 \Phi = 0$ and must be finite (regular) throughout the interior, including at the origin $r = 0$ and along the polar axis $\theta = 0, \pi$.

#### Step 1: Axisymmetric separation of variables

Since both the spherical domain and the boundary data are independent of the azimuthal angle $\phi$, the unique bounded harmonic solution inside $r < a$ is axisymmetric: $\Phi = \Phi(r, \theta)$. Laplace's equation in spherical polar coordinates then becomes
\[
\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2 \frac{\partial \Phi}{\partial r}\right)
+ \frac{1}{r^2 \sin\theta}\frac{\partial}{\partial \theta}\left(\sin\theta \frac{\partial \Phi}{\partial \theta}\right) = 0.
\]
Seeking separable solutions $\Phi(r, \theta) = R(r)\Theta(\theta)$ and multiplying by $r^2/(R\Theta)$ gives
\[
\frac{1}{R}\frac{d}{dr}\left(r^2 \frac{dR}{dr}\right)
= -\frac{1}{\Theta \sin\theta}\frac{d}{d\theta}\left(\sin\theta \frac{d\Theta}{d\theta}\right) = \lambda,
\]
where $\lambda$ is a separation constant.

- **Angular equation:** Set $x = \cos\theta \in [-1, 1]$ and $y(x) = \Theta(\theta)$. Since $\dfrac{d}{d\theta} = -\sin\theta\dfrac{d}{dx}$ and $\sin^2\theta = 1 - x^2$, the angular equation becomes Legendre's differential equation:
  \[
  \frac{d}{dx}\left((1 - x^2)\frac{dy}{dx}\right) + \lambda y = 0, \qquad -1 < x < 1.
  \]
  Regularity (boundedness) at the singular endpoints $x = \pm 1$ ($\theta = 0, \pi$) requires $\lambda = \ell(\ell + 1)$ for $\ell \in \{0, 1, 2, \dots\}$, with regular solutions proportional to the Legendre polynomials $P_\ell(x) = P_\ell(\cos\theta)$ (the second linearly independent solutions $Q_\ell(x)$ diverge logarithmically as $x \to \pm 1$).

- **Radial equation:** With $\lambda = \ell(\ell + 1)$,
  \[
  r^2 R'' + 2r R' - \ell(\ell + 1) R = 0.
  \]
  Trying $R(r) = r^\alpha$ gives $\alpha(\alpha + 1) = \ell(\ell + 1)$, with distinct roots $\alpha = \ell$ and $\alpha = -(\ell + 1)$, so $R_\ell(r) = A_\ell r^\ell + B_\ell r^{-(\ell + 1)}$. Regularity at the origin $r = 0$ inside the sphere forces $B_\ell = 0$ for every $\ell \ge 0$.

*(Connection with Part 1: Expanding $1/|\mathbf{r} - \mathbf{r}'| = \sum_{\ell=0}^\infty r^\ell (r')^{-(\ell+1)} P_\ell(\cos\theta)$ for $\mathbf{r}' = r'\hat{\mathbf{z}}$ with $r' > r$ immediately reflects why each solid harmonic $r^\ell P_\ell(\cos\theta)$ satisfies Laplace's equation.)*

Superposing all regular axisymmetric modes, we write the general interior potential as
\[
\Phi(r, \theta) = \sum_{\ell = 0}^{\infty} c_\ell \left(\frac{r}{a}\right)^\ell P_\ell(\cos\theta), \qquad 0 \le r < a.
\]

#### Step 2: Matching the boundary condition via Legendre orthogonality

At $r = a$, setting $x = \cos\theta \in [-1, 1]$, the boundary data is
\[
g(x) = \begin{cases}
+V, & 0 < x \le 1, \\
-V, & -1 \le x < 0,
\end{cases}
\qquad
\sum_{\ell = 0}^{\infty} c_\ell P_\ell(x) = g(x).
\]
Multiplying both sides by $P_m(x)$, integrating over $[-1, 1]$, and applying the orthogonality relation $\int_{-1}^{1} P_\ell(x) P_m(x)\,dx = \dfrac{2}{2\ell + 1}\delta_{\ell m}$ yields
\[
\frac{2}{2\ell + 1} c_\ell = \int_{-1}^{1} g(x) P_\ell(x)\,dx = V\int_{0}^{1} P_\ell(x)\,dx - V\int_{-1}^{0} P_\ell(x)\,dx.
\]
Using the parity of the Legendre polynomials, $P_\ell(-x) = (-1)^\ell P_\ell(x)$, the substitution $x \mapsto -x$ gives $\int_{-1}^{0} P_\ell(x)\,dx = (-1)^\ell \int_{0}^{1} P_\ell(x)\,dx$. Therefore,
\[
\frac{2}{2\ell + 1} c_\ell = V\bigl(1 - (-1)^\ell\bigr)\int_{0}^{1} P_\ell(x)\,dx.
\]
- **Even $\ell = 2k$ ($k = 0, 1, 2, \dots$):** $1 - (-1)^{2k} = 0$, so
  \[
  c_{2k} = 0.
  \]
  (This also follows directly from the antisymmetry of $g(x)$ under $x \mapsto -x$.)
- **Odd $\ell = 2k + 1$ ($k = 0, 1, 2, \dots$):** $1 - (-1)^{2k+1} = 2$, so
  \[
  c_{2k+1} = (4k + 3) V \int_{0}^{1} P_{2k+1}(x)\,dx.
  \]

#### Step 3: Exact evaluation of $\int_0^1 P_{2k+1}(x)\,dx$

We evaluate the integral using Rodrigues' formula,
\[
P_\ell(x) = \frac{1}{2^\ell \ell!}\frac{d^\ell}{dx^\ell}\bigl((x^2 - 1)^\ell\bigr).
\]
For any $\ell \ge 1$:
\[
\int_{0}^{1} P_\ell(x)\,dx = \frac{1}{2^\ell \ell!} \left[ \frac{d^{\ell - 1}}{dx^{\ell - 1}}\bigl((x^2 - 1)^\ell\bigr) \right]_{x = 0}^{x = 1}.
\]
- At $x = 1$, $(x^2 - 1)^\ell = (x - 1)^\ell (x + 1)^\ell$ has a zero of order $\ell$, so its $(\ell - 1)$-th derivative vanishes at $x = 1$.
- At $x = 0$, expand $(x^2 - 1)^\ell$ by the binomial theorem:
  \[
  (x^2 - 1)^\ell = \sum_{j = 0}^{\ell} \binom{\ell}{j} (-1)^{\ell - j} x^{2j}.
  \]
  For $\ell = 2k + 1$, we take the $(2k)$-th derivative at $x = 0$, so only the $j = k$ term ($x^{2k}$) survives, contributing
  \[
  (2k)! \binom{2k + 1}{k} (-1)^{(2k + 1) - k} = (-1)^{k + 1} (2k)! \frac{(2k + 1)!}{k! (k + 1)!}.
  \]

Remembering the minus sign from the lower integration limit $x = 0$:
\[
\int_{0}^{1} P_{2k + 1}(x)\,dx
= -\frac{1}{2^{2k + 1} (2k + 1)!} \cdot (-1)^{k + 1} (2k)! \frac{(2k + 1)!}{k! (k + 1)!}
= \frac{(-1)^k (2k)!}{2^{2k + 1} k! (k + 1)!}.
\]
*(Equivalently, using $(2\ell+1)P_\ell(x) = P'_{\ell+1}(x) - P'_{\ell-1}(x)$ and $P_n(1) = 1$, we get $(2\ell+1)\int_0^1 P_\ell(x)\,dx = P_{\ell-1}(0) - P_{\ell+1}(0) = (-1)^k \frac{(4k+3)(2k)!}{2^{2k+1}k!(k+1)!}$, matching identically.)*

Therefore, for $k = 0, 1, 2, \dots$:
\[
c_{2k + 1} = V (-1)^k \frac{(4k + 3)(2k)!}{2^{2k + 1} k! (k + 1)!}
= V (-1)^k (4k + 3) \frac{(2k - 1)!!}{(2k + 2)!!}
= V \bigl(P_{2k}(0) - P_{2k+2}(0)\bigr),
\]
with the convention $(-1)!! = 1$. The first few non-zero coefficients are:
\[
c_1 = \frac{3}{2}V, \qquad c_3 = -\frac{7}{8}V, \qquad c_5 = \frac{11}{16}V, \qquad c_7 = -\frac{75}{128}V.
\]

#### Step 4: Final potential and convergence check

Substituting $c_\ell$ into the interior expansion, the potential inside the sphere $r < a$ is
\[
\Phi(r, \theta) = V \sum_{k = 0}^{\infty} (-1)^k \frac{(4k + 3)(2k)!}{2^{2k + 1} k! (k + 1)!} \left(\frac{r}{a}\right)^{2k + 1} P_{2k + 1}(\cos\theta).
\]
Explicitly, the leading terms are
\[
\Phi(r, \theta) = V\left[\frac{3}{2}\frac{r}{a} P_1(\cos\theta) - \frac{7}{8}\left(\frac{r}{a}\right)^3 P_3(\cos\theta) + \frac{11}{16}\left(\frac{r}{a}\right)^5 P_5(\cos\theta) - \frac{75}{128}\left(\frac{r}{a}\right)^7 P_7(\cos\theta) + \cdots\right].
\]
Since $|P_{2k+1}(\cos\theta)| \le 1$ and $|c_{2k+1}| \sim \dfrac{2V}{\sqrt{\pi k}}$ by Stirling's approximation, the factor $(r/a)^{2k+1}$ ensures absolute and uniform convergence (together with all derivatives) on every closed ball $r \le r_0 < a$, so $\Phi$ is smooth and harmonic throughout $r < a$. On the equatorial plane $\theta = \pi/2$, $P_{2k+1}(0) = 0$ so $\Phi(r, \pi/2) = 0$ (in particular $\Phi = 0$ at the centre $r = 0$, as required by antisymmetry and the mean value property). Along the polar axis ($\theta = 0$, where $P_{2k+1}(1) = 1$), the series sums in closed form via $\sum_{m=0}^\infty P_{2m}(0) t^{2m} = (1 + t^2)^{-1/2}$ ($t = r/a$) to
\[
\Phi(r, 0) = V\left[\frac{a}{r} - \frac{a^2 - r^2}{r\sqrt{a^2 + r^2}}\right] \xrightarrow{r \to a^-} +V,
\]
providing an explicit confirmation of the boundary limit.

## Final Answers

1. For $\mathbf{r} \neq \mathbf{r}'$ in $\mathbb{R}^3$,
   \[
   \nabla^2\left(\frac{1}{|\mathbf{r} - \mathbf{r}'|}\right) = -\frac{3}{|\mathbf{r} - \mathbf{r}'|^3} + \frac{3|\mathbf{r} - \mathbf{r}'|}{|\mathbf{r} - \mathbf{r}'|^4} = 0.
   \]
2. Inside a sphere of radius $a$ bounded by hemispheres at $+V$ ($0 \le \theta < \pi/2$) and $-V$ ($\pi/2 < \theta \le \pi$),
   \[
   \Phi(r, \theta) = V \sum_{k = 0}^{\infty} (-1)^k \frac{(4k + 3)(2k)!}{2^{2k + 1} k! (k + 1)!} \left(\frac{r}{a}\right)^{2k + 1} P_{2k + 1}(\cos\theta)
   \]
   \[
   = V\left[\frac{3}{2}\frac{r}{a} P_1(\cos\theta) - \frac{7}{8}\left(\frac{r}{a}\right)^3 P_3(\cos\theta) + \frac{11}{16}\left(\frac{r}{a}\right)^5 P_5(\cos\theta) - \cdots\right].
   \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-2_qx03_check.py` (executed with `python3`, exit code 0).

```text
IB_methods_sheet-2_qx03_check.py: ALL CHECKS PASSED
  Part 1 max 4th-order FD |nabla^2 (1/|r-r'|)|:       9.869e-11
  Part 2 exact rational orthogonality & c_l (l<=15):  EXACT
  Part 2 first 4 odd coefficients c_{2k+1}/V:        3/2, -7/8, 11/16, -75/128
  Part 2 max spherical-polar Laplace residual:        1.719e-09
  Part 2 max z-axis series-vs-closed-form error:      1.377e-14
  Part 2 max off-axis series-vs-3D-Poisson error:     2.927e-11
  Part 2 near-boundary r=0.999a (theta=0.6, 2.5):     3.69804, -3.69797 (target +/-3.7)
```

## Confidence & Edge-Case Notes

- **Confidence:** High.
- **Phase B watchpoints:**
  1. **Sphere radius symbol:** The printed question does not name the sphere's radius (saying only "inside a spherical region bounded by two (conducting) hemispheres"). We introduced radius $a > 0$; the stored solution might use $a$, $R$, or take unit radius $a = 1$ (check whether it states its choice clearly!).
  2. **Sign and factorial factors in $c_{2k+1}$:** Very common slips when applying Rodrigues' formula are dropping the minus sign from the lower limit $x = 0$, getting the wrong power of $2$ ($2^{2k+1}$ vs $2^{2k}$), or getting $k!(k+1)!$ vs $(k!)^2$. Verify the first four terms: $\frac{3}{2}, -\frac{7}{8}, \frac{11}{16}, -\frac{75}{128}$ (times $V$). Also check equivalent forms such as $V\bigl(P_{2k}(0) - P_{2k+2}(0)\bigr)$ or $V(-1)^k(4k+3)\frac{(2k-1)!!}{(2k+2)!!}$ or $\frac{(4k+3)}{2(k+1)}V P_{2k}(0)$.
  3. **Sign convention for hemispheres:** Swapping which hemisphere is $+V$ and which is $-V$ flips the overall sign of $\Phi$; verify consistency with the chosen coordinate convention.
  4. **Regularity justifications:** Check that both $B_\ell = 0$ (regularity at $r = 0$) and the restriction to Legendre polynomials $P_\ell$ with $\lambda = \ell(\ell+1)$ (regularity at $\theta = 0, \pi$) are justified, and that $c_0 = 0$ and all even modes $c_{2k} = 0$ are clearly established.
