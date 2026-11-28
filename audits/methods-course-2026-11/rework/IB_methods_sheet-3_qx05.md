---
id: IB/methods/sheet-3/qx05
slug: IB_methods_sheet-3_qx05
source_archive: private/raw/sheets/IB/methods/2001-02/sheet-3.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-3_qx05_p1.png` (Michaelmas term 2001, Dr E.P.S. Shellard, *Methods — Examples III*, Question 5, under *Extrema and Lagrange multipliers*). Question 5 finishes cleanly in the middle of page 1.

- **Status:** `ok` (exact match for every word, symbol, sign, exponent, coefficient, and punctuation mark).
- **Purely typographic note:** The printed heading *Stationary points.* is both italicised and underlined (consistent with the sheet's title style for all questions), rendered in `question.tex` as `\emph{Stationary points.}`. The symbol `$R^3$` faithfully matches the printed italic $R^3$.

## Full Solution

We analyse the smooth function $\phi : \mathbb{R}^3 \to \mathbb{R}$ defined by
\[
\phi(x_1, x_2, x_3) = \frac{1}{4}(x_1^4 + x_2^4 + x_3^4) - x_2 x_3 - x_3 x_1 - x_1 x_2.
\]

### 1. All stationary points of $\phi$

The condition $\nabla \phi(x_1, x_2, x_3) = \mathbf{0}$ yields the system
\[
\frac{\partial \phi}{\partial x_1} = x_1^3 - x_2 - x_3 = 0, \qquad
\frac{\partial \phi}{\partial x_2} = x_2^3 - x_3 - x_1 = 0, \qquad
\frac{\partial \phi}{\partial x_3} = x_3^3 - x_1 - x_2 = 0.
\]
Setting $S = x_1 + x_2 + x_3$ and adding $x_i$ to the $i$-th equation gives
\[
x_1^3 + x_1 = S, \qquad x_2^3 + x_2 = S, \qquad x_3^3 + x_3 = S.
\]
Consider $g(t) = t^3 + t$ for $t \in \mathbb{R}$. Since $g'(t) = 3t^2 + 1 \ge 1 > 0$ for all $t \in \mathbb{R}$, $g$ is strictly increasing and therefore injective on $\mathbb{R}$. (Equivalently, subtracting any two equations gives $(x_i - x_j)(x_i^2 + x_i x_j + x_j^2 + 1) = 0$, and $x_i^2 + x_i x_j + x_j^2 + 1 = (x_i + \tfrac{1}{2}x_j)^2 + \tfrac{3}{4}x_j^2 + 1 \ge 1 > 0$, forcing $x_i = x_j$.)

Hence every stationary point satisfies $x_1 = x_2 = x_3 = t$ for some $t \in \mathbb{R}$. Substituting back into $x_1^3 - x_2 - x_3 = 0$ gives
\[
t^3 - 2t = t(t^2 - 2) = 0 \implies t \in \{0, \sqrt{2}, -\sqrt{2}\}.
\]
Thus $\phi$ has **exactly three** stationary points in $\mathbb{R}^3$:
\[
\mathbf{x}_+ = (\sqrt{2}, \sqrt{2}, \sqrt{2}), \qquad
\mathbf{x}_- = (-\sqrt{2}, -\sqrt{2}, -\sqrt{2}), \qquad
\mathbf{x}_0 = (0, 0, 0).
\]

### 2. Global minimum value and number of minimisers

We prove directly via a sum-of-squares identity that $\phi(\mathbf{x}) \ge -3$ for all $\mathbf{x} \in \mathbb{R}^3$, with equality precisely at $\mathbf{x}_+$ and $\mathbf{x}_-$. Observe that
\[
\frac{1}{2}\bigl[(x_1 - x_2)^2 + (x_2 - x_3)^2 + (x_3 - x_1)^2\bigr]
= x_1^2 + x_2^2 + x_3^2 - x_2 x_3 - x_3 x_1 - x_1 x_2,
\]
and
\[
\frac{1}{4}\sum_{i=1}^3 (x_i^2 - 2)^2
= \frac{1}{4}(x_1^4 + x_2^4 + x_3^4) - (x_1^2 + x_2^2 + x_3^2) + 3.
\]
Adding these two identities cancels $x_1^2 + x_2^2 + x_3^2$ and yields
\[
\phi(x_1, x_2, x_3) = -3 + \frac{1}{4}\sum_{i=1}^3 (x_i^2 - 2)^2 + \frac{1}{2}\bigl[(x_1 - x_2)^2 + (x_2 - x_3)^2 + (x_3 - x_1)^2\bigr].
\]
Every term after $-3$ is non-negative, so $\phi(x_1, x_2, x_3) \ge -3$ everywhere in $\mathbb{R}^3$. Equality $\phi = -3$ holds if and only if every square vanishes simultaneously:
\[
x_1 = x_2 = x_3 \quad \text{and} \quad x_1^2 = x_2^2 = x_3^2 = 2,
\]
which holds at **exactly $2$ points**: $\mathbf{x}_+ = (\sqrt{2}, \sqrt{2}, \sqrt{2})$ and $\mathbf{x}_- = (-\sqrt{2}, -\sqrt{2}, -\sqrt{2})$.

*(Alternative classification via the Hessian and coercivity.)* The Hessian matrix is
\[
H(\mathbf{x}) = \begin{pmatrix} 3x_1^2 & -1 & -1 \\ -1 & 3x_2^2 & -1 \\ -1 & -1 & 3x_3^2 \end{pmatrix}.
\]
At $\mathbf{x}_\pm$, $3x_i^2 = 6$, so $H(\mathbf{x}_\pm) = 7I - J$, where $J$ is the $3\times 3$ all-ones matrix (which has eigenvalues $3$ on $\operatorname{span}\{(1,1,1)^T\}$ and $0$ on the orthogonal plane $x_1 + x_2 + x_3 = 0$). Hence $H(\mathbf{x}_\pm)$ has eigenvalues $4, 7, 7 > 0$, confirming that $\mathbf{x}_\pm$ are strict local minima with $\phi(\mathbf{x}_\pm) = \frac{3}{4}(4) - 3(2) = -3$. Moreover, since $x_2 x_3 + x_3 x_1 + x_1 x_2 \le x_1^2 + x_2^2 + x_3^2$, we have $\phi(\mathbf{x}) \ge \frac{1}{4}\sum_{i=1}^3 x_i^4 - \sum_{i=1}^3 x_i^2 \to +\infty$ as $\|\mathbf{x}\| \to \infty$, so $\phi$ attains a global minimum on $\mathbb{R}^3$, which must occur at a stationary point; comparing $\phi(\mathbf{x}_\pm) = -3 < 0 = \phi(\mathbf{x}_0)$ again shows that the least value is $-3$, attained at exactly $2$ points.

### 3. The unique saddle point and the tangent double cone to $\phi = 0$

The only remaining stationary point is $\mathbf{x}_0 = (0,0,0)$, where $\phi(\mathbf{0}) = 0$ and the Hessian is
\[
H(\mathbf{0}) = \begin{pmatrix} 0 & -1 & -1 \\ -1 & 0 & -1 \\ -1 & -1 & 0 \end{pmatrix} = I - J.
\]
Its eigenvalues and an orthonormal eigenbasis are:
- $\lambda_1 = 1 - 3 = -2$, with unit eigenvector $\mathbf{e}_1 = \frac{1}{\sqrt{3}}(1, 1, 1)^T$;
- $\lambda_2 = \lambda_3 = 1 - 0 = +1$ (multiplicity $2$), spanning the plane $x_1 + x_2 + x_3 = 0$, with orthonormal eigenvectors e.g. $\mathbf{e}_2 = \frac{1}{\sqrt{2}}(1, -1, 0)^T$ and $\mathbf{e}_3 = \frac{1}{\sqrt{6}}(1, 1, -2)^T$.

Because $H(\mathbf{0})$ has both positive and negative eigenvalues ($-2, +1, +1$), it is indefinite, so $\mathbf{x}_0 = (0,0,0)$ is a **saddle point**. Since $\phi$ has only three stationary points in total and the other two are strict minima, $\phi$ has **exactly one** saddle point.

To determine the local geometry of the surface $\phi(\mathbf{x}) = 0$ near $\mathbf{x} = \mathbf{0}$, decompose $\phi$ into its quadratic and quartic parts:
\[
\phi(\mathbf{x}) = Q(\mathbf{x}) + \frac{1}{4}(x_1^4 + x_2^4 + x_3^4), \qquad Q(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T H(\mathbf{0}) \mathbf{x} = -x_2 x_3 - x_3 x_1 - x_1 x_2.
\]
Writing $\mathbf{x} = u_1 \mathbf{e}_1 + u_2 \mathbf{e}_2 + u_3 \mathbf{e}_3$ in the orthonormal eigenbasis, the quadratic form diagonalises as
\[
Q(\mathbf{x}) = \frac{1}{2}\bigl(\lambda_1 u_1^2 + \lambda_2 u_2^2 + \lambda_3 u_3^2\bigr) = -u_1^2 + \frac{1}{2}(u_2^2 + u_3^2).
\]
Equivalently, introducing spherical polar coordinates $(r, \theta, \psi)$ about the axis $\mathbf{e}_1 = \frac{1}{\sqrt{3}}(1,1,1)^T$ so that $u_1 = r\cos\theta$ and $\sqrt{u_2^2 + u_3^2} = r\sin\theta$ ($\theta \in [0, \pi]$ is the angle between $\mathbf{x}$ and $\mathbf{e}_1$), we have
\[
Q(\mathbf{x}) = \frac{1}{2}r^2\bigl(\sin^2\theta - 2\cos^2\theta\bigr).
\]
The zero set $Q(\mathbf{x}) = 0$ is therefore the circular **double cone**
\[
u_2^2 + u_3^2 = 2 u_1^2 \iff \tan^2\theta = 2 \iff \theta = \alpha \text{ or } \theta = \pi - \alpha,
\]
with apex at the origin, axis along the line $x_1 = x_2 = x_3$, and semi-vertical angle (semi-angle)
\[
\alpha = \tan^{-1}(\sqrt{2}).
\]
(As a geometric check, each coordinate axis lies on this cone since $Q = 0$ on the axes and $\cos\theta = \mathbf{e}_1 \cdot (1,0,0)^T = 1/\sqrt{3} \implies \tan\theta = \sqrt{2}$.)

Finally, we verify rigorously that the surface $\phi(\mathbf{x}) = 0$ is tangent to this double cone at the origin:
- For any non-zero point $\mathbf{x} = r\hat{\mathbf{n}}$ ($r > 0$, $\|\hat{\mathbf{n}}\| = 1$) on $\phi(\mathbf{x}) = 0$, dividing $\phi(r\hat{\mathbf{n}}) = 0$ by $r^2 > 0$ gives
  \[
  Q(\hat{\mathbf{n}}) = \frac{1}{2}(\sin^2\theta - 2\cos^2\theta) = -\frac{1}{4}r^2(\hat{n}_1^4 + \hat{n}_2^4 + \hat{n}_3^4).
  \]
  Hence as $r \to 0^+$ along $\phi = 0$, $\sin^2\theta - 2\cos^2\theta = O(r^2) \to 0$, so every limiting tangent ray lies on the double cone $\tan^2\theta = 2$ (approaching from the interior side $\tan\theta < \sqrt{2}$ where $Q < 0$ balances the positive quartic term).
- Conversely, for every azimuthal angle $\psi$ and every $\theta$ close to $\alpha$ with $\theta < \alpha$ (or close to $\pi - \alpha$ with $\theta > \pi - \alpha$), we have $Q(\hat{\mathbf{n}}) < 0$ while $\sum \hat{n}_i^4 \ge \frac{1}{3} > 0$, giving a unique $r(\theta, \psi) = \sqrt{-4Q(\hat{\mathbf{n}})/\sum\hat{n}_i^4} > 0$ on $\phi = 0$ with $r \to 0^+$ as $\theta \to \alpha^-$ or $(\pi - \alpha)^+$. Thus both nappes of the double cone are realised as the exact tangent cone to $\phi = 0$ at $\mathbf{0}$, and the surface normal $\nabla\phi(\mathbf{x}) = H(\mathbf{0})\mathbf{x} + O(r^3)$ converges in direction to the cone normal $H(\mathbf{0})\mathbf{x}$ (which has norm $\sqrt{2}\,r + O(r^3)$) with angular error $O(r^2)$ as $r \to 0^+$.

## Final Answers

- **Number of points where $\phi$ takes its minimum value:** $2$, namely $(\sqrt{2}, \sqrt{2}, \sqrt{2})$ and $(-\sqrt{2}, -\sqrt{2}, -\sqrt{2})$.
- **Least value of $\phi$:** $-3$.
- **Saddle point:** Unique saddle point at $(0, 0, 0)$, where $H(\mathbf{0})$ has eigenvalues $-2, +1, +1$.
- **Tangent cone to $\phi = 0$ at $(0,0,0)$:** Circular double cone $u_2^2 + u_3^2 = 2u_1^2$ (i.e. $x_2 x_3 + x_3 x_1 + x_1 x_2 = 0$) with axis along $x_1 = x_2 = x_3$ and semi-angle $\tan^{-1}(\sqrt{2})$.

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-3_qx05_check.py`

```text
[PASS] Sum-of-squares identity verified on 5000 random points (max err = 1.14e-13)
[PASS] Stationary points, values, gradients, and Hessians verified
[PASS] No extraneous stationary points on [-2,2]^3 grid (min |grad| outside 0.15-balls = 0.1655)
[PASS] Hessian eigenpairs verified: minima eigenvalues (4, 7, 7), saddle eigenvalues (-2, 1, 1)
[PASS] Surface phi=0 approaches both nappes of cone (semi-angle atan(sqrt(2))=0.955317 rad): max |phi|=6.94e-18, max normal angle/r^2=0.7071
ALL CHECKS PASSED for IB_methods_sheet-3_qx05.
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  1. Ensure the solution actually proves **all** stationary points satisfy $x_1 = x_2 = x_3$ (via injectivity of $t^3 + t$ or factoring $(x_i - x_j)(x_i^2 + x_i x_j + x_j^2 + 1) = 0$) rather than merely guessing symmetric stationary points $x_1 = x_2 = x_3$.
  2. Ensure the solution proves $-3$ is the **global** minimum (either via the exact sum-of-squares decomposition or via coercivity $\phi(\mathbf{x}) \to +\infty$ as $\|\mathbf{x}\| \to \infty$), not just that $\pm(\sqrt{2},\sqrt{2},\sqrt{2})$ are local minima.
  3. Watch out for the factor of $\frac{1}{2}$ in the quadratic form $Q = -u_1^2 + \frac{1}{2}(u_2^2 + u_3^2)$, which gives $\sqrt{u_2^2 + u_3^2}/|u_1| = \sqrt{2}$, i.e., semi-angle $\tan^{-1}(\sqrt{2})$ (not $\tan^{-1}(1/\sqrt{2})$—it is easy to invert $\tan\theta$ if confusing the axis $\mathbf{e}_1$ with the transverse directions!).
