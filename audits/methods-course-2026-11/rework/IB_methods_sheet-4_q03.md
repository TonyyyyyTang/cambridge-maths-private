---
id: IB/methods/sheet-4/q03
slug: IB_methods_sheet-4_q03
source_archive: private/raw/sheets/IB/methods/2025-26/sheet-4.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against the 200 DPI rendered page `IB_methods_sheet-4_q03_p1.png` (`private/raw/sheets/IB/methods/2025-26/sheet-4.pdf`, page 1, Question 3). Question 3 sits near the upper third of page 1 and is immediately followed by Question 4 on the same page, so there is no spillover onto page 2.

- Equation `$u_{xx}+yu_{yy}+\frac12 u_y=0$` matches printed $u_{xx}+yu_{yy}+\tfrac{1}{2}u_y=0$.
- Canonical form `$U_{\eta\xi}=0$` matches printed $U_{\eta\xi}=0$ (including capital $U$ and subscript order $\eta\xi$).
- Target general solution `$u(x,y)=f(x+2\sqrt{-y})+g(x-2\sqrt{-y})$` and surrounding wording match verbatim.

Verdict: **`transcription: ok`** (exact match).

## Full Solution

### 1. Identification of the hyperbolic region

Consider the second-order linear PDE in $(x,y) \in \mathbb{R}^2$,
\[
u_{xx} + y u_{yy} + \frac{1}{2} u_y = 0.
\]
Comparing the principal part with the standard form $a(x,y) u_{xx} + b(x,y) u_{xy} + c(x,y) u_{yy}$, we have
\[
a(x,y) = 1, \qquad b(x,y) = 0, \qquad c(x,y) = y.
\]
The discriminant is
\[
\Delta(x,y) = b^2 - 4ac = -4y.
\]
Hence the PDE is:
- **elliptic** for $y > 0$ ($\Delta < 0$),
- **parabolic** on the line $y = 0$ ($\Delta = 0$),
- **hyperbolic** for $y < 0$ ($\Delta > 0$).

We therefore work on the open lower half-plane $\Omega_- = \{(x,y) \in \mathbb{R}^2 : y < 0\}$.

### 2. Characteristic curves and characteristic coordinates

A curve $\phi(x,y) = \text{const}$ is characteristic if $\phi$ satisfies the characteristic equation
\[
a \phi_x^2 + b \phi_x \phi_y + c \phi_y^2 = \phi_x^2 + y \phi_y^2 = 0.
\]
Equivalently, along a characteristic curve $y = y(x)$ (where $\mathrm{d}y/\mathrm{d}x = -\phi_x/\phi_y$),
\[
a\left(\frac{\mathrm{d}y}{\mathrm{d}x}\right)^2 - b\frac{\mathrm{d}y}{\mathrm{d}x} + c = \left(\frac{\mathrm{d}y}{\mathrm{d}x}\right)^2 + y = 0.
\]
In the hyperbolic region $y < 0$, we have $-y > 0$ and thus two distinct real families of characteristic slopes:
\[
\frac{\mathrm{d}y}{\mathrm{d}x} = \pm\sqrt{-y}.
\]
Separating variables (valid since $\sqrt{-y} > 0$ on $\Omega_-$) and letting $s = -y > 0$ (so $\mathrm{d}y = -\mathrm{d}s$),
\[
\int \frac{\mathrm{d}y}{\sqrt{-y}} = -\int s^{-1/2}\,\mathrm{d}s = -2\sqrt{-y} + \text{const} = \pm x + \text{const},
\]
so the two families of characteristic curves are
\[
x + 2\sqrt{-y} = \text{const} \qquad \text{and} \qquad x - 2\sqrt{-y} = \text{const}.
\]

Accordingly, introduce the characteristic coordinates on $\Omega_- = \{(x,y) : y < 0\}$:
\[
\xi(x,y) = x + 2\sqrt{-y}, \qquad \eta(x,y) = x - 2\sqrt{-y}.
\]
The first partial derivatives are
\[
\xi_x = 1, \qquad \xi_y = -\frac{1}{\sqrt{-y}}, \qquad \eta_x = 1, \qquad \eta_y = \frac{1}{\sqrt{-y}}.
\]
Let us verify non-degeneracy and global invertibility of the coordinate map $(x,y) \mapsto (\xi,\eta)$:
- The Jacobian determinant is
  \[
  J = \frac{\partial(\xi,\eta)}{\partial(x,y)} = \begin{vmatrix} 1 & -(-y)^{-1/2} \\ 1 & (-y)^{-1/2} \end{vmatrix} = \frac{2}{\sqrt{-y}} > 0 \qquad \text{for all } y < 0.
  \]
- Moreover, $\xi - \eta = 4\sqrt{-y} > 0$, and the map is a smooth bijection (diffeomorphism) from $\Omega_- = \{(x,y) : y < 0\}$ onto the open half-plane
  \[
  D = \{(\xi,\eta) \in \mathbb{R}^2 : \xi > \eta\},
  \]
  with explicit smooth inverse
  \[
  x = \frac{\xi + \eta}{2}, \qquad y = -\frac{(\xi - \eta)^2}{16}.
  \]

### 3. Transformation to canonical form

Assume $u \in C^2(\Omega_-)$, and define $U \in C^2(D)$ by $U(\xi,\eta) = u(x(\xi,\eta), y(\xi,\eta))$, so that $u(x,y) = U(\xi(x,y),\eta(x,y))$ and $U_{\xi\eta} = U_{\eta\xi}$ by Clairaut's theorem.

First compute the second partial derivatives of $\xi$ and $\eta$ for $y < 0$:
\[
\xi_{xx} = 0, \qquad \xi_{yy} = \frac{\partial}{\partial y}\!\left(-(-y)^{-1/2}\right) = -\frac{1}{2}(-y)^{-3/2},
\]
\[
\eta_{xx} = 0, \qquad \eta_{yy} = \frac{\partial}{\partial y}\!\left((-y)^{-1/2}\right) = \frac{1}{2}(-y)^{-3/2}.
\]
Applying the chain rule to $u(x,y) = U(\xi,\eta)$:
\[
u_x = U_\xi \xi_x + U_\eta \eta_x = U_\xi + U_\eta,
\]
\[
u_{xx} = U_{\xi\xi} + 2 U_{\eta\xi} + U_{\eta\eta},
\]
\[
u_y = U_\xi \xi_y + U_\eta \eta_y = \frac{-U_\xi + U_\eta}{\sqrt{-y}},
\]
and
\[
\begin{aligned}
u_{yy} &= U_{\xi\xi}\xi_y^2 + 2 U_{\eta\xi}\xi_y\eta_y + U_{\eta\eta}\eta_y^2 + U_\xi \xi_{yy} + U_\eta \eta_{yy} \\
&= \frac{1}{-y} U_{\xi\xi} - \frac{2}{-y} U_{\eta\xi} + \frac{1}{-y} U_{\eta\eta} - \frac{1}{2(-y)^{3/2}} U_\xi + \frac{1}{2(-y)^{3/2}} U_\eta.
\end{aligned}
\]
Multiplying $u_{yy}$ by $y = -(-y)$ yields
\[
y u_{yy} = -U_{\xi\xi} + 2 U_{\eta\xi} - U_{\eta\eta} + \frac{1}{2\sqrt{-y}} U_\xi - \frac{1}{2\sqrt{-y}} U_\eta.
\]
Also,
\[
\frac{1}{2} u_y = -\frac{1}{2\sqrt{-y}} U_\xi + \frac{1}{2\sqrt{-y}} U_\eta.
\]
Summing the three contributions gives the exact operator identity
\[
\begin{aligned}
u_{xx} + y u_{yy} + \frac{1}{2} u_y
&= (1 - 1) U_{\xi\xi} + (2 + 2) U_{\eta\xi} + (1 - 1) U_{\eta\eta} \\
&\quad + \left(\frac{1}{2\sqrt{-y}} - \frac{1}{2\sqrt{-y}}\right) U_\xi + \left(-\frac{1}{2\sqrt{-y}} + \frac{1}{2\sqrt{-y}}\right) U_\eta \\
&= 4 U_{\eta\xi}.
\end{aligned}
\]
Thus the PDE reduces in the hyperbolic region $y < 0$ to the canonical form
\[
4 U_{\eta\xi} = 0 \iff U_{\eta\xi} = 0 \qquad \text{on } D = \{(\xi,\eta) \in \mathbb{R}^2 : \xi > \eta\}.
\]

### 4. Deduction of the general solution

We solve $\frac{\partial}{\partial \xi}(U_\eta) = 0$ on $D = \{(\xi,\eta) : \xi > \eta\}$. Notice that:
- For every fixed $\eta \in \mathbb{R}$, the horizontal slice $\{\xi \in \mathbb{R} : (\xi,\eta) \in D\} = (\eta,\infty)$ is a connected open interval. Hence $\frac{\partial}{\partial\xi}(U_\eta) = 0$ implies that $U_\eta(\xi,\eta)$ is independent of $\xi$, i.e.
  \[
  U_\eta(\xi,\eta) = h(\eta)
  \]
  for some $C^1$ function $h : \mathbb{R} \to \mathbb{R}$.
- Let $g : \mathbb{R} \to \mathbb{R}$ be an antiderivative of $h$, so $g'(\eta) = h(\eta)$ (with $g \in C^2(\mathbb{R})$). Then on $D$,
  \[
  \frac{\partial}{\partial \eta}\bigl(U(\xi,\eta) - g(\eta)\bigr) = U_\eta(\xi,\eta) - h(\eta) = 0.
  \]
- For every fixed $\xi \in \mathbb{R}$, the vertical slice $\{\eta \in \mathbb{R} : (\xi,\eta) \in D\} = (-\infty,\xi)$ is also a connected open interval. Hence $U(\xi,\eta) - g(\eta)$ is independent of $\eta$, so
  \[
  U(\xi,\eta) = f(\xi) + g(\eta)
  \]
  for some $C^2$ function $f : \mathbb{R} \to \mathbb{R}$.

Conversely, for any $C^2$ functions $f, g : \mathbb{R} \to \mathbb{R}$, $U(\xi,\eta) = f(\xi) + g(\eta)$ satisfies $U_{\eta\xi} = 0$ on $D$. Transforming back to $(x,y) \in \Omega_-$ gives the general solution in the hyperbolic region $y < 0$:
\[
u(x,y) = f\!\left(x + 2\sqrt{-y}\right) + g\!\left(x - 2\sqrt{-y}\right),
\]
for arbitrary (twice-differentiable, for classical $C^2$ solutions) functions $f, g : \mathbb{R} \to \mathbb{R}$ (unique up to shifting $f \mapsto f + c$, $g \mapsto g - c$ by a constant $c$).

## Final Answers

- **Hyperbolic region**: $y < 0$ (discriminant $\Delta = -4y > 0$).
- **Characteristic coordinates**: $\xi = x + 2\sqrt{-y}$, $\eta = x - 2\sqrt{-y}$ (diffeomorphism from $\{y < 0\}$ onto $\{\xi > \eta\}$ with Jacobian $J = 2/\sqrt{-y} > 0$).
- **Transformed PDE / canonical form**: $u_{xx} + y u_{yy} + \frac{1}{2}u_y = 4 U_{\eta\xi} = 0 \iff U_{\eta\xi} = 0$.
- **General solution in $y < 0$**: $u(x,y) = f(x + 2\sqrt{-y}) + g(x - 2\sqrt{-y})$ for arbitrary ($C^2$) functions $f, g$.

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-4_q03_check.py`

```text
IB_methods_sheet-4_q03 verification PASSED:
  max inverse map error      = 8.882e-16
  max characteristic eq err  = 2.237e-11
  max Jacobian error         = 5.542e-11
  max operator identity err  = 1.325e-07 (rel 2.039e-09)  [u_xx + y*u_yy + 0.5*u_y == 4*U_{eta,xi}]
  max general sol PDE resid  = 2.159e-08
```

## Confidence & Edge-Case Notes

- **Confidence**: high.
- **Phase B watchpoints**:
  1. Check that the existing solution explicitly determines the hyperbolic region ($y < 0$) and verifies the Jacobian $J = \frac{\partial(\xi,\eta)}{\partial(x,y)} = 2/\sqrt{-y} \neq 0$ (invertibility of the coordinate change).
  2. Watch the sign of $\xi_{yy} = -\frac{1}{2}(-y)^{-3/2}$ and $\eta_{yy} = +\frac{1}{2}(-y)^{-3/2}$: a sign mistake here causes $y\xi_{yy} + \frac{1}{2}\xi_y$ to fail to cancel (in fact $y\xi_{yy} = +\frac{1}{2}(-y)^{-1/2}$ cancels $\frac{1}{2}\xi_y = -\frac{1}{2}(-y)^{-1/2}$).
  3. Watch that the image domain in $(\xi,\eta)$-space is the half-plane $\xi > \eta$ (whose horizontal and vertical slices are connected intervals), which justifies integrating $U_{\eta\xi} = 0$ globally to $f(\xi) + g(\eta)$.
  4. Watch the factor $4$ in $u_{xx} + y u_{yy} + \frac{1}{2}u_y = 4 U_{\eta\xi}$ (comes from $2 U_{\eta\xi}$ in $u_{xx}$ plus $+2 U_{\eta\xi}$ in $y u_{yy}$, since $y \cdot (-2/(-y)) = +2$).
