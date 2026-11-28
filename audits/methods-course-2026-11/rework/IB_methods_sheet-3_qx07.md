---
id: IB/methods/sheet-3/qx07
slug: IB_methods_sheet-3_qx07
source_archive: private/raw/sheets/IB/methods/2001-02/sheet-3.pdf
pages_checked: [1, 2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-3_qx07_p1.png` (Michaelmas term 2001, Dr E.P.S. Shellard, *Methods — Examples III*, Question 7, under *Calculus of variations and the Euler-Lagrange equations*). Also inspected `IB_methods_sheet-3_qx07_p2.png` because Question 7 sits at the bottom of page 1; page 2 begins directly with Question 8, confirming nothing spills over.

- **Status:** `ok` (exact match for every word, symbol, italic $l$, exponent `^{\frac12}`, punctuation, and bracketed note).
- **Purely typographic note:** The printed heading *Brachistochrone.* is both italicised and underlined, matching the sheet's title convention, encoded as `\emph{Brachistochrone.}`.

## Full Solution

### 1. Formulation of the travel-time functional

Choose Cartesian coordinates $(x, y, z)$ with departure point $A = (0, 0, 0)$ and destination $B = (l, 0, 0)$ ($l > 0$), with $x$ horizontal from $A$ to $B$, $y$ measured **vertically downwards** (so underground tunnels satisfy $y \ge 0$), and $z$ horizontal perpendicular to $AB$.

Because the tunnel is frictionless, gravity $g > 0$ is uniform, and the carriage is released from rest at $y = 0$, conservation of mechanical energy gives
\[
\frac{1}{2} m v^2 - m g y = 0 \implies v = \sqrt{2 g y}.
\]
In particular, $y \ge 0$ along the entire tunnel (and $y > 0$ on the interior so the carriage moves and reaches $B$ in finite time). The element of arc length is $ds = \sqrt{dx^2 + dy^2 + dz^2} \ge \sqrt{dx^2 + dy^2}$, with equality if and only if $dz \equiv 0$. Since $v = \sqrt{2gy}$ is independent of $z$, projecting any 3D curve onto the vertical plane $z = 0$ strictly reduces the travel time unless $z \equiv 0$. Hence any minimising tunnel lies in the vertical plane $z = 0$.

For a smooth curve in the $(x, y)$-plane connecting $(0, 0)$ to $(l, 0)$, the travel time is
\[
T = \int_A^B \frac{ds}{v} = \int_0^l F(y, y')\,dx, \qquad F(y, y') = \frac{\sqrt{1 + y'^2}}{\sqrt{2gy}},
\]
where $y' = dy/dx$. (Equivalently, in parametric form $(x(u), y(u))$ avoiding any coordinate singularity where the tangent becomes vertical,
\[
T[x, y] = \frac{1}{\sqrt{2g}} \int_{u_0}^{u_1} L(y, \dot{x}, \dot{y})\,du, \qquad L(y, \dot{x}, \dot{y}) = \sqrt{\frac{\dot{x}^2 + \dot{y}^2}{y}},
\]
which does not even assume *a priori* that the curve is a single-valued graph $y = y(x)$.)

### 2. Euler–Lagrange equation and Beltrami first integral

Since $F(y, y')$ has no explicit dependence on $x$, any smooth stationary curve satisfies the Euler–Lagrange equation $\frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) - \frac{\partial F}{\partial y} = 0$, which admits the **Beltrami first integral**:
\[
\frac{d}{dx}\!\left(F - y'\frac{\partial F}{\partial y'}\right)
= \frac{\partial F}{\partial y}y' + \frac{\partial F}{\partial y'}y'' - y''\frac{\partial F}{\partial y'} - y'\frac{d}{dx}\!\left(\frac{\partial F}{\partial y'}\right)
= -y'\!\left[\frac{d}{dx}\!\left(\frac{\partial F}{\partial y'}\right) - \frac{\partial F}{\partial y}\right] = 0.
\]
Hence $F - y'\frac{\partial F}{\partial y'} = \text{const}$ along the extremal. Computing this directly:
\[
\frac{\partial F}{\partial y'} = \frac{y'}{\sqrt{2gy}\sqrt{1 + y'^2}}
\implies
F - y'\frac{\partial F}{\partial y'} = \frac{(1 + y'^2) - y'^2}{\sqrt{2gy}\sqrt{1 + y'^2}} = \frac{1}{\sqrt{2gy(1 + y'^2)}} = \text{const}.
\]
(Notice that the parametric Euler–Lagrange equation for the cyclic coordinate $x$ gives $\frac{\partial L}{\partial \dot{x}} = \frac{\dot{x}}{\sqrt{y(\dot{x}^2 + \dot{y}^2)}} = C = \text{const}$ directly, with no singularity at vertical tangents $\dot{x} = 0$. Since $x$ increases from $0$ to $l > 0$, we must have $C > 0$, forcing $\dot{x} > 0$ whenever $y > 0$, so every extremal automatically has $x$ strictly increasing and is indeed a graph $y = y(x)$!)

Write the positive constant in $y(1 + y'^2) = \text{const}$ as $2a$ with $a > 0$ (physically, $2a = y_{\max}$ is the maximum depth, attained where $y' = 0$). Then
\[
y(1 + y'^2) = 2a \implies \left(\frac{dy}{dx}\right)^2 = \frac{2a - y}{y}.
\]
(Note that differentiating $y(1 + y'^2) = 2a$ gives $y'(1 + y'^2 + 2y y'') = 0$; our cycloid below has $y' = 0$ only at the isolated lowest point $x = l/2$ and satisfies the full second-order Euler–Lagrange equation $2yy'' + 1 + y'^2 = 0$ everywhere on $(0, l)$.)

### 3. Parametric cycloid solution and boundary conditions

Set
\[
y = a(1 - \cos\theta) = 2a\sin^2\!\frac{\theta}{2},
\]
so that $y$ increases from $0$ to $2a$ as $\theta$ runs from $0$ to $\pi$, and decreases back from $2a$ to $0$ as $\theta$ runs from $\pi$ to $2\pi$. Let us check signs carefully on both halves:
- For $\theta \in (0, \pi)$ (descending from $A$ to the bottom): $dy/dx = +\sqrt{(2a - y)/y} = +\cot(\theta/2) > 0$.
- For $\theta \in (\pi, 2\pi)$ (ascending from the bottom to $B$): $dy/dx = -\sqrt{(2a - y)/y} = -|\cot(\theta/2)| = +\cot(\theta/2) < 0$, since $\theta/2 \in (\pi/2, \pi)$ where $\cot(\theta/2) < 0$.

Thus $dy/dx = \cot(\theta/2)$ holds smoothly across the entire interval $\theta \in (0, 2\pi)$, including at $\theta = \pi$ where $dy/dx = 0$. Using $dy/d\theta = a\sin\theta = 2a\sin(\theta/2)\cos(\theta/2)$, we obtain
\[
\frac{dx}{d\theta} = \frac{dx}{dy}\frac{dy}{d\theta} = \tan\frac{\theta}{2} \cdot 2a\sin\frac{\theta}{2}\cos\frac{\theta}{2} = 2a\sin^2\!\frac{\theta}{2} = a(1 - \cos\theta).
\]
Integrating and imposing $(x, y) = (0, 0)$ at $\theta = 0$ gives the inverted cycloid
\[
x(\theta) = a(\theta - \sin\theta), \qquad y(\theta) = a(1 - \cos\theta).
\]
To reach $B = (l, 0)$ with $y(\theta_B) = a(1 - \cos\theta_B) = 0$ and $x(\theta_B) = l > 0$, we must have $\theta_B = 2\pi N$ for some integer $N \ge 1$, giving $a = l/(2\pi N)$. As shown immediately below, $N$ arches take time $\sqrt{N}\,(2\pi l/g)^{1/2}$, which is strictly minimised by the single arch $N = 1$. Hence
\[
\theta_B = 2\pi, \qquad a = \frac{l}{2\pi}.
\]
The maximum depth is $y(\pi) = 2a = l/\pi$, and the slope is $dy/dx = \cot(\theta/2)$: as $\theta \to 0^+$ (at $A$), $dy/dx \to +\infty$, and as $\theta \to 2\pi^-$ (at $B$), $dy/dx \to -\infty$, confirming that the tunnel is **vertical at both $A$ and $B$**.

### 4. Minimum travel time

Along the cycloid for $\theta \in [0, 2\pi]$ (where $\sin(\theta/2) \ge 0$):
\[
ds = \sqrt{\left(\frac{dx}{d\theta}\right)^2 + \left(\frac{dy}{d\theta}\right)^2}\,d\theta
= a\sqrt{(1 - \cos\theta)^2 + \sin^2\theta}\,d\theta
= a\sqrt{2(1 - \cos\theta)}\,d\theta
= 2a\sin\frac{\theta}{2}\,d\theta,
\]
and
\[
v = \sqrt{2gy} = \sqrt{2ga(1 - \cos\theta)} = \sqrt{4ga\sin^2\!\frac{\theta}{2}} = 2\sqrt{ga}\sin\frac{\theta}{2}.
\]
The factors of $2\sin(\theta/2)$ cancel identically, yielding a completely regular integrand:
\[
dt = \frac{ds}{v} = \sqrt{\frac{a}{g}}\,d\theta.
\]
Integrating from $\theta = 0$ to $\theta = 2\pi N$ with $a = l/(2\pi N)$ gives
\[
T_N = \int_0^{2\pi N} \sqrt{\frac{a}{g}}\,d\theta = 2\pi N\sqrt{\frac{l}{2\pi N g}} = \sqrt{N}\left(\frac{2\pi l}{g}\right)^{\frac{1}{2}}.
\]
Taking the single arch $N = 1$ gives the travel time
\[
T = \left(\frac{2\pi l}{g}\right)^{\frac{1}{2}}.
\]

### 5. Proof of global minimality (Hamilton–Jacobi calibration)

The Euler–Lagrange equation only shows stationarity; to prove rigorously that $T = (2\pi l/g)^{1/2}$ is the **strict global minimum** among all piecewise-$C^1$ curves from $A$ to $B$ in $\{y \ge 0\}$, consider the family of cycloids emanating from $(0, 0)$ with parameters $a > 0$, $\theta \in (0, 2\pi)$. Since
\[
\frac{d}{d\theta}\!\left(\frac{\theta - \sin\theta}{1 - \cos\theta}\right)
= \frac{4\sin(\theta/2)\bigl(\sin(\theta/2) - (\theta/2)\cos(\theta/2)\bigr)}{(1 - \cos\theta)^2} > 0
\quad \text{on } (0, 2\pi)
\]
(as $p(u) = \sin u - u\cos u$ satisfies $p(0) = 0$ and $p'(u) = u\sin u > 0$ on $(0, \pi)$), the ratio $x/y = (\theta - \sin\theta)/(1 - \cos\theta)$ increases strictly from $0$ to $+\infty$ on $(0, 2\pi)$. Thus every point $(x, y)$ in the open quadrant $\{x > 0, y > 0\}$ lies on a unique cycloid $(a, \theta) \in (0, \infty) \times (0, 2\pi)$.

Define the travel-time field $S : \mathbb{R} \times [0, \infty) \to [0, \infty)$ by
\[
S(x, y) = \begin{cases}
\theta\sqrt{a/g}, & x > 0,\; y > 0, \\
\sqrt{2y/g}, & x \le 0,\; y \ge 0, \\
\sqrt{2\pi x/g}, & x > 0,\; y = 0.
\end{cases}
\]
Direct differentiation on $\{x > 0, y > 0\}$ using $dx = (\theta - \sin\theta)\,da + a(1-\cos\theta)\,d\theta$ and $dy = (1-\cos\theta)\,da + a\sin\theta\,d\theta$ gives
\[
dx + \cot\frac{\theta}{2}\,dy
= \bigl(\theta - \sin\theta + 2\sin\frac{\theta}{2}\cos\frac{\theta}{2}\bigr)\,da + 2a\bigl(\sin^2\frac{\theta}{2} + \cos^2\frac{\theta}{2}\bigr)\,d\theta
= \theta\,da + 2a\,d\theta,
\]
and dividing by $2\sqrt{ga}$ yields
\[
dS = d\!\left(\theta\sqrt{\frac{a}{g}}\right) = \frac{1}{2\sqrt{ga}}\,dx + \frac{\cot(\theta/2)}{2\sqrt{ga}}\,dy.
\]
Hence on $\{x > 0, y > 0\}$,
\[
|\nabla S|^2 = \frac{1 + \cot^2(\theta/2)}{4ga} = \frac{1}{4ga\sin^2(\theta/2)} = \frac{1}{2gy}.
\]
As $x \to 0^+$ with $y > 0$, $\theta \to 0^+$ and $a \sim 2y/\theta^2 \to \infty$, so $S \to \sqrt{2y/g}$ and $\nabla S \to (0, 1/\sqrt{2gy})$, matching the formula on $\{x \le 0, y > 0\}$; thus $S$ is globally $C^1$ on $\{y > 0\}$ and continuous on $\{y \ge 0\}$ (with $S(0, 0) = 0$ and $S(l, 0) = (2\pi l/g)^{1/2}$).

By Cauchy–Schwarz, for **any** piecewise-$C^1$ curve $\gamma(u) = (x(u), y(u))$ from $(0,0)$ to $(l,0)$ with $y > 0$ on $(0, 1)$,
\[
\frac{\sqrt{\dot{x}^2 + \dot{y}^2}}{\sqrt{2gy}} = |\nabla S|\,|\dot{\gamma}| \ge \nabla S \cdot \dot{\gamma} = \frac{d}{du} S(\gamma(u)).
\]
Integrating from $u = 0$ to $u = 1$ gives
\[
T[\gamma] \ge S(l, 0) - S(0, 0) = \left(\frac{2\pi l}{g}\right)^{\frac{1}{2}},
\]
with equality if and only if $\dot{\gamma}$ is a positive multiple of $\nabla S$ everywhere, i.e., $\gamma$ is the unique single-arch cycloid with $a = l/(2\pi)$ in the vertical plane $z = 0$.

## Final Answers

- **Optimal tunnel curve (in the vertical plane $z = 0$, $y$ downwards):** the single-arch inverted cycloid
  \[
  x(\theta) = \frac{l}{2\pi}(\theta - \sin\theta), \qquad y(\theta) = \frac{l}{2\pi}(1 - \cos\theta), \qquad \theta \in [0, 2\pi],
  \]
  which reaches maximum depth $l/\pi$ at $x = l/2$ and has vertical tangents ($dy/dx \to \pm\infty$) at $A$ ($\theta = 0$) and $B$ ($\theta = 2\pi$).
- **Minimum travel time:** $T_{\min} = \left(\dfrac{2\pi l}{g}\right)^{\frac{1}{2}}$.

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-3_qx07_check.py`

```text
[PASS] Cycloid satisfies Beltrami (max err=6.66e-15) and full EL ODE 2yy''+1+y'^2=0 (max err=1.05e-11)
[PASS] Endpoints (0,0), (l,0) and vertical tangents at A, B verified
[PASS] Travel time quadrature: exact=1.5394170507, param=1.5394170506, x-int=1.5394170507
[PASS] Hamilton-Jacobi calibration S(x,y)=theta*sqrt(a/g) verified (max |grad S|^2 err=1.10e-09)
[PASS] Strictly minimal vs competitors (times in units of sqrt(l/g)): cycloid=2.50663, semicircle=2.62262, parabola=2.61414, V-tunnel=2.82843
ALL CHECKS PASSED for IB_methods_sheet-3_qx07.
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  1. Check the sign of $dy/dx$ on the ascending half $\theta \in (\pi, 2\pi)$: careless taking of the positive square root $+\sqrt{(2a-y)/y}$ everywhere only describes the descending half $0 \le \theta \le \pi$. Using $\cot(\theta/2)$ handles both halves smoothly across the bottom $\theta = \pi$.
  2. Check the constant $a = l/(2\pi)$ vs $l/\pi$: a classic slip is forgetting the factor of $2$ in $x(2\pi) = 2\pi a$ (confusing half-arch $x(\pi) = \pi a$ with full arch $x(2\pi) = 2\pi a$), which gives the wrong factor inside $(2\pi l/g)^{1/2}$.
  3. Watch whether the existing solution rules out multiple cycloid arches ($N \ge 2$ gives $T_N = \sqrt{N}\,T_1 > T_1$) and justifies reduction to the vertical plane $z = 0$ and vertical tangents at $A, B$.
  4. Note that Beltrami's first integral only shows stationarity, not minimality; including the global minimality calibration (or at least explicitly remarking on minimality/convexity) strengthens the proof substantially.
