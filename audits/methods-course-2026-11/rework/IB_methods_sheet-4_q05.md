---
id: IB/methods/sheet-4/q05
slug: IB_methods_sheet-4_q05
source_archive: private/raw/sheets/IB/methods/2025-26/sheet-4.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against the 200 DPI rendered page `IB_methods_sheet-4_q05_p1.png` (`private/raw/sheets/IB/methods/2025-26/sheet-4.pdf`, page 1, Question 5). Question 5 sits in the middle of page 1 and is followed by Questions 6–8 on the same page, so there is no page spillover.

- Initial value problem $(\dagger)$ matches printed layout, bold $\mathbf{R}^n$ and $\mathbf{x}$, $\kappa\Delta u$, $F(\mathbf{x},t)$, $f(\mathbf{x})$, and domains $\mathbf{R}^n \times (0,\infty)$ and $\mathbf{R}^n$.
- Inequality $\|u_1(\cdot,t) - u_2(\cdot,t)\| \le C\|f_1 - f_2\|$, norm definition $\|u(\cdot,t)\|^2 = \int |u(\mathbf{x},t)|^2\,\mathrm{d}^n\mathbf{x}\text{ etc.}$, and all surrounding wording ("for some constant $C>0$ you should determine. Deduce that $(\dagger)$ is well-posed with respect to the $\|\cdot\|$ norm.") match verbatim.

Verdict: **`transcription: ok`** (exact match).

## Full Solution

### 1. Setup, assumptions, and Fourier transform convention

We assume the thermal diffusivity satisfies $\kappa > 0$ (standard for the forward heat equation; see edge-case note on $\kappa < 0$ below) and work in the class of smooth solutions decaying sufficiently rapidly as $|\mathbf{x}| \to \infty$ (together with their spatial derivatives) so that spatial Fourier transforms exist, differentiation under the integral sign in $t$ is valid, and boundary terms vanish upon integration by parts.

We adopt the $n$-dimensional Fourier transform convention used on this sheet (Question 4): for a rapidly decaying function $g : \mathbb{R}^n \to \mathbb{C}$,
\[
\hat{g}(\boldsymbol{\lambda}) = \int_{\mathbb{R}^n} g(\mathbf{x}) e^{-\mathrm{i}\boldsymbol{\lambda}\cdot\mathbf{x}}\,\mathrm{d}^n\mathbf{x}, \qquad
g(\mathbf{x}) = \frac{1}{(2\pi)^n}\int_{\mathbb{R}^n} \hat{g}(\boldsymbol{\lambda}) e^{\mathrm{i}\boldsymbol{\lambda}\cdot\mathbf{x}}\,\mathrm{d}^n\boldsymbol{\lambda}.
\]
Under this convention, Parseval's theorem (Plancherel's identity) states that
\[
\|g\|^2 = \int_{\mathbb{R}^n} |g(\mathbf{x})|^2\,\mathrm{d}^n\mathbf{x} = \frac{1}{(2\pi)^n}\int_{\mathbb{R}^n} |\hat{g}(\boldsymbol{\lambda})|^2\,\mathrm{d}^n\boldsymbol{\lambda}.
\]
(Note: any other Fourier transform normalisation rescales Parseval's identity by a common positive constant, which cancels identically in the inequality below, so the resulting constant $C$ is convention-independent.)

### 2. Linearity and reduction to the homogeneous problem

Let $u_1$ and $u_2$ be solutions to $(\dagger)$ with the same forcing $F(\mathbf{x},t)$ and initial data $f_1$ and $f_2$, respectively. Define
\[
v(\mathbf{x},t) = u_1(\mathbf{x},t) - u_2(\mathbf{x},t), \qquad g(\mathbf{x}) = f_1(\mathbf{x}) - f_2(\mathbf{x}).
\]
By linearity of the operator $\partial_t - \kappa\Delta$, the forcing $F(\mathbf{x},t)$ cancels out and $v$ satisfies the homogeneous initial value problem
\[
\begin{cases}
v_t - \kappa \Delta v = 0, & (\mathbf{x},t) \in \mathbb{R}^n \times (0,\infty), \\
v(\mathbf{x},0) = g(\mathbf{x}), & \mathbf{x} \in \mathbb{R}^n.
\end{cases}
\]

### 3. Spatial Fourier transform and solution in frequency space

Let $\hat{v}(\boldsymbol{\lambda},t) = \int_{\mathbb{R}^n} v(\mathbf{x},t) e^{-\mathrm{i}\boldsymbol{\lambda}\cdot\mathbf{x}}\,\mathrm{d}^n\mathbf{x}$. For each $j \in \{1,\dots,n\}$, integrating by parts twice in $x_j$ (with boundary terms at infinity vanishing by the decay hypothesis) gives
\[
\widehat{\partial_{x_j}^2 v}(\boldsymbol{\lambda},t)
= \int_{\mathbb{R}^n} (\partial_{x_j}^2 v) e^{-\mathrm{i}\boldsymbol{\lambda}\cdot\mathbf{x}}\,\mathrm{d}^n\mathbf{x}
= \int_{\mathbb{R}^n} v \,\partial_{x_j}^2\!\left(e^{-\mathrm{i}\boldsymbol{\lambda}\cdot\mathbf{x}}\right)\mathrm{d}^n\mathbf{x}
= (-\mathrm{i}\lambda_j)^2 \hat{v}(\boldsymbol{\lambda},t)
= -\lambda_j^2 \hat{v}(\boldsymbol{\lambda},t).
\]
Summing over $j = 1, \dots, n$ yields
\[
\widehat{\Delta v}(\boldsymbol{\lambda},t) = -|\boldsymbol{\lambda}|^2 \hat{v}(\boldsymbol{\lambda},t), \qquad |\boldsymbol{\lambda}|^2 = \sum_{j=1}^n \lambda_j^2.
\]
Also $\widehat{v_t}(\boldsymbol{\lambda},t) = \frac{\partial \hat{v}}{\partial t}(\boldsymbol{\lambda},t)$. Taking the spatial Fourier transform of the homogeneous PDE and initial condition therefore gives, for each fixed $\boldsymbol{\lambda} \in \mathbb{R}^n$, the linear first-order ODE
\[
\frac{\partial \hat{v}}{\partial t}(\boldsymbol{\lambda},t) + \kappa |\boldsymbol{\lambda}|^2 \hat{v}(\boldsymbol{\lambda},t) = 0, \qquad \hat{v}(\boldsymbol{\lambda},0) = \hat{g}(\boldsymbol{\lambda}).
\]
Multiplying by the integrating factor $e^{\kappa |\boldsymbol{\lambda}|^2 t}$ gives
\[
\frac{\partial}{\partial t}\!\left(e^{\kappa |\boldsymbol{\lambda}|^2 t} \hat{v}(\boldsymbol{\lambda},t)\right) = 0
\implies \hat{v}(\boldsymbol{\lambda},t) = \hat{g}(\boldsymbol{\lambda}) e^{-\kappa |\boldsymbol{\lambda}|^2 t} \qquad \text{for all } t \ge 0.
\]

### 4. Application of Parseval's theorem and determination of $C$

Since $\kappa > 0$, $t > 0$, and $|\boldsymbol{\lambda}|^2 \ge 0$, we have $\kappa |\boldsymbol{\lambda}|^2 t \ge 0$ and hence
\[
0 < e^{-2\kappa |\boldsymbol{\lambda}|^2 t} \le 1 \qquad \text{for all } \boldsymbol{\lambda} \in \mathbb{R}^n,\; t > 0.
\]
Applying Parseval's theorem twice gives, for every $t > 0$,
\[
\begin{aligned}
\|u_1(\cdot,t) - u_2(\cdot,t)\|^2
= \|v(\cdot,t)\|^2
&= \frac{1}{(2\pi)^n}\int_{\mathbb{R}^n} |\hat{v}(\boldsymbol{\lambda},t)|^2\,\mathrm{d}^n\boldsymbol{\lambda} \\
&= \frac{1}{(2\pi)^n}\int_{\mathbb{R}^n} |\hat{g}(\boldsymbol{\lambda})|^2 e^{-2\kappa |\boldsymbol{\lambda}|^2 t}\,\mathrm{d}^n\boldsymbol{\lambda} \\
&\le \frac{1}{(2\pi)^n}\int_{\mathbb{R}^n} |\hat{g}(\boldsymbol{\lambda})|^2\,\mathrm{d}^n\boldsymbol{\lambda} \\
&= \|g\|^2 = \|f_1 - f_2\|^2.
\end{aligned}
\]
Taking non-negative square roots gives
\[
\|u_1(\cdot,t) - u_2(\cdot,t)\| \le \|f_1 - f_2\| \qquad \text{for all } t > 0.
\]
Thus the stated inequality holds with
\[
C = 1
\]
(uniformly for all $t > 0$).

**Sharpness of $C = 1$**: No smaller constant $C < 1$ can hold for all initial data $f_1, f_2$, even at a fixed $t > 0$. Indeed, taking $g_\varepsilon(\mathbf{x}) = e^{-\varepsilon |\mathbf{x}|^2}$ ($\varepsilon > 0$) gives $\hat{g}_\varepsilon(\boldsymbol{\lambda}) = (\pi/\varepsilon)^{n/2} e^{-|\boldsymbol{\lambda}|^2/(4\varepsilon)}$, and an exact Gaussian integral shows
\[
\frac{\|v_\varepsilon(\cdot,t)\|}{\|g_\varepsilon\|} = (1 + 4\kappa\varepsilon t)^{-n/4} \longrightarrow 1 \qquad \text{as } \varepsilon \to 0^+.
\]
Hence $C = 1$ is the optimal constant.

### 5. Well-posedness of $(\dagger)$ with respect to $\|\cdot\|$

Recall that a problem is **well-posed** (in the sense of Hadamard) if (1) a solution exists, (2) the solution is unique, and (3) the solution depends continuously on the initial data (stability) in the chosen norm.

1. **Existence**: Taking the spatial Fourier transform directly of the inhomogeneous problem $(\dagger)$ yields
   \[
   \frac{\partial \hat{u}}{\partial t}(\boldsymbol{\lambda},t) + \kappa|\boldsymbol{\lambda}|^2 \hat{u}(\boldsymbol{\lambda},t) = \hat{F}(\boldsymbol{\lambda},t), \qquad \hat{u}(\boldsymbol{\lambda},0) = \hat{f}(\boldsymbol{\lambda}),
   \]
   whose solution via the integrating factor $e^{\kappa|\boldsymbol{\lambda}|^2 t}$ is
   \[
   \hat{u}(\boldsymbol{\lambda},t) = \hat{f}(\boldsymbol{\lambda}) e^{-\kappa|\boldsymbol{\lambda}|^2 t} + \int_0^t \hat{F}(\boldsymbol{\lambda},s) e^{-\kappa|\boldsymbol{\lambda}|^2(t-s)}\,\mathrm{d}s.
   \]
   Since $\widehat{K_t}(\boldsymbol{\lambda}) = e^{-\kappa|\boldsymbol{\lambda}|^2 t}$ is the Fourier transform of the $n$-dimensional heat kernel $K_t(\mathbf{x}) = (4\pi\kappa t)^{-n/2}\exp(-|\mathbf{x}|^2/(4\kappa t))$, inverting the Fourier transform via the convolution theorem provides the explicit solution (Duhamel's formula)
   \[
   u(\mathbf{x},t) = (K_t * f)(\mathbf{x}) + \int_0^t \bigl(K_{t-s} * F(\cdot,s)\bigr)(\mathbf{x})\,\mathrm{d}s,
   \]
   which exists (for suitably regular/decaying $f, F$, e.g. $f \in L^2(\mathbb{R}^n)$ and $\int_0^t \|F(\cdot,s)\|\,\mathrm{d}s < \infty$, with $\|u(\cdot,t)\| \le \|f\| + \int_0^t \|F(\cdot,s)\|\,\mathrm{d}s$).

2. **Uniqueness**: Suppose $u_1$ and $u_2$ are two solutions in this decaying/$L^2$ class corresponding to the same initial data $f_1 = f_2 = f$ (and same $F$). Then $\|f_1 - f_2\| = 0$, so the inequality gives
   \[
   0 \le \|u_1(\cdot,t) - u_2(\cdot,t)\| \le \|f_1 - f_2\| = 0 \implies \int_{\mathbb{R}^n} |u_1(\mathbf{x},t) - u_2(\mathbf{x},t)|^2\,\mathrm{d}^n\mathbf{x} = 0
   \]
   for every $t > 0$. For continuous solutions, this forces $u_1(\mathbf{x},t) = u_2(\mathbf{x},t)$ for all $(\mathbf{x},t) \in \mathbb{R}^n \times (0,\infty)$ (and in $L^2(\mathbb{R}^n)$, equality holds almost everywhere). Hence the solution is unique.

3. **Continuous dependence on the initial data (stability)**: For any $\varepsilon > 0$, choosing $\delta = \varepsilon$ ensures that whenever $\|f_1 - f_2\| < \delta$, we have
   \[
   \|u_1(\cdot,t) - u_2(\cdot,t)\| \le \|f_1 - f_2\| < \varepsilon \qquad \text{for all } t > 0.
   \]
   Thus the solution map $f \mapsto u(\cdot,t)$ is (Lipschitz) continuous with respect to the $L^2$ norm $\|\cdot\|$, uniformly in $t > 0$.

Since existence, uniqueness, and continuous dependence on the initial data all hold, $(\dagger)$ is well-posed with respect to the $\|\cdot\|$ norm.

## Final Answers

- **Stability bound**: $\|u_1(\cdot,t) - u_2(\cdot,t)\| \le \|f_1 - f_2\|$ for all $t > 0$ (assuming $\kappa > 0$), so **$C = 1$** (sharp).
- **Well-posedness**: $(\dagger)$ is well-posed in the $L^2(\mathbb{R}^n)$ norm $\|\cdot\|$: existence follows from the Fourier/Duhamel representation, uniqueness follows by setting $f_1 = f_2$, and continuous dependence (uniform Lipschitz stability) follows immediately from $\|u_1(\cdot,t) - u_2(\cdot,t)\| \le \|f_1 - f_2\|$.

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-4_q05_check.py`

```text
IB_methods_sheet-4_q05 verification PASSED:
  Parseval error (1D)             = 8.882e-16
  Fourier vs convolution diff     = 2.776e-16
  1D heat PDE FD residual         = 5.944e-08
  1D contraction ratios ||v(t)||/||g||: [(0.05, 0.938664), (0.45, 0.766769), (2.0, 0.597372)]
  2D ||f1-f2|| = 0.436351, ||u1(0.5)-u2(0.5)|| = 0.127976 <= ||f1-f2||
  Sharpness of C=1 (eps=1e-5, t=1): {1: 0.999993, 2: 0.999986, 3: 0.999979, 4: 0.999972}
  Backwards heat (kappa=-0.1, k=10, t=0.1) amplification = 2.7884 > 1
```

## Confidence & Edge-Case Notes

- **Confidence**: high.
- **Phase B watchpoints**:
  1. Check that $\kappa > 0$ is explicitly noted (if $\kappa < 0$, the factor $e^{-2\kappa|\boldsymbol{\lambda}|^2 t}$ grows unboundedly as $|\boldsymbol{\lambda}| \to \infty$ and the backwards heat equation is severely ill-posed).
  2. Check that the inhomogeneous term $F(\mathbf{x},t)$ is explicitly cancelled by considering $v = u_1 - u_2$ (a common oversight is ignoring $F$ without mentioning why it disappears).
  3. Check that all three Hadamard criteria (**existence**, **uniqueness**, and **continuous dependence on data**) are addressed when deducing well-posedness — not just stability alone!
  4. Verify that $C = 1$ (any $C \ge 1$ is technically true, so check the existing solution obtains $C = 1$ and does not accidentally leave a spurious factor of $(2\pi)^n$ or $(2\pi)^{n/2}$ from inconsistent Parseval normalisation!).
