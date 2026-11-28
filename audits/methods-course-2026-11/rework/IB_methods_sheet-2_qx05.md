---
id: IB/methods/sheet-2/qx05
slug: IB_methods_sheet-2_qx05
source_archive: private/raw/sheets/IB/methods/2026-27/sheet-2.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-2_qx05_p1.png` (Question 5 on page 1 of *Mathematical Tripos Part IB, Methods, Example Sheet 2*, Michaelmas 2026, Prof. A.C.L. Ashton). Question 6 follows immediately on page 1, so there is no spillover onto page 2.

Every word, symbol, boundary condition, sign, constant ($A$, $L$, $c$, $\lambda$, $\frac{1}{2}L$), derivative ($\partial y/\partial t$), and delta function matches the printed PDF exactly. No discrepancies found (`transcription: ok`).

## Full Solution

### 1. Governing PDE, boundary conditions, and initial conditions for $t > 0$

For $t > 0$, the transverse displacement $y(x,t)$ of the uniform string satisfies the one-dimensional wave equation
\[
\frac{\partial^2 y}{\partial t^2} = c^2 \frac{\partial^2 y}{\partial x^2}, \qquad 0 < x < L, \quad t > 0,
\]
where the wave speed $c > 0$ is identified from the pre-strike fundamental mode $y(x,t) = A\sin(\pi x/L)\sin(\pi ct/L)$ for $t < 0$ (which indeed satisfies $y_{tt} = -(\pi c/L)^2 y = c^2 y_{xx}$). The fixed-end boundary conditions are
\[
y(0,t) = y(L,t) = 0 \qquad (t > 0).
\]

To obtain the initial data at $t = 0^+$, evaluate the pre-strike motion and its time derivative as $t \to 0^-$:
\[
y(x,0^-) = \lim_{t \to 0^-} A\sin\!\left(\frac{\pi x}{L}\right)\sin\!\left(\frac{\pi ct}{L}\right) = 0,
\]
\[
\frac{\partial y}{\partial t}(x,0^-) = \lim_{t \to 0^-} \frac{\pi c A}{L}\sin\!\left(\frac{\pi x}{L}\right)\cos\!\left(\frac{\pi ct}{L}\right) = \frac{\pi c A}{L}\sin\!\left(\frac{\pi x}{L}\right).
\]
An impulsive strike at $t = 0$ produces an instantaneous finite change in velocity without an instantaneous jump in displacement (mathematically, the forced equation $y_{tt} - c^2 y_{xx} = \lambda\delta(x-\frac{1}{2}L)\delta(t)$ gives continuity of $y$ across $t=0$ and a jump $\left[y_t\right]_{0^-}^{0^+} = \lambda\delta(x-\frac{1}{2}L)$). Therefore the initial conditions at $t = 0^+$ are
\[
y(x,0^+) = 0, \qquad \frac{\partial y}{\partial t}(x,0^+) = \frac{\pi c A}{L}\sin\!\left(\frac{\pi x}{L}\right) + \lambda\delta\!\left(x - \frac{L}{2}\right), \qquad 0 < x < L.
\]

### 2. Separation of variables for $t > 0$

Look for non-trivial separable solutions $y(x,t) = X(x)T(t)$ of $y_{tt} = c^2 y_{xx}$ satisfying $X(0) = X(L) = 0$. Substituting into the wave equation gives
\[
\frac{T''(t)}{c^2 T(t)} = \frac{X''(x)}{X(x)} = \mu,
\]
where $\mu$ is a separation constant. Checking all signs of $\mu$:
- If $\mu = \alpha^2 > 0$ ($\alpha > 0$), $X(x) = C_1 \cosh(\alpha x) + C_2 \sinh(\alpha x)$; then $X(0)=0 \implies C_1=0$ and $X(L)=0 \implies C_2\sinh(\alpha L)=0 \implies C_2=0$ (trivial).
- If $\mu = 0$, $X(x) = C_1 x + C_2$; then $X(0)=X(L)=0 \implies C_1=C_2=0$ (trivial).
- If $\mu = -k^2 < 0$ ($k > 0$), $X(x) = C_1\cos(kx) + C_2\sin(kx)$; then $X(0)=0 \implies C_1=0$, and $X(L)=0$ with $C_2 \neq 0$ requires $\sin(kL)=0$, i.e.
\[
k = k_n = \frac{n\pi}{L}, \qquad n = 1, 2, 3, \dots
\]
The corresponding temporal ODE is
\[
T_n''(t) + \left(\frac{n\pi c}{L}\right)^2 T_n(t) = 0 \implies T_n(t) = a_n\cos\!\left(\frac{n\pi ct}{L}\right) + b_n\sin\!\left(\frac{n\pi ct}{L}\right).
\]
By linearity and superposition, the general solution for $t > 0$ is
\[
y(x,t) = \sum_{n=1}^\infty \left[ a_n\cos\!\left(\frac{n\pi ct}{L}\right) + b_n\sin\!\left(\frac{n\pi ct}{L}\right) \right]\sin\!\left(\frac{n\pi x}{L}\right).
\]

### 3. Imposing the initial conditions at $t = 0^+$

Recall the orthogonality relation of the Fourier sine basis on $[0,L]$:
\[
\int_0^L \sin\!\left(\frac{n\pi x}{L}\right)\sin\!\left(\frac{m\pi x}{L}\right)\,dx = \frac{L}{2}\delta_{nm} \qquad (n, m \ge 1).
\]
First, $y(x,0^+) = 0$ gives
\[
\sum_{n=1}^\infty a_n \sin\!\left(\frac{n\pi x}{L}\right) = 0 \implies a_n = 0 \quad \text{for all } n \ge 1.
\]
Thus
\[
y(x,t) = \sum_{n=1}^\infty b_n \sin\!\left(\frac{n\pi x}{L}\right)\sin\!\left(\frac{n\pi ct}{L}\right).
\]
Differentiating with respect to $t$ and setting $t = 0^+$ gives (in the Fourier/distributional sense)
\[
\sum_{n=1}^\infty \frac{n\pi c}{L} b_n \sin\!\left(\frac{n\pi x}{L}\right) = \frac{\pi c A}{L}\sin\!\left(\frac{\pi x}{L}\right) + \lambda\delta\!\left(x - \frac{L}{2}\right).
\]
Multiplying both sides by $\sin(m\pi x/L)$ and integrating from $x=0$ to $x=L$ (valid since $x = \frac{1}{2}L$ lies strictly inside $(0,L)$):
\[
\frac{L}{2} \cdot \frac{m\pi c}{L} b_m = \frac{\pi c A}{L}\cdot\frac{L}{2}\delta_{m1} + \lambda \int_0^L \delta\!\left(x - \frac{L}{2}\right)\sin\!\left(\frac{m\pi x}{L}\right)\,dx,
\]
so
\[
\frac{m\pi c}{2} b_m = \frac{\pi c A}{2}\delta_{m1} + \lambda\sin\!\left(\frac{m\pi}{2}\right).
\]
Hence
\[
b_m = A\delta_{m1} + \frac{2\lambda}{m\pi c}\sin\!\left(\frac{m\pi}{2}\right).
\]
Evaluating $\sin(m\pi/2)$:
- For even $m = 2k$, $\sin(k\pi) = 0$, so $b_{2k} = 0$ (as expected by symmetry: the midpoint $x = \frac{1}{2}L$ is a node of every even harmonic, so striking the midpoint cannot excite even modes).
- For odd $m = 2k+1$ ($k = 0, 1, 2, \dots$), $\sin\!\left(\frac{(2k+1)\pi}{2}\right) = (-1)^k$, giving
\[
b_1 = A + \frac{2\lambda}{\pi c}, \qquad b_{2k+1} = \frac{2\lambda(-1)^k}{(2k+1)\pi c} \quad (k \ge 1).
\]

Therefore, for $t > 0$,
\[
y(x,t) = A\sin\!\left(\frac{\pi x}{L}\right)\sin\!\left(\frac{\pi ct}{L}\right) + \frac{2\lambda}{\pi c}\sum_{k=0}^\infty \frac{(-1)^k}{2k+1}\sin\!\left(\frac{(2k+1)\pi x}{L}\right)\sin\!\left(\frac{(2k+1)\pi ct}{L}\right).
\]

### 4. Equivalent d'Alembert (piecewise-constant) form of the impulse response

Using $\sin\alpha\sin\beta = \frac{1}{2}[\cos(\alpha-\beta)-\cos(\alpha+\beta)]$, the impulse contribution $y_1(x,t)$ can be rewritten as
\[
y_1(x,t) = \frac{\lambda}{\pi c}\sum_{k=0}^\infty \frac{(-1)^k}{2k+1}\left[\cos\!\left(\frac{(2k+1)\pi(x-ct)}{L}\right) - \cos\!\left(\frac{(2k+1)\pi(x+ct)}{L}\right)\right].
\]
Since the $2L$-periodic square wave $S(u) = \operatorname{sgn}(\cos(\pi u/L))$ (equal to $+1$ on $|u| < L/2$ and $-1$ on $L/2 < |u| < L$, and $0$ at jump points) has Fourier series $S(u) = \frac{4}{\pi}\sum_{k=0}^\infty \frac{(-1)^k}{2k+1}\cos\!\left(\frac{(2k+1)\pi u}{L}\right)$, we obtain the exact closed-form d'Alembert representation (away from the characteristic jump lines $x \pm ct \in \frac{1}{2}L + L\mathbb{Z}$):
\[
y(x,t) = A\sin\!\left(\frac{\pi x}{L}\right)\sin\!\left(\frac{\pi ct}{L}\right) + \frac{\lambda}{4c}\left[S(x-ct) - S(x+ct)\right].
\]
For $0 < ct < \frac{1}{2}L$, the impulsive term is a flat-topped rectangular pulse of height $\frac{\lambda}{2c}$ on $|x - \frac{1}{2}L| < ct$ and $0$ outside, expanding outward from the midpoint at speed $c$ before reflecting with inverted sign off the fixed ends $x=0$ and $x=L$.

## Final Answers

For $t > 0$ and $0 \le x \le L$,
\[
y(x,t) = A\sin\!\left(\frac{\pi x}{L}\right)\sin\!\left(\frac{\pi ct}{L}\right) + \frac{2\lambda}{\pi c}\sum_{k=0}^\infty \frac{(-1)^k}{2k+1}\sin\!\left(\frac{(2k+1)\pi x}{L}\right)\sin\!\left(\frac{(2k+1)\pi ct}{L}\right),
\]
equivalently
\[
y(x,t) = \left(A + \frac{2\lambda}{\pi c}\right)\sin\!\left(\frac{\pi x}{L}\right)\sin\!\left(\frac{\pi ct}{L}\right) + \frac{2\lambda}{\pi c}\sum_{k=1}^\infty \frac{(-1)^k}{2k+1}\sin\!\left(\frac{(2k+1)\pi x}{L}\right)\sin\!\left(\frac{(2k+1)\pi ct}{L}\right),
\]
or in d'Alembert form (away from characteristics $x \pm ct \in \frac{1}{2}L + L\mathbb{Z}$),
\[
y(x,t) = A\sin\!\left(\frac{\pi x}{L}\right)\sin\!\left(\frac{\pi ct}{L}\right) + \frac{\lambda}{4c}\left[\operatorname{sgn}\!\left(\cos\frac{\pi(x-ct)}{L}\right) - \operatorname{sgn}\!\left(\cos\frac{\pi(x+ct)}{L}\right)\right].
\]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-2_qx05_check.py`

Verbatim stdout (exit code 0):
```text
IB_methods_sheet-2_qx05_check: ALL PASSED
  max PDE finite-diff error (K=30): 2.487e-06
  max series vs d'Alembert error (K=6000): 1.968e-04
  distributional initial velocity error: 1.217e-09
  max smoothed-delta coefficient error (n=1..9): 3.571e-07
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`. Both separation of variables and d'Alembert's formula agree identically and pass all numerical checks (PDE, BCs, distributional velocity pairing, smoothed-delta projection, characteristic jump midpoint).
- **Phase B watchpoints:**
  - Check the factor of $2$ and $c$ in $\frac{2\lambda}{\pi c}$ (a common error is dropping $L$ cancellation or getting $\frac{2\lambda}{L}$ instead of $\frac{2\lambda}{\pi c}$; dimensionally $[\lambda] = L^2 T^{-1}$ since $[\lambda\delta(x-L/2)] = L T^{-1}$, so $[\lambda/c] = L$ is independent of $L$).
  - Check that $y(x,0^+) = 0$ (no jump in displacement from a finite velocity impulse), so all cosine-in-$t$ coefficients $a_n = 0$.
  - Check the sign of $(-1)^k = \sin\!\left(\frac{(2k+1)\pi}{2}\right)$ and the complete absence of even harmonics ($n = 2, 4, \dots$).
  - Check whether the existing solution correctly keeps the pre-existing $A\sin(\pi x/L)\sin(\pi ct/L)$ mode (i.e., $b_1 = A + \frac{2\lambda}{\pi c}$) rather than forgetting the initial velocity $\frac{\pi c A}{L}\sin(\pi x/L)$ from $t < 0$.
