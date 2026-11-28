---
id: IB/methods/sheet-3/qx16
slug: IB_methods_sheet-3_qx16
source_archive: private/raw/sheets/IB/methods/2006-07/sheet-3.ps
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-3_qx16_p2.png` (2006–07 IB Methods Sheet 3, Question 6, under *The one-dimensional wave equation*). Question 6 sits at the top of page 2 and finishes well above Question 7.

- **Status:** `ok` (exact match for every word, symbol, curly `\ell`, `\equiv`, `4V/\ell^2`, punctuation, and the roman subscript `\text{n odd}` under `\sum`).
- **Purely typographic note:** The printed heading *Modes on a string.* is italicised and underlined (sheet house style), encoded as `\emph{Modes on a string.}`.

## Full Solution

### 1. Wave equation and normal-mode expansion

Small transverse vibrations $y(x, t)$ of a uniform string of line density $\mu > 0$ and tension $T > 0$ satisfy the one-dimensional wave equation
\[
\mu \frac{\partial^2 y}{\partial t^2} = T \frac{\partial^2 y}{\partial x^2}
\iff
\frac{\partial^2 y}{\partial t^2} = c^2 \frac{\partial^2 y}{\partial x^2}, \qquad c = \sqrt{\frac{T}{\mu}} > 0,
\]
on $0 < x < \ell$, $t > 0$, subject to the fixed-end boundary conditions
\[
y(0, t) = 0, \qquad y(\ell, t) = 0 \qquad (t \ge 0).
\]

Seeking separated solutions $y(x, t) = X(x)\Theta(t)$ gives
\[
\frac{\Theta''(t)}{c^2 \Theta(t)} = \frac{X''(x)}{X(x)} = -\lambda,
\]
where $\lambda$ is a separation constant:
- If $\lambda < 0$ (write $\lambda = -\kappa^2$, $\kappa > 0$), then $X(x) = \alpha\cosh(\kappa x) + \beta\sinh(\kappa x)$; $X(0) = 0 \implies \alpha = 0$, and $X(\ell) = \beta\sinh(\kappa\ell) = 0 \implies \beta = 0$.
- If $\lambda = 0$, then $X(x) = \alpha x + \beta$; $X(0) = X(\ell) = 0 \implies \alpha = \beta = 0$.
- If $\lambda > 0$ (write $\lambda = k^2$, $k > 0$), then $X(x) = \alpha\cos(kx) + \beta\sin(kx)$; $X(0) = 0 \implies \alpha = 0$, and $X(\ell) = \beta\sin(k\ell) = 0$ yields non-trivial solutions if and only if
  \[
  k = k_n = \frac{n\pi}{\ell}, \qquad n = 1, 2, 3, \dots,
  \]
  with spatial eigenfunctions $X_n(x) = \sin(k_n x) = \sin(n\pi x/\ell)$ and corresponding angular frequencies
  \[
  \omega_n = c k_n = \frac{n\pi c}{\ell} = \frac{n\pi}{\ell}\sqrt{\frac{T}{\mu}}.
  \]

Superposing all normal modes gives the general solution satisfying the wave equation and boundary conditions:
\[
y(x, t) = \sum_{n=1}^{\infty} \bigl[A_n \cos(\omega_n t) + B_n \sin(\omega_n t)\bigr] \sin\!\left(\frac{n\pi x}{\ell}\right).
\]

### 2. Normal-mode amplitudes from the initial conditions

Applying the first initial condition $y(x, 0) = 0$ on $(0, \ell)$ yields
\[
\sum_{n=1}^{\infty} A_n \sin\!\left(\frac{n\pi x}{\ell}\right) = 0 \quad (0 < x < \ell).
\]
By the orthogonality relation
\[
\int_0^{\ell} \sin\!\left(\frac{m\pi x}{\ell}\right) \sin\!\left(\frac{n\pi x}{\ell}\right) dx = \frac{\ell}{2}\,\delta_{mn} \qquad (m, n \ge 1),
\]
we obtain $A_n = 0$ for all $n \ge 1$. Thus
\[
y(x, t) = \sum_{n=1}^{\infty} B_n \sin(\omega_n t) \sin\!\left(\frac{n\pi x}{\ell}\right).
\]

Differentiating with respect to $t$ and setting $t = 0$ (justified *a posteriori* since $B_n = O(n^{-4})$, so the series for $y$, $y_t$, $y_x$, and even $y_{tt}, y_{xx}$ converge uniformly and absolutely) gives
\[
y_t(x, 0) = \sum_{n=1}^{\infty} \omega_n B_n \sin\!\left(\frac{n\pi x}{\ell}\right) = \frac{4V}{\ell^2}\,x(\ell - x), \qquad 0 < x < \ell.
\]
Projecting onto $\sin(n\pi x/\ell)$ using orthogonality:
\[
\omega_n B_n = \frac{2}{\ell} \int_0^{\ell} \frac{4V}{\ell^2}\,x(\ell - x)\sin\!\left(\frac{n\pi x}{\ell}\right) dx.
\]
Substituting $u = x/\ell$ ($dx = \ell\,du$):
\[
\omega_n B_n = 8V \int_0^1 (u - u^2)\sin(n\pi u)\,du.
\]
Integrating by parts twice:
\begin{align*}
\int_0^1 (u - u^2)\sin(n\pi u)\,du
&= \underbrace{\left[-(u - u^2)\frac{\cos(n\pi u)}{n\pi}\right]_0^1}_{=\,0} + \frac{1}{n\pi}\int_0^1 (1 - 2u)\cos(n\pi u)\,du \\
&= \underbrace{\left[(1 - 2u)\frac{\sin(n\pi u)}{n^2\pi^2}\right]_0^1}_{=\,0} + \frac{2}{n^2\pi^2}\int_0^1 \sin(n\pi u)\,du \\
&= \frac{2}{n^2\pi^2}\left[-\frac{\cos(n\pi u)}{n\pi}\right]_0^1
= \frac{2\bigl(1 - (-1)^n\bigr)}{n^3\pi^3}
= \begin{cases}
\dfrac{4}{n^3\pi^3}, & n \text{ odd}, \\[8pt]
0, & n \text{ even}.
\end{cases}
\end{align*}
(Vanishing for even $n$ also follows immediately from symmetry: $x(\ell - x)$ is symmetric about $x = \ell/2$, whereas $\sin(n\pi x/\ell)$ is antisymmetric about $x = \ell/2$ for even $n$.)

Therefore
\[
\omega_n B_n = \begin{cases}
\dfrac{32V}{n^3\pi^3}, & n \text{ odd}, \\[8pt]
0, & n \text{ even},
\end{cases}
\]
and dividing by $\omega_n = n\pi c/\ell = (n\pi/\ell)\sqrt{T/\mu}$ gives the normal-mode displacement amplitudes:
\[
A_n = 0 \quad (\forall\, n \ge 1), \qquad
B_n = \begin{cases}
\dfrac{32 V \ell}{n^4 \pi^4 c} = \dfrac{32 V \ell}{n^4 \pi^4}\sqrt{\dfrac{\mu}{T}}, & n \text{ odd}, \\[10pt]
0, & n \text{ even}.
\end{cases}
\]
The displacement of the string is therefore
\[
y(x, t) = \frac{32 V \ell}{\pi^4 c} \sum_{n\text{ odd}} \frac{1}{n^4}\sin(\omega_n t)\sin\!\left(\frac{n\pi x}{\ell}\right).
\]

### 3. Kinetic and potential energies at time $t$

To avoid clashing with tension $T$ and initial midpoint speed $V$, denote the kinetic and potential energies at time $t$ by $K(t)$ and $U(t)$ respectively. For small transverse vibrations ($|y_x| \ll 1$, so arc length is $\int_0^\ell \sqrt{1 + y_x^2}\,dx \approx \ell + \frac{1}{2}\int_0^\ell y_x^2\,dx$):
\[
K(t) = \frac{1}{2}\mu \int_0^{\ell} y_t(x, t)^2\,dx, \qquad
U(t) = \frac{1}{2}T \int_0^{\ell} y_x(x, t)^2\,dx.
\]

The uniformly convergent series for $y_t$ and $y_x$ are
\[
y_t(x, t) = \sum_{n\text{ odd}} \omega_n B_n \cos(\omega_n t)\sin\!\left(\frac{n\pi x}{\ell}\right), \qquad
y_x(x, t) = \sum_{n\text{ odd}} k_n B_n \sin(\omega_n t)\cos\!\left(\frac{n\pi x}{\ell}\right).
\]
Using the orthogonality relations on $[0, \ell]$ for both sines and cosines,
\[
\int_0^{\ell} \sin\!\left(\frac{m\pi x}{\ell}\right)\sin\!\left(\frac{n\pi x}{\ell}\right) dx
= \int_0^{\ell} \cos\!\left(\frac{m\pi x}{\ell}\right)\cos\!\left(\frac{n\pi x}{\ell}\right) dx
= \frac{\ell}{2}\,\delta_{mn} \qquad (m, n \ge 1),
\]
all cross-terms vanish and we obtain:
- **Kinetic energy:**
  \[
  K(t) = \frac{1}{2}\mu \cdot \frac{\ell}{2} \sum_{n\text{ odd}} \omega_n^2 B_n^2 \cos^2(\omega_n t)
  = \frac{\mu \ell}{4} \sum_{n\text{ odd}} \left(\frac{32V}{n^3\pi^3}\right)^2 \cos^2(\omega_n t)
  = \frac{256\,\mu \ell V^2}{\pi^6} \sum_{n\text{ odd}} \frac{\cos^2(\omega_n t)}{n^6}.
  \]
- **Potential energy:** Since $T k_n^2 = (\mu c^2) k_n^2 = \mu (c k_n)^2 = \mu \omega_n^2$,
  \[
  U(t) = \frac{1}{2}T \cdot \frac{\ell}{2} \sum_{n\text{ odd}} k_n^2 B_n^2 \sin^2(\omega_n t)
  = \frac{\mu \ell}{4} \sum_{n\text{ odd}} \omega_n^2 B_n^2 \sin^2(\omega_n t)
  = \frac{256\,\mu \ell V^2}{\pi^6} \sum_{n\text{ odd}} \frac{\sin^2(\omega_n t)}{n^6}.
  \]

Adding $K(t)$ and $U(t)$ and using $\cos^2(\omega_n t) + \sin^2(\omega_n t) = 1$, the total energy is conserved in time:
\[
E(t) = K(t) + U(t) = \frac{256\,\mu \ell V^2}{\pi^6} \sum_{n\text{ odd}} \frac{1}{n^6}.
\]
(Energy conservation also follows directly from the PDE and fixed-end BCs: $\frac{dE}{dt} = \int_0^\ell (\mu y_t y_{tt} + T y_x y_{xt})\,dx = [T y_x y_t]_0^\ell + \int_0^\ell y_t(\mu y_{tt} - T y_{xx})\,dx = 0$.)

### 4. Comparison with the initial energy and evaluation of $\sum_{n\text{ odd}} n^{-6}$

At $t = 0$, $y(x, 0) = 0 \implies y_x(x, 0) = 0$, so $U(0) = 0$, and the initial total energy is purely kinetic:
\[
E(0) = K(0) = \frac{1}{2}\mu \int_0^{\ell} y_t(x, 0)^2\,dx
= \frac{1}{2}\mu \int_0^{\ell} \frac{16V^2}{\ell^4}\,x^2(\ell - x)^2\,dx
= \frac{8\mu V^2}{\ell^4} \int_0^{\ell} x^2(\ell - x)^2\,dx.
\]
Substituting $u = x/\ell$:
\[
\int_0^{\ell} x^2(\ell - x)^2\,dx
= \ell^5 \int_0^1 (u^2 - 2u^3 + u^4)\,du
= \ell^5 \left[\frac{u^3}{3} - \frac{u^4}{2} + \frac{u^5}{5}\right]_0^1
= \ell^5\left(\frac{1}{3} - \frac{1}{2} + \frac{1}{5}\right)
= \frac{\ell^5}{30}.
\]
Hence the initial energy of the string is
\[
E(0) = \frac{8\mu V^2}{\ell^4} \cdot \frac{\ell^5}{30} = \frac{4}{15}\,\mu \ell V^2.
\]
Equating $E(t) = E(0)$ (equivalently applying Parseval's identity to $y_t(x, 0)$) gives
\[
\frac{256\,\mu \ell V^2}{\pi^6} \sum_{n\text{ odd}} \frac{1}{n^6} = \frac{4}{15}\,\mu \ell V^2.
\]
Cancelling the non-zero factor $\mu \ell V^2$ yields
\[
\sum_{n\text{ odd}} \frac{1}{n^6} = \frac{4\pi^6}{15 \times 256} = \frac{\pi^6}{15 \times 64} = \frac{\pi^6}{960}.
\]

## Final Answers

- **Normal-mode displacement amplitudes** in $y(x, t) = \sum_{n=1}^\infty \bigl[A_n \cos(\omega_n t) + B_n \sin(\omega_n t)\bigr]\sin(n\pi x/\ell)$ (with $\omega_n = \frac{n\pi}{\ell}\sqrt{T/\mu}$):
  \[
  A_n = 0 \quad (\forall\, n \ge 1), \qquad
  B_n = \begin{cases}
  \dfrac{32 V \ell}{n^4 \pi^4}\sqrt{\dfrac{\mu}{T}}, & n \text{ odd}, \\[8pt]
  0, & n \text{ even}.
  \end{cases}
  \]
- **Kinetic and potential energies at time $t$:**
  \[
  K(t) = \frac{256\,\mu \ell V^2}{\pi^6} \sum_{n\text{ odd}} \frac{\cos^2(\omega_n t)}{n^6}, \qquad
  U(t) = \frac{256\,\mu \ell V^2}{\pi^6} \sum_{n\text{ odd}} \frac{\sin^2(\omega_n t)}{n^6}.
  \]
- **Initial energy and series identity:** $E(0) = \frac{4}{15}\,\mu \ell V^2 = E(t)$, giving $\displaystyle \sum_{n\text{ odd}} \frac{1}{n^6} = \frac{\pi^6}{960}$.

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-3_qx16_check.py`

```text
[PASS] Fourier sine coefficients omega_n*B_n verified for n=1..12 (max err=2.07e-12)
[PASS] BCs and ICs verified (max |y(x,0)|=0.00e+00, max |y_t(x,0)-f(x)|=1.73e-06)
[PASS] Wave PDE mu*y_tt = T*y_xx verified by finite differences (max res=3.50e-06)
[PASS] Modal K(t), U(t) match spatial quadrature (max K err=8.88e-16, max U err=5.55e-16) and E(t)=E(0)=4/15*mu*ell*V^2 (err=2.36e-13)
[PASS] Series identity verified: sum_(n odd) 1/n^6 = 1.001447076640942, pi^6/960 = 1.001447076640942 (diff=0.00e+00)
ALL CHECKS PASSED for IB_methods_sheet-3_qx16.
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  1. Watch for notation clashes: tension is $T$ and initial midpoint speed is $V$, so avoid calling kinetic energy $T(t)$ or potential energy $V(t)$!
  2. Watch for the factor of $\omega_n = n\pi c/\ell$ when extracting the displacement amplitude $B_n$ from the Fourier coefficient $\omega_n B_n = 32V/(n^3\pi^3)$ of $y_t(x, 0)$ — forgetting to divide by $\omega_n$ leaves $O(n^{-3})$ instead of $O(n^{-4})$.
  3. Watch for the factor of $\frac{1}{2}$ from $\int_0^\ell \sin^2(n\pi x/\ell)\,dx = \frac{\ell}{2}$ combined with the $\frac{1}{2}$ in $\frac{1}{2}\mu\int y_t^2 dx$ (giving $\frac{1}{4}\mu\ell \cdot 32^2 = 256\mu\ell$).
  4. Verify both even-$n$ vanishing and the exact arithmetic $\int_0^1 u^2(1-u)^2 du = 1/30 \implies E(0) = \frac{4}{15}\mu\ell V^2 \implies \frac{4/15}{256} = \frac{1}{960}$.
