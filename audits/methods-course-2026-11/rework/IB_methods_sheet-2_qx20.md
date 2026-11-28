---
id: IB/methods/sheet-2/qx20
slug: IB_methods_sheet-2_qx20
source_archive: private/raw/sheets/IB/methods/2007-08/sheet-2.pdf
pages_checked: [1, 2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-2_qx20_p1.png` (Question 4 on page 1 of Michaelmas 2007 Methods Examples II, Dr C. P. Caulfield). Also inspected page 2 (`IB_methods_sheet-2_qx24_p2.png`), which begins cleanly with Section *The Dirac delta function*, Question 5, confirming no spillover onto page 2.

**Verdict: `ok`.** All words, signs, exponents, derivatives ($y'''' + k^2 y'' = f(t)$), initial conditions ($y(0)=y'(0)=y''(0)=y'''(0)=0$), piecewise formula for $G(t,\tau)$, forcing $f(t)=e^{-t}$, and the bracketed hint match verbatim. Purely typographic differences only: the printed title *Higher order initial value problem$^*$.* is italicised and underlined (encoded as `\emph{...}`), and "is given by" is placed in the text line before the display rather than inside the display line.

## Full Solution

We consider the fourth-order linear initial value problem (with $' \equiv \frac{d}{dt}$ and non-zero constant $k$; we take $k \in \mathbb{R}\setminus\{0\}$ so that $k^{-2}, k^{-3}$, and $(1+k^2)^{-1}$ are well-defined):
\[
\mathcal{L}y \equiv y'''' + k^2 y'' = f(t), \qquad t > 0, \qquad y(0) = y'(0) = y''(0) = y'''(0) = 0.
\]

### 1. Construction of the Green's function $G(t,\tau)$

For fixed $\tau > 0$, the Green's function $G(t,\tau)$ satisfies
\[
\frac{\partial^4 G}{\partial t^4} + k^2 \frac{\partial^2 G}{\partial t^2} = \delta(t - \tau), \qquad t > 0,
\]
subject to the homogeneous initial conditions
\[
G(0,\tau) = G_t(0,\tau) = G_{tt}(0,\tau) = G_{ttt}(0,\tau) = 0,
\]
where subscripts denote partial derivatives with respect to $t$.

#### Region $0 < t < \tau$
On $0 < t < \tau$, $G(t,\tau)$ satisfies the homogeneous ODE $G_{tttt} + k^2 G_{tt} = 0$. Its auxiliary equation is
\[
m^4 + k^2 m^2 = m^2(m^2 + k^2) = 0,
\]
which for $k \neq 0$ has roots $m = 0$ (double) and $m = \pm ik$. Hence the general solution on $0 \le t < \tau$ is
\[
G(t,\tau) = A + Bt + C\cos(kt) + D\sin(kt).
\]
Differentiating and applying the initial conditions at $t = 0$:
\[
\begin{aligned}
G(0,\tau) &= A + C = 0, \\
G_t(0,\tau) &= B + kD = 0, \\
G_{tt}(0,\tau) &= -k^2 C = 0, \\
G_{ttt}(0,\tau) &= -k^3 D = 0.
\end{aligned}
\]
Since $k \neq 0$, the last two equations give $C = D = 0$, and the first two then give $A = B = 0$. Therefore,
\[
G(t,\tau) = 0 \qquad \text{for } 0 \le t < \tau.
\]

#### Continuity and jump conditions at $t = \tau$
In the differential equation $G_{tttt} + k^2 G_{tt} = \delta(t-\tau)$, the highest derivative $G_{tttt}$ must contain $\delta(t-\tau)$ and no higher derivatives of $\delta$. Consequently, $G_{ttt}(t,\tau)$ has a finite jump discontinuity at $t = \tau$, and is otherwise smooth and bounded near $t = \tau$. Integrating successively shows that $G_{tt}$, $G_t$, and $G$ are all continuous at $t = \tau$:
\[
\big[G\big]_{t=\tau^-}^{t=\tau^+} = 0, \qquad \big[G_t\big]_{t=\tau^-}^{t=\tau^+} = 0, \qquad \big[G_{tt}\big]_{t=\tau^-}^{t=\tau^+} = 0.
\]
Now integrate the ODE from $t = \tau - \varepsilon$ to $t = \tau + \varepsilon$ ($\varepsilon > 0$):
\[
\big[G_{ttt}(t,\tau)\big]_{\tau-\varepsilon}^{\tau+\varepsilon} + k^2 \big[G_t(t,\tau)\big]_{\tau-\varepsilon}^{\tau+\varepsilon} = \int_{\tau-\varepsilon}^{\tau+\varepsilon} \delta(t-\tau)\,dt = 1.
\]
Letting $\varepsilon \to 0^+$ and using the continuity of $G_t$ at $t = \tau$, we obtain the jump condition
\[
\big[G_{ttt}(t,\tau)\big]_{t=\tau^-}^{t=\tau^+} = 1.
\]
Since $G(t,\tau) = 0$ for $0 < t < \tau$, all left-hand limits at $t = \tau^-$ vanish, so the right-hand limits at $t = \tau^+$ are
\[
G(\tau^+,\tau) = 0, \qquad G_t(\tau^+,\tau) = 0, \qquad G_{tt}(\tau^+,\tau) = 0, \qquad G_{ttt}(\tau^+,\tau) = 1.
\]

#### Region $t > \tau$
Since the differential operator $\mathcal{L} = \frac{d^4}{dt^4} + k^2\frac{d^2}{dt^2}$ has constant coefficients (time-translation invariant) and the matching data at $t = \tau^+$ are independent of $\tau$, we write $G(t,\tau) = g(s)$ with $s = t - \tau > 0$. Then $g(s)$ solves
\[
g''''(s) + k^2 g''(s) = 0 \quad (s > 0), \qquad g(0) = g'(0) = g''(0) = 0, \quad g'''(0) = 1.
\]
Writing the general homogeneous solution as
\[
g(s) = A + Bs + C\cos(ks) + D\sin(ks),
\]
its derivatives are
\[
\begin{aligned}
g'(s) &= B - kC\sin(ks) + kD\cos(ks), \\
g''(s) &= -k^2 C\cos(ks) - k^2 D\sin(ks), \\
g'''(s) &= k^3 C\sin(ks) - k^3 D\cos(ks).
\end{aligned}
\]
Imposing the conditions at $s = 0^+$:
- $g''(0) = -k^2 C = 0 \implies C = 0$ (since $k \neq 0$),
- $g(0) = A + C = 0 \implies A = 0$,
- $g'''(0) = -k^3 D = 1 \implies D = -k^{-3}$,
- $g'(0) = B + kD = 0 \implies B = -kD = k^{-2}$.

Hence $g(s) = k^{-2} s - k^{-3}\sin(ks)$ for $s > 0$, giving
\[
G(t,\tau) = \begin{cases}
0, & 0 < t < \tau, \\
k^{-2}(t-\tau) - k^{-3}\sin k(t-\tau), & t > \tau,
\end{cases}
\]
as required.

### 2. Solution for $f(t) = e^{-t}$

For $t \ge 0$, the solution to the IVP is given by
\[
y(t) = \int_0^\infty G(t,\tau) f(\tau)\,d\tau = \int_0^t \left[k^{-2}(t-\tau) - k^{-3}\sin k(t-\tau)\right] e^{-\tau}\,d\tau.
\]
Substituting $s = t - \tau$ (so $\tau = t - s$, $d\tau = -ds$, and $s$ runs from $0$ to $t$):
\[
y(t) = e^{-t} \left[ k^{-2} \int_0^t s e^s\,ds - k^{-3} \int_0^t e^s \sin(ks)\,ds \right].
\]
Evaluating the two elementary integrals:
1. Integration by parts gives
   \[
   \int_0^t s e^s\,ds = \big[(s - 1)e^s\big]_0^t = (t - 1)e^t + 1.
   \]
2. Integrating by parts twice (or noting $\frac{d}{ds}\big[\frac{e^s(\sin ks - k\cos ks)}{1+k^2}\big] = e^s \sin ks$, valid whenever $1 + k^2 \neq 0$, in particular for all real $k$):
   \[
   I_2 = \int_0^t e^s \sin(ks)\,ds = \big[e^s\sin(ks)\big]_0^t - k\big[e^s\cos(ks)\big]_0^t - k^2 I_2
   \implies (1+k^2) I_2 = e^t\big(\sin(kt) - k\cos(kt)\big) + k,
   \]
   so
   \[
   \int_0^t e^s \sin(ks)\,ds = \frac{e^t\big(\sin(kt) - k\cos(kt)\big) + k}{1 + k^2}.
   \]

Substituting both results back into $y(t)$:
\[
\begin{aligned}
y(t) &= k^{-2}\big(t - 1 + e^{-t}\big) - k^{-3}\frac{\sin(kt) - k\cos(kt) + k e^{-t}}{1 + k^2} \\
&= \frac{t - 1}{k^2} + \left(\frac{1}{k^2} - \frac{1}{k^2(1+k^2)}\right)e^{-t} + \frac{k\cos(kt) - \sin(kt)}{k^3(1+k^2)} \\
&= \frac{e^{-t}}{1 + k^2} + \frac{t - 1}{k^2} + \frac{k\cos(kt) - \sin(kt)}{k^3(1 + k^2)}.
\end{aligned}
\]

### 3. Direct verification of the ODE and initial conditions

Differentiating $y(t)$ four times with respect to $t$:
\[
\begin{aligned}
y(t) &= \frac{e^{-t}}{1+k^2} + \frac{t-1}{k^2} + \frac{\cos(kt)}{k^2(1+k^2)} - \frac{\sin(kt)}{k^3(1+k^2)}, \\
y'(t) &= -\frac{e^{-t}}{1+k^2} + \frac{1}{k^2} - \frac{\sin(kt)}{k(1+k^2)} - \frac{\cos(kt)}{k^2(1+k^2)}, \\
y''(t) &= \frac{e^{-t}}{1+k^2} - \frac{\cos(kt)}{1+k^2} + \frac{\sin(kt)}{k(1+k^2)}, \\
y'''(t) &= -\frac{e^{-t}}{1+k^2} + \frac{k\sin(kt)}{1+k^2} + \frac{\cos(kt)}{1+k^2}, \\
y''''(t) &= \frac{e^{-t}}{1+k^2} + \frac{k^2\cos(kt)}{1+k^2} - \frac{k\sin(kt)}{1+k^2}.
\end{aligned}
\]

**Differential equation:**
\[
\begin{aligned}
y''''(t) + k^2 y''(t) &= \left(\frac{e^{-t}}{1+k^2} + \frac{k^2\cos(kt) - k\sin(kt)}{1+k^2}\right) + k^2\left(\frac{e^{-t}}{1+k^2} - \frac{\cos(kt)}{1+k^2} + \frac{\sin(kt)}{k(1+k^2)}\right) \\
&= \frac{1 + k^2}{1 + k^2} e^{-t} = e^{-t}.
\end{aligned}
\]
So the ODE is satisfied identically for all $t > 0$.

**Initial conditions at $t = 0$:**
\[
\begin{aligned}
y(0) &= \frac{1}{1+k^2} - \frac{1}{k^2} + \frac{1}{k^2(1+k^2)} - 0 = \frac{k^2 - (1+k^2) + 1}{k^2(1+k^2)} = 0, \\
y'(0) &= -\frac{1}{1+k^2} + \frac{1}{k^2} - 0 - \frac{1}{k^2(1+k^2)} = \frac{-k^2 + (1+k^2) - 1}{k^2(1+k^2)} = 0, \\
y''(0) &= \frac{1}{1+k^2} - \frac{1}{1+k^2} + 0 = 0, \\
y'''(0) &= -\frac{1}{1+k^2} + 0 + \frac{1}{1+k^2} = 0.
\end{aligned}
\]
Hence all four initial conditions are satisfied.

## Final Answers

- **Green's function** (for $k \neq 0$):
  \[
  G(t,\tau) = \begin{cases}
  0, & 0 < t < \tau, \\
  k^{-2}(t-\tau) - k^{-3}\sin k(t-\tau), & t > \tau.
  \end{cases}
  \]
- **Solution for $f(t) = e^{-t}$** (for $k \neq 0$, $k^2 \neq -1$):
  \[
  y(t) = \frac{e^{-t}}{1+k^2} + \frac{t-1}{k^2} + \frac{k\cos(kt) - \sin(kt)}{k^3(1+k^2)}.
  \]
- Directly verified that $y'''' + k^2 y'' = e^{-t}$ and $y(0) = y'(0) = y''(0) = y'''(0) = 0$.

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-2_qx20_check.py`

Verbatim stdout:
```text
IB_methods_sheet-2_qx20_check.py: ALL CHECKS PASSED
  max quadrature vs closed-form error (f=e^-t): 1.776e-15
  max analytic ODE residual:                    1.388e-16
  max IC error at t=0:                          4.441e-16
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  - Verify that the continuity of $G, G_t, G_{tt}$ at $t = \tau$ and the jump $[G_{ttt}]_{t=\tau^-}^{t=\tau^+} = 1$ (specifically showing why the $k^2 G_{tt}$ term contributes $k^2 [G_t]_{\tau^-}^{\tau^+} = 0$ upon integration across $\tau$) are properly justified, not merely asserted.
  - Check algebraic signs in the final solution for $f(t) = e^{-t}$: especially $+\frac{t-1}{k^2}$ (not $1-t$) and $+\frac{k\cos(kt)-\sin(kt)}{k^3(1+k^2)}$; also confirm the coefficient of $e^{-t}$ simplifies to $\frac{1}{1+k^2}$ (matching the particular integral).
  - Note implicitly $k \neq 0$ (and $k^2 \neq -1$, automatic for real $k$). As a limiting sanity check, as $k \to 0$, $G(t,\tau) \to \frac{1}{6}(t-\tau)^3$ and $y(t) \to e^{-t} - 1 + t - \frac{t^2}{2} + \frac{t^3}{6}$.
