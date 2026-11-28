---
id: IB/methods/sheet-3/q03
slug: IB_methods_sheet-3_q03
source_archive: private/raw/sheets/IB/methods/2025-26/sheet-3.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-3_q03_p1.png` (Question 3 on page 1 of Michaelmas 2025 Methods Example Sheet 3, Prof. A. C. L. Ashton; cleanly followed on the same page by Questions 4 and 5).

**Verdict: `ok`.** Every word, symbol, sign, domain (`0 < x < 1`), boundary condition (`y(0) = y(1) = 0,`), and even the omitted article ("to boundary value problem") and trailing comma after `y(0) = y(1) = 0,` match verbatim. (The minor horizontal alignment of `y(0) = y(1) = 0,` inside `\begin{aligned}` is purely typographic.)

## Full Solution

We consider the boundary value problem
\[
\mathcal{L}y \equiv y'' - 2y' + y = f(x), \qquad 0 < x < 1, \qquad y(0) = y(1) = 0,
\]
for a general (say continuous) forcing function $f$ on $[0,1]$.

### 1. Homogeneous problem and existence/uniqueness

The homogeneous equation is $y'' - 2y' + y = 0$, with auxiliary equation
\[
m^2 - 2m + 1 = (m - 1)^2 = 0,
\]
having a repeated root $m = 1$. Its general solution is
\[
y_h(x) = (A + Bx)e^x.
\]
Imposing the boundary conditions:
- $y_h(0) = A = 0$,
- $y_h(1) = B e = 0 \implies B = 0$.

Thus the homogeneous BVP admits only the trivial solution $y \equiv 0$. Consequently, the BVP has a unique Green's function and a unique solution for any continuous forcing $f$.

### 2. Construction of the Green's function $G(x,\xi)$

For fixed $\xi \in (0,1)$, we seek $G(x,\xi)$ satisfying
\[
\frac{\partial^2 G}{\partial x^2} - 2\frac{\partial G}{\partial x} + G = \delta(x - \xi), \qquad 0 < x < 1, \qquad G(0,\xi) = G(1,\xi) = 0.
\]
(Note: $\mathcal{L}$ is **not** self-adjoint as written, so $G(x,\xi)$ will **not** be symmetric in $x$ and $\xi$; we also note the equivalent self-adjoint formulation below.)

For $x \neq \xi$, $G(x,\xi)$ satisfies the homogeneous ODE $G_{xx} - 2G_x + G = 0$:
- On $0 \le x < \xi$, the homogeneous solution satisfying $G(0,\xi) = 0$ is proportional to $y_1(x) = x e^x$.
- On $\xi < x \le 1$, the homogeneous solution satisfying $G(1,\xi) = 0$ is proportional to $y_2(x) = (1 - x)e^x$.

Hence we take
\[
G(x,\xi) = \begin{cases}
A(\xi)\, x e^x, & 0 \le x < \xi, \\
B(\xi)\, (1 - x) e^x, & \xi < x \le 1.
\end{cases}
\]

#### Continuity and jump conditions at $x = \xi$
1. **Continuity of $G$ at $x = \xi$:** Since $G_{xx}$ contains only $\delta(x-\xi)$ (and no $\delta'$), $G_x$ has a finite jump at $x = \xi$, and therefore $G$ is continuous at $x = \xi$:
   \[
   G(\xi^+,\xi) = G(\xi^-,\xi) \implies B(\xi)(1-\xi)e^\xi = A(\xi)\xi e^\xi \implies B(\xi)(1-\xi) = A(\xi)\xi.
   \]
   We may therefore write $A(\xi) = C(\xi)(1-\xi)$ and $B(\xi) = C(\xi)\xi$ for some function $C(\xi)$, giving
   \[
   G(x,\xi) = \begin{cases}
   C(\xi)\,(1-\xi) x e^x, & 0 \le x \le \xi, \\
   C(\xi)\,\xi (1-x) e^x, & \xi \le x \le 1.
   \end{cases}
   \]

2. **Jump in $G_x$ at $x = \xi$:** Integrating $G_{xx} - 2G_x + G = \delta(x-\xi)$ over $[\xi-\varepsilon, \xi+\varepsilon]$ and letting $\varepsilon \to 0^+$:
   \[
   \big[G_x\big]_{x=\xi^-}^{x=\xi^+} - 2\big[G\big]_{x=\xi^-}^{x=\xi^+} + \lim_{\varepsilon\to 0^+}\int_{\xi-\varepsilon}^{\xi+\varepsilon} G(x,\xi)\,dx = 1.
   \]
   Since $G$ is continuous (and bounded) at $x = \xi$, the last two terms vanish, yielding
   \[
   G_x(\xi^+,\xi) - G_x(\xi^-,\xi) = 1.
   \]
   Computing the one-sided derivatives:
   - For $x < \xi$: $\dfrac{\partial}{\partial x}(x e^x) = (1+x)e^x$, so $G_x(\xi^-,\xi) = C(\xi)(1-\xi)(1+\xi)e^\xi = C(\xi)(1-\xi^2)e^\xi$.
   - For $x > \xi$: $\dfrac{\partial}{\partial x}\big((1-x)e^x\big) = -x e^x$, so $G_x(\xi^+,\xi) = -C(\xi)\xi^2 e^\xi$.

   Imposing the unit jump condition:
   \[
   G_x(\xi^+,\xi) - G_x(\xi^-,\xi) = C(\xi)e^\xi\big[-\xi^2 - (1 - \xi^2)\big] = -C(\xi)e^\xi = 1
   \implies C(\xi) = -e^{-\xi}.
   \]

Substituting $C(\xi) = -e^{-\xi}$ gives the Green's function:
\[
G(x,\xi) = \begin{cases}
-x(1-\xi)e^{x-\xi}, & 0 \le x \le \xi, \\
-\xi(1-x)e^{x-\xi}, & \xi \le x \le 1.
\end{cases}
\]

*Remark (self-adjoint cross-check):* Multiplying the ODE by the integrating factor $e^{-2x}$ puts it in Sturm–Liouville form $\tilde{\mathcal{L}}y \equiv (e^{-2x}y')' + e^{-2x}y = e^{-2x}f(x)$. Its self-adjoint Green's function $\tilde{G}(x,\xi)$ satisfies $\tilde{\mathcal{L}}_x \tilde{G} = \delta(x-\xi)$ with jump $[\tilde{G}_x]_{\xi^-}^{\xi^+} = 1/p(\xi) = e^{2\xi}$, giving the manifestly symmetric kernel
\[
\tilde{G}(x,\xi) = e^{2\xi} G(x,\xi) = \begin{cases}
-x(1-\xi)e^{x+\xi}, & 0 \le x \le \xi, \\
-\xi(1-x)e^{x+\xi}, & \xi \le x \le 1,
\end{cases}
\]
and $y(x) = \int_0^1 \tilde{G}(x,\xi)\,e^{-2\xi}f(\xi)\,d\xi = \int_0^1 G(x,\xi) f(\xi)\,d\xi$, in exact agreement.

### 3. General solution of the BVP and direct verification

The unique solution for a general continuous forcing $f(x)$ is
\[
\begin{aligned}
y(x) &= \int_0^1 G(x,\xi) f(\xi)\,d\xi \\
&= \int_0^x \big(-\xi(1-x)e^{x-\xi}\big) f(\xi)\,d\xi + \int_x^1 \big(-x(1-\xi)e^{x-\xi}\big) f(\xi)\,d\xi \\
&= -e^x\left[(1-x)\int_0^x \xi e^{-\xi} f(\xi)\,d\xi + x\int_x^1 (1-\xi)e^{-\xi} f(\xi)\,d\xi\right].
\end{aligned}
\]

**Direct verification:** Write $y(x) = e^x u(x)$, where
\[
u(x) = -(1-x)\int_0^x \xi e^{-\xi} f(\xi)\,d\xi - x\int_x^1 (1-\xi)e^{-\xi} f(\xi)\,d\xi.
\]
Clearly $u(0) = 0$ and $u(1) = 0$, so $y(0) = y(1) = 0$. Differentiating $u(x)$ via the Leibniz rule:
\[
\begin{aligned}
u'(x) &= \int_0^x \xi e^{-\xi} f(\xi)\,d\xi - (1-x)x e^{-x} f(x) - \int_x^1 (1-\xi)e^{-\xi} f(\xi)\,d\xi + x(1-x)e^{-x} f(x) \\
&= \int_0^x \xi e^{-\xi} f(\xi)\,d\xi - \int_x^1 (1-\xi)e^{-\xi} f(\xi)\,d\xi,
\end{aligned}
\]
and differentiating once more:
\[
u''(x) = x e^{-x} f(x) + (1-x)e^{-x} f(x) = e^{-x} f(x).
\]
Since $y' = e^x(u'+u)$ and $y'' = e^x(u''+2u'+u)$, we obtain
\[
y'' - 2y' + y = e^x u'' = e^x\big(e^{-x} f(x)\big) = f(x),
\]
confirming both the ODE and boundary conditions identically.

## Final Answers

- **Green's function** (for $\mathcal{L}G = \delta(x-\xi)$ with $G(0,\xi)=G(1,\xi)=0$):
  \[
  G(x,\xi) = \begin{cases}
  -x(1-\xi)e^{x-\xi}, & 0 \le x \le \xi, \\
  -\xi(1-x)e^{x-\xi}, & \xi \le x \le 1.
  \end{cases}
  \]
- **Solution to the BVP** (unique for any continuous $f$):
  \[
  y(x) = -e^x\left[(1-x)\int_0^x \xi e^{-\xi} f(\xi)\,d\xi + x\int_x^1 (1-\xi)e^{-\xi} f(\xi)\,d\xi\right].
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-3_q03_check.py`

Verbatim stdout:
```text
IB_methods_sheet-3_q03_check.py: ALL CHECKS PASSED
  max exact vs Green-integral error:  1.110e-16
  max FD ODE residual (cos(2x)+x^2):  1.916e-09
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  - **Non-self-adjoint trap:** The operator $y'' - 2y' + y$ is NOT self-adjoint, so $G(x,\xi)$ solving $G_{xx}-2G_x+G = \delta(x-\xi)$ is NOT symmetric (it satisfies $G(x,\xi) = e^{2(x-\xi)} G(\xi,x)$). If a solution converts to Sturm–Liouville form $(e^{-2x}y')' + e^{-2x}y = e^{-2x}f(x)$ first, it MUST include both the jump $[\tilde{G}_x]_{\xi^-}^{\xi^+} = e^{2\xi}$ AND the factor $e^{-2\xi}$ multiplying $f(\xi)$ inside the solution integral!
  - **Sign check:** $G(x,\xi) \le 0$ on $[0,1]^2$ (for $f > 0$, $u'' = e^{-x}f > 0$ with $u(0)=u(1)=0$ forces $u < 0$, hence $y = e^x u < 0$). Check that the overall minus sign is present.
  - **Branch swap in the integral:** In $\int_0^x$, $\xi < x$ (so $x > \xi$, using $-\xi(1-x)e^{x-\xi}$); in $\int_x^1$, $\xi > x$ (so $x < \xi$, using $-x(1-\xi)e^{x-\xi}$). Swapping these two branches is a common error.
