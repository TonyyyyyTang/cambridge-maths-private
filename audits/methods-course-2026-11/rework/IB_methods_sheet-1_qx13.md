---
id: IB/methods/sheet-1/qx13
slug: IB_methods_sheet-1_qx13
source_archive: private/raw/sheets/IB/methods/2001-02/sheet-1.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against Question 7 on page 2 of `2001-02/sheet-1.pdf` (pre-rendered 200 DPI PNG `IB_methods_sheet-1_qx13_p2.png`).

- **Result:** `transcription: ok`.
- **Every word, symbol, sign, constant, initial condition, and forcing term matches exactly**, including:
  - Tension $T = \mu c^2$, resistive force per unit length $-2k\mu y_t$, constant $k = \pi c/\ell$.
  - Target PDE $c^2\frac{\partial^2 y}{\partial x^2} = \frac{\partial^2 y}{\partial t^2} + 2k\frac{\partial y}{\partial t}$.
  - Initial conditions $y(x, 0) = A\sin(\pi x/\ell)$ and $y_t(x, 0) = 0$.
  - Extra transverse force $F_0\sin(\pi x/\ell)\cos(\pi c t/\ell)$ per unit length in part (ii).
- **Purely typographic notes (not discrepancies):** the printed heading *Damped string motion.* is underlined as well as italicised, whereas `question.tex` uses `\emph{Damped string motion.}` and formats `(i)`/`(ii)` using an `enumerate` list.

---

## Full Solution

### Part (i)

#### 1. Derivation of the equation of motion

Consider an infinitesimal element of the stretched string between $x$ and $x + \delta x$ undergoing small transverse displacements $y(x, t)$ ($|y_x| \ll 1$), so that horizontal motion is negligible and the tension $T$ remains uniform along the string.
- The mass of the element is $\delta m = \mu\,\delta x$.
- If $\theta(x, t)$ is the angle of the string with the $x$-axis, then $\tan\theta = \frac{\partial y}{\partial x}$, and for small slopes $\sin\theta \approx \tan\theta = \frac{\partial y}{\partial x}$ (to leading order, neglecting $O(y_x^3)$).
- The net transverse restoring force due to the tension acting at $x + \delta x$ and $x$ is
  $$T\sin\theta(x + \delta x, t) - T\sin\theta(x, t) \approx T\left[\frac{\partial y}{\partial x}(x + \delta x, t) - \frac{\partial y}{\partial x}(x, t)\right] \approx T\frac{\partial^2 y}{\partial x^2}\delta x.$$
- The resistive force from the surrounding medium on the element of length $\delta x$ is $(-2k\mu y_t)\,\delta x$.

By Newton's second law for the transverse motion of the element:
$$(\mu\,\delta x)\frac{\partial^2 y}{\partial t^2} = T\frac{\partial^2 y}{\partial x^2}\delta x - 2k\mu\frac{\partial y}{\partial t}\delta x.$$
Dividing by $\mu\,\delta x$, taking the limit $\delta x \to 0$, and substituting $T/\mu = c^2$ gives
$$\frac{\partial^2 y}{\partial t^2} = c^2\frac{\partial^2 y}{\partial x^2} - 2k\frac{\partial y}{\partial t},$$
that is,
$$c^2\frac{\partial^2 y}{\partial x^2} = \frac{\partial^2 y}{\partial t^2} + 2k\frac{\partial y}{\partial t}.$$

#### 2. Solution of the initial-boundary value problem

The fixed ends impose the Dirichlet boundary conditions $y(0, t) = y(\ell, t) = 0$. Expanding $y(x, t)$ in the complete orthogonal spatial sine basis on $[0, \ell]$:
$$y(x, t) = \sum_{n=1}^\infty q_n(t)\sin\left(\frac{n\pi x}{\ell}\right).$$
Substituting into the PDE and using orthogonality of $\sin(n\pi x/\ell)$ on $(0, \ell)$, each mode amplitude $q_n(t)$ satisfies the damped harmonic oscillator equation
$$\ddot{q}_n + 2k\dot{q}_n + \omega_n^2 q_n = 0, \qquad \omega_n = \frac{n\pi c}{\ell} = nk,$$
where we used the given value $k = \frac{\pi c}{\ell}$.

The initial conditions $y(x, 0) = A\sin(\pi x/\ell)$ and $y_t(x, 0) = 0$ yield
$$q_1(0) = A, \quad \dot{q}_1(0) = 0; \qquad q_n(0) = 0, \quad \dot{q}_n(0) = 0 \quad (n \ge 2).$$
For every $n \ge 2$, the homogeneous linear ODE with zero initial data has the unique solution $q_n(t) \equiv 0$. (Uniqueness of the full PDE solution also follows directly from the energy functional $E(t) = \frac{1}{2}\int_0^\ell (w_t^2 + c^2 w_x^2)\,dx$ for the difference $w = y_1 - y_2$ of any two solutions, since $\dot{E}(t) = -2k\int_0^\ell w_t^2\,dx \le 0$ and $E(0) = 0$ imply $w \equiv 0$.)

For the fundamental mode $n = 1$, $\omega_1 = \frac{\pi c}{\ell} = k$, so the ODE is **critically damped**:
$$\ddot{q}_1 + 2k\dot{q}_1 + k^2 q_1 = 0.$$
Its auxiliary equation is $r^2 + 2kr + k^2 = (r + k)^2 = 0$, with repeated root $r = -k$. Hence the general solution is
$$q_1(t) = (C + Dt)e^{-kt}.$$
Applying the initial conditions:
- $q_1(0) = A \implies C = A$;
- $\dot{q}_1(t) = \bigl[D - k(C + Dt)\bigr]e^{-kt}$, so $\dot{q}_1(0) = D - kC = 0 \implies D = kA$.

Thus $q_1(t) = A(1 + kt)e^{-kt}$, and the displacement is
$$y(x, t) = A(1 + kt)e^{-kt}\sin\left(\frac{\pi x}{\ell}\right) = A\left(1 + \frac{\pi c t}{\ell}\right)e^{-\pi c t/\ell}\sin\left(\frac{\pi x}{\ell}\right).$$

---

### Part (ii)

With an additional transverse driving force $F_0\sin(\pi x/\ell)\cos(\pi c t/\ell)$ **per unit length**, Newton's second law on the element $\delta x$ becomes
$$\mu\frac{\partial^2 y}{\partial t^2} = T\frac{\partial^2 y}{\partial x^2} - 2k\mu\frac{\partial y}{\partial t} + F_0\sin\left(\frac{\pi x}{\ell}\right)\cos\left(\frac{\pi c t}{\ell}\right).$$
Dividing by $\mu$ and recalling $T/\mu = c^2$ and $k = \pi c/\ell$:
$$\frac{\partial^2 y}{\partial t^2} + 2k\frac{\partial y}{\partial t} - c^2\frac{\partial^2 y}{\partial x^2} = \frac{F_0}{\mu}\sin\left(\frac{\pi x}{\ell}\right)\cos(kt).$$

Since the forcing has spatial dependence purely in the fundamental mode $\sin(\pi x/\ell)$, we seek a steady-state forced (particular) oscillation of the form
$$y_p(x, t) = q_p(t)\sin\left(\frac{\pi x}{\ell}\right),$$
which satisfies the boundary conditions $y_p(0, t) = y_p(\ell, t) = 0$. Substituting into the forced PDE (with $c^2(\pi/\ell)^2 = k^2$) gives
$$\ddot{q}_p + 2k\dot{q}_p + k^2 q_p = \frac{F_0}{\mu}\cos(kt).$$

To solve this, complexify: let $z(t)$ satisfy $\ddot{z} + 2k\dot{z} + k^2 z = \frac{F_0}{\mu}e^{ikt}$ with $q_p(t) = \operatorname{Re}\bigl(z(t)\bigr)$, and try $z(t) = Z e^{ikt}$:
$$\bigl[(ik)^2 + 2k(ik) + k^2\bigr]Z = \frac{F_0}{\mu} \implies (-k^2 + 2ik^2 + k^2)Z = \frac{F_0}{\mu} \implies 2ik^2 Z = \frac{F_0}{\mu}.$$
(Because the string is driven at the undamped fundamental resonance frequency $\omega = \pi c/\ell = \omega_1 = k$, the inertia and elastic restoring terms cancel identically and the driving force is balanced purely by the damping term.) Hence
$$Z = \frac{F_0}{2ik^2\mu} = -\frac{iF_0}{2\mu k^2},$$
and taking the real part of $Z e^{ikt} = -\frac{iF_0}{2\mu k^2}(\cos kt + i\sin kt)$ gives
$$q_p(t) = \operatorname{Re}(Z e^{ikt}) = \frac{F_0}{2\mu k^2}\sin(kt).$$

Using $\mu k^2 = \mu\frac{\pi^2 c^2}{\ell^2} = \frac{\pi^2 T}{\ell^2}$, the resulting steady-state forced oscillation is
$$y_p(x, t) = \frac{F_0}{2\mu k^2}\sin\left(\frac{\pi x}{\ell}\right)\sin(kt) = \frac{F_0 \ell^2}{2\pi^2 \mu c^2}\sin\left(\frac{\pi x}{\ell}\right)\sin\left(\frac{\pi c t}{\ell}\right) = \frac{F_0 \ell^2}{2\pi^2 T}\sin\left(\frac{\pi x}{\ell}\right)\sin\left(\frac{\pi c t}{\ell}\right),$$
which has finite amplitude (regulated by damping) and lags the forcing $\cos(\pi ct/\ell)$ in phase by $\pi/2$.

*Remark on transients:* Every homogeneous mode $q_n^{(h)}(t)$ decays exponentially as $e^{-kt} \to 0$ as $t \to \infty$ (for $n = 1$, $(C_1 + D_1 t)e^{-kt}$; for $n \ge 2$, $e^{-kt}[C_n\cos(\sqrt{n^2-1}\,kt) + D_n\sin(\sqrt{n^2-1}\,kt)]$), so regardless of the initial conditions the motion approaches $y_p(x,t)$ exponentially fast. In particular, if the Part (i) initial conditions $y(x,0) = A\sin(\pi x/\ell)$, $y_t(x,0) = 0$ are also imposed simultaneously with the forcing, the complete transient + forced solution is
$$y(x, t) = \left[\left(A + \left(kA - \frac{F_0}{2\mu k}\right)t\right)e^{-kt} + \frac{F_0}{2\mu k^2}\sin(kt)\right]\sin\left(\frac{\pi x}{\ell}\right).$$

---

## Final Answers

- **(i)** Equation of motion: $c^2\dfrac{\partial^2 y}{\partial x^2} = \dfrac{\partial^2 y}{\partial t^2} + 2k\dfrac{\partial y}{\partial t}$ (with $c^2 = T/\mu$).
  Displacement for $k = \pi c/\ell$, $y(x,0) = A\sin(\pi x/\ell)$, $y_t(x,0) = 0$:
  $$y(x, t) = A\left(1 + \frac{\pi c t}{\ell}\right)e^{-\pi c t/\ell}\sin\left(\frac{\pi x}{\ell}\right) = A(1 + kt)e^{-kt}\sin\left(\frac{\pi x}{\ell}\right).$$
- **(ii)** Resulting (steady-state) forced oscillation:
  $$y_p(x, t) = \frac{F_0}{2\mu k^2}\sin\left(\frac{\pi x}{\ell}\right)\sin(kt) = \frac{F_0 \ell^2}{2\pi^2 \mu c^2}\sin\left(\frac{\pi x}{\ell}\right)\sin\left(\frac{\pi c t}{\ell}\right) = \frac{F_0 \ell^2}{2\pi^2 T}\sin\left(\frac{\pi x}{\ell}\right)\sin\left(\frac{\pi c t}{\ell}\right).$$

---

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_qx13_check.py` (exit code 0).

```text
Params (l=1.0, mu=1.0, c=1.0, A=1.0, F0=1.0): PDE1_res=3.38e-07, RK4_1=1.97e-14, PDE2_res=2.11e-08, RK4_2=9.14e-14, steady_diff(t>=20/k)=4.12e-08
Params (l=2.5, mu=1.7, c=3.2, A=0.8, F0=4.5): PDE1_res=2.91e-07, RK4_1=1.32e-14, PDE2_res=4.57e-08, RK4_2=4.75e-14, steady_diff(t>=20/k)=3.12e-08
Params (l=0.7, mu=0.4, c=1.9, A=-1.3, F0=2.1): PDE1_res=1.34e-06, RK4_1=4.24e-14, PDE2_res=4.44e-08, RK4_2=3.38e-14, steady_diff(t>=20/k)=5.78e-08
ALL CHECKS PASSED for IB/methods/sheet-1/qx13.
```

---

## Confidence & Edge-Case Notes

- **Confidence:** high.
- **Phase B watchpoints:**
  1. In part (i), verify recognition of **critical damping** for the $n = 1$ mode ($\omega_1 = \pi c/\ell = k$), yielding $(C + Dt)e^{-kt}$ rather than oscillatory or two distinct real exponentials, and verify $D = kA$ from $y_t(x,0) = 0$.
  2. In part (ii), check the factor of $\mu$! The forcing $F_0\sin(\pi x/\ell)\cos(\pi ct/\ell)$ is a force *per unit length*, so dividing Newton's second law by $\mu$ puts $F_0/\mu$ on the RHS (leading to $\frac{F_0}{2\mu k^2} = \frac{F_0\ell^2}{2\pi^2 T}$, NOT $\frac{F_0}{2k^2}$).
  3. Check the sign of the forced oscillation: $+\frac{F_0}{2\mu k^2}\sin(kt)$ (NOT minus!). Easy sign check: $2k\dot{q}_p = \frac{F_0}{\mu}\cos(kt)$ requires $\dot{q}_p > 0$ at $t = 0$, hence $+\sin(kt)$.
  4. Watch whether the existing solution also adds the transient solution for part (ii) with Part (i)'s initial conditions — both the steady-state forced oscillation and the full transient + forced formula are recorded above for clean comparison.
