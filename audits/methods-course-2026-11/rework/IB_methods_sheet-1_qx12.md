---
id: IB/methods/sheet-1/qx12
slug: IB_methods_sheet-1_qx12
source_archive: private/raw/sheets/IB/methods/2001-02/sheet-1.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against Question 6 on page 2 of `2001-02/sheet-1.pdf` (pre-rendered 200 DPI PNG `IB_methods_sheet-1_qx12_p2.png`).

- **Result:** `transcription: ok`.
- **Every word, symbol, sign, coefficient, boundary/initial condition, and displayed formula matches exactly**, including:
  - $y(x,0) = 0$, $y_t(x,0) = \frac{4V}{\ell^2}x(\ell - x)$ for $0 < x < \ell$, where $y_t \equiv \partial y/\partial t$.
  - Target identity $\sum_{\text{n odd}} \frac{1}{n^6} = \frac{\pi^6}{960}$.
- **Purely typographic note (not a discrepancy):** the printed heading *Modes on a string.* is underlined as well as italicised, whereas `question.tex` uses `\emph{Modes on a string.}`.

---

## Full Solution

### 1. Wave equation and normal modes

For small transverse vibrations of a uniform stretched string of mass per unit length $\mu$ and constant tension $T$, the displacement $y(x,t)$ satisfies the one-dimensional wave equation
$$\frac{\partial^2 y}{\partial t^2} = c^2 \frac{\partial^2 y}{\partial x^2}, \qquad c = \sqrt{\frac{T}{\mu}},$$
on $0 < x < \ell$, subject to the fixed-end Dirichlet boundary conditions
$$y(0, t) = 0, \qquad y(\ell, t) = 0 \qquad (t \ge 0).$$

Seek separated solutions $y(x,t) = X(x)\Theta(t)$:
$$\frac{\Theta''(t)}{c^2 \Theta(t)} = \frac{X''(x)}{X(x)} = -\lambda,$$
where $\lambda$ is a separation constant. For the spatial boundary-value problem $X'' + \lambda X = 0$ with $X(0) = X(\ell) = 0$:
- If $\lambda < 0$, write $\lambda = -\kappa^2$ ($\kappa > 0$): $X(x) = P\cosh(\kappa x) + Q\sinh(\kappa x)$; $X(0) = 0 \implies P = 0$, and $X(\ell) = Q\sinh(\kappa\ell) = 0 \implies Q = 0$ (trivial).
- If $\lambda = 0$: $X(x) = P + Qx$; $X(0) = X(\ell) = 0 \implies P = Q = 0$ (trivial).
- If $\lambda > 0$, write $\lambda = k^2$ ($k > 0$): $X(x) = P\cos(kx) + Q\sin(kx)$; $X(0) = 0 \implies P = 0$, and $X(\ell) = Q\sin(k\ell) = 0$ yields non-trivial solutions if and only if $k\ell = n\pi$ for $n \in \mathbb{Z}_{>0}$, giving eigenvalues and angular frequencies
  $$k_n = \frac{n\pi}{\ell}, \qquad \omega_n = c k_n = \frac{n\pi c}{\ell} = \frac{n\pi}{\ell}\sqrt{\frac{T}{\mu}}, \qquad n = 1, 2, 3, \dots$$
  with spatial modes $\sin\left(\frac{n\pi x}{\ell}\right)$.

The temporal equation $\Theta'' + \omega_n^2 \Theta = 0$ gives $\Theta_n(t) = A_n \cos(\omega_n t) + B_n \sin(\omega_n t)$. By linearity, the general solution is
$$y(x, t) = \sum_{n=1}^\infty \bigl[A_n \cos(\omega_n t) + B_n \sin(\omega_n t)\bigr]\sin\left(\frac{n\pi x}{\ell}\right).$$

### 2. Amplitudes of the normal modes

Impose the initial condition $y(x, 0) = 0$ for $0 < x < \ell$:
$$\sum_{n=1}^\infty A_n \sin\left(\frac{n\pi x}{\ell}\right) = 0.$$
By the orthogonality relation
$$\int_0^\ell \sin\left(\frac{m\pi x}{\ell}\right)\sin\left(\frac{n\pi x}{\ell}\right)\,dx = \frac{\ell}{2}\delta_{mn} \qquad (m, n \ge 1),$$
we obtain $A_n = 0$ for all $n \ge 1$. Hence
$$y(x, t) = \sum_{n=1}^\infty B_n \sin(\omega_n t)\sin\left(\frac{n\pi x}{\ell}\right).$$

Now impose the initial velocity condition $y_t(x, 0) = \frac{4V}{\ell^2}x(\ell - x)$ for $0 < x < \ell$:
$$y_t(x, 0) = \sum_{n=1}^\infty \omega_n B_n \sin\left(\frac{n\pi x}{\ell}\right) = \frac{4V}{\ell^2}x(\ell - x).$$
Projecting onto $\sin(n\pi x/\ell)$ gives
$$\omega_n B_n = \frac{2}{\ell}\int_0^\ell \frac{4V}{\ell^2}x(\ell - x)\sin\left(\frac{n\pi x}{\ell}\right)\,dx.$$
Substituting $u = x/\ell$ ($dx = \ell\,du$, $x(\ell - x) = \ell^2(u - u^2)$):
$$\omega_n B_n = 8V \int_0^1 (u - u^2)\sin(n\pi u)\,du.$$
Integrate by parts twice (noting that $u - u^2 = 0$ at $u = 0, 1$ and $\sin(n\pi u) = 0$ at $u = 0, 1$):
$$\begin{aligned}
\int_0^1 (u - u^2)\sin(n\pi u)\,du
&= \underbrace{\left[-(u - u^2)\frac{\cos(n\pi u)}{n\pi}\right]_0^1}_{=\,0} + \frac{1}{n\pi}\int_0^1 (1 - 2u)\cos(n\pi u)\,du \\
&= \frac{1}{n\pi}\left(\underbrace{\left[(1 - 2u)\frac{\sin(n\pi u)}{n\pi}\right]_0^1}_{=\,0} + \frac{2}{n\pi}\int_0^1 \sin(n\pi u)\,du\right) \\
&= \frac{2}{(n\pi)^2}\left[-\frac{\cos(n\pi u)}{n\pi}\right]_0^1
= \frac{2\bigl(1 - (-1)^n\bigr)}{n^3\pi^3}.
\end{aligned}$$
Therefore the velocity-mode coefficients are
$$\omega_n B_n = \frac{16V\bigl(1 - (-1)^n\bigr)}{n^3\pi^3} = \begin{cases} \dfrac{32V}{n^3\pi^3}, & n \text{ odd}, \\[8pt] 0, & n \text{ even} \end{cases}$$
(the vanishing of all even modes reflects the symmetry of $y_t(x, 0)$ about $x = \ell/2$). Dividing by $\omega_n = \frac{n\pi c}{\ell}$ gives the normal-mode displacement amplitudes:
$$A_n = 0 \quad (\forall n \ge 1), \qquad B_n = \begin{cases} \dfrac{32V\ell}{n^4\pi^4 c} = \dfrac{32V\ell}{n^4\pi^4}\sqrt{\dfrac{\mu}{T}}, & n \text{ odd}, \\[10pt] 0, & n \text{ even}. \end{cases}$$
Thus the motion of the string is
$$y(x, t) = \frac{32V\ell}{\pi^4 c}\sum_{n\text{ odd}} \frac{1}{n^4}\sin\left(\frac{n\pi c t}{\ell}\right)\sin\left(\frac{n\pi x}{\ell}\right).$$

### 3. Kinetic and potential energies at time $t$

Since $B_n = O(n^{-4})$, the term-by-term differentiated series for $y_t$ and $y_x$ have coefficients of order $O(n^{-3})$, so both converge absolutely and uniformly for all $(x, t)$:
$$y_t(x, t) = \sum_{n\text{ odd}} \omega_n B_n \cos(\omega_n t)\sin\left(\frac{n\pi x}{\ell}\right), \qquad y_x(x, t) = \sum_{n\text{ odd}} \frac{n\pi}{\ell}B_n \sin(\omega_n t)\cos\left(\frac{n\pi x}{\ell}\right).$$

Using the orthogonality relations on $[0, \ell]$ for $m, n \ge 1$:
$$\int_0^\ell \sin\left(\frac{m\pi x}{\ell}\right)\sin\left(\frac{n\pi x}{\ell}\right)\,dx = \frac{\ell}{2}\delta_{mn}, \qquad \int_0^\ell \cos\left(\frac{m\pi x}{\ell}\right)\cos\left(\frac{n\pi x}{\ell}\right)\,dx = \frac{\ell}{2}\delta_{mn},$$
we evaluate the kinetic energy $K(t)$ and potential energy $U(t)$ (for small slopes, stretching work is $T\int_0^\ell (\sqrt{1+y_x^2}-1)\,dx \approx \frac{1}{2}T\int_0^\ell y_x^2\,dx$):

- **Kinetic energy:**
  $$\begin{aligned}
  K(t) &= \frac{1}{2}\mu \int_0^\ell \left(\frac{\partial y}{\partial t}\right)^2 dx
  = \frac{1}{2}\mu \cdot \frac{\ell}{2}\sum_{n\text{ odd}} (\omega_n B_n)^2 \cos^2(\omega_n t) \\
  &= \frac{\mu\ell}{4}\sum_{n\text{ odd}} \left(\frac{32V}{n^3\pi^3}\right)^2 \cos^2\left(\frac{n\pi c t}{\ell}\right)
  = \frac{256\mu\ell V^2}{\pi^6}\sum_{n\text{ odd}} \frac{1}{n^6}\cos^2\left(\frac{n\pi c t}{\ell}\right).
  \end{aligned}$$

- **Potential energy:**
  Since $T = \mu c^2$ and $\frac{n\pi}{\ell}c = \omega_n$, we have $T\left(\frac{n\pi}{\ell}B_n\right)^2 = \mu(\omega_n B_n)^2$. Hence
  $$\begin{aligned}
  U(t) &= \frac{1}{2}T \int_0^\ell \left(\frac{\partial y}{\partial x}\right)^2 dx
  = \frac{1}{2}T \cdot \frac{\ell}{2}\sum_{n\text{ odd}} \left(\frac{n\pi}{\ell}B_n\right)^2 \sin^2(\omega_n t) \\
  &= \frac{\mu\ell}{4}\sum_{n\text{ odd}} (\omega_n B_n)^2 \sin^2(\omega_n t)
  = \frac{256\mu\ell V^2}{\pi^6}\sum_{n\text{ odd}} \frac{1}{n^6}\sin^2\left(\frac{n\pi c t}{\ell}\right),
  \end{aligned}$$
  where $c = \sqrt{T/\mu}$.

- **Total energy:**
  Adding $K(t)$ and $U(t)$ and using $\cos^2(\omega_n t) + \sin^2(\omega_n t) = 1$ for every $n$:
  $$E(t) = K(t) + U(t) = \frac{256\mu\ell V^2}{\pi^6}\sum_{n\text{ odd}} \frac{1}{n^6},$$
  which is constant in time (energy conservation).

### 4. Comparison with the initial energy and evaluation of the series

At $t = 0$, $y(x, 0) = 0$ implies $y_x(x, 0) = 0$, so $U(0) = 0$. The initial energy is therefore purely kinetic:
$$E(0) = K(0) = \frac{1}{2}\mu\int_0^\ell [y_t(x, 0)]^2\,dx = \frac{1}{2}\mu \cdot \frac{16V^2}{\ell^4}\int_0^\ell x^2(\ell - x)^2\,dx.$$
Setting $u = x/\ell$:
$$\int_0^\ell x^2(\ell - x)^2\,dx = \ell^5 \int_0^1 (u^2 - 2u^3 + u^4)\,du = \ell^5\left[\frac{u^3}{3} - \frac{u^4}{2} + \frac{u^5}{5}\right]_0^1 = \ell^5\left(\frac{1}{3} - \frac{1}{2} + \frac{1}{5}\right) = \frac{\ell^5}{30}.$$
Hence
$$E(0) = 8\mu\frac{V^2}{\ell^4}\cdot\frac{\ell^5}{30} = \frac{4}{15}\mu\ell V^2.$$

Equating $E(t) = E(0)$ (equivalently Parseval's theorem applied to $y_t(x, 0)$ at $t = 0$):
$$\frac{256\mu\ell V^2}{\pi^6}\sum_{n\text{ odd}} \frac{1}{n^6} = \frac{4}{15}\mu\ell V^2.$$
Dividing through by $\frac{256\mu\ell V^2}{\pi^6}$ (for non-zero $\mu, \ell, V$) yields
$$\sum_{n\text{ odd}} \frac{1}{n^6} = \frac{4\pi^6}{15 \times 256} = \frac{\pi^6}{15 \times 64} = \frac{\pi^6}{960}.$$

---

## Final Answers

- **Normal-mode amplitudes** (in $y(x,t) = \sum_{n=1}^\infty [A_n\cos(\omega_n t) + B_n\sin(\omega_n t)]\sin(n\pi x/\ell)$ with $\omega_n = n\pi c/\ell$, $c = \sqrt{T/\mu}$):
  $$A_n = 0 \quad (\forall n \ge 1), \qquad B_n = \begin{cases} \dfrac{32V\ell}{n^4\pi^4 c} = \dfrac{32V\ell}{n^4\pi^4}\sqrt{\dfrac{\mu}{T}}, & n \text{ odd}, \\[10pt] 0, & n \text{ even}. \end{cases}$$
- **Kinetic energy at time $t$:**
  $$K(t) = \frac{256\mu\ell V^2}{\pi^6}\sum_{n\text{ odd}} \frac{1}{n^6}\cos^2\left(\frac{n\pi c t}{\ell}\right).$$
- **Potential energy at time $t$:**
  $$U(t) = \frac{256\mu\ell V^2}{\pi^6}\sum_{n\text{ odd}} \frac{1}{n^6}\sin^2\left(\frac{n\pi c t}{\ell}\right).$$
- **Total/initial energy:** $E(t) = E(0) = \dfrac{4}{15}\mu\ell V^2$, giving $\displaystyle\sum_{n\text{ odd}} \frac{1}{n^6} = \frac{\pi^6}{960}$.

---

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-1_qx12_check.py` (exit code 0).

```text
Exact rational check: I=1/30, E0_coeff=4/15, prefactor=256, ratio=1/960
Numerical sum_{n odd} 1/n^6 = 1.001447076640942, pi^6/960 = 1.001447076640942, err = 0.00e+00
Params (mu=1.0, T=1.0, l=1.0, V=1.0): modes, PDE (res=8.33e-07), K(t), U(t), E(t)=E(0)=0.266667 verified.
Params (mu=2.3, T=7.5, l=1.8, V=3.1): modes, PDE (res=3.35e-06), K(t), U(t), E(t)=E(0)=10.609440 verified.
Params (mu=0.4, T=12.0, l=2.5, V=0.7): modes, PDE (res=3.66e-06), K(t), U(t), E(t)=E(0)=0.130667 verified.
ALL CHECKS PASSED for IB/methods/sheet-1/qx12.
```

---

## Confidence & Edge-Case Notes

- **Confidence:** high.
- **Phase B watchpoints:**
  1. Check the exact factor of $\omega_n = n\pi c/\ell$ when passing from velocity-mode coefficients $\omega_n B_n = \frac{32V}{n^3\pi^3}$ to displacement-mode amplitudes $B_n = \frac{32V\ell}{n^4\pi^4 c}$ (common error: forgetting to divide by $\omega_n$, or dropping $\ell/c$).
  2. Check the $1/4$ prefactor in $K(t)$ and $U(t)$: $\frac{1}{2}\mu$ from kinetic energy times $\frac{\ell}{2}$ from $\int_0^\ell \sin^2(n\pi x/\ell)\,dx = \ell/2$ gives $\frac{\mu\ell}{4}$, and $32^2/4 = 256$.
  3. Check that term-by-term differentiation/integration for $K(t), U(t)$ is valid (here $B_n = O(n^{-4})$, so $y_t, y_x$ series converge uniformly) and that $\cos(m\pi x/\ell)$ orthogonality for $m, n \ge 1$ also has weight $\ell/2$ (no $n = 0$ mode).
