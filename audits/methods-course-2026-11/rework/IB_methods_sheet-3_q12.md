---
id: IB/methods/sheet-3/q12
slug: IB_methods_sheet-3_q12
source_archive: private/raw/sheets/IB/methods/2025-26/sheet-3.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-3_q12_p2.png` (Question 12 at the bottom of page 2 of Michaelmas 2025 Methods Example Sheet 3, Prof. A. C. L. Ashton; verified that page 3 does not exist).

**Verdict: `ok`.** Every word, symbol, derivative order, domain (`t > 0`), initial condition (`y(0) = \cdots = y^{(n-1)}(0) = 0.`), and even the printed sheet's exact phrasing ("for a $n$-times continuously differentiable functions.") match verbatim.

## Full Solution

Let $n \ge 1$ be a positive integer.

### 1. Construction of the Green's function and solution of the IVP

We consider the $n$-th order initial value problem
\[
y^{(n)}(t) = f(t), \qquad t > 0, \qquad y(0) = y'(0) = \cdots = y^{(n-1)}(0) = 0,
\]
where $f$ is continuous on $[0,\infty)$.

For fixed $\tau > 0$, we seek the Green's function $G(t,\tau)$ satisfying
\[
\frac{\partial^n G}{\partial t^n}(t,\tau) = \delta(t - \tau), \qquad t > 0, \qquad \left.\frac{\partial^k G}{\partial t^k}(t,\tau)\right|_{t=0} = 0 \quad (k = 0, 1, \dots, n-1).
\]

#### Region $0 \le t < \tau$
For $0 < t < \tau$, $\frac{\partial^n G}{\partial t^n} = 0$, so $G(t,\tau)$ is a polynomial in $t$ of degree at most $n-1$:
\[
G(t,\tau) = \sum_{k=0}^{n-1} c_k t^k.
\]
The initial conditions at $t = 0$ give $k!\,c_k = \left.\frac{\partial^k G}{\partial t^k}\right|_{t=0} = 0$ for all $k = 0, 1, \dots, n-1$, so $c_k = 0$ for all $k$. Thus
\[
G(t,\tau) = 0 \qquad \text{for } 0 \le t < \tau.
\]

#### Continuity and jump conditions at $t = \tau$
Integrating $\frac{\partial^n G}{\partial t^n} = \delta(t-\tau)$ across $t = \tau$ from $\tau - \varepsilon$ to $\tau + \varepsilon$ and letting $\varepsilon \to 0^+$ gives a unit jump in the $(n-1)$-th derivative:
\[
\left[\frac{\partial^{n-1} G}{\partial t^{n-1}}(t,\tau)\right]_{t=\tau^-}^{t=\tau^+} = 1.
\]
Since $\frac{\partial^{n-1} G}{\partial t^{n-1}}$ is piecewise smooth and bounded near $t = \tau$ (and no derivatives of $\delta(t-\tau)$ appear on the right-hand side), all lower derivatives (when $n \ge 2$) are continuous at $t = \tau$:
\[
\left[\frac{\partial^k G}{\partial t^k}(t,\tau)\right]_{t=\tau^-}^{t=\tau^+} = 0 \qquad (k = 0, 1, \dots, n-2).
\]
(For $n = 1$, there are no lower derivatives, only the unit jump $[G]_{\tau^-}^{\tau^+} = 1$.)

Since all left-hand limits at $t = \tau^-$ vanish, the right-hand limits at $t = \tau^+$ are
\[
\left.\frac{\partial^k G}{\partial t^k}(t,\tau)\right|_{t=\tau^+} = \begin{cases}
0, & k = 0, 1, \dots, n-2, \\
1, & k = n-1.
\end{cases}
\]

#### Region $t > \tau$
For $t > \tau$, $\frac{\partial^n G}{\partial t^n} = 0$, so $G(t,\tau)$ is a polynomial of degree at most $n-1$, which we expand in powers of $(t-\tau)$:
\[
G(t,\tau) = \sum_{k=0}^{n-1} a_k (t - \tau)^k.
\]
Since $\left.\frac{\partial^m G}{\partial t^m}(t,\tau)\right|_{t=\tau^+} = m!\, a_m$ for $0 \le m \le n-1$, matching with the jump and continuity conditions yields $a_0 = a_1 = \cdots = a_{n-2} = 0$ and $(n-1)!\,a_{n-1} = 1$, i.e. $a_{n-1} = \frac{1}{(n-1)!}$.

Combining both regions, for every integer $n \ge 1$ (with the convention $0! = 1$ and $(t-\tau)^0 = 1$ for $n = 1$),
\[
G(t,\tau) = \begin{cases}
0, & 0 \le t < \tau, \\
\dfrac{(t - \tau)^{n-1}}{(n-1)!}, & t > \tau,
\end{cases}
= \frac{(t - \tau)^{n-1}}{(n-1)!} H(t - \tau),
\]
where $H$ is the Heaviside step function.

#### Solution and rigorous verification of existence & uniqueness
Using the Green's function, the solution to the IVP for $t \ge 0$ is
\[
y(t) = \int_0^\infty G(t,\tau) f(\tau)\,d\tau = \frac{1}{(n-1)!}\int_0^t (t - \tau)^{n-1} f(\tau)\,d\tau.
\]

To make the subsequent proof of Taylor's theorem completely rigorous, we verify directly that this function is the **unique** solution for any continuous $f$:
1. **Existence (verification):** By induction on $k = 0, 1, \dots, n-1$, we claim
   \[
   y^{(k)}(t) = \frac{1}{(n-1-k)!}\int_0^t (t - \tau)^{n-1-k} f(\tau)\,d\tau.
   \]
   The base case $k = 0$ is the definition of $y(t)$. If the formula holds for some $0 \le k \le n-2$, then $n - 1 - k \ge 1$, so the integrand $(t-\tau)^{n-1-k} f(\tau)$ vanishes at $\tau = t$; by the Leibniz integral rule,
   \[
   y^{(k+1)}(t) = \frac{1}{(n-1-k)!}\left[0 + \int_0^t (n-1-k)(t-\tau)^{n-2-k} f(\tau)\,d\tau\right] = \frac{1}{(n-2-k)!}\int_0^t (t-\tau)^{n-2-k} f(\tau)\,d\tau,
   \]
   completing the induction. Setting $t = 0$ yields $y^{(k)}(0) = 0$ for all $k = 0, 1, \dots, n-1$. Setting $k = n-1$ gives $y^{(n-1)}(t) = \int_0^t f(\tau)\,d\tau$, and differentiating once more by the Fundamental Theorem of Calculus (using continuity of $f$) gives $y^{(n)}(t) = f(t)$ for $t > 0$.
2. **Uniqueness:** Suppose $z$ is any $n$-times differentiable function on $[0,\infty)$ solving the same IVP, and put $w = z - y$. Then $w^{(n)}(t) = 0$ for all $t > 0$ and $w(0) = w'(0) = \cdots = w^{(n-1)}(0) = 0$. Applying the Mean Value Theorem successively $n$ times: $w^{(n)} \equiv 0 \implies w^{(n-1)}$ is constant and equals $w^{(n-1)}(0) = 0$; inductively $w^{(n-2)} \equiv 0, \dots, w \equiv 0$. Thus the solution $y$ is unique.

### 2. Proof of Taylor's theorem with integral remainder

**Theorem.** *Let $g$ be $n$-times continuously differentiable ($g \in C^n(I)$) on an interval $I$, and let $a, x \in I$. Then*
\[
g(x) = \sum_{k=0}^{n-1} \frac{g^{(k)}(a)}{k!}(x - a)^k + \frac{1}{(n-1)!}\int_a^x (x - s)^{n-1} g^{(n)}(s)\,ds.
\]

**Proof.** Define the remainder function (with shift $t = x - a$) by
\[
y(t) := g(a + t) - \sum_{k=0}^{n-1} \frac{g^{(k)}(a)}{k!} t^k,
\]
for all $t$ such that $a + t \in I$. For each $0 \le m \le n-1$, differentiating $m$ times with respect to $t$ gives
\[
y^{(m)}(t) = g^{(m)}(a + t) - \sum_{k=m}^{n-1} \frac{g^{(k)}(a)}{(k - m)!} t^{k - m}.
\]
Evaluating at $t = 0$, only the $k = m$ term survives in the sum, so
\[
y^{(m)}(0) = g^{(m)}(a) - g^{(m)}(a) = 0 \qquad (m = 0, 1, \dots, n-1).
\]
Differentiating $n$ times annihilates the polynomial of degree at most $n-1$, yielding
\[
y^{(n)}(t) = g^{(n)}(a + t) =: f(t).
\]
Because $g \in C^n(I)$, $f$ is continuous.

- **Case $x > a$ ($t = x - a > 0$):** By the unique solution to the IVP proved in Part 1,
  \[
  y(t) = \frac{1}{(n-1)!}\int_0^t (t - \tau)^{n-1} g^{(n)}(a + \tau)\,d\tau.
  \]
- **Case $x = a$ ($t = 0$):** Both sides vanish identically ($0 = 0$).
- **Case $x < a$ ($t = x - a < 0$):** Define $h(u) := y(-u)$ for $u \ge 0$ with $a - u \in I$. Then $h^{(m)}(0) = (-1)^m y^{(m)}(0) = 0$ for $0 \le m \le n-1$, and $h^{(n)}(u) = (-1)^n y^{(n)}(-u) = (-1)^n g^{(n)}(a - u)$. Applying Part 1 at $u = -t > 0$ gives
  \[
  y(t) = h(-t) = \frac{(-1)^n}{(n-1)!}\int_0^{-t} (-t - \sigma)^{n-1} g^{(n)}(a - \sigma)\,d\sigma.
  \]
  Substituting $\tau = -\sigma$ ($d\sigma = -d\tau$, limits $\sigma: 0 \to -t$ become $\tau: 0 \to t$) and noting $(-t - \sigma)^{n-1} = (-1)^{n-1}(t - \tau)^{n-1}$ and $(-1)^n(-1)^{n-1}(-1) = (-1)^{2n} = 1$, we obtain the exact same formula:
  \[
  y(t) = \frac{1}{(n-1)!}\int_0^t (t - \tau)^{n-1} g^{(n)}(a + \tau)\,d\tau.
  \]

Finally, setting $t = x - a$ and substituting $s = a + \tau$ ($ds = d\tau$, $t - \tau = x - s$, limits $\tau: 0 \to x-a$ become $s: a \to x$) yields
\[
g(x) - \sum_{k=0}^{n-1} \frac{g^{(k)}(a)}{k!}(x - a)^k = \frac{1}{(n-1)!}\int_a^x (x - s)^{n-1} g^{(n)}(s)\,ds,
\]
and rearranging gives Taylor's theorem with integral remainder:
\[
g(x) = \sum_{k=0}^{n-1} \frac{g^{(k)}(a)}{k!}(x - a)^k + \frac{1}{(n-1)!}\int_a^x (x - s)^{n-1} g^{(n)}(s)\,ds. \qquad \square
\]
(In particular, for $a = 0$, $g(x) = \sum_{k=0}^{n-1} \frac{g^{(k)}(0)}{k!}x^k + \frac{1}{(n-1)!}\int_0^x (x-s)^{n-1} g^{(n)}(s)\,ds$.)

## Final Answers

- **Green's function** (for $n \ge 1$):
  \[
  G(t,\tau) = \begin{cases}
  0, & 0 \le t < \tau, \\
  \dfrac{(t - \tau)^{n-1}}{(n-1)!}, & t > \tau.
  \end{cases}
  \]
- **IVP solution** (for $t \ge 0$):
  \[
  y(t) = \frac{1}{(n-1)!}\int_0^t (t - \tau)^{n-1} f(\tau)\,d\tau.
  \]
- **Taylor's theorem with integral remainder** (for $g \in C^n(I)$ and $a, x \in I$):
  \[
  g(x) = \sum_{k=0}^{n-1} \frac{g^{(k)}(a)}{k!}(x - a)^k + \frac{1}{(n-1)!}\int_a^x (x - s)^{n-1} g^{(n)}(s)\,ds.
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-3_q12_check.py`

Verbatim stdout:
```text
IB_methods_sheet-3_q12_check.py: ALL CHECKS PASSED
  exact rational monomial IVP checks: 42/42 exact
  max IVP quadrature error (f=e^t):   1.776e-15
  max Taylor integral-remainder err:  1.421e-14
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Phase B watchpoints:**
  - **Index consistency with $C^n$:** Since $g$ is $n$-times continuously differentiable ($g \in C^n$), the Taylor polynomial must stop at degree $n-1$ ($\sum_{k=0}^{n-1}$) and the integral remainder must be $\frac{1}{(n-1)!}\int_a^x (x-s)^{n-1} g^{(n)}(s)\,ds$ (NOT $\sum_{k=0}^n$ with $g^{(n+1)}$, which would require $g \in C^{n+1}$).
  - **Uniqueness step:** To deduce that the remainder $y(t) = g(a+t) - \sum_{k=0}^{n-1}\frac{g^{(k)}(a)}{k!}t^k$ *equals* $\frac{1}{(n-1)!}\int_0^t (t-\tau)^{n-1} g^{(n)}(a+\tau)\,d\tau$, one needs uniqueness of solutions to the IVP (proved cleanly via the Mean Value Theorem).
  - **Edge cases:** $n = 1$ gives $G(t,\tau) = H(t-\tau)$ and reduces Taylor's theorem to the Fundamental Theorem of Calculus $g(x) = g(a) + \int_a^x g'(s)\,ds$; handling $x < a$ (in addition to $x > a$) via reflection $h(u) = y(-u)$ makes the theorem hold on an arbitrary interval $I$.
