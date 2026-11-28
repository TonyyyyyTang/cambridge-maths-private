---
id: IB/methods/sheet-2/qx19
slug: IB_methods_sheet-2_qx19
source_archive: private/raw/sheets/IB/methods/2007-08/sheet-2.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-2_qx19_p1.png` (Question 3, *Finite asymptotics*, on page 1 of *Mathematical Tripos Part IB, Michaelmas term 2007, Dr C. P. Caulfield, METHODS — EXAMPLES II*). Question 4 follows immediately on page 1, so there is no spillover onto page 2.

Every word, symbol, sign, operator definition, piecewise condition, constant ($A$, $B$, $R$, $a$), and punctuation mark matches the printed PDF exactly (even the omission of a comma after $B$ in "for suitable constants $A, B$"). Purely typographic difference only: underlined italic heading in the PDF vs `\emph{Finite asymptotics.}` in `question.tex`. Verdict: `transcription: ok`.

## Full Solution

Throughout we work on the half-line $x > 0$ (with $x = 0$ understood by continuous extension, and parameters $a > 0$, $R > 0$).

### 1. Homogeneous solutions of $L_x y = 0$

We are given
\[
L_x y = -\frac{1}{x^2}\frac{d}{dx}\!\left(x^2\frac{dy}{dx}\right) + y.
\]
Make the substitution $y(x) = z(x)/x$. Then
\[
\frac{dy}{dx} = \frac{x z' - z}{x^2} \implies x^2\frac{dy}{dx} = x z' - z,
\]
and differentiating once more gives
\[
\frac{d}{dx}\!\left(x^2\frac{dy}{dx}\right) = x z'' + z' - z' = x z''.
\]
Hence
\[
L_x y = -\frac{1}{x^2}(x z'') + \frac{z}{x} = \frac{-z'' + z}{x}.
\]
For $x > 0$, $L_x y = 0$ is therefore equivalent to
\[
z'' - z = 0,
\]
whose general solution is $z(x) = C_1 e^x + C_2 e^{-x} = D_1 \sinh x + D_2 \cosh x$. Thus the general solution of $L_x y = 0$ on $(0, \infty)$ is
\[
y(x) = \frac{C_1 e^x + C_2 e^{-x}}{x} = \frac{D_1 \sinh x + D_2 \cosh x}{x}.
\]

Now impose each asymptotic condition in turn:
- **(a) Bounded as $x \to 0^+$:** As $x \to 0^+$, $\frac{\sinh x}{x} = 1 + \frac{x^2}{6} + O(x^4) \to 1$, whereas $\frac{\cosh x}{x} = \frac{1}{x} + \frac{x}{2} + O(x^3) \to +\infty$. Hence $y(x)$ is bounded as $x \to 0^+$ if and only if $D_2 = 0$, giving
  \[
  y_1(x) = C \frac{\sinh x}{x} \qquad (C \text{ arbitrary constant}).
  \]
- **(b) Bounded as $x \to \infty$:** As $x \to \infty$, $\frac{e^x}{x} \to +\infty$, whereas $\frac{e^{-x}}{x} \to 0$. Hence $y(x)$ is bounded as $x \to \infty$ if and only if $C_1 = 0$, giving
  \[
  y_2(x) = K \frac{e^{-x}}{x} \qquad (K \text{ arbitrary constant}).
  \]

Notice that the only solution of $L_x y = 0$ satisfying *both* (a) and (b) is $y \equiv 0$, so the Green's function exists uniquely.

### 2. Construction of the Green's function $G(x, a)$

We seek $G(x, a)$ ($x > 0$, $a > 0$) satisfying
\[
L_x G(x, a) = -\frac{1}{x^2}\frac{\partial}{\partial x}\!\left(x^2\frac{\partial G}{\partial x}\right) + G = \delta(x - a),
\]
subject to conditions (a) ($G(x,a)$ bounded as $x \to 0^+$) and (b) ($G(x,a)$ bounded as $x \to \infty$).

For $x \neq a$, $L_x G(x, a) = 0$. From Section 1:
- For $0 < x < a$, condition (a) gives $G(x, a) = \alpha(a)\dfrac{\sinh x}{x}$.
- For $x > a$, condition (b) gives $G(x, a) = \beta(a)\dfrac{e^{-x}}{x}$.

Now apply the two matching conditions at $x = a$:

1. **Continuity of $G(x, a)$ at $x = a$:**
   \[
   \alpha(a)\frac{\sinh a}{a} = \beta(a)\frac{e^{-a}}{a} \implies \alpha(a)\sinh a = \beta(a)e^{-a}.
   \]
   Writing $\alpha(a) = C(a)e^{-a}$ and $\beta(a) = C(a)\sinh a$ for some function $C(a)$, we have
   \[
   G(x, a) = \begin{cases}
   C(a) e^{-a}\dfrac{\sinh x}{x}, & 0 < x \le a, \\[6pt]
   C(a) \sinh a\dfrac{e^{-x}}{x}, & x \ge a.
   \end{cases}
   \]

2. **Derivative jump condition at $x = a$:**
   Multiply $L_x G(x,a) = \delta(x - a)$ by $-x^2$ to put it into self-adjoint Sturm–Liouville form:
   \[
   \frac{\partial}{\partial x}\!\left(x^2\frac{\partial G}{\partial x}\right) - x^2 G = -x^2 \delta(x - a) = -a^2 \delta(x - a).
   \]
   Integrating from $x = a - \varepsilon$ to $x = a + \varepsilon$ and letting $\varepsilon \to 0^+$ (using continuity of $G$ at $x = a$, so $\int_{a-\varepsilon}^{a+\varepsilon} x^2 G\,dx \to 0$) yields
   \[
   \left[x^2\frac{\partial G}{\partial x}\right]_{x=a^-}^{x=a^+} = -a^2 \iff \left[\frac{\partial G}{\partial x}\right]_{x=a^-}^{x=a^+} = -1.
   \]
   (Note the crucial minus sign from the leading coefficient $p(x) = -1$ of $y''$ in $L_x y = -y'' - \frac{2}{x}y' + y$.)

   Compute $x^2 \frac{\partial G}{\partial x}$ on either side of $x = a$:
   - For $0 < x < a$:
     \[
     x^2\frac{\partial G}{\partial x} = C(a)e^{-a}(x\cosh x - \sinh x).
     \]
   - For $x > a$:
     \[
     x^2\frac{\partial G}{\partial x} = C(a)\sinh a\,(-x e^{-x} - e^{-x}) = -C(a)\sinh a\,(x + 1)e^{-x}.
     \]
   Taking the jump across $x = a$:
   \[
   \left[x^2\frac{\partial G}{\partial x}\right]_{a^-}^{a^+} = C(a)e^{-a}\bigl[-(a+1)\sinh a - (a\cosh a - \sinh a)\bigr] = -a C(a)e^{-a}(\sinh a + \cosh a) = -a C(a),
   \]
   since $\sinh a + \cosh a = e^a$. Equating this to $-a^2$ gives
   \[
   -a C(a) = -a^2 \implies C(a) = a.
   \]

Therefore the Green's function is
\[
G(x, a) = \begin{cases}
\dfrac{a e^{-a}\sinh x}{x}, & 0 < x \le a, \\[8pt]
\dfrac{a \sinh a\, e^{-x}}{x}, & x \ge a,
\end{cases}
\]
with continuous limit $G(0, a) = a e^{-a}$ at $x = 0$. (As a consistency check, $G(x,a) > 0$ for all $x,a>0$ as expected for the positive operator $-\nabla^2 + 1$, and $x^2 G(x,a) = a^2 G(a,x)$ reflects self-adjointness of $L_x$ with respect to weight $x^2$.)

### 3. Solution of the inhomogeneous problem using $G(x, a)$

We solve $L_x y(x) = f(x)$ subject to (a) and (b), where
\[
f(x) = \begin{cases} 1, & 0 \le x \le R, \\ 0, & x > R. \end{cases}
\]
By the defining property $L_x G(x, a) = \delta(x - a)$,
\[
y(x) = \int_0^\infty G(x, a) f(a)\,da = \int_0^R G(x, a)\,da.
\]
(Note: there is no extra factor of $a^2$ in the integral because $L_x G(x,a) = \delta(x-a)$ already incorporates the factor $a^2$ into $G(x,a)$.)

We evaluate the two elementary integrals needed via integration by parts:
\[
\int_0^u a\sinh a\,da = \bigl[a\cosh a\bigr]_0^u - \int_0^u \cosh a\,da = u\cosh u - \sinh u,
\]
\[
\int_u^R a e^{-a}\,da = \bigl[-a e^{-a}\bigr]_u^R + \int_u^R e^{-a}\,da = \bigl[-(a+1)e^{-a}\bigr]_u^R = (u + 1)e^{-u} - (R + 1)e^{-R}.
\]

**Region 1: $0 < x \le R$ (and $x = 0$ by continuity)**

Splitting the integral at $a = x$ (where $G(x,a)$ changes branch):
\[
y(x) = \int_0^x \frac{a\sinh a\,e^{-x}}{x}\,da + \int_x^R \frac{a e^{-a}\sinh x}{x}\,da
= \frac{e^{-x}}{x}\int_0^x a\sinh a\,da + \frac{\sinh x}{x}\int_x^R a e^{-a}\,da.
\]
Substituting the two integrals:
\[
\begin{aligned}
y(x) &= \frac{e^{-x}}{x}(x\cosh x - \sinh x) + \frac{\sinh x}{x}\bigl[(x + 1)e^{-x} - (R + 1)e^{-R}\bigr] \\
&= \frac{e^{-x}}{x}\bigl[x\cosh x - \sinh x + x\sinh x + \sinh x\bigr] - (R + 1)e^{-R}\frac{\sinh x}{x} \\
&= \frac{e^{-x}}{x}\bigl[x(\cosh x + \sinh x)\bigr] - (R + 1)e^{-R} x^{-1}\sinh x \\
&= 1 - (R + 1)e^{-R} x^{-1}\sinh x,
\end{aligned}
\]
since $\cosh x + \sinh x = e^x$. At $x = 0$, by continuous extension $\lim_{x\to 0^+} x^{-1}\sinh x = 1$, so $y(0) = 1 - (R+1)e^{-R}$, which is finite.

**Region 2: $x > R$**

For all $a \in (0, R)$ we have $x > a$, so only the $x \ge a$ branch of $G(x, a)$ enters:
\[
y(x) = \int_0^R \frac{a\sinh a\,e^{-x}}{x}\,da = (R\cosh R - \sinh R)x^{-1}e^{-x}.
\]

Thus $y(x)$ has the exact required form
\[
y(x) = \begin{cases}
1 + A x^{-1}\sinh x, & 0 \le x \le R, \\
B x^{-1}e^{-x}, & x > R,
\end{cases}
\]
with the explicit constants
\[
A = -(R + 1)e^{-R}, \qquad B = R\cosh R - \sinh R.
\]

**Direct consistency check at $x = R$:**
Both $y$ and $y'$ must be continuous across the finite jump of $f$ at $x = R$. Indeed, letting $z(x) = x y(x)$:
- $z(R^-) = R - (R+1)e^{-R}\sinh R = \frac{1}{2}[(R-1) + (R+1)e^{-2R}] = (R\cosh R - \sinh R)e^{-R} = z(R^+)$, so $y(R^-) = y(R^+)$.
- $z'(R^-) = 1 - (R+1)e^{-R}\cosh R = \frac{1}{2}[1 - R - (R+1)e^{-2R}] = -(R\cosh R - \sinh R)e^{-R} = z'(R^+)$, so $y'(R^-) = y'(R^+)$.

## Final Answers

- **Homogeneous solutions of $L_x y = 0$:**
  - (a) Bounded as $x \to 0$: $y(x) = C x^{-1}\sinh x$ ($C$ constant).
  - (b) Bounded as $x \to \infty$: $y(x) = K x^{-1}e^{-x}$ ($K$ constant).
- **Green's function ($a > 0$):**
  \[
  G(x, a) = \begin{cases}
  a e^{-a} x^{-1}\sinh x, & 0 < x \le a, \\
  a \sinh a\, x^{-1}e^{-x}, & x \ge a.
  \end{cases}
  \]
- **Solution for the step-function RHS:**
  \[
  y(x) = \begin{cases}
  1 + A x^{-1}\sinh x, & 0 \le x \le R, \\
  B x^{-1}e^{-x}, & x > R,
  \end{cases}
  \qquad \text{with} \qquad A = -(R + 1)e^{-R}, \quad B = R\cosh R - \sinh R.
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-2_qx19_check.py`

Verbatim stdout (exit code 0):
```text
IB_methods_sheet-2_qx19_check: ALL PASSED
  homogeneous ODE err: 7.559e-07
  Green's fn: ODE err=3.221e-07, cont err=9.274e-13, jump err=2.594e-06, sym err=2.220e-16
  Green reproduction err (psi=exp(-x^2)): 3.331e-16
  step-RHS: quad vs closed-form err=4.441e-16, ODE err=5.952e-07
  matching at x=R: C^0 err=4.513e-14, C^1 err=1.760e-07, y'' jump(+1) err=3.045e-05
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`. Verified via three independent derivations of $G(x,a)$ (direct jump, Wronskian formula, and reduction to $H(x,a) = x G(x,a)$), exact algebraic $C^1$ matching at $x = R$, and machine-precision numerical quadrature and Green reproduction checks.
- **Phase B watchpoints:**
  - **Major IB Methods pitfall (weight/symmetry & sign of $G(x,a)$):** Because $L_x = -\frac{1}{x^2}\frac{d}{dx}(x^2\frac{d}{dx}) + 1$ has leading coefficient $p(x) = -1$ (negative!) and is NOT yet multiplied by $x^2$, the jump condition is $[\partial G/\partial x]_{a^-}^{a^+} = -1$, equivalently $[x^2 \partial G/\partial x]_{a^-}^{a^+} = -a^2$. Consequently $G(x,a)$ is **not** symmetric in $(x,a)$ — it carries an extra factor of $a$ (specifically $x^2 G(x,a) = a^2 G(a,x)$), and $y(x) = \int_0^R G(x,a)\,da$ has NO extra $a^2$ weight. Watch carefully for solutions that erroneously drop $a^2$ in the jump condition or flip the sign of $G(x,a)$!
  - Check that the constants $A = -(R+1)e^{-R}$ and $B = R\cosh R - \sinh R$ are found explicitly (not merely left implicit) and that the cancellation of $e^{-x}(\cosh x + \sinh x) = 1$ in $0 \le x \le R$ is derived cleanly from the integral without fudging signs.
