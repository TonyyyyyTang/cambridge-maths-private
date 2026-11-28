---
id: IB/methods/sheet-4/q09
slug: IB_methods_sheet-4_q09
source_archive: private/raw/sheets/IB/methods/2015-16/all-sheets.pdf
pages_checked: [8]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against the 200 DPI rendered page `IB_methods_sheet-4_q09_p8.png` (`private/raw/sheets/IB/methods/2015-16/all-sheets.pdf`, page 8, Question 10). Question 10 sits in the middle of page 8 and is followed by Question 11 on the same page, so there is no page spillover.

- Title: printed as underlined italic *Dirichlet Green's function for the sphere*$^*$. (stored as `\emph{Dirichlet Green’s function for the sphere}$^*$.` — purely typographic omission of underline in restricted LaTeX).
- Part `i)`: exact wording including bold **interior**, displayed formula
  \[
  G(\mathbf{x};\mathbf{x}_0) = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{a}{|\mathbf{x}_0|}\frac{1}{4\pi|\mathbf{x}-\mathbf{x}_0^*|}, \qquad \text{where} \qquad \mathbf{x}_0^* = \frac{a^2\mathbf{x}_0}{|\mathbf{x}_0|^2}.
  \]
  Every sign, fraction, bold vector, subscript, superscript, and punctuation mark (comma after first display formula, period at end) matches verbatim.
- Part `ii)`: exact wording including bold **exterior**: `Derive the Dirichlet Green’s function for the Laplacian for the \textbf{exterior} of a spherical domain of radius $a$.`

Verdict: **`transcription: ok`** (exact match; only purely typographic omission of title underline).

## Full Solution

We work in $\mathbb{R}^3$ with the sphere of radius $a > 0$ centred at the origin, $S_a = \{\mathbf{x} \in \mathbb{R}^3 : |\mathbf{x}| = a\}$. Consistent with the singular term $\frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|}$ given in part (i), the Dirichlet Green's function $G(\mathbf{x};\mathbf{x}_0)$ for the Laplacian on a domain $\Omega \subset \mathbb{R}^3$ (for fixed source $\mathbf{x}_0 \in \Omega$) satisfies
\[
\begin{cases}
\nabla_{\mathbf{x}}^2 G(\mathbf{x};\mathbf{x}_0) = \delta(\mathbf{x} - \mathbf{x}_0), & \mathbf{x} \in \Omega, \\
G(\mathbf{x};\mathbf{x}_0) = 0, & \mathbf{x} \in \partial\Omega,
\end{cases}
\]
together with $G(\mathbf{x};\mathbf{x}_0) \to 0$ as $|\mathbf{x}| \to \infty$ when $\Omega$ is unbounded (part (ii)).

Recall that the free-space fundamental solution in $\mathbb{R}^3$ is $\Phi(\mathbf{x} - \mathbf{p}) = \frac{-1}{4\pi|\mathbf{x} - \mathbf{p}|}$, which is smooth and harmonic for $\mathbf{x} \neq \mathbf{p}$ and satisfies $\nabla_{\mathbf{x}}^2 \Phi(\mathbf{x} - \mathbf{p}) = \delta(\mathbf{x} - \mathbf{p})$ in $\mathbb{R}^3$ (since by the divergence theorem on $B_\varepsilon(\mathbf{p})$, $\int_{\partial B_\varepsilon(\mathbf{p})} \frac{\partial \Phi}{\partial r}\,\mathrm{d}S = \frac{1}{4\pi\varepsilon^2}\cdot 4\pi\varepsilon^2 = 1$).

---

### Part (i): Interior of the sphere ($B_a = \{\mathbf{x} \in \mathbb{R}^3 : |\mathbf{x}| < a\}$)

Let $\mathbf{x}_0 \in B_a$ with $0 < |\mathbf{x}_0| < a$ (we handle the removable special case $\mathbf{x}_0 = \mathbf{0}$ explicitly below). Consider
\[
G(\mathbf{x};\mathbf{x}_0) = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{a}{|\mathbf{x}_0|}\frac{1}{4\pi|\mathbf{x}-\mathbf{x}_0^*|}, \qquad \mathbf{x}_0^* = \frac{a^2 \mathbf{x}_0}{|\mathbf{x}_0|^2}.
\]

1. **Location of the image point and PDE in $B_a$**:
   The distance of the image point $\mathbf{x}_0^*$ (the inversion of $\mathbf{x}_0$ in the sphere $S_a$) from the origin is
   \[
   |\mathbf{x}_0^*| = \frac{a^2 |\mathbf{x}_0|}{|\mathbf{x}_0|^2} = \frac{a^2}{|\mathbf{x}_0|} > a,
   \]
   because $0 < |\mathbf{x}_0| < a$. Hence $\mathbf{x}_0^* \notin \overline{B_a}$: for every $\mathbf{x} \in \overline{B_a}$, $|\mathbf{x} - \mathbf{x}_0^*| \ge |\mathbf{x}_0^*| - |\mathbf{x}| \ge \frac{a^2}{|\mathbf{x}_0|} - a > 0$. Therefore the image term $\frac{a}{|\mathbf{x}_0|}\frac{1}{4\pi|\mathbf{x}-\mathbf{x}_0^*|}$ is smooth and harmonic on an open neighbourhood of $\overline{B_a}$:
   \[
   \nabla_{\mathbf{x}}^2\left(\frac{a}{|\mathbf{x}_0|}\frac{1}{4\pi|\mathbf{x}-\mathbf{x}_0^*|}\right) = 0 \qquad \text{for all } \mathbf{x} \in B_a.
   \]
   By linearity,
   \[
   \nabla_{\mathbf{x}}^2 G(\mathbf{x};\mathbf{x}_0) = \nabla_{\mathbf{x}}^2\left(\frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|}\right) + 0 = \delta(\mathbf{x} - \mathbf{x}_0) \qquad \text{for } \mathbf{x} \in B_a.
   \]

2. **Homogeneous Dirichlet boundary condition on $|\mathbf{x}| = a$**:
   For any $\mathbf{x} \in S_a$ (so $|\mathbf{x}|^2 = a^2$), compute $|\mathbf{x} - \mathbf{x}_0^*|^2$:
   \[
   \begin{aligned}
   |\mathbf{x} - \mathbf{x}_0^*|^2
   &= |\mathbf{x}|^2 - 2\mathbf{x}\cdot\mathbf{x}_0^* + |\mathbf{x}_0^*|^2 \\
   &= a^2 - \frac{2a^2}{|\mathbf{x}_0|^2}(\mathbf{x}\cdot\mathbf{x}_0) + \frac{a^4}{|\mathbf{x}_0|^2} \\
   &= \frac{a^2}{|\mathbf{x}_0|^2}\left(|\mathbf{x}_0|^2 - 2\mathbf{x}\cdot\mathbf{x}_0 + a^2\right) \\
   &= \frac{a^2}{|\mathbf{x}_0|^2}\left(|\mathbf{x}_0|^2 - 2\mathbf{x}\cdot\mathbf{x}_0 + |\mathbf{x}|^2\right) \\
   &= \frac{a^2}{|\mathbf{x}_0|^2}|\mathbf{x} - \mathbf{x}_0|^2.
   \end{aligned}
   \]
   Taking positive square roots (valid since $|\mathbf{x}| = a > |\mathbf{x}_0| > 0$, so both distances are strictly positive):
   \[
   |\mathbf{x} - \mathbf{x}_0^*| = \frac{a}{|\mathbf{x}_0|}|\mathbf{x} - \mathbf{x}_0| \qquad \text{for all } |\mathbf{x}| = a.
   \]
   Substituting this identity into $G(\mathbf{x};\mathbf{x}_0)$ on $|\mathbf{x}| = a$ gives
   \[
   G(\mathbf{x};\mathbf{x}_0) = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{a}{|\mathbf{x}_0|}\cdot \frac{1}{4\pi \frac{a}{|\mathbf{x}_0|}|\mathbf{x}-\mathbf{x}_0|} = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{1}{4\pi|\mathbf{x}-\mathbf{x}_0|} = 0.
   \]

3. **Uniqueness**:
   If $G_1$ and $G_2$ are two Dirichlet Green's functions for $B_a$, then $h(\mathbf{x}) = G_1(\mathbf{x};\mathbf{x}_0) - G_2(\mathbf{x};\mathbf{x}_0)$ has no singularity in $B_a$, satisfies $\nabla^2 h = 0$ in $B_a$, and vanishes on $\partial B_a$. By Green's first identity (or the maximum principle),
   \[
   \int_{B_a} |\nabla h|^2\,\mathrm{d}V = \int_{\partial B_a} h\frac{\partial h}{\partial n}\,\mathrm{d}S - \int_{B_a} h\nabla^2 h\,\mathrm{d}V = 0,
   \]
   so $h$ is constant in $B_a$ and hence $h \equiv 0$ by the boundary condition. Thus $G$ is the unique Dirichlet Green's function for $B_a$.

4. **Special case $\mathbf{x}_0 = \mathbf{0}$ (and symmetric form)**:
   For any $\mathbf{x}_0 \neq \mathbf{0}$, bringing $|\mathbf{x}_0|$ inside the square root gives
   \[
   |\mathbf{x}_0|\,|\mathbf{x} - \mathbf{x}_0^*| = \sqrt{|\mathbf{x}|^2|\mathbf{x}_0|^2 - 2a^2 (\mathbf{x}\cdot\mathbf{x}_0) + a^4},
   \]
   so we can write
   \[
   G(\mathbf{x};\mathbf{x}_0) = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{a}{4\pi\sqrt{|\mathbf{x}|^2|\mathbf{x}_0|^2 - 2a^2(\mathbf{x}\cdot\mathbf{x}_0) + a^4}},
   \]
   which is manifestly symmetric in $\mathbf{x}$ and $\mathbf{x}_0$ and has a smooth limit as $\mathbf{x}_0 \to \mathbf{0}$:
   \[
   G(\mathbf{x};\mathbf{0}) = \frac{-1}{4\pi|\mathbf{x}|} + \frac{1}{4\pi a}.
   \]
   One verifies directly that $\nabla_{\mathbf{x}}^2 G(\mathbf{x};\mathbf{0}) = \delta(\mathbf{x})$ in $B_a$ and $G(\mathbf{x};\mathbf{0}) = 0$ on $|\mathbf{x}| = a$.

---

### Part (ii): Exterior of the sphere ($E_a = \{\mathbf{x} \in \mathbb{R}^3 : |\mathbf{x}| > a\}$)

Now let $\mathbf{x}_0 \in E_a$, so $|\mathbf{x}_0| > a > 0$. We seek $G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0)$ satisfying:
1. $\nabla_{\mathbf{x}}^2 G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) = \delta(\mathbf{x} - \mathbf{x}_0)$ for $|\mathbf{x}| > a$,
2. $G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) = 0$ for $|\mathbf{x}| = a$,
3. $G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) \to 0$ as $|\mathbf{x}| \to \infty$.

**Derivation via the method of images**:
By rotational symmetry about the line through $\mathbf{0}$ and $\mathbf{x}_0$, we try placing a single image charge of strength $q$ at a point $\mathbf{x}_0^* = \lambda \mathbf{x}_0$ along the ray from the origin towards $\mathbf{x}_0$, with $\lambda > 0$ chosen so that $\mathbf{x}_0^*$ lies strictly inside the sphere ($|\mathbf{x}_0^*| = \lambda |\mathbf{x}_0| < a$, hence outside $\overline{E_a}$):
\[
G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{q}{4\pi|\mathbf{x}-\mathbf{x}_0^*|}.
\]
Because $\mathbf{x}_0^* \notin \overline{E_a}$, the second term is smooth and harmonic on $E_a$, so $\nabla_{\mathbf{x}}^2 G_{\mathrm{ext}} = \delta(\mathbf{x}-\mathbf{x}_0)$ holds automatically in $E_a$.

Imposing the Dirichlet boundary condition $G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) = 0$ on $|\mathbf{x}| = a$ requires $q > 0$ and
\[
|\mathbf{x} - \mathbf{x}_0^*| = q|\mathbf{x} - \mathbf{x}_0| \qquad \text{for all } \mathbf{x} \text{ with } |\mathbf{x}| = a.
\]
Squaring both sides and using $|\mathbf{x}|^2 = a^2$ and $\mathbf{x}_0^* = \lambda \mathbf{x}_0$:
\[
a^2 - 2\lambda(\mathbf{x}\cdot\mathbf{x}_0) + \lambda^2|\mathbf{x}_0|^2 = q^2\bigl(a^2 - 2(\mathbf{x}\cdot\mathbf{x}_0) + |\mathbf{x}_0|^2\bigr).
\]
Since $\mathbf{x}\cdot\mathbf{x}_0$ varies over the entire interval $[-a|\mathbf{x}_0|, +a|\mathbf{x}_0|]$ as $\mathbf{x}$ ranges over $S_a$, this identity holds on all of $S_a$ if and only if the coefficients of $\mathbf{x}\cdot\mathbf{x}_0$ and the constant terms match separately:
\[
\begin{cases}
\lambda = q^2, \\
a^2 + \lambda^2|\mathbf{x}_0|^2 = q^2(a^2 + |\mathbf{x}_0|^2).
\end{cases}
\]
Substituting $q^2 = \lambda$ into the second equation gives
\[
a^2 + \lambda^2|\mathbf{x}_0|^2 - \lambda a^2 - \lambda|\mathbf{x}_0|^2 = (\lambda - 1)\bigl(\lambda|\mathbf{x}_0|^2 - a^2\bigr) = 0.
\]
The root $\lambda = 1$ gives $\mathbf{x}_0^* = \mathbf{x}_0 \in E_a$ (yielding the trivial zero function, not a Green's function). The unique admissible root is
\[
\lambda = \frac{a^2}{|\mathbf{x}_0|^2}, \qquad q = \sqrt{\lambda} = \frac{a}{|\mathbf{x}_0|}, \qquad \mathbf{x}_0^* = \frac{a^2 \mathbf{x}_0}{|\mathbf{x}_0|^2}.
\]

**Verification of all conditions for $E_a$**:
1. **Image lies outside $\overline{E_a}$**: Since $|\mathbf{x}_0| > a$, we have
   \[
   |\mathbf{x}_0^*| = \frac{a^2}{|\mathbf{x}_0|} < a,
   \]
   so $\mathbf{x}_0^*$ lies strictly inside the open ball $B_a$ (and hence strictly outside $\overline{E_a}$). Therefore $\frac{a}{|\mathbf{x}_0|}\frac{1}{4\pi|\mathbf{x}-\mathbf{x}_0^*|}$ is smooth and harmonic on $E_a$, and $\nabla_{\mathbf{x}}^2 G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) = \delta(\mathbf{x}-\mathbf{x}_0)$ in $E_a$.
2. **Dirichlet boundary condition**: By construction (and by the exact same algebra as in part (i)), $|\mathbf{x}-\mathbf{x}_0^*| = \frac{a}{|\mathbf{x}_0|}|\mathbf{x}-\mathbf{x}_0|$ on $|\mathbf{x}| = a$, so $G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) = 0$ on $S_a$.
3. **Decay at infinity**: For fixed $\mathbf{x}_0 \in E_a$, as $|\mathbf{x}| \to \infty$,
   \[
   G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) = -\frac{1}{4\pi|\mathbf{x}|}\left(1 - \frac{a}{|\mathbf{x}_0|}\right) + O\!\left(\frac{1}{|\mathbf{x}|^2}\right) \longrightarrow 0.
   \]
4. **Uniqueness**: If $G_1$ and $G_2$ both satisfy $\nabla^2 G_i = \delta(\mathbf{x}-\mathbf{x}_0)$ in $E_a$, $G_i = 0$ on $S_a$, and $G_i(\mathbf{x};\mathbf{x}_0) \to 0$ as $|\mathbf{x}| \to \infty$, then $h = G_1 - G_2$ is harmonic on $E_a$, vanishes on $S_a$, and tends to $0$ as $|\mathbf{x}| \to \infty$. Given any $\varepsilon > 0$, choose $R > a$ sufficiently large that $|h(\mathbf{x})| < \varepsilon$ for $|\mathbf{x}| = R$. By the maximum and minimum principles on the bounded spherical shell $a < |\mathbf{x}| < R$, $|h(\mathbf{x})| \le \varepsilon$ throughout $a < |\mathbf{x}| < R$. Letting $R \to \infty$ and $\varepsilon \to 0^+$ yields $h \equiv 0$ on $E_a$. (Note that without the decay condition at infinity, one could add any multiple of $1 - a/|\mathbf{x}|$, which is harmonic on $E_a$ and vanishes on $S_a$.)

Thus the Dirichlet Green's function for the exterior of the sphere of radius $a$ has the **same algebraic expression** as for the interior (now evaluated for $|\mathbf{x}| > a$ and $|\mathbf{x}_0| > a$, with the image point $\mathbf{x}_0^*$ lying inside the sphere rather than outside):
\[
G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{a}{|\mathbf{x}_0|}\frac{1}{4\pi|\mathbf{x}-\mathbf{x}_0^*|}, \qquad \text{where} \qquad \mathbf{x}_0^* = \frac{a^2\mathbf{x}_0}{|\mathbf{x}_0|^2}.
\]

## Final Answers

- **(i) Interior ($|\mathbf{x}| < a$, $|\mathbf{x}_0| < a$)**:
  \[
  G(\mathbf{x};\mathbf{x}_0) = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{a}{|\mathbf{x}_0|}\frac{1}{4\pi|\mathbf{x}-\mathbf{x}_0^*|}, \qquad \mathbf{x}_0^* = \frac{a^2\mathbf{x}_0}{|\mathbf{x}_0|^2} \quad (\text{for } \mathbf{x}_0 \neq \mathbf{0}),
  \]
  with continuous limit $G(\mathbf{x};\mathbf{0}) = \frac{-1}{4\pi|\mathbf{x}|} + \frac{1}{4\pi a}$ at $\mathbf{x}_0 = \mathbf{0}$ (equivalently $G(\mathbf{x};\mathbf{x}_0) = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{a}{4\pi\sqrt{|\mathbf{x}|^2|\mathbf{x}_0|^2 - 2a^2\mathbf{x}\cdot\mathbf{x}_0 + a^4}}$ for all $\mathbf{x}_0 \in B_a$).
- **(ii) Exterior ($|\mathbf{x}| > a$, $|\mathbf{x}_0| > a$)**:
  \[
  G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) = \frac{-1}{4\pi|\mathbf{x}-\mathbf{x}_0|} + \frac{a}{|\mathbf{x}_0|}\frac{1}{4\pi|\mathbf{x}-\mathbf{x}_0^*|}, \qquad \mathbf{x}_0^* = \frac{a^2\mathbf{x}_0}{|\mathbf{x}_0|^2},
  \]
  which has the same algebraic form, with $|\mathbf{x}_0^*| = a^2/|\mathbf{x}_0| < a$ lying inside the sphere and $G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) \to 0$ as $|\mathbf{x}| \to \infty$.

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-4_q09_check.py`

```text
IB_methods_sheet-4_q09 verification PASSED:
  Interior max BC error          = 1.388e-17
  Interior Laplacian residual    = 3.469e-10
  Interior delta flux            = 1.00000000 (target 1.0)
  Interior Poisson err (poly)    = 6.119e-06
  Exterior max BC error          = 6.939e-18
  Exterior Laplacian residual    = 5.204e-10
  Exterior delta flux            = 1.00000000 (target 1.0)
  Exterior far-field rel err     = 5.207e-03
  Exterior Poisson err (a/|x|)   = 1.661e-07
  Exterior Poisson err (1/|x-p|) = 5.129e-08
```

## Confidence & Edge-Case Notes

- **Confidence**: high.
- **Phase B watchpoints**:
  1. **Special point $\mathbf{x}_0 = \mathbf{0}$ in part (i)**: The printed formula has $|\mathbf{x}_0|$ in denominators; check whether the existing solution remarks on the limit $\mathbf{x}_0 \to \mathbf{0}$ giving $G(\mathbf{x};\mathbf{0}) = -\frac{1}{4\pi|\mathbf{x}|} + \frac{1}{4\pi a}$.
  2. **Exterior image location & decay at infinity in part (ii)**: Although the algebraic formula is identical in parts (i) and (ii), part (ii) specifically asks to *derive* $G_{\mathrm{ext}}$, and requires checking BOTH that $|\mathbf{x}_0^*| = a^2/|\mathbf{x}_0| < a$ (so the image charge now lies strictly *inside* the sphere, outside the exterior domain!) AND that $G_{\mathrm{ext}}(\mathbf{x};\mathbf{x}_0) \to 0$ as $|\mathbf{x}| \to \infty$ (without which $1 - a/|\mathbf{x}|$ could be added arbitrarily!).
  3. **Sign of image charge & PDE convention**: Watch for sign errors in $\frac{a}{|\mathbf{x}_0|}$ vs $\frac{|\mathbf{x}_0|}{a}$! A very common mistake in part (ii) is flipping the factor to $\frac{|\mathbf{x}_0|}{a}$ when moving from interior to exterior, whereas in fact $q = \frac{a}{|\mathbf{x}_0|}$ in BOTH cases (for exterior, $q < 1$ because the image charge is smaller than the external source charge!).
  4. **Verification of both $\lambda = 1$ and $\lambda = a^2/|\mathbf{x}_0|^2$**: In deriving part (ii), solving $(\lambda - 1)(\lambda|\mathbf{x}_0|^2 - a^2) = 0$ cleanly rejects $\lambda = 1$ (which lies inside $E_a$ and gives $G \equiv 0$).
