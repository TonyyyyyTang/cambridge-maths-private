---
id: IB/methods/sheet-2/qx10
slug: IB_methods_sheet-2_qx10
source_archive: private/raw/sheets/IB/methods/2026-27/sheet-2.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word-for-word and symbol-for-symbol against `IB_methods_sheet-2_qx10_p2.png` (Question 10 on page 2 of *Mathematical Tripos Part IB, Methods, Example Sheet 2*, Michaelmas 2026, Prof. A.C.L. Ashton). Page 2 is the final page of the sheet (no page 3 exists) and Question 10 terminates cleanly with a full stop at the bottom of page 2.

Every word, symbol, equation, boundary/initial condition, assumption (a)–(d), and subpart (i)–(iv) matches the printed PDF exactly. Purely typographic differences only: systems $(A)$ and $(B)$ are printed side-by-side on one line in the PDF versus stacked in consecutive display blocks in `question.tex`, and the printed sheet uses upright `\mathrm{d}V` whereas `question.tex` uses `dV`. Verdict: `transcription: ok`.

## Full Solution

Throughout, $\Omega \subset \mathbf{R}^3$ is a bounded (open, connected) domain with sufficiently smooth boundary $\partial\Omega$ for the divergence theorem (Green's identities) to apply, and each eigenfunction $\psi_n$ is real-valued, smooth enough for Green's identities (e.g. $C^2(\Omega)\cap C^1(\overline{\Omega})$), and non-identically zero in $\Omega$. We assume $\kappa > 0$ and $c > 0$.

### Part (i): Positivity of $\lambda_n$ and orthogonality for distinct eigenvalues

**1. Proof that $\lambda_n > 0$ for each $n$:**

By the product rule for divergence,
\[
\nabla \cdot (\psi_n \nabla\psi_n) = |\nabla\psi_n|^2 + \psi_n \Delta\psi_n.
\]
Integrating over $\Omega$ and applying the divergence theorem (Green's first identity),
\[
\int_{\partial\Omega} \psi_n (\nabla\psi_n \cdot \mathbf{n})\,dS = \int_\Omega |\nabla\psi_n|^2\,dV + \int_\Omega \psi_n \Delta\psi_n\,dV,
\]
where $\mathbf{n}$ is the outward unit normal to $\partial\Omega$. By assumption (b), $\psi_n = 0$ on $\partial\Omega$, so the boundary integral vanishes. Substituting $-\Delta\psi_n = \lambda_n \psi_n$ from assumption (a) gives
\[
\lambda_n \int_\Omega \psi_n^2\,dV = -\int_\Omega \psi_n \Delta\psi_n\,dV = \int_\Omega |\nabla\psi_n|^2\,dV.
\]
Since $\psi_n$ is real, continuous, and not identically zero on $\Omega$, we have $\int_\Omega \psi_n^2\,dV > 0$, while $|\nabla\psi_n|^2 \ge 0$ everywhere in $\Omega$. Therefore
\[
\lambda_n = \frac{\int_\Omega |\nabla\psi_n|^2\,dV}{\int_\Omega \psi_n^2\,dV} \ge 0.
\]
Suppose for contradiction that $\lambda_n = 0$. Then $\int_\Omega |\nabla\psi_n|^2\,dV = 0$. Since $|\nabla\psi_n|^2 \ge 0$ is continuous, this forces $\nabla\psi_n \equiv \mathbf{0}$ everywhere in $\Omega$. Hence $\psi_n$ is constant on $\Omega$ (or on each connected component of $\Omega$). Because $\Omega$ is bounded, its boundary $\partial\Omega$ is non-empty (and meets the closure of every connected component), and $\psi_n = 0$ on $\partial\Omega$. By continuity up to $\partial\Omega$, the constant must be $0$, so $\psi_n \equiv 0$ in $\Omega$, contradicting that $\psi_n$ is an eigenfunction. Thus $\lambda_n > 0$ for every $n$.

**2. Proof that $\int_\Omega \psi_n \psi_m\,dV = 0$ for $\lambda_n \neq \lambda_m$:**

Applying the divergence theorem to $\psi_n \nabla\psi_m - \psi_m \nabla\psi_n$ (Green's second identity),
\[
\int_\Omega (\psi_n \Delta\psi_m - \psi_m \Delta\psi_n)\,dV = \int_{\partial\Omega} \left(\psi_n \frac{\partial\psi_m}{\partial n} - \psi_m \frac{\partial\psi_n}{\partial n}\right)\,dS.
\]
Since $\psi_n = \psi_m = 0$ on $\partial\Omega$, the surface integral vanishes. Using $\Delta\psi_n = -\lambda_n\psi_n$ and $\Delta\psi_m = -\lambda_m\psi_m$, the left-hand side becomes
\[
\int_\Omega \left(\psi_n(-\lambda_m\psi_m) - \psi_m(-\lambda_n\psi_n)\right)\,dV = (\lambda_n - \lambda_m)\int_\Omega \psi_n\psi_m\,dV = 0.
\]
Hence, whenever $\lambda_n \neq \lambda_m$, we must have
\[
\int_\Omega \psi_n \psi_m\,dV = 0.
\]

### Part (ii): Orthogonality for all $n \neq m$ without loss of generality

Define the real $L^2(\Omega)$ inner product
\[
\langle u, v \rangle = \int_\Omega u(\mathbf{x})v(\mathbf{x})\,dV.
\]
This is symmetric, bilinear, and positive definite on real continuous functions on $\Omega$.

For each distinct eigenvalue $\lambda$, let $E_\lambda$ be the corresponding eigenspace, i.e. the set of solutions of $-\Delta\psi = \lambda\psi$ in $\Omega$ with $\psi = 0$ on $\partial\Omega$. By linearity of $-\Delta$ and homogeneity of the Dirichlet boundary condition, $E_\lambda$ is a real vector space: any linear combination of eigenfunctions with eigenvalue $\lambda$ is again an eigenfunction with eigenvalue $\lambda$ (or identically zero).

By assumption (c), each eigenspace $E_\lambda$ has finite dimension $d_\lambda = \dim E_\lambda < \infty$. Let $\{\phi_1, \dots, \phi_{d_\lambda}\}$ be a basis of $E_\lambda$ extracted from the given eigenfunctions with eigenvalue $\lambda$. Applying the Gram–Schmidt orthogonalisation process to $\{\phi_1, \dots, \phi_{d_\lambda}\}$ with respect to the inner product $\langle\cdot,\cdot\rangle$, we obtain a new basis $\{\tilde\psi_1, \dots, \tilde\psi_{d_\lambda}\}$ of $E_\lambda$ satisfying
\[
\langle \tilde\psi_j, \tilde\psi_k \rangle = \int_\Omega \tilde\psi_j \tilde\psi_k\,dV = 0 \qquad (j \neq k).
\]
Replacing the eigenfunctions in each eigenspace $E_\lambda$ by this orthogonal basis preserves all four properties (a)–(d):
- (a) and (b) hold because each $\tilde\psi_j$ is a linear combination of elements of $E_\lambda$;
- (c) still holds since $\dim E_\lambda = d_\lambda$ is unchanged;
- (d) completeness holds because each original eigenfunction in $E_\lambda$ is a finite linear combination of the new basis of $E_\lambda$, so the span of the system is unchanged.

Finally, by part (i), any two eigenfunctions belonging to *different* eigenspaces ($\lambda_n \neq \lambda_m$) are automatically orthogonal. Combining orthogonality within each eigenspace and between distinct eigenspaces, we may assume without loss of generality that
\[
\int_\Omega \psi_n \psi_m\,dV = 0 \qquad \text{for all } n \neq m.
\]

### Part (iii): Separation of variables for the heat problem $(A)$ and formal exponential operator

**1. Solution of $(A)$ by separation of variables:**

Seek non-trivial separable solutions $\varphi(\mathbf{x}, t) = X(\mathbf{x})T(t)$ of $\varphi_t - \kappa\Delta\varphi = 0$ in $\Omega\times(0,\infty)$ with $\varphi = 0$ on $\partial\Omega\times[0,\infty)$. Substituting into the PDE and dividing by $\kappa X(\mathbf{x})T(t)$ yields
\[
\frac{T'(t)}{\kappa T(t)} = \frac{\Delta X(\mathbf{x})}{X(\mathbf{x})} = -\lambda,
\]
where $\lambda$ must be constant since the left-hand side depends only on $t$ and the right-hand side depends only on $\mathbf{x}$.

The spatial part and boundary condition ($T \not\equiv 0$) yield the Dirichlet eigenvalue problem
\[
-\Delta X = \lambda X \quad \text{in } \Omega, \qquad X = 0 \quad \text{on } \partial\Omega.
\]
Its eigenpairs are $\{(\psi_n, \lambda_n)\}_{n=1}^\infty$ (no other eigenvalue can exist, because by part (i) any eigenfunction for $\lambda \notin \{\lambda_n\}$ would be orthogonal to every $\psi_n$, forcing it to vanish identically by completeness (d)). For $\lambda = \lambda_n$, the temporal equation is
\[
T_n'(t) = -\lambda_n \kappa T_n(t) \implies T_n(t) = \alpha_n e^{-\lambda_n \kappa t},
\]
for an arbitrary constant $\alpha_n \in \mathbf{R}$.

By linearity and homogeneity of the PDE and boundary condition, superposition gives
\[
\varphi(\mathbf{x}, t) = \sum_{n=1}^\infty \alpha_n e^{-\lambda_n \kappa t}\psi_n(\mathbf{x}).
\]
Imposing the initial condition $\varphi(\mathbf{x}, 0) = f(\mathbf{x})$ in $\Omega$:
\[
f(\mathbf{x}) = \sum_{n=1}^\infty \alpha_n \psi_n(\mathbf{x}),
\]
which exists by completeness (assumption (d), assuming $f = 0$ on $\partial\Omega$). Multiplying by $\psi_m(\mathbf{x})$, integrating over $\Omega$ (term-by-term), and applying orthogonality from part (ii),
\[
\int_\Omega f(\mathbf{x})\psi_m(\mathbf{x})\,dV = \sum_{n=1}^\infty \alpha_n \int_\Omega \psi_n(\mathbf{x})\psi_m(\mathbf{x})\,dV = \alpha_m \int_\Omega \psi_m(\mathbf{x})^2\,dV.
\]
Since $\psi_m \not\equiv 0$, $\int_\Omega \psi_m^2\,dV > 0$, so relabelling $m \to n$ gives
\[
\varphi(\mathbf{x}, t) = \sum_{n=1}^\infty \alpha_n e^{-\lambda_n \kappa t}\psi_n(\mathbf{x}), \qquad \alpha_n = \frac{\int_\Omega f\psi_n\,dV}{\int_\Omega \psi_n^2\,dV}.
\]

**2. Formal interpretation as $\varphi(\mathbf{x},t) = e^{\kappa t\Delta}\varphi(\mathbf{x},0)$:**

Since $\Delta\psi_n = -\lambda_n \psi_n$, induction on $p \in \mathbf{Z}_{\ge 0}$ gives
\[
\Delta^p \psi_n = (-\lambda_n)^p \psi_n.
\]
Applying the formal operator series $e^{\kappa t\Delta} = \sum_{p=0}^\infty \frac{(\kappa t)^p}{p!}\Delta^p$ to an eigenfunction $\psi_n$ yields
\[
e^{\kappa t\Delta}\psi_n = \sum_{p=0}^\infty \frac{(\kappa t)^p}{p!}\Delta^p \psi_n = \sum_{p=0}^\infty \frac{(-\lambda_n \kappa t)^p}{p!}\psi_n = e^{-\lambda_n \kappa t}\psi_n.
\]
Writing the initial condition as $\varphi(\mathbf{x}, 0) = f(\mathbf{x}) = \sum_{n=1}^\infty \alpha_n \psi_n(\mathbf{x})$ and applying $e^{\kappa t\Delta}$ formally term-by-term (by linearity, interchanging the sums over $p$ and $n$):
\[
e^{\kappa t\Delta}\varphi(\mathbf{x}, 0) = e^{\kappa t\Delta}\!\left(\sum_{n=1}^\infty \alpha_n \psi_n(\mathbf{x})\right) = \sum_{n=1}^\infty \alpha_n \left(e^{\kappa t\Delta}\psi_n(\mathbf{x})\right) = \sum_{n=1}^\infty \alpha_n e^{-\lambda_n \kappa t}\psi_n(\mathbf{x}) = \varphi(\mathbf{x}, t).
\]
(Equivalently, differentiating the PDE $\varphi_t = \kappa\Delta\varphi$ repeatedly gives $\partial_t^p \varphi = (\kappa\Delta)^p\varphi$, so a formal Taylor series in $t$ about $t=0$ reads $\varphi(\mathbf{x},t) = \sum_{p=0}^\infty \frac{t^p}{p!}\partial_t^p \varphi(\mathbf{x},0) = \sum_{p=0}^\infty \frac{(\kappa t)^p}{p!}\Delta^p\varphi(\mathbf{x},0)$. The interpretation is *formal* because for general $f$, $\Delta^p f$ need not exist or vanish on $\partial\Omega$, and the power series need not converge, whereas the eigenfunction series converges rapidly for $t>0$ since $\lambda_n > 0$.)

### Part (iv): Separation of variables for the wave problem $(B)$ and formal trigonometric operators

**1. Solution of $(B)$ by separation of variables:**

Seek non-trivial separable solutions $\varphi(\mathbf{x}, t) = X(\mathbf{x})T(t)$ of $\varphi_{tt} - c^2\Delta\varphi = 0$ in $\Omega\times(0,\infty)$ with $\varphi = 0$ on $\partial\Omega\times[0,\infty)$. Substituting into the PDE and dividing by $c^2 X(\mathbf{x})T(t)$ gives
\[
\frac{T''(t)}{c^2 T(t)} = \frac{\Delta X(\mathbf{x})}{X(\mathbf{x})} = -\lambda.
\]
Again, the spatial boundary value problem is $-\Delta X = \lambda X$ in $\Omega$, $X = 0$ on $\partial\Omega$, with eigenpairs $\{(\psi_n, \lambda_n)\}_{n=1}^\infty$. Crucially, by part (i), $\lambda_n > 0$ for every $n$, so $\sqrt{\lambda_n} > 0$ is real and strictly positive.

For each $n$, the temporal ODE is
\[
T_n''(t) + c^2 \lambda_n T_n(t) = 0,
\]
with general solution
\[
T_n(t) = \beta_n \cos\!\left(ct\sqrt{\lambda_n}\right) + \tilde\gamma_n \sin\!\left(ct\sqrt{\lambda_n}\right),
\]
for real constants $\beta_n, \tilde\gamma_n$. By superposition,
\[
\varphi(\mathbf{x}, t) = \sum_{n=1}^\infty \left[\beta_n \cos\!\left(ct\sqrt{\lambda_n}\right) + \tilde\gamma_n \sin\!\left(ct\sqrt{\lambda_n}\right)\right]\psi_n(\mathbf{x}).
\]
Differentiating term-by-term in $t$:
\[
\varphi_t(\mathbf{x}, t) = \sum_{n=1}^\infty c\sqrt{\lambda_n}\left[-\beta_n \sin\!\left(ct\sqrt{\lambda_n}\right) + \tilde\gamma_n \cos\!\left(ct\sqrt{\lambda_n}\right)\right]\psi_n(\mathbf{x}).
\]
Imposing the initial conditions $\varphi(\mathbf{x}, 0) = g(\mathbf{x})$ and $\varphi_t(\mathbf{x}, 0) = h(\mathbf{x})$:
\[
g(\mathbf{x}) = \sum_{n=1}^\infty \beta_n \psi_n(\mathbf{x}), \qquad h(\mathbf{x}) = \sum_{n=1}^\infty c\sqrt{\lambda_n}\,\tilde\gamma_n \psi_n(\mathbf{x}).
\]
Using orthogonality from part (ii) (multiplying by $\psi_m$ and integrating over $\Omega$, noting $c\sqrt{\lambda_m} > 0$), we obtain $\tilde\gamma_n = \gamma_n / (c\sqrt{\lambda_n})$ where
\[
\beta_n = \frac{\int_\Omega g\psi_n\,dV}{\int_\Omega \psi_n^2\,dV}, \qquad \gamma_n = \frac{\int_\Omega h\psi_n\,dV}{\int_\Omega \psi_n^2\,dV}.
\]
Thus the solution to $(B)$ is
\[
\varphi(\mathbf{x}, t) = \sum_{n=1}^\infty \left[ \gamma_n \frac{\sin\!\left(ct\sqrt{\lambda_n}\right)}{c\sqrt{\lambda_n}} + \beta_n \cos\!\left(ct\sqrt{\lambda_n}\right) \right]\psi_n(\mathbf{x}).
\]

**2. Relation to the formal operator expression:**

Because $\cos(z)$ and $\frac{\sin(z)}{z}$ are even entire functions of $z$, their Maclaurin series involve only *integer* powers of $z^2 = -c^2 t^2 \Delta$ (so no fractional powers of $-\Delta$ are required in the formal power series):
\[
\cos\!\left(ct\sqrt{-\Delta}\right) = \sum_{p=0}^\infty \frac{(-1)^p (ct)^{2p}}{(2p)!}(-\Delta)^p,
\]
\[
\frac{\sin\!\left(ct\sqrt{-\Delta}\right)}{c\sqrt{-\Delta}} = \sum_{p=0}^\infty \frac{(-1)^p c^{2p} t^{2p+1}}{(2p+1)!}(-\Delta)^p.
\]
(Equivalently, since $-\Delta\psi_n = \lambda_n \psi_n$ with $\lambda_n > 0$, we may also define $\sqrt{-\Delta}\,\psi_n = \sqrt{\lambda_n}\,\psi_n$ directly on the eigenbasis.) Using $(-\Delta)^p \psi_n = \lambda_n^p \psi_n$, we find on each eigenfunction $\psi_n$:
\[
\cos\!\left(ct\sqrt{-\Delta}\right)\psi_n = \sum_{p=0}^\infty \frac{(-1)^p (ct)^{2p}\lambda_n^p}{(2p)!}\psi_n = \cos\!\left(ct\sqrt{\lambda_n}\right)\psi_n,
\]
\[
\frac{\sin\!\left(ct\sqrt{-\Delta}\right)}{c\sqrt{-\Delta}}\psi_n = \sum_{p=0}^\infty \frac{(-1)^p c^{2p} t^{2p+1}\lambda_n^p}{(2p+1)!}\psi_n = \frac{\sin\!\left(ct\sqrt{\lambda_n}\right)}{c\sqrt{\lambda_n}}\psi_n.
\]
Substituting the eigenfunction expansions of the initial data $\varphi(\mathbf{x}, 0) = g(\mathbf{x}) = \sum_{n=1}^\infty \beta_n \psi_n(\mathbf{x})$ and $\varphi_t(\mathbf{x}, 0) = h(\mathbf{x}) = \sum_{n=1}^\infty \gamma_n \psi_n(\mathbf{x})$ and applying the formal operators term-by-term:
\[
\frac{\sin\!\left(ct\sqrt{-\Delta}\right)}{c\sqrt{-\Delta}}\varphi_t(\mathbf{x},0) + \cos\!\left(ct\sqrt{-\Delta}\right)\varphi(\mathbf{x},0)
= \sum_{n=1}^\infty \left[\gamma_n \frac{\sin\!\left(ct\sqrt{\lambda_n}\right)}{c\sqrt{\lambda_n}} + \beta_n \cos\!\left(ct\sqrt{\lambda_n}\right)\right]\psi_n(\mathbf{x}),
\]
which coincides identically with the separated-variables solution of $(B)$.

## Final Answers

- **(i)** $\lambda_n = \dfrac{\int_\Omega |\nabla\psi_n|^2\,dV}{\int_\Omega \psi_n^2\,dV} > 0$ (equality $\lambda_n = 0$ would force $\nabla\psi_n \equiv \mathbf{0}$ and hence $\psi_n \equiv 0$ by $\psi_n|_{\partial\Omega}=0$); $(\lambda_n - \lambda_m)\int_\Omega \psi_n\psi_m\,dV = 0$ by Green's second identity, so $\int_\Omega \psi_n\psi_m\,dV = 0$ for $\lambda_n \neq \lambda_m$.
- **(ii)** Each eigenspace $E_\lambda$ is a finite-dimensional real vector space (assumption (c)); applying Gram–Schmidt with respect to $\langle u,v\rangle = \int_\Omega uv\,dV$ within each $E_\lambda$ yields an orthogonal basis preserving (a)–(d), while distinct eigenspaces are already orthogonal by (i).
- **(iii)** $\varphi(\mathbf{x},t) = \sum_{n=1}^\infty \alpha_n e^{-\lambda_n\kappa t}\psi_n(\mathbf{x})$ with $\alpha_n = \dfrac{\int_\Omega f\psi_n\,dV}{\int_\Omega \psi_n^2\,dV}$; since $\Delta^p\psi_n = (-\lambda_n)^p\psi_n$, formally $e^{\kappa t\Delta}\psi_n = e^{-\lambda_n\kappa t}\psi_n$, giving $\varphi(\mathbf{x},t) = e^{\kappa t\Delta}\varphi(\mathbf{x},0)$ term-by-term.
- **(iv)** The solution of $(B)$ is
  \[
  \varphi(\mathbf{x},t) = \sum_{n=1}^\infty \left[\gamma_n \frac{\sin\!\left(ct\sqrt{\lambda_n}\right)}{c\sqrt{\lambda_n}} + \beta_n \cos\!\left(ct\sqrt{\lambda_n}\right)\right]\psi_n(\mathbf{x}),
  \]
  where $\beta_n = \dfrac{\int_\Omega g\psi_n\,dV}{\int_\Omega \psi_n^2\,dV}$ and $\gamma_n = \dfrac{\int_\Omega h\psi_n\,dV}{\int_\Omega \psi_n^2\,dV}$. Expanding $\cos\!\left(ct\sqrt{-\Delta}\right) = \sum_{p=0}^\infty \frac{(-1)^p (ct)^{2p}}{(2p)!}(-\Delta)^p$ and $\frac{\sin(ct\sqrt{-\Delta})}{c\sqrt{-\Delta}} = \sum_{p=0}^\infty \frac{(-1)^p c^{2p}t^{2p+1}}{(2p+1)!}(-\Delta)^p$ and using $(-\Delta)^p\psi_n = \lambda_n^p\psi_n$ reproduces this series term-by-term.

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-2_qx10_check.py`

Verbatim stdout (exit code 0):
```text
IB_methods_sheet-2_qx10_check: ALL PASSED
  part (i) max Rayleigh quotient error: 5.613e-13
  part (i) max distinct-eigenvalue orthogonality error: 1.069e-16
  part (ii) degenerate E_6 Gram-Schmidt orth error: 7.924e-14, eig eq error: 6.613e-08
  part (iii) heat IC err: 4.582e-06, PDE err: 3.954e-07, operator series err: 5.551e-17
  part (iv) wave IC err (g,h): (0.000e+00, 6.979e-11), PDE err: 1.025e-05, operator series err: 7.286e-17
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`. Every proof step and operator identity is completely rigorous and verified numerically on concrete 3D domains (including a degenerate eigenspace on $(0,\pi)^3$).
- **Phase B watchpoints:**
  - In part (i), check that the proof explicitly rules out $\lambda_n = 0$ (not just $\lambda_n \ge 0$!) by showing $\nabla\psi_n \equiv \mathbf{0}$ implies $\psi_n$ is constant and hence $0$ via the Dirichlet boundary condition on the bounded domain $\Omega$, contradicting $\psi_n \not\equiv 0$.
  - In part (ii), check that finite multiplicity (assumption (c)), linearity of $-\Delta$ and the homogeneous BC, and preservation of completeness (assumption (d)) are all addressed when invoking Gram–Schmidt.
  - In part (iv), check that positivity $\lambda_n > 0$ from (i) is used to ensure $\sqrt{\lambda_n} > 0$ (no division by zero in $\frac{\sin(ct\sqrt{\lambda_n})}{c\sqrt{\lambda_n}}$ and oscillatory rather than hyperbolic time dependence), and check the exact powers of $c$ and $t$ in the Maclaurin series of $\frac{\sin(ct\sqrt{-\Delta})}{c\sqrt{-\Delta}} = \sum_{p=0}^\infty \frac{(-1)^p c^{2p} t^{2p+1}}{(2p+1)!}(-\Delta)^p$.
