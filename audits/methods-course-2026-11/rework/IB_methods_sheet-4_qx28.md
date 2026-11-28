---
id: IB/methods/sheet-4/qx28
slug: IB_methods_sheet-4_qx28
source_archive: private/raw/sheets/IB/methods/2015-16/all-sheets.pdf
pages_checked: [7]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word for word and symbol for symbol against `IB_methods_sheet-4_qx28_p7.png` (Michaelmas term 2015, *METHODS — EXAMPLES IV*, Dr. D. Skinner, Course C10—No. B8, page 7, Question 4 under *Properties of Green's functions*). Question 4 sits cleanly near the middle of page 7.

- **Result**: `ok` (exact mathematical and verbal match).
- **Purely typographic notes**:
  - The underlined italic title `Symmetry.` is rendered as `\emph{Symmetry.}`.
  - All mathematical symbols ($G(\mathbf{r};\mathbf{r}_0)$, $\mathcal{D}$, $\mathbf{r}\ne\mathbf{r}_0$) and prose match verbatim.

## Full Solution

### Definition and Conventions

We adopt the standard Cambridge IB Methods sign convention in which the Dirichlet Green's function $G(\mathbf{r}; \mathbf{r}_0)$ for the Laplacian in a three-dimensional domain $\mathcal{D}$ (with piecewise smooth boundary $\partial\mathcal{D}$) satisfies
\[
\nabla^2_{\mathbf{r}} G(\mathbf{r}; \mathbf{r}_0) = \delta^{(3)}(\mathbf{r} - \mathbf{r}_0) \quad (\mathbf{r}, \mathbf{r}_0 \in \mathcal{D}),
\qquad
G(\mathbf{r}; \mathbf{r}_0) = 0 \quad (\mathbf{r} \in \partial\mathcal{D}).
\]
(Note: if the opposite convention $\nabla^2_{\mathbf{r}} G = -\delta^{(3)}(\mathbf{r}-\mathbf{r}_0)$ is used, $G$ simply reverses sign and the proof and symmetry result are unchanged.)

More precisely, for each fixed interior point $\mathbf{r}_0 \in \mathcal{D}$, $G(\mathbf{r}; \mathbf{r}_0)$ has the decomposition
\[
G(\mathbf{r}; \mathbf{r}_0) = -\frac{1}{4\pi |\mathbf{r} - \mathbf{r}_0|} + H(\mathbf{r}; \mathbf{r}_0),
\]
where the regular part $H(\cdot; \mathbf{r}_0)$ is harmonic ($\nabla^2_{\mathbf{r}} H = 0$) and $C^2$ throughout $\mathcal{D}$ (and continuous up to $\partial\mathcal{D}$, with $H(\mathbf{r};\mathbf{r}_0) = \frac{1}{4\pi|\mathbf{r}-\mathbf{r}_0|}$ on $\partial\mathcal{D}$). In particular, $\nabla^2_{\mathbf{r}} G(\mathbf{r}; \mathbf{r}_0) = 0$ for all $\mathbf{r} \in \mathcal{D}\setminus\{\mathbf{r}_0\}$. If $\mathcal{D}$ is unbounded, $G$ also satisfies the standard 3D decay conditions at infinity: $G(\mathbf{r};\mathbf{r}_0) = O(|\mathbf{r}|^{-1})$ and $|\nabla_{\mathbf{r}} G(\mathbf{r};\mathbf{r}_0)| = O(|\mathbf{r}|^{-2})$ as $|\mathbf{r}| \to \infty$.

---

### Green's Second Identity

For any bounded region $\Omega \subset \mathbb{R}^3$ with piecewise smooth boundary $\partial\Omega$ (and outward unit normal $\hat{\mathbf{n}}$) and any two functions $u, v \in C^2(\Omega) \cap C^1(\overline{\Omega})$, applying the Divergence Theorem to the vector field $u\nabla v - v\nabla u$ yields Green's second identity:
\[
\int_{\Omega} \bigl(u \nabla^2 v - v \nabla^2 u\bigr)\,\mathrm{d}V
= \int_{\Omega} \nabla \cdot (u\nabla v - v\nabla u)\,\mathrm{d}V
= \int_{\partial\Omega} \left(u \frac{\partial v}{\partial n} - v \frac{\partial u}{\partial n}\right)\mathrm{d}S.
\]

---

### Rigorous Proof via Excised Balls

Fix any two distinct interior points $\mathbf{r}_1, \mathbf{r}_2 \in \mathcal{D}$ with $\mathbf{r}_1 \neq \mathbf{r}_2$, and define
\[
u(\mathbf{r}) = G(\mathbf{r}; \mathbf{r}_1), \qquad v(\mathbf{r}) = G(\mathbf{r}; \mathbf{r}_2).
\]
Because $u$ is singular at $\mathbf{r}_1$ and $v$ is singular at $\mathbf{r}_2$, Green's second identity cannot be applied directly on all of $\mathcal{D}$ without first excising small balls around both singularities.

#### Case 1: $\mathcal{D}$ is bounded

Choose $\varepsilon_0 > 0$ small enough that the closed balls $\overline{B_{\varepsilon_0}(\mathbf{r}_1)}$ and $\overline{B_{\varepsilon_0}(\mathbf{r}_2)}$ are disjoint and contained strictly within $\mathcal{D}$, i.e.
\[
0 < \varepsilon_0 < \min\!\left(\tfrac{1}{2}|\mathbf{r}_1 - \mathbf{r}_2|,\; \operatorname{dist}(\mathbf{r}_1, \partial\mathcal{D}),\; \operatorname{dist}(\mathbf{r}_2, \partial\mathcal{D})\right).
\]
For any $0 < \varepsilon \le \varepsilon_0$, define the excised domain
\[
\mathcal{D}_\varepsilon = \mathcal{D} \setminus \bigl(\overline{B_\varepsilon(\mathbf{r}_1)} \cup \overline{B_\varepsilon(\mathbf{r}_2)}\bigr).
\]
Its boundary is the disjoint union of three surfaces,
\[
\partial\mathcal{D}_\varepsilon = \partial\mathcal{D} \cup S_\varepsilon(\mathbf{r}_1) \cup S_\varepsilon(\mathbf{r}_2),
\]
where $S_\varepsilon(\mathbf{r}_k) = \{\mathbf{r} : |\mathbf{r} - \mathbf{r}_k| = \varepsilon\}$ ($k = 1, 2$). Crucially, on each excised sphere $S_\varepsilon(\mathbf{r}_k)$, the outward unit normal $\hat{\mathbf{n}}$ to $\mathcal{D}_\varepsilon$ points **inward** toward $\mathbf{r}_k$: introducing spherical radius $\rho_k = |\mathbf{r} - \mathbf{r}_k|$ and outward radial unit vector $\hat{\mathbf{e}}_{\rho_k} = (\mathbf{r}-\mathbf{r}_k)/\rho_k$, we have
\[
\hat{\mathbf{n}} = -\hat{\mathbf{e}}_{\rho_k}, \qquad \frac{\partial}{\partial n} = -\frac{\partial}{\partial \rho_k} \qquad \text{on } S_\varepsilon(\mathbf{r}_k).
\]

On $\mathcal{D}_\varepsilon$, both $u$ and $v$ are $C^2$ and harmonic ($\nabla^2 u = \nabla^2 v = 0$). Applying Green's second identity on $\mathcal{D}_\varepsilon$ therefore gives
\[
0 = \int_{\mathcal{D}_\varepsilon} \bigl(u \nabla^2 v - v \nabla^2 u\bigr)\,\mathrm{d}V
= \int_{\partial\mathcal{D}} \left(u \frac{\partial v}{\partial n} - v \frac{\partial u}{\partial n}\right)\mathrm{d}S + I_1(\varepsilon) + I_2(\varepsilon),
\]
where
\[
I_k(\varepsilon) = \int_{S_\varepsilon(\mathbf{r}_k)} \left(u \frac{\partial v}{\partial n} - v \frac{\partial u}{\partial n}\right)\mathrm{d}S \qquad (k = 1, 2).
\]
By the homogeneous Dirichlet boundary condition, $u = 0$ and $v = 0$ everywhere on $\partial\mathcal{D}$, so the integral over $\partial\mathcal{D}$ vanishes identically, leaving
\[
I_1(\varepsilon) + I_2(\varepsilon) = 0 \qquad \text{for all } 0 < \varepsilon \le \varepsilon_0.
\]

Now we evaluate $\lim_{\varepsilon \to 0^+} I_1(\varepsilon)$:
- On $\overline{B_{\varepsilon_0}(\mathbf{r}_1)}$, $u(\mathbf{r}) = -\dfrac{1}{4\pi \rho_1} + H(\mathbf{r}; \mathbf{r}_1)$ where $H(\cdot; \mathbf{r}_1)$ is $C^2$, so $|H|$ and $|\nabla H|$ are bounded by some constant $M_1$. On $S_\varepsilon(\mathbf{r}_1)$,
  \[
  \frac{\partial u}{\partial n} = -\frac{\partial}{\partial \rho_1}\left(-\frac{1}{4\pi \rho_1} + H(\mathbf{r};\mathbf{r}_1)\right)\Bigg|_{\rho_1=\varepsilon}
  = -\frac{1}{4\pi \varepsilon^2} - \frac{\partial H}{\partial \rho_1}(\mathbf{r}; \mathbf{r}_1).
  \]
- Since $\mathbf{r}_1 \neq \mathbf{r}_2$, $v(\mathbf{r}) = G(\mathbf{r}; \mathbf{r}_2)$ is $C^2$ on $\overline{B_{\varepsilon_0}(\mathbf{r}_1)}$, so $|v| \le M_2$ and $|\nabla v| \le M_3$ there for constants $M_2, M_3$.

Splitting $I_1(\varepsilon)$ into two integrals:
1. For the first integral, since $\operatorname{Area}(S_\varepsilon(\mathbf{r}_1)) = 4\pi \varepsilon^2$,
   \[
   \left|\int_{S_\varepsilon(\mathbf{r}_1)} u \frac{\partial v}{\partial n}\,\mathrm{d}S\right|
   \le 4\pi\varepsilon^2 \left(\frac{1}{4\pi\varepsilon} + M_1\right) M_3
   = \varepsilon M_3 + 4\pi\varepsilon^2 M_1 M_3 \xrightarrow[\varepsilon\to 0^+]{} 0.
   \]
2. For the second integral,
   \[
   -\int_{S_\varepsilon(\mathbf{r}_1)} v \frac{\partial u}{\partial n}\,\mathrm{d}S
   = \frac{1}{4\pi\varepsilon^2} \int_{S_\varepsilon(\mathbf{r}_1)} v(\mathbf{r})\,\mathrm{d}S
   + \int_{S_\varepsilon(\mathbf{r}_1)} v(\mathbf{r})\,\frac{\partial H(\mathbf{r};\mathbf{r}_1)}{\partial \rho_1}\,\mathrm{d}S.
   \]
   The second term is bounded in magnitude by $4\pi\varepsilon^2 M_2 M_1 = O(\varepsilon^2) \to 0$. The first term is the spherical average of the continuous function $v$ over $S_\varepsilon(\mathbf{r}_1)$:
   \[
   \left|\frac{1}{4\pi\varepsilon^2} \int_{S_\varepsilon(\mathbf{r}_1)} v(\mathbf{r})\,\mathrm{d}S - v(\mathbf{r}_1)\right|
   \le \sup_{\mathbf{r}\in S_\varepsilon(\mathbf{r}_1)} |v(\mathbf{r}) - v(\mathbf{r}_1)| \xrightarrow[\varepsilon\to 0^+]{} 0.
   \]
   (In fact, since $v$ is harmonic on $\overline{B_{\varepsilon_0}(\mathbf{r}_1)}$, the spherical mean equals $v(\mathbf{r}_1)$ identically for every $\varepsilon$, and similarly $I_1(\varepsilon) = v(\mathbf{r}_1)$ holds for every finite $\varepsilon \le \varepsilon_0$.)

Thus
\[
\lim_{\varepsilon \to 0^+} I_1(\varepsilon) = v(\mathbf{r}_1) = G(\mathbf{r}_1; \mathbf{r}_2).
\]

By the exact same argument on $S_\varepsilon(\mathbf{r}_2)$ (or directly by exchanging $u \leftrightarrow v$ and $\mathbf{r}_1 \leftrightarrow \mathbf{r}_2$, noting that the integrand $u\frac{\partial v}{\partial n} - v\frac{\partial u}{\partial n}$ is antisymmetric under $u \leftrightarrow v$):
- $-\int_{S_\varepsilon(\mathbf{r}_2)} v \frac{\partial u}{\partial n}\,\mathrm{d}S = O(\varepsilon) \to 0$, while
- $\int_{S_\varepsilon(\mathbf{r}_2)} u \frac{\partial v}{\partial n}\,\mathrm{d}S = -\frac{1}{4\pi\varepsilon^2}\int_{S_\varepsilon(\mathbf{r}_2)} u(\mathbf{r})\,\mathrm{d}S + O(\varepsilon^2) \to -u(\mathbf{r}_2) = -G(\mathbf{r}_2; \mathbf{r}_1)$.

Hence
\[
\lim_{\varepsilon \to 0^+} I_2(\varepsilon) = -u(\mathbf{r}_2) = -G(\mathbf{r}_2; \mathbf{r}_1).
\]

Taking $\varepsilon \to 0^+$ in $I_1(\varepsilon) + I_2(\varepsilon) = 0$ gives
\[
G(\mathbf{r}_1; \mathbf{r}_2) - G(\mathbf{r}_2; \mathbf{r}_1) = 0 \implies G(\mathbf{r}_1; \mathbf{r}_2) = G(\mathbf{r}_2; \mathbf{r}_1).
\]

#### Case 2: $\mathcal{D}$ is unbounded

If $\mathcal{D}$ is unbounded, take a large ball $B_R(\mathbf{0})$ of radius $R > \max(|\mathbf{r}_1|, |\mathbf{r}_2|) + \varepsilon_0$ and apply Green's second identity on
\[
\mathcal{D}_{\varepsilon, R} = \bigl(\mathcal{D} \cap B_R(\mathbf{0})\bigr) \setminus \bigl(\overline{B_\varepsilon(\mathbf{r}_1)} \cup \overline{B_\varepsilon(\mathbf{r}_2)}\bigr).
\]
This introduces one extra surface integral over $\Sigma_R = \mathcal{D} \cap \partial B_R(\mathbf{0})$:
\[
I_1(\varepsilon) + I_2(\varepsilon) + \int_{\Sigma_R} \left(u \frac{\partial v}{\partial n} - v \frac{\partial u}{\partial n}\right)\mathrm{d}S = 0.
\]
Using the 3D decay conditions $u, v = O(R^{-1})$ and $|\nabla u|, |\nabla v| = O(R^{-2})$ on $\Sigma_R$, together with $\operatorname{Area}(\Sigma_R) \le 4\pi R^2$:
\[
\left|\int_{\Sigma_R} \left(u \frac{\partial v}{\partial n} - v \frac{\partial u}{\partial n}\right)\mathrm{d}S\right|
\le 4\pi R^2 \cdot O(R^{-1} \cdot R^{-2}) = O(R^{-1}) \xrightarrow[R\to\infty]{} 0.
\]
Letting $R \to \infty$ and $\varepsilon \to 0^+$ again gives $G(\mathbf{r}_1; \mathbf{r}_2) = G(\mathbf{r}_2; \mathbf{r}_1)$.

Relabelling $\mathbf{r}_1 = \mathbf{r}$ and $\mathbf{r}_2 = \mathbf{r}_0$, we conclude that
\[
G(\mathbf{r}; \mathbf{r}_0) = G(\mathbf{r}_0; \mathbf{r}) \qquad \text{for all } \mathbf{r} \neq \mathbf{r}_0 \text{ in } \mathcal{D}.
\]

*(Remark: In formal distributional notation, applying Green's second identity directly on $\mathcal{D}$ gives $0 = \int_{\mathcal{D}} [G(\mathbf{r};\mathbf{r}_1)\delta^{(3)}(\mathbf{r}-\mathbf{r}_2) - G(\mathbf{r};\mathbf{r}_2)\delta^{(3)}(\mathbf{r}-\mathbf{r}_1)]\,\mathrm{d}V = G(\mathbf{r}_2;\mathbf{r}_1) - G(\mathbf{r}_1;\mathbf{r}_2)$; the excised-spheres argument above is the rigorous justification of this formal manipulation.)*

## Final Answers

For any Dirichlet Green's function $G(\mathbf{r};\mathbf{r}_0)$ of the 3D Laplacian on a domain $\mathcal{D}$ (bounded or unbounded, with standard decay $G = O(R^{-1})$, $|\nabla G| = O(R^{-2})$ at infinity if unbounded), applying Green's second identity to $u(\mathbf{r}) = G(\mathbf{r};\mathbf{r}_1)$ and $v(\mathbf{r}) = G(\mathbf{r};\mathbf{r}_2)$ on $\mathcal{D} \setminus (\overline{B_\varepsilon(\mathbf{r}_1)} \cup \overline{B_\varepsilon(\mathbf{r}_2)})$ (intersected with $B_R(\mathbf{0})$ if unbounded) and taking $\varepsilon \to 0^+$ (and $R \to \infty$) yields
\[
G(\mathbf{r}; \mathbf{r}_0) = G(\mathbf{r}_0; \mathbf{r}) \qquad \text{for all } \mathbf{r} \neq \mathbf{r}_0 \text{ in } \mathcal{D}.
\]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-4_qx28_check.py` (exit code 0).

Tests both a bounded 3D domain (ball $|\mathbf{r}| < a$, Kelvin inversion Green's function) and an unbounded 3D domain (half-space $z > 0$, method of images), directly integrating the excised-sphere surface integrals $I_1(\varepsilon), I_2(\varepsilon)$ and outer hemisphere integral $I_R$ via Gauss–Legendre/trapezoidal spherical quadrature.

```text
OK: Ball (a=2.5) — max symmetry err=2.78e-17, max BC err=1.39e-17, |laplacian G|=2.78e-09.
OK: Excised-sphere Green's 2nd identity on ball (eps=0.1..0.005): I1=-0.0120536373 (G(r1;r2)=-0.0120536373), I2=0.0120536373 (-G(r2;r1)=0.0120536373), |I1+I2|=3.30e-17, int|u dv/dn|=3.88e-05.
OK: Half-space z>0 verified — |I1+I2|=2.60e-17, outer hemisphere |I_R|: R=10 -> 2.24e-21, R=50 -> 2.57e-22, R=200 -> 7.98e-24.
```

## Confidence & Edge-Case Notes

- **Confidence**: `high`.
- **Phase B watchpoints**:
  - Check whether the existing solution relies solely on the formal delta-function substitution $\int_{\mathcal{D}} (u \nabla^2 v - v \nabla^2 u)\,\mathrm{d}V$ without justifying the excision of small balls around $\mathbf{r}_1$ and $\mathbf{r}_2$ (where $u$ and $v$ are not $C^2$). Excising $B_\varepsilon(\mathbf{r}_1)$ and $B_\varepsilon(\mathbf{r}_2)$ is standard rigorous IB Methods technique.
  - Check the sign of the outward normal on the excised spheres $S_\varepsilon(\mathbf{r}_k)$: $\hat{\mathbf{n}} = -\hat{\mathbf{e}}_{\rho_k}$, so $\frac{\partial}{\partial n}\left(-\frac{1}{4\pi\rho_k}\right) = -\frac{1}{4\pi\varepsilon^2}$. A sign mistake on $\hat{\mathbf{n}}$ flips both $I_1$ and $I_2$ simultaneously (so the final conclusion still comes out right, masking the intermediate sign error!).
  - Check whether unbounded domains $\mathcal{D}$ are addressed (since the question states "an *arbitrary* three-dimensional domain $\mathcal{D}$"), including the decay rate $O(R^2 \cdot R^{-1} \cdot R^{-2}) = O(R^{-1}) \to 0$.
