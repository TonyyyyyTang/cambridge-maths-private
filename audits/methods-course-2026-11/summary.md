# Whole-Course Cross-Model Audit Summary: IB Methods

- **Auditor:** Gemini 4 via Antigravity (high-effort reasoning), cross-model to the Codex (GPT) solver
- **Date:** 2026-11-28 | **Seed:** `20261201` | **Sample:** 68 / 336 `ai-draft` units (20.2%: 31 / 154 sheet, 37 / 182 Tripos)
- **Blindness:** Clean — all 68 blind reworks (`rework/<slug>.md`), all 68 Python verification scripts (`rework/<slug>_check.py`, all passing exit code 0), and all 68 visual PDF/PS transcription checks were completed and frozen before opening any stored `solution.tex`, unit metadata, or prior audit file.

## Error Rates & Verdict Breakdown

| Metric | Overall (`n = 68`) | Example Sheets (`n = 31`) | Tripos (`n = 37`) |
|---|---:|---:|---:|
| **Substantive error rate** | **0 / 68 (0.0%)** | 0 / 31 (0.0%) | 0 / 37 (0.0%) |
| **Transcription error rate** | **0 / 68 (0.0%)** | 0 / 31 (0.0%) | 0 / 37 (0.0%) |
| **Minor issue rate** | **21 / 68 (30.9%)** | 7 / 31 (22.6%) | 14 / 37 (37.8%) |
| **Unsure** | **0 / 68 (0.0%)** | 0 / 31 (0.0%) | 0 / 37 (0.0%) |
| **Agree (exact results & valid argument)** | **47 / 68 (69.1%)** | 24 / 31 (77.4%) | 23 / 37 (62.2%) |

- **Full-course review triggered?** **No** (`0 / 68` substantive errors). Mathematical accuracy across the entire sample is exceptionally strong: every formula, sign, constant, asymptotic, edge case (including notorious traps such as the hypergeometric sign on $(0,1)$, non-self-adjoint Green's functions, Bessel denominator sign $\tilde\lambda - j_n^2$, exterior sphere image weight $a/|\mathbf{x}_0|$ and far-field decay, $n=1,2,\ge 3$ octant images, the $a>0$ vs $a>1$ defect + $a=1$ counterexample in `2010/P4/q17B`, and unrestricted $a_0$ in `2023/P2/q14A`), verbatim question text, and all 3 stored solution diagrams (`quartic-and-extensions.png`, `gaussian-pair.png`, `robin-roots.png`) verified accurately.
- **Structural pattern behind the `minor` rate:** Ten of the 14 Tripos `minor` verdicts occur on older **2002–2015 Section II** solutions drafted as ultra-compressed 1–7 line outlines that assert rather than prove explicit "Show that" / "Derive" / "State and prove" parts. By contrast, **9 of 10 sampled 2016–2024 Tripos solutions** are clean `agree`.

## Questions Needing Attention (Ordered by Severity)

All 21 `minor` units have exact results; verbatim quotes and validated drop-in restricted-LaTeX replacements are in `report.yaml`.

### 1. Severe Minor — Missing Required Statement/Object, Factual Slip, or Unproved Key Theorem (7 units)

1. **`IB/2011/P1/q14A`** *(Tripos)*: Asked to *"State and prove Parseval's theorem"*, never states Parseval's theorem (neither hypotheses nor the $N\to\infty$ identity); invokes off-schedule periodic-$L^2$ density and Stone–Weierstrass without checking hypotheses.
2. **`IB/2010/P4/q17B`** *(Tripos)*: Asked to *"Construct an appropriate Green's function $G(\mathbf{r};\mathbf{r}_0)$ satisfying $\partial G/\partial z=0$ at $z=0$"*, describes it in words without writing $G(\mathbf{r};\mathbf{r}_0)$ or verifying $\partial G/\partial z|_{z=0}=0$, compresses Green's third identity, and includes `\textbf{Faculty correction.}` meta-commentary.
3. **`IB/2006/P1/q14A`** *(Tripos)*: Single 7-sentence paragraph for a 20-mark Section II question; skips all three "show that" proofs and misstates *"the two off-diagonal elements involving the third direction"* (there are four: $M_{13}, M_{23}, M_{31}, M_{32}$).
4. **`IB/2009/P4/q5B`** *(Tripos)*: Mislabels the inequality $4 = r^2 + xy \ge r^2/2$ as an *"identity"*, gives no justification for the bound or equality condition $y=-x$, and bypasses the Lagrange multiplier / determinant method set up in Parts 1–2.
5. **`IB/methods/sheet-3/qx07`** *(sheet)*: Cites Wirtinger's inequality $\int_0^T u^2\,dt \le (T/\pi)^2\int_0^T \dot u^2\,dt$ as an unproved black box carrying the entire lower bound (provable in 3 lines via IB Fourier sine series and Parseval), bypassing Euler–Lagrange.
6. **`IB/2011/P2/q5A`** *(Tripos)*: Never defines $\langle u,v\rangle$, skips the integration-by-parts derivation, replaces IB regular-endpoint self-adjointness with unproved Part II/III maximal-domain claims, skips checking $\lambda_n\neq\lambda_m$ via $m,n\ge 0$, omits Rodrigues' formula and the Beta integral, and ends with an expression rather than an equation.
7. **`IB/2019/P3/q15D`** *(Tripos)*: Jumps straight from $\frac{H(T)}{4\pi c R}[\delta(R-cT)-\delta(R+cT)]$ to $\frac{\delta(T-R/c)}{4\pi c^2 R}$ without stating $c>0$ or justifying discarding $\delta(R+cT)$, dropping $H(T)$, and the $1/c$ Jacobian; skips Part 1 differentiation, Part 2 angular-integral substitution, and Part 3 radial integration.

### 2. Standard Minor — Explicitly Requested Derivations Skipped (8 units)

8. **`IB/methods/sheet-2/qx10`** *(sheet)*: Instructed to solve (A) and (B) *"using separation of variables"* and show $\alpha_n$, skips separation of variables and orthogonality projection, omits $\lambda_n>0\implies\sqrt{\lambda_n}>0$, and omits the $(-\Delta)^p$ power series.
9. **`IB/methods/sheet-2/qx25`** *(sheet)*: Instructed *"Using integration by parts ..., show that"*, skips the IBP step and the $\beta = j_n + \epsilon \to j_n$ numerator/denominator cancellation; leaves *"Simplicity follows by ODE uniqueness"* unnamed.
10. **`IB/2012/P1/q14C`** *(Tripos)*: Never defines unweighted $\langle u,v\rangle$, skips both requested IBP derivations (self-adjointness and Rayleigh quotient), substitutes Part II adjoint-domain jargon, and omits $q\ge 0, k>0, y(0)=0$ in strict positivity.
11. **`IB/2015/P1/q14C`** *(Tripos)*: Skips the requested IBP derivation in (ii), renames $y_1, y_2 \to u, v$, drops integral limits in (iii), skips deriving $G(x,\xi)$ in (iv), and uses Part II *"maximal domain"*.
12. **`IB/2005/P4/q5H`** *(Tripos)*: Skips all three Section I "Show how / Show that" derivations (chain rule for $4\partial_\xi\partial_\eta$, $y_t/y_x$ cross-term cancellation and reflection, $2L$-periodicity of $g'$ and $dE/dt=0$) and introduces undefined $u$.
13. **`IB/methods/sheet-2/qx20`** *(sheet)*: Skips writing and evaluating the convolution integral for $f(t)=e^{-t}$, skips justifying $G=0$ ($s<0$) and integration constants, omits $y',y'',y''',y''''$ in the verification, and renames $\tau \to t'$.
14. **`IB/methods/sheet-4/q09`** *(sheet)*: Skips the core algebraic expansion proving $|\mathbf{x}-\mathbf{x}_0^*|^2=(a^2/|\mathbf{x}_0|^2)|\mathbf{x}-\mathbf{x}_0|^2$ on $|\mathbf{x}|=a$ and the $|\mathbf{x}_0^*|$ check, renames all variables, drops `i)`/`ii)`, and never writes out the exterior formula in ii).
15. **`IB/methods/sheet-2/qx03`** *(sheet)*: Invokes $(2n+1)P_n=P'_{n+1}-P'_{n-1}$ without proof (not in hint), leaves $a_n = V[P_{n-1}(0)-P_{n+1}(0)]$ unevaluated with no closed form or leading terms, and skips spherical separation of variables.

### 3. Compression Minor — Multi-Part Section II / Sheet Proofs Asserted Without Algebra (6 units)

16. **`IB/2007/P4/q16E`** *(Tripos)*: 3-line Section II answer skips the Beltrami chain/product-rule cancellation and `= C`, omits deriving $T$ (leaving $r'$ undefined), skips $u'=-r'/r^2$ and $\cos\theta+\sin\theta > 0$, and bypasses Part 1 Beltrami.
17. **`IB/2003/P2/q11C`** *(Tripos)*: 3-line Section II answer skips defining isotropy and expanding $c'_{ijkl}$, skips $\beta e_{ij}+\gamma e_{ji}$, skips proving $p$ is a scalar and $d_{ij}$ traceless ($\delta_{kk}=3$), skips $e_{ij}e_{ij}$ and $E$, and gives no explicit test strain tensors.
18. **`IB/2002/P2/q11C`** *(Tripos)*: 3-line Section II answer asserts tensor transformation/symmetry without checking orthogonality contractions, skips $T_{ij}v_j$, dimension count $1+2=3$, and explicit eigenvalue $0$ at $\mathbf{B}=\mathbf{0}$, and skips the Levi-Civita contraction.
19. **`IB/2008/P1/q14D`** *(Tripos)*: Single-paragraph, zero-display Section II answer skips all derivations in both "Show that" parts (never checks $C\ne 0$ or $k>c(x)$, skips solving for $y'$, never writes $y'=x/\sqrt{x_2^2-x^2}$ or integrates it).
20. **`IB/2004/P2/q17B`** *(Tripos)*: 7-line Section II answer skips the parallel-axis expansion and $\int r_k\,dm=0$, asserts cube average $a^2/3$ without integrating or defining $M=8\rho a^3$, and skips $|X|^2=3a^2$ and $I(P)e, I(P)v$.
21. **`IB/methods/sheet-4/qx02`** *(sheet)*: 6-line proof omits orthogonality ($QQ^T=I$) in (i), skips axis-coincidence verification and tensoriality of $u_iu_j$ in (ii), omits uniqueness in the arbitrary frame from Frame 1 components, and drops `(i)`/`(ii)` structure.

### 4. High-Priority Style & Schedule Observations (for a later style pass)

- **Off-schedule Part II/III methods** *(against `private/raw/schedules/schedules-2025-26.pdf`)*: Highest priority are **`IB/methods/sheet-4/q05`** (Bochner $L^1([0,T];L^2)$ contraction semigroup $S(t)$ instead of heat-kernel convolution), **`IB/2011/P2/q5A`** (unbounded-operator maximal/adjoint domain proof), **`IB/2018/P1/q14C`** (Part II $L^1$/Fubini/distributional-pairing derivation bypassing Part 1 dual convolution), **`IB/2011/P1/q14A`** (periodic-$L^2$ density + Stone–Weierstrass), and **`IB/2012/P1/q14C`** (adjoint-domain jargon). Lighter Part II phrasing (`boundary trace`, `maximal domain`, `dominated convergence`, `almost everywhere`, `maximum principle`) appears across 13 additional units.
- **Historical 2000–2009 syllabus context (not off-schedule solution errors):** 13 sampled questions come from older IB Methods topics later moved to IA Vector Calculus or IB Variational Principles (6 Cartesian tensors, 5 Calculus of Variations, 2 multivariable extrema).
- **Standalone question context (`IB/methods/sheet-1/qx24`, `qx25`, `qx27`):** Transcriptions match `2007-08/sheet-1.pdf` p. 2 verbatim, but the shared section header defining $\tilde f(k)=\int_{-\infty}^\infty f(x)e^{-ikx}\,dx$ sits above Q6–Q11; `qx25` does not restate it.
- **Systematic `templates/STYLE.md` drift:** Widespread imperative voice (`Write...`, `Put...`); dropped part labels (7 units) and `\textbf{(i)}`/`\textbf{(a)}` bold labels (9 units); American `-ize`/`disk's` spellings (12 units) and unhyphenated `nonzero`/`nonnegative` (27 units); editorial meta-commentary (4 units); `"Green function"` without `'s` (6 units); omitted integral limits/`dx` (7 units).
