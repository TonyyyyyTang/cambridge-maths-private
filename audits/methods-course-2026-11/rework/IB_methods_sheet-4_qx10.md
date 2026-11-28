---
id: IB/methods/sheet-4/qx10
slug: IB_methods_sheet-4_qx10
source_archive: private/raw/sheets/IB/methods/2000-01/sheet-4.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word for word and symbol for symbol against `IB_methods_sheet-4_qx10_p2.png` (Michaelmas term 2000, *METHODS — EXAMPLES IV*, Dr E.P.S. Shellard, page 2, Question 10). Question 10 sits cleanly in the upper-middle section of page 2.

- **Result**: `ok` (exact mathematical and verbal match).
- **Purely typographic notes**:
  - The underlined italic heading `Symmetric tensors.` is rendered as `\emph{Symmetric tensors.}`.
  - Every index (`T_{ii}`, `T_{ik}T_{ki}`, `T_{ik}\,T_{km}\,T_{mi}`), fraction (`\tfrac12`, `\tfrac16`), sign, and exponent in the cubic equation matches the printed sheet exactly.

## Full Solution

### Part (i): Direct proof of rotational invariance

Consider an arbitrary rank-$2$ Cartesian tensor $T_{ik}$ in $\mathbb{R}^3$ (not assumed symmetric). Under a rotation of the Cartesian coordinate axes represented by an orthogonal matrix $L = (\ell_{ij}) \in SO(3)$, the components transform according to
\[
T'_{ik} = \ell_{ip}\,\ell_{kq}\,T_{pq},
\]
where the direction-cosine matrix satisfies the orthogonality relations
\[
\ell_{ip}\,\ell_{iq} = \delta_{pq} \qquad (\text{i.e. } L^T L = I)
\]
(and similarly $\ell_{pi}\,\ell_{qi} = \delta_{pq}$). Note that we use the Einstein summation convention over repeated indices from $1$ to $3$.

We check each of the three quantities directly in the rotated frame:

1. **Invariance of $\alpha = T_{ii}$**:
   \[
   \alpha' = T'_{ii} = \ell_{ip}\,\ell_{iq}\,T_{pq} = \delta_{pq}\,T_{pq} = T_{pp} = T_{ii} = \alpha.
   \]

2. **Invariance of $\beta = T_{ik}T_{ki}$**:
   \[
   \beta' = T'_{ik}\,T'_{ki}
   = (\ell_{ip}\,\ell_{kq}\,T_{pq})(\ell_{kr}\,\ell_{is}\,T_{rs})
   = (\ell_{ip}\,\ell_{is})(\ell_{kq}\,\ell_{kr})\,T_{pq}\,T_{rs}.
   \]
   Applying $\ell_{ip}\,\ell_{is} = \delta_{ps}$ and $\ell_{kq}\,\ell_{kr} = \delta_{qr}$:
   \[
   \beta' = \delta_{ps}\,\delta_{qr}\,T_{pq}\,T_{rs} = T_{pq}\,T_{qp} = T_{ik}\,T_{ki} = \beta.
   \]

3. **Invariance of $\gamma = T_{ik}\,T_{km}\,T_{mi}$**:
   \[
   \begin{aligned}
   \gamma' = T'_{ik}\,T'_{km}\,T'_{mi}
   &= (\ell_{ip}\,\ell_{kq}\,T_{pq})(\ell_{kr}\,\ell_{ms}\,T_{rs})(\ell_{mt}\,\ell_{iu}\,T_{tu}) \\
   &= (\ell_{ip}\,\ell_{iu})(\ell_{kq}\,\ell_{kr})(\ell_{ms}\,\ell_{mt})\,T_{pq}\,T_{rs}\,T_{tu}.
   \end{aligned}
   \]
   Applying $\ell_{ip}\,\ell_{iu} = \delta_{pu}$, $\ell_{kq}\,\ell_{kr} = \delta_{qr}$, and $\ell_{ms}\,\ell_{mt} = \delta_{st}$:
   \[
   \gamma' = \delta_{pu}\,\delta_{qr}\,\delta_{st}\,T_{pq}\,T_{rs}\,T_{tu}
   = T_{pq}\,T_{qs}\,T_{sp}
   = T_{ik}\,T_{km}\,T_{mi} = \gamma.
   \]

Hence $\alpha$, $\beta$, and $\gamma$ are all scalar invariants under any rotation of the coordinate axes (indeed under any orthogonal transformation $L \in O(3)$, since only $L^T L = I$ was used, and without requiring $T_{ik}$ to be symmetric).

---

### Part (ii): Symmetric tensor, eigenvalues, and the characteristic cubic

Now suppose $T_{ik}$ is a real symmetric rank-$2$ tensor ($T_{ik} = T_{ki}$). By the spectral theorem (principal-axes theorem) for real symmetric $3\times 3$ matrices, $T_{ik}$ has three real eigenvalues $\lambda_1, \lambda_2, \lambda_3$ (counted with multiplicity) and an orthonormal basis of corresponding eigenvectors. By reversing the sign of one eigenvector if necessary, this orthonormal eigenbasis can be chosen to be right-handed, so there exists a proper rotation $L \in SO(3)$ transforming to the principal-axis frame in which $T_{ik}$ is diagonal:
\[
T'_{ik} = \begin{pmatrix}
\lambda_1 & 0 & 0 \\
0 & \lambda_2 & 0 \\
0 & 0 & \lambda_3
\end{pmatrix}.
\]
Because $\alpha, \beta, \gamma$ are invariant under rotations by part (i), we may evaluate them directly in this diagonal frame (where $T'_{ik} = 0$ for $i \neq k$, so only diagonal terms survive in the sums):
\[
\begin{aligned}
\alpha &= T'_{ii} = \lambda_1 + \lambda_2 + \lambda_3, \\
\beta &= T'_{ik}\,T'_{ki} = \sum_{i=1}^3 (T'_{ii})^2 = \lambda_1^2 + \lambda_2^2 + \lambda_3^2, \\
\gamma &= T'_{ik}\,T'_{km}\,T'_{mi} = \sum_{i=1}^3 (T'_{ii})^3 = \lambda_1^3 + \lambda_2^3 + \lambda_3^3.
\end{aligned}
\]

#### Deduction of the characteristic cubic

The eigenvalues $\lambda_1, \lambda_2, \lambda_3$ are the roots of the characteristic equation $\det(\lambda I - T) = 0$, which factorises as
\[
(\lambda - \lambda_1)(\lambda - \lambda_2)(\lambda - \lambda_3) = \lambda^3 - e_1 \lambda^2 + e_2 \lambda - e_3 = 0,
\]
where the elementary symmetric polynomials of the eigenvalues are
\[
e_1 = \lambda_1 + \lambda_2 + \lambda_3, \qquad
e_2 = \lambda_1\lambda_2 + \lambda_2\lambda_3 + \lambda_3\lambda_1, \qquad
e_3 = \lambda_1\lambda_2\lambda_3.
\]
We express $e_1, e_2, e_3$ in terms of $\alpha, \beta, \gamma$:

1. Immediately,
   \[
   e_1 = \lambda_1 + \lambda_2 + \lambda_3 = \alpha.
   \]

2. Squaring $\alpha = \lambda_1 + \lambda_2 + \lambda_3$ gives
   \[
   \alpha^2 = (\lambda_1^2 + \lambda_2^2 + \lambda_3^2) + 2(\lambda_1\lambda_2 + \lambda_2\lambda_3 + \lambda_3\lambda_1) = \beta + 2e_2,
   \]
   and hence
   \[
   e_2 = \tfrac{1}{2}(\alpha^2 - \beta).
   \]

3. Each eigenvalue $\lambda_r$ ($r = 1, 2, 3$) satisfies $\lambda_r^3 - e_1\lambda_r^2 + e_2\lambda_r - e_3 = 0$. Summing this equation over $r = 1, 2, 3$ yields
   \[
   \sum_{r=1}^3 \lambda_r^3 - e_1 \sum_{r=1}^3 \lambda_r^2 + e_2 \sum_{r=1}^3 \lambda_r - 3 e_3 = 0,
   \]
   that is,
   \[
   \gamma - e_1 \beta + e_2 \alpha - 3 e_3 = 0.
   \]
   Substituting $e_1 = \alpha$ and $e_2 = \tfrac{1}{2}(\alpha^2 - \beta)$:
   \[
   3 e_3 = \gamma - \alpha\beta + \tfrac{1}{2}\alpha(\alpha^2 - \beta) = \tfrac{1}{2}(\alpha^3 - 3\alpha\beta + 2\gamma),
   \]
   so
   \[
   e_3 = \tfrac{1}{6}(\alpha^3 - 3\alpha\beta + 2\gamma).
   \]

Substituting these expressions for $e_1, e_2, e_3$ back into $\lambda^3 - e_1\lambda^2 + e_2\lambda - e_3 = 0$ gives the required cubic equation for the eigenvalues:
\[
\lambda^3 - \alpha\lambda^2 + \tfrac{1}{2}(\alpha^2 - \beta)\lambda - \tfrac{1}{6}(\alpha^3 - 3\alpha\beta + 2\gamma) = 0.
\]

## Final Answers

- **(i)** Direct substitution of $T'_{ik} = \ell_{ip}\ell_{kq}T_{pq}$ and orthogonality $\ell_{ip}\ell_{iq} = \delta_{pq}$ yields $\alpha' = \alpha$, $\beta' = \beta$, and $\gamma' = \gamma$ for any rank-$2$ Cartesian tensor $T_{ik}$.
- **(ii)** In terms of the eigenvalues $\lambda_1, \lambda_2, \lambda_3$ of a real symmetric tensor $T_{ik}$:
  \[
  \alpha = \lambda_1 + \lambda_2 + \lambda_3, \qquad
  \beta = \lambda_1^2 + \lambda_2^2 + \lambda_3^2, \qquad
  \gamma = \lambda_1^3 + \lambda_2^3 + \lambda_3^3.
  \]
  With $e_1 = \alpha$, $e_2 = \tfrac{1}{2}(\alpha^2-\beta)$, and $e_3 = \tfrac{1}{6}(\alpha^3 - 3\alpha\beta + 2\gamma)$, the characteristic cubic for the eigenvalues is
  \[
  \lambda^3 - \alpha\lambda^2 + \tfrac{1}{2}(\alpha^2 - \beta)\lambda - \tfrac{1}{6}(\alpha^3 - 3\alpha\beta + 2\gamma) = 0.
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-4_qx10_check.py` (exit code 0).

```text
OK: Part (i) verified on 60 arbitrary non-symmetric tensors (50 SO(3), 10 O(3)\SO(3)) — max |da|=2.58e-14, |db|=2.34e-13, |dg|=3.30e-12.
OK: Part (ii) exact rational identities verified on 7 eigenvalue triples (incl. repeated/zero).
OK: 30 dense rotated symmetric tensors verified — max root err=6.39e-14, max det err=1.14e-13, max Cayley-Hamilton err=7.82e-14.
```

## Confidence & Edge-Case Notes

- **Confidence**: `high`.
- **Phase B watchpoints**:
  - In part (i), check that the solution does NOT mistakenly assume symmetry of $T_{ik}$ (and does NOT confuse $\beta = T_{ik}T_{ki} = \operatorname{tr}(T^2)$ with $T_{ik}T_{ik} = \operatorname{tr}(T T^T)$). Also check that distinct dummy indices are used cleanly when expanding $T'_{ik}T'_{ki}$ and $T'_{ik}T'_{km}T'_{mi}$ (no index appearing more than twice!).
  - In part (ii), check that diagonalisation via a *proper* rotation (choosing a right-handed orthonormal eigenbasis) is justified before appealing to part (i)'s rotational invariance.
  - Check the derivation of $e_3 = \frac{1}{6}(\alpha^3 - 3\alpha\beta + 2\gamma)$ is actually proved rather than merely asserted by quoting Newton's identities.
