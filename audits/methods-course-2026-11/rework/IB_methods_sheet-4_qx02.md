---
id: IB/methods/sheet-4/qx02
slug: IB_methods_sheet-4_qx02
source_archive: private/raw/sheets/IB/methods/2000-01/sheet-4.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

Checked word for word and symbol for symbol against `IB_methods_sheet-4_qx02_p1.png` (Michaelmas term 2000, *METHODS — EXAMPLES IV*, Dr E.P.S. Shellard, Course C10—No. B9, page 1, Question 2). Question 2 sits well above the bottom of page 1, so no text spills onto page 2.

- **Result**: `ok` (exact mathematical and verbal match).
- **Purely typographic notes**:
  - The underlined italic title `Tensor transformations.` is encoded as `\emph{Tensor transformations.}` (standard restricted-LaTeX convention).
  - Inside the bracketed `[Answer: ...]`, the words `or` and `in the standard notation.` are italicised on the printed sheet and encoded with `\text{...}` in `question.tex`. All matrix entries, commas, subscripts (`\ell_{i1}\ell_{j1}`), and wording match verbatim.

## Full Solution

### Setup: Frames and Rotation Matrices

Let the three right-handed Cartesian frames of reference in $\mathbb{R}^3$ be denoted by Frame $1$ ($n=0$), Frame $2$ ($n=1$), and Frame $3$ ($n=2$), with orthonormal basis vectors $\{\mathbf{e}^{(n+1)}_1, \mathbf{e}^{(n+1)}_2, \mathbf{e}^{(n+1)}_3\}$ for $n=0,1,2$.

We are given that the $i^{\text{th}}$ axis of Frame $1$ coincides with the $(i+n)^{\text{th}}$ axis (modulo $3$, taking residues in $\{1,2,3\}$) of Frame $(n+1)$:
\[
\mathbf{e}^{(1)}_i = \mathbf{e}^{(n+1)}_{(i+n)\bmod 3} \qquad (i=1,2,3;\; n=0,1,2).
\]
Explicitly:
- **Frame $1$ ($n=0$)**: $\mathbf{e}^{(1)}_1, \mathbf{e}^{(1)}_2, \mathbf{e}^{(1)}_3$.
- **Frame $2$ ($n=1$)**: $\mathbf{e}^{(2)}_2 = \mathbf{e}^{(1)}_1$, $\mathbf{e}^{(2)}_3 = \mathbf{e}^{(1)}_2$, $\mathbf{e}^{(2)}_1 = \mathbf{e}^{(1)}_3$.
- **Frame $3$ ($n=2$)**: $\mathbf{e}^{(3)}_3 = \mathbf{e}^{(1)}_1$, $\mathbf{e}^{(3)}_1 = \mathbf{e}^{(1)}_2$, $\mathbf{e}^{(3)}_2 = \mathbf{e}^{(1)}_3$.

In the standard direction-cosine convention for transforming from Frame $1$ (unprimed) to a primed Cartesian frame $\{\mathbf{e}'_i\}$, the transformation matrix $L = (\ell_{ij})$ has entries
\[
\ell_{ij} = \mathbf{e}'_i \cdot \mathbf{e}^{(1)}_j.
\]
Thus $L$ is orthogonal ($L L^T = L^T L = I$, i.e. $\ell_{ip}\ell_{jp} = \delta_{ij}$ and $\ell_{pi}\ell_{pj} = \delta_{ij}$), vectors transform as $v'_i = \ell_{ij} v_j$, and a rank-$2$ Cartesian tensor transforms as
\[
T'_{ij} = \ell_{ip}\,\ell_{jq}\,T_{pq}, \qquad \text{or in matrix form} \qquad T' = L T L^T.
\]

For Frame $2$ and Frame $3$ relative to Frame $1$, the direction-cosine matrices $L^{(2)}_{ij} = \mathbf{e}^{(2)}_i \cdot \mathbf{e}^{(1)}_j$ and $L^{(3)}_{ij} = \mathbf{e}^{(3)}_i \cdot \mathbf{e}^{(1)}_j$ are:
\[
L^{(2)} = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}, \qquad
L^{(3)} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}.
\]
Both are proper rotations: $L^{(2)} (L^{(2)})^T = L^{(3)} (L^{(3)})^T = I$ and $\det L^{(2)} = \det L^{(3)} = +1$ (cyclic permutations of axes preserve right-handedness, corresponding to rotations by $\pm 2\pi/3$ about $(1,1,1)$). Note also that $L^{(3)} = (L^{(2)})^2 = (L^{(2)})^T$.

---

### Part (i)

The given $3\times 3$ arrays in the three frames are
\[
A^{(1)} = \begin{pmatrix}1&0&0\\0&1&0\\0&0&1\end{pmatrix}, \qquad
A^{(2)} = \begin{pmatrix}1&0&0\\0&1&0\\0&0&0\end{pmatrix}, \qquad
A^{(3)} = \begin{pmatrix}1&0&0\\0&0&0\\0&0&0\end{pmatrix}.
\]
Suppose for contradiction that this physical entity is a rank-$2$ Cartesian tensor. In Frame $1$, its components are $A^{(1)}_{pq} = \delta_{pq}$ (the identity matrix $I$). Under the proper rotation $L^{(2)}$ to Frame $2$, the tensor transformation law requires its components in Frame $2$ to be
\[
\bigl(A^{(2)}_{\text{tensor}}\bigr)_{ij} = \ell^{(2)}_{ip}\,\ell^{(2)}_{jq}\,\delta_{pq} = \ell^{(2)}_{ip}\,\ell^{(2)}_{jp} = \delta_{ij},
\]
by orthogonality of $L^{(2)}$ (equivalently $L^{(2)} I (L^{(2)})^T = L^{(2)} (L^{(2)})^T = I$). Likewise, in Frame $3$ its components would have to be $\delta_{ij} = I$ (indeed $\delta_{ij}$ is isotropic, taking the form $I$ in every Cartesian frame).

However, the observed components are $A^{(2)} = \operatorname{diag}(1,1,0) \neq I$ and $A^{(3)} = \operatorname{diag}(1,0,0) \neq I$. (Equivalently, the trace $T_{ii}$ of any rank-$2$ Cartesian tensor is invariant under orthogonal transformations, since $T'_{ii} = \ell_{ip}\ell_{iq}T_{pq} = \delta_{pq}T_{pq} = T_{pp}$, whereas here $\operatorname{tr}(A^{(1)}) = 3$, $\operatorname{tr}(A^{(2)}) = 2$, and $\operatorname{tr}(A^{(3)}) = 1$ are unequal.)

Furthermore, because $\det L^{(2)} = \det L^{(3)} = +1$, even a rank-$2$ pseudotensor would obey the identical transformation law under these proper rotations. Hence the entity cannot be a tensor.

---

### Part (ii)

Now consider the entity with respective components
\[
T^{(1)} = \begin{pmatrix}1&0&0\\0&0&0\\0&0&0\end{pmatrix}, \qquad
T^{(2)} = \begin{pmatrix}0&0&0\\0&1&0\\0&0&0\end{pmatrix}, \qquad
T^{(3)} = \begin{pmatrix}0&0&0\\0&0&0\\0&0&1\end{pmatrix}.
\]
In Frame $1$, the components are
\[
T^{(1)}_{pq} = \delta_{p1}\delta_{q1}.
\]
Under the rank-$2$ tensor transformation law from Frame $1$ to Frame $2$, using $\ell^{(2)}_{i1} = \mathbf{e}^{(2)}_i \cdot \mathbf{e}^{(1)}_1 = \delta_{i2}$ (the first column of $L^{(2)}$ is $(0,1,0)^T$):
\[
\bigl(L^{(2)} T^{(1)} (L^{(2)})^T\bigr)_{ij}
= \ell^{(2)}_{ip}\,\ell^{(2)}_{jq}\,\delta_{p1}\delta_{q1}
= \ell^{(2)}_{i1}\,\ell^{(2)}_{j1}
= \delta_{i2}\delta_{j2}
= \begin{pmatrix}0&0&0\\0&1&0\\0&0&0\end{pmatrix}_{ij}
= T^{(2)}_{ij}.
\]
Similarly, from Frame $1$ to Frame $3$, using $\ell^{(3)}_{i1} = \mathbf{e}^{(3)}_i \cdot \mathbf{e}^{(1)}_1 = \delta_{i3}$ (the first column of $L^{(3)}$ is $(0,0,1)^T$):
\[
\bigl(L^{(3)} T^{(1)} (L^{(3)})^T\bigr)_{ij}
= \ell^{(3)}_{ip}\,\ell^{(3)}_{jq}\,\delta_{p1}\delta_{q1}
= \ell^{(3)}_{i1}\,\ell^{(3)}_{j1}
= \delta_{i3}\delta_{j3}
= \begin{pmatrix}0&0&0\\0&0&0\\0&0&1\end{pmatrix}_{ij}
= T^{(3)}_{ij}.
\]
Because the transformation matrices form a group, consistency between Frame $2$ and Frame $3$ follows immediately: the rotation from Frame $2$ to Frame $3$ is $L^{(3\leftarrow 2)} = L^{(3)}(L^{(2)})^T = L^{(2)}$, and
\[
L^{(3\leftarrow 2)} T^{(2)} (L^{(3\leftarrow 2)})^T
= L^{(3)}(L^{(2)})^T \bigl(L^{(2)} T^{(1)} (L^{(2)})^T\bigr) L^{(2)} (L^{(3)})^T
= L^{(3)} T^{(1)} (L^{(3)})^T = T^{(3)}.
\]
Thus the components in all three frames are mutually consistent with the rank-$2$ Cartesian tensor transformation law. Since verifying the transformation law for only three frames is necessary but not sufficient (a tensor must transform correctly under *all* rotations $L \in SO(3)$), this shows that the entity **could** be a rank-$2$ tensor (specifically, it is consistent with the dyadic/outer-product tensor $T = \mathbf{e}^{(1)}_1 \otimes \mathbf{e}^{(1)}_1$, whose components in any frame are $u_i u_j$ where $\mathbf{u} = \mathbf{e}^{(1)}_1$).

#### Components in an arbitrary Cartesian frame

Assume now that the entity is indeed a rank-$2$ Cartesian tensor, and consider an arbitrary Cartesian frame $\{\mathbf{e}'_1, \mathbf{e}'_2, \mathbf{e}'_3\}$ whose axes have direction cosines $(\lambda, \mu, \nu)$ with respect to the $1$-axis of Frame $1$, i.e.
\[
\ell_{11} = \mathbf{e}'_1 \cdot \mathbf{e}^{(1)}_1 = \lambda, \qquad
\ell_{21} = \mathbf{e}'_2 \cdot \mathbf{e}^{(1)}_1 = \mu, \qquad
\ell_{31} = \mathbf{e}'_3 \cdot \mathbf{e}^{(1)}_1 = \nu,
\]
with $\lambda^2 + \mu^2 + \nu^2 = \ell_{i1}\ell_{i1} = 1$.

Applying the rank-$2$ tensor transformation law from Frame $1$ to this arbitrary frame:
\[
T'_{ij} = \ell_{ip}\,\ell_{jq}\,T^{(1)}_{pq} = \ell_{ip}\,\ell_{jq}\,\delta_{p1}\delta_{q1} = \ell_{i1}\,\ell_{j1}.
\]
Notice that the result depends solely on the first column $(\ell_{11}, \ell_{21}, \ell_{31}) = (\lambda, \mu, \nu)$ of $L$ (the remaining two columns $\ell_{i2}, \ell_{i3}$ do not enter because $T^{(1)}_{pq} = 0$ unless $p=q=1$). Writing out the $3\times 3$ matrix explicitly:
\[
T' = \begin{pmatrix}
\lambda^2 & \lambda\mu & \lambda\nu \\
\lambda\mu & \mu^2 & \mu\nu \\
\lambda\nu & \mu\nu & \nu^2
\end{pmatrix}.
\]
As consistency checks: $T'$ is symmetric, $\operatorname{tr}(T') = \lambda^2+\mu^2+\nu^2 = 1 = \operatorname{tr}(T^{(1)})$, $(T')^2 = T'$ (rank-$1$ orthogonal projection onto the $1$-axis of Frame $1$), and specialising $(\lambda,\mu,\nu)$ to $(1,0,0)$, $(0,1,0)$, $(0,0,1)$ recovers $T^{(1)}, T^{(2)}, T^{(3)}$ respectively.

## Final Answers

- **(i)** The entity cannot be a rank-$2$ Cartesian tensor: transforming $A^{(1)} = I$ by the proper rotations $L^{(2)}$ and $L^{(3)}$ yields $I$ in both Frame $2$ and Frame $3$, contradicting $A^{(2)} = \operatorname{diag}(1,1,0)$ and $A^{(3)} = \operatorname{diag}(1,0,0)$ (and the traces $3, 2, 1$ fail to be invariant).
- **(ii)** Transforming $T^{(1)}_{pq} = \delta_{p1}\delta_{q1}$ by $L^{(2)}$ and $L^{(3)}$ reproduces $T^{(2)}$ and $T^{(3)}$ exactly, so the entity *could* be a rank-$2$ tensor (consistency across three frames is necessary but not sufficient to establish tensoriality for all frames). Assuming it is a tensor, its components in an arbitrary Cartesian frame whose axes have direction cosines $(\lambda,\mu,\nu)$ relative to the $1$-axis of Frame $1$ are
  \[
  T'_{ij} = \ell_{i1}\ell_{j1}
  = \begin{pmatrix}
  \lambda^2 & \lambda\mu & \lambda\nu \\
  \lambda\mu & \mu^2 & \mu\nu \\
  \lambda\nu & \mu\nu & \nu^2
  \end{pmatrix}
  \qquad (\lambda^2 + \mu^2 + \nu^2 = 1).
  \]

## Computational Verification

Script: `private/audits/methods-course-2026-11/rework/IB_methods_sheet-4_qx02_check.py` (exit code 0).

```text
OK: Frame rotation matrices L^(1), L^(2), L^(3) are exact SO(3) cyclic permutations.
OK: Part (i) verified — L A^(1) L^T = I != A^(2), A^(3); traces 3, 2, 1 disagree.
OK: Part (ii) consistency across frames 1, 2, 3 verified exactly.
OK: 50 random SO(3) frames verified — max formula err=0.00e+00, max trace err=2.22e-16, max idempotence err=2.22e-16.
```

## Confidence & Edge-Case Notes

- **Confidence**: `high`.
- **Phase B watchpoints**:
  - Verify the index convention for $L^{(2)}$ and $L^{(3)}$! A common error is transposing $L^{(2)}$ and $L^{(3)}$ (swapping Frame $2$ and Frame $3$): here "the $i^{\text{th}}$ axis of the first frame coincides with the $(i+n)^{\text{th}}$ axis of the $(n+1)^{\text{th}}$ frame" means $\mathbf{e}^{(1)}_1 = \mathbf{e}^{(2)}_2$ (so $L^{(2)}_{21} = 1$, first column $(0,1,0)^T$) and $\mathbf{e}^{(1)}_1 = \mathbf{e}^{(3)}_3$ (so $L^{(3)}_{31} = 1$, first column $(0,0,1)^T$). Getting the transpose wrong would swap $T^{(2)}$ and $T^{(3)}$ in part (ii).
  - Verify that the existing solution explains *why* part (ii) only shows the entity *could* be a tensor (three frames are necessary, not sufficient; tensor transformation must hold for all rotations in $SO(3)$).
  - Verify that $\ell_{ij} = \mathbf{e}'_i \cdot \mathbf{e}_j$ is stated clearly so $\ell_{i1}\ell_{j1}$ matches the direction-cosine matrix convention (if someone uses the inverse/transposed convention $R_{ij} = \mathbf{e}_i \cdot \mathbf{e}'_j$, the formula would be $R_{1i}R_{1j}$ instead of $\ell_{i1}\ell_{j1}$).
