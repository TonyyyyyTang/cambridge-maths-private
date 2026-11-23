---
id: IB/markov-chains/sheet-2/q11
slug: IB_markov-chains_sheet-2_q11
source_archive: private/raw/sheets/IB/markov-chains/2025-26/sheet-2.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

# Blind Rework: `IB/markov-chains/sheet-2/q11`

## Transcription Check

**Status:** `ok` (exact match).

Compared `question.tex` word-for-word and symbol-for-symbol against the rendered screenshot of page 1 of `private/raw/sheets/IB/markov-chains/2025-26/sheet-2.pdf` (Exercise 11, followed by Exercise 12 on the same page, so no truncation occurs):

- Printed PDF:
  *"11. A fair die is thrown repeatedly. Let $X_n$ denote the sum of the first $n$ throws. Find*
  $$\lim_{n\to\infty} \mathbb{P}[X_n \text{ is a multiple of } 13]$$
  *quoting carefully any general theorems that you use."*
- `question.tex`: exact match in every word, symbol, display equation, and punctuation (including no comma before "quoting").

---

## Full Solution

### 1. Construction of the Markov chain modulo $13$

Let $D_1, D_2, \dots$ be i.i.d. random variables representing the outcomes of successive throws of the fair six-sided die, so that
\[
\mathbb{P}(D_k = r) = \frac{1}{6} \qquad \text{for } r \in \{1, 2, 3, 4, 5, 6\}.
\]
Set $X_0 := 0$ and $X_n := \sum_{k=1}^n D_k$ for $n \ge 1$. For each $n \ge 0$, define $Y_n \in I := \{0, 1, 2, \dots, 12\}$ to be the unique remainder of $X_n$ upon division by $13$ (i.e. $Y_n \equiv X_n \pmod{13}$). Then $Y_0 = 0$, and for every $n \ge 0$,
\[
Y_{n+1} = (Y_n + D_{n+1}) \bmod 13.
\]
Since $D_{n+1}$ is independent of $(D_1, \dots, D_n)$ and hence of $(Y_0, \dots, Y_n)$, for any $i_0, \dots, i_n, j \in I$ with $\mathbb{P}(Y_0 = i_0, \dots, Y_n = i_n) > 0$,
\[
\mathbb{P}(Y_{n+1} = j \mid Y_0 = i_0, \dots, Y_n = i_n) = \mathbb{P}((i_n + D_{n+1}) \bmod 13 = j) = p_{i_n, j},
\]
where the $13 \times 13$ transition matrix $P = (p_{ij})_{i, j \in I}$ is given by
\[
p_{ij} = \begin{cases}
\frac{1}{6}, & \text{if } (j - i) \bmod 13 \in \{1, 2, 3, 4, 5, 6\}, \\
0, & \text{otherwise.}
\end{cases}
\]
Because the residues $1, 2, 3, 4, 5, 6$ are distinct modulo $13$, each row of $P$ contains six entries equal to $\frac{1}{6}$ and seven entries equal to $0$, so $P$ is a stochastic matrix and $(Y_n)_{n \ge 0}$ is $\mathrm{Markov}(\delta_0, P)$. Moreover,
\[
\mathbb{P}[X_n \text{ is a multiple of } 13] = \mathbb{P}(Y_n = 0 \mid Y_0 = 0) = p_{00}^{(n)}.
\]

### 2. Statement of the Convergence to Equilibrium Theorem

We shall use the following standard theorem of IB Markov Chains:

> **Theorem (Convergence to equilibrium):** Let $P$ be an **irreducible** and **aperiodic** stochastic matrix on a countable state space $I$, and suppose that $P$ has an **invariant distribution** $\pi$ (equivalently, $P$ is irreducible, aperiodic, and positive recurrent). Let $(Y_n)_{n \ge 0}$ be $\mathrm{Markov}(\lambda, P)$ for any initial distribution $\lambda$ on $I$. Then for every $j \in I$,
> \[
> \mathbb{P}(Y_n = j) \to \pi_j \qquad \text{as } n \to \infty.
> \]
> In particular, $p_{ij}^{(n)} \to \pi_j$ as $n \to \infty$ for all $i, j \in I$.

(Notice that irreducibility also implies that the invariant distribution $\pi$, if it exists, is unique.)

### 3. Verification of all hypotheses

1. **Irreducibility:**
   For every $i \in I$, we have $p_{i, (i+1) \bmod 13} = \frac{1}{6} > 0$ (corresponding to throwing a $1$). Thus for any $i, j \in I$, setting $k := (j - i) \bmod 13 \in \{0, 1, \dots, 12\}$ gives
   \[
   p_{ij}^{(k)} \ge \left(\frac{1}{6}\right)^k > 0.
   \]
   (Alternatively, any integer in $\{3, 4, \dots, 18\}$ is a sum of three numbers in $\{1, \dots, 6\}$, and $\{3, \dots, 18\}$ covers all $13$ residue classes modulo $13$, so in fact $p_{ij}^{(3)} > 0$ for all $i, j \in I$.) Hence $i \to j$ for all $i, j \in I$, so $P$ is **irreducible**.

2. **Aperiodicity:**
   Recall that the period of a state $i \in I$ is $d_i := \gcd\{n \ge 1 : p_{ii}^{(n)} > 0\}$. Consider state $0$ (in fact the same argument applies directly to every $i \in I$):
   - Throwing $(6, 6, 1)$ gives sum $13 \equiv 0 \pmod{13}$ in $3$ steps, so $p_{00}^{(3)} \ge \left(\frac{1}{6}\right)^3 > 0$.
   - Throwing $(6, 5, 1, 1)$ gives sum $13 \equiv 0 \pmod{13}$ in $4$ steps, so $p_{00}^{(4)} \ge \left(\frac{1}{6}\right)^4 > 0$.
   (Note: $p_{00}^{(1)} = p_{00}^{(2)} = 0$ since the sum of one or two die throws lies in $\{1, \dots, 12\}$.)
   Since $\gcd(3, 4) = 1$, we have $d_0 = 1$. Since $P$ is irreducible and periodicity is a class property (or directly by translation invariance of $P$, giving $d_i = 1$ for every $i \in I$), $P$ is **aperiodic**.

3. **Existence of an invariant distribution:**
   We observe that $P$ is **doubly stochastic**: for every $j \in I$,
   \[
   \sum_{i \in I} p_{ij} = \sum_{r=1}^6 p_{(j-r)\bmod 13,\, j} = \sum_{r=1}^6 \frac{1}{6} = 1,
   \]
   since for each $r \in \{1, \dots, 6\}$ there is a unique $i = (j - r) \bmod 13 \in I$ with $(j - i) \bmod 13 = r$, and these six states $i$ are distinct.

   Define the uniform probability distribution $\pi$ on $I$ by $\pi_i := \frac{1}{13}$ for all $i \in I$. Clearly $\pi_i \ge 0$ and $\sum_{i \in I} \pi_i = 1$. For every $j \in I$,
   \[
   (\pi P)_j = \sum_{i \in I} \pi_i p_{ij} = \frac{1}{13} \sum_{i \in I} p_{ij} = \frac{1}{13} \cdot 1 = \frac{1}{13} = \pi_j.
   \]
   Hence $\pi = \left(\frac{1}{13}, \dots, \frac{1}{13}\right)$ is an (and by irreducibility, the unique) **invariant distribution** for $P$.

### 4. Conclusion

Since $P$ is irreducible and aperiodic and has invariant distribution $\pi$, the Convergence to Equilibrium Theorem applies and yields
\[
\lim_{n \to \infty} \mathbb{P}[X_n \text{ is a multiple of } 13] = \lim_{n \to \infty} p_{00}^{(n)} = \pi_0 = \frac{1}{13}. \qquad \square
\]

---

## Final Answers

\[
\lim_{n\to\infty}\mathbb{P}[X_n \text{ is a multiple of } 13] = \frac{1}{13}.
\]

---

## Computational Verification

Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_q11_check.py`

Command: `python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_q11_check.py"` (exit code 0)

Verbatim stdout:
```text
[PASS] 13x13 transition matrix P is doubly stochastic.
[PASS] Irreducibility & aperiodicity verified: p_00^(1)=0, p_00^(2)=0, p_00^(3)=7/72, p_00^(4)=35/324 (gcd(3,4)=1; all entries of P^3 and P^4 are strictly positive).
[PASS] pi = (1/13,...,1/13) satisfies pi*P = pi and rank(P^T - I) = 12 = 12 (unique).
[PASS] Brute-force die enumeration matches p_00^(n) exactly for n = 1..5.
       Convergence of exact p_00^(n) to 1/13 = 0.076923076923...:
         n= 3: p_00^(n) = 0.097222222222222, |err| = 2.030e-02
         n= 4: p_00^(n) = 0.108024691358025, |err| = 3.110e-02
         n= 5: p_00^(n) = 0.063014403292181, |err| = 1.391e-02
         n=10: p_00^(n) = 0.075562244915240, |err| = 1.361e-03
         n=20: p_00^(n) = 0.076851417181357, |err| = 7.166e-05
         n=30: p_00^(n) = 0.076925191577829, |err| = 2.115e-06
         n=50: p_00^(n) = 0.076923075480110, |err| = 1.443e-09
         n=80: p_00^(n) = 0.076923076923055, |err| = 2.240e-14
[PASS] Fourier cross-check: max_{k!=0} |lambda_k| = 0.691352 < 1; DFT matches exact P^20 to 0.00e+00.
ALL CHECKS PASSED.
```

---

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Subtleties / pitfalls to watch for in Phase B:**
  1. **Quoting the right theorem with all hypotheses:** The question explicitly asks for *"quoting carefully any general theorems that you use"*. One must quote the **Convergence to Equilibrium Theorem** (NOT merely the Ergodic Theorem, which only gives time-average/Cesàro convergence) and explicitly verify **all three hypotheses**: irreducibility, **aperiodicity**, and existence of an invariant distribution $\pi$.
  2. **Aperiodicity proof pitfall:** $p_{00}^{(1)} = 0$ and $p_{00}^{(2)} = 0$ (since $2 \le D_1 + D_2 \le 12 < 13$), so one must use e.g. $p_{00}^{(3)} > 0$ and $p_{00}^{(4)} > 0$ (with $\gcd(3, 4) = 1$), or $p_{00}^{(3)} > 0$ and $p_{00}^{(13)} > 0$. Watch out for any careless claim that $p_{00}^{(2)} > 0$ or $p_{ii} > 0$.
  3. **Doubly stochastic justification:** Showing $\sum_i p_{ij} = 1$ explicitly verifies that the uniform distribution $\pi_i = \frac{1}{13}$ is invariant without guessing.
