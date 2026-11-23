---
id: IB/markov-chains/sheet-2/q15
slug: IB_markov-chains_sheet-2_q15
source_archive: private/raw/sheets/IB/markov-chains/2025-26/sheet-2.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

# Blind Rework: `IB/markov-chains/sheet-2/q15`

## Transcription Check

**Status:** `ok` (exact match).

Compared `question.tex` word-for-word and symbol-for-symbol against the rendered screenshot of page 2 of `private/raw/sheets/IB/markov-chains/2025-26/sheet-2.pdf` (Exercise 15, followed by Exercise 16 on the same page):

- Printed PDF: *"15. A random walk on the set $\{0, 1, 2, \dots, b\}$ has transition matrix given by $p_{00} = 1 - \lambda_0$, $p_{bb} = 1 - \mu_b$, $p_{i,i+1} = \lambda_i$ and $p_{i+1,i} = \mu_{i+1}$ for $0 \le i < b$, where $0 < \lambda_i, \mu_i < 1$ for all $i$, and $\lambda_i + \mu_i = 1$ for $1 \le i < b$. Show that this process is time-reversible in equilibrium."*
- `question.tex`: exact match in every word, symbol, subscript, inequality, and punctuation mark.

---

## Full Solution

### 1. Verification of stochasticity, irreducibility, and uniqueness of equilibrium

Let $I := \{0, 1, 2, \dots, b\}$, where $b \in \mathbb{Z}_{\ge 1}$ (note that $b \ge 1$ is implicit since $0$ and $b$ are distinct boundary states with $p_{00} = 1 - \lambda_0 < 1$ and $p_{bb} = 1 - \mu_b < 1$; the boundary case $b = 1$ is valid, with the interior condition $1 \le i < b$ holding vacuously).

All unspecified entries of $P = (p_{ij})_{i,j \in I}$ are $0$. Checking row sums:
- Row $0$: $p_{00} + p_{01} = (1 - \lambda_0) + \lambda_0 = 1$, with both entries in $(0, 1)$.
- Row $i$ ($1 \le i < b$): $p_{i, i-1} + p_{i, i+1} = \mu_i + \lambda_i = 1$, with both entries in $(0, 1)$ (and $p_{ii} = 0$).
- Row $b$: $p_{b, b-1} + p_{bb} = \mu_b + (1 - \mu_b) = 1$, with both entries in $(0, 1)$.

Hence $P$ is a stochastic matrix on the finite state space $I$.

**Irreducibility (and aperiodicity):** Since $\lambda_i > 0$ for all $0 \le i < b$ and $\mu_{i+1} > 0$ for all $0 \le i < b$, for any $0 \le j < k \le b$ we have
\[
p_{jk}^{(k-j)} \ge \lambda_j \lambda_{j+1} \cdots \lambda_{k-1} > 0, \qquad p_{kj}^{(k-j)} \ge \mu_k \mu_{k-1} \cdots \mu_{j+1} > 0.
\]
Hence every pair of states communicates, so $P$ is **irreducible**. *(In addition, $p_{00} = 1 - \lambda_0 > 0$ and $p_{bb} = 1 - \mu_b > 0$, so $P$ is also aperiodic, though aperiodicity is not needed for time-reversibility in equilibrium.)*

Because $I$ is finite and $P$ is irreducible, the Markov chain is positive recurrent and possesses a **unique** invariant distribution $\pi$ (with $\pi_i > 0$ for all $i \in I$).

### 2. Constructing a distribution in detailed balance with $P$

We seek a probability distribution $\pi = (\pi_0, \pi_1, \dots, \pi_b)$ on $I$ satisfying the **detailed balance equations**
\[
\pi_i p_{ij} = \pi_j p_{ji} \qquad \text{for all } i, j \in I.
\]
- For $i = j$, $\pi_i p_{ii} = \pi_i p_{ii}$ holds identically.
- For $|i - j| \ge 2$, $p_{ij} = 0 = p_{ji}$, so $\pi_i p_{ij} = 0 = \pi_j p_{ji}$ holds identically.
- For adjacent states $\{i, i+1\}$ ($0 \le i < b$), detailed balance requires
  \[
  \pi_i p_{i, i+1} = \pi_{i+1} p_{i+1, i} \iff \pi_i \lambda_i = \pi_{i+1} \mu_{i+1}.
  \]

Since $\mu_{i+1} > 0$ for all $0 \le i < b$, this is equivalent to the recurrence
\[
\pi_{i+1} = \frac{\lambda_i}{\mu_{i+1}} \pi_i \qquad (0 \le i < b).
\]
By induction on $i$, setting $w_0 := 1$ and
\[
w_i := \prod_{k=0}^{i-1} \frac{\lambda_k}{\mu_{k+1}} = \frac{\lambda_0 \lambda_1 \cdots \lambda_{i-1}}{\mu_1 \mu_2 \cdots \mu_i} \qquad (1 \le i \le b),
\]
we have $\pi_i = \pi_0 w_i$ for all $0 \le i \le b$. Because $0 < \lambda_k, \mu_{k+1} < 1$ for every $0 \le k < b$, each $w_i > 0$, and because $I$ is finite, the normalising constant
\[
Z := \sum_{i=0}^b w_i = 1 + \sum_{i=1}^b \prod_{k=0}^{i-1} \frac{\lambda_k}{\mu_{k+1}}
\]
is finite and strictly positive ($Z > 1$). Defining
\[
\pi_i := \frac{w_i}{Z} = \frac{1}{Z} \prod_{k=0}^{i-1} \frac{\lambda_k}{\mu_{k+1}} \qquad (0 \le i \le b)
\]
(with the empty product for $i = 0$ understood as $1$), we obtain a probability distribution on $I$ ($\pi_i > 0$ for all $i$, $\sum_{i=0}^b \pi_i = 1$).

Let us verify explicitly that $\pi$ satisfies detailed balance: for every $0 \le i < b$,
\[
\pi_i p_{i, i+1} = \frac{w_i}{Z} \lambda_i = \frac{w_i \frac{\lambda_i}{\mu_{i+1}}}{Z} \mu_{i+1} = \frac{w_{i+1}}{Z} \mu_{i+1} = \pi_{i+1} p_{i+1, i},
\]
and as noted above $\pi_i p_{ij} = \pi_j p_{ji}$ holds trivially whenever $i = j$ or $|i - j| \ge 2$.

### 3. Invariance, uniqueness, and time-reversibility in equilibrium

Summing the detailed balance equations over $i \in I$ shows that $\pi$ is an invariant distribution for $P$: for every $j \in I$,
\[
(\pi P)_j = \sum_{i=0}^b \pi_i p_{ij} = \sum_{i=0}^b \pi_j p_{ji} = \pi_j \sum_{i=0}^b p_{ji} = \pi_j.
\]
By irreducibility of $P$, $\pi$ is the **unique** invariant distribution of the chain, so the process in equilibrium is uniquely $(X_n)_{n \ge 0} \sim \mathrm{Markov}(\pi, P)$.

Finally, we verify directly that $\mathrm{Markov}(\pi, P)$ is time-reversible. For any integer $N \ge 1$ and any states $i_0, i_1, \dots, i_N \in I$, repeated application of $\pi_i p_{ij} = \pi_j p_{ji}$ yields
\[
\begin{aligned}
\mathbb{P}_\pi(X_0 = i_0, X_1 = i_1, \dots, X_N = i_N)
&= \pi_{i_0} p_{i_0 i_1} p_{i_1 i_2} \cdots p_{i_{N-1} i_N} \\
&= p_{i_1 i_0} \pi_{i_1} p_{i_1 i_2} \cdots p_{i_{N-1} i_N} \\
&= \cdots \\
&= p_{i_1 i_0} p_{i_2 i_1} \cdots p_{i_N i_{N-1}} \pi_{i_N} \\
&= \pi_{i_N} p_{i_N i_{N-1}} \cdots p_{i_1 i_0} \\
&= \mathbb{P}_\pi(X_0 = i_N, X_1 = i_{N-1}, \dots, X_N = i_0).
\end{aligned}
\]
Thus, for every $N \ge 1$, the time-reversed process $(X_{N-n})_{n=0}^N$ has the same distribution as $(X_n)_{n=0}^N$ (equivalently, its transition matrix $\hat{p}_{ji} = \frac{\pi_i}{\pi_j} p_{ij}$ coincides with $P$). Therefore the process is time-reversible in equilibrium. $\square$

---

## Final Answers

- The chain is irreducible (and aperiodic) on $\{0, 1, \dots, b\}$ ($b \ge 1$), with **unique** invariant distribution
  \[
  \pi_i = \frac{1}{Z} \prod_{k=0}^{i-1} \frac{\lambda_k}{\mu_{k+1}} \qquad (0 \le i \le b), \qquad Z = 1 + \sum_{i=1}^b \prod_{k=0}^{i-1} \frac{\lambda_k}{\mu_{k+1}}
  \]
  (with empty product $1$ at $i = 0$).
- This unique equilibrium distribution $\pi$ satisfies detailed balance $\pi_i p_{ij} = \pi_j p_{ji}$ for all $i, j \in \{0, \dots, b\}$, so the process is **time-reversible in equilibrium**.

---

## Computational Verification

Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_q15_check.py`

Command: `python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_q15_check.py"` (exit code 0)

Verbatim stdout:
```text
[PASS] Boundary b=1 (b=1): P^1 > 0 (irreducible & aperiodic), unique pi matches formula, detailed balance holds, all paths N<=4 reversible.
[PASS] Asymmetric b=2 (b=2): P^2 > 0 (irreducible & aperiodic), unique pi matches formula, detailed balance holds, all paths N<=4 reversible.
[PASS] Heterogeneous b=3 (b=3): P^3 > 0 (irreducible & aperiodic), unique pi matches formula, detailed balance holds, all paths N<=3 reversible.
[PASS] Biased b=4 (b=4): P^4 > 0 (irreducible & aperiodic), unique pi matches formula, detailed balance holds, all paths N<=3 reversible.
[PASS] Larger b=7 (b=7): P^7 > 0 (irreducible & aperiodic), unique pi matches formula, detailed balance holds, detailed balance verified (full path check skipped for large b).
ALL CHECKS PASSED.
```

---

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Subtleties / pitfalls to watch for in Phase B:**
  1. **Uniqueness of invariant distribution (irreducibility):** Exhibiting *one* distribution $\pi$ in detailed balance shows that the chain started in $\pi$ is reversible, and proving irreducibility (using $\lambda_i, \mu_{i+1} > 0$) guarantees that $\pi$ is the **unique** invariant distribution — so there is no ambiguity about "in equilibrium".
  2. **Boundary case $b = 1$:** Handled cleanly via the empty-product convention (interior condition $\lambda_i + \mu_i = 1$ for $1 \le i < b$ is vacuous; $b = 0$ is excluded by $0 < \lambda_0, \mu_b < 1$ and stochasticity).
  3. **Index range of $\lambda_i$ and $\mu_i$:** Notice carefully that $\pi_{i+1} = \frac{\lambda_i}{\mu_{i+1}}\pi_i$ uses $\lambda_k$ for $0 \le k \le i-1$ and $\mu_{k+1}$ for $1 \le k+1 \le i$ (watch out for off-by-one index slips like writing $\lambda_k / \mu_k$!).
  4. **Strict positivity $0 < \lambda_i, \mu_i < 1$:** Ensures no division by zero in $\lambda_k / \mu_{k+1}$, guarantees irreducibility ($\lambda_i > 0, \mu_{i+1} > 0$), and guarantees $p_{00} = 1 - \lambda_0 > 0, p_{bb} = 1 - \mu_b > 0$ (hence also aperiodicity, and $P^b$ strictly positive).
