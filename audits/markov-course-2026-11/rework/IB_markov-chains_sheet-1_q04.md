---
id: IB/markov-chains/sheet-1/q04
slug: IB_markov-chains_sheet-1_q04
source_archive: private/raw/sheets/IB/markov-chains/2025-26/sheet-1.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

- Checked word-for-word and symbol-for-symbol against Question 4 on page 1 of `private/raw/sheets/IB/markov-chains/2025-26/sheet-1.pdf` (and the 200 DPI render `IB_markov-chains_sheet-1_q03_p1.png`).
- Question 4 lies mid-page and is followed immediately by Question 5, so nothing is truncated across pages.
- **Result:** Exact match (`transcription: ok`). Purely typographic note: `\Z` in `question.tex` is the project's shared macro for $\mathbb{Z}$.

## Full Solution

Let $S_0 = 0$ and write $S_n = \sum_{m=1}^{n} \xi_m$ for $n \ge 1$, where $(\xi_m)_{m\ge 1}$ are i.i.d. random variables with
$$
\mathbb{P}(\xi_m = +1) = p, \qquad \mathbb{P}(\xi_m = -1) = q := 1 - p, \qquad p \in [0, 1].
$$
Both $X_n = |S_n|$ and $Y_n = M_n - S_n$ take values in the countable state space $\mathbb{Z}_{\ge 0} = \{0, 1, 2, \dots\}$ with $X_0 = 0$ and $Y_0 = 0$.

---

### Part 1: $X_n = |S_n|$ is a Markov chain, and its transition probabilities

We will show that $(X_n)_{n\ge 0}$ is a **time-homogeneous** Markov chain on $\mathbb{Z}_{\ge 0}$ with initial state $X_0 = 0$ and transition probabilities
$$
p^X_{0, 1} = 1, \qquad
p^X_{k, k+1} = \frac{p^{k+1} + q^{k+1}}{p^k + q^k}, \qquad
p^X_{k, k-1} = \frac{pq(p^{k-1} + q^{k-1})}{p^k + q^k} \qquad (k \ge 1),
$$
and $p^X_{k, j} = 0$ for all other $j \in \mathbb{Z}_{\ge 0}$.

First verify that $P^X$ is a stochastic matrix: for $k = 0$, $p^X_{0,1} = 1$; for every $k \ge 1$, $p^k + q^k > 0$ (since $p, q \ge 0$ and $p + q = 1$), the entries are non-negative, and
$$
p^X_{k, k+1} + p^X_{k, k-1}
= \frac{p^{k+1} + q^{k+1} + p^k q + p q^k}{p^k + q^k}
= \frac{(p + q)(p^k + q^k)}{p^k + q^k} = 1.
$$

#### Degenerate cases $p = 1$ ($q = 0$) and $p = 0$ ($q = 1$)
- If $p = 1$ ($q = 0$), then $S_n = n$ almost surely for all $n \ge 0$, so $X_n = n$ almost surely.
- If $p = 0$ ($q = 1$), then $S_n = -n$ almost surely for all $n \ge 0$, so $X_n = n$ almost surely.

In both cases $(X_n)_{n\ge 0}$ is trivially a homogeneous Markov chain with $p^X_{k, k+1} = 1$ for all $k \ge 0$ (and $0$ otherwise), which agrees with the formula above (since for $k \ge 1$, $\frac{1^{k+1}+0}{1^k+0} = 1$ and $\frac{0}{1} = 0$).

#### Main case $p \in (0, 1)$ (so $pq > 0$)
Fix $n \ge 0$ and any history $(x_0, x_1, \dots, x_n) \in \mathbb{Z}_{\ge 0}^{n+1}$ such that the event
$$
A := \{X_0 = x_0, X_1 = x_1, \dots, X_n = x_n\}
$$
has $\mathbb{P}(A) > 0$. (In particular $x_0 = 0$ since $S_0 = 0$.) Set $k := x_n \in \mathbb{Z}_{\ge 0}$.

**Subcase $k = 0$:** On the event $A$, $|S_n| = 0 \implies S_n = 0$, so $X_{n+1} = |S_n + \xi_{n+1}| = |\xi_{n+1}| = 1$ almost surely. Hence
$$
\mathbb{P}(X_{n+1} = 1 \mid A) = 1 = p^X_{0, 1}, \qquad \mathbb{P}(X_{n+1} = j \mid A) = 0 = p^X_{0, j} \quad (j \ne 1).
$$

**Subcase $k \ge 1$:** On the event $A$, $S_n \in \{+k, -k\}$. Define the finite sets of sample paths of $S$:
$$
\mathcal{P}_{\pm}(x_0, \dots, x_n) := \bigl\{ (s_0, \dots, s_n) \in \mathbb{Z}^{n+1} : s_0 = 0,\ |s_m| = x_m \ (0 \le m \le n),\ |s_m - s_{m-1}| = 1 \ (1 \le m \le n),\ s_n = \pm k \bigr\}.
$$
Observe two crucial properties:
1. **Equal path probabilities within each sign:** Let $(s_0, \dots, s_n) \in \mathcal{P}_+(x_0, \dots, x_n)$, and let $u$ and $d$ be the numbers of $+1$ and $-1$ steps among $s_1 - s_0, \dots, s_n - s_{n-1}$. Then $u + d = n$ and $u - d = s_n - s_0 = k$, so
   $$
   u = \frac{n + k}{2}, \qquad d = \frac{n - k}{2}.
   $$
   (In particular, if $\mathcal{P}_+$ is non-empty, $(n+k)/2$ and $(n-k)/2$ are non-negative integers.) Because $\xi_1, \dots, \xi_n$ are independent, every path in $\mathcal{P}_+(x_0, \dots, x_n)$ has the exact same probability:
   $$
   \mathbb{P}(S_0 = s_0, \dots, S_n = s_n) = p^{(n+k)/2} q^{(n-k)/2} = (pq)^{(n-k)/2} p^k.
   $$
   Similarly, every path in $\mathcal{P}_-(x_0, \dots, x_n)$ has $u = (n-k)/2$ and $d = (n+k)/2$, and therefore has probability
   $$
   \mathbb{P}(S_0 = s_0, \dots, S_n = s_n) = p^{(n-k)/2} q^{(n+k)/2} = (pq)^{(n-k)/2} q^k.
   $$
2. **Equal path counts by sign reflection:** The map $\Phi(s_0, s_1, \dots, s_n) = (-s_0, -s_1, \dots, -s_n)$ is a bijection between $\mathcal{P}_+(x_0, \dots, x_n)$ and $\mathcal{P}_-(x_0, \dots, x_n)$: it uses $s_0 = 0 \implies -s_0 = 0$, preserves $|-s_m| = |s_m| = x_m$ and $|(-s_m) - (-s_{m-1})| = |s_m - s_{m-1}| = 1$, sends $s_n = +k$ to $-s_n = -k$, and is its own inverse. Let
   $$
   N := |\mathcal{P}_+(x_0, \dots, x_n)| = |\mathcal{P}_-(x_0, \dots, x_n)|.
   $$

Summing over the disjoint path events in $\mathcal{P}_+$ and $\mathcal{P}_-$ yields
$$
\mathbb{P}(A, S_n = +k) = N (pq)^{(n-k)/2} p^k, \qquad
\mathbb{P}(A, S_n = -k) = N (pq)^{(n-k)/2} q^k,
$$
and therefore
$$
\mathbb{P}(A) = \mathbb{P}(A, S_n = +k) + \mathbb{P}(A, S_n = -k) = N (pq)^{(n-k)/2} (p^k + q^k).
$$
Since $\mathbb{P}(A) > 0$ and $pq > 0$, we have $N(pq)^{(n-k)/2} > 0$, so it cancels identically to give the conditional distribution of the sign of $S_n$:
$$
\mathbb{P}(S_n = +k \mid A) = \frac{p^k}{p^k + q^k}, \qquad \mathbb{P}(S_n = -k \mid A) = \frac{q^k}{p^k + q^k}.
$$
Notably, this depends **only** on $k = x_n$, not on $n$ or $x_0, \dots, x_{n-1}$.

Now compute $\mathbb{P}(X_{n+1} = x_{n+1} \mid A)$. Both events $A \cap \{S_n = +k\}$ and $A \cap \{S_n = -k\}$ are determined by $(S_0, \dots, S_n)$ and are therefore independent of $\xi_{n+1}$. Because $k \ge 1$ and $\xi_{n+1} \in \{-1, +1\}$:
- On $A \cap \{S_n = +k\}$, $X_{n+1} = |k + \xi_{n+1}| = k + \xi_{n+1}$, so $X_{n+1} = k+1 \iff \xi_{n+1} = +1$ and $X_{n+1} = k-1 \iff \xi_{n+1} = -1$.
- On $A \cap \{S_n = -k\}$, $X_{n+1} = |-k + \xi_{n+1}| = k - \xi_{n+1}$, so $X_{n+1} = k+1 \iff \xi_{n+1} = -1$ and $X_{n+1} = k-1 \iff \xi_{n+1} = +1$.

(Check $k = 1$: $|1 + (-1)| = 0 = k-1$ and $|-1 + 1| = 0 = k-1$, so these identities hold for $k = 1$ as well as $k \ge 2$.)

Hence:
$$
\begin{aligned}
\mathbb{P}(X_{n+1} = k+1 \mid A)
&= \frac{\mathbb{P}(A, S_n = +k, \xi_{n+1} = +1) + \mathbb{P}(A, S_n = -k, \xi_{n+1} = -1)}{\mathbb{P}(A)} \\
&= \mathbb{P}(\xi_{n+1} = +1)\,\mathbb{P}(S_n = +k \mid A) + \mathbb{P}(\xi_{n+1} = -1)\,\mathbb{P}(S_n = -k \mid A) \\
&= p \cdot \frac{p^k}{p^k + q^k} + q \cdot \frac{q^k}{p^k + q^k} = \frac{p^{k+1} + q^{k+1}}{p^k + q^k} = p^X_{k, k+1},
\end{aligned}
$$
$$
\begin{aligned}
\mathbb{P}(X_{n+1} = k-1 \mid A)
&= \frac{\mathbb{P}(A, S_n = +k, \xi_{n+1} = -1) + \mathbb{P}(A, S_n = -k, \xi_{n+1} = +1)}{\mathbb{P}(A)} \\
&= q \cdot \frac{p^k}{p^k + q^k} + p \cdot \frac{q^k}{p^k + q^k} = \frac{pq(p^{k-1} + q^{k-1})}{p^k + q^k} = p^X_{k, k-1},
\end{aligned}
$$
and $\mathbb{P}(X_{n+1} = j \mid A) = 0 = p^X_{k, j}$ for all $j \notin \{k-1, k+1\}$.

Since $\mathbb{P}(X_{n+1} = j \mid X_0 = x_0, \dots, X_n = x_n) = p^X_{x_n, j}$ depends only on $x_n$ and $j$ (and not on $n$ or $x_0, \dots, x_{n-1}$), $(X_n)_{n\ge 0}$ is a homogeneous Markov chain with transition matrix $P^X$.

*Special case (symmetric walk $p = q = 1/2$):* $p^X_{0, 1} = 1$, and for all $k \ge 1$, $p^X_{k, k+1} = p^X_{k, k-1} = 1/2$.

*(Note why $S_0 = 0$ is crucial when $p \ne 1/2$: if $S_0 = 1$ deterministically and $p \ne 1/2$, then $\mathbb{P}(X_1 = 2 \mid X_0 = 1) = p$, whereas $\mathbb{P}(X_3 = 2 \mid X_0 = 1, X_1 = 0, X_2 = 1) = p^2 + q^2 \ne p$, so $|S_n|$ is not a Markov chain.)*

---

### Part 2: $Y_n = M_n - S_n$ is a Markov chain

Recall $M_n = \max_{0 \le k \le n} S_k$ and $Y_n = M_n - S_n \in \mathbb{Z}_{\ge 0}$, with $Y_0 = 0$.

For every $n \ge 0$, we have $M_{n+1} = \max(M_n, S_{n+1})$ and $S_{n+1} = S_n + \xi_{n+1}$. Therefore, pointwise:
$$
\begin{aligned}
Y_{n+1}
&= M_{n+1} - S_{n+1} \\
&= \max(M_n, S_{n+1}) - S_{n+1} \\
&= \max(M_n - S_{n+1}, 0) \\
&= \max(M_n - S_n - \xi_{n+1}, 0) \\
&= \max(Y_n - \xi_{n+1}, 0) = (Y_n - \xi_{n+1})^+.
\end{aligned}
$$
Since $S_m = \sum_{r=1}^m \xi_r$ for $0 \le m \le n$, the entire history $(Y_0, Y_1, \dots, Y_n)$ is a deterministic function of $(\xi_1, \dots, \xi_n)$. Because $\xi_1, \xi_2, \dots$ are independent, $\xi_{n+1}$ is independent of $(\xi_1, \dots, \xi_n)$ and hence independent of $(Y_0, Y_1, \dots, Y_n)$.

Consequently, for any $n \ge 0$ and any $y_0, \dots, y_n, y_{n+1} \in \mathbb{Z}_{\ge 0}$ with $\mathbb{P}(Y_0 = y_0, \dots, Y_n = y_n) > 0$,
$$
\begin{aligned}
\mathbb{P}(Y_{n+1} = y_{n+1} \mid Y_0 = y_0, \dots, Y_n = y_n)
&= \mathbb{P}\bigl((Y_n - \xi_{n+1})^+ = y_{n+1} \mid Y_0 = y_0, \dots, Y_n = y_n\bigr) \\
&= \mathbb{P}\bigl((y_n - \xi_{n+1})^+ = y_{n+1} \mid Y_0 = y_0, \dots, Y_n = y_n\bigr) \\
&= \mathbb{P}\bigl((y_n - \xi_{n+1})^+ = y_{n+1}\bigr),
\end{aligned}
$$
where the final step uses the independence of $\xi_{n+1}$ from $(Y_0, \dots, Y_n)$.

Because the right-hand side depends only on $y_n$ and $y_{n+1}$ (and neither on $n$ nor on $y_0, \dots, y_{n-1}$), $(Y_n)_{n\ge 0}$ is a **time-homogeneous** Markov chain on $\mathbb{Z}_{\ge 0}$ with initial state $Y_0 = 0$ and transition probabilities $p^Y_{i, j} = \mathbb{P}\bigl((i - \xi_1)^+ = j\bigr)$ given explicitly by:
- For $i = 0$: since $(0 - 1)^+ = 0$ (prob $p$) and $(0 - (-1))^+ = 1$ (prob $q$),
  $$
  p^Y_{0, 0} = p, \qquad p^Y_{0, 1} = q = 1 - p, \qquad p^Y_{0, j} = 0 \quad (j \notin \{0, 1\}).
  $$
- For $i \ge 1$: since $(i - 1)^+ = i - 1$ (prob $p$) and $(i - (-1))^+ = i + 1$ (prob $q$),
  $$
  p^Y_{i, i-1} = p, \qquad p^Y_{i, i+1} = q = 1 - p, \qquad p^Y_{i, j} = 0 \quad (j \notin \{i-1, i+1\}).
  $$
This holds for every $p \in [0, 1]$ (including the boundary cases $p = 0$ and $p = 1$).

## Final Answers

1. $X_n = |S_n|$ is a homogeneous Markov chain on $\mathbb{Z}_{\ge 0}$ (starting at $X_0 = 0$) with transition probabilities (where $q = 1 - p$, $p \in [0, 1]$):
   $$
   p^X_{0, 1} = 1, \qquad
   p^X_{k, k+1} = \frac{p^{k+1} + q^{k+1}}{p^k + q^k}, \qquad
   p^X_{k, k-1} = \frac{pq(p^{k-1} + q^{k-1})}{p^k + q^k} \qquad (k \ge 1),
   $$
   and $p^X_{k, j} = 0$ otherwise. (For symmetric $p = 1/2$: $p^X_{0, 1} = 1$ and $p^X_{k, k+1} = p^X_{k, k-1} = 1/2$ for $k \ge 1$.)
2. $Y_n = M_n - S_n$ satisfies the recursion $Y_{n+1} = \max(Y_n - \xi_{n+1}, 0)$ with $\xi_{n+1}$ independent of $(Y_0, \dots, Y_n)$, and is therefore a homogeneous Markov chain on $\mathbb{Z}_{\ge 0}$ (starting at $Y_0 = 0$) with transition probabilities:
   $$
   p^Y_{0, 0} = p, \quad p^Y_{0, 1} = q, \qquad p^Y_{i, i-1} = p, \quad p^Y_{i, i+1} = q \quad (i \ge 1),
   $$
   and $0$ otherwise.

## Computational Verification

- Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-1_q04_check.py`
- Command: `python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-1_q04_check.py"` (exit code 0)
- Verbatim stdout:
```text
PASS: X_n = |S_n| Markov property & transition probabilities verified across 328 histories (N=8, 6 values of p).
PASS: Y_n = M_n - S_n Markov property, Lindley recursion & transition probabilities verified across 1036 histories.
PASS: Counterexample with S_0=1, p=1/3 confirmed: P(X_1=2|X_0=1)=1/3 != P(X_3=2|X_0=1,X_1=0,X_2=1)=5/9.
```

## Confidence & Edge-Case Notes

- **Confidence:** High.
- **Subtleties / Phase B comparison watchpoints:**
  1. **Common fatal trap in Part 1:** Many solutions accidentally assume $p = 1/2$ (where $p^X_{k, k\pm 1} = 1/2$) or falsely claim $\mathbb{P}(S_n = +k \mid X_0 = x_0, \dots, X_n = k) = p$ for all $k \ge 1$ (it is actually $\frac{p^k}{p^k + q^k}$, which equals $p$ only at $k = 1$!). Check that the existing solution handles asymmetric $p \ne 1/2$ with the exact formula $\frac{p^{k+1}+q^{k+1}}{p^k+q^k}$ and $\frac{pq(p^{k-1}+q^{k-1})}{p^k+q^k}$.
  2. **Role of $S_0 = 0$:** Without $S_0 = 0$, the reflection bijection $(s_0, \dots, s_n) \mapsto (-s_0, \dots, -s_n)$ fails to preserve the starting state and $|S_n|$ is NOT Markov when $p \ne 1/2$. Watch that the proof uses $S_0 = 0$ explicitly.
  3. **Boundary cases:** Check $k = 0$ (where $p^X_{0, 1} = 1$, NOT $\frac{p+q}{1+1} = 1/2$), $k = 1$ (making sure $|1 + \xi|$ and $|-1 + \xi|$ are verified), and $p \in \{0, 1\}$ (avoiding division by $pq = 0$).
  4. **Direction of transitions in Part 2:** Since $Y_n = M_n - S_n$ *decreases* when $S$ steps up ($\xi_{n+1} = +1$, prob $p$) and *increases* when $S$ steps down ($\xi_{n+1} = -1$, prob $q$), $p^Y_{i, i-1} = p$ and $p^Y_{i, i+1} = q$ (and $p^Y_{0,0} = p$, $p^Y_{0,1} = q$) — easy to flip $p$ and $q$ by accident!
