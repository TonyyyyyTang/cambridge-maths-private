---
id: IB/markov-chains/sheet-1/q03
slug: IB_markov-chains_sheet-1_q03
source_archive: private/raw/sheets/IB/markov-chains/2025-26/sheet-1.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

- Checked word-for-word and symbol-for-symbol against Question 3 on page 1 of `private/raw/sheets/IB/markov-chains/2025-26/sheet-1.pdf` (and the 200 DPI render `IB_markov-chains_sheet-1_q03_p1.png`).
- Question 3 lies mid-page and is followed immediately by Question 4, so nothing is truncated across pages.
- **Result:** Exact match (`transcription: ok`). No missing, extra, or altered words, symbols, subscripts, or punctuation.
- *Mild note on the original printed wording:* The printed sheet writes "an unbounded increasing sequence of positive integers" and then asks for the case $n_r = 2r$ (which has $n_0 = 0 \in \mathbb{Z}_{\ge 0}$, non-negative rather than strictly positive). Our proof handles any non-decreasing sequence $(n_r)_{r\ge 0}$ in $\mathbb{Z}_{\ge 0}$, seamlessly covering both strictly positive increasing sequences and $n_r = 2r$.

## Full Solution

### Part 1: $Y_r = X_{n_r}$ is a (possibly inhomogeneous) Markov chain

Let $I$ be the countable state space of the Markov chain $X = (X_n)_{n\ge 0}$. Assume first (the standard IB definition) that $X$ is a time-homogeneous Markov chain with initial distribution $\lambda = (\lambda_i)_{i\in I}$ and transition matrix $P = (p_{ij})_{i,j\in I}$. (We note at the end how the identical argument works if $X$ is itself inhomogeneous.)

Let $(n_r)_{r\ge 0}$ be any non-decreasing sequence of non-negative integers $0 \le n_0 \le n_1 \le n_2 \le \dots$ (in particular, any unbounded strictly increasing sequence of positive integers, as well as $n_r = 2r$). Set $Y_r = X_{n_r}$ for $r \ge 0$, and define
$$
\mu = \lambda P^{n_0}, \qquad Q(r) = P^{n_{r+1} - n_r} \quad (r \ge 0),
$$
where $P^0 = I_{\mathrm{id}}$ (the identity matrix with entries $\delta_{ij}$).

First observe that each $Q(r)$ is a stochastic matrix and $\mu$ is a probability distribution on $I$:
- For any $d \ge 0$, every entry $(P^d)_{ij}$ is non-negative, and by induction on $d$ each row sums to $1$: for $d = 0$ this is immediate; assuming $\sum_{j\in I} (P^d)_{ij} = 1$ for all $i \in I$, Tonelli's theorem (for sums of non-negative terms) gives
  $$
  \sum_{j\in I} (P^{d+1})_{ij} = \sum_{j\in I} \sum_{l\in I} (P^d)_{il} p_{lj} = \sum_{l\in I} (P^d)_{il} \sum_{j\in I} p_{lj} = \sum_{l\in I} (P^d)_{il} = 1.
  $$
- Similarly, $\mu_i = \sum_{j\in I} \lambda_j (P^{n_0})_{ji} \ge 0$ and $\sum_{i\in I} \mu_i = \sum_{j\in I} \lambda_j \sum_{i\in I} (P^{n_0})_{ji} = \sum_{j\in I} \lambda_j = 1$.

We now compute the finite-dimensional distributions of $Y = (Y_r)_{r\ge 0}$. Recall that because $X$ is $\mathrm{Markov}(\lambda, P)$, for every $m \ge 0$ and all $j_0, j_1, \dots, j_m \in I$,
$$
\mathbb{P}(X_0 = j_0, X_1 = j_1, \dots, X_m = j_m) = \lambda_{j_0} \prod_{t=0}^{m-1} p_{j_t j_{t+1}}.
$$
Fix $r \ge 0$ and $i_0, i_1, \dots, i_r \in I$. Assume first that $0 \le n_0 < n_1 < \dots < n_r$ (the strictly increasing case). The event $\{Y_0 = i_0, Y_1 = i_1, \dots, Y_r = i_r\} = \{X_{n_0} = i_0, X_{n_1} = i_1, \dots, X_{n_r} = i_r\}$ is the countable disjoint union of the events $\{X_0 = j_0, X_1 = j_1, \dots, X_{n_r} = j_{n_r}\}$ over all choices of the intermediate indices $j_t \in I$ for $t \in \{0, 1, \dots, n_r\} \setminus \{n_0, n_1, \dots, n_r\}$, with $j_{n_k} = i_k$ fixed for each $k = 0, \dots, r$.

By countable additivity of $\mathbb{P}$ and Tonelli's theorem (factoring a multiple sum of non-negative terms over disjoint sets of summation indices):
$$
\begin{aligned}
\mathbb{P}(Y_0 = i_0, \dots, Y_r = i_r)
&= \sum_{\substack{j_t \in I \text{ for } t \notin \{n_0,\dots,n_r\} \\ j_{n_k} = i_k \text{ for } 0 \le k \le r}} \lambda_{j_0} \prod_{t=0}^{n_r - 1} p_{j_t j_{t+1}} \\
&= \left( \sum_{j_0, \dots, j_{n_0 - 1} \in I} \lambda_{j_0} p_{j_0 j_1} \cdots p_{j_{n_0 - 1} i_0} \right)
   \prod_{k=0}^{r-1} \left( \sum_{j_{n_k + 1}, \dots, j_{n_{k+1} - 1} \in I} p_{i_k j_{n_k + 1}} \cdots p_{j_{n_{k+1} - 1} i_{k+1}} \right).
\end{aligned}
$$
Examining each factor:
- The leading factor equals $\lambda_{i_0} = (\lambda P^0)_{i_0} = \mu_{i_0}$ if $n_0 = 0$, and equals $(\lambda P^{n_0})_{i_0} = \mu_{i_0}$ if $n_0 \ge 1$ by the definition of matrix multiplication.
- For each $k \in \{0, \dots, r-1\}$, with $d_k = n_{k+1} - n_k \ge 1$, the $k$-th factor is $p_{i_k i_{k+1}} = (P^1)_{i_k i_{k+1}}$ when $d_k = 1$ (empty sum of intermediate indices), and by matrix multiplication is $(P^{d_k})_{i_k i_{k+1}} = Q(k)_{i_k i_{k+1}}$ for every $d_k \ge 1$.
- (Note: if weakly increasing $n_{k+1} = n_k$ is also allowed, then $\{X_{n_k}=i_k, X_{n_{k+1}}=i_{k+1}\}$ is empty unless $i_k = i_{k+1}$, contributing a factor $\delta_{i_k i_{k+1}} = (P^0)_{i_k i_{k+1}} = Q(k)_{i_k i_{k+1}}$, so the same formula holds.)

Hence, for all $r \ge 0$ and all $i_0, \dots, i_r \in I$,
$$
\mathbb{P}(Y_0 = i_0, Y_1 = i_1, \dots, Y_r = i_r) = \mu_{i_0} \prod_{k=0}^{r-1} Q(k)_{i_k i_{k+1}}.
$$
From this joint distribution identity, whenever $r \ge 0$ and $i_0, \dots, i_{r+1} \in I$ satisfy $\mathbb{P}(Y_0 = i_0, \dots, Y_r = i_r) > 0$, we obtain
$$
\begin{aligned}
\mathbb{P}(Y_{r+1} = i_{r+1} \mid Y_0 = i_0, \dots, Y_r = i_r)
&= \frac{\mathbb{P}(Y_0 = i_0, \dots, Y_r = i_r, Y_{r+1} = i_{r+1})}{\mathbb{P}(Y_0 = i_0, \dots, Y_r = i_r)} \\
&= \frac{\mu_{i_0} \left(\prod_{k=0}^{r-1} Q(k)_{i_k i_{k+1}}\right) Q(r)_{i_r i_{r+1}}}{\mu_{i_0} \prod_{k=0}^{r-1} Q(k)_{i_k i_{k+1}}} \\
&= Q(r)_{i_r i_{r+1}} = (P^{n_{r+1} - n_r})_{i_r i_{r+1}}.
\end{aligned}
$$
Moreover, summing $\mathbb{P}(Y_0 = i_0, \dots, Y_r = i_r, Y_{r+1} = i_{r+1}) = \mathbb{P}(Y_0 = i_0, \dots, Y_r = i_r)\, Q(r)_{i_r i_{r+1}}$ over all $i_0, \dots, i_{r-1} \in I$ yields $\mathbb{P}(Y_r = i_r, Y_{r+1} = i_{r+1}) = \mathbb{P}(Y_r = i_r)\, Q(r)_{i_r i_{r+1}}$, so whenever $\mathbb{P}(Y_r = i_r) > 0$,
$$
\mathbb{P}(Y_{r+1} = i_{r+1} \mid Y_r = i_r) = Q(r)_{i_r i_{r+1}}.
$$
Thus $(Y_r)_{r\ge 0}$ satisfies the Markov property with initial distribution $\mu = \lambda P^{n_0}$ and step-$r$ transition matrix $Q(r) = P^{n_{r+1} - n_r}$. Because $n_{r+1} - n_r$ may depend on $r$, $Q(r)$ may depend on $r$, so $Y$ is in general an inhomogeneous Markov chain.

*Remark (inhomogeneous $X$):* If $X$ is itself an inhomogeneous Markov chain with step-$t$ transition matrices $P(t)$, the exact same block factorisation gives $\mu = \lambda P(0)\cdots P(n_0 - 1)$ and $Q(r) = P(n_r) P(n_r + 1) \cdots P(n_{r+1} - 1)$ (empty products being the identity matrix), so $Y$ is still a (possibly inhomogeneous) Markov chain.

---

### Part 2: Transition probabilities of $Y_r = X_{2r}$ when $X$ is a simple random walk

Let $X = (X_n)_{n\ge 0}$ be a simple random walk on $\mathbb{Z}$ with right-step probability $p \in [0, 1]$ and left-step probability $q = 1 - p$, i.e.
$$
p_{i, i+1} = p, \qquad p_{i, i-1} = q = 1 - p, \qquad p_{i, j} = 0 \quad (j \notin \{i-1, i+1\})
$$
for all $i \in \mathbb{Z}$.

When $n_r = 2r$ (for $r \ge 0$), the step increments are constant:
$$
n_{r+1} - n_r = 2(r+1) - 2r = 2 \qquad \text{for all } r \ge 0.
$$
Hence $Y = (X_{2r})_{r\ge 0}$ is a **time-homogeneous** Markov chain on $\mathbb{Z}$ (with initial distribution $\mu = \lambda P^0 = \lambda$) and one-step transition matrix $Q = P^2$.

Computing $q_{ij} = (P^2)_{ij} = \sum_{k\in\mathbb{Z}} p_{ik} p_{kj}$ for any $i, j \in \mathbb{Z}$:
- Only $k = i+1$ and $k = i-1$ satisfy $p_{ik} > 0$, so $q_{ij} = p_{i, i+1} p_{i+1, j} + p_{i, i-1} p_{i-1, j}$.
- For $j = i + 2$: only $k = i + 1$ contributes, giving $q_{i, i+2} = p \cdot p = p^2$.
- For $j = i$: both $k = i + 1$ and $k = i - 1$ contribute, giving $q_{i, i} = p \cdot q + q \cdot p = 2pq = 2p(1-p)$.
- For $j = i - 2$: only $k = i - 1$ contributes, giving $q_{i, i-2} = q \cdot q = q^2 = (1-p)^2$.
- For all other $j \in \mathbb{Z}$ (i.e. $j \notin \{i-2, i, i+2\}$): $q_{i, j} = 0$.

As a sanity check, $p^2 + 2p(1-p) + (1-p)^2 = (p + (1-p))^2 = 1$, and the boundary cases $p = 0$ (deterministic step $-2$) and $p = 1$ (deterministic step $+2$) are correctly captured.

In particular, if $X$ is a **simple symmetric random walk** ($p = q = 1/2$), the transition probabilities of $Y$ are:
$$
q_{i, i+2} = \frac{1}{4}, \qquad q_{i, i} = \frac{1}{2}, \qquad q_{i, i-2} = \frac{1}{4}, \qquad q_{i, j} = 0 \text{ otherwise.}
$$
(Note that $Y$ preserves parity, so $2\mathbb{Z}$ and $2\mathbb{Z}+1$ do not communicate; if $X_0 = 0$, $Y$ stays in $2\mathbb{Z}$.)

## Final Answers

1. $Y_r = X_{n_r}$ is a (possibly inhomogeneous) Markov chain with initial distribution $\mu = \lambda P^{n_0}$ and step-$r$ transition matrix $Q(r) = P^{n_{r+1} - n_r}$ (for homogeneous $X$ with transition matrix $P$; more generally $Q(r) = P(n_r)\cdots P(n_{r+1}-1)$ if $X$ is inhomogeneous).
2. For $n_r = 2r$ and a simple random walk $X$ on $\mathbb{Z}$ with step probabilities $p_{i, i+1} = p$, $p_{i, i-1} = 1 - p$ ($p \in [0, 1]$), $Y_r = X_{2r}$ is a **homogeneous** Markov chain with transition probabilities:
   $$
   q_{i, i+2} = p^2, \qquad q_{i, i} = 2p(1-p), \qquad q_{i, i-2} = (1-p)^2, \qquad q_{i, j} = 0 \text{ otherwise.}
   $$
   In the symmetric case ($p = 1/2$): $q_{i, i+2} = 1/4$, $q_{i, i} = 1/2$, $q_{i, i-2} = 1/4$, and $0$ otherwise.

## Computational Verification

- Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-1_q03_check.py`
- Command: `python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-1_q03_check.py"` (exit code 0)
- Verbatim stdout:
```text
PASS: Subsampled homogeneous chain Y_r = X_{n_r} verified for sequences [1,2,5,7], [0,1,3,6], [1,3,3,5].
PASS: Subsampled inhomogeneous chain verified with Q(r) = P(n_r)...P(n_{r+1}-1).
PASS: Simple random walk Y_r = X_{2r} transition probabilities verified for p in {0, 1/5, 1/3, 1/2, 3/4, 1}.
  General p: q_{i,i+2} = p^2, q_{i,i} = 2p(1-p), q_{i,i-2} = (1-p)^2.
  Symmetric p = 1/2: q_{i,i+2} = 1/4, q_{i,i} = 1/2, q_{i,i-2} = 1/4.
```

## Confidence & Edge-Case Notes

- **Confidence:** High.
- **Subtleties / Phase B comparison watchpoints:**
  1. **Convention for "simple random walk":** Check whether the existing solution assumes a general asymmetric simple random walk ($p, 1-p$) on $\mathbb{Z}$ or only the symmetric case ($p = 1/2$). Both are stated cleanly above; note Question 4 on the same sheet explicitly says "simple (possibly asymmetric) random walk on $\mathbb{Z}$", suggesting general $(p, 1-p)$ is the natural level of generality (with $p = 1/2$ as the symmetric special case).
  2. **$n_0 = 0$ vs "positive integers":** The sheet's statement says "positive integers" yet sets $n_r = 2r$ (giving $n_0 = 0$). Make sure the proof works cleanly for $n_0 = 0$ as well as $n_0 \ge 1$.
  3. **Inhomogeneity:** Ensure the proof explicitly identifies the step-$r$ transition matrix $P^{n_{r+1} - n_r}$ and explains why $Y$ can be inhomogeneous (varying $n_{r+1} - n_r$), and does not divide by zero when conditioning on $\mathbb{P}(Y_0=i_0, \dots, Y_r=i_r)$.
