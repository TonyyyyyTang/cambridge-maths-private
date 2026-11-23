---
id: IB/markov-chains/sheet-2/q14
slug: IB_markov-chains_sheet-2_q14
source_archive: private/raw/sheets/IB/markov-chains/2025-26/sheet-2.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

# Blind Rework: `IB/markov-chains/sheet-2/q14`

## Transcription Check

**Status:** `ok` (exact match).

Compared `question.tex` word-for-word and symbol-for-symbol against the rendered screenshot of page 2 of `private/raw/sheets/IB/markov-chains/2025-26/sheet-2.pdf` (Exercise 14, followed by Exercises 15 and 16 on the same page):

- Intro: *"In each of the following cases determine whether the stochastic matrix $P$ corresponds to a chain which is reversible in equilibrium:"* — exact match.
- (a): $\begin{pmatrix} p & 1-p \\ q & 1-q \end{pmatrix}$; — exact match (including trailing semicolon).
- (b): $\begin{pmatrix} 0 & p & 1-p \\ 1-p & 0 & p \\ p & 1-p & 0 \end{pmatrix}$; — exact match (all 9 entries and trailing semicolon).
- (c): $I = \{0, 1, 2, \dots\}$ and $p_{01} = 1$, $p_{i,i+1} = p$, $p_{i,i-1} = 1-p$ for $i \ge 1$. — exact match (including trailing period).
- (d): $p_{ij} = p_{ji}$ for all $i, j \in I$. — exact match (including trailing period).

---

## Full Solution

### Preliminaries: Definition and Criterion for Reversibility in Equilibrium

Let $P = (p_{ij})_{i,j \in I}$ be a stochastic matrix on a countable state space $I$.
For a Markov chain with transition matrix $P$ to be *in equilibrium*, there must exist an **invariant (probability) distribution** $\pi = (\pi_i)_{i \in I}$, meaning
\[
\pi_i \ge 0 \quad (\forall i \in I), \qquad \sum_{i \in I} \pi_i = 1, \qquad \text{and} \qquad \pi P = \pi.
\]
Given an invariant distribution $\pi$, the equilibrium chain $(X_n)_{n \ge 0} \sim \mathrm{Markov}(\pi, P)$ is **time-reversible** if for every integer $N \ge 1$ and all $i_0, i_1, \dots, i_N \in I$,
\[
\mathbb{P}_\pi(X_0 = i_0, X_1 = i_1, \dots, X_N = i_N) = \mathbb{P}_\pi(X_0 = i_N, X_1 = i_{N-1}, \dots, X_N = i_0).
\]

**Lemma (Detailed balance characterisation):** Let $\pi$ be a probability distribution on $I$.
1. If $\pi$ and $P$ satisfy the **detailed balance equations**
   \[
   \pi_i p_{ij} = \pi_j p_{ji} \qquad \text{for all } i, j \in I,
   \]
   then $\pi$ is automatically an invariant distribution for $P$, and $\mathrm{Markov}(\pi, P)$ is time-reversible.
2. Conversely, if $\pi$ is an invariant distribution for $P$ and $\mathrm{Markov}(\pi, P)$ is time-reversible, then $\pi$ and $P$ satisfy detailed balance.

*Proof.*
(1) Summing detailed balance over $i \in I$ gives $(\pi P)_j = \sum_{i \in I} \pi_i p_{ij} = \sum_{i \in I} \pi_j p_{ji} = \pi_j \sum_{i \in I} p_{ji} = \pi_j$ for all $j \in I$, so $\pi P = \pi$. For any $N \ge 1$ and $i_0, \dots, i_N \in I$, repeatedly substituting $\pi_{i_k} p_{i_k i_{k+1}} = p_{i_{k+1} i_k} \pi_{i_{k+1}}$ gives
\[
\begin{aligned}
\mathbb{P}_\pi(X_0 = i_0, \dots, X_N = i_N)
&= \pi_{i_0} p_{i_0 i_1} p_{i_1 i_2} \cdots p_{i_{N-1} i_N} \\
&= p_{i_1 i_0} \pi_{i_1} p_{i_1 i_2} \cdots p_{i_{N-1} i_N} \\
&= \cdots = p_{i_1 i_0} p_{i_2 i_1} \cdots p_{i_N i_{N-1}} \pi_{i_N} \\
&= \mathbb{P}_\pi(X_0 = i_N, \dots, X_N = i_0).
\end{aligned}
\]
(2) Taking $N = 1$ in the definition of time-reversibility directly yields $\pi_i p_{ij} = \mathbb{P}_\pi(X_0 = i, X_1 = j) = \mathbb{P}_\pi(X_0 = j, X_1 = i) = \pi_j p_{ji}$ for all $i, j \in I$. $\square$

*(Remark: Note that this equivalence holds for arbitrary $P$ and invariant distribution $\pi$, without requiring irreducibility, aperiodicity, or $\pi_i > 0$. However, when $P$ has no invariant distribution at all, no equilibrium exists and the chain cannot be reversible in equilibrium; and when $P$ is reducible with multiple invariant distributions, one should check whether detailed balance holds for all invariant distributions.)*

We now analyse parts (a)–(d) in full, including all boundary and degenerate parameter values.

---

### Part (a)

\[
P = \begin{pmatrix} p & 1-p \\ q & 1-q \end{pmatrix}, \qquad I = \{1, 2\}.
\]
Since $P$ is a stochastic matrix, $p, q \in [0, 1]$.

Let $\pi = (\pi_1, \pi_2)$ be a probability distribution ($\pi_1, \pi_2 \ge 0$, $\pi_1 + \pi_2 = 1$). The invariance equations $\pi P = \pi$ are
\[
\pi_1 p + \pi_2 q = \pi_1 \iff \pi_1(1 - p) = \pi_2 q,
\]
and the second equation $\pi_1(1 - p) + \pi_2(1 - q) = \pi_2$ is equivalent to the same relation $\pi_1(1 - p) = \pi_2 q$.

Meanwhile, the detailed balance equations $\pi_i p_{ij} = \pi_j p_{ji}$ are trivial for $i = j$, and for $\{i, j\} = \{1, 2\}$ reduce to the single equation
\[
\pi_1 p_{12} = \pi_2 p_{21} \iff \pi_1(1 - p) = \pi_2 q,
\]
which is **identical** to the invariance equation $\pi P = \pi$.

Checking existence of $\pi$ across all $(p, q) \in [0, 1]^2$:
- **Case 1: $(p, q) \ne (1, 0)$** (so $1 - p + q > 0$). There is a **unique** invariant distribution
  \[
  \pi = \left(\frac{q}{1 - p + q},\; \frac{1 - p}{1 - p + q}\right),
  \]
  and it satisfies detailed balance. (Specifically: for $0 \le p < 1$ and $0 < q \le 1$, $P$ is irreducible — aperiodic except at $(p, q) = (0, 1)$, where it has period $2$ and is still reversible in equilibrium with $\pi = (\frac{1}{2}, \frac{1}{2})$; for $p = 1, q > 0$, state $1$ is absorbing and $\pi = (1, 0)$; for $p < 1, q = 0$, state $2$ is absorbing and $\pi = (0, 1)$.)
- **Case 2: $(p, q) = (1, 0)$.** Then $P = I_2$ (reducible, both states absorbing). **Every** distribution $\pi = (\alpha, 1 - \alpha)$ with $\alpha \in [0, 1]$ is invariant, and every such $\pi$ satisfies detailed balance ($\pi_1 \cdot 0 = 0 = \pi_2 \cdot 0$), since $X_n = X_0$ a.s. for all $n$.

**Conclusion for (a):** **Yes — reversible in equilibrium for all $p, q \in [0, 1]$** (every invariant distribution satisfies detailed balance, and at least one always exists).

---

### Part (b)

\[
P = \begin{pmatrix} 0 & p & 1-p \\ 1-p & 0 & p \\ p & 1-p & 0 \end{pmatrix}, \qquad I = \{1, 2, 3\}.
\]
Since $P$ is a stochastic matrix, $p \in [0, 1]$.

First, $P$ is **irreducible for every $p \in [0, 1]$**:
- For $0 < p < 1$, all off-diagonal entries are strictly positive ($p_{ij} > 0$ for all $i \ne j$).
- For $p = 1$, $P$ is the deterministic $3$-cycle $1 \to 2 \to 3 \to 1$, which is irreducible (of period $3$).
- For $p = 0$, $P$ is the deterministic $3$-cycle $1 \to 3 \to 2 \to 1$, which is irreducible (of period $3$).

Since $I$ is finite and $P$ is irreducible, $P$ has a **unique** invariant distribution $\pi$ for every $p \in [0, 1]$. Moreover, every column of $P$ sums to $0 + p + (1 - p) = 1$, so $P$ is doubly stochastic and hence the unique invariant distribution is
\[
\pi = \left(\frac{1}{3},\; \frac{1}{3},\; \frac{1}{3}\right)
\]
for all $p \in [0, 1]$.

Now check detailed balance $\pi_i p_{ij} = \pi_j p_{ji}$. Since $\pi_1 = \pi_2 = \pi_3 = \frac{1}{3} > 0$, detailed balance holds if and only if $p_{ij} = p_{ji}$ for all $i, j \in \{1, 2, 3\}$. Looking at $(i, j) = (1, 2)$,
\[
p_{12} = p, \qquad p_{21} = 1 - p,
\]
so detailed balance requires $p = 1 - p$, i.e. $p = \frac{1}{2}$. Conversely, when $p = \frac{1}{2}$, we have $p_{12} = p_{21} = p_{23} = p_{32} = p_{31} = p_{13} = \frac{1}{2}$, so $P$ is symmetric and $\pi_i p_{ij} = \pi_j p_{ji}$ holds for all $i, j \in \{1, 2, 3\}$.

**Conclusion for (b):** **Reversible in equilibrium if and only if $p = \frac{1}{2}$.**

---

### Part (c)

\[
I = \{0, 1, 2, \dots\}, \qquad p_{01} = 1, \qquad p_{i, i+1} = p, \qquad p_{i, i-1} = 1 - p \quad (i \ge 1),
\]
with all other entries $p_{ij} = 0$, and $p \in [0, 1]$ so that $P$ is a stochastic matrix.

The invariance equations $(\pi P)_j = \pi_j$ are:
\[
\begin{aligned}
j = 0: &\qquad (1 - p)\pi_1 = \pi_0, \\
j = 1: &\qquad \pi_0 + (1 - p)\pi_2 = \pi_1, \\
j \ge 2: &\qquad p\pi_{j-1} + (1 - p)\pi_{j+1} = \pi_j.
\end{aligned}
\]
We analyse all $p \in [0, 1]$:

1. **Case $p = 1$:**
   Here $p_{i, i+1} = 1$ for all $i \ge 0$. The invariance equations give $\pi_0 = 0$ and $\pi_j = \pi_{j-1}$ for all $j \ge 1$, hence $\pi_j = 0$ for all $j \ge 0$. Thus **no invariant distribution exists** (all states are transient), so the chain is **not** reversible in equilibrium.

2. **Case $0 < p < 1$:**
   Since $p > 0$ and $1 - p > 0$, $P$ is **irreducible** (and periodic of period $2$, as every transition changes the parity of the state). Solving $(\pi P)_j = \pi_j$ recursively:
   - From $j = 0$: $(1 - p)\pi_1 = \pi_0$, i.e. $\pi_1 = \frac{1}{1-p}\pi_0$.
   - From $j = 1$: $(1 - p)\pi_2 = \pi_1 - \pi_0 = \pi_1 - (1 - p)\pi_1 = p\pi_1$.
   - Inductively, for $j \ge 2$, assuming $(1 - p)\pi_j = p\pi_{j-1}$, the $j$-th invariance equation gives
     \[
     (1 - p)\pi_{j+1} = \pi_j - p\pi_{j-1} = \pi_j - (1 - p)\pi_j = p\pi_j.
     \]
   Hence every invariant measure $\pi$ satisfies the detailed balance equations $\pi_0 \cdot 1 = (1-p)\pi_1$ and $p\pi_i = (1-p)\pi_{i+1}$ ($i \ge 1$), and has the exact form
   \[
   \pi_i = \frac{1}{1 - p}\left(\frac{p}{1 - p}\right)^{i-1} \pi_0 \qquad (i \ge 1).
   \]
   Summing over $i \ge 0$, with $r := \frac{p}{1-p} > 0$:
   \[
   \sum_{i=0}^\infty \pi_i = \pi_0 \left(1 + \frac{1}{1 - p} \sum_{k=0}^\infty r^k\right).
   \]
   - **Subcase $0 < p < \frac{1}{2}$:** Here $0 < r = \frac{p}{1-p} < 1$, so the geometric series converges to $\frac{1}{1 - r} = \frac{1-p}{1-2p}$, giving
     \[
     \sum_{i=0}^\infty \pi_i = \pi_0 \left(1 + \frac{1}{1 - 2p}\right) = \pi_0 \frac{2(1 - p)}{1 - 2p}.
     \]
     Setting this equal to $1$ yields the **unique invariant distribution**
     \[
     \pi_0 = \frac{1 - 2p}{2(1 - p)}, \qquad \pi_i = \frac{1 - 2p}{2(1 - p)^2}\left(\frac{p}{1 - p}\right)^{i-1} \quad (i \ge 1).
     \]
     Since $\pi$ satisfies the detailed balance equations $\pi_0 p_{01} = \pi_1 p_{10}$ and $\pi_i p_{i,i+1} = \pi_{i+1} p_{i+1,i}$ ($i \ge 1$), and $p_{ij} = p_{ji} = 0$ whenever $|i - j| \ne 1$, the chain **is reversible in equilibrium**.
   - **Subcase $\frac{1}{2} \le p < 1$:** Here $r = \frac{p}{1-p} \ge 1$. If $\pi_0 > 0$, then $\pi_i \ge \frac{\pi_0}{1-p} > 0$ for all $i \ge 1$, so $\sum_{i=0}^\infty \pi_i = \infty$; if $\pi_0 = 0$, then $\pi_i = 0$ for all $i \ge 0$. Hence **no invariant distribution exists** (the chain is null recurrent for $p = \frac{1}{2}$ and transient for $\frac{1}{2} < p < 1$; although there is a non-zero invariant *measure* in detailed balance, it has infinite total mass and cannot be normalised to a probability distribution). Consequently, the chain has no equilibrium and is **not** reversible in equilibrium.

3. **Boundary case $p = 0$:**
   Here $p_{01} = 1$ and $p_{i, i-1} = 1$ for all $i \ge 1$, while $p_{i, i+1} = 0$ for all $i \ge 1$. The chain is reducible: $\{0, 1\}$ is a closed communicating class (period $2$) and $\{2, 3, \dots\}$ are transient states. The invariance equations become $\pi_1 = \pi_0$, $\pi_0 + \pi_2 = \pi_1$ (hence $\pi_2 = 0$), and $\pi_{j+1} = \pi_j$ for $j \ge 2$ (hence $\pi_j = 0$ for all $j \ge 2$). Thus there is a **unique invariant distribution**
   \[
   \pi = \left(\frac{1}{2},\; \frac{1}{2},\; 0,\; 0,\; \dots\right).
   \]
   Checking detailed balance: $\pi_0 p_{01} = \frac{1}{2} \cdot 1 = \pi_1 p_{10}$, and for all $i \ge 1$, $\pi_i p_{i,i+1} = \pi_i \cdot 0 = 0 = 0 \cdot 1 = \pi_{i+1} p_{i+1, i}$. Hence detailed balance holds and the chain **is reversible in equilibrium** at $p = 0$ as well! *(Note: if one implicitly assumes the chain is irreducible, i.e. $p \in (0, 1)$, then the condition is $0 < p < \frac{1}{2}$.)*

**Conclusion for (c):** **Reversible in equilibrium if and only if $0 \le p < \frac{1}{2}$** (or $0 < p < \frac{1}{2}$ if restricted to irreducible chains). For $\frac{1}{2} \le p \le 1$, no invariant distribution exists.

---

### Part (d)

\[
p_{ij} = p_{ji} \qquad \text{for all } i, j \in I.
\]
Here $P$ is a symmetric stochastic matrix on a countable state space $I$. We determine the exact necessary and sufficient condition for $P$ to have an equilibrium and be reversible in equilibrium — treating both finite and infinite $I$, and both irreducible and reducible $P$.

**Step 1: Every communicating class is closed.**
Because $p_{ij} = p_{ji}$ for all $i, j \in I$, for every $n \ge 1$ and all $i, j \in I$,
\[
p_{ij}^{(n)} = \sum_{k_1, \dots, k_{n-1} \in I} p_{i k_1} p_{k_1 k_2} \cdots p_{k_{n-1} j} = \sum_{k_1, \dots, k_{n-1} \in I} p_{j k_{n-1}} \cdots p_{k_2 k_1} p_{k_1 i} = p_{ji}^{(n)}.
\]
Thus $i \to j \iff j \to i$. Consequently, every communicating class $C \subseteq I$ is **closed**, and there are no transitions between distinct communicating classes ($p_{ij} = 0$ whenever $i, j$ lie in different classes).

**Step 2: Structure of any invariant distribution $\pi$ (maximum principle).**
Suppose $\pi$ is any invariant distribution for $P$ ($\pi_i \ge 0$, $\sum_{i \in I} \pi_i = 1$, $\pi P = \pi$), and let $C \subseteq I$ be any communicating class. Set $w_C := \sum_{i \in C} \pi_i \in [0, 1]$.

Suppose $w_C > 0$. Pick $i^* \in C$ with $\pi_{i^*} > 0$. Because $\sum_{i \in C} \pi_i = w_C \le 1$, the set $\{i \in C : \pi_i \ge \pi_{i^*}\}$ is finite and non-empty, so $\pi$ attains its maximum on $C$,
\[
M_C := \max_{i \in C} \pi_i > 0,
\]
at some state $j_0 \in C$. Since no state outside $C$ can transition into $C$ and $p_{i j_0} = p_{j_0 i}$ with $\sum_{i \in C} p_{j_0 i} = 1$, the invariance equation $(\pi P)_{j_0} = \pi_{j_0}$ gives
\[
M_C = \pi_{j_0} = \sum_{i \in C} \pi_i p_{i j_0} = \sum_{i \in C} \pi_i p_{j_0 i},
\]
and therefore
\[
\sum_{i \in C} (M_C - \pi_i) p_{j_0 i} = 0.
\]
Since every term is non-negative, whenever $p_{j_0 i} > 0$ we must have $\pi_i = M_C$. Repeating inductively along paths, any $k \in C$ reachable from $j_0$ in $n$ positive-probability steps satisfies $\pi_k = M_C$. Because $C$ is a communicating class, $\pi_k = M_C$ for **all** $k \in C$.

Hence $w_C = \sum_{k \in C} \pi_k = M_C |C|$.
- If $|C| = \infty$, this contradicts $w_C \le 1$ (as $M_C > 0$). Therefore we must have $w_C = 0$, i.e. $\pi_k = 0$ for all $k \in C$.
- If $|C| < \infty$, then $\pi_k = \frac{w_C}{|C|}$ is **constant** on $C$.

**Step 3: Necessary and sufficient condition for existence and reversibility.**
1. **Existence of an equilibrium:** By Step 2, an invariant distribution $\pi$ exists **if and only if $P$ has at least one finite communicating class**!
   - If $P$ has at least one finite communicating class $C_0$, then the uniform distribution on $C_0$ ($\pi_i = \frac{1}{|C_0|}$ for $i \in C_0$, $\pi_i = 0$ for $i \notin C_0$) is an invariant distribution (and the general invariant distribution is any convex combination $\sum_{C \text{ finite}} w_C u^{(C)}$ of uniform distributions on finite communicating classes).
   - If **no** communicating class is finite, then Step 2 forces $\pi_i = 0$ for all $i \in I$, contradicting $\sum_{i \in I} \pi_i = 1$; thus **no invariant distribution exists**!
2. **Automatic detailed balance whenever an equilibrium exists:** Let $\pi$ be *any* invariant distribution of $P$. For any $i, j \in I$:
   - If $i, j$ lie in different communicating classes, then $p_{ij} = 0 = p_{ji}$, so $\pi_i p_{ij} = 0 = \pi_j p_{ji}$.
   - If $i, j$ lie in the same communicating class $C$, then by Step 2 we have $\pi_i = \pi_j$ (either both $\frac{w_C}{|C|}$ if $|C| < \infty$, or both $0$ if $|C| = \infty$). Together with $p_{ij} = p_{ji}$, this gives $\pi_i p_{ij} = \pi_j p_{ji}$.
   Hence **every** invariant distribution $\pi$ satisfies detailed balance, so the chain is reversible in equilibrium.

**Special cases to highlight:**
- **If $I$ is finite:** Every communicating class is finite, so an invariant distribution always exists (for irreducible $P$, uniquely the uniform distribution $\pi_i = \frac{1}{|I|}$; for reducible $P$, any mixture of uniform distributions on its closed classes) and every equilibrium is **reversible**.
- **If $I$ is infinite and $P$ is irreducible:** The single communicating class $I$ is infinite, so **no invariant distribution exists** (for example, the simple symmetric random walk on $\mathbb{Z}$ with $p_{i,i\pm 1} = \frac{1}{2}$ is symmetric and null recurrent, with only the unnormalisable counting measure $\mu_i \equiv 1$); hence it is **not** reversible in equilibrium.
- **If $I$ is infinite and $P$ is reducible:** It **is** reversible in equilibrium **if and only if** $P$ has at least one finite communicating class (for example, $P = I_\infty$ on $\mathbb{Z}_{\ge 0}$, or $p_{01} = p_{10} = 1$ combined with reflected simple symmetric random walk on $\{2, 3, \dots\}$); if all communicating classes are infinite, no equilibrium exists and it is **not** reversible in equilibrium.

**Conclusion for (d):** **Reversible in equilibrium if and only if $P$ has at least one finite communicating class** (in particular: **always** when $I$ is finite; and for irreducible $P$, **if and only if $I$ is finite**). Whenever an invariant distribution exists, every invariant distribution is constant on each communicating class and satisfies detailed balance. $\square$

---

## Final Answers

- **(a)** **Yes, always** (for all $p, q \in [0, 1]$). Unique $\pi = \left(\frac{q}{1-p+q}, \frac{1-p}{1-p+q}\right)$ when $(p, q) \ne (1, 0)$; every distribution $\pi = (\alpha, 1-\alpha)$ when $(p, q) = (1, 0)$.
- **(b)** **Yes if and only if $p = \frac{1}{2}$** (for all $p \in [0, 1]$, $P$ is irreducible with unique invariant distribution $\pi = (\frac{1}{3}, \frac{1}{3}, \frac{1}{3})$).
- **(c)** **Yes if and only if $0 \le p < \frac{1}{2}$** (unique invariant distribution $\pi_0 = \frac{1-2p}{2(1-p)}$, $\pi_i = \frac{1-2p}{2(1-p)^2}\left(\frac{p}{1-p}\right)^{i-1}$ for $i \ge 1$, reducing to $\pi = (\frac{1}{2}, \frac{1}{2}, 0, \dots)$ at the reducible boundary $p = 0$; for $\frac{1}{2} \le p \le 1$, no invariant distribution exists; if restricted to irreducible chains, $0 < p < \frac{1}{2}$).
- **(d)** **Yes if and only if $P$ has at least one finite communicating class** — in particular, **always if $I$ is finite**, and **for irreducible $P$, if and only if $I$ is finite** (when all communicating classes are infinite, e.g. simple symmetric random walk on $\mathbb{Z}$, no invariant distribution exists, though the counting measure $\mu_i \equiv 1$ is an unnormalisable invariant measure in detailed balance).

---

## Computational Verification

Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_q14_check.py`

Command: `python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_q14_check.py"` (exit code 0)

Verbatim stdout:
```text
[PASS] Part (a): all 25 parameter pairs (p,q) in grid (including all 4 corners) have invariant distribution(s), all satisfying detailed balance & path reversal (N<=4).
[PASS] Part (b): unique invariant dist (1/3,1/3,1/3) for all p in [0,1] (incl. p=0,1); detailed balance & path reversal hold iff p = 1/2.
[PASS] Part (c): p=0 (reducible, unique pi=(1/2,1/2,0,...)) and 0<p<1/2 (irreducible, exact pi) satisfy invariance, detailed balance & path reversal; 1/2<=p<=1 have no invariant distribution.
[PASS] Part (d): verified finite irreducible symmetric (unique uniform pi), finite reducible symmetric (all class-wise constant mixtures reversible), and infinite reducible symmetric with finite closed class (reversible in equilibrium).
ALL CHECKS PASSED.
```

---

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Subtleties / pitfalls to watch for in Phase B:**
  1. **Part (a) degenerate cases:** Check whether $(p, q) = (1, 0)$ (where $1 - p + q = 0$, so the formula $\left(\frac{q}{1-p+q}, \frac{1-p}{1-p+q}\right)$ divides by $0$!) and absorbing boundaries ($p = 1$ or $q = 0$) and periodic boundary $(p, q) = (0, 1)$ are handled accurately.
  2. **Part (b) uniqueness at $p = 0, 1$:** One must verify that $\pi = (\frac{1}{3}, \frac{1}{3}, \frac{1}{3})$ is the *unique* invariant distribution (including at $p = 0$ and $p = 1$, where $P$ is a deterministic 3-cycle and still irreducible); otherwise showing detailed balance fails for $(\frac{1}{3}, \frac{1}{3}, \frac{1}{3})$ alone wouldn't rule out another invariant distribution.
  3. **Part (c) boundary $p = 0$ and threshold $p = \frac{1}{2}$:**
     - Watch out for $p = 0$: it is reducible ($\{0, 1\}$ closed, $\{2, 3, \dots\}$ transient), yet has the unique invariant distribution $(\frac{1}{2}, \frac{1}{2}, 0, \dots)$ satisfying detailed balance! Many solutions forget $p = 0$ or write $0 < p < \frac{1}{2}$ without mentioning $p = 0$ (and $p = 1$).
     - Watch out for $\frac{1}{2} \le p < 1$: detailed balance has a non-zero solution $\mu$, which is only an infinite *measure*, NOT an invariant *distribution*, so there is no equilibrium.
     - Note also the chain has period $2$ for all $0 \le p < 1$, which does **not** prevent reversibility in equilibrium.
  4. **Part (d) finite vs infinite AND irreducible vs reducible:**
     - A classic error is claiming "always reversible since $\pi_i = \text{const}$ satisfies detailed balance" without noticing that on an infinite state space $I$ (e.g. irreducible symmetric chain like SSRW on $\mathbb{Z}$) no invariant *distribution* exists!
     - A subtler second error is claiming "reversible iff $I$ is finite" without noticing that a **reducible** symmetric chain on an infinite state space $I$ (such as $P = I_\infty$, or $p_{01} = p_{10} = 1$ plus reflected SSRW on $\{2, 3, \dots\}$) *does* have invariant distributions (on its finite communicating classes), all of which satisfy detailed balance! The exact necessary and sufficient condition is that $P$ possesses at least one finite communicating class (which reduces to $|I| < \infty$ when $P$ is irreducible).
     - Also watch out on finite reducible symmetric chains: one must prove that *every* invariant distribution (not just the uniform distribution on all of $I$) satisfies detailed balance, which our maximum-principle proof establishes cleanly.
