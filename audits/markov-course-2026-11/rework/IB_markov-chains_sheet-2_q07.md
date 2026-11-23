---
id: IB/markov-chains/sheet-2/q07
slug: IB_markov-chains_sheet-2_q07
source_archive: private/raw/sheets/IB/markov-chains/2025-26/sheet-2.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

# Blind Rework: `IB/markov-chains/sheet-2/q07`

## Transcription Check

**Status:** `ok` (exact match).

Compared `question.tex` word-for-word and symbol-for-symbol against the rendered screenshot of page 1 of `private/raw/sheets/IB/markov-chains/2025-26/sheet-2.pdf` (Exercise 7, followed by Exercise 8 on the same page, so no truncation occurs):

- Printed PDF: *"7. The rooted binary tree is an infinite graph $T$ with one distinguished vertex $R$ from which comes a single edge; at every other vertex there are three edges and there are no closed loops. The random walk on $T$ jumps from a vertex along each available edge with equal probability. Show that the random walk is transient."*
- `question.tex`: exact match in every word, symbol, and punctuation mark.

---

## Full Solution

Let $V$ denote the countable vertex set of the rooted binary tree $T$, and let $(X_n)_{n \ge 0}$ be the simple random walk on $T$.

### 1. Structure of $T$ and irreducibility of $(X_n)$

Because $T$ is a tree (a connected graph with no closed loops), for every vertex $v \in V$ there is a unique simple path in $T$ from the root $R$ to $v$. Define the **depth** (graph distance from $R$) of $v \in V$ by
\[
d(v) := \text{length (number of edges) of the unique simple path in } T \text{ from } R \text{ to } v \in \mathbb{Z}_{\ge 0}.
\]
Then:
- $d(v) = 0$ if and only if $v = R$.
- Since $R$ has degree $1$, $R$ has a unique neighbour $v_1 \in V$, which satisfies $d(v_1) = 1$, and $p_{R, v_1} = 1$.
- For every vertex $v \in V \setminus \{R\}$, we have $d(v) \ge 1$ and $\deg(v) = 3$. One incident edge is the final edge of the unique simple path from $R$ to $v$, connecting $v$ to its unique parent $u$ with $d(u) = d(v) - 1$. The other two incident edges cannot lead to any vertex $w$ with $d(w) \le d(v)$ (nor to the same vertex), because otherwise $T$ would contain a closed loop. Hence those two edges lead to two distinct neighbours (children) of $v$ at depth $d(v) + 1$. Since the walk jumps along each incident edge with equal probability, $p_{vw} = \frac{1}{3}$ for each of the three neighbours $w$ of $v$.

For any two vertices $x, y \in V$, connectedness of $T$ yields a finite path $x = w_0, w_1, \dots, w_m = y$ of adjacent vertices in $T$. Along this path every transition probability is at least $\frac{1}{3}$, so
\[
p_{xy}^{(m)} \ge p_{w_0 w_1} p_{w_1 w_2} \cdots p_{w_{m-1} w_m} \ge \left(\frac{1}{3}\right)^m > 0.
\]
Hence $x \to y$ for all $x, y \in V$, so the Markov chain $(X_n)_{n \ge 0}$ is **irreducible**.

By the class property of recurrence and transience (for an irreducible Markov chain, either all states are recurrent or all states are transient), to show that $(X_n)_{n \ge 0}$ is transient it suffices to show that the single state $R$ is transient.

### 2. Reduction to the distance Markov chain $Y_n = d(X_n)$

Define $Y_n := d(X_n) \in \mathbb{Z}_{\ge 0}$ for $n \ge 0$. We first verify rigorously that $(Y_n)_{n \ge 0}$ is a Markov chain on $\mathbb{Z}_{\ge 0}$ with transition matrix $Q = (q_{ij})_{i,j \ge 0}$ given by
\[
q_{0,1} = 1, \qquad q_{i, i+1} = \frac{2}{3}, \qquad q_{i, i-1} = \frac{1}{3} \quad (i \ge 1),
\]
and $q_{ij} = 0$ otherwise.

Let $n \ge 0$ and let $y_0, y_1, \dots, y_n, y_{n+1} \in \mathbb{Z}_{\ge 0}$ be such that $\mathbb{P}(Y_0 = y_0, \dots, Y_n = y_n) > 0$. Write $A := \{Y_0 = y_0, \dots, Y_{n-1} = y_{n-1}\}$. Note that $A$ is determined by $(X_0, \dots, X_{n-1})$. For any vertex $x \in V$ with $d(x) = y_n$, by the Markov property of $(X_n)$,
\[
\mathbb{P}(Y_{n+1} = y_{n+1} \mid X_n = x, A) = \mathbb{P}(d(X_{n+1}) = y_{n+1} \mid X_n = x).
\]
From the tree structure above:
- If $y_n = 0$, then $x = R$ and $d(X_{n+1}) = 1$ with probability $1$, so $\mathbb{P}(d(X_{n+1}) = y_{n+1} \mid X_n = x) = q_{0, y_{n+1}}$.
- If $y_n \ge 1$, then $x \ne R$ has $1$ neighbour at depth $y_n - 1$ and $2$ neighbours at depth $y_n + 1$, each chosen with probability $\frac{1}{3}$, so $\mathbb{P}(d(X_{n+1}) = y_{n+1} \mid X_n = x) = q_{y_n, y_{n+1}}$.

Crucially, this probability depends on $x$ only through $d(x) = y_n$. Summing over the disjoint events $\{X_n = x\}$ with $d(x) = y_n$, we obtain
\[
\begin{aligned}
\mathbb{P}(Y_{n+1} = y_{n+1}, Y_n = y_n, A)
&= \sum_{x \in V:\, d(x) = y_n} \mathbb{P}(Y_{n+1} = y_{n+1} \mid X_n = x, A)\,\mathbb{P}(X_n = x, A) \\
&= \sum_{x \in V:\, d(x) = y_n} q_{y_n, y_{n+1}}\,\mathbb{P}(X_n = x, A) \\
&= q_{y_n, y_{n+1}}\,\mathbb{P}(Y_n = y_n, A).
\end{aligned}
\]
Dividing by $\mathbb{P}(Y_n = y_n, A) > 0$ gives
\[
\mathbb{P}(Y_{n+1} = y_{n+1} \mid Y_0 = y_0, \dots, Y_n = y_n) = q_{y_n, y_{n+1}}.
\]
Thus $(Y_n)_{n \ge 0}$ is indeed a Markov chain on $\mathbb{Z}_{\ge 0}$ with transition matrix $Q$.

Since $R$ is the unique vertex at depth $0$, we have $\{X_n = R\} = \{Y_n = 0\}$ for every $n \ge 0$. Consequently, starting from $X_0 = R$ (so $Y_0 = 0$), the first return time to $R$,
\[
T_R := \inf\{n \ge 1 : X_n = R\},
\]
coincides with the first return time of $(Y_n)$ to $0$,
\[
T_0^Y := \inf\{n \ge 1 : Y_n = 0\}.
\]
In particular, $\mathbb{P}_R(T_R < \infty) = \mathbb{P}_0(T_0^Y < \infty)$.

### 3. Computing the return probability and proving transience

For each $k \in \mathbb{Z}_{\ge 0}$, define the hitting time $H_0 := \inf\{n \ge 0 : Y_n = 0\}$ and hitting probability
\[
h_k := \mathbb{P}_k(H_0 < \infty).
\]
Starting from $Y_0 = 0$, the first jump goes deterministically to $Y_1 = 1$, so by the Markov property at time $1$,
\[
\mathbb{P}_0(T_0^Y < \infty) = \mathbb{P}_1(H_0 < \infty) = h_1.
\]
We determine $h_k$ via the **Hitting Probability Theorem** (and also verify the upper bound directly by induction so the proof is completely self-contained).

**Theorem (Hitting probabilities):** For a Markov chain $(Y_n)$ with transition matrix $Q$ and target set $A \subseteq I$, the vector of hitting probabilities $h_k = \mathbb{P}_k(H_A < \infty)$ is the minimal non-negative solution to
\[
h_k = 1 \quad (k \in A), \qquad h_k = \sum_{j \in I} q_{kj} h_j \quad (k \notin A).
\]

For $A = \{0\}$, the equations are
\[
h_0 = 1, \qquad h_k = \frac{2}{3} h_{k+1} + \frac{1}{3} h_{k-1} \quad \text{for all } k \ge 1.
\]
Rewriting the recurrence as
\[
2(h_{k+1} - h_k) = h_k - h_{k-1} \quad (k \ge 1),
\]
we obtain by induction
\[
h_{k+1} - h_k = \left(\frac{1}{2}\right)^k (h_1 - h_0) \quad (k \ge 0).
\]
Summing from $0$ to $m - 1$ yields, for every $m \ge 0$,
\[
h_m - h_0 = (h_1 - h_0) \sum_{k=0}^{m-1} 2^{-k} = 2(1 - 2^{-m})(h_1 - h_0).
\]
Since $h_0 = 1$, putting $B := 2(1 - h_1)$ gives the general solution
\[
h_m = (1 - B) + B \cdot 2^{-m} \qquad (m \ge 0).
\]
For $(h_m)_{m \ge 0}$ to be non-negative for all $m \ge 0$, taking $m \to \infty$ requires $1 - B \ge 0$, i.e. $B \le 1$. Since $1 - 2^{-m} \ge 0$ for every $m \ge 0$, writing $h_m = 2^{-m} + (1 - B)(1 - 2^{-m})$ shows that every non-negative solution satisfies $h_m \ge 2^{-m}$ for all $m \ge 0$, with equality attained by the non-negative solution $x_m^* := 2^{-m}$ (taking $B = 1$). Hence the minimal non-negative solution is $h_m = 2^{-m}$ for all $m \ge 0$.

*(Self-contained direct verification of $h_k \le 2^{-k}$ without quoting minimality:* let $h_k^{(n)} := \mathbb{P}_k(H_0 \le n)$ for $n \ge 0, k \ge 0$. For $n = 0$, $h_0^{(0)} = 1 = 2^0$ and $h_k^{(0)} = 0 \le 2^{-k}$ for $k \ge 1$. Assuming $h_k^{(n)} \le 2^{-k}$ for all $k \ge 0$, we have $h_0^{(n+1)} = 1 = 2^0$ and, for $k \ge 1$, conditioning on the first step gives
\[
h_k^{(n+1)} = \frac{2}{3} h_{k+1}^{(n)} + \frac{1}{3} h_{k-1}^{(n)} \le \frac{2}{3} \cdot 2^{-(k+1)} + \frac{1}{3} \cdot 2^{-(k-1)} = 2^{-k}\left(\frac{1}{3} + \frac{2}{3}\right) = 2^{-k}.
\]
By induction $h_k^{(n)} \le 2^{-k}$ for all $n, k \ge 0$. Since $\{H_0 \le n\} \uparrow \{H_0 < \infty\}$ as $n \to \infty$, continuity of probability from below gives $h_k = \lim_{n \to \infty} h_k^{(n)} \le 2^{-k}$. Together with $h_k \ge 2^{-k}$, we obtain $h_k = 2^{-k}$.)

In particular,
\[
\mathbb{P}_R(T_R < \infty) = \mathbb{P}_0(T_0^Y < \infty) = h_1 = \frac{1}{2} < 1.
\]
Therefore $R$ is a transient state of $(X_n)$. Because $(X_n)$ is irreducible and transience is a class property, **every** state of $(X_n)$ is transient; that is, the random walk on $T$ is transient. $\square$

---

## Final Answers

- The random walk $(X_n)_{n \ge 0}$ on the rooted binary tree $T$ is **transient** (and irreducible).
- Explicitly, the return probability to the root $R$ is $\mathbb{P}_R(T_R < \infty) = \frac{1}{2} < 1$ (and more generally $\mathbb{P}_v(H_R < \infty) = 2^{-d(v)}$, and $\mathbb{E}_R[\text{total visits to } R] = \sum_{n=0}^\infty p_{RR}^{(n)} = 2 < \infty$).

---

## Computational Verification

Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_q07_check.py`

Command: `python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_q07_check.py"` (exit code 0)

Verbatim stdout:
```text
[PASS] Truncated tree T_6 (64 vertices): lumpability verified and exact hitting prob h_v = (2^-d - 2^-6)/(1 - 2^-6) for all vertices.
       h_1^(L=6) = 31/63 = 0.49206349 (limit as L->inf: 1/2)
[PASS] Exact Catalan formula f_00^(2m) = (1/3)*C_(m-1)*(2/9)^(m-1) verified for m=1..15.
       Partial sum P_0(T_0 <= 30) = 0.4970538634 (< 0.5, exact limit 1/2)
       Partial sum sum_{n=0}^30 q^(n)(0,0) = 1.9697218185 (< 2.0, exact limit 2)
[PASS] Direct induction bound h_k^(n) <= 2^(-k) verified for 0 <= n,k <= 30; h_1^(30) = 0.4970538634.
ALL CHECKS PASSED.
```

---

## Confidence & Edge-Case Notes

- **Confidence:** `high`.
- **Subtleties / pitfalls to watch for in Phase B:**
  1. **Lumpability justification:** Claiming $Y_n = d(X_n)$ is a Markov chain without justifying why the transition probabilities of $d(X_{n+1})$ depend only on $d(X_n)$ (and not on the specific vertex $X_n$ or past trajectory) is a common gap.
  2. **Irreducibility and class property:** Showing $R$ is transient only proves transience of the entire walk once one notes that $(X_n)$ is irreducible (since $T$ is connected) and transience is a class property.
  3. **Minimal non-negative solution:** When solving $h_k = \frac{2}{3}h_{k+1} + \frac{1}{3}h_{k-1}$ to get $h_k = (1-B) + B \cdot 2^{-k}$, one must explicitly invoke the **minimal** non-negative solution (which forces $B=1$, NOT $B=0$, since $1 - 2^{-k} > 0$), rather than guessing $h_k \to 0$ as $k \to \infty$ without proof! (Our direct induction bound $h_k^{(n)} \le 2^{-k}$ also verifies this from first principles.)
  4. **Bijection at depth $0$:** Transience of $0$ for $Y_n$ implies transience of $R$ for $X_n$ specifically because $R$ is the **unique** vertex with $d(v) = 0$ (so $\{X_n = R\} = \{Y_n = 0\}$). (In contrast, for a vertex $v$ at depth $k \ge 1$, $\{Y_n = k\}$ does not imply $\{X_n = v\}$, though of course transience of $R$ + irreducibility yields transience of every vertex $v$.)
