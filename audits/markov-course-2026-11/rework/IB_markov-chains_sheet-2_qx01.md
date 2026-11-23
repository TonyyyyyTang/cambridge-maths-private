---
id: IB/markov-chains/sheet-2/qx01
slug: IB_markov-chains_sheet-2_qx01
source_archive: private/raw/sheets/IB/markov-chains/2012-13/sheet-2.pdf
pages_checked: [1, 2]
transcription: ok
blindness: clean
---

## Transcription Check

**Status:** `ok` (exact match).

Visual inspection of Page 1 (Question 6) of `private/raw/sheets/IB/markov-chains/2012-13/sheet-2.pdf` (and confirmation that Question 7 follows immediately on Page 1, with nothing spilling onto Page 2) shows a word-for-word, symbol-for-symbol, punctuation-for-punctuation match with `question.tex`:
- Opening sentence: "A particle moves on the eight vertices of a cube in the following way: at each step the particle is equally likely to move to each of the three adjacent vertices, independently of its past motion. Let $i$ be the initial vertex occupied by the particle, $o$ the vertex opposite $i$. Calculate each of the following quantities:" — exact match.
- Subparts `(i)`, `(ii)`, `(iii)`: exact match, including semicolons after `(i)` and `(ii)` and a full stop after `(iii)`.
- No figures, hints, or extra conditions are omitted.

## Full Solution

### Setup and Properties of the Markov Chain

Identify the eight vertices of the cube with $V = \{0,1\}^3$, where two vertices $u, v \in V$ are adjacent if and only if their Hamming distance $d(u, v) = \sum_{r=1}^3 |u_r - v_r|$ equals $1$. Without loss of generality, let $i = (0,0,0)$ and let the opposite vertex be $o = (1,1,1)$ (so $d(i,o) = 3$).

Let $(X_n)_{n \ge 0}$ be the position of the particle at step $n$, starting from $X_0 = i$. Since each step is chosen uniformly from the three adjacent vertices independently of the past, $(X_n)_{n \ge 0}$ is a discrete-time Markov chain on the finite state space $V$ ($|V| = 8$) with transition matrix $P = (p_{uv})_{u,v \in V}$ given by
$$p_{uv} = \begin{cases} \frac{1}{3} & \text{if } d(u,v) = 1, \\ 0 & \text{otherwise.} \end{cases}$$

We verify the basic structural properties of $(X_n)$:
1. **Irreducibility:** For any $u, v \in V$, flipping the $d(u,v) \le 3$ differing coordinates one by one gives a path of length $d(u,v)$ along edges of the cube from $u$ to $v$, so $p_{uv}^{(d(u,v))} \ge (1/3)^3 > 0$. Thus every state communicates with every other state, so the chain is irreducible.
2. **Positive recurrence and uniqueness of the invariant distribution:** An irreducible Markov chain on a finite state space is positive recurrent and possesses a unique invariant distribution $\pi$.
3. **Explicit invariant distribution:** Since the adjacency relation on the cube is symmetric, $p_{uv} = p_{vu}$ for all $u, v \in V$ (so $P$ is symmetric, hence doubly stochastic). Taking the uniform distribution $\pi_v = \frac{1}{8}$ for all $v \in V$, we have $\sum_{v \in V} \pi_v = 1$ and, for every $v \in V$,
   $$(\pi P)_v = \sum_{u \in V} \frac{1}{8} p_{uv} = \frac{1}{8} \sum_{u \in V} p_{vu} = \frac{1}{8} = \pi_v$$
   (equivalently, $\pi_u p_{uv} = \pi_v p_{vu}$ for all $u, v \in V$, so the chain is reversible with respect to $\pi$). By uniqueness, $\pi_v = \frac{1}{8}$ ($v \in V$) is the unique invariant distribution.

*(Remark on periodicity: because the cube graph is bipartite — the parity of $u_1 + u_2 + u_3$ flips at every step — the chain has period $2$. However, the standard theorems for expected return times and expected visits between returns require only irreducibility and (positive) recurrence, not aperiodicity.)*

Define the first return time to $i$ by
$$T_i = \inf\{n \ge 1 : X_n = i\}$$
and the first hitting time of $o$ by
$$H_o = \inf\{n \ge 0 : X_n = o\}.$$
Since $X_0 = i \neq o$, we have $H_o = \inf\{n \ge 1 : X_n = o\} \ge 3$ almost surely.

---

### (i) Expected number of steps until the particle returns to $i$

Since $(X_n)$ is irreducible and positive recurrent with unique invariant distribution $\pi$, the standard theorem relating invariant distributions to expected return times gives
$$m_v = \mathbb{E}_v[T_v] = \frac{1}{\pi_v} \qquad \text{for all } v \in V.$$
Evaluating at $v = i$, where $\pi_i = \frac{1}{8}$, yields
$$\mathbb{E}_i[T_i] = \frac{1}{1/8} = 8.$$

*(Independent check via the distance chain in part (iii) below: starting from $i$, the first step moves deterministically to a vertex at distance $1$ from $i$. By symmetry of the cube, the expected hitting time of $i$ from a vertex at distance $1$ from $i$ equals $e_1 = 7$, so $\mathbb{E}_i[T_i] = 1 + e_1 = 8$.)*

---

### (ii) Expected number of visits to $o$ until the first return to $i$

For each $v \in V$, let $\gamma_v^i$ denote the expected number of visits to $v$ strictly before the first return to $i$:
$$\gamma_v^i = \mathbb{E}_i\!\left[\sum_{n=0}^{T_i - 1} \mathbf{1}_{\{X_n = v\}}\right].$$
Note that for $v = o \neq i$, we have $X_0 = i \neq o$ and $X_{T_i} = i \neq o$ on $\{T_i < \infty\}$ (which has probability $1$ by recurrence). Therefore
$$\sum_{n=0}^{T_i - 1} \mathbf{1}_{\{X_n = o\}} = \sum_{n=1}^{T_i - 1} \mathbf{1}_{\{X_n = o\}} = \sum_{n=1}^{T_i} \mathbf{1}_{\{X_n = o\}},$$
so there is no endpoint ambiguity in "until the first return to $i$".

By the standard IB theorem on invariant measures for an irreducible recurrent Markov chain:
- The vector $\gamma^i = (\gamma_v^i : v \in V)$ is an invariant measure of $P$ with $\gamma_i^i = 1$ and $0 < \gamma_v^i < \infty$ for all $v \in V$.
- Every invariant measure of an irreducible recurrent Markov chain is a positive scalar multiple of $\gamma^i$.

Since $\pi$ is also an invariant measure of $P$, there exists a constant $c > 0$ such that $\gamma_v^i = c \pi_v$ for all $v \in V$. Evaluating at $v = i$ gives $1 = \gamma_i^i = c \pi_i$, so $c = 1/\pi_i$ and
$$\gamma_v^i = \frac{\pi_v}{\pi_i} \qquad \text{for all } v \in V.$$
In particular, at $v = o$:
$$\gamma_o^i = \frac{\pi_o}{\pi_i} = \frac{1/8}{1/8} = 1.$$

*(Independent check via hitting probabilities on the distance chain $(Y_n)$ where $Y_n = d(X_n, i) \in \{0,1,2,3\}$: letting $h_d = \mathbb{P}(H_3 < H_0 \mid Y_0 = d)$ for $d \in \{1,2\}$, first-step analysis gives $h_1 = \frac{2}{3} h_2$ and $h_2 = \frac{2}{3} h_1 + \frac{1}{3}$, hence $h_1 = \frac{2}{5}$ and $h_2 = \frac{3}{5}$. Starting at $Y_0 = 0$, the particle jumps to $1$, reaches $3$ (vertex $o$) before returning to $0$ with probability $h_1 = \frac{2}{5}$, and from each visit to $3$ jumps to $2$ and returns to $3$ before hitting $0$ with probability $h_2 = \frac{3}{5}$. Thus, given at least one visit to $o$, the total number of visits to $o$ before returning to $i$ is $\mathrm{Geom}(1 - h_2) = \mathrm{Geom}(2/5)$ with mean $\frac{5}{2}$, giving expected visits $\frac{2}{5} \cdot \frac{5}{2} = 1$.)*

---

### (iii) Expected number of steps until the first visit to $o$

Let $Z_n = d(X_n, o) \in \{0, 1, 2, 3\}$ be the Hamming distance from $X_n$ to $o$. Then $X_n = o \iff Z_n = 0$, so
$$H_o = \inf\{n \ge 0 : X_n = o\} = \inf\{n \ge 0 : Z_n = 0\},$$
and $X_0 = i \iff Z_0 = 3$.

We first verify explicitly that $(Z_n)_{n \ge 0}$ is a Markov chain on $\{0, 1, 2, 3\}$ (lumpability):
- If $d(u, o) = 0$ (so $u = o$), all $3$ neighbours of $u$ have distance $1$ from $o$.
- If $d(u, o) = 1$, $u$ differs from $o$ in $1$ coordinate: flipping that coordinate moves to distance $0$ (probability $1/3$), while flipping either of the other $2$ coordinates moves to distance $2$ (probability $2/3$).
- If $d(u, o) = 2$, $u$ differs from $o$ in $2$ coordinates: flipping either of those moves to distance $1$ (probability $2/3$), while flipping the remaining coordinate moves to distance $3$ (probability $1/3$).
- If $d(u, o) = 3$ (so $u = i$), all $3$ neighbours of $u$ have distance $2$ from $o$.

Because $\mathbb{P}(Z_{n+1} = s \mid X_n = u)$ depends on $u$ only through $r = d(u, o)$ (call this value $\tilde{p}_{rs}$), for any $z_0, \dots, z_{n+1} \in \{0,1,2,3\}$ with $\mathbb{P}(Z_0 = z_0, \dots, Z_n = z_n) > 0$, conditioning on $X_n = u$ (with $d(u,o) = z_n$) and using the Markov property of $(X_n)$ gives
$$\mathbb{P}(Z_{n+1} = z_{n+1} \mid Z_0 = z_0, \dots, Z_n = z_n) = \sum_{u : d(u,o) = z_n} \tilde{p}_{z_n, z_{n+1}} \mathbb{P}(X_n = u \mid Z_0 = z_0, \dots, Z_n = z_n) = \tilde{p}_{z_n, z_{n+1}}.$$
Hence $(Z_n)_{n \ge 0}$ is a Markov chain on $\{0, 1, 2, 3\}$ with transition matrix
$$\tilde{P} = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 1/3 & 0 & 2/3 & 0 \\ 0 & 2/3 & 0 & 1/3 \\ 0 & 0 & 1 & 0 \end{pmatrix}.$$

Let $e_r = \mathbb{E}[H_o \mid Z_0 = r]$ for $r \in \{0, 1, 2, 3\}$. First note that $e_r < \infty$ for every $r$: from any state $r \in \{0,1,2,3\}$, $\mathbb{P}(H_o \le 3 \mid Z_0 = r) \ge (1/3)^3 = 1/27$, so by the Markov property $\mathbb{P}(H_o > 3m \mid Z_0 = r) \le (26/27)^m$ for all $m \ge 0$, and hence $e_r = \sum_{n=0}^\infty \mathbb{P}(H_o > n \mid Z_0 = r) \le 3 \sum_{m=0}^\infty (26/27)^m = 81 < \infty$.

By the theorem on expected hitting times (first-step analysis), $(e_0, e_1, e_2, e_3)$ is the minimal non-negative solution to
$$e_0 = 0,$$
$$e_1 = 1 + \frac{1}{3} e_0 + \frac{2}{3} e_2,$$
$$e_2 = 1 + \frac{2}{3} e_1 + \frac{1}{3} e_3,$$
$$e_3 = 1 + e_2.$$
We solve this linear system directly:
1. Substituting $e_3 = 1 + e_2$ into the equation for $e_2$:
   $$e_2 = 1 + \frac{2}{3} e_1 + \frac{1}{3}(1 + e_2) = \frac{4}{3} + \frac{2}{3} e_1 + \frac{1}{3} e_2 \implies \frac{2}{3} e_2 = \frac{4}{3} + \frac{2}{3} e_1 \implies e_2 = e_1 + 2.$$
2. Substituting $e_0 = 0$ and $e_2 = e_1 + 2$ into the equation for $e_1$:
   $$e_1 = 1 + \frac{2}{3}(e_1 + 2) = \frac{7}{3} + \frac{2}{3} e_1 \implies \frac{1}{3} e_1 = \frac{7}{3} \implies e_1 = 7.$$
3. Back-substituting gives $e_2 = 9$ and $e_3 = 10$.

Because this $3 \times 3$ linear system has a unique solution $(e_1, e_2, e_3) = (7, 9, 10)$ (which is strictly positive), it is necessarily the minimal non-negative solution. (Equivalently, on the original $8$-vertex chain, $h_v = e_{d(v,o)}$ solves $h_o = 0$, $h_v = 1 + \sum_{w} p_{vw} h_w$ for $v \neq o$, and any difference $u$ of two solutions is harmonic on $V \setminus \{o\}$ with $u_o = 0$, hence identically $0$ by the maximum principle on the finite connected cube graph.)

Therefore, the expected number of steps until the first visit to $o$ starting from $i$ is
$$\mathbb{E}_i[H_o] = e_3 = 10.$$

## Final Answers

- **(i)** Expected number of steps until the particle returns to $i$: **$8$**.
- **(ii)** Expected number of visits to $o$ until the first return to $i$: **$1$**.
- **(iii)** Expected number of steps until the first visit to $o$: **$10$**.

## Computational Verification

Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_qx01_check.py`

Verbatim stdout (`python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_qx01_check.py"`, exit code 0):
```text
IB/markov-chains/sheet-2/qx01 verification PASSED
Unique invariant distribution: pi_v = 1/8 for all 8 vertices (rank(I-P) = 7)
(i)   E_i[T_i] (exact 8-vertex solve) = 8 (1/pi_i = 8)
(ii)  gamma_o^i (exact 8-vertex solve) = 1 (pi_o/pi_i = 1)
(iii) E_i[H_o] (exact 8-vertex solve) = 10
      E_v[H_o] by distance to o: d=0 -> 0, d=1 -> 7, d=2 -> 9, d=3 -> 10
Hitting probabilities on distance chain: P_1(H_3 < H_0) = 2/5, P_2(H_3 < H_0) = 3/5
200-step exact partial sums: E_i[T_i] ~ 7.9999999978, gamma_o^i ~ 1.0000000000, E_i[H_o] ~ 9.9999999967
Exact partial + Markov tail: E_i[T_i] = 8, gamma_o^i = 1, E_i[H_o] = 10
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`. All three quantities were derived theoretically and verified independently in exact rational arithmetic both on the full 8-vertex transition matrix (without assuming symmetry) and via the lumped 4-state distance chain, hitting probabilities, and 200-step exact Markov-tail expansion.
- **Subtleties / Phase B checkpoints:**
  - Verify that the existing solution checks irreducibility and finiteness (giving positive recurrence and uniqueness of $\pi$) before invoking $m_i = 1/\pi_i$ and $\gamma_o^i = \pi_o/\pi_i$.
  - Watch out for any false claim that the chain is aperiodic (it has period $2$, being a random walk on a bipartite graph; aperiodicity is NOT needed for (i)–(iii)).
  - In (iii), verify that symmetry/lumpability (or uniqueness of the linear system for expected hitting times) is justified rather than asserted without explanation.
