---
id: IB/markov-chains/sheet-2/qx06
slug: IB_markov-chains_sheet-2_qx06
source_archive: private/raw/sheets/IB/markov-chains/2012-13/sheet-2.pdf
pages_checked: [2, 3]
transcription: ok
blindness: clean
---

## Transcription Check

**Status:** `ok` (exact match).

Visual inspection of Page 2 (Question 14, under the **Extras** section) and Page 3 (which begins with Question 15, confirming nothing of Question 14 spills over) of `private/raw/sheets/IB/markov-chains/2012-13/sheet-2.pdf` confirms a word-for-word, symbol-for-symbol, punctuation-for-punctuation match with `question.tex`:
- "We wish to simulate a random variable $X$ for which we know $P(X=i)\propto\alpha_i$, $i\in I$, but we do not know the constant of proportionality because $\sum_{i\in I}\alpha_i$ is too difficult calculate." — exact match! Notably, the printed sheet itself contains the author typo `"too difficult calculate"` (missing `"to"`), which `question.tex` faithfully preserves verbatim.
- "Let $P=(p_{ij})$ be the (known) transition matrix of some irreducible Markov chain on $I$. Define" — exact match.
- Displayed equation:
  `q_{ij}=\min\{p_{ij},\alpha_jp_{ji}/\alpha_i\},\ i\ne j,\quad\text{and}\quad q_{ii}=1-\sum_{j\ne i}q_{ij}.` — exact match, including the trailing full stop.
- "Prove that $Q=(q_{ij})$ is the transition matrix of a reversible Markov chain with invariant distribution proportional to $\alpha$. How can this help with the problem of simulating $X$?" — exact match.
- Note: the **Extras** heading and its introductory note ("The following questions are optional, ...") are section-level framing for Questions 14–17, not part of Question 14 itself.

## Full Solution

### Preliminary Hypotheses

Since $X$ is a discrete random variable taking values in the countable state space $I$ with $P(X = i) \propto \alpha_i$ for $i \in I$, there exists a normalising constant $Z = \sum_{k \in I} \alpha_k \in (0, \infty)$ such that
$$\pi_i := P(X = i) = \frac{\alpha_i}{Z} \qquad (i \in I),$$
with $\alpha_i \ge 0$ for all $i \in I$. Moreover, for the ratio $\alpha_j p_{ji} / \alpha_i$ in the definition of $q_{ij}$ to be well-defined for every $i \in I$, we implicitly assume that $\alpha_i > 0$ for all $i \in I$ (equivalently, $I$ is taken to be the support $\{i : P(X = i) > 0\}$ of $X$). All arguments below apply equally whether $I$ is finite or countably infinite.

---

### Part 1: $Q = (q_{ij})$ is a Transition (Stochastic) Matrix

We must show that for every $i \in I$: (a) $q_{ij} \ge 0$ for all $j \in I$, and (b) $\sum_{j \in I} q_{ij} = 1$.

1. **Non-negativity of off-diagonal entries:** For any $i, j \in I$ with $i \neq j$, since $P$ is a stochastic matrix we have $p_{ij} \ge 0$ and $p_{ji} \ge 0$, and since $\alpha_i, \alpha_j > 0$ we have $\alpha_j p_{ji} / \alpha_i \ge 0$. Hence
   $$q_{ij} = \min\!\left\{p_{ij}, \frac{\alpha_j p_{ji}}{\alpha_i}\right\} \ge 0.$$
2. **Upper bound on off-diagonal row sum and non-negativity of $q_{ii}$:** By definition of the minimum,
   $$0 \le q_{ij} \le p_{ij} \qquad \text{for all } j \neq i.$$
   Summing over all $j \neq i$ (all terms are non-negative, so the series converges in $[0, \infty]$ and comparison applies even when $I$ is countably infinite):
   $$0 \le \sum_{j \neq i} q_{ij} \le \sum_{j \neq i} p_{ij} = 1 - p_{ii} \le 1,$$
   where we used that $P$ is a stochastic matrix ($\sum_{j \in I} p_{ij} = 1$ and $p_{ii} \ge 0$). Consequently, the sum $\sum_{j \neq i} q_{ij}$ is finite and lies in $[0, 1]$, so $q_{ii} = 1 - \sum_{j \neq i} q_{ij}$ is well-defined and satisfies
   $$q_{ii} = 1 - \sum_{j \neq i} q_{ij} \ge 1 - (1 - p_{ii}) = p_{ii} \ge 0$$
   (and $q_{ii} \le 1$).
3. **Unit row sums:** For every $i \in I$,
   $$\sum_{j \in I} q_{ij} = q_{ii} + \sum_{j \neq i} q_{ij} = \left(1 - \sum_{j \neq i} q_{ij}\right) + \sum_{j \neq i} q_{ij} = 1.$$

Therefore $Q = (q_{ij})_{i,j \in I}$ is a stochastic (transition) matrix on $I$.

---

### Part 2: Detailed Balance, Reversibility, and Invariant Distribution $\pi \propto \alpha$

Define the probability distribution $\pi = (\pi_i : i \in I)$ by $\pi_i = \alpha_i / Z$, where $Z = \sum_{k \in I} \alpha_k \in (0, \infty)$. Then $\pi_i > 0$ for all $i \in I$, $\sum_{i \in I} \pi_i = 1$, and $\pi$ is proportional to $\alpha$ (indeed $\pi_i = P(X = i)$).

We verify the **detailed balance equations** $\pi_i q_{ij} = \pi_j q_{ji}$ for all $i, j \in I$:
- If $i = j$, $\pi_i q_{ii} = \pi_i q_{ii}$ holds trivially.
- If $i \neq j$, since $\alpha_i > 0$ we can multiply inside the minimum:
  $$\alpha_i q_{ij} = \alpha_i \min\!\left\{p_{ij}, \frac{\alpha_j p_{ji}}{\alpha_i}\right\} = \min\{\alpha_i p_{ij}, \alpha_j p_{ji}\}.$$
  Because $\min\{\alpha_i p_{ij}, \alpha_j p_{ji}\}$ is symmetric in $i$ and $j$, it follows immediately that
  $$\alpha_i q_{ij} = \min\{\alpha_i p_{ij}, \alpha_j p_{ji}\} = \min\{\alpha_j p_{ji}, \alpha_i p_{ij}\} = \alpha_j q_{ji}.$$
  Dividing both sides by $Z > 0$ yields
  $$\pi_i q_{ij} = \pi_j q_{ji} \qquad \text{for all } i, j \in I.$$

We deduce directly from detailed balance that $\pi$ is an **invariant distribution** of $Q$: for every $j \in I$, summing over $i \in I$ (valid for countable $I$ as all terms are non-negative),
$$(\pi Q)_j = \sum_{i \in I} \pi_i q_{ij} = \sum_{i \in I} \pi_j q_{ji} = \pi_j \sum_{i \in I} q_{ji} = \pi_j \cdot 1 = \pi_j,$$
since $Q$ is stochastic. Thus $\pi Q = \pi$.

Furthermore, because $\pi$ is an invariant distribution satisfying detailed balance with $Q$, the Markov chain $(Y_n)_{n \ge 0}$ with transition matrix $Q$ in equilibrium ($Y_0 \sim \pi$) is **reversible**: for any $N \ge 1$ and $i_0, i_1, \dots, i_N \in I$,
$$\mathbb{P}(Y_0 = i_0, Y_1 = i_1, \dots, Y_N = i_N) = \pi_{i_0} q_{i_0 i_1} \cdots q_{i_{N-1} i_N} = \pi_{i_N} q_{i_N i_{N-1}} \cdots q_{i_1 i_0} = \mathbb{P}(Y_0 = i_N, Y_1 = i_{N-1}, \dots, Y_N = i_0),$$
so $(Y_{N-n})_{0 \le n \le N}$ has the same joint distribution as $(Y_n)_{0 \le n \le N}$. Hence $Q$ is the transition matrix of a reversible Markov chain with invariant distribution $\pi \propto \alpha$.

---

### Part 3: How This Helps Simulate $X$ (Metropolis–Hastings / MCMC)

#### 1. Simulating transitions of $Q$ without knowing $Z = \sum_{i \in I} \alpha_i$

The key advantage of $Q$ is that its off-diagonal transition probabilities depend on $\alpha$ **only through the ratios** $\alpha_j / \alpha_i$, so the intractable normalising constant $Z = \sum_{k \in I} \alpha_k$ cancels completely. Moreover, we do not even need to evaluate the sum $\sum_{j \neq i} q_{ij}$ to compute $q_{ii}$!

Instead, we can simulate one step of a Markov chain $(Y_n)_{n \ge 0}$ with transition matrix $Q$ from current state $Y_n = i$ via the **Metropolis–Hastings** proposal/accept–reject rule:
1. **Propose** a candidate state $j \in I$ according to the known proposal distribution $(p_{ij})_{j \in I}$.
2. **Accept or reject:**
   - If $j = i$, set $Y_{n+1} = i$.
   - If $j \neq i$ (note $p_{ij} > 0$ almost surely when $j$ is proposed from $p_{i\cdot}$), compute the acceptance probability
     $$a_{ij} = \min\!\left\{1, \frac{\alpha_j p_{ji}}{\alpha_i p_{ij}}\right\} \in [0, 1].$$
     Independently draw $U \sim \mathrm{Uniform}(0, 1)$: if $U \le a_{ij}$, **accept** and set $Y_{n+1} = j$; otherwise **reject** and set $Y_{n+1} = i$.

Let us verify that this two-step mechanism has exact transition matrix $Q$:
- For $j \neq i$: if $p_{ij} = 0$, the transition probability is $0 = q_{ij}$; if $p_{ij} > 0$, the probability of moving from $i$ to $j$ is
  $$p_{ij} a_{ij} = p_{ij} \min\!\left\{1, \frac{\alpha_j p_{ji}}{\alpha_i p_{ij}}\right\} = \min\!\left\{p_{ij}, \frac{\alpha_j p_{ji}}{\alpha_i}\right\} = q_{ij}.$$
- For $j = i$: the probability of remaining at $i$ is the complement of moving to any $j \neq i$, namely
  $$p_{ii} + \sum_{j \neq i} p_{ij}(1 - a_{ij}) = p_{ii} + \sum_{j \neq i}(p_{ij} - q_{ij}) = 1 - \sum_{j \neq i} q_{ij} = q_{ii}.$$

Thus we can simulate $(Y_n)_{n \ge 0}$ using only samples from $P$ and evaluations of $\alpha_j / \alpha_i$.

#### 2. Convergence to $P(X = i) = \pi_i$ — Crucial Conditions on Irreducibility and Aperiodicity of $Q$

To use $(Y_n)$ to approximately sample $X$, we want $\mathbb{P}(Y_n = j) \to \pi_j = P(X = j)$ as $n \to \infty$ regardless of the initial state $Y_0$. By the standard **convergence to equilibrium theorem** for Markov chains on a countable state space $I$, since $Q$ already has the invariant distribution $\pi$, it suffices that **$Q$ is irreducible and aperiodic** (irreducibility + existence of an invariant distribution $\pi$ implies $Q$ is positive recurrent and $\pi$ is the **unique** invariant distribution; adding aperiodicity implies $\mathbb{P}(Y_n = j) \to \pi_j$ for all $j \in I$ and any initial distribution).

**Crucial caveat:** Merely assuming that $P$ is irreducible does **not** automatically guarantee either irreducibility or aperiodicity of $Q$:

1. **Irreducibility of $Q$:** For $i \neq j$, since $\alpha_i, \alpha_j > 0$, we have
   $$q_{ij} > 0 \iff p_{ij} > 0 \text{ and } p_{ji} > 0.$$
   Hence $Q$ is irreducible if and only if the undirected graph of **bidirectional** transitions $\{(i,j) : p_{ij} > 0 \text{ and } p_{ji} > 0\}$ is connected.
   - *Counterexample if $P$ is only irreducible:* On $I = \{1, 2, 3\}$, the directed 3-cycle
     $$P = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$$
     is irreducible, yet $p_{ij} p_{ji} = 0$ for all $i \neq j$, so $q_{ij} = 0$ for all $i \neq j$ and $Q = I_3$ is completely reducible (every state is absorbing and $Q^n = I_3$ never converges to $\pi$).
   - *How to ensure irreducibility of $Q$:* Choose $P$ to be irreducible with **symmetric support**, i.e. $p_{ij} > 0 \iff p_{ji} > 0$ for all $i, j \in I$ (for instance, a symmetric proposal $p_{ij} = p_{ji}$ such as a random walk on a connected undirected graph on $I$). Then $Q$ has the same positive off-diagonal transitions as $P$, so $Q$ is irreducible.

2. **Aperiodicity of $Q$:** Even if $P$ is symmetric and irreducible (so $Q$ is irreducible), $Q$ can still be **periodic**.
   - *Counterexample:* On $I = \{1, 2\}$, take $P = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ and $\alpha = (1, 1)$. Then $Q = P$ is irreducible and reversible with unique invariant distribution $\pi = (\frac{1}{2}, \frac{1}{2})$, but has period $2$, so $\mathbb{P}(Y_n = 1)$ alternates between $0$ and $1$ and fails to converge.
   - *Characterisation and how to ensure aperiodicity of $Q$:* Suppose $Q$ is irreducible. Since $q_{ii} = p_{ii} + \sum_{j \neq i}(p_{ij} - q_{ij}) \ge p_{ii}$, if $q_{ii} > 0$ for even a single state $i \in I$, then $Q$ is aperiodic. Notice that $q_{ii} = 0$ for all $i$ occurs if and only if $p_{ii} = 0$ for all $i$ AND $q_{ij} = p_{ij}$ for all $i \neq j$, i.e., $\alpha_i p_{ij} = \alpha_j p_{ji}$ for all $i, j$ (so $P$ is already $\alpha$-reversible and $Q = P$). Thus:
     - Whenever $P$ is **not** already $\alpha$-reversible, an irreducible $Q$ automatically acquires at least one self-loop $q_{ii} > 0$ (from rejected proposals) and is therefore aperiodic!
     - In general, choosing $P$ with $p_{ii} > 0$ for at least one $i \in I$ (for example, using a lazy proposal $\tilde{P} = \frac{1}{2}(I + P)$, which gives $\tilde{q}_{ii} \ge \frac{1}{2} > 0$ for all $i$) guarantees that $Q$ is aperiodic unconditionally.

#### 3. Using the Simulated Chain in Practice

With $P$ chosen so that $Q$ is irreducible and aperiodic:
- **Approximate simulation of $X$:** Start $(Y_n)$ from any convenient initial state $Y_0 = i_0$ and run the Metropolis–Hastings chain for a large ("burn-in") number of steps $N$. Since $\mathbb{P}(Y_N = j) \to \pi_j = P(X = j)$ as $N \to \infty$, the state $Y_N$ has distribution close to that of $X$ and serves as an approximate sample of $X$. (To obtain multiple approximately independent samples, one may run independent copies of the chain or take states spaced far apart.)
- **Estimating expectations $\mathbb{E}[f(X)]$ (ergodic averages):** Even without aperiodicity (requiring only irreducibility and positive recurrence of $Q$), by the **ergodic theorem**, for any bounded function $f: I \to \mathbb{R}$ (or any $f$ with $\sum_{i \in I} |f(i)| \pi_i < \infty$),
  $$\frac{1}{n} \sum_{k=0}^{n-1} f(Y_k) \xrightarrow{\text{a.s.}} \sum_{i \in I} f(i) \pi_i = \mathbb{E}[f(X)] \qquad \text{as } n \to \infty,$$
  and in particular the empirical proportions $\frac{1}{n}\sum_{k=0}^{n-1} \mathbf{1}_{\{Y_k = j\}} \xrightarrow{\text{a.s.}} P(X = j)$.

## Final Answers

- $Q = (q_{ij})$ is a valid stochastic matrix ($q_{ij} \ge 0$ for $j \neq i$, $0 \le \sum_{j \neq i} q_{ij} \le 1 - p_{ii} \le 1 \implies q_{ii} \ge p_{ii} \ge 0$, and $\sum_j q_{ij} = 1$).
- Taking $\pi_i = \alpha_i / Z$ with $Z = \sum_{k \in I} \alpha_k \in (0, \infty)$ (assuming $\alpha_i > 0$ on $I$), detailed balance holds since $\alpha_i q_{ij} = \min\{\alpha_i p_{ij}, \alpha_j p_{ji}\} = \alpha_j q_{ji}$, so $\pi_i q_{ij} = \pi_j q_{ji}$ for all $i, j \in I$; summing over $i$ gives $\pi Q = \pi$, so $Q$ is the transition matrix of a reversible Markov chain with invariant distribution $\pi \propto \alpha$.
- **Application to simulating $X$:** One simulates transitions of $Q$ via the Metropolis–Hastings proposal ($j \sim p_{i\cdot}$) and accept–reject step ($a_{ij} = \min\{1, \frac{\alpha_j p_{ji}}{\alpha_i p_{ij}}\}$), which requires only ratios $\alpha_j / \alpha_i$ and never needs $Z = \sum_{i \in I} \alpha_i$. Provided $P$ is chosen so that $Q$ is **irreducible** (requires connected bidirectional edges of $P$, e.g. symmetric support $p_{ij} > 0 \iff p_{ji} > 0$; irreducibility of $P$ alone is *insufficient*, e.g. a directed 3-cycle gives $Q = I_3$) and **aperiodic** (e.g. $p_{ii} > 0$ for some $i$, or a lazy proposal; otherwise e.g. $P = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ with $\alpha = (1,1)$ gives period-2 $Q$), $Q$ is positive recurrent with unique invariant distribution $\pi$ and $\mathbb{P}(Y_n = j) \to \pi_j = P(X = j)$, so $Y_N$ for large $N$ approximates a draw from $X$, and ergodic averages $\frac{1}{n}\sum_{k=0}^{n-1} f(Y_k) \to \mathbb{E}[f(X)]$ a.s.

## Computational Verification

Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_qx06_check.py`

Verbatim stdout (`python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_qx06_check.py"`, exit code 0):
```text
IB/markov-chains/sheet-2/qx06 verification PASSED
Example 1 (4-state asymmetric P, sym support): irreducible=True, period=1, TV(10)=5.808495e-04, TV(30)=3.477403e-10
Example 2 (directed 3-cycle P): P irreducible=True, Q=I_3 irreducible=False, TV(50)=5/6
Example 3 (2 sym pairs + directed bridge): P irreducible=True, Q irreducible=False (classes {0,1}, {2,3})
Example 4 (2-cycle, alpha=(1,1)): Q=P irreducible=True, period=2, TV(50)=1/2
Example 5 (2-cycle, alpha=(1,2)): Q irreducible=True, period=1, TV(20)=1/1572864 (6.357829e-07)
Example 6 (lazy 2-cycle, alpha=(1,1)): Q irreducible=True, period=1, TV(1)=0
Example 7 (truncated Geom(2/3), m=7, periodic P): Q irreducible=True, period=1, TV(50)=5.771300e-03
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`. Every step of the proof holds for both finite and countably infinite state spaces $I$, and every subtle hypothesis and counterexample was verified in exact rational arithmetic.
- **Subtleties / Phase B checkpoints:**
  - **Transcription note:** The original printed sheet has the typo `"too difficult calculate."` (missing `"to"`), which `question.tex` faithfully reproduces.
  - **Major pitfall #1 (Irreducibility of $Q$):** Watch very carefully for solutions that claim "$Q$ is irreducible because $P$ is irreducible" — that is **FALSE** (counterexample: directed 3-cycle $P$ gives $Q = I_3$!). One needs bidirectional connectivity ($p_{ij} > 0 \text{ and } p_{ji} > 0$, e.g. symmetric support $p_{ij} > 0 \iff p_{ji} > 0$).
  - **Major pitfall #2 (Aperiodicity of $Q$):** Watch for solutions that claim $\mathbb{P}(Y_n = j) \to \pi_j$ without checking aperiodicity of $Q$ (counterexample: symmetric 2-cycle with $\alpha = (1, 1)$ gives periodic $Q = P$). Also note the sharp characterisation: an irreducible $Q$ can be periodic *only* if $p_{ii} = 0$ for all $i$ and $P$ already satisfies detailed balance with $\alpha$.
  - **Major pitfall #3 (Implicit assumptions on $\alpha$):** Check whether $\alpha_i > 0$ for all $i \in I$ (needed to avoid division by zero in $\alpha_j p_{ji}/\alpha_i$) and $\sum_{i \in I} \alpha_i < \infty$ (needed so $\alpha$ normalises to a probability distribution, implied by $P(X = i) \propto \alpha_i$) are noted.
  - **Major pitfall #4 (How to simulate $Q$ without computing $q_{ii}$):** Watch whether the solution explains the two-step proposal + accept–reject mechanism ($j \sim p_{i\cdot}$, accept with probability $\min\{1, \alpha_j p_{ji}/(\alpha_i p_{ij})\}$) rather than pretending one can compute $q_{ii} = 1 - \sum_{j \neq i} q_{ij}$ explicitly when $I$ is large/infinite.
