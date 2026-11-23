---
id: IB/markov-chains/sheet-2/qx05
slug: IB_markov-chains_sheet-2_qx05
source_archive: private/raw/sheets/IB/markov-chains/2012-13/sheet-2.pdf
pages_checked: [2]
transcription: ok
blindness: clean
---

## Transcription Check

**Status:** `ok` (exact match).

Visual inspection of Page 2 (Question 13) of `private/raw/sheets/IB/markov-chains/2012-13/sheet-2.pdf` (and confirmation that the **Extras** heading and Question 14 follow immediately on Page 2) shows a word-for-word, symbol-for-symbol, punctuation-for-punctuation match with `question.tex`:
- "An opera singer is due to perform a long series of concerts. Having a fine artistic temperament, she is liable to pull out each night with probability $1/2$. Once this has happened she will not sing again until the promoter convinces her of his high regard. This he does by sending flowers every day until she returns. Flowers costing $x$ thousand pounds, $0\le x\le1$, bring about a reconciliation with probability $\sqrt{x}$. The promoter stands to make £750 from each successful concert. How much should he spend on flowers?" — exact match.
- Purely typographic note: the printed sheet renders the pound sign in Computer Modern italic (`\pounds` = *£*750), while `question.tex` uses UTF-8 `£750` (identical meaning).

## Full Solution

### Markov Chain Formulation

Measure all monetary amounts in thousands of pounds, so a successful concert brings revenue $£750 = \frac{3}{4}$ (thousand pounds) and each daily bouquet of flowers costs $x \in [0, 1]$ (thousand pounds). Suppose the promoter chooses a fixed daily flower expenditure $x \in [0, 1]$ whenever the singer is refusing to sing.

Model the concert series as a discrete-time Markov chain $(X_n)_{n \ge 0}$ on the state space $I = \{0, 1\}$, where:
- $X_n = 1$ if the singer performs on night $n$ (a successful concert);
- $X_n = 0$ if the singer does not perform on night $n$ (she has pulled out and has not yet reconciled).

The day-to-day evolution is as follows:
- If $X_n = 1$ (she sings on night $n$), the promoter makes $\frac{3}{4}$ on night $n$ and sends no flowers the following day. On night $n+1$, she pulls out with probability $\frac{1}{2}$ (so $X_{n+1} = 0$) and performs again with probability $\frac{1}{2}$ (so $X_{n+1} = 1$).
- If $X_n = 0$ (she does not sing on night $n$), the promoter makes $0$ on night $n$ and sends flowers costing $x$ the following day in an attempt to reconcile her for night $n+1$. With probability $\sqrt{x}$ a reconciliation occurs and she returns to sing on night $n+1$ ($X_{n+1} = 1$); with probability $1 - \sqrt{x}$ she remains unreconciled and does not sing on night $n+1$ ($X_{n+1} = 0$).

Thus $(X_n)_{n \ge 0}$ is a Markov chain on $I = \{0, 1\}$ with transition matrix
$$P = \begin{pmatrix} p_{00} & p_{01} \\ p_{10} & p_{11} \end{pmatrix} = \begin{pmatrix} 1 - \sqrt{x} & \sqrt{x} \\ 1/2 & 1/2 \end{pmatrix}.$$

Let $R(1) = \frac{3}{4}$ and $R(0) = -x$ be the net profit (in thousands of pounds) associated with each visit to state $1$ or state $0$ respectively (noting that any boundary effect of at most one flower delivery at $n=0$ or $n=N$ changes the average over $N$ nights by at most $x/N \to 0$). Letting $V_i(N) = \sum_{n=0}^{N-1} \mathbf{1}_{\{X_n = i\}}$ for $i \in \{0,1\}$, the promoter's average net profit per night over the first $N$ nights is
$$\bar{R}_N = \frac{1}{N} \sum_{n=0}^{N-1} R(X_n) = \frac{3}{4} \frac{V_1(N)}{N} - x \frac{V_0(N)}{N}.$$

*(Renewal-cycle sanity check: starting from a concert she sings, she sings a run of $K \sim \mathrm{Geom}(1/2)$ consecutive nights with $\mathbb{E}[K] = 2$, pulls out on night $K+1$, and then receives flowers for $M \sim \mathrm{Geom}(\sqrt{x})$ consecutive days with $\mathbb{E}[M] = 1/\sqrt{x}$ before returning to sing on night $K+M+1$. Thus each cycle consists of $K + M$ nights containing $K$ successful concerts and $M$ flower deliveries, which matches associating $+3/4$ with each visit to state $1$ and $-x$ with each visit to state $0$.)*

---

### Long-Run Average Profit: Case Analysis in $x \in [0, 1]$

#### Case 1: $x = 0$ (degenerate / reducible case)

When $x = 0$, we have $\sqrt{x} = 0$ and
$$P = \begin{pmatrix} 1 & 0 \\ 1/2 & 1/2 \end{pmatrix}.$$
The chain is **not irreducible**: state $0$ is absorbing ($p_{00} = 1$) and state $1$ is transient (from $X_0 = 1$, she sings a finite number $K \sim \mathrm{Geom}(1/2)$ of nights and then remains in state $0$ forever). Solving $\pi P = \pi$ still gives the unique invariant distribution $\pi = (1, 0)$.

Directly from the sample paths:
- If $X_0 = 0$, then $X_n = 0$ for all $n \ge 0$, so $\bar{R}_N = 0$ for all $N \ge 1$.
- If $X_0 = 1$, then $V_1(N) \le K < \infty$ almost surely and the flower cost is $x = 0$, so $0 \le \bar{R}_N \le \frac{3K}{4N} \to 0$ almost surely (and $\mathbb{E}[\bar{R}_N] = \frac{3}{4N}\sum_{n=0}^{N-1} 2^{-n} = \frac{3}{4N}(2 - 2^{1-N}) \to 0$).

Hence for $x = 0$, the long-run average profit per night is $f(0) = 0$.

#### Case 2: $0 < x \le 1$ (irreducible case)

When $0 < x \le 1$, we have $\sqrt{x} \in (0, 1]$, so $p_{01} = \sqrt{x} > 0$ and $p_{10} = \frac{1}{2} > 0$. Therefore:
- The Markov chain $(X_n)$ is **irreducible**.
- Since the state space $I = \{0, 1\}$ is finite, the chain is **positive recurrent** and has a **unique invariant distribution** $\pi = (\pi_0, \pi_1)$.
- (Note also that $p_{11} = \frac{1}{2} > 0$, so the chain is aperiodic for all $x \in (0, 1]$, including the endpoint $x = 1$ where $p_{00} = 0$.)

The invariant distribution satisfies $\pi P = \pi$ and $\pi_0 + \pi_1 = 1$. From the first component of $\pi P = \pi$:
$$(1 - \sqrt{x})\pi_0 + \frac{1}{2}\pi_1 = \pi_0 \iff \sqrt{x}\,\pi_0 = \frac{1}{2}\pi_1 \iff \pi_1 = 2\sqrt{x}\,\pi_0.$$
Normalising $\pi_0 + \pi_1 = (1 + 2\sqrt{x})\pi_0 = 1$ gives
$$\pi_0 = \frac{1}{1 + 2\sqrt{x}}, \qquad \pi_1 = \frac{2\sqrt{x}}{1 + 2\sqrt{x}}.$$

By the **ergodic theorem** for irreducible positive recurrent Markov chains, for any initial distribution of $X_0$,
$$\frac{V_i(N)}{N} \xrightarrow{\text{a.s.}} \pi_i \qquad \text{as } N \to \infty, \quad i \in \{0, 1\}.$$
Hence, as $N \to \infty$, the average net profit per night converges almost surely (and in expectation, by bounded convergence since $|R(X_n)| \le 1$) to
$$f(x) = \frac{3}{4}\pi_1 - x\pi_0 = \frac{\frac{3}{2}\sqrt{x} - x}{1 + 2\sqrt{x}} \qquad \text{(in thousands of pounds per night)}.$$
Notice that this formula also gives $f(0) = 0$, so $f(x) = \frac{\frac{3}{2}\sqrt{x} - x}{1 + 2\sqrt{x}}$ is valid for all $x \in [0, 1]$.

---

### Maximising the Long-Run Average Profit

Set $u = \sqrt{x}$. As $x$ ranges over $[0, 1]$, $u$ ranges bijectively over $[0, 1]$ with $x = u^2$. Define
$$g(u) = f(u^2) = \frac{\frac{3}{2} u - u^2}{1 + 2u}, \qquad u \in [0, 1].$$

Since $1 + 2u \ge 1 > 0$ on $[0, 1]$, $g$ is smooth on $[0, 1]$. Differentiating via the quotient rule:
$$g'(u) = \frac{\left(\frac{3}{2} - 2u\right)(1 + 2u) - \left(\frac{3}{2} u - u^2\right) \cdot 2}{(1 + 2u)^2} = \frac{\frac{3}{2} + u - 4u^2 - 3u + 2u^2}{(1 + 2u)^2} = \frac{\frac{3}{2} - 2u - 2u^2}{(1 + 2u)^2}.$$
Factorising the numerator:
$$g'(u) = \frac{3 - 4u - 4u^2}{2(1 + 2u)^2} = \frac{(3 + 2u)(1 - 2u)}{2(1 + 2u)^2}.$$

Since $3 + 2u > 0$ and $2(1 + 2u)^2 > 0$ for all $u \in [0, 1]$, the sign of $g'(u)$ is the sign of $1 - 2u$:
- $g'(u) > 0$ for $u \in [0, 1/2)$ ($g$ is strictly increasing);
- $g'(1/2) = 0$;
- $g'(u) < 0$ for $u \in (1/2, 1]$ ($g$ is strictly decreasing).

Hence $g$ attains its **unique global maximum** on $[0, 1]$ at $u = \frac{1}{2}$.

*(Elementary global-maximum verification without calculus: for every $u \in [0, 1]$,*
$$\frac{1}{4} - g(u) = \frac{(1 + 2u) - 4\left(\frac{3}{2} u - u^2\right)}{4(1 + 2u)} = \frac{4u^2 - 4u + 1}{4(1 + 2u)} = \frac{(2u - 1)^2}{4(1 + 2u)} \ge 0,$$
*with equality if and only if $u = \frac{1}{2}$.)*

Comparing key values:
- At $u = 0$ ($x = 0$, £0): $g(0) = 0$ (£0 per night).
- At $u = \frac{1}{2}$ ($x = \frac{1}{4}$, £250): $g(1/2) = \frac{1}{4}$ (£250 per night, with $\pi_0 = \pi_1 = \frac{1}{2}$).
- At $u = 1$ ($x = 1$, £1000): $g(1) = \frac{1}{6}$ (£166.67 per night, with $\pi_0 = \frac{1}{3}$, $\pi_1 = \frac{2}{3}$).

Therefore, the optimal spend on flowers is
$$x = u^2 = \left(\frac{1}{2}\right)^2 = \frac{1}{4} \text{ thousand pounds} = £250$$
each day that he sends flowers.

*(Remark on dynamic policies: $x = \frac{1}{4}$ is optimal even among all time-varying or history-dependent spending policies $x_n \in [0, 1]$. Indeed, setting $V(0) = 0$ and $V(1) = 1$, we have*
$$R(1) + \mathbb{E}[V(X_{n+1}) \mid X_n = 1] - V(1) = \frac{3}{4} + \frac{1}{2} - 1 = \frac{1}{4},$$
*and for any choice $x \in [0, 1]$ in state $0$,*
$$R(0, x) + \mathbb{E}[V(X_{n+1}) \mid X_n = 0, x] - V(0) = -x + \sqrt{x} = \frac{1}{4} - \left(\sqrt{x} - \frac{1}{2}\right)^2 \le \frac{1}{4},$$
*with equality iff $x = \frac{1}{4}$. Summing over $n = 0, \dots, N-1$ and taking expectations gives $\frac{1}{N}\mathbb{E}[\sum_{n=0}^{N-1} R_n] \le \frac{1}{4} + \frac{V(X_0) - \mathbb{E}[V(X_N)]}{N} \to \frac{1}{4}$, so no policy can beat the long-run average profit of £250 per night achieved by constant $x = \frac{1}{4}$.)*

## Final Answers

- The promoter should spend **$x = \frac{1}{4}$ thousand pounds, i.e., £250**, on flowers each day (until she returns).
- At this optimal spend, $\sqrt{x} = \frac{1}{2}$, the invariant distribution is $\pi = (\frac{1}{2}, \frac{1}{2})$, and the maximum long-run average profit is **£250 per night** (with long-run average flower spend £125 per night).

## Computational Verification

Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_qx05_check.py`

Verbatim stdout (`python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-2_qx05_check.py"`, exit code 0):
```text
IB/markov-chains/sheet-2/qx05 verification PASSED
Optimal u = sqrt(x) = 1/2, optimal x = 1/4 thousand pounds = £250
At optimum: pi = (1/2, 1/2), long-run profit = 1/4 thousand pounds/night = £250/night
Endpoint values: g(0) = 0 (£0/night), g(1) = 1/6 (£166.67/night)
Global-maximum identity 1/4 - g(u) = (2u-1)^2 / (4(1+2u)) and Bellman inequality verified on 104 rational points.
```

## Confidence & Edge-Case Notes

- **Confidence:** `high`. The Markov chain model matches the renewal-cycle accounting identically, the optimisation factorises cleanly to $(3+2u)(1-2u)$ with exact algebraic identity $\frac{1}{4} - g(u) = \frac{(2u-1)^2}{4(1+2u)}$, and the Bellman remark proves optimality even over non-stationary strategies.
- **Subtleties / Phase B checkpoints:**
  - **Edge case $x = 0$:** Watch for solutions that claim $(X_n)$ is irreducible for all $x \in [0, 1]$ without separating out $x = 0$ (where state $0$ is absorbing and state $1$ is transient, though $\pi = (1, 0)$ is still the unique invariant distribution and the long-run profit is still $0$).
  - **Maximisation justification:** Verify that the stationary point $u = 1/2$ is proved to be a global maximum on $[0, 1]$ (via sign of $g'$, second derivative + endpoints, or the completing-the-square identity), not merely identified from $g'(u) = 0$.
  - **Units:** Check that $x = 1/4$ is converted back to **£250** (since $x$ is in thousands of pounds, and concert profit is £750 $= 3/4$ thousand pounds).
