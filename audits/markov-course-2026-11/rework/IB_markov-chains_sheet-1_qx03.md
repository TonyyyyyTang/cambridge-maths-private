---
id: IB/markov-chains/sheet-1/qx03
slug: IB_markov-chains_sheet-1_qx03
source_archive: private/raw/sheets/IB/markov-chains/2012-13/sheet-1.pdf
pages_checked: [1]
transcription: ok
blindness: clean
---

## Transcription Check

- Checked word-for-word and symbol-for-symbol against Question 1 on page 1 of `private/raw/sheets/IB/markov-chains/2012-13/sheet-1.pdf` (and the 200 DPI render `IB_markov-chains_sheet-1_qx03_p1.png`).
- Question 1 is at the top of page 1 and followed immediately by Question 2, so nothing is truncated across pages.
- **Result:** Exact match (`transcription: ok`).
- *Purely typographic note:* In the printed PDF, forced line breaks appear before "Show" and before "Deduce" (without paragraph indentation), whereas `question.tex` sets the three sentences in a single paragraph; all words, symbols ($P$ rather than $\mathbb{P}$, `\mid`, `\ldots`, `\bigcup_{n=1}^{\infty}`), labels `(a)`, `(b)`, and punctuation match identically.

## Full Solution

### Part 1: Constant conditional probability on a countable partition

Let $B_1, B_2, \ldots$ be pairwise disjoint events with $\bigcup_{n=1}^{\infty} B_n = \Omega$, and let $A$ be another event such that $P(A \mid B_n) = p$ for all $n$ (understood, in the standard way, for every $n \ge 1$ with $P(B_n) > 0$, as $P(A \mid B_n)$ is only defined when $P(B_n) > 0$).

First verify that for **every** $n \ge 1$,
$$
P(A \cap B_n) = p\, P(B_n).
$$
- If $P(B_n) > 0$, then by the definition of conditional probability,
  $$
  P(A \cap B_n) = P(A \mid B_n)\, P(B_n) = p\, P(B_n).
  $$
- If $P(B_n) = 0$, then $A \cap B_n \subseteq B_n$ implies $0 \le P(A \cap B_n) \le P(B_n) = 0$, so $P(A \cap B_n) = 0 = p\, P(B_n)$ holds automatically.

Now, because $B_1, B_2, \ldots$ are pairwise disjoint and $\bigcup_{n=1}^{\infty} B_n = \Omega$, the events $A \cap B_1, A \cap B_2, \ldots$ are pairwise disjoint and satisfy
$$
A = A \cap \Omega = A \cap \bigcup_{n=1}^{\infty} B_n = \bigcup_{n=1}^{\infty} (A \cap B_n).
$$
Applying countable additivity of the probability measure $P$ twice (first to $A = \bigsqcup_{n=1}^{\infty} (A \cap B_n)$, and then to $\Omega = \bigsqcup_{n=1}^{\infty} B_n$), we obtain:
$$
P(A) = \sum_{n=1}^{\infty} P(A \cap B_n) = \sum_{n=1}^{\infty} p\, P(B_n) = p \sum_{n=1}^{\infty} P(B_n) = p\, P(\Omega) = p \cdot 1 = p.
$$
*(Note: The same conclusion holds for a finite partition $B_1, \dots, B_m$ of $\Omega$, which is the special case obtained by setting $B_n = \emptyset$ for all $n > m$, since $P(\emptyset) = 0$.)*

---

### Part 2: Equivalence of (a) and (b) for discrete random variables

Let $X$ and $Y$ be discrete random variables, so there exist countable (finite or countably infinite) sets $S_X$ and $S_Y$ such that $P(X \in S_X) = 1$ and $P(Y \in S_Y) = 1$. Let
$$
S_Y^+ := \{y \in S_Y : P(Y = y) > 0\}
$$
be the support of $Y$. Then $S_Y^+$ is a non-empty countable set with $\sum_{y \in S_Y^+} P(Y = y) = 1$.

Recall the precise meanings of the two statements:
- **(a)** $X$ and $Y$ are independent: for all $x \in S_X$ and $y \in S_Y$,
  $$
  P(X = x, Y = y) = P(X = x)\, P(Y = y)
  $$
  (equivalently, by countable additivity, $P(X \in E, Y \in F) = P(X \in E)\, P(Y \in F)$ for all subsets $E \subseteq S_X$ and $F \subseteq S_Y$).
- **(b)** The conditional distribution of $X$ given $Y = y$ is independent of $y$: for each $x \in S_X$, the conditional probability $P(X = x \mid Y = y)$ takes the same value $p_x$ for all $y \in S_Y^+$ (i.e., all $y$ for which conditioning on $\{Y = y\}$ is well-defined).

We prove both implications.

#### Proof of (a) $\implies$ (b)
Assume $X$ and $Y$ are independent. Let $x \in S_X$ and let $y \in S_Y^+$ (so $P(Y = y) > 0$). By the definition of conditional probability and independence,
$$
P(X = x \mid Y = y) = \frac{P(X = x, Y = y)}{P(Y = y)} = \frac{P(X = x)\, P(Y = y)}{P(Y = y)} = P(X = x).
$$
The right-hand side does not depend on $y \in S_Y^+$. (Likewise, for any $E \subseteq S_X$, $P(X \in E \mid Y = y) = \sum_{x \in E} P(X = x \mid Y = y) = P(X \in E)$ is independent of $y$.) Hence the conditional distribution of $X$ given $Y = y$ is independent of $y$, proving (b).

#### Proof of (b) $\implies$ (a) (deduced from Part 1)
Assume (b): for each $x \in S_X$, there exists a constant $p_x \in [0, 1]$ such that
$$
P(X = x \mid Y = y) = p_x \qquad \text{for all } y \in S_Y^+.
$$
Construct a countable partition $B_1, B_2, \ldots$ of $\Omega$ as follows:
- Enumerate the countable set $S_Y^+$ as $\{y_1, y_2, \dots\}$ (finite or countably infinite) and set $B_n = \{Y = y_n\}$ for each element of $S_Y^+$.
- If $S_Y^+ = \{y_1, \dots, y_m\}$ is finite, define $B_{m+1} = \Omega \setminus \bigcup_{k=1}^m \{Y = y_k\}$ and $B_n = \emptyset$ for all $n > m + 1$.
- If $S_Y^+ = \{y_1, y_2, \dots\}$ is countably infinite, include $B_0 = \Omega \setminus \bigcup_{k=1}^{\infty} \{Y = y_k\}$ (re-indexing $B_0, B_1, B_2, \dots$ as a single sequence indexed by $\mathbb{N}$).

In every case, the events $B_1, B_2, \ldots$ are pairwise disjoint with $\bigcup_{n=1}^{\infty} B_n = \Omega$; every $B_n$ of the form $\{Y = y\}$ ($y \in S_Y^+$) satisfies $P(B_n) > 0$, while the remaining event(s) (the null complement and any empty padding events) have probability $0$ (since $P(\bigcup_{y \in S_Y^+} \{Y = y\}) = \sum_{y \in S_Y^+} P(Y = y) = 1$).

Now fix any $x \in S_X$ and set $A := \{X = x\}$. For every $n$ with $P(B_n) > 0$, $B_n = \{Y = y\}$ for some $y \in S_Y^+$, so by hypothesis (b),
$$
P(A \mid B_n) = P(X = x \mid Y = y) = p_x.
$$
By Part 1, we deduce that
$$
P(X = x) = P(A) = p_x.
$$
Finally, let $x \in S_X$ and $y \in S_Y$ be arbitrary:
- If $P(Y = y) > 0$ (i.e. $y \in S_Y^+$), then
  $$
  P(X = x, Y = y) = P(X = x \mid Y = y)\, P(Y = y) = p_x\, P(Y = y) = P(X = x)\, P(Y = y).
  $$
- If $P(Y = y) = 0$, then $\{X = x, Y = y\} \subseteq \{Y = y\}$ implies $0 \le P(X = x, Y = y) \le P(Y = y) = 0$, so
  $$
  P(X = x, Y = y) = 0 = P(X = x) \cdot 0 = P(X = x)\, P(Y = y).
  $$

Thus $P(X = x, Y = y) = P(X = x)\, P(Y = y)$ for all $x \in S_X$ and $y \in S_Y$. By countable additivity, for all $E \subseteq S_X$ and $F \subseteq S_Y$,
$$
P(X \in E, Y \in F) = \sum_{x \in E} \sum_{y \in F} P(X = x, Y = y) = \left(\sum_{x \in E} P(X = x)\right)\left(\sum_{y \in F} P(Y = y)\right) = P(X \in E)\, P(Y \in F),
$$
so $X$ and $Y$ are independent, proving (a).

## Final Answers

- **Part 1:** Proved $P(A) = \sum_{n=1}^{\infty} P(A \cap B_n) = p \sum_{n=1}^{\infty} P(B_n) = p\, P(\Omega) = p$ via countable additivity.
- **Part 2:** Proved (a) $\iff$ (b), using Part 1 directly in (b) $\implies$ (a) to deduce that the common conditional probability $p_x = P(X = x \mid Y = y)$ equals the marginal probability $P(X = x)$, whence $P(X = x, Y = y) = P(X = x) P(Y = y)$ for all $x, y$.

## Computational Verification

- Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-1_qx03_check.py`
- Command: `python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-1_qx03_check.py"` (exit code 0)
- Verbatim stdout:
```text
PASS: Part 1 verified on finite (with zero-measure event) and geometric countably infinite partitions.
PASS: Part 2 equivalence (a) <=> (b) verified on independent (with zero-prob y, degenerate Y) and dependent distributions.
```

## Confidence & Edge-Case Notes

- **Confidence:** High.
- **Subtleties / Phase B comparison watchpoints:**
  1. **Zero-probability events ($P(B_n) = 0$ and $P(Y = y) = 0$):** Conditioning on $\{Y = y\}$ is only defined for $y \in S_Y^+ = \{y : P(Y = y) > 0\}$. A rigorous proof should explicitly handle $P(B_n) = 0$ (where $P(A \cap B_n) = 0 = p P(B_n)$ automatically) and $P(Y = y) = 0$ (where $P(X = x, Y = y) = 0 = P(X = x) P(Y = y)$), rather than dividing by $P(Y = y)$ unconditionally.
  2. **Finite vs countably infinite support of $Y$:** Part 1 is stated for an infinite sequence $B_1, B_2, \ldots$, whereas a discrete random variable $Y$ may take only finitely many values. Padding with empty sets (or handling finite/infinite partitions uniformly via the $P(B_n) = 0$ remark) closes this gap cleanly.
  3. **Direct deduction:** Ensure (b) $\implies$ (a) genuinely invokes Part 1 (rather than merely reproving it without referencing Part 1).
