---
id: IB/markov-chains/sheet-1/qx12
slug: IB_markov-chains_sheet-1_qx12
source_archive: private/raw/sheets/IB/markov-chains/2010-11/sheet-1.pdf
pages_checked: [2, 3]
transcription: ok
blindness: clean
---

## Transcription Check

- Checked word-for-word against Question 13 on page 3 of `private/raw/sheets/IB/markov-chains/2010-11/sheet-1.pdf` (and the 200 DPI render `IB_markov-chains_sheet-1_qx12_p3.png`).
- Question 13 is the final question on page 3; nothing follows it.
- **Result:** Exact verbatim match to Question 13 (`transcription: ok`).
- **Important self-containedness / metadata flag:** Question 13 reads *"In Example 8, which states are recurrent and which are transient?"* and depends directly on Question 8 ("Example 8") printed on **page 2** of `private/raw/sheets/IB/markov-chains/2010-11/sheet-1.pdf`:
  > **8.** Identify the communicating classes of the transition matrix
  > $$P = \begin{pmatrix} \frac{1}{2} & 0 & 0 & 0 & \frac{1}{2} \\ 0 & \frac{1}{2} & 0 & \frac{1}{2} & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 0 & \frac{1}{4} & \frac{1}{4} & \frac{1}{4} & \frac{1}{4} \\ \frac{1}{2} & 0 & 0 & 0 & \frac{1}{2} \end{pmatrix}.$$
  > Which of the classes are closed?
  Note that `question.tex` reproduces Question 13 verbatim without quoting $P$ from Question 8, and the packet locator lists only `pages: [3]` rather than `[2, 3]`. If rendered as a standalone unit without linking/quoting Example 8, a reader would not see the matrix $P$.

## Full Solution

Label the state space $I = \{1, 2, 3, 4, 5\}$ corresponding to the rows and columns of the $5 \times 5$ transition matrix $P$ from Example 8 (page 2 of `2010-11/sheet-1.pdf`):
$$
P = \begin{pmatrix}
\frac{1}{2} & 0 & 0 & 0 & \frac{1}{2} \\
0 & \frac{1}{2} & 0 & \frac{1}{2} & 0 \\
0 & 0 & 1 & 0 & 0 \\
0 & \frac{1}{4} & \frac{1}{4} & \frac{1}{4} & \frac{1}{4} \\
\frac{1}{2} & 0 & 0 & 0 & \frac{1}{2}
\end{pmatrix}.
$$
(All entries are non-negative and every row sums to $1$.)

### Step 1: Communicating classes and closedness

Listing the positive one-step transitions ($p_{ij} > 0$):
- From $1$: $1 \to 1$ ($\frac{1}{2}$), $1 \to 5$ ($\frac{1}{2}$).
- From $2$: $2 \to 2$ ($\frac{1}{2}$), $2 \to 4$ ($\frac{1}{2}$).
- From $3$: $3 \to 3$ ($1$).
- From $4$: $4 \to 2$ ($\frac{1}{4}$), $4 \to 3$ ($\frac{1}{4}$), $4 \to 4$ ($\frac{1}{4}$), $4 \to 5$ ($\frac{1}{4}$).
- From $5$: $5 \to 1$ ($\frac{1}{2}$), $5 \to 5$ ($\frac{1}{2}$).

From these transitions we identify the communicating classes ($i \leftrightarrow j$ iff $i$ leads to $j$ and $j$ leads to $i$):
1. **$\{1, 5\}$:** Since $p_{15} = \frac{1}{2} > 0$ and $p_{51} = \frac{1}{2} > 0$, we have $1 \leftrightarrow 5$. Moreover, $p_{1j} = p_{5j} = 0$ for all $j \notin \{1, 5\}$, so no state in $\{1, 5\}$ leads to any state outside $\{1, 5\}$. Thus $\{1, 5\}$ is a **closed communicating class** (finite, size $2$).
2. **$\{3\}$:** Since $p_{33} = 1$ and $p_{3j} = 0$ for all $j \ne 3$, state $3$ is absorbing and leads only to itself. Thus $\{3\}$ is a **closed communicating class** (finite, size $1$).
3. **$\{2, 4\}$:** Since $p_{24} = \frac{1}{2} > 0$ and $p_{42} = \frac{1}{4} > 0$, we have $2 \leftrightarrow 4$. Neither $2$ nor $4$ can communicate with $1, 3,$ or $5$ because $\{1, 5\}$ and $\{3\}$ are closed. However, $p_{43} = \frac{1}{4} > 0$ and $p_{45} = \frac{1}{4} > 0$ (and $p^{(2)}_{23} = p_{24} p_{43} = \frac{1}{8} > 0$), so $\{2, 4\}$ can escape into $\{3\}$ and $\{1, 5\}$. Thus $\{2, 4\}$ is a **non-closed (open) communicating class**.

### Step 2: Classification of recurrent and transient states

#### Method 1: Via standard Part IB theorems
Recall the two standard results for discrete-time Markov chains:
- **Theorem 1:** Recurrence and transience are class properties, and every recurrent communicating class is closed (equivalently, every state in a non-closed communicating class is transient).
- **Theorem 2:** Every finite closed communicating class is recurrent.

Checking the hypotheses on each communicating class:
- $\{1, 5\}$ is closed and finite ($|\{1, 5\}| = 2 < \infty$) $\implies$ **states $1$ and $5$ are recurrent**.
- $\{3\}$ is closed and finite ($|\{3\}| = 1 < \infty$) $\implies$ **state $3$ is recurrent**.
- $\{2, 4\}$ is not closed $\implies$ **states $2$ and $4$ are transient**.

*(Moreover, since each recurrent class is finite, states $1, 3, 5$ are in fact **positive recurrent**, with mean return times $m_3 = 1$ and $m_1 = m_5 = 2$.)*

#### Method 2: Direct verification from the definition of recurrence and transience
Let $T_i = \inf\{n \ge 1 : X_n = i\}$ and $f_i = \mathbb{P}_i(T_i < \infty)$; state $i$ is recurrent iff $f_i = 1$ (equivalently $\sum_{n=0}^{\infty} p^{(n)}_{ii} = \infty$) and transient iff $f_i < 1$ (equivalently $\sum_{n=0}^{\infty} p^{(n)}_{ii} < \infty$).

1. **State $3$:** Since $p_{33} = 1$, $\mathbb{P}_3(T_3 = 1) = 1$, so $f_3 = 1$ (and $\sum_{n=0}^{\infty} p^{(n)}_{33} = \sum_{n=0}^{\infty} 1 = \infty$). Hence **$3$ is recurrent**.
2. **States $1$ and $5$:** Starting from $1$, the chain stays in $\{1, 5\}$ forever and transitions to $1$ or $5$ each with probability $\frac{1}{2}$ at every step, independently of which of $\{1, 5\}$ it occupies. Thus, for every $n \ge 1$,
   $$
   \mathbb{P}_1(T_1 > n) = \mathbb{P}_1(X_1 = 5, X_2 = 5, \dots, X_n = 5) = \left(\frac{1}{2}\right)^n \longrightarrow 0 \quad (n \to \infty),
   $$
   so $\mathbb{P}_1(T_1 = \infty) = \lim_{n\to\infty} \mathbb{P}_1(T_1 > n) = 0$ and $f_1 = 1$. (Equivalently, $p^{(n)}_{11} = \frac{1}{2}$ for all $n \ge 1$, so $\sum_{n=0}^{\infty} p^{(n)}_{11} = \infty$.) By symmetry between rows/columns $1$ and $5$, $f_5 = 1$ as well. Hence **$1$ and $5$ are recurrent**.
3. **States $2$ and $4$:** Because $\{3\}$ and $\{1, 5\}$ are closed, once the chain enters $\{3, 5\}$ it can never visit $\{2, 4\}$ again. Hence
   $$
   \mathbb{P}_4(T_4 = \infty) \ge \mathbb{P}_4(X_1 \in \{3, 5\}) = p_{43} + p_{45} = \frac{1}{4} + \frac{1}{4} = \frac{1}{2} > 0,
   $$
   $$
   \mathbb{P}_2(T_2 = \infty) \ge \mathbb{P}_2(X_1 = 4, X_2 \in \{3, 5\}) = p_{24}(p_{43} + p_{45}) = \frac{1}{2} \cdot \frac{1}{2} = \frac{1}{4} > 0,
   $$
   which already gives $f_4 \le \frac{1}{2} < 1$ and $f_2 \le \frac{3}{4} < 1$, proving directly that **$2$ and $4$ are transient**.

In fact, we can compute the exact return probabilities and expected numbers of visits on $\{2, 4\}$ cleanly:
- The restriction of $P$ to $\{2, 4\}$ is the rank-1 matrix
  $$
  Q = \begin{pmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{4} & \frac{1}{4} \end{pmatrix},
  $$
  which satisfies $Q^2 = \frac{3}{4} Q$, and therefore by induction $Q^n = \left(\frac{3}{4}\right)^{n-1} Q$ for all $n \ge 1$.
- Since $\{1, 3, 5\}$ cannot lead back to $\{2, 4\}$, the $\{2, 4\}$-block of $P^n$ is precisely $Q^n$ for all $n \ge 0$. Summing over $n \ge 0$ gives the exact expected-visit (fundamental) matrix:
  $$
  \sum_{n=0}^{\infty} Q^n = (I - Q)^{-1} = I + \sum_{n=1}^{\infty} \left(\frac{3}{4}\right)^{n-1} Q = I + 4Q = \begin{pmatrix} 3 & 2 \\ 1 & 2 \end{pmatrix}.
  $$
- Therefore:
  $$
  \sum_{n=0}^{\infty} p^{(n)}_{22} = 3 < \infty \implies f_2 = 1 - \frac{1}{3} = \frac{2}{3} < 1,
  $$
  $$
  \sum_{n=0}^{\infty} p^{(n)}_{44} = 2 < \infty \implies f_4 = 1 - \frac{1}{2} = \frac{1}{2} < 1,
  $$
  confirming once again that states $2$ and $4$ are transient.

## Final Answers

With states labelled $1, 2, 3, 4, 5$ corresponding to the rows/columns of $P$ in Example 8:
- **Recurrent states:** $1, 3, 5$ (forming the finite closed communicating classes $\{1, 5\}$ and $\{3\}$; all three are positive recurrent, with $f_1 = f_3 = f_5 = 1$).
- **Transient states:** $2, 4$ (forming the non-closed communicating class $\{2, 4\}$, with exact return probabilities $f_2 = \frac{2}{3}$ and $f_4 = \frac{1}{2}$).

## Computational Verification

- Script: `private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-1_qx12_check.py`
- Command: `python3 "private/audits/markov-course-2026-11/rework/IB_markov-chains_sheet-1_qx12_check.py"` (exit code 0)
- Verbatim stdout:
```text
PASS: Communicating classes = {1, 5} (closed), {3} (closed), {2, 4} (open).
PASS: Exact rank-1 identity Q^m = (3/4)^{m-1} Q verified for m = 1..30.
PASS: Exact first-return probabilities f_i^{(m)} verified for all states i in {1..5}, m = 1..30.
PASS: Recurrent states = {1, 3, 5} (f_1 = f_3 = f_5 = 1); transient states = {2, 4} (f_2 = 2/3, f_4 = 1/2).
PASS: Fundamental matrix on {2, 4} is (I - Q)^(-1) = [['3', '2'], ['1', '2']].
```

## Confidence & Edge-Case Notes

- **Confidence:** High.
- **Subtleties / Phase B comparison watchpoints:**
  1. **Cross-reference to Example 8:** `question.tex` says "In Example 8" without reproducing the matrix $P$ (which is on page 2 of `2010-11/sheet-1.pdf`), and the locator only lists `pages: [3]`. Check whether the existing solution quotes the matrix $P$ from Question 8 and verifies all entries accurately (especially row 4 being $(0, \frac{1}{4}, \frac{1}{4}, \frac{1}{4}, \frac{1}{4})$, NOT $(\frac{1}{4}, \frac{1}{4}, \frac{1}{4}, \frac{1}{4}, 0)$).
  2. **State indexing:** States are conventionally indexed $\{1, 2, 3, 4, 5\}$ (matching Norris Example 1.2.2 and standard IB sheet conventions); check that the existing solution explicitly states its indexing convention.
  3. **Justification of recurrence/transience:** Check that the solution justifies *both* why $\{1, 5\}$ and $\{3\}$ are recurrent (finite + closed communicating classes, not merely closed — closed infinite classes can be transient!) and why $\{2, 4\}$ is transient (not closed, or direct bound $\mathbb{P}_i(T_i = \infty) > 0$).
