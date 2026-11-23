# Whole-Course Cross-Model Audit Summary: IB Markov Chains

- **Auditor:** Gemini 4 via Antigravity (high-effort reasoning), cross-model to the Codex (GPT) solver
- **Date:** 2026-11-22 | **Seed:** `20261122` | **Sample:** 29 / 141 units (20.6%: 11 / 51 sheet, 18 / 90 Tripos)
- **Blindness:** Clean — all 29 reworks and 29 exact-rational Python check scripts were written and verified passing before opening any stored solution, audit report, or unit metadata.

## Error Rates & Verdict Breakdown

| Metric | Overall (`n = 29`) | Example Sheets (`n = 11`) | Tripos (`n = 18`) |
|---|---:|---:|---:|
| **Substantive error rate** | **0 / 29 (0.0%)** | 0 / 11 (0.0%) | 0 / 18 (0.0%) |
| **Transcription error rate** | **0 / 29 (0.0%)** | 0 / 11 (0.0%) | 0 / 18 (0.0%) |
| **Minor issue rate** | **4 / 29 (13.8%)** | 3 / 11 (27.3%) | 1 / 18 (5.6%) |
| **Unsure** | **0 / 29 (0.0%)** | 0 / 11 (0.0%) | 0 / 18 (0.0%) |
| **Agree (exact results & valid argument)** | **25 / 29 (86.2%)** | 8 / 11 (72.7%) | 17 / 18 (94.4%) |

- **Full-course review triggered?** **No** (`0/29` substantive errors). All 29 stored answers, classifications, edge cases (including subtle boundary cases in `sheet-1/q04`, `sheet-2/q14`, `sheet-2/qx05`, `sheet-2/qx06`, `2013/P1/q20H`, `2014/P4/q9H`, `2019/P1/q20H`), verbatim question texts, and both stored question diagrams (`random-walk-tree.png`, `unopoly-board.png`) are mathematically sound and faithful to the archived PDFs.

## Questions Needing Attention (Ordered by Severity)

### 1. Minor Justification / Completeness Issues (4 units; full quoted fixes in `report.yaml`)

1. **`IB/markov-chains/sheet-2/q15`** *(minor — unchecked uniqueness hypothesis)*: Constructs a reversible distribution $\pi$ in detailed balance and concludes reversibility "in this equilibrium", but **never checks irreducibility** (from $0 < \lambda_i, \mu_{i+1}$) to prove that $\pi$ is the **unique** invariant distribution.
2. **`IB/2016/P3/q9H`** *(minor — missing stopping-time check)*: Invokes the strong Markov property on "the $k$th visit when it occurs" without defining the $k$th visit time $T_i^{(k)}$, verifying it is a stopping time (with $X_{T_i^{(k)}} = i$ on $\{T_i^{(k)} < \infty\}$), or writing out the induction and $r = 0$ conditioning case.
3. **`IB/markov-chains/sheet-2/q11`** *(minor — missing required theorem statement)*: The question explicitly instructs *"quoting carefully any general theorems that you use"*, yet the solution never quotes the Convergence to Equilibrium Theorem (alluding to it in half a sentence), omits the class property of periodicity when passing from residue $0$ to the chain, and asserts uniqueness of the uniform distribution without citing irreducibility.
4. **`IB/markov-chains/sheet-2/q07`** *(minor — over-compressed proof sketch)*: Skips stating/justifying that distance $D_n$ is a Markov chain (lumpability), quotes the gambler's ruin formula $(1 - 1/2)/(1 - 2^{-M})$ without writing down or solving the recurrence $u_k = \frac{2}{3}u_{k+1} + \frac{1}{3}u_{k-1}$, and asserts almost-sure finite-interval exit without justification.

### 2. High-Priority Style & Schedule Observations (for a later style pass)

- **Off-schedule Part II methods** *(against `private/raw/schedules/schedules-2025-26.pdf`)*:
  - **`IB/2021/P1/q19H`** *(highest style priority)*: Proves part (a) using $\sigma$-algebra filtrations, an exponential martingale $r^{S_j}$, Optional Stopping, Strong Law of Large Numbers, and Dominated Convergence — all Part II Probability & Measure — rather than the native IB minimal-non-negative-solution proof.
  - **`IB/markov-chains/sheet-1/q03`** (and lightly **`IB/2014/P4/q9H`**): Uses filtrations and the tower property of conditional expectation instead of elementary IB finite-dimensional summation over intermediate states.
  - *(Note)*: `sheet-2/qx05`, `sheet-2/qx06`, `2013/P3/q9H` use the ergodic theorem (starred on the 2025–26 schedule, though prompted by those 2012–13 questions), and `2016/P3/q9H`, `2019/P4/q9H`, `2021/P2/q18H` cite Monotone Convergence by name for non-negative sums.
- **Metadata / self-containedness note (`IB/markov-chains/sheet-1/qx12`)**: `question.tex` verbatim matches Question 13 on p. 3 (*"In Example 8, which states are recurrent and which are transient?"*), so transcription is `ok`, but the $5\times 5$ matrix $P$ sits on p. 2 of `2010-11/sheet-1.pdf` (locator lists only `[3]`) and is not reproduced in `solution.tex`.
- **Section II Tripos compression**: Seven 20-mark Tripos solutions (`2008/P2/q20H`, `2009/P2/q20H`, `2011/P2/q20H`, `2012/P2/q20H`, `2013/P1/q20H`, `2019/P1/q20H`, `2021/P1/q19H`) compress full derivations into dense prose paragraphs without key intermediate displayed equations.
- **Notation & `templates/STYLE.md` drift**: Missing part labels (`sheet-2/qx01`, `2011/P2/q20H`, `2024/P3/q8H`); `\textbf{...}` headings/labels (`sheet-1/q04`, `sheet-2/q14`, `2014/P4/q9H`); $h_i$ used for expected hitting *times* instead of $k_i$ (`sheet-2/qx01`, `2013/P1/q20H`); `\Pr` vs `\mathbb{P}` (`sheet-2/q11`, `2021/P1/q19H`); non-Norris terms (`kernel`, `stationary mass`, `cut fluxes` in `sheet-2/q14`, `2009/P2/q20H`); American spellings (`normalization`/`normalizable`/`realizes` in `sheet-2/q14`, `sheet-2/qx06`, `2012/P2/q20H`, `2021/P2/q18H`); source meta-commentary (`sheet-1/q03`).
