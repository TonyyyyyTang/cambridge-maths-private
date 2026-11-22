# Task: whole-course audit of IB Markov Chains (Antigravity)

You are a fresh, independent auditor. The solutions in this repository were written and double-checked by Codex (GPT models). You run on a different model family, so your job is a **cross-model audit**. Use the highest reasoning effort available.

## Before you start
1. Read `AGENTS.md` at the repository root (the project rules) and this file. Where they differ, this file decides the scope of your task.
2. **Write only inside `private/audits/markov-course-2026-11/`.** Do not modify, move or delete any existing file anywhere in the repository. Do not create branches, commit or push: the main agent (Claude Code) reviews your output and commits it.
3. `private/` is a separate private git repository. Nothing from it may be published.

## The sample (already drawn; do not redraw)
`sample.yaml` lists 29 question units. That is 20% (rounded up) of each stratum of the course's 141 units: 11 of 51 example-sheet units and 18 of 90 Tripos units, drawn with seed 20261122. Each entry gives a blind packet in `packets/` (the exact question text and source locators only) and the path where your rework goes.

## Phase A: blind rework (all 29 before any comparison)
Until **every** rework file is written, you must not open anything in the units' own folders: `solution.tex`, `independent-solution.md`, `comparison.md`, `check.py`, `meta.yaml`, `transcription-review.md`, `audit.md`, `history/`. The same goes for any existing `private/*audit*.yaml`. All of these sit in the same repository, so blindness depends on your discipline. If you open one by accident, say exactly what you saw in that question's rework file.

For each packet:
1. **Transcription check.** Open the archived PDF at the page given in the packet (`private/raw/...`; run `private/storage.sh fetch` first) and compare the packet's question text word for word against the printed page: every symbol, condition, constant, subpart and mark. Record `transcription: ok`, or describe each discrepancy.
2. **Solve it yourself**, from the statement alone, at the level of a full Tripos answer. Bookwork must be proved, not cited as "standard".
3. **Check computable results** (hitting probabilities, invariant distributions, expected times, matrix powers, generating functions) with SymPy or exact rational arithmetic. Save the script as `rework/<slug>_check.py`, run it, and paste its output into the rework.
4. Write `rework/<slug>.md`: the full solution, the final answers stated clearly, the check output, and your confidence.

## Phase B: comparison
Only once all 29 reworks exist, compare each against the stored `solution.tex` in its unit folder. Give one verdict per question:
- `agree`: same results, valid argument.
- `minor`: results correct, but a step lacks justification, a hypothesis goes unchecked, or the presentation is unclear or wrong in a way a supervisor would mark down.
- `substantive`: a wrong result, an invalid proof step, a missing case, or a part not answered.
- `transcription`: the stored question text differs from the printed source.
- `unsure`: you cannot decide. Say what would settle it.

For every verdict other than `agree`, quote the exact step, explain the problem, and give a corrected version. **Do not edit the stored solution.** If you find a substantive error, AGENTS.md requires a full-course review. Do not start it: flag it prominently, and the main agent will decide.

Note separately, even for `agree` verdicts, these style observations for a later pass:
- methods outside the IB Markov Chains schedule (see `private/raw/schedules/schedules-2025-26.pdf`);
- steps a Cambridge supervisor would expect but the solution skips;
- notation or naming that differs from the standard course usage.

## Deliverables (in this folder)
1. `report.yaml`, following the format of existing reports such as `private/markov-sheet-2-union-audit.yaml`. Read that file only **after** Phase A is complete. Include at least:
   - `scope`, `auditor` (state the model, e.g. "Claude Opus 5.5 via Antigravity, cross-model to the Codex solver"), `seed`, `population`, `sample_size`, `sample_fraction`, `strata`
   - `independence`: an honest account of anything you saw before freezing your reworks
   - `results`: one entry per question with `id`, `verdict`, `transcription`, `details`, `suggested_fix`
   - counts per verdict, the substantive error rate with its sample size (e.g. 0/29), and the transcription error rate
   - `style_observations`
2. `summary.md`: a one-page human-readable summary. Start with the error rates, then list the questions that need attention, ordered by severity.

When you finish, reply with the summary and the list of files you created.
