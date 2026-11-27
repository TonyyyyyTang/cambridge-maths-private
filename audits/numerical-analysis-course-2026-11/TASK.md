# Task: whole-course audit of IB Numerical Analysis (Antigravity)

You are a fresh, independent auditor. The solutions in this repository were written and double-checked by Codex (GPT models); you run on a different model family, so this is a **cross-model audit**. Use the strongest model and the highest reasoning effort available.

## Before you start
1. Read `AGENTS.md` at the repository root and this file. Where they differ, this file decides the scope of your task.
2. **Write only inside `private/audits/numerical-analysis-course-2026-11/`.** Do not modify, move or delete any existing file anywhere else. Do not create branches, commit or push: the main agent (Claude Code) reviews your output and commits it. Other agents may be editing other parts of the repository at the same time; ignore their changes.
3. `private/` is a separate private git repository. Nothing from it may be published.

## The sample (already drawn; do not redraw)
`sample.yaml` lists 24 question units: 20% (rounded up) of each stratum of the course's 119 ai-draft units (24 of 119 tripos units), drawn with seed 20261208. Each entry gives a blind packet in `packets/` (exact question text, source locators and the question figures you may open) and the path where your rework goes.

## Phase A: blind rework (all 24 before any comparison)
Until **every** rework file is written, do not open anything in the units' own folders except the question figures listed in the packet. In particular, do not open `solution.tex`, `independent-solution.md`, `comparison.md`, `check.py`, `meta.yaml`, `transcription-review.md`, `audit.md`, `variants/` or `history/`, and do not open any `private/*audit*` file or other audit folder. Blindness depends on your discipline: if you open anything by accident, record exactly what you saw in that question's rework file.

For each packet:
1. **Transcription check.** Run `private/storage.sh fetch` once. Open the archived PDF at the given page(s) (render with e.g. `gs -q -dNOPAUSE -dBATCH -sDEVICE=png16m -r150 -dFirstPage=P -dLastPage=P -sOutputFile=page.png FILE.pdf`, keeping images outside git, e.g. in /tmp). Compare the packet's question text word for word with the printed page: every symbol, condition, constant, subpart and mark. Record `transcription: ok` or each discrepancy.
2. **Solve it yourself** from the statement alone, at the level of a full Tripos answer. Bookwork must be proved, not cited as "standard".
3. **Check computable results** with SymPy, exact arithmetic or careful numerics. Save the script as `rework/<slug>_check.py`, run it, and paste its output into the rework.
4. Write `rework/<slug>.md` with the full solution, the final answers stated clearly, the check output and your confidence.

This is long. You may work across several sessions: a rework file that already exists is finished, so skip it. Phase A rules still apply until all 24 exist.

## Phase B: comparison
Only once all 24 reworks exist, compare each against the stored `solution.tex` in its unit folder. Give one verdict per question:
- `agree`: same results, valid argument.
- `minor`: results correct, but a step lacks justification, a hypothesis is unchecked, or the presentation is unclear or wrong in a way a supervisor would mark down.
- `substantive`: a wrong result, an invalid proof step, a missing case, or a part not answered.
- `transcription`: the stored question text differs from the printed source.
- `unsure`: you cannot decide. Say what would settle it.

For every verdict other than `agree`, quote the exact step, explain the problem, and give a corrected version. **Do not edit the stored solution.** A substantive error triggers a full-course review under AGENTS.md; do not start it, but flag it prominently.

Also note style observations (even for `agree`): methods outside the IB Numerical Analysis schedule (`private/raw/schedules/schedules-2025-26.pdf`), steps a supervisor would expect but the solution skips, and notation or spelling that departs from `templates/STYLE.md`.

## Deliverables (in this folder)
1. `report.yaml`, in the format of `private/audits/markov-course-2026-11/report.yaml` (read it only after Phase A). Include `scope`, `date`, `auditor` (state the exact model, e.g. "Gemini ... via Antigravity, cross-model to the Codex (GPT) solver"), `seed`, `selector`, `population`, `sample_size`, `sample_fraction`, `strata`, `independence` (an honest account of anything seen before the reworks were frozen), `verdict_counts`, `substantive_error_rate`, `minor_error_rate`, `transcription_error_rate`, `unsure_rate`, `full_course_review_triggered`, `results` (one entry per question: `id`, `verdict`, `transcription`, `details`, `suggested_fix`) and `style_observations`.
2. `summary.md`: one page; error rates first, then the questions needing attention by severity.

When finished, reply with the summary and the list of files you created.
