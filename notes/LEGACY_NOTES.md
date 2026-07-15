# Origins-DNA legacy notes — append-only

Existing entries in this file must never be edited or removed. Corrections are appended as later entries.

## NOTE-0001 — 2026-07-15 — process invariant

Adrian Lipa requires the project to preserve both legacy change logs and legacy notes while maintaining clean code and clean LaTeX. The operational interpretation is:

- change provenance remains in `logs/changes.jsonl`;
- contextual, conceptual and methodological notes remain in this append-only file;
- source code is refactored forward without rewriting historical records;
- LaTeX is maintained as reviewable source with explicit notation, provenance and epistemic boundaries.

## NOTE-0002 — 2026-07-15 — first controlled comparative result

The v0.3 composition-matched synthetic motif benchmark did not show an advantage for the current fixed-size orbital/atomic summary. The dinucleotide-frequency baseline reached mean balanced accuracy 0.6354, while the full orbital summary reached 0.5365. This negative/neutral result is retained as a constraint on future development rather than hidden or reinterpreted as validation.

The result applies only to the declared synthetic task, classifier and summary map. It does not disprove richer local orbital representations, but it blocks any claim that the present summary is already superior to conventional sequence-order features.

## NOTE-0003 — 2026-07-15 — pre-publication entrypoint correction

During the final v0.3 validation, direct execution of `python scripts/run_controlled_benchmark.py` failed when the package had not been installed in editable mode. The branch still pointed to v0.2 and no v0.3 commit had been published. The script was corrected to resolve the repository `src/` directory explicitly, and a subprocess regression test now verifies this exact invocation.

This note preserves the failure as process evidence rather than presenting the final passing state as if the defect had never existed.
