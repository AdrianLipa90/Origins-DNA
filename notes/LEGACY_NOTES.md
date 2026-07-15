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

## NOTE-0004 — 2026-07-15 — local-order diagnostic result

The v0.4 benchmark tested the explicit diagnosis that the v0.3 global orbital summary lost useful local order. The 48-dimensional lags 1--3 orbital tensor reached mean balanced accuracy 0.6771, compared with 0.5365 for the global orbital summary. The paired improvement was 0.1406 with a 95% bootstrap interval of [0.1021, 0.1771].

The result is only partially supportive. The same orbital tensor exceeded the 2-mer baseline and was close to the 3-mer baseline, but the task-aligned 4-mer baseline reached 0.8521 and exceeded the orbital tensor by 0.1750 [0.1479, 0.1969]. This failure against 4-mer is retained without reinterpretation. The comparison is not capacity-matched because the 4-mer vector has 256 dimensions and the orbital multilag tensor has 48.

The next scientifically useful gate is a capacity-matched comparison followed by a provenance-controlled biological dataset. No synthetic result is to be promoted as biological mechanism evidence.

## NOTE-0005 — 2026-07-15 — GitHub blob transport correction

Immediate readback of the first v0.4 branch commit showed corrupted bytes inside the GitHub blob for `logs/changes.jsonl` sequence 5, although the validated local source was valid. `main` was not affected and no pull request had been opened. A later corrective commit restores the intended sequence-5 record and appends sequence 6. The faulty commit remains in Git history; it is not rewritten or hidden.

## NOTE-0006 — 2026-07-15 — v0.5 preregistration

V0.5 predeclares a dimension-matched capacity proxy before observing its results. Orbital multilag tensors and tetranucleotide frequencies are projected independently through the same label-blind Rademacher projection family at output dimensions 12, 24 and 48, while the v0.4 task, splits, classifier, controls and metrics remain frozen.

Equal dimension is not treated as proof of equal information content or equal effective model capacity. Raw 48-dimensional orbital and raw 256-dimensional tetranucleotide representations remain diagnostic anchors only. Result interpretation will be appended after CI execution and will be retained whether favourable, neutral or unfavourable.

## NOTE-0007 — 2026-07-15 — pre-result deterministic-bundle test correction

The first isolated full-suite run of the preregistered v0.5 head produced 65 passes and one failed test. The benchmark implementation did not fail. The test incorrectly required summaries written under different output basenames (`first.json` and `second.json`) to be byte-identical even though each summary intentionally records its own archive filename.

Before executing the full v0.5 benchmark, the test was corrected to write the same basename in separate directories and to compare the resulting full JSON, deterministic gzip archive, summary and paired-comparison artifact. The corrected isolated regression passed 66 of 66 tests. No budget, projection seed, task, split, classifier, metric or primary estimand was changed.

## NOTE-0008 — 2026-07-15 — pre-result execution-efficiency correction

The first full v0.5 benchmark attempt exceeded the execution window before producing any artifact. Inspection showed redundant implementation work: every projected representation repeatedly converted and validated the same fixed matrix for every sequence, while the same raw orbital and tetranucleotide vectors were recomputed independently for each dimension budget.

Before any v0.5 result was observed, matrix validation was moved to construction time and source vectors were shared through deterministic in-memory caches. The Rademacher matrices, projection seeds, dimensions, data, splits, classifier, metrics and estimand remain unchanged. A smoke run retained 66 of 66 passing tests and completed normally.
