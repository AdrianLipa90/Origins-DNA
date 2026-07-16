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

## NOTE-0009 — 2026-07-15 — v0.5 dimension-matched result

The preregistered result is budget-dependent. At 12 output dimensions, the orbital projection reached mean balanced accuracy 0.6646 while the projected tetranucleotide representation reached 0.5479. The paired tetranucleotide-minus-orbital difference was -0.1167 with a 95% interval of [-0.1552, -0.0823], favouring the orbital source under severe compression.

At 24 dimensions, the projected tetranucleotide representation reached 0.7198 and the orbital projection reached 0.6677; the paired difference was +0.0521 [0.0042, 0.1000]. At 48 dimensions, tetranucleotide reached 0.7115 and orbital reached 0.6542; the paired difference was +0.0573 [0.0260, 0.0885]. The raw anchors remain 0.8521 for 256-dimensional tetranucleotide and 0.6771 for the 48-dimensional orbital tensor.

No global winner is inferred. The pattern shows representation-specific information loss under one preregistered projection seed. Equal output dimension remains only a capacity proxy, and projection-seed sensitivity has not yet been tested. All label-permutation and sequence-shuffle aggregates remain near chance. Two independent executions produced byte-identical result bundles.

Validation was performed in an isolated checkout reconstructed from the GitHub head: 66 tests passed, `compileall` passed, the direct script passed, and the LaTeX report compiled twice without warnings or layout defects. The new workflow did not produce an observable GitHub Actions run because it was not present in the stacked base branch; no CI PASS is claimed.

## NOTE-0010 — 2026-07-15 — exploratory projection-seed pilot and independent confirmation freeze

A 16-projection-seed sensitivity pilot was executed after its outline had been stated in chat but before the schedule was committed to the repository. It is therefore explicitly exploratory, not preregistered, and is excluded from the v0.6 confirmatory aggregates. The pilot used the original v0.5 seed `50501` and 15 seeds derived from master seed `60601`.

The exploratory paired difference `projected 4-mer minus projected orbital` was approximately +0.0014 at 12 dimensions with hierarchical interval [-0.0229, +0.0238], +0.0467 at 24 dimensions with [+0.0270, +0.0650], and +0.0632 at 48 dimensions with [+0.0452, +0.0820]. These observations motivated a larger independent schedule but are not pooled with it.

The confirmatory schedule is frozen before execution: 32 unique seeds derived from master seed `70701`, with all 16 pilot seeds excluded; budgets 12, 24 and 48; frozen v0.5 data seeds, task, splits, classifier, controls and metrics; and a hierarchical bootstrap that resamples projection seeds and data seeds. The full design and interpretation gates are recorded in `docs/PROJECTION_SEED_SENSITIVITY_PROTOCOL.md`. No confirmatory result had been executed when this note was appended.

## NOTE-0011 — 2026-07-15 — independent confirmation result and null-bias block

The independent 32-projection-seed confirmation found a positive observed-label contrast `projected 4-mer minus projected orbital` at every declared dimension. At 12 dimensions the mean was +0.0233 with hierarchical interval [0.0046, 0.0425], but only 68.75% of projection-seed means were positive, below the preregistered 75% sign-stability gate. The 12-dimensional result is therefore projection-sensitive or unresolved.

At 24 dimensions the mean was +0.0444 [0.0301, 0.0590] with 84.38% positive projection seeds. At 48 dimensions the mean was +0.0775 [0.0645, 0.0908] with 100% positive projection seeds. These two observed-label comparisons pass the preregistered observed gate for a 4-mer advantage.

The substantive promotion gate nevertheless fails. The label-permutation null excludes zero at 24 and 48 dimensions, while the sequence-shuffle null excludes zero at 12, 24 and 48 dimensions. The observed advantages are reproducible software results, but they cannot be promoted as clean representation effects while representation-family differences remain present under null conditions. The overall verdict is `WARN_REVIEW_REQUIRED`, with `promotion_allowed=false`.

The next gate is a separately preregistered exchangeable, null-calibrated and bias-corrected observed-minus-null contrast. Potential explanations such as finite-sample geometry, rejection-sampling effects or residual composition after shuffling remain hypotheses only. Biological dataset work and any DNA-decryption claim remain blocked.

## NOTE-0012 — 2026-07-15 — direct validation for main

The post-merge audit found that `main` contained the same file tree as the final validated v0.6 head, but the merge commit itself had no GitHub Actions run. The existing workflow triggered pushes only on `research/origins-dna-v0.5-dimension-matched`, so later merges into `main` were not independently checked.

The workflow is generalized forward-only: pushes to `main`, `research/**` and `agent/**` are validated; pull-request path coverage includes source, tests, scripts, LaTeX, protocols, results, README and package metadata; and every LaTeX report is compiled twice and checked for warnings or layout defects. The historical merge graph is retained unchanged rather than rewritten.


## NOTE-0013 — 2026-07-16 — accidental staging merge cleanup

While work on v0.7 was beginning, PR #18 merged a temporary staging branch into `main`. That merge added `.github/workflows/bootstrap-ci-main-validation.yml` and replaced the validated read-only validation workflow with a staging variant using `contents: write` and an explicit head-ref checkout. PRs #19–#22 then merged `main` and helper branches in both directions, making the graph noisier.

A direct tree comparison against the clean PR-17 merge found no changes to source code, tests, scientific results, package version, `logs/changes.jsonl` or prior legacy notes. The contamination is limited to the two workflow files. The correction removes the bootstrap workflow and restores the exact validated workflow forward-only. Existing merge commits remain immutable process evidence; no force-push or history rewrite is used.
