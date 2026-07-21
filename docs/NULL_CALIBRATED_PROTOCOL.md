# Origins-DNA v0.7 — independent null-calibrated contrast protocol

## Status

This document freezes the v0.7 analysis before execution of its default result.
The aggregate v0.6 outcome is already known; therefore v0.7 is not presented as
a blind reinterpretation of the same seeds. It uses projection and data seed
schedules that are disjoint from every v0.5/v0.6 exploratory or confirmatory
schedule.

The work remains a synthetic software experiment. It cannot establish a
biological mechanism or DNA decryption.

## Frozen task

The following v0.6 elements remain unchanged:

- balanced composition-matched synthetic motif task;
- sample count 400 and sequence length 64;
- motif `CGCG`;
- train fraction 0.7;
- nearest-centroid classifier;
- balanced accuracy;
- output-dimension budgets 12, 24 and 48;
- deterministic label-blind Rademacher projections;
- label-permutation and sequence-shuffle controls.

## Independent schedules

Projection seeds:

- master seed: `80801`;
- count: 32;
- all 16 exploratory pilot seeds and all 32 v0.6 confirmatory seeds excluded.

Data seeds:

- master seed: `90901`;
- count: 8;
- prior data seeds `11, 23, 37, 41, 53, 67, 79, 83` excluded.

Schedules are derived deterministically with the existing `derived_seed`
operator and are committed before result execution.

## Primary estimands

Let

\[
\Delta_c = \mathrm{BA}_{4mer,c} - \mathrm{BA}_{orbital,c}
\]

for condition \(c\), paired within projection seed, data seed, output dimension
and held-out split.

V0.7 preserves two distinct calibrated estimands:

\[
C_{\mathrm{perm}} = \Delta_{\mathrm{observed}} -
                     \Delta_{\mathrm{label\ permutation}},
\]

\[
C_{\mathrm{shuffle}} = \Delta_{\mathrm{observed}} -
                        \Delta_{\mathrm{sequence\ shuffle}}.
\]

The two nulls are never averaged or selected after observation.

## Uncertainty and exchangeability

For each null and budget:

1. a hierarchical bootstrap resamples projection seeds and then data seeds
   within selected projection seeds;
2. a two-sided cluster sign-flip randomization test flips complete projection
   seed clusters, preserving within-cluster dependence;
3. projection-seed sign stability is reported.

Default gates:

- hierarchical 95% interval excludes zero;
- at least 75% of projection-seed means have the interval direction;
- cluster sign-flip \(p < 0.05\).

## Decision rule

A budget passes only when both null-calibrated estimands pass all three gates
and point in the same direction. The overall synthetic promotion gate passes
only when every declared budget passes.

No pooled-null fallback is allowed. A disagreement between nulls is
`WARN_REVIEW_REQUIRED`, not an invitation to select the favourable null.

## Control diagnostics

Sequence shuffling must preserve the exact A/C/G/T count of every sequence.
The maximum per-sequence residual is recorded and must equal zero.

Negative examples are still generated through motif-avoidance rejection
sampling. V0.7 records this as an unresolved generator limitation rather than
claiming that observed-minus-null subtraction causally removes it. Biological
dataset progression remains blocked regardless of the synthetic decision.

## Reproducibility

The result bundle will contain:

- full per-projection/per-data calibrated records;
- deterministic gzip archive;
- compact summary;
- separate decisions for every budget;
- shuffle-composition diagnostics;
- all seeds, thresholds and resample counts.

Results must be retained whether favourable, neutral or unfavourable.
