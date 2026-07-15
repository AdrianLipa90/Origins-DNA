# Origins-DNA v0.6 projection-seed sensitivity protocol

## Status and chronology

This document freezes the independent confirmation design before the 32-seed confirmatory run.

A 16-seed exploratory pilot was executed after the design had been stated in chat but before it was committed to the repository. It is therefore **not preregistered** and is excluded from all confirmatory aggregates. The pilot used the original v0.5 seed `50501` plus 15 seeds derived from master seed `60601`.

The exploratory pilot produced the following non-canonical indications for the paired difference `projected 4-mer minus projected orbital`:

| Output dimension | Pilot mean difference | Hierarchical 95% interval | Positive projection means |
|---:|---:|---:|---:|
| 12 | +0.0014 | [-0.0229, +0.0238] | 9/16 |
| 24 | +0.0467 | [+0.0270, +0.0650] | 14/16 |
| 48 | +0.0632 | [+0.0452, +0.0820] | 16/16 |

The pilot motivates a larger independent schedule; it is not pooled with the confirmatory result.

## Frozen confirmatory schedule

- projection master seed: `70701`;
- projection-seed count: `32`;
- schedule: `derived_seed(70701, ordinal)`, skipping duplicates and every exploratory pilot seed;
- output-dimension budgets: `12`, `24`, `48`;
- data seeds: `11, 23, 37, 41, 53, 67, 79, 83`;
- samples per data seed: `400`;
- sequence length: `64`;
- positive motif: `CGCG`;
- negative insert: composition-matched cyclic decoy;
- train fraction: `0.7`;
- classifier: standardized nearest centroid;
- projection family: deterministic label-blind Rademacher matrices;
- family contrast: projected tetranucleotide frequency minus projected orbital multilag tensor.

The v0.5 raw feature maps, matrix construction, splits, classifier and metrics remain frozen.

## Conditions

Every confirmatory projection seed is evaluated under:

1. observed labels;
2. training-label permutation null;
3. within-sequence shuffle null.

No projection seed may be removed because of its result. No output dimension may be selected post hoc as the sole primary result.

## Estimands

For each output dimension and condition:

1. calculate the paired balanced-accuracy difference for every data seed;
2. calculate each projection seed's mean difference across data seeds;
3. calculate the mean and median across projection seeds;
4. calculate the fraction of projection-seed means above and below zero;
5. calculate the minimum and maximum projection-seed mean;
6. calculate a hierarchical 95% bootstrap interval by resampling projection seeds and then data seeds within each selected projection.

The hierarchical bootstrap uses `5000` deterministic resamples. Within-seed v0.5 intervals retain `2000` resamples.

## Interpretation gates

For observed labels at a given dimension:

- **stable 4-mer advantage**: hierarchical interval lies above zero and at least 75% of projection-seed means are positive;
- **stable orbital advantage**: hierarchical interval lies below zero and at least 75% of projection-seed means are negative;
- **projection-sensitive or unresolved**: all other outcomes.

Null conditions should remain centred near zero. A non-zero null result blocks substantive interpretation and requires review.

These gates classify only this synthetic software task. Equal output dimension is not equal information content, equal geometry or equal effective statistical capacity.

## Integrity and promotion boundary

- The exploratory pilot remains visible and separate.
- The 32 confirmatory seeds must be committed before execution.
- Negative, mixed and contradictory outcomes are retained without reinterpretation.
- No synthetic outcome may be promoted as a biological mechanism, general superiority or DNA decryption.
- Biological data work remains blocked until this sensitivity gate is complete and provenance-controlled dataset rules are separately frozen.
