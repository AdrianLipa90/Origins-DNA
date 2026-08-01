# Origins-DNA v0.9.1 — preregistration amendment 001

## Timing

This clarification is recorded before the first v0.9.1 workflow or predictive benchmark run.

## Projection matching

The `phase_34` and `phase_mean` source vectors have identical coordinate order and dimension. At each output budget they use the same deterministic Rademacher projection matrix. Their comparison therefore changes only the frozen aggregation exponent:

\[
L^{-3/4}\sum_k x_k
\quad\text{versus}\quad
L^{-1}\sum_k x_k.
\]

The reciprocal-even channel is projected independently. In the split representation, phase and even coordinates are projected in separate blocks and concatenated; no projection row mixes the two channels.

## Primary inference across budgets

To avoid selecting a favourable budget after observation, the primary gate aggregates budgets before inference.

For each seed and representation family, calculate the arithmetic mean balanced accuracy over the fixed output budgets

\[
12,24,48.
\]

Then form paired seed-level differences. A primary contrast passes only when its preregistered 95% bootstrap confidence interval has a lower bound strictly above zero.

The three-quarter phase signal status is `PASS` only if all three primary contrasts pass:

1. `phase_34 observed - even_34 observed`;
2. `phase_34 observed - phase_mean observed`;
3. `phase_34 observed - phase_34 sequence_shuffle_null`.

Budget-specific results remain secondary and all are reported.

## Integrity boundary

This amendment does not change the operators, data labels, motif, decoy, exponent, classifier, seeds, budgets or decision threshold. It only makes the already declared matched-projection and across-budget inference rules explicit before execution.
