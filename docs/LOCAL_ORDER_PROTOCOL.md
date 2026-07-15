# Origins-DNA v0.4 local-order diagnostic protocol

## Question

Did the v0.3 fixed-size orbital summary fail mainly because global averaging discarded local sequence order?

The diagnostic compares local orbital tensors with explicit k-mer controls on the same deterministic synthetic motif task used in v0.3. It does not test a biological mechanism.

## Frozen task

- binary, balanced synthetic classification;
- 400 sequences per seed;
- sequence length 64;
- positive insertion: `CGCG`;
- negative insertion: composition-matched cyclic decoy `GCGC`;
- train fraction 0.70, stratified within each seed;
- seeds: 11, 23, 37, 41, 53, 67, 79, 83;
- classifier: standardized nearest centroid;
- metrics: balanced accuracy and Matthews correlation coefficient;
- confidence interval: deterministic 95% bootstrap interval over seeds, 2,000 resamples.

## Representations

Conventional controls:

1. dinucleotide frequency, 16 dimensions;
2. trinucleotide frequency, 64 dimensions;
3. tetranucleotide frequency, 256 dimensions.

Orbital diagnostics:

1. global summary from v0.3, 6 dimensions;
2. lag-1 ordered orbital outer-product tensor, 16 dimensions;
3. lags 1--3 ordered orbital outer-product tensor, 48 dimensions;
4. phase-only lags 1--3 ablation, 12 dimensions;
5. class-only lags 1--3 ablation, 12 dimensions.

For state vector \(s_i\in\mathbb{R}^4\), lag \(\ell\), and channels \(a,b\), the local tensor component is

\[
T^{(\ell)}_{ab}=\frac{1}{n-\ell}\sum_{i=1}^{n-\ell}s_{i,a}s_{i+\ell,b}.
\]

The representation uses no label, motif identity or learned task-specific template.

## Null conditions

Every representation is evaluated on identical held-out indices under:

1. observed labels;
2. training-label permutation;
3. independent within-sequence shuffling.

## Declared comparisons

Paired, per-seed balanced-accuracy differences are reported for:

- orbital multilag minus orbital global;
- orbital multilag minus orbital pair;
- orbital multilag minus 2-mer;
- orbital multilag minus 3-mer;
- 4-mer minus orbital multilag;
- orbital multilag minus phase-only;
- orbital multilag minus classes-only.

Separate bootstrap confidence intervals are computed from the paired per-seed differences.

## Interpretation gates

- Recovery above the global orbital summary supports the narrow diagnosis that local-order information mattered.
- Beating 2-mer or 3-mer on this task does not establish general superiority.
- Failure against 4-mer must be retained and reported.
- The 4-mer baseline has greater dimensionality than the 48-dimensional orbital multilag tensor; this benchmark is task-aligned, not capacity-matched.
- Null performance near chance is required for a minimally credible observed-label result.
- No synthetic result is biological evidence or DNA decryption.
