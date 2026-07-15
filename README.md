# Origins-DNA

A reproducible research scaffold for testing information-theoretic representations of DNA.

## Status

Version `0.6.0` is a **projection-seed sensitivity diagnostic with a null-bias block**, not a claim that DNA has been decrypted. It provides:

- strict DNA normalization, reversible 2-bit encoding and deterministic mutation simulation;
- an experimental orbital/atomic map with explicit complement symmetries;
- frozen v0.3-v0.5 benchmarks with negative, mixed and partially supportive results retained;
- local orbital tensors and conventional k-mer controls;
- equal-output-dimension comparisons at 12, 24 and 48 dimensions;
- an explicitly disclosed 16-seed exploratory pilot;
- an independent 32-projection-seed confirmation with hierarchical bootstrap intervals;
- observed-label, label-permutation and sequence-shuffle conditions;
- append-only change logs, append-only legacy notes and clean LaTeX reports.

## Current controlled result

The paired contrast is `projected tetranucleotide - projected orbital` balanced accuracy. Positive values favour projected 4-mer features.

### Independent observed-label confirmation

| Output dimension | Mean difference | Hierarchical 95% CI | Positive projection seeds | Preregistered observed gate |
|---:|---:|---:|---:|---|
| 12 | +0.0233 | [0.0046, 0.0425] | 68.75% | unresolved: sign agreement below 75% |
| 24 | +0.0444 | [0.0301, 0.0590] | 84.38% | observed 4-mer advantage |
| 48 | +0.0775 | [0.0645, 0.0908] | 100.00% | observed 4-mer advantage |

### Null gate

The substantive promotion gate **fails**:

- label-permutation null excludes zero at 24 and 48 dimensions;
- sequence-shuffle null excludes zero at 12, 24 and 48 dimensions.

Therefore the observed-label advantages at 24 and 48 dimensions are reproducible but are not promoted as clean representation effects. The overall verdict is `WARN_REVIEW_REQUIRED`: a null-calibrated and bias-corrected contrast must be preregistered next.

The 16-seed pilot is retained separately because it was observed before repository preregistration. It is not pooled with the 32-seed confirmation.

## Quick start

```bash
python -m pip install -e ".[test]"
pytest
python scripts/run_controlled_benchmark.py
python scripts/run_local_order_benchmark.py
python scripts/run_dimension_matched_benchmark.py
python scripts/run_projection_seed_sensitivity.py
```

The runners write full JSON, deterministic gzip archives and compact summaries.

## Minimal example

```python
from origins_dna import SyntheticMotifConfig, run_projection_seed_sensitivity

report = run_projection_seed_sensitivity(
    SyntheticMotifConfig(
        sample_count=40,
        sequence_length=24,
        seeds=(3, 5),
        bootstrap_resamples=50,
    ),
    budgets=(12,),
    projection_seeds=(101, 202),
    hierarchical_bootstrap_resamples=50,
)
assert report.benchmark_status == "synthetic_projection_seed_sensitivity_result"
```

## Research discipline

1. Biological, mathematical and software claims remain separate.
2. Every source-code change requires a new append-only entry in `logs/changes.jsonl`.
3. `notes/LEGACY_NOTES.md` is append-only.
4. Existing logs and notes are never deleted, reordered or silently rewritten.
5. Comparisons use conventional controls, null models, ablations and identical held-out indices.
6. Negative, neutral, mixed and blocked results are preserved plainly.
7. A passing test or synthetic benchmark is not evidence of a biological mechanism.
8. Code and LaTeX must satisfy `docs/DEVELOPMENT_POLICY.md`.

See the protocols under `docs/`, including `PROJECTION_SEED_SENSITIVITY_PROTOCOL.md`, and the corresponding sources under `latex/`.
