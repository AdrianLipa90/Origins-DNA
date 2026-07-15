# Origins-DNA

A reproducible research scaffold for testing information-theoretic representations of DNA.

## Status

Version `0.5.0` is a **dimension-matched synthetic benchmark**, not a claim that DNA has been decrypted. It provides:

- strict DNA normalization, reversible 2-bit encoding and deterministic mutation simulation;
- an experimental orbital/atomic map with explicit complement symmetries;
- frozen v0.3 and v0.4 benchmarks with negative and partially supportive results retained;
- local orbital tensors and conventional k-mer controls;
- a preregistered equal-output-dimension capacity proxy at 12, 24 and 48 dimensions;
- deterministic label-blind Rademacher projections, identical held-out splits and paired confidence intervals;
- observed-label, label-permutation and sequence-shuffle conditions;
- append-only change logs, append-only legacy notes and clean LaTeX reports.

## Current controlled result

V0.5 compresses the 48-dimensional orbital multilag tensor and the 256-dimensional tetranucleotide vector independently through the same declared random-projection family. Equal dimension is only a capacity proxy; it does not equalize information content.

| Representation | Output dimension | Mean balanced accuracy | 95% bootstrap CI |
|---|---:|---:|---:|
| Tetranucleotide raw anchor | 256 | 0.8521 | [0.8344, 0.8698] |
| Tetranucleotide projected | 24 | 0.7198 | [0.6865, 0.7531] |
| Tetranucleotide projected | 48 | 0.7115 | [0.6906, 0.7323] |
| Orbital multilag raw anchor | 48 | 0.6771 | [0.6562, 0.7021] |
| Orbital projected | 24 | 0.6677 | [0.6396, 0.6958] |
| Orbital projected | 12 | 0.6646 | [0.6385, 0.6896] |
| Orbital projected | 48 | 0.6542 | [0.6177, 0.6844] |
| Tetranucleotide projected | 12 | 0.5479 | [0.5156, 0.5823] |

The preregistered paired estimand is `projected tetranucleotide - projected orbital`:

- 12 dimensions: `-0.1167`, CI `[-0.1552, -0.0823]`;
- 24 dimensions: `+0.0521`, CI `[0.0042, 0.1000]`;
- 48 dimensions: `+0.0573`, CI `[0.0260, 0.0885]`.

Thus the result is budget-dependent. The orbital representation is stronger after severe compression to 12 dimensions, while projected 4-mer features are modestly stronger at 24 and 48 dimensions. This is not general superiority for either family: the comparison uses one preregistered projection seed, and equal output dimension does not imply equal retained information. All label-permutation and sequence-shuffle aggregates remain near chance.

## Quick start

```bash
python -m pip install -e ".[test]"
pytest
python scripts/run_controlled_benchmark.py
python scripts/run_local_order_benchmark.py
python scripts/run_dimension_matched_benchmark.py
```

The runners write full JSON, deterministic gzip archives and compact summaries. V0.4 and v0.5 also write paired-comparison JSON.

## Minimal example

```python
from origins_dna import SyntheticMotifConfig, run_dimension_matched_benchmark

report = run_dimension_matched_benchmark(
    SyntheticMotifConfig(sample_count=40, sequence_length=24, seeds=(3, 5)),
    budgets=(12, 24),
)
assert report.benchmark_status == "synthetic_dimension_matched_capacity_proxy_result"
```

## Research discipline

1. Biological, mathematical and software claims remain separate.
2. Every source-code change requires a new append-only entry in `logs/changes.jsonl`.
3. `notes/LEGACY_NOTES.md` is append-only.
4. Existing logs and notes are never deleted, reordered or silently rewritten.
5. Comparisons use conventional controls, null models, ablations and identical held-out indices.
6. Negative, neutral, mixed and partially supportive results are preserved plainly.
7. A passing test or synthetic benchmark is not evidence of a biological mechanism.
8. Code and LaTeX must satisfy `docs/DEVELOPMENT_POLICY.md`.

See `docs/BENCHMARK_PROTOCOL.md` for v0.3, `docs/LOCAL_ORDER_PROTOCOL.md` for v0.4, `docs/DIMENSION_MATCHED_PROTOCOL.md` for v0.5, and the corresponding sources under `latex/`.
