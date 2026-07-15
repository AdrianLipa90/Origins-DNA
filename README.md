# Origins-DNA

A reproducible research scaffold for testing information-theoretic representations of DNA.

## Status

Version `0.4.0` is a **local-order diagnostic benchmark**, not a claim that DNA has been decrypted. It provides:

- strict DNA normalization, reversible 2-bit encoding and deterministic mutation simulation;
- an experimental orbital/atomic feature map with explicit complement symmetries;
- frozen v0.3 global-summary benchmarks and their negative result;
- directional orbital lag tensors that retain local order without using labels or motif templates;
- explicit 2-mer, 3-mer and 4-mer controls plus phase-only and class-only ablations;
- observed-label, label-permutation and sequence-shuffle conditions;
- paired per-seed differences with deterministic bootstrap confidence intervals;
- append-only change logs, append-only legacy notes and clean LaTeX reports.

## Current controlled result

The v0.4 diagnostic confirms that global averaging discarded substantial local information, while also preserving the stronger conventional result:

| Representation | Dimension | Mean balanced accuracy | 95% bootstrap CI |
|---|---:|---:|---:|
| Tetranucleotide frequency | 256 | 0.8521 | [0.8333, 0.8688] |
| Orbital multilag tensor | 48 | 0.6771 | [0.6562, 0.7021] |
| Trinucleotide frequency | 64 | 0.6469 | [0.6250, 0.6646] |
| Orbital multilag classes only | 12 | 0.6448 | [0.6167, 0.6729] |
| Dinucleotide frequency | 16 | 0.6354 | [0.6198, 0.6500] |
| Orbital pair tensor | 16 | 0.6219 | [0.6031, 0.6427] |
| Orbital multilag phase only | 12 | 0.5656 | [0.5406, 0.5854] |
| Orbital global summary | 6 | 0.5365 | [0.5198, 0.5573] |

Paired over identical seeds and held-out splits:

- multilag orbital minus global orbital: `+0.1406`, CI `[0.1021, 0.1771]`;
- multilag orbital minus dinucleotide: `+0.0417`, CI `[0.0188, 0.0667]`;
- multilag orbital minus trinucleotide: `+0.0302`, CI `[-0.0010, 0.0646]`;
- tetranucleotide minus multilag orbital: `+0.1750`, CI `[0.1479, 0.1969]`.

Thus local-order retention materially improves the orbital map, but the task-aligned 4-mer remains decisively stronger. The dimensionalities are not matched, so v0.4 is a diagnostic rather than a capacity-controlled superiority test.

## Quick start

```bash
python -m pip install -e ".[test]"
pytest
python scripts/run_controlled_benchmark.py
python scripts/run_local_order_benchmark.py
```

The runners write full JSON, deterministic gzip archives and compact summaries. The v0.4 runner also writes paired-comparison JSON.

## Minimal example

```python
from origins_dna import SyntheticMotifConfig, run_local_order_benchmark

report = run_local_order_benchmark(
    SyntheticMotifConfig(sample_count=40, sequence_length=24, seeds=(3, 5))
)
assert report.benchmark_status == "synthetic_local_order_diagnostic_result"
```

## Research discipline

1. Biological, mathematical and software claims remain separate.
2. Every source-code change requires a new append-only entry in `logs/changes.jsonl`.
3. `notes/LEGACY_NOTES.md` is append-only.
4. Existing logs and notes are never deleted, reordered or silently rewritten.
5. Comparisons use conventional controls, null models, ablations and identical held-out indices.
6. Negative, neutral and partially supportive results are preserved plainly.
7. A passing test or synthetic benchmark is not evidence of a biological mechanism.
8. Code and LaTeX must satisfy `docs/DEVELOPMENT_POLICY.md`.

See `docs/BENCHMARK_PROTOCOL.md` for v0.3, `docs/LOCAL_ORDER_PROTOCOL.md` for v0.4, and the corresponding sources under `latex/`.
