# Origins-DNA

A reproducible research scaffold for testing information-theoretic representations of DNA.

## Status

Version `0.3.0` is a **controlled software benchmark**, not a claim that DNA has been decrypted. It provides:

- strict DNA sequence normalization and reversible 2-bit encoding;
- deterministic mutation simulation and conventional sequence statistics;
- an experimental orbital/atomic feature map with explicit complement symmetries;
- dependency-free classifiers, stratified splits and deterministic bootstrap intervals;
- observed-label, label-permutation and sequence-shuffle conditions;
- orbital channel ablations against base-frequency, dinucleotide and positional one-hot controls;
- append-only change logs and append-only legacy notes;
- clean, compilable LaTeX notes for the representation and benchmark.

## Current controlled result

On the v0.3 composition-matched synthetic motif task, the conventional dinucleotide-frequency baseline performed best:

| Representation | Mean balanced accuracy | 95% bootstrap CI |
|---|---:|---:|
| Dinucleotide frequency | 0.6354 | [0.6188, 0.6490] |
| Orbital full | 0.5365 | [0.5188, 0.5583] |
| Orbital phase only | 0.5354 | [0.5167, 0.5594] |
| Orbital classes only | 0.5375 | [0.5146, 0.5563] |
| Positional one-hot | 0.5312 | [0.5094, 0.5531] |
| Base frequency | 0.5250 | [0.5010, 0.5479] |

This result does **not** support an advantage for the current fixed-size orbital summary on this task. The result is retained because negative and neutral evidence is part of the project record.

## Quick start

```bash
python -m pip install -e ".[test]"
pytest
python scripts/run_controlled_benchmark.py
```

The runner writes full JSON plus a deterministic gzip archive and compact summary. The repository retains the compressed full result and the summary.

## Minimal examples

```python
from origins_dna import decode_2bit, encode_2bit

encoded = encode_2bit("ACGTAC")
assert decode_2bit(encoded) == "ACGTAC"
```

```python
from origins_dna import SyntheticMotifConfig, run_controlled_benchmark

report = run_controlled_benchmark(
    SyntheticMotifConfig(sample_count=40, sequence_length=24, seeds=(3, 5, 7))
)
assert report.benchmark_status == "synthetic_controlled_software_result"
```

## Research discipline

1. Biological, mathematical and software claims must remain separate.
2. Every source-code change requires a new append-only entry in `logs/changes.jsonl`.
3. `notes/LEGACY_NOTES.md` is also append-only.
4. Existing logs and notes must never be deleted, reordered or silently rewritten.
5. Hypotheses must be tested against conventional baselines, null models and ablations.
6. Negative and neutral results must be preserved and reported plainly.
7. A passing software test or synthetic benchmark is not evidence of a biological mechanism.
8. Code and LaTeX must satisfy the gates in `docs/DEVELOPMENT_POLICY.md`.

See `docs/BENCHMARK_PROTOCOL.md` for the v0.3 protocol and
`latex/controlled_benchmark_note.tex` for the clean LaTeX report.
