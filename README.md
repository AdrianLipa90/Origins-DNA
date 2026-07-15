# Origins-DNA

A reproducible research scaffold for testing information-theoretic representations of DNA.

## Status

Version `0.1.0` is a **software baseline**, not a claim that DNA has been decrypted. It provides:

- strict DNA sequence normalization and validation;
- reversible 2-bit encoding (`A=00`, `C=01`, `G=10`, `T=11`);
- deterministic substitution-mutation simulation;
- basic sequence statistics and k-mer counts;
- tests that define the first stable behaviour of the project;
- an append-only change log.

The intended research direction is to compare conventional sequence representations against
NOEMA-inspired, orbital and atomic discretizations under explicit null models, held-out data,
ablations and reproducible metrics.

## Quick start

```bash
python -m pip install -e ".[test]"
pytest
```

## Minimal example

```python
from origins_dna import decode_2bit, encode_2bit

encoded = encode_2bit("ACGTAC")
assert decode_2bit(encoded) == "ACGTAC"
```

## Research discipline

1. Biological and software claims must be separated.
2. Every source-code change requires a new append-only entry in `logs/changes.jsonl`.
3. Existing log entries must never be deleted, reordered or rewritten.
4. Hypotheses must be tested against conventional baselines and null models.
5. A passing software test is not evidence of a biological mechanism.

See `docs/RESEARCH_PROTOCOL.md` for the first validation plan.
