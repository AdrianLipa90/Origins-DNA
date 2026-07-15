# Origins-DNA

A reproducible research scaffold for testing information-theoretic representations of DNA.

## Status

Version `0.2.0` is a **software and representation baseline**, not a claim that DNA has been decrypted. It provides:

- strict DNA sequence normalization and validation;
- reversible 2-bit encoding (`A=00`, `C=01`, `G=10`, `T=11`);
- deterministic substitution-mutation simulation;
- basic sequence statistics and k-mer counts;
- an experimental orbital/atomic feature map with explicit complement symmetries;
- deterministic local window consolidation and fixed-size summaries;
- tests that define stable software behaviour;
- append-only change logs and append-only legacy notes;
- a clean, compilable LaTeX technical note for the v0.2 representation.

The intended research direction is to compare conventional sequence representations against
NOEMA-inspired, orbital and atomic discretizations under explicit null models, held-out data,
ablations and reproducible metrics.

## Quick start

```bash
python -m pip install -e ".[test]"
pytest
```

## Minimal examples

```python
from origins_dna import decode_2bit, encode_2bit

encoded = encode_2bit("ACGTAC")
assert decode_2bit(encoded) == "ACGTAC"
```

```python
from origins_dna import orbital_atomic_encode, orbital_atomic_summary

states = orbital_atomic_encode("ACGT")
summary = orbital_atomic_summary("ACGT")
assert len(states) == 4
assert summary["phase_coherence"] == 0.0
```

## Research discipline

1. Biological and software claims must be separated.
2. Every source-code change requires a new append-only entry in `logs/changes.jsonl`.
3. `notes/LEGACY_NOTES.md` is also append-only.
4. Existing logs and notes must never be deleted, reordered or silently rewritten.
5. Hypotheses must be tested against conventional baselines and null models.
6. A passing software test is not evidence of a biological mechanism.
7. Code and LaTeX must satisfy the gates in `docs/DEVELOPMENT_POLICY.md`.

See `docs/RESEARCH_PROTOCOL.md` for the validation plan and
`latex/orbital_atomic_note.tex` for the mathematical definition of the v0.2 map.
