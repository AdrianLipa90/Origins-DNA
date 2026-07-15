# Origins-DNA development policy

## Non-destructive history

- `logs/changes.jsonl` is append-only.
- `notes/LEGACY_NOTES.md` is append-only.
- Existing legacy records are not deleted, reordered, silently corrected, or rewritten.
- Corrections are new entries that identify the earlier record and explain the resolution.

## Clean code gate

Every source change should have:

1. a narrow public API and explicit type contracts;
2. deterministic behaviour unless randomness is an explicit, seeded input;
3. tests for invariants, validation errors, edge cases and declared symmetries;
4. no biological conclusion inferred from software correctness;
5. no hidden destructive repository operation.

## Controlled benchmark gate

Every comparative benchmark should:

1. pre-declare the task, representations, classifier, splits and metrics;
2. use identical held-out indices for compared representations within a seed;
3. include conventional controls, relevant null models and ablations;
4. preserve all results, including negative and neutral outcomes;
5. write a machine-readable result with reproducible configuration and checksum;
6. distinguish synthetic software performance from biological evidence.

## Clean LaTeX gate

Every LaTeX source should:

1. compile without errors and avoid avoidable warnings;
2. keep notation defined before use;
3. separate definitions, hypotheses, methods, results and limitations;
4. use labels and cross-references instead of hard-coded numbering;
5. distinguish theorem, model assumption, simulation result and interpretation;
6. avoid generated prose or equations that are not traceable to code, data or derivation.

## Review boundary

Passing unit tests promotes code only to **software-verified** status. Scientific promotion additionally requires pre-declared tasks, conventional baselines, null models, ablations, held-out evaluation and reproducibility.
