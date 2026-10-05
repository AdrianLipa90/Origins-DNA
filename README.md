# Origins-DNA

A reproducible research scaffold for testing information-theoretic, orbital and relational representations of DNA.

## Current status

- **research lineage:** v0.9.1
- **package version:** 0.9.1
- **current reconciliation base:** `main@447359b69c62ba9fab6902d2b059e9c6fcbb3a1e`
- **evidence class:** mathematical/software validation plus controlled synthetic experiments
- **biological mechanism:** not established
- **DNA decryption claim:** not established

The repository has accumulated several research layers. Version numbers in older files describe the historical stage at which those files were written; they are not all current package-version declarations.

## Research chronology

| Stage | Main role | Current interpretation |
|---|---|---|
| v0.1 | reversible normalization / 2-bit encoding | software baseline |
| v0.2 | orbital/atomic nucleotide map | exact representation symmetries; experimental map |
| v0.3 | first controlled benchmark | global orbital summary underperformed a dinucleotide baseline |
| v0.4 | local-order tensors | local order recovered signal; task-aligned 4-mers remained stronger |
| v0.5 | dimension-matched comparison | budget-dependent result |
| v0.6 | projection-seed sensitivity | reproducible differences, but null bias blocked promotion |
| v0.7 | independent null-calibrated design | implementation and preregistration merged; default result remains pending |
| v0.8 | exact DHSE projective bridge | technical exactness; composition-first predictive result near chance |
| v0.9 | ordered noncommutative DHSE trajectory | trajectory retained more order information; primary causal gate did not pass |
| v0.9.1 | reciprocal-odd phase / reciprocal-even split | partial synthetic signal; the declared 3/4 aggregation did not beat ordinary averaging |

Negative, mixed and blocked outcomes remain part of the record.

## Version axes

Research-study numbering and non-research subsystem numbering are now kept separate.

- `v0.x` refers to the Origins-DNA research lineage.
- The legacy file `docs/CRYPTOGRAPHIC_BOUNDARY_V0_9.md` is a **design-only historical crypto document**. Its filename collides with the research v0.9 label and must not be read as an implemented Origins-DNA v0.9 cryptographic release.
- Any future authenticated-encryption implementation should use its own profile namespace (for example `AEAD-PROFILE-0.1`) and dedicated implementation/tests.

At this reconciliation point there is **no `origins_dna.crypto` implementation on main**.

## v0.7 scientific gate

The independent null-calibrated analysis is implemented and preregistered, but the default v0.7 result is still a distinct pending scientific execution. The existing v0.6 conclusion remains:

`WARN_REVIEW_REQUIRED`

The result must not be silently promoted by later code or documentation work.

## Cross-domain mathematics

Two current repositories were checked as mathematical references during the 2026-10-05 reconciliation:

- **On-Primes** supplies exact modular/phase constraints, including the centered-hexagonal channel and a rank-two two-target frontier. Its unconditional rank-two global propagation lemma remains **open**.
- **Relational-Dynamic-State-Space (RDSS)** supplies the current structural relational algebra. Its structural algebra is closed at the recorded 407-claim ledger, while physical binding remains a separate external evidence frontier.

These are consistency/reference layers. They are **not biological evidence for DNA** and are not imported as such.

The placeholder repository `relational-loop-algebra` currently contains only its initial README and is not used as the algebraic source of truth here.

## Local validation policy

GitHub Actions is not treated as available validation authority for this reconciliation. The 2026-10-05 maintenance branch records local checks instead.

The recorded local pass includes:

- source `compileall`;
- 2,056 randomized reversible 2-bit encoding round-trips;
- complement and reverse-complement involution checks;
- deterministic mutation checks;
- shuffle composition preservation;
- representation-dimension checks for 4D base, 16D dinucleotide, 256D tetranucleotide and 48D orbital-multilag features;
- independent checks of the centered-hexagonal prime channel and rank-two controls;
- independent checks of the RDSS κ normalization, pair-force antisymmetry, frame second difference, gauge invariance and Lyapunov sign.

The DNA source check is explicitly a **targeted local behavioral reconstruction**, not a byte-identical full-repository checkout and not a claim of a complete upstream `pytest` run.

See `docs/REPOSITORY_STATUS_2026_10_05.md` and `receipts/LOCAL_RECONCILIATION_VALIDATION_20261005.json`.

## Quick start

```bash
python -m pip install -e ".[test]"
python -m pytest
```

Benchmark entry points include:

```bash
python scripts/run_controlled_benchmark.py
python scripts/run_local_order_benchmark.py
python scripts/run_dimension_matched_benchmark.py
python scripts/run_projection_seed_sensitivity.py
```

The v0.7 default null-calibrated execution is intentionally not listed as an ordinary smoke command because it is a preregistered scientific result-producing action.

## Research discipline

1. Biological, mathematical, cryptographic and software claims remain separate.
2. Every source-code or research-status change receives a new append-only entry in `logs/changes.jsonl`.
3. `notes/LEGACY_NOTES.md` is append-only.
4. Existing logs and notes are never deleted, reordered or silently rewritten.
5. Comparisons use conventional controls, null models, ablations and identical held-out indices where declared.
6. Negative, neutral, mixed and blocked results are preserved plainly.
7. A passing test or synthetic benchmark is not evidence of a biological mechanism.
8. Open results imported from another repository remain open here.
9. Historical branch notes may be marked superseded, but their original content is retained.

See the protocols under `docs/` and the corresponding sources under `latex/`.
