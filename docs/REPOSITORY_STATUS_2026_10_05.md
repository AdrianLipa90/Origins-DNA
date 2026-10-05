# Origins-DNA Repository Status — 2026-10-05

Status: **RECONCILIATION / LOCAL VALIDATION PASS WITH EXPLICIT SCOPE BOUNDARY**

## Source identity

- repository: `AdrianLipa90/Origins-DNA`
- reconciliation base: `main@447359b69c62ba9fab6902d2b059e9c6fcbb3a1e`
- maintenance branch: `maintenance/origins-dna-reconcile-20261005`
- latest merged research implementation inherited from: v0.9.1 line
- GitHub Actions used for this reconciliation: **no**

## Repository debt found

The base main branch carried four conflicting status surfaces:

1. `README.md` advertised v0.6.0;
2. `pyproject.toml` advertised v0.8.0;
3. source modules and merged research history include v0.9/v0.9.1;
4. `docs/CRYPTOGRAPHIC_BOUNDARY_V0_9.md` reused the label v0.9 for an AEAD design document.

The AEAD merge did not add an `origins_dna.crypto` implementation or crypto tests. Therefore the cryptographic document is classified here as **design-only legacy material**, not implemented functionality.

## Scientific state

### v0.7

The independent null-calibrated analysis is present in code and preregistered, but no default v0.7 result is promoted by this maintenance pass. Its execution remains a distinct scientific action.

### v0.8-v0.9.1

The repository contains the exact DHSE bridge, ordered causal trajectory and reciprocal odd/even phase split. The current evidence boundary remains synthetic. Historical negative and mixed results are retained; the 3/4 readout is not promoted as uniquely privileged.

### PhaseNav Natural

A plan exists outside main, but the planned native artifacts are not present on current main. This maintenance pass does not fabricate them.

### Genomic grammar V17

The local Library contains a later V17 candidate package with self-test PASS / mechanics-only status. It is explicitly `CANDIDATE_ONLY`, `biological_evidence=false`, and is not imported by this reconciliation.

## On-Primes cross-check

Reference head:

`AdrianLipa90/On-Primes@372125962a42e1b225337bf31416233da4ca436b`

Local independent checks reproduced the current exact centered-hexagonal lane:

[
H_n=3n^2+3n+1.
]

Checked:

- `H_n ≡ 1 (mod 6)` for every `0 <= n <= 10000`;
- every prime divisor found in that range lies in the `1 mod 6` channel;
- for every prime `5 <= q <= 5000`, the modular obstruction has exactly two roots for `q ≡ 1 mod 6` and zero roots for `q ≡ 5 mod 6`;
- Eisenstein norm identity reproduced.

Rank-two controls also reproduced:

- same-exponent norm: `2371644030976`;
- prime support: `{2, 11, 71, 181}`;
- `q=11`: exponent logs `(1,1)`;
- `q=61`: exponent logs `(9,11)`, genuinely off diagonal;
- the global torsion relation remains satisfied.

The unconditional rank-two global propagation theorem remains **OPEN**. Nothing in Origins-DNA promotes it.

## RDSS algebra cross-check

Reference head:

`AdrianLipa90/Relational-Dynamic-State-Space@d116bc0b60ebadc65d13b0306427857e75ec0f5c`

Reference closure receipt reports:

- `PASS_STRUCTURAL_ALGEBRA_CLOSED`;
- 407 claims;
- 0 derivation gates;
- 1 evidence gate: physical OFHM binding.

Independent local numerical/algebraic checks reproduced:

- (kappa=ln2/(24pi));
- (pikappa=ln2/24);
- twelve maximally mismatched pair potentials sum to (ln2);
- pair forces are antisymmetric and sum to zero;
- the symmetric second difference annihilates affine trajectories;
- Moire potential is invariant under a global phase shift;
- overdamped gradient flow gives (dot Vle0).

These statements are used as structural consistency references only. RDSS physical binding remains an external evidence frontier and is not turned into a DNA biology claim.

## Local DNA validation

Because a byte-identical Git checkout/archive was not available through the current connector path, this maintenance pass does **not** claim a complete upstream `pytest` execution.

A targeted behavioral reconstruction of current main source surfaces passed:

- `compileall`: PASS;
- 2,056 randomized 2-bit encode/decode round trips;
- strict normalization rejection controls;
- deterministic mutation replay;
- complement and reverse-complement involutions;
- finite orbital summaries;
- deterministic synthetic dataset generation/splitting;
- sequence-shuffle composition preservation;
- vector dimensions: base 4, dinucleotide 16, tetranucleotide 256, orbital-multilag 48;
- held-out nearest-centroid smoke evaluation.

Machine-readable details are in `receipts/LOCAL_RECONCILIATION_VALIDATION_20261005.json`.

## Version policy after reconciliation

- package/research implementation line: **0.9.1**;
- historical study labels stay unchanged in their own frozen documents;
- crypto work receives a separate profile namespace instead of reusing research v0.9;
- open cross-repository mathematics stays open;
- biological interpretation is not inferred from mathematical/software PASS states.

## Next gates

1. obtain or materialize a byte-identical full checkout and run the full local test suite;
2. execute the preregistered v0.7 default result only as an explicit scientific action;
3. bridge V17 as a candidate layer under its existing evidence firewall;
4. implement PhaseNav Natural separately;
5. implement AEAD under a separate crypto-profile namespace with dedicated tests.

No merge to main is authorized by this document.
