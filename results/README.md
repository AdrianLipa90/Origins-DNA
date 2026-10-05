# Current results boundary

> **Reconciled 2026-10-05.** The committed deterministic result bundles in this directory currently cover v0.3-v0.6. Later v0.8-v0.9.1 source/protocol layers exist in the repository, but this directory does not contain equivalent committed machine-readable bundles for those stages. The preregistered v0.7 default result remains pending. Do not infer a missing result artifact from the presence of later source code.

# Results provenance

The repository retains deterministic machine-readable artifacts for each controlled benchmark.

## v0.3 global-summary benchmark

- summary: `v0.3_synthetic_motif_benchmark_summary.json`
- complete deterministic archive: `v0.3_synthetic_motif_benchmark.json.gz`
- schema: `origins-dna.benchmark.v1`
- generator: `scripts/run_controlled_benchmark.py`
- uncompressed JSON SHA-256: `07ed3bc941821785131b61c667ca9a81957fc187d5ac7c33d1875a5f83b6ff21`
- deterministic gzip SHA-256: `01f44457c8932b7d799ef0bc09926a76121017084216d1fc01edd9f9d7300545`
- compact summary SHA-256: `ff169e3fb076b8db8e512f78b8758e239b4339b89efce77c525d6524f864d3ab`

## v0.4 local-order diagnostic

- summary: `v0.4_local_order_benchmark_summary.json`
- paired comparisons: `v0.4_local_order_benchmark_comparisons.json`
- complete deterministic archive: `v0.4_local_order_benchmark.json.gz`
- schema: `origins-dna.local-order-benchmark.v1`
- comparison schema: `origins-dna.local-order-comparisons.v1`
- generator: `scripts/run_local_order_benchmark.py`
- protocol: `docs/LOCAL_ORDER_PROTOCOL.md`
- uncompressed JSON SHA-256: `57affa6729a304d1b5d4d14ba7a7f1e4b111bc59f8303d0324d785d7db12e980`
- deterministic gzip SHA-256: `dea0b159e2254e6a349a45300d553f1452dc40f691cbd2ec345aaefdd0045235`
- compact summary SHA-256: `abc9295bf9f4c07e126f702e9d076cc6104087a68ff942d31cdd72f6731e93d9`
- paired comparisons SHA-256: `9ec6ca1189723efbe2e3f811164ea6b65170dab4dd3b0491d11269cc7c0af6aa`

## v0.5 dimension-matched capacity proxy

- summary: `v0.5_dimension_matched_benchmark_summary.json`
- paired comparisons: `v0.5_dimension_matched_benchmark_comparisons.json`
- complete deterministic archive: `v0.5_dimension_matched_benchmark.json.gz`
- schema: `origins-dna.dimension-matched-benchmark.v1`
- comparison schema: `origins-dna.dimension-matched-comparisons.v1`
- generator: `scripts/run_dimension_matched_benchmark.py`
- protocol: `docs/DIMENSION_MATCHED_PROTOCOL.md`
- projection seed: `50501`
- matched dimensions: `12, 24, 48`
- uncompressed JSON SHA-256: `d62c6c130fc4186ef72d073a955efa3afb386f37ed98811bb75c8ef033c4ca3c`
- deterministic gzip SHA-256: `a62272a186b5dd31d0f09a85a28eaf52cd5b805013c3b399122a1fab210f724d`
- compact summary SHA-256: `614479f8e5b4da4de46f7cced1d04f8d9bde3ebb685fb49ac6b30b8d51d50ef1`
- paired comparisons SHA-256: `c200c50407961d92257f7cb827c85d4ed7237a7a7cc36422920bcea720ec4104`

## v0.6 independent projection-seed confirmation

- summary: `v0.6_projection_seed_sensitivity_summary.json`
- decision record: `v0.6_projection_seed_sensitivity_decision.json`
- complete deterministic archive: `v0.6_projection_seed_sensitivity.json.gz`
- schema: `origins-dna.projection-seed-sensitivity.v1`
- decision schema: `origins-dna.projection-seed-sensitivity-decision.v1`
- generator: `scripts/run_projection_seed_sensitivity.py`
- protocol: `docs/PROJECTION_SEED_SENSITIVITY_PROTOCOL.md`
- projection master seed: `70701`
- independent projection seeds: `32`
- excluded exploratory pilot seeds: `16`
- matched dimensions: `12, 24, 48`
- uncompressed JSON SHA-256: `feeb5350584d9328474faa7886017b197268f6ccb0655d908ffa1f88f2f9802c`
- deterministic gzip SHA-256: `2e4710b9b2d05fe81157da40ad1a38dfe0a84ecca2d87543d50903defe5de4d4`
- compact summary SHA-256: `dc27c033028b17ed50f517715b1e74b226dd5914248cbbc57818d583c329599c`
- decision record SHA-256: `938f7ae53fbef40db022c008ec69c4eb66682ca2af9b7a982e7a098680f4dc00`
- overall verdict: `WARN_REVIEW_REQUIRED`
- promotion allowed: `false`, because preregistered null contrasts are displaced from zero

Independent executions at each validated source state produced deterministic artifacts. The uncompressed full JSON files are generated locally but not retained in Git; the deterministic compressed archives contain every seed-level record.

These artifacts are synthetic software results, not biological evidence.
