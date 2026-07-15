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

Independent executions at the validated source state produced byte-identical full JSON and gzip artifacts. The uncompressed full JSON files are generated locally but not retained in Git; the deterministic compressed archives contain every seed-level record.

These artifacts are synthetic software results, not biological evidence.
