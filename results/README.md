# Results provenance

The repository retains two deterministic artifacts for the v0.3 synthetic software benchmark:

- `v0.3_synthetic_motif_benchmark_summary.json`: aggregate, review-friendly JSON;
- `v0.3_synthetic_motif_benchmark.json.gz`: complete result including every seed-level record.

Provenance:

- schema: `origins-dna.benchmark.v1`
- generator: `scripts/run_controlled_benchmark.py`
- default configuration: declared in `docs/BENCHMARK_PROTOCOL.md`
- uncompressed full JSON SHA-256: `07ed3bc941821785131b61c667ca9a81957fc187d5ac7c33d1875a5f83b6ff21`
- deterministic gzip SHA-256: `01f44457c8932b7d799ef0bc09926a76121017084216d1fc01edd9f9d7300545`
- compact summary SHA-256: `ff169e3fb076b8db8e512f78b8758e239b4339b89efce77c525d6524f864d3ab`
- reproducibility check: two independent executions in the same validated source state produced byte-identical artifacts

The artifacts are not biological evidence. They record the performance of defined software representations on a synthetic task.
