#!/usr/bin/env python3
"""Run the frozen Origins-DNA v0.9.1 three-quarter causal benchmark."""
from __future__ import annotations

import argparse
from pathlib import Path

from origins_dna.benchmark_types import SyntheticMotifConfig
from origins_dna.causal_three_quarter_benchmark import (
    run_three_quarter_benchmark,
    write_three_quarter_bundle,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/v0.9.1_three_quarter_causal_benchmark.json"),
    )
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--length", type=int, default=64)
    parser.add_argument("--motif", default="CGCG")
    parser.add_argument("--projection-seed", type=int, default=134091)
    parser.add_argument("--bootstrap-resamples", type=int, default=2000)
    return parser


def main() -> None:
    args = _parser().parse_args()
    config = SyntheticMotifConfig(
        sample_count=args.samples,
        sequence_length=args.length,
        motif=args.motif,
        bootstrap_resamples=args.bootstrap_resamples,
    )
    report = run_three_quarter_benchmark(
        config,
        projection_seed=args.projection_seed,
    )
    outputs = write_three_quarter_bundle(report, args.output)
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
