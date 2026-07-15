#!/usr/bin/env python3
"""Run the deterministic Origins-DNA v0.3 synthetic benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = _REPOSITORY_ROOT / "src"
if str(_SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(_SOURCE_ROOT))

from origins_dna import SyntheticMotifConfig, run_controlled_benchmark


def _parse_seeds(value: str) -> tuple[int, ...]:
    try:
        seeds = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError("seeds must be comma-separated integers") from error
    if not seeds:
        raise argparse.ArgumentTypeError("at least one seed is required")
    return seeds


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/v0.3_synthetic_motif_benchmark.json"),
    )
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--length", type=int, default=64)
    parser.add_argument("--motif", default="CGCG")
    parser.add_argument(
        "--seeds",
        type=_parse_seeds,
        default=(11, 23, 37, 41, 53, 67, 79, 83),
    )
    parser.add_argument("--bootstrap-resamples", type=int, default=2000)
    arguments = parser.parse_args()

    config = SyntheticMotifConfig(
        sample_count=arguments.samples,
        sequence_length=arguments.length,
        motif=arguments.motif,
        seeds=arguments.seeds,
        bootstrap_resamples=arguments.bootstrap_resamples,
    )
    report = run_controlled_benchmark(config)
    full_json, archive, summary = report.write_bundle(arguments.output)
    print(full_json)
    print(archive)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
