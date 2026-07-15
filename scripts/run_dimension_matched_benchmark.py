#!/usr/bin/env python3
"""Run the deterministic Origins-DNA v0.5 dimension-matched benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = _REPOSITORY_ROOT / "src"
if str(_SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(_SOURCE_ROOT))

from origins_dna.capacity_matched import (
    DEFAULT_DIMENSION_BUDGETS,
    run_dimension_matched_benchmark,
    write_dimension_matched_bundle,
)
from origins_dna.benchmark_types import SyntheticMotifConfig


def _parse_integer_tuple(value: str, *, label: str) -> tuple[int, ...]:
    try:
        values = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"{label} must be comma-separated integers") from error
    if not values:
        raise argparse.ArgumentTypeError(f"at least one {label[:-1]} is required")
    return values


def _parse_seeds(value: str) -> tuple[int, ...]:
    return _parse_integer_tuple(value, label="seeds")


def _parse_budgets(value: str) -> tuple[int, ...]:
    return _parse_integer_tuple(value, label="budgets")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/v0.5_dimension_matched_benchmark.json"),
    )
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--length", type=int, default=64)
    parser.add_argument("--motif", default="CGCG")
    parser.add_argument(
        "--seeds",
        type=_parse_seeds,
        default=(11, 23, 37, 41, 53, 67, 79, 83),
    )
    parser.add_argument(
        "--budgets",
        type=_parse_budgets,
        default=DEFAULT_DIMENSION_BUDGETS,
    )
    parser.add_argument("--projection-seed", type=int, default=50_501)
    parser.add_argument("--bootstrap-resamples", type=int, default=2000)
    arguments = parser.parse_args()

    report = run_dimension_matched_benchmark(
        SyntheticMotifConfig(
            sample_count=arguments.samples,
            sequence_length=arguments.length,
            motif=arguments.motif,
            seeds=arguments.seeds,
            bootstrap_resamples=arguments.bootstrap_resamples,
        ),
        budgets=arguments.budgets,
        projection_seed=arguments.projection_seed,
    )
    outputs = write_dimension_matched_bundle(report, arguments.output)
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
