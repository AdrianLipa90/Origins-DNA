#!/usr/bin/env python3
"""Run the independent Origins-DNA v0.6 projection-seed sensitivity analysis."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = _REPOSITORY_ROOT / "src"
if str(_SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(_SOURCE_ROOT))

from origins_dna.benchmark_types import SyntheticMotifConfig
from origins_dna.projection_sensitivity import run_projection_seed_sensitivity


def _parse_ints(value: str) -> tuple[int, ...]:
    try:
        values = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError("values must be comma-separated integers") from error
    if not values:
        raise argparse.ArgumentTypeError("at least one integer is required")
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/v0.6_projection_seed_sensitivity.json"),
    )
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--length", type=int, default=64)
    parser.add_argument("--motif", default="CGCG")
    parser.add_argument(
        "--data-seeds",
        type=_parse_ints,
        default=(11, 23, 37, 41, 53, 67, 79, 83),
    )
    parser.add_argument("--budgets", type=_parse_ints, default=(12, 24, 48))
    parser.add_argument("--projection-master-seed", type=int, default=70_701)
    parser.add_argument("--projection-seed-count", type=int, default=32)
    parser.add_argument("--projection-seeds", type=_parse_ints)
    parser.add_argument("--bootstrap-resamples", type=int, default=2_000)
    parser.add_argument("--hierarchical-bootstrap-resamples", type=int, default=5_000)
    arguments = parser.parse_args()

    report = run_projection_seed_sensitivity(
        SyntheticMotifConfig(
            sample_count=arguments.samples,
            sequence_length=arguments.length,
            motif=arguments.motif,
            seeds=arguments.data_seeds,
            bootstrap_resamples=arguments.bootstrap_resamples,
        ),
        budgets=arguments.budgets,
        projection_master_seed=arguments.projection_master_seed,
        projection_seed_count=arguments.projection_seed_count,
        projection_seeds=arguments.projection_seeds,
        hierarchical_bootstrap_resamples=arguments.hierarchical_bootstrap_resamples,
    )
    full_json, archive, summary = report.write_bundle(arguments.output)
    print(full_json)
    print(archive)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
