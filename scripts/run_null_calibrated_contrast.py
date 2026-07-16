#!/usr/bin/env python3
"""Run the preregistered Origins-DNA v0.7 null-calibrated contrast."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = _REPOSITORY_ROOT / "src"
if str(_SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(_SOURCE_ROOT))

from origins_dna.null_calibration import run_null_calibrated_analysis


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
        default=Path("results/v0.7_null_calibrated_contrast.json"),
    )
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--length", type=int, default=64)
    parser.add_argument("--motif", default="CGCG")
    parser.add_argument("--budgets", type=_parse_ints, default=(12, 24, 48))
    parser.add_argument("--projection-master-seed", type=int, default=80_801)
    parser.add_argument("--projection-seed-count", type=int, default=32)
    parser.add_argument("--data-master-seed", type=int, default=90_901)
    parser.add_argument("--data-seed-count", type=int, default=8)
    parser.add_argument("--bootstrap-resamples", type=int, default=2_000)
    parser.add_argument("--hierarchical-bootstrap-resamples", type=int, default=5_000)
    parser.add_argument("--cluster-sign-flip-resamples", type=int, default=20_000)
    parser.add_argument("--significance-level", type=float, default=0.05)
    parser.add_argument("--sign-stability-threshold", type=float, default=0.75)
    arguments = parser.parse_args()

    report = run_null_calibrated_analysis(
        sample_count=arguments.samples,
        sequence_length=arguments.length,
        motif=arguments.motif,
        budgets=arguments.budgets,
        projection_master_seed=arguments.projection_master_seed,
        projection_seed_count=arguments.projection_seed_count,
        data_master_seed=arguments.data_master_seed,
        data_seed_count=arguments.data_seed_count,
        bootstrap_resamples=arguments.bootstrap_resamples,
        hierarchical_bootstrap_resamples=arguments.hierarchical_bootstrap_resamples,
        cluster_sign_flip_resamples=arguments.cluster_sign_flip_resamples,
        significance_level=arguments.significance_level,
        sign_stability_threshold=arguments.sign_stability_threshold,
    )
    full_json, archive, summary = report.write_bundle(arguments.output)
    print(full_json)
    print(archive)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
