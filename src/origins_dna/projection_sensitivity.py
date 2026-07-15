"""Projection-seed sensitivity analysis introduced in Origins-DNA v0.6.

The confirmatory schedule is independent of the exploratory 16-seed pilot. The
pilot is retained as process evidence and is not included in confirmatory
aggregates. Results remain synthetic software evidence only.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import gzip
import json
import math
from pathlib import Path
import random
import statistics
from typing import Sequence

from .benchmark_data import derived_seed
from .benchmark_types import SyntheticMotifConfig
from .capacity_matched import (
    DEFAULT_DIMENSION_BUDGETS,
    paired_dimension_difference,
    run_dimension_matched_benchmark,
)

DEFAULT_PROJECTION_MASTER_SEED = 70_701
DEFAULT_PROJECTION_SEED_COUNT = 32
DEFAULT_HIERARCHICAL_BOOTSTRAP_RESAMPLES = 5_000
EXPLORATORY_PILOT_PROJECTION_SEEDS = (
    50_501,
    471_639_659,
    471_737_068,
    471_834_477,
    471_931_886,
    472_029_295,
    472_126_704,
    472_224_113,
    472_321_522,
    472_418_931,
    472_516_340,
    472_613_749,
    472_711_158,
    472_808_567,
    472_905_976,
    473_003_385,
)
_CONDITIONS = (
    "observed_labels",
    "label_permutation_null",
    "sequence_shuffle_null",
)


@dataclass(frozen=True, slots=True)
class ProjectionSeedResult:
    projection_seed: int
    budget: int
    condition: str
    mean_balanced_accuracy_difference: float
    per_data_seed_differences: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class ProjectionSensitivityAggregate:
    budget: int
    condition: str
    projection_seed_count: int
    data_seed_count: int
    mean_balanced_accuracy_difference: float
    median_projection_mean_difference: float
    hierarchical_ci_low: float
    hierarchical_ci_high: float
    fraction_projection_means_positive: float
    fraction_projection_means_negative: float
    minimum_projection_mean_difference: float
    maximum_projection_mean_difference: float


@dataclass(frozen=True, slots=True)
class ProjectionSensitivityReport:
    schema: str
    benchmark_status: str
    config: SyntheticMotifConfig
    budgets: tuple[int, ...]
    projection_master_seed: int
    projection_seeds: tuple[int, ...]
    excluded_exploratory_projection_seeds: tuple[int, ...]
    hierarchical_bootstrap_resamples: int
    seed_results: tuple[ProjectionSeedResult, ...]
    aggregate_results: tuple[ProjectionSensitivityAggregate, ...]
    interpretation_boundary: str

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": self.schema,
            "benchmark_status": self.benchmark_status,
            "config": asdict(self.config),
            "budgets": list(self.budgets),
            "projection_master_seed": self.projection_master_seed,
            "projection_seeds": list(self.projection_seeds),
            "excluded_exploratory_projection_seeds": list(
                self.excluded_exploratory_projection_seeds
            ),
            "hierarchical_bootstrap_resamples": self.hierarchical_bootstrap_resamples,
            "seed_results": [asdict(result) for result in self.seed_results],
            "aggregate_results": [asdict(result) for result in self.aggregate_results],
            "interpretation_boundary": self.interpretation_boundary,
        }

    def summary_dict(self, *, full_result_archive: str | None = None) -> dict[str, object]:
        payload = self.to_dict()
        payload.pop("seed_results")
        if full_result_archive is not None:
            payload["full_result_archive"] = full_result_archive
        return payload

    def write_bundle(self, path: str | Path) -> tuple[Path, Path, Path]:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )
        archive = output.with_name(output.name + ".gz")
        with archive.open("wb") as raw_file:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_file, mtime=0) as compressed:
                compressed.write(output.read_bytes())
        summary = output.with_name(output.stem + "_summary.json")
        summary.write_text(
            json.dumps(
                self.summary_dict(full_result_archive=archive.name),
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
        return output, archive, summary


def projection_seed_schedule(
    *,
    master_seed: int = DEFAULT_PROJECTION_MASTER_SEED,
    count: int = DEFAULT_PROJECTION_SEED_COUNT,
    excluded: Sequence[int] = EXPLORATORY_PILOT_PROJECTION_SEEDS,
) -> tuple[int, ...]:
    """Return a deterministic unique schedule excluding exploratory seeds."""

    if not isinstance(master_seed, int):
        raise TypeError("master_seed must be an integer")
    if not isinstance(count, int):
        raise TypeError("count must be an integer")
    if count <= 0:
        raise ValueError("count must be positive")
    excluded_values = tuple(excluded)
    if any(not isinstance(value, int) for value in excluded_values):
        raise TypeError("excluded must contain integers")
    excluded_set = set(excluded_values)
    output: list[int] = []
    ordinal = 0
    while len(output) < count:
        candidate = derived_seed(master_seed, ordinal)
        ordinal += 1
        if candidate in excluded_set or candidate in output:
            continue
        output.append(candidate)
    return tuple(output)


def run_projection_seed_sensitivity(
    config: SyntheticMotifConfig,
    *,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
    projection_master_seed: int = DEFAULT_PROJECTION_MASTER_SEED,
    projection_seed_count: int = DEFAULT_PROJECTION_SEED_COUNT,
    projection_seeds: Sequence[int] | None = None,
    hierarchical_bootstrap_resamples: int = DEFAULT_HIERARCHICAL_BOOTSTRAP_RESAMPLES,
) -> ProjectionSensitivityReport:
    """Run the frozen v0.5 comparison over independent projection seeds."""

    materialized_budgets = _validate_positive_unique_ints(budgets, name="budgets")
    if projection_seeds is None:
        materialized_projection_seeds = projection_seed_schedule(
            master_seed=projection_master_seed,
            count=projection_seed_count,
        )
    else:
        materialized_projection_seeds = _validate_positive_unique_ints(
            projection_seeds,
            name="projection_seeds",
        )
        exploratory_overlap = set(materialized_projection_seeds).intersection(
            EXPLORATORY_PILOT_PROJECTION_SEEDS
        )
        if exploratory_overlap:
            raise ValueError(
                "confirmatory projection_seeds overlap exploratory pilot seeds: "
                f"{sorted(exploratory_overlap)}"
            )
    if hierarchical_bootstrap_resamples <= 0:
        raise ValueError("hierarchical_bootstrap_resamples must be positive")

    seed_results: list[ProjectionSeedResult] = []
    for projection_seed in materialized_projection_seeds:
        report = run_dimension_matched_benchmark(
            config,
            budgets=materialized_budgets,
            projection_seed=projection_seed,
        )
        for condition in _CONDITIONS:
            for budget in materialized_budgets:
                comparison = paired_dimension_difference(
                    report,
                    budget,
                    condition=condition,
                    seed=derived_seed(projection_seed, budget),
                )
                seed_results.append(
                    ProjectionSeedResult(
                        projection_seed=projection_seed,
                        budget=budget,
                        condition=condition,
                        mean_balanced_accuracy_difference=float(
                            comparison["mean_balanced_accuracy_difference"]
                        ),
                        per_data_seed_differences=tuple(
                            float(value) for value in comparison["per_seed_differences"]
                        ),
                    )
                )

    aggregates = _aggregate_projection_seed_results(
        seed_results,
        budgets=materialized_budgets,
        projection_seeds=materialized_projection_seeds,
        data_seed_count=len(config.seeds),
        confidence=config.confidence,
        resamples=hierarchical_bootstrap_resamples,
        seed=derived_seed(projection_master_seed, 60_699),
    )
    return ProjectionSensitivityReport(
        schema="origins-dna.projection-seed-sensitivity.v1",
        benchmark_status="synthetic_projection_seed_sensitivity_result",
        config=config,
        budgets=materialized_budgets,
        projection_master_seed=projection_master_seed,
        projection_seeds=materialized_projection_seeds,
        excluded_exploratory_projection_seeds=EXPLORATORY_PILOT_PROJECTION_SEEDS,
        hierarchical_bootstrap_resamples=hierarchical_bootstrap_resamples,
        seed_results=tuple(seed_results),
        aggregate_results=aggregates,
        interpretation_boundary=(
            "This independent projection-seed confirmation follows an exploratory pilot that was "
            "observed before repository preregistration. It quantifies synthetic software sensitivity "
            "only and does not establish equal information content, a biological mechanism, general "
            "representational superiority or DNA decryption."
        ),
    )


def _validate_positive_unique_ints(values: Sequence[int], *, name: str) -> tuple[int, ...]:
    materialized = tuple(values)
    if not materialized:
        raise ValueError(f"{name} must not be empty")
    if any(not isinstance(value, int) for value in materialized):
        raise TypeError(f"{name} must contain integers")
    if any(value <= 0 for value in materialized):
        raise ValueError(f"{name} must be positive")
    if len(set(materialized)) != len(materialized):
        raise ValueError(f"{name} must be unique")
    return tuple(sorted(materialized))


def _aggregate_projection_seed_results(
    results: Sequence[ProjectionSeedResult],
    *,
    budgets: Sequence[int],
    projection_seeds: Sequence[int],
    data_seed_count: int,
    confidence: float,
    resamples: int,
    seed: int,
) -> tuple[ProjectionSensitivityAggregate, ...]:
    by_key = {
        (result.projection_seed, result.budget, result.condition): result
        for result in results
    }
    output: list[ProjectionSensitivityAggregate] = []
    ordinal = 0
    for condition in _CONDITIONS:
        for budget in budgets:
            group = tuple(
                by_key[(projection_seed, budget, condition)]
                for projection_seed in projection_seeds
            )
            projection_means = tuple(
                result.mean_balanced_accuracy_difference for result in group
            )
            ci_low, ci_high = hierarchical_bootstrap_confidence_interval(
                tuple(result.per_data_seed_differences for result in group),
                confidence=confidence,
                resamples=resamples,
                seed=derived_seed(seed, ordinal),
            )
            ordinal += 1
            output.append(
                ProjectionSensitivityAggregate(
                    budget=budget,
                    condition=condition,
                    projection_seed_count=len(group),
                    data_seed_count=data_seed_count,
                    mean_balanced_accuracy_difference=(
                        math.fsum(projection_means) / len(projection_means)
                    ),
                    median_projection_mean_difference=statistics.median(projection_means),
                    hierarchical_ci_low=ci_low,
                    hierarchical_ci_high=ci_high,
                    fraction_projection_means_positive=(
                        sum(value > 0.0 for value in projection_means) / len(projection_means)
                    ),
                    fraction_projection_means_negative=(
                        sum(value < 0.0 for value in projection_means) / len(projection_means)
                    ),
                    minimum_projection_mean_difference=min(projection_means),
                    maximum_projection_mean_difference=max(projection_means),
                )
            )
    return tuple(output)


def hierarchical_bootstrap_confidence_interval(
    nested_values: Sequence[Sequence[float]],
    *,
    confidence: float,
    resamples: int,
    seed: int,
) -> tuple[float, float]:
    """Bootstrap projection seeds, then data seeds within each selected projection."""

    groups = tuple(tuple(float(value) for value in group) for group in nested_values)
    if not groups or any(not group for group in groups):
        raise ValueError("nested_values must contain non-empty groups")
    if any(not math.isfinite(value) for group in groups for value in group):
        raise ValueError("nested_values must be finite")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie strictly between 0 and 1")
    if resamples <= 0:
        raise ValueError("resamples must be positive")

    rng = random.Random(seed)
    projection_count = len(groups)
    data_count = len(groups[0])
    if any(len(group) != data_count for group in groups):
        raise ValueError("all projection groups must contain equal data-seed counts")
    means: list[float] = []
    for _ in range(resamples):
        total = 0.0
        for _ in range(projection_count):
            group = groups[rng.randrange(projection_count)]
            total += math.fsum(group[rng.randrange(data_count)] for _ in range(data_count))
        means.append(total / (projection_count * data_count))
    means.sort()
    tail = (1.0 - confidence) / 2.0
    low_index = min(resamples - 1, max(0, math.floor(tail * resamples)))
    high_index = min(resamples - 1, max(0, math.ceil((1.0 - tail) * resamples) - 1))
    return means[low_index], means[high_index]
