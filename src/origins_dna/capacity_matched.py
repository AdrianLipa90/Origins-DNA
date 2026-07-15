"""Dimension-matched benchmark introduced in Origins-DNA v0.5.

The comparison uses deterministic, label-blind signed random projections to
place orbital multilag tensors and tetranucleotide frequencies in the same
output dimensions. Equal dimensionality is a capacity proxy, not a proof of
equal information content or equal effective model complexity.
"""

from __future__ import annotations

from functools import lru_cache
import json
import math
import random
from pathlib import Path
from typing import Sequence

from .benchmark_data import derived_seed
from .benchmark_features import (
    Representation,
    kmer_frequency_vector,
    orbital_multilag_tensor_vector,
)
from .benchmark_runner import bootstrap_mean_confidence_interval, run_registered_benchmark
from .benchmark_types import BenchmarkReport, SyntheticMotifConfig, Vector

DEFAULT_DIMENSION_BUDGETS = (12, 24, 48)
_ORBITAL_INPUT_DIMENSION = 48
_TETRANUCLEOTIDE_INPUT_DIMENSION = 256
_RAW_ORBITAL_NAME = "orbital_multilag_raw_48"
_RAW_TETRANUCLEOTIDE_NAME = "tetranucleotide_raw_256"


def make_signed_projection(
    input_dimension: int,
    output_dimension: int,
    *,
    seed: int,
) -> tuple[tuple[float, ...], ...]:
    """Create a deterministic Rademacher projection matrix.

    Matrix entries are +/-1/sqrt(output_dimension). The matrix is generated
    without labels, train/test data or task outcomes.
    """

    if not isinstance(input_dimension, int) or not isinstance(output_dimension, int):
        raise TypeError("projection dimensions must be integers")
    if input_dimension <= 0 or output_dimension <= 0:
        raise ValueError("projection dimensions must be positive")
    if not isinstance(seed, int):
        raise TypeError("seed must be an integer")

    rng = random.Random(seed)
    scale = 1.0 / math.sqrt(output_dimension)
    return tuple(
        tuple(scale if rng.getrandbits(1) else -scale for _ in range(input_dimension))
        for _ in range(output_dimension)
    )


def _apply_validated_projection(vector: Vector, matrix: tuple[tuple[float, ...], ...]) -> Vector:
    """Apply an already validated matrix without repeating structural scans."""

    return tuple(sum(weight * value for weight, value in zip(row, vector)) for row in matrix)


def apply_projection(vector: Sequence[float], matrix: Sequence[Sequence[float]]) -> Vector:
    """Apply a rectangular projection matrix to one finite numeric vector."""

    materialized = tuple(float(value) for value in vector)
    rows = tuple(tuple(float(value) for value in row) for row in matrix)
    if not rows:
        raise ValueError("matrix must not be empty")
    input_dimension = len(rows[0])
    if input_dimension <= 0:
        raise ValueError("matrix rows must not be empty")
    if any(len(row) != input_dimension for row in rows):
        raise ValueError("matrix rows must have equal length")
    if len(materialized) != input_dimension:
        raise ValueError(
            f"vector has dimension {len(materialized)}; expected {input_dimension}"
        )
    if any(not math.isfinite(value) for value in materialized):
        raise ValueError("vector values must be finite")
    if any(not math.isfinite(value) for row in rows for value in row):
        raise ValueError("matrix values must be finite")
    return _apply_validated_projection(materialized, rows)


def projected_representation(
    representation: Representation,
    *,
    input_dimension: int,
    output_dimension: int,
    seed: int,
) -> Representation:
    """Wrap a representation in one precomputed label-blind projection."""

    matrix = make_signed_projection(input_dimension, output_dimension, seed=seed)

    @lru_cache(maxsize=None)
    def transform(sequence: str) -> Vector:
        vector = representation(sequence)
        if len(vector) != input_dimension:
            raise ValueError(
                f"representation produced dimension {len(vector)}; expected {input_dimension}"
            )
        return _apply_validated_projection(vector, matrix)

    return transform


def _validate_budgets(budgets: Sequence[int]) -> tuple[int, ...]:
    materialized = tuple(budgets)
    if not materialized:
        raise ValueError("budgets must not be empty")
    if any(not isinstance(value, int) for value in materialized):
        raise TypeError("budgets must contain integers")
    if any(value <= 0 for value in materialized):
        raise ValueError("budgets must be positive")
    if len(set(materialized)) != len(materialized):
        raise ValueError("budgets must be unique")
    return tuple(sorted(materialized))


def dimension_matched_representation_registry(
    sequence_length: int,
    *,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
    projection_seed: int = 50_501,
) -> dict[str, Representation]:
    """Return raw anchors plus matched orbital/4-mer projection pairs."""

    if sequence_length <= 3:
        raise ValueError("sequence_length must exceed 3")
    materialized_budgets = _validate_budgets(budgets)
    if not isinstance(projection_seed, int):
        raise TypeError("projection_seed must be an integer")

    @lru_cache(maxsize=None)
    def orbital_source(sequence: str) -> Vector:
        return orbital_multilag_tensor_vector(sequence)

    @lru_cache(maxsize=None)
    def tetranucleotide_source(sequence: str) -> Vector:
        return kmer_frequency_vector(sequence, 4)

    registry: dict[str, Representation] = {
        _RAW_ORBITAL_NAME: orbital_source,
        _RAW_TETRANUCLEOTIDE_NAME: tetranucleotide_source,
    }
    for ordinal, budget in enumerate(materialized_budgets):
        family_seed = derived_seed(projection_seed, ordinal)
        registry[f"orbital_projected_{budget}"] = projected_representation(
            orbital_source,
            input_dimension=_ORBITAL_INPUT_DIMENSION,
            output_dimension=budget,
            seed=derived_seed(family_seed, 1),
        )
        registry[f"tetranucleotide_projected_{budget}"] = projected_representation(
            tetranucleotide_source,
            input_dimension=_TETRANUCLEOTIDE_INPUT_DIMENSION,
            output_dimension=budget,
            seed=derived_seed(family_seed, 2),
        )
    return registry


def run_dimension_matched_benchmark(
    config: SyntheticMotifConfig,
    *,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
    projection_seed: int = 50_501,
) -> BenchmarkReport:
    """Run equal-output-dimension comparisons under the frozen v0.4 protocol."""

    return run_registered_benchmark(
        config,
        dimension_matched_representation_registry(
            config.sequence_length,
            budgets=budgets,
            projection_seed=projection_seed,
        ),
        schema="origins-dna.dimension-matched-benchmark.v1",
        benchmark_status="synthetic_dimension_matched_capacity_proxy_result",
        interpretation_boundary=(
            "Equal output dimension and an identical classifier are controlled here, but this is "
            "only a capacity proxy. Results do not establish equal information content, a biological "
            "mechanism, general superiority or DNA decryption."
        ),
    )


def paired_dimension_difference(
    report: BenchmarkReport,
    budget: int,
    *,
    condition: str = "observed_labels",
    seed: int = 50_901,
) -> dict[str, object]:
    """Compare projected 4-mer minus projected orbital accuracy for one budget."""

    if not isinstance(budget, int):
        raise TypeError("budget must be an integer")
    if budget <= 0:
        raise ValueError("budget must be positive")
    left = f"tetranucleotide_projected_{budget}"
    right = f"orbital_projected_{budget}"
    by_key = {
        (result.seed, result.representation, result.condition): result.balanced_accuracy
        for result in report.seed_results
    }
    differences: list[float] = []
    for benchmark_seed in report.config.seeds:
        try:
            differences.append(
                by_key[(benchmark_seed, left, condition)]
                - by_key[(benchmark_seed, right, condition)]
            )
        except KeyError as error:
            raise ValueError(
                f"missing matched result for budget={budget}, seed={benchmark_seed}, "
                f"condition={condition}"
            ) from error
    ci_low, ci_high = bootstrap_mean_confidence_interval(
        differences,
        confidence=report.config.confidence,
        resamples=report.config.bootstrap_resamples,
        seed=seed,
    )
    return {
        "budget": budget,
        "left": left,
        "right": right,
        "condition": condition,
        "seed_count": len(differences),
        "mean_balanced_accuracy_difference": math.fsum(differences) / len(differences),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "per_seed_differences": differences,
    }


def _budgets_from_report(report: BenchmarkReport) -> tuple[int, ...]:
    prefix = "orbital_projected_"
    budgets = tuple(
        sorted(
            int(name.removeprefix(prefix))
            for name in report.representations
            if name.startswith(prefix)
        )
    )
    if not budgets:
        raise ValueError("report contains no matched orbital representations")
    for budget in budgets:
        if f"tetranucleotide_projected_{budget}" not in report.representations:
            raise ValueError(f"report is missing tetranucleotide pair for budget {budget}")
    return budgets


def dimension_matched_comparison_summary(report: BenchmarkReport) -> dict[str, object]:
    """Build paired comparisons for every declared output dimension."""

    if report.schema != "origins-dna.dimension-matched-benchmark.v1":
        raise ValueError("report must be a dimension-matched benchmark")
    budgets = _budgets_from_report(report)
    comparisons = [
        paired_dimension_difference(
            report,
            budget,
            seed=derived_seed(50_901, ordinal),
        )
        for ordinal, budget in enumerate(budgets)
    ]
    return {
        "schema": "origins-dna.dimension-matched-comparisons.v1",
        "source_report_schema": report.schema,
        "metric": "balanced_accuracy",
        "capacity_control": "equal output dimension, identical classifier and identical split",
        "projection": "deterministic label-blind Rademacher random projection",
        "budgets": list(budgets),
        "comparisons": comparisons,
        "interpretation_boundary": report.interpretation_boundary,
    }


def write_dimension_matched_bundle(
    report: BenchmarkReport,
    path: str | Path,
) -> tuple[Path, Path, Path, Path]:
    """Write deterministic full, gzip, summary and paired-comparison artifacts."""

    full_json, archive, summary = report.write_bundle(path)
    summary_payload = json.loads(summary.read_text(encoding="utf-8"))
    summary.write_text(
        json.dumps(summary_payload, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    comparisons = full_json.with_name(full_json.stem + "_comparisons.json")
    comparisons.write_text(
        json.dumps(
            dimension_matched_comparison_summary(report),
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n",
        encoding="utf-8",
    )
    return full_json, archive, summary, comparisons
