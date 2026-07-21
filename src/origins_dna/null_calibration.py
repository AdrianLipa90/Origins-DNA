"""Null-calibrated observed-minus-null contrasts for Origins-DNA v0.7.

The v0.7 analysis uses projection and data seeds disjoint from the exploratory
v0.5/v0.6 schedules. It preserves label-permutation and sequence-shuffle nulls
as separate estimands; they are never pooled into one convenient baseline.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import gzip
import json
import math
from pathlib import Path
import random
import statistics
from typing import Sequence

from .benchmark_data import (
    derived_seed,
    generate_synthetic_motif_dataset,
    shuffle_sequence,
)
from .benchmark_types import SyntheticMotifConfig
from .capacity_matched import DEFAULT_DIMENSION_BUDGETS
from .projection_sensitivity import (
    EXPLORATORY_PILOT_PROJECTION_SEEDS,
    ProjectionSensitivityReport,
    hierarchical_bootstrap_confidence_interval,
    projection_seed_schedule,
    run_projection_seed_sensitivity,
)

DEFAULT_NULL_CALIBRATION_PROJECTION_MASTER_SEED = 80_801
DEFAULT_NULL_CALIBRATION_PROJECTION_SEED_COUNT = 32
DEFAULT_NULL_CALIBRATION_DATA_MASTER_SEED = 90_901
DEFAULT_NULL_CALIBRATION_DATA_SEED_COUNT = 8
DEFAULT_CLUSTER_SIGN_FLIP_RESAMPLES = 20_000
DEFAULT_SIGN_STABILITY_THRESHOLD = 0.75
DEFAULT_SIGNIFICANCE_LEVEL = 0.05
_PRIOR_DATA_SEEDS = (11, 23, 37, 41, 53, 67, 79, 83)
_NULL_CONDITIONS = ("label_permutation_null", "sequence_shuffle_null")


@dataclass(frozen=True, slots=True)
class NullCalibratedSeedResult:
    projection_seed: int
    budget: int
    null_condition: str
    mean_calibrated_difference: float
    per_data_seed_calibrated_differences: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class NullCalibratedAggregate:
    budget: int
    null_condition: str
    projection_seed_count: int
    data_seed_count: int
    mean_calibrated_difference: float
    median_projection_mean_difference: float
    hierarchical_ci_low: float
    hierarchical_ci_high: float
    fraction_projection_means_positive: float
    fraction_projection_means_negative: float
    minimum_projection_mean_difference: float
    maximum_projection_mean_difference: float
    cluster_sign_flip_p_value: float
    direction: str
    interval_gate_passed: bool
    sign_stability_gate_passed: bool
    randomization_gate_passed: bool
    calibrated_gate_passed: bool


@dataclass(frozen=True, slots=True)
class NullCalibratedDecision:
    budget: int
    label_permutation_direction: str
    sequence_shuffle_direction: str
    consensus_direction: str
    both_nulls_passed: bool
    promotion_gate_passed: bool


@dataclass(frozen=True, slots=True)
class NullControlDiagnostics:
    data_seed_count: int
    sequence_count: int
    maximum_per_sequence_base_count_residual_after_shuffle: int
    shuffle_composition_gate_passed: bool
    rejection_sampling_status: str


@dataclass(frozen=True, slots=True)
class NullCalibratedReport:
    schema: str
    benchmark_status: str
    source_report_schema: str
    config: SyntheticMotifConfig
    budgets: tuple[int, ...]
    projection_master_seed: int
    projection_seeds: tuple[int, ...]
    data_master_seed: int
    excluded_prior_projection_seeds: tuple[int, ...]
    excluded_prior_data_seeds: tuple[int, ...]
    hierarchical_bootstrap_resamples: int
    cluster_sign_flip_resamples: int
    significance_level: float
    sign_stability_threshold: float
    seed_results: tuple[NullCalibratedSeedResult, ...]
    aggregate_results: tuple[NullCalibratedAggregate, ...]
    decisions: tuple[NullCalibratedDecision, ...]
    diagnostics: NullControlDiagnostics
    overall_verdict: str
    promotion_allowed: bool
    biological_progression_allowed: bool
    interpretation_boundary: str

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": self.schema,
            "benchmark_status": self.benchmark_status,
            "source_report_schema": self.source_report_schema,
            "config": asdict(self.config),
            "budgets": list(self.budgets),
            "projection_master_seed": self.projection_master_seed,
            "projection_seeds": list(self.projection_seeds),
            "data_master_seed": self.data_master_seed,
            "excluded_prior_projection_seeds": list(self.excluded_prior_projection_seeds),
            "excluded_prior_data_seeds": list(self.excluded_prior_data_seeds),
            "hierarchical_bootstrap_resamples": self.hierarchical_bootstrap_resamples,
            "cluster_sign_flip_resamples": self.cluster_sign_flip_resamples,
            "significance_level": self.significance_level,
            "sign_stability_threshold": self.sign_stability_threshold,
            "seed_results": [asdict(result) for result in self.seed_results],
            "aggregate_results": [asdict(result) for result in self.aggregate_results],
            "decisions": [asdict(result) for result in self.decisions],
            "diagnostics": asdict(self.diagnostics),
            "overall_verdict": self.overall_verdict,
            "promotion_allowed": self.promotion_allowed,
            "biological_progression_allowed": self.biological_progression_allowed,
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


def _unique_derived_schedule(
    *,
    master_seed: int,
    count: int,
    excluded: Sequence[int],
) -> tuple[int, ...]:
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
        if candidate <= 0 or candidate in excluded_set or candidate in output:
            continue
        output.append(candidate)
    return tuple(output)


def prior_projection_seeds() -> tuple[int, ...]:
    """Return every projection seed observed before the v0.7 confirmation."""

    v06_confirmatory = projection_seed_schedule(master_seed=70_701, count=32)
    return tuple(sorted(set(EXPLORATORY_PILOT_PROJECTION_SEEDS).union(v06_confirmatory)))


def null_calibration_projection_seed_schedule(
    *,
    master_seed: int = DEFAULT_NULL_CALIBRATION_PROJECTION_MASTER_SEED,
    count: int = DEFAULT_NULL_CALIBRATION_PROJECTION_SEED_COUNT,
) -> tuple[int, ...]:
    return _unique_derived_schedule(
        master_seed=master_seed,
        count=count,
        excluded=prior_projection_seeds(),
    )


def null_calibration_data_seed_schedule(
    *,
    master_seed: int = DEFAULT_NULL_CALIBRATION_DATA_MASTER_SEED,
    count: int = DEFAULT_NULL_CALIBRATION_DATA_SEED_COUNT,
) -> tuple[int, ...]:
    return _unique_derived_schedule(
        master_seed=master_seed,
        count=count,
        excluded=_PRIOR_DATA_SEEDS,
    )


def cluster_sign_flip_p_value(
    projection_means: Sequence[float],
    *,
    resamples: int = DEFAULT_CLUSTER_SIGN_FLIP_RESAMPLES,
    seed: int,
) -> float:
    """Two-sided projection-cluster sign-flip randomization p-value."""

    values = tuple(float(value) for value in projection_means)
    if not values:
        raise ValueError("projection_means must not be empty")
    if any(not math.isfinite(value) for value in values):
        raise ValueError("projection_means must be finite")
    if resamples <= 0:
        raise ValueError("resamples must be positive")
    observed = abs(math.fsum(values) / len(values))
    if observed == 0.0 and all(value == 0.0 for value in values):
        return 1.0
    rng = random.Random(seed)
    exceedances = 0
    for _ in range(resamples):
        permuted = math.fsum(
            value if rng.getrandbits(1) else -value for value in values
        ) / len(values)
        if abs(permuted) >= observed - 1e-15:
            exceedances += 1
    return (exceedances + 1.0) / (resamples + 1.0)


def sequence_shuffle_composition_diagnostics(
    config: SyntheticMotifConfig,
) -> NullControlDiagnostics:
    """Verify the exact per-sequence base-count invariant of the shuffle null."""

    maximum_residual = 0
    sequence_count = 0
    for data_seed in config.seeds:
        dataset = generate_synthetic_motif_dataset(config, seed=data_seed)
        for index, sequence in enumerate(dataset.sequences):
            shuffled = shuffle_sequence(
                sequence,
                seed=derived_seed(data_seed, 10_000 + index),
            )
            for base in "ACGT":
                maximum_residual = max(
                    maximum_residual,
                    abs(sequence.count(base) - shuffled.count(base)),
                )
            sequence_count += 1
    return NullControlDiagnostics(
        data_seed_count=len(config.seeds),
        sequence_count=sequence_count,
        maximum_per_sequence_base_count_residual_after_shuffle=maximum_residual,
        shuffle_composition_gate_passed=maximum_residual == 0,
        rejection_sampling_status=(
            "unresolved_design_limitation: negative sequences are generated by motif-avoidance "
            "rejection sampling; v0.7 records but does not causally remove that generator effect"
        ),
    )


def calibrate_projection_sensitivity_report(
    source_report: ProjectionSensitivityReport,
    *,
    cluster_sign_flip_resamples: int = DEFAULT_CLUSTER_SIGN_FLIP_RESAMPLES,
    significance_level: float = DEFAULT_SIGNIFICANCE_LEVEL,
    sign_stability_threshold: float = DEFAULT_SIGN_STABILITY_THRESHOLD,
) -> NullCalibratedReport:
    """Subtract each null from observed labels within every paired seed cell."""

    if source_report.schema != "origins-dna.projection-seed-sensitivity.v1":
        raise ValueError("source_report must be a projection-seed-sensitivity report")
    if not 0.0 < significance_level < 1.0:
        raise ValueError("significance_level must lie strictly between 0 and 1")
    if not 0.5 < sign_stability_threshold <= 1.0:
        raise ValueError("sign_stability_threshold must lie in (0.5, 1]")
    if cluster_sign_flip_resamples <= 0:
        raise ValueError("cluster_sign_flip_resamples must be positive")

    by_key = {
        (result.projection_seed, result.budget, result.condition): result
        for result in source_report.seed_results
    }
    seed_results: list[NullCalibratedSeedResult] = []
    for projection_seed in source_report.projection_seeds:
        for budget in source_report.budgets:
            observed = by_key[(projection_seed, budget, "observed_labels")]
            for null_condition in _NULL_CONDITIONS:
                null = by_key[(projection_seed, budget, null_condition)]
                if len(observed.per_data_seed_differences) != len(
                    null.per_data_seed_differences
                ):
                    raise ValueError("observed and null data-seed counts differ")
                calibrated = tuple(
                    observed_value - null_value
                    for observed_value, null_value in zip(
                        observed.per_data_seed_differences,
                        null.per_data_seed_differences,
                    )
                )
                seed_results.append(
                    NullCalibratedSeedResult(
                        projection_seed=projection_seed,
                        budget=budget,
                        null_condition=null_condition,
                        mean_calibrated_difference=math.fsum(calibrated) / len(calibrated),
                        per_data_seed_calibrated_differences=calibrated,
                    )
                )

    aggregates = _aggregate_calibrated_results(
        seed_results,
        budgets=source_report.budgets,
        projection_seeds=source_report.projection_seeds,
        data_seed_count=len(source_report.config.seeds),
        confidence=source_report.config.confidence,
        hierarchical_bootstrap_resamples=source_report.hierarchical_bootstrap_resamples,
        cluster_sign_flip_resamples=cluster_sign_flip_resamples,
        significance_level=significance_level,
        sign_stability_threshold=sign_stability_threshold,
        seed=derived_seed(source_report.projection_master_seed, 70_701),
    )
    decisions = build_null_calibrated_decisions(aggregates, budgets=source_report.budgets)
    diagnostics = sequence_shuffle_composition_diagnostics(source_report.config)
    all_passed = all(decision.promotion_gate_passed for decision in decisions)
    promotion_allowed = all_passed and diagnostics.shuffle_composition_gate_passed
    return NullCalibratedReport(
        schema="origins-dna.null-calibrated-contrast.v1",
        benchmark_status="synthetic_independent_null_calibrated_confirmation",
        source_report_schema=source_report.schema,
        config=source_report.config,
        budgets=source_report.budgets,
        projection_master_seed=source_report.projection_master_seed,
        projection_seeds=source_report.projection_seeds,
        data_master_seed=DEFAULT_NULL_CALIBRATION_DATA_MASTER_SEED,
        excluded_prior_projection_seeds=prior_projection_seeds(),
        excluded_prior_data_seeds=_PRIOR_DATA_SEEDS,
        hierarchical_bootstrap_resamples=source_report.hierarchical_bootstrap_resamples,
        cluster_sign_flip_resamples=cluster_sign_flip_resamples,
        significance_level=significance_level,
        sign_stability_threshold=sign_stability_threshold,
        seed_results=tuple(seed_results),
        aggregate_results=aggregates,
        decisions=decisions,
        diagnostics=diagnostics,
        overall_verdict=(
            "CONTINUE_RESEARCH"
            if promotion_allowed
            else "WARN_REVIEW_REQUIRED"
        ),
        promotion_allowed=promotion_allowed,
        biological_progression_allowed=False,
        interpretation_boundary=(
            "This is a synthetic independent-seed null-calibrated software analysis. "
            "Label-permutation and sequence-shuffle nulls remain separate and both must pass "
            "in the same direction. Even a passing result would not establish a biological "
            "mechanism, equal information content or DNA decryption. Motif-avoidance rejection "
            "sampling remains an unresolved generator limitation."
        ),
    )


def run_null_calibrated_analysis(
    *,
    sample_count: int = 400,
    sequence_length: int = 64,
    motif: str = "CGCG",
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
    projection_master_seed: int = DEFAULT_NULL_CALIBRATION_PROJECTION_MASTER_SEED,
    projection_seed_count: int = DEFAULT_NULL_CALIBRATION_PROJECTION_SEED_COUNT,
    data_master_seed: int = DEFAULT_NULL_CALIBRATION_DATA_MASTER_SEED,
    data_seed_count: int = DEFAULT_NULL_CALIBRATION_DATA_SEED_COUNT,
    bootstrap_resamples: int = 2_000,
    hierarchical_bootstrap_resamples: int = 5_000,
    cluster_sign_flip_resamples: int = DEFAULT_CLUSTER_SIGN_FLIP_RESAMPLES,
    significance_level: float = DEFAULT_SIGNIFICANCE_LEVEL,
    sign_stability_threshold: float = DEFAULT_SIGN_STABILITY_THRESHOLD,
) -> NullCalibratedReport:
    """Run the preregistered independent v0.7 schedule and calibrate both nulls."""

    projection_seeds = null_calibration_projection_seed_schedule(
        master_seed=projection_master_seed,
        count=projection_seed_count,
    )
    data_seeds = null_calibration_data_seed_schedule(
        master_seed=data_master_seed,
        count=data_seed_count,
    )
    source = run_projection_seed_sensitivity(
        SyntheticMotifConfig(
            sample_count=sample_count,
            sequence_length=sequence_length,
            motif=motif,
            seeds=data_seeds,
            bootstrap_resamples=bootstrap_resamples,
        ),
        budgets=budgets,
        projection_master_seed=projection_master_seed,
        projection_seed_count=projection_seed_count,
        projection_seeds=projection_seeds,
        hierarchical_bootstrap_resamples=hierarchical_bootstrap_resamples,
    )
    report = calibrate_projection_sensitivity_report(
        source,
        cluster_sign_flip_resamples=cluster_sign_flip_resamples,
        significance_level=significance_level,
        sign_stability_threshold=sign_stability_threshold,
    )
    return replace(report, data_master_seed=data_master_seed)


def _aggregate_calibrated_results(
    results: Sequence[NullCalibratedSeedResult],
    *,
    budgets: Sequence[int],
    projection_seeds: Sequence[int],
    data_seed_count: int,
    confidence: float,
    hierarchical_bootstrap_resamples: int,
    cluster_sign_flip_resamples: int,
    significance_level: float,
    sign_stability_threshold: float,
    seed: int,
) -> tuple[NullCalibratedAggregate, ...]:
    by_key = {
        (result.projection_seed, result.budget, result.null_condition): result
        for result in results
    }
    output: list[NullCalibratedAggregate] = []
    ordinal = 0
    for null_condition in _NULL_CONDITIONS:
        for budget in budgets:
            group = tuple(
                by_key[(projection_seed, budget, null_condition)]
                for projection_seed in projection_seeds
            )
            projection_means = tuple(
                result.mean_calibrated_difference for result in group
            )
            ci_low, ci_high = hierarchical_bootstrap_confidence_interval(
                tuple(
                    result.per_data_seed_calibrated_differences
                    for result in group
                ),
                confidence=confidence,
                resamples=hierarchical_bootstrap_resamples,
                seed=derived_seed(seed, ordinal),
            )
            p_value = cluster_sign_flip_p_value(
                projection_means,
                resamples=cluster_sign_flip_resamples,
                seed=derived_seed(seed, 10_000 + ordinal),
            )
            ordinal += 1
            positive_fraction = (
                sum(value > 0.0 for value in projection_means) / len(projection_means)
            )
            negative_fraction = (
                sum(value < 0.0 for value in projection_means) / len(projection_means)
            )
            interval_direction = (
                "tetranucleotide_advantage"
                if ci_low > 0.0
                else "orbital_advantage"
                if ci_high < 0.0
                else "unresolved"
            )
            sign_passed = (
                positive_fraction >= sign_stability_threshold
                if interval_direction == "tetranucleotide_advantage"
                else negative_fraction >= sign_stability_threshold
                if interval_direction == "orbital_advantage"
                else False
            )
            interval_passed = interval_direction != "unresolved"
            randomization_passed = p_value < significance_level
            gate_passed = interval_passed and sign_passed and randomization_passed
            output.append(
                NullCalibratedAggregate(
                    budget=budget,
                    null_condition=null_condition,
                    projection_seed_count=len(group),
                    data_seed_count=data_seed_count,
                    mean_calibrated_difference=(
                        math.fsum(projection_means) / len(projection_means)
                    ),
                    median_projection_mean_difference=statistics.median(
                        projection_means
                    ),
                    hierarchical_ci_low=ci_low,
                    hierarchical_ci_high=ci_high,
                    fraction_projection_means_positive=positive_fraction,
                    fraction_projection_means_negative=negative_fraction,
                    minimum_projection_mean_difference=min(projection_means),
                    maximum_projection_mean_difference=max(projection_means),
                    cluster_sign_flip_p_value=p_value,
                    direction=interval_direction,
                    interval_gate_passed=interval_passed,
                    sign_stability_gate_passed=sign_passed,
                    randomization_gate_passed=randomization_passed,
                    calibrated_gate_passed=gate_passed,
                )
            )
    return tuple(output)


def build_null_calibrated_decisions(
    aggregates: Sequence[NullCalibratedAggregate],
    *,
    budgets: Sequence[int],
) -> tuple[NullCalibratedDecision, ...]:
    by_key = {
        (result.budget, result.null_condition): result
        for result in aggregates
    }
    output: list[NullCalibratedDecision] = []
    for budget in budgets:
        label = by_key[(budget, "label_permutation_null")]
        shuffle = by_key[(budget, "sequence_shuffle_null")]
        same_direction = (
            label.calibrated_gate_passed
            and shuffle.calibrated_gate_passed
            and label.direction == shuffle.direction
            and label.direction != "unresolved"
        )
        output.append(
            NullCalibratedDecision(
                budget=budget,
                label_permutation_direction=label.direction,
                sequence_shuffle_direction=shuffle.direction,
                consensus_direction=label.direction if same_direction else "unresolved",
                both_nulls_passed=same_direction,
                promotion_gate_passed=same_direction,
            )
        )
    return tuple(output)
