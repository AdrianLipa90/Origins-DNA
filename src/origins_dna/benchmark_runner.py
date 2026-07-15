"""Controlled benchmark execution, aggregation and confidence intervals."""

from __future__ import annotations

import math
import random
from typing import Sequence

from .benchmark_data import (
    derived_seed,
    generate_synthetic_motif_dataset,
    permute_labels,
    shuffle_sequence,
    stratified_split_indices,
)
from .benchmark_features import Representation, representation_registry
from .benchmark_model import confusion_matrix, fit_nearest_centroid
from .benchmark_types import (
    AggregateResult,
    BenchmarkReport,
    SeedResult,
    SyntheticDataset,
    SyntheticMotifConfig,
)


def run_controlled_benchmark(config: SyntheticMotifConfig) -> BenchmarkReport:
    return run_registered_benchmark(
        config,
        representation_registry(config.sequence_length),
        schema="origins-dna.benchmark.v1",
        benchmark_status="synthetic_controlled_software_result",
        interpretation_boundary=(
            "Results quantify a synthetic software task only. They neither validate a biological "
            "mechanism nor establish that DNA has been decrypted."
        ),
    )


def run_registered_benchmark(
    config: SyntheticMotifConfig,
    registry: dict[str, Representation],
    *,
    schema: str,
    benchmark_status: str,
    interpretation_boundary: str,
) -> BenchmarkReport:
    """Execute one deterministic registry under shared splits and controls."""

    if not registry:
        raise ValueError("registry must not be empty")
    if not schema.strip():
        raise ValueError("schema must not be empty")
    if not benchmark_status.strip():
        raise ValueError("benchmark_status must not be empty")
    if not interpretation_boundary.strip():
        raise ValueError("interpretation_boundary must not be empty")
    results: list[SeedResult] = []
    for seed in config.seeds:
        dataset = generate_synthetic_motif_dataset(config, seed=seed)
        train_indices, test_indices = stratified_split_indices(
            dataset.labels,
            train_fraction=config.train_fraction,
            seed=derived_seed(seed, 1),
        )
        results.extend(
            _evaluate_condition(
                dataset,
                registry,
                train_indices,
                test_indices,
                seed=seed,
                condition="observed_labels",
            )
        )
        results.extend(
            _evaluate_condition(
                dataset,
                registry,
                train_indices,
                test_indices,
                seed=seed,
                condition="label_permutation_null",
                permute_training_labels=True,
            )
        )
        shuffled = SyntheticDataset(
            sequences=tuple(
                shuffle_sequence(sequence, seed=derived_seed(seed, 10_000 + index))
                for index, sequence in enumerate(dataset.sequences)
            ),
            labels=dataset.labels,
            motif=dataset.motif,
            decoy=dataset.decoy,
        )
        results.extend(
            _evaluate_condition(
                shuffled,
                registry,
                train_indices,
                test_indices,
                seed=seed,
                condition="sequence_shuffle_null",
            )
        )
    return BenchmarkReport(
        schema=schema,
        benchmark_status=benchmark_status,
        config=config,
        representations=tuple(registry),
        seed_results=tuple(results),
        aggregate_results=_aggregate_results(results, config),
        interpretation_boundary=interpretation_boundary,
    )


def _evaluate_condition(
    dataset: SyntheticDataset,
    registry: dict[str, Representation],
    train_indices: Sequence[int],
    test_indices: Sequence[int],
    *,
    seed: int,
    condition: str,
    permute_training_labels: bool = False,
) -> tuple[SeedResult, ...]:
    train_labels = tuple(dataset.labels[index] for index in train_indices)
    if permute_training_labels:
        train_labels = permute_labels(train_labels, seed=derived_seed(seed, 2))
    test_labels = tuple(dataset.labels[index] for index in test_indices)
    output: list[SeedResult] = []
    for name, representation in registry.items():
        vectors = tuple(representation(sequence) for sequence in dataset.sequences)
        model = fit_nearest_centroid(
            tuple(vectors[index] for index in train_indices),
            train_labels,
        )
        predictions = tuple(model.predict(vectors[index]) for index in test_indices)
        metrics = confusion_matrix(test_labels, predictions)
        output.append(
            SeedResult(
                seed=seed,
                representation=name,
                condition=condition,
                train_size=len(train_indices),
                test_size=len(test_indices),
                accuracy=metrics.accuracy,
                balanced_accuracy=metrics.balanced_accuracy,
                matthews_correlation=metrics.matthews_correlation,
            )
        )
    return tuple(output)


def _aggregate_results(
    results: Sequence[SeedResult], config: SyntheticMotifConfig
) -> tuple[AggregateResult, ...]:
    groups: dict[tuple[str, str], list[SeedResult]] = {}
    for result in results:
        groups.setdefault((result.representation, result.condition), []).append(result)
    aggregates: list[AggregateResult] = []
    for ordinal, ((representation, condition), group) in enumerate(sorted(groups.items())):
        balanced = tuple(item.balanced_accuracy for item in group)
        mcc = tuple(item.matthews_correlation for item in group)
        ci_low, ci_high = bootstrap_mean_confidence_interval(
            balanced,
            confidence=config.confidence,
            resamples=config.bootstrap_resamples,
            seed=derived_seed(17_071, ordinal),
        )
        aggregates.append(
            AggregateResult(
                representation=representation,
                condition=condition,
                seed_count=len(group),
                mean_balanced_accuracy=math.fsum(balanced) / len(balanced),
                ci_low=ci_low,
                ci_high=ci_high,
                mean_matthews_correlation=math.fsum(mcc) / len(mcc),
            )
        )
    return tuple(aggregates)


def bootstrap_mean_confidence_interval(
    values: Sequence[float],
    *,
    confidence: float,
    resamples: int,
    seed: int,
) -> tuple[float, float]:
    if not values:
        raise ValueError("values must not be empty")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie strictly between 0 and 1")
    if resamples <= 0:
        raise ValueError("resamples must be positive")
    materialized = tuple(float(value) for value in values)
    if any(not math.isfinite(value) for value in materialized):
        raise ValueError("values must be finite")
    if len(materialized) == 1:
        return materialized[0], materialized[0]
    rng = random.Random(seed)
    sample_size = len(materialized)
    means = sorted(
        math.fsum(materialized[rng.randrange(sample_size)] for _ in range(sample_size))
        / sample_size
        for _ in range(resamples)
    )
    tail = (1.0 - confidence) / 2.0
    low_index = min(resamples - 1, max(0, math.floor(tail * resamples)))
    high_index = min(resamples - 1, max(0, math.ceil((1.0 - tail) * resamples) - 1))
    return means[low_index], means[high_index]
