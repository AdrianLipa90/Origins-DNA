"""Paired-background v0.9.1 benchmark for the three-quarter causal readout."""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
import math
from pathlib import Path
import random
from typing import Callable, Sequence

from .benchmark_data import derived_seed, permute_labels, shuffle_sequence
from .benchmark_features import (
    Representation,
    kmer_frequency_vector,
    orbital_multilag_tensor_vector,
)
from .benchmark_model import confusion_matrix, fit_nearest_centroid
from .benchmark_runner import bootstrap_mean_confidence_interval
from .benchmark_types import (
    AggregateResult,
    BenchmarkReport,
    SeedResult,
    SyntheticMotifConfig,
)
from .capacity_matched import DEFAULT_DIMENSION_BUDGETS, projected_representation
from .causal_dhse import causal_dhse_vector
from .causal_three_quarter import (
    even_three_quarter_vector,
    exact_three_quarter_gate_receipt,
    phase_ordinary_mean_vector,
    phase_three_quarter_vector,
)
from .encoding import normalize_dna

_BASES = "ACGT"


@dataclass(frozen=True, slots=True)
class PairedSyntheticDataset:
    sequences: tuple[str, ...]
    labels: tuple[int, ...]
    pair_ids: tuple[int, ...]
    motif: str
    decoy: str

    def __post_init__(self) -> None:
        if not (len(self.sequences) == len(self.labels) == len(self.pair_ids)):
            raise ValueError("paired dataset fields must have equal length")
        if set(self.labels) != {0, 1}:
            raise ValueError("paired dataset must contain both classes")
        for pair_id in set(self.pair_ids):
            indices = [i for i, value in enumerate(self.pair_ids) if value == pair_id]
            if len(indices) != 2 or {self.labels[i] for i in indices} != {0, 1}:
                raise ValueError("each pair must contain one sequence per class")


def _cyclic_decoy(motif: str) -> str:
    normalized = normalize_dna(motif)
    for shift in range(1, len(normalized)):
        candidate = normalized[shift:] + normalized[:shift]
        if candidate != normalized:
            return candidate
    raise ValueError("motif must admit a distinct cyclic decoy")


def _occurrences(sequence: str, pattern: str) -> tuple[int, ...]:
    return tuple(
        index
        for index in range(0, len(sequence) - len(pattern) + 1)
        if sequence[index : index + len(pattern)] == pattern
    )


def _insert(background: str, start: int, word: str) -> str:
    return background[:start] + word + background[start + len(word) :]


def generate_paired_motif_dataset(
    config: SyntheticMotifConfig, *, seed: int
) -> PairedSyntheticDataset:
    motif = normalize_dna(config.motif)
    decoy = _cyclic_decoy(motif)
    rng = random.Random(seed)
    pair_count = config.sample_count // 2
    records: list[tuple[str, int, int]] = []
    for pair_id in range(pair_count):
        for _ in range(100_000):
            background = "".join(rng.choice(_BASES) for _ in range(config.sequence_length))
            start = rng.randrange(config.sequence_length - len(motif) + 1)
            positive = _insert(background, start, motif)
            negative = _insert(background, start, decoy)
            if (
                _occurrences(positive, motif) == (start,)
                and not _occurrences(positive, decoy)
                and _occurrences(negative, decoy) == (start,)
                and not _occurrences(negative, motif)
            ):
                records.append((positive, 1, pair_id))
                records.append((negative, 0, pair_id))
                break
        else:
            raise RuntimeError("failed to generate an uncontaminated paired background")
    rng.shuffle(records)
    sequences, labels, pair_ids = zip(*records)
    return PairedSyntheticDataset(
        tuple(sequences), tuple(labels), tuple(pair_ids), motif, decoy
    )


def paired_split_indices(
    pair_ids: Sequence[int], *, train_fraction: float, seed: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    unique = sorted(set(int(value) for value in pair_ids))
    if len(unique) < 2:
        raise ValueError("at least two pairs are required")
    if not 0.0 < train_fraction < 1.0:
        raise ValueError("train_fraction must lie strictly between zero and one")
    rng = random.Random(seed)
    rng.shuffle(unique)
    train_count = min(len(unique) - 1, max(1, round(len(unique) * train_fraction)))
    train_pairs = set(unique[:train_count])
    train = tuple(index for index, pair_id in enumerate(pair_ids) if pair_id in train_pairs)
    test = tuple(index for index, pair_id in enumerate(pair_ids) if pair_id not in train_pairs)
    if set(pair_ids[index] for index in train) & set(pair_ids[index] for index in test):
        raise AssertionError("pair leakage detected")
    return train, test


def _dimension(representation: Representation, sequence_length: int) -> int:
    probe = ("ACGT" * ((sequence_length + 3) // 4))[:sequence_length]
    return len(representation(probe))


def _concatenated_projection(
    left: Representation,
    right: Representation,
    *,
    left_dimension: int,
    right_dimension: int,
    output_dimension: int,
    seed: int,
) -> Representation:
    if output_dimension % 2:
        raise ValueError("split output dimension must be even")
    half = output_dimension // 2
    left_projected = projected_representation(
        left,
        input_dimension=left_dimension,
        output_dimension=half,
        seed=derived_seed(seed, 1),
    )
    right_projected = projected_representation(
        right,
        input_dimension=right_dimension,
        output_dimension=half,
        seed=derived_seed(seed, 2),
    )

    @lru_cache(maxsize=None)
    def transform(sequence: str) -> tuple[float, ...]:
        return left_projected(sequence) + right_projected(sequence)

    return transform


def three_quarter_representation_registry(
    sequence_length: int,
    *,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
    projection_seed: int = 134_091,
) -> dict[str, Representation]:
    @lru_cache(maxsize=None)
    def v09(sequence: str) -> tuple[float, ...]:
        return causal_dhse_vector(sequence)

    @lru_cache(maxsize=None)
    def phase_34(sequence: str) -> tuple[float, ...]:
        return phase_three_quarter_vector(sequence)

    @lru_cache(maxsize=None)
    def even_34(sequence: str) -> tuple[float, ...]:
        return even_three_quarter_vector(sequence)

    @lru_cache(maxsize=None)
    def phase_mean(sequence: str) -> tuple[float, ...]:
        return phase_ordinary_mean_vector(sequence)

    @lru_cache(maxsize=None)
    def tetranucleotide(sequence: str) -> tuple[float, ...]:
        return kmer_frequency_vector(sequence, 4)

    @lru_cache(maxsize=None)
    def orbital(sequence: str) -> tuple[float, ...]:
        return orbital_multilag_tensor_vector(sequence)

    sources: dict[str, Representation] = {
        "causal_v09": v09,
        "phase_34": phase_34,
        "even_34": even_34,
        "phase_mean": phase_mean,
        "tetranucleotide": tetranucleotide,
        "orbital": orbital,
    }
    dimensions = {
        name: _dimension(representation, sequence_length)
        for name, representation in sources.items()
    }
    registry: dict[str, Representation] = {}
    for budget_ordinal, budget in enumerate(sorted(tuple(budgets))):
        if budget <= 0 or budget % 2:
            raise ValueError("budgets must be positive even integers")
        family_seed = derived_seed(projection_seed, budget_ordinal)
        projection_salts = {
            "causal_v09": 1,
            "phase_34": 2,
            "phase_mean": 2,
            "even_34": 3,
            "tetranucleotide": 4,
            "orbital": 5,
        }
        if dimensions["phase_34"] != dimensions["phase_mean"]:
            raise AssertionError("phase exponent ablations must have equal source dimension")
        for name, representation in sources.items():
            registry[f"{name}_projected_{budget}"] = projected_representation(
                representation,
                input_dimension=dimensions[name],
                output_dimension=budget,
                seed=derived_seed(family_seed, projection_salts[name]),
            )
        registry[f"split_phase_dhse_34_projected_{budget}"] = _concatenated_projection(
            phase_34,
            even_34,
            left_dimension=dimensions["phase_34"],
            right_dimension=dimensions["even_34"],
            output_dimension=budget,
            seed=derived_seed(family_seed, 101),
        )
    return registry


def _evaluate(
    dataset: PairedSyntheticDataset,
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
            tuple(vectors[index] for index in train_indices), train_labels
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


def _transformed_dataset(
    dataset: PairedSyntheticDataset,
    transform: Callable[[str, int], str],
) -> PairedSyntheticDataset:
    sequences = tuple(
        transform(sequence, index) for index, sequence in enumerate(dataset.sequences)
    )
    return PairedSyntheticDataset(
        sequences, dataset.labels, dataset.pair_ids, dataset.motif, dataset.decoy
    )


def _aggregate_results(
    results: Sequence[SeedResult], config: SyntheticMotifConfig
) -> tuple[AggregateResult, ...]:
    groups: dict[tuple[str, str], list[SeedResult]] = {}
    for result in results:
        groups.setdefault((result.representation, result.condition), []).append(result)
    output: list[AggregateResult] = []
    for ordinal, ((representation, condition), group) in enumerate(sorted(groups.items())):
        balanced = tuple(item.balanced_accuracy for item in group)
        mcc = tuple(item.matthews_correlation for item in group)
        ci_low, ci_high = bootstrap_mean_confidence_interval(
            balanced,
            confidence=config.confidence,
            resamples=config.bootstrap_resamples,
            seed=derived_seed(34_091, ordinal),
        )
        output.append(
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
    return tuple(output)


def run_three_quarter_benchmark(
    config: SyntheticMotifConfig,
    *,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
    projection_seed: int = 134_091,
) -> BenchmarkReport:
    receipt = exact_three_quarter_gate_receipt(config.motif, _cyclic_decoy(config.motif))
    if receipt["technical_status"] != "PASS":
        raise AssertionError("three-quarter exact gates failed")
    registry = three_quarter_representation_registry(
        config.sequence_length,
        budgets=budgets,
        projection_seed=projection_seed,
    )
    results: list[SeedResult] = []
    for seed in config.seeds:
        dataset = generate_paired_motif_dataset(config, seed=seed)
        train_indices, test_indices = paired_split_indices(
            dataset.pair_ids,
            train_fraction=config.train_fraction,
            seed=derived_seed(seed, 1),
        )
        results.extend(
            _evaluate(
                dataset,
                registry,
                train_indices,
                test_indices,
                seed=seed,
                condition="observed_labels",
            )
        )
        results.extend(
            _evaluate(
                dataset,
                registry,
                train_indices,
                test_indices,
                seed=seed,
                condition="label_permutation_null",
                permute_training_labels=True,
            )
        )
        shuffled = _transformed_dataset(
            dataset,
            lambda sequence, index: shuffle_sequence(
                sequence, seed=derived_seed(seed, 10_000 + index)
            ),
        )
        results.extend(
            _evaluate(
                shuffled,
                registry,
                train_indices,
                test_indices,
                seed=seed,
                condition="sequence_shuffle_null",
            )
        )
        reversed_dataset = _transformed_dataset(
            dataset, lambda sequence, _index: sequence[::-1]
        )
        results.extend(
            _evaluate(
                reversed_dataset,
                registry,
                train_indices,
                test_indices,
                seed=seed,
                condition="sequence_reversal",
            )
        )
    return BenchmarkReport(
        schema="origins-dna.causal-three-quarter-benchmark.v1",
        benchmark_status="synthetic_paired_three_quarter_result",
        config=config,
        representations=tuple(registry),
        seed_results=tuple(results),
        aggregate_results=_aggregate_results(results, config),
        interpretation_boundary=(
            "This paired synthetic benchmark tests a fixed three-quarter-scaled "
            "projective readout. It does not establish biological causation, "
            "abiogenesis, DNA decryption or literal Möbius dynamics in DNA."
        ),
    )


def _seed_budget_means(
    report: BenchmarkReport,
    family: str,
    *,
    condition: str,
    budgets: Sequence[int],
) -> dict[int, float]:
    by_key = {
        (row.seed, row.representation, row.condition): row.balanced_accuracy
        for row in report.seed_results
    }
    output: dict[int, float] = {}
    for seed in report.config.seeds:
        values = tuple(
            by_key[(seed, f"{family}_projected_{budget}", condition)]
            for budget in budgets
        )
        output[seed] = math.fsum(values) / len(values)
    return output


def _paired_contrast(
    report: BenchmarkReport,
    *,
    left_family: str,
    right_family: str,
    left_condition: str = "observed_labels",
    right_condition: str = "observed_labels",
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
    seed: int,
) -> dict[str, object]:
    left = _seed_budget_means(
        report, left_family, condition=left_condition, budgets=budgets
    )
    right = _seed_budget_means(
        report, right_family, condition=right_condition, budgets=budgets
    )
    differences = tuple(left[item] - right[item] for item in report.config.seeds)
    ci_low, ci_high = bootstrap_mean_confidence_interval(
        differences,
        confidence=report.config.confidence,
        resamples=report.config.bootstrap_resamples,
        seed=seed,
    )
    return {
        "left_family": left_family,
        "right_family": right_family,
        "left_condition": left_condition,
        "right_condition": right_condition,
        "budgets": list(budgets),
        "seed_count": len(differences),
        "mean_balanced_accuracy_difference": math.fsum(differences) / len(differences),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "per_seed_differences": list(differences),
        "passed": ci_low > 0.0,
    }


def three_quarter_comparison_summary(
    report: BenchmarkReport,
    *,
    budgets: Sequence[int] = DEFAULT_DIMENSION_BUDGETS,
) -> dict[str, object]:
    contrasts = (
        _paired_contrast(
            report,
            left_family="phase_34",
            right_family="even_34",
            budgets=budgets,
            seed=134_501,
        ),
        _paired_contrast(
            report,
            left_family="phase_34",
            right_family="phase_mean",
            budgets=budgets,
            seed=134_502,
        ),
        _paired_contrast(
            report,
            left_family="phase_34",
            right_family="phase_34",
            right_condition="sequence_shuffle_null",
            budgets=budgets,
            seed=134_503,
        ),
        _paired_contrast(
            report,
            left_family="split_phase_dhse_34",
            right_family="causal_v09",
            budgets=budgets,
            seed=134_504,
        ),
    )
    primary = contrasts[:3]
    return {
        "schema": "origins-dna.causal-three-quarter-comparisons.v1",
        "source_schema": report.schema,
        "primary_budget_aggregation": "per-seed arithmetic mean over budgets 12,24,48",
        "contrasts": list(contrasts),
        "three_quarter_phase_signal_status": (
            "PASS" if all(bool(row["passed"]) for row in primary) else "FAIL"
        ),
    }


def write_three_quarter_bundle(
    report: BenchmarkReport, output: str | Path
) -> tuple[Path, Path, Path, Path, Path]:
    full_json, archive, summary = report.write_bundle(output)
    comparisons = full_json.with_name(full_json.stem + "_comparisons.json")
    comparisons.write_text(
        json.dumps(
            three_quarter_comparison_summary(report),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    gates = full_json.with_name(full_json.stem + "_gates.json")
    gates.write_text(
        json.dumps(exact_three_quarter_gate_receipt(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return full_json, archive, summary, comparisons, gates
