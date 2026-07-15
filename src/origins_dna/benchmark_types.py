"""Typed records and deterministic report serialization for benchmarks."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import gzip
import json
from pathlib import Path
from typing import Sequence

from .encoding import normalize_dna

Vector = tuple[float, ...]


@dataclass(frozen=True, slots=True)
class SyntheticMotifConfig:
    """Configuration for a balanced composition-matched motif task."""

    sample_count: int = 400
    sequence_length: int = 64
    motif: str = "CGCG"
    train_fraction: float = 0.7
    seeds: tuple[int, ...] = (11, 23, 37, 41, 53, 67, 79, 83)
    bootstrap_resamples: int = 2000
    confidence: float = 0.95

    def __post_init__(self) -> None:
        if self.sample_count < 4 or self.sample_count % 2:
            raise ValueError("sample_count must be even and at least 4")
        if self.sequence_length <= 0:
            raise ValueError("sequence_length must be positive")
        motif = normalize_dna(self.motif)
        if not motif:
            raise ValueError("motif must not be empty")
        if len(motif) > self.sequence_length:
            raise ValueError("motif must fit inside sequence_length")
        if not 0.0 < self.train_fraction < 1.0:
            raise ValueError("train_fraction must lie strictly between 0 and 1")
        if not self.seeds:
            raise ValueError("at least one seed is required")
        if len(set(self.seeds)) != len(self.seeds):
            raise ValueError("seeds must be unique")
        if self.bootstrap_resamples <= 0:
            raise ValueError("bootstrap_resamples must be positive")
        if not 0.0 < self.confidence < 1.0:
            raise ValueError("confidence must lie strictly between 0 and 1")


@dataclass(frozen=True, slots=True)
class SyntheticDataset:
    """Balanced sequence dataset and binary labels."""

    sequences: tuple[str, ...]
    labels: tuple[int, ...]
    motif: str
    decoy: str

    def __post_init__(self) -> None:
        if len(self.sequences) != len(self.labels):
            raise ValueError("sequences and labels must have equal length")
        if set(self.labels) != {0, 1}:
            raise ValueError("dataset must contain both binary classes")


@dataclass(frozen=True, slots=True)
class ConfusionMatrix:
    """Binary confusion matrix with derived metrics."""

    true_negative: int
    false_positive: int
    false_negative: int
    true_positive: int

    @property
    def total(self) -> int:
        return self.true_negative + self.false_positive + self.false_negative + self.true_positive

    @property
    def accuracy(self) -> float:
        return (self.true_negative + self.true_positive) / self.total

    @property
    def balanced_accuracy(self) -> float:
        specificity = self.true_negative / (self.true_negative + self.false_positive)
        sensitivity = self.true_positive / (self.true_positive + self.false_negative)
        return 0.5 * (specificity + sensitivity)

    @property
    def matthews_correlation(self) -> float:
        import math

        numerator = self.true_positive * self.true_negative - self.false_positive * self.false_negative
        denominator = math.sqrt(
            (self.true_positive + self.false_positive)
            * (self.true_positive + self.false_negative)
            * (self.true_negative + self.false_positive)
            * (self.true_negative + self.false_negative)
        )
        return 0.0 if denominator == 0.0 else numerator / denominator


@dataclass(frozen=True, slots=True)
class SeedResult:
    seed: int
    representation: str
    condition: str
    train_size: int
    test_size: int
    accuracy: float
    balanced_accuracy: float
    matthews_correlation: float


@dataclass(frozen=True, slots=True)
class AggregateResult:
    representation: str
    condition: str
    seed_count: int
    mean_balanced_accuracy: float
    ci_low: float
    ci_high: float
    mean_matthews_correlation: float


@dataclass(frozen=True, slots=True)
class BenchmarkReport:
    schema: str
    benchmark_status: str
    config: SyntheticMotifConfig
    representations: tuple[str, ...]
    seed_results: tuple[SeedResult, ...]
    aggregate_results: tuple[AggregateResult, ...]
    interpretation_boundary: str

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": self.schema,
            "benchmark_status": self.benchmark_status,
            "config": asdict(self.config),
            "representations": list(self.representations),
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

    def write_json(self, path: str | Path) -> None:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    def write_bundle(self, path: str | Path) -> tuple[Path, Path, Path]:
        output = Path(path)
        self.write_json(output)
        archive = output.with_name(output.name + ".gz")
        with archive.open("wb") as raw_file:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_file, mtime=0) as compressed:
                compressed.write(output.read_bytes())
        summary = output.with_name(output.stem + "_summary.json")
        summary.write_text(
            json.dumps(self.summary_dict(full_result_archive=archive.name), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return output, archive, summary


@dataclass(frozen=True, slots=True)
class NearestCentroidModel:
    means: Vector
    scales: Vector
    negative_centroid: Vector
    positive_centroid: Vector

    def predict(self, vector: Sequence[float]) -> int:
        from .benchmark_model import predict_nearest_centroid

        return predict_nearest_centroid(self, vector)
