"""Public facade for controlled Origins-DNA benchmarks."""

from .benchmark_data import (
    generate_synthetic_motif_dataset,
    permute_labels,
    shuffle_sequence,
    stratified_split_indices,
)
from .benchmark_features import (
    base_frequency_vector,
    dinucleotide_frequency_vector,
    orbital_classes_vector,
    orbital_full_vector,
    orbital_phase_vector,
    positional_one_hot_vector,
    representation_registry,
)
from .benchmark_model import confusion_matrix, fit_nearest_centroid
from .benchmark_runner import bootstrap_mean_confidence_interval, run_controlled_benchmark
from .benchmark_types import (
    AggregateResult,
    BenchmarkReport,
    ConfusionMatrix,
    NearestCentroidModel,
    SeedResult,
    SyntheticDataset,
    SyntheticMotifConfig,
)

__all__ = [
    "AggregateResult",
    "BenchmarkReport",
    "ConfusionMatrix",
    "NearestCentroidModel",
    "SeedResult",
    "SyntheticDataset",
    "SyntheticMotifConfig",
    "base_frequency_vector",
    "bootstrap_mean_confidence_interval",
    "confusion_matrix",
    "dinucleotide_frequency_vector",
    "fit_nearest_centroid",
    "generate_synthetic_motif_dataset",
    "orbital_classes_vector",
    "orbital_full_vector",
    "orbital_phase_vector",
    "permute_labels",
    "positional_one_hot_vector",
    "representation_registry",
    "run_controlled_benchmark",
    "shuffle_sequence",
    "stratified_split_indices",
]
