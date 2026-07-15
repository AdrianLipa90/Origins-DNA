"""Conventional sequence features and orbital ablations."""

from __future__ import annotations

from typing import Callable

from .benchmark_types import Vector
from .encoding import normalize_dna
from .features import base_frequencies, kmer_counts
from .orbital_atomic import orbital_atomic_summary

_BASES = "ACGT"
_DINUCLEOTIDES = tuple(left + right for left in _BASES for right in _BASES)
_ORBITAL_KEYS = (
    "mean_phase_cos",
    "mean_phase_sin",
    "phase_coherence",
    "mean_ring_class",
    "mean_bond_class",
    "signed_phase_flow_per_step",
)
Representation = Callable[[str], Vector]


def representation_registry(sequence_length: int) -> dict[str, Representation]:
    if sequence_length <= 0:
        raise ValueError("sequence_length must be positive")
    return {
        "base_frequency": base_frequency_vector,
        "dinucleotide_frequency": dinucleotide_frequency_vector,
        "positional_one_hot": lambda sequence: positional_one_hot_vector(sequence, sequence_length),
        "orbital_full": orbital_full_vector,
        "orbital_phase_only": orbital_phase_vector,
        "orbital_classes_only": orbital_classes_vector,
    }


def base_frequency_vector(sequence: str) -> Vector:
    frequencies = base_frequencies(sequence)
    return tuple(frequencies[base] for base in _BASES)


def dinucleotide_frequency_vector(sequence: str) -> Vector:
    normalized = normalize_dna(sequence)
    counts = kmer_counts(normalized, 2)
    denominator = max(1, len(normalized) - 1)
    return tuple(counts.get(kmer, 0) / denominator for kmer in _DINUCLEOTIDES)


def positional_one_hot_vector(sequence: str, expected_length: int) -> Vector:
    normalized = normalize_dna(sequence)
    if len(normalized) != expected_length:
        raise ValueError(
            f"sequence length {len(normalized)} does not match expected_length {expected_length}"
        )
    return tuple(
        1.0 if observed == candidate else 0.0
        for observed in normalized
        for candidate in _BASES
    )


def orbital_full_vector(sequence: str) -> Vector:
    summary = orbital_atomic_summary(sequence)
    return tuple(summary[key] for key in _ORBITAL_KEYS)


def orbital_phase_vector(sequence: str) -> Vector:
    summary = orbital_atomic_summary(sequence)
    keys = ("mean_phase_cos", "mean_phase_sin", "phase_coherence", "signed_phase_flow_per_step")
    return tuple(summary[key] for key in keys)


def orbital_classes_vector(sequence: str) -> Vector:
    summary = orbital_atomic_summary(sequence)
    return (summary["mean_ring_class"], summary["mean_bond_class"])
