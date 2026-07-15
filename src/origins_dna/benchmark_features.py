"""Conventional sequence features and orbital ablations."""

from __future__ import annotations

from itertools import product
from typing import Callable, Sequence

from .benchmark_types import Vector
from .encoding import normalize_dna
from .features import base_frequencies, kmer_counts
from .orbital_atomic import orbital_atomic_encode, orbital_atomic_summary

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
    """Return the frozen v0.3 representation registry."""

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


def local_order_representation_registry(sequence_length: int) -> dict[str, Representation]:
    """Return the v0.4 registry for diagnosing local-order information loss."""

    if sequence_length <= 3:
        raise ValueError("sequence_length must exceed the maximum lag of 3")
    return {
        "dinucleotide_frequency": dinucleotide_frequency_vector,
        "trinucleotide_frequency": lambda sequence: kmer_frequency_vector(sequence, 3),
        "tetranucleotide_frequency": lambda sequence: kmer_frequency_vector(sequence, 4),
        "orbital_global_summary": orbital_full_vector,
        "orbital_pair_tensor": orbital_pair_tensor_vector,
        "orbital_multilag_tensor": orbital_multilag_tensor_vector,
        "orbital_multilag_phase_only": orbital_multilag_phase_vector,
        "orbital_multilag_classes_only": orbital_multilag_classes_vector,
    }


def base_frequency_vector(sequence: str) -> Vector:
    frequencies = base_frequencies(sequence)
    return tuple(frequencies[base] for base in _BASES)


def kmer_frequency_vector(sequence: str, k: int) -> Vector:
    """Return normalized frequencies for all lexicographic A/C/G/T k-mers."""

    if not isinstance(k, int):
        raise TypeError("k must be an integer")
    if k <= 0:
        raise ValueError("k must be positive")
    normalized = normalize_dna(sequence)
    vocabulary = tuple("".join(symbols) for symbols in product(_BASES, repeat=k))
    counts = kmer_counts(normalized, k)
    denominator = max(1, len(normalized) - k + 1)
    return tuple(counts.get(kmer, 0) / denominator for kmer in vocabulary)


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


def orbital_lag_tensor_vector(
    sequence: str,
    *,
    lags: Sequence[int],
    channels: Sequence[int] = (0, 1, 2, 3),
) -> Vector:
    """Average ordered outer products of orbital channels at declared lags.

    The vector preserves directional local order without using task labels or the
    planted motif. Channels follow ``OrbitalAtomicState.as_tuple``.
    """

    normalized = normalize_dna(sequence)
    materialized_lags = tuple(lags)
    materialized_channels = tuple(channels)
    if not materialized_lags:
        raise ValueError("lags must not be empty")
    if any(not isinstance(lag, int) for lag in materialized_lags):
        raise TypeError("lags must contain integers")
    if any(lag <= 0 for lag in materialized_lags):
        raise ValueError("lags must be positive")
    if len(set(materialized_lags)) != len(materialized_lags):
        raise ValueError("lags must be unique")
    if not materialized_channels:
        raise ValueError("channels must not be empty")
    if any(not isinstance(channel, int) for channel in materialized_channels):
        raise TypeError("channels must contain integers")
    if any(channel < 0 or channel > 3 for channel in materialized_channels):
        raise ValueError("channels must lie between 0 and 3")
    if len(set(materialized_channels)) != len(materialized_channels):
        raise ValueError("channels must be unique")
    if normalized and max(materialized_lags) >= len(normalized):
        raise ValueError("each lag must be smaller than the sequence length")
    if not normalized:
        return (0.0,) * (
            len(materialized_lags) * len(materialized_channels) * len(materialized_channels)
        )

    states = tuple(state.as_tuple() for state in orbital_atomic_encode(normalized))
    output: list[float] = []
    for lag in materialized_lags:
        pair_count = len(states) - lag
        for left_channel in materialized_channels:
            for right_channel in materialized_channels:
                output.append(
                    sum(
                        states[index][left_channel] * states[index + lag][right_channel]
                        for index in range(pair_count)
                    )
                    / pair_count
                )
    return tuple(output)


def orbital_pair_tensor_vector(sequence: str) -> Vector:
    return orbital_lag_tensor_vector(sequence, lags=(1,))


def orbital_multilag_tensor_vector(sequence: str) -> Vector:
    return orbital_lag_tensor_vector(sequence, lags=(1, 2, 3))


def orbital_multilag_phase_vector(sequence: str) -> Vector:
    return orbital_lag_tensor_vector(sequence, lags=(1, 2, 3), channels=(0, 1))


def orbital_multilag_classes_vector(sequence: str) -> Vector:
    return orbital_lag_tensor_vector(sequence, lags=(1, 2, 3), channels=(2, 3))
