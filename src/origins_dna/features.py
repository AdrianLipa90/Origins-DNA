"""Small, dependency-free sequence feature baselines."""

from __future__ import annotations

from collections import Counter
import math

from .encoding import normalize_dna


def base_frequencies(sequence: str) -> dict[str, float]:
    """Return normalized A/C/G/T frequencies.

    The empty sequence maps every base to 0.0.
    """

    normalized = normalize_dna(sequence)
    if not normalized:
        return {base: 0.0 for base in "ACGT"}

    counts = Counter(normalized)
    length = len(normalized)
    return {base: counts.get(base, 0) / length for base in "ACGT"}


def kmer_counts(sequence: str, k: int) -> dict[str, int]:
    """Count overlapping k-mers in lexicographic insertion order."""

    if not isinstance(k, int):
        raise TypeError("k must be an integer")
    if k <= 0:
        raise ValueError("k must be positive")

    normalized = normalize_dna(sequence)
    counts: dict[str, int] = {}
    for start in range(max(0, len(normalized) - k + 1)):
        kmer = normalized[start : start + k]
        counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def shannon_entropy_bits(sequence: str) -> float:
    """Return zero-order Shannon entropy in bits per nucleotide."""

    frequencies = base_frequencies(sequence)
    return -sum(p * math.log2(p) for p in frequencies.values() if p > 0.0)
