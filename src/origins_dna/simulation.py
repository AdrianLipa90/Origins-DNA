"""Deterministic mutation simulation for controlled software experiments."""

from __future__ import annotations

from dataclasses import dataclass
import random

from .encoding import normalize_dna

_BASES = "ACGT"


@dataclass(frozen=True, slots=True)
class MutationRecord:
    """One substitution event."""

    position: int
    original: str
    mutated: str


def mutate_substitutions(
    sequence: str,
    rate: float,
    *,
    seed: int,
) -> tuple[str, tuple[MutationRecord, ...]]:
    """Apply independent substitutions with a deterministic pseudo-random seed."""

    if not isinstance(rate, (int, float)):
        raise TypeError("rate must be numeric")
    rate = float(rate)
    if not 0.0 <= rate <= 1.0:
        raise ValueError("rate must be between 0 and 1")

    normalized = normalize_dna(sequence)
    rng = random.Random(seed)
    output = list(normalized)
    records: list[MutationRecord] = []

    for position, original in enumerate(normalized):
        if rng.random() < rate:
            alternatives = [base for base in _BASES if base != original]
            mutated = rng.choice(alternatives)
            output[position] = mutated
            records.append(MutationRecord(position, original, mutated))

    return "".join(output), tuple(records)


def hamming_distance(left: str, right: str) -> int:
    """Count substitutions between equal-length normalized sequences."""

    left_normalized = normalize_dna(left)
    right_normalized = normalize_dna(right)
    if len(left_normalized) != len(right_normalized):
        raise ValueError("sequences must have equal length")
    return sum(a != b for a, b in zip(left_normalized, right_normalized))
