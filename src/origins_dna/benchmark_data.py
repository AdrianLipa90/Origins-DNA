"""Deterministic synthetic data and control transformations."""

from __future__ import annotations

import random
from typing import Sequence

from .benchmark_types import SyntheticDataset, SyntheticMotifConfig
from .encoding import normalize_dna

_BASES = "ACGT"


def generate_synthetic_motif_dataset(config: SyntheticMotifConfig, *, seed: int) -> SyntheticDataset:
    motif = normalize_dna(config.motif)
    decoy = _composition_matched_decoy(motif)
    rng = random.Random(seed)
    records: list[tuple[str, int]] = []
    for _ in range(config.sample_count // 2):
        records.append((_sequence_with_insert(config.sequence_length, motif, rng), 1))
        records.append((_negative_sequence(config.sequence_length, motif, decoy, rng), 0))
    rng.shuffle(records)
    sequences, labels = zip(*records)
    return SyntheticDataset(tuple(sequences), tuple(labels), motif, decoy)


def _composition_matched_decoy(motif: str) -> str:
    for shift in range(1, len(motif)):
        candidate = motif[shift:] + motif[:shift]
        if candidate != motif:
            return candidate
    raise ValueError("motif must admit a distinct cyclic decoy")


def _sequence_with_insert(length: int, insert: str, rng: random.Random) -> str:
    background = [rng.choice(_BASES) for _ in range(length)]
    start = rng.randrange(length - len(insert) + 1)
    background[start : start + len(insert)] = insert
    return "".join(background)


def _negative_sequence(length: int, motif: str, decoy: str, rng: random.Random) -> str:
    for _ in range(10_000):
        sequence = _sequence_with_insert(length, decoy, rng)
        if motif not in sequence:
            return sequence
    raise RuntimeError("failed to generate a motif-free negative sequence")


def shuffle_sequence(sequence: str, *, seed: int) -> str:
    bases = list(normalize_dna(sequence))
    random.Random(seed).shuffle(bases)
    return "".join(bases)


def permute_labels(labels: Sequence[int], *, seed: int) -> tuple[int, ...]:
    permuted = list(labels)
    random.Random(seed).shuffle(permuted)
    return tuple(permuted)


def stratified_split_indices(
    labels: Sequence[int], *, train_fraction: float, seed: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if not 0.0 < train_fraction < 1.0:
        raise ValueError("train_fraction must lie strictly between 0 and 1")
    rng = random.Random(seed)
    train: list[int] = []
    test: list[int] = []
    for class_label in (0, 1):
        indices = [index for index, label in enumerate(labels) if label == class_label]
        if len(indices) < 2:
            raise ValueError("each class requires at least two examples")
        rng.shuffle(indices)
        train_count = min(len(indices) - 1, max(1, round(len(indices) * train_fraction)))
        train.extend(indices[:train_count])
        test.extend(indices[train_count:])
    rng.shuffle(train)
    rng.shuffle(test)
    return tuple(train), tuple(test)


def derived_seed(seed: int, salt: int) -> int:
    return (seed * 1_000_003 + salt * 97_409) & 0xFFFFFFFF
