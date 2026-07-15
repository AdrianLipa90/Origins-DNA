import pytest

from origins_dna import hamming_distance, mutate_substitutions


def test_mutation_is_deterministic_for_seed() -> None:
    first = mutate_substitutions("ACGTACGT", 0.5, seed=7)
    second = mutate_substitutions("ACGTACGT", 0.5, seed=7)
    assert first == second


def test_zero_rate_preserves_sequence() -> None:
    mutated, records = mutate_substitutions("ACGT", 0.0, seed=1)
    assert mutated == "ACGT"
    assert records == ()


def test_full_rate_mutates_every_position() -> None:
    mutated, records = mutate_substitutions("ACGT", 1.0, seed=1)
    assert hamming_distance("ACGT", mutated) == 4
    assert len(records) == 4


def test_invalid_rate_is_rejected() -> None:
    with pytest.raises(ValueError):
        mutate_substitutions("ACGT", 1.01, seed=1)


def test_hamming_distance_requires_equal_lengths() -> None:
    with pytest.raises(ValueError):
        hamming_distance("A", "AC")
